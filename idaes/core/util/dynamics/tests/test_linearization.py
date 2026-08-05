#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
import pytest
from itertools import combinations

from numpy import array, eye, logspace, ndarray, zeros
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    exp,
    Param,
    Reference,
    Set,
    sin,
    Suffix,
    TransformationFactory,
    units as pyunits,
    value,
    Var,
)
from pyomo.common.collections import ComponentSet
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP

from idaes.core.util.constants import Constants
from idaes.core.util.dynamics.linearization import (
    _validate_vardata_collections,
    _validate_steady_state,
    c2d,
    linearize_system,
)
from idaes.core.util.exceptions import BurntToast

__author__ = "Douglas Allan"


def _component_index(lst, comp):
    """
    Function to emulate the function of list.index for collections of
    Pyomo components. Since Pyomo variables cannot be compared using ==,
    we check for component identity instead.

    Args:
        lst (iter): Iterable containing Pyomo components
        comp (ComponentData): Pyomo component which we want to find

    Returns:
        Index at which comp occurs in lst

    Raises:
        ValueError if comp does not occur in lst
    """
    for i, v in enumerate(lst):
        if v is comp:
            return i
    raise ValueError(f"Component {comp} is not present in collection.")


@pytest.mark.unit
class TestValidateVardataCollections(object):
    """
    Tests for the basic features of the linearize_system function.

    Tests in this file were prepared with the assistance of Google
    Gemini 3.5 Flash.
    """

    @pytest.fixture
    @classmethod
    def base_setup(cls):
        """Sets up a Pyomo model and base variables for partitioning."""
        m = ConcreteModel()
        # Create variables representing each distinct category
        m.deriv_v = Var()
        m.diff_v = Var()
        m.alg_v = Var()
        m.input_v = Var()
        m.dist_v = Var()

        # Setup baseline dictionaries
        vardata_lists = {
            "deriv": [m.deriv_v],
            "diff": [m.diff_v],
            "alg": [m.alg_v],
            "input": [m.input_v],
            "dist": [m.dist_v],
            "total": [m.deriv_v, m.diff_v, m.alg_v, m.input_v, m.dist_v],
            "output": [m.deriv_v, m.alg_v, m.dist_v],  # subset of total
        }

        vardata_sets = {key: ComponentSet(lst) for key, lst in vardata_lists.items()}

        return m, vardata_lists, vardata_sets

    def test_validation_success(self, base_setup):
        """Verifies valid configurations pass validation and output lists are correctly partitioned."""
        _, vardata_lists, vardata_sets = base_setup

        output_partition = _validate_vardata_collections(vardata_lists, vardata_sets)

        # Check that outputs are mapped and original list ordering is preserved
        assert output_partition["deriv"] == ComponentSet([base_setup[0].deriv_v])
        assert output_partition["alg"] == ComponentSet([base_setup[0].alg_v])
        assert output_partition["dist"] == ComponentSet([base_setup[0].dist_v])
        assert output_partition["diff"] == ComponentSet()
        assert output_partition["input"] == ComponentSet()

    def test_output_partition_assignment(self, base_setup):
        """Verifies that input and differential variables that show up as outputs
        are correctly assigned to their respective partitioned ComponentSets."""
        m, vardata_lists, vardata_sets = base_setup

        # Configure outputs to explicitly include a differential and an input variable
        # alongside an algebraic variable
        vardata_lists["output"] = [m.input_v, m.alg_v, m.diff_v]
        vardata_sets["output"] = ComponentSet(vardata_lists["output"])

        output_partition = _validate_vardata_collections(vardata_lists, vardata_sets)

        assert output_partition["input"] == ComponentSet([m.input_v])
        assert output_partition["alg"] == ComponentSet([m.alg_v])
        assert output_partition["diff"] == ComponentSet([m.diff_v])

        # Ensure other categories remain empty
        assert output_partition["deriv"] == ComponentSet()
        assert output_partition["dist"] == ComponentSet()

    @pytest.mark.parametrize(
        "category", ["total", "deriv", "diff", "alg", "input", "output", "dist"]
    )
    def test_duplicate_detection(self, base_setup, category):
        """Verifies ValueError is raised when duplicate variables are present in any category."""
        m, vardata_lists, vardata_sets = base_setup

        # Duplicate the first element in the target category list
        duplicated_var = vardata_lists[category][0]
        vardata_lists[category].append(duplicated_var)
        vardata_sets[category] = ComponentSet(
            vardata_lists[category]
        )  # sets inherently drop duplicates

        with pytest.raises(ValueError) as excinfo:
            _validate_vardata_collections(vardata_lists, vardata_sets)

        assert f"Duplicate variables found in the '{category}' variable list" in str(
            excinfo.value
        )

    # Generate all 6 unique pairs of disjoint categories:
    # ('diff', 'deriv'), ('diff', 'input'), ('diff', 'dist'), ('deriv', 'input'), ('deriv', 'dist'), ('input', 'dist')
    disjoint_pairs = list(combinations(["diff", "deriv", "input", "dist"], 2))

    @pytest.mark.parametrize("cat_a, cat_b", disjoint_pairs)
    def test_disjoint_categories(self, base_setup, cat_a, cat_b):
        """Verifies ValueError is raised when any pair of disjoint categories overlaps."""
        m, vardata_lists, vardata_sets = base_setup

        # Share a variable between cat_a and cat_b
        shared_var = vardata_lists[cat_a][0]
        vardata_lists[cat_b].append(shared_var)

        vardata_sets[cat_a] = ComponentSet(vardata_lists[cat_a])
        vardata_sets[cat_b] = ComponentSet(vardata_lists[cat_b])

        with pytest.raises(ValueError) as excinfo:
            _validate_vardata_collections(vardata_lists, vardata_sets)

        assert "are not disjoint" in str(excinfo.value)

    @pytest.mark.parametrize("category", ["input", "dist", "output"])
    def test_missing_from_total(self, base_setup, category):
        """Verifies ValueError is raised if input, dist, or output is missing from total."""
        m, vardata_lists, vardata_sets = base_setup

        # Create an external variable not present in the 'total' set
        m.external_v = Var()

        vardata_lists[category].append(m.external_v)
        vardata_sets[category] = ComponentSet(vardata_lists[category])

        with pytest.raises(ValueError) as excinfo:
            _validate_vardata_collections(vardata_lists, vardata_sets)

        assert "is not present in the total variables set" in str(excinfo.value)


@pytest.mark.unit
class TestSteadyStateValidation(object):
    @pytest.fixture
    @classmethod
    def pyomo_system(cls):
        """Fixture containing a real Pyomo model with initialized vars and constraints."""
        m = ConcreteModel()
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)

        m.v_ok = Var(initialize=1e-6)
        m.v_fail = Var(initialize=0.5)

        # Clean equality constraint (v_ok == 0) -> body is v_ok, lb = ub = 0
        m.con_ok = Constraint(expr=m.v_ok == 0)

        # Failing equality constraint (v_fail == 0) -> body is v_fail, lb = ub = 0
        m.con_fail = Constraint(expr=m.v_fail == 0)

        # Inequality constraint (0 <= v_ok <= 10) -> lb = 0, ub = 10
        m.con_ineq = Constraint(expr=(0.0, m.v_ok, 10.0))

        return m

    def test_steady_state_validation_success(self, pyomo_system):
        m = pyomo_system
        m.scaling_factor[m.v_ok] = 1.0
        m.scaling_factor[m.con_ok] = 1.0

        # Ensure only the passing constraint is checked
        m.con_fail.deactivate()
        m.con_ineq.deactivate()

        # Should run without error
        _validate_steady_state(
            t_block=m,
            deriv_var_set=ComponentSet([m.v_ok]),
            steady_state_derivative_tolerance=1e-5,
            constraint_tolerance=1e-5,
        )

    def test_steady_state_validation_derivative_failure(self, pyomo_system):
        m = pyomo_system
        m.scaling_factor[m.v_fail] = 2.0  # Scaled value = 1.0

        m.con_ok.deactivate()
        m.con_fail.deactivate()
        m.con_ineq.deactivate()

        with pytest.raises(ValueError) as exc_info:
            _validate_steady_state(
                t_block=m,
                deriv_var_set=ComponentSet([m.v_fail]),
                steady_state_derivative_tolerance=0.1,
                constraint_tolerance=1e-5,
            )

        assert "Nonzero derivative variables" in str(exc_info.value)
        assert "v_fail: 5.00e-01 (scaling factor=2.00e+00)" in str(exc_info.value)

    def test_steady_state_validation_scaled_derivative_success(self, pyomo_system):
        m = pyomo_system
        m.scaling_factor[m.v_fail] = 0.1  # Scaled value = 5e-2

        m.con_ok.deactivate()
        m.con_fail.deactivate()
        m.con_ineq.deactivate()

        _validate_steady_state(
            t_block=m,
            deriv_var_set=ComponentSet([m.v_fail]),
            steady_state_derivative_tolerance=0.07,
            constraint_tolerance=1e-5,
        )

    def test_steady_state_validation_scaled_derivative_failure(self, pyomo_system):
        m = pyomo_system
        m.scaling_factor[m.v_fail] = 3.0  # Scaled value = 1.5

        m.con_ok.deactivate()
        m.con_fail.deactivate()
        m.con_ineq.deactivate()

        with pytest.raises(ValueError) as exc_info:
            _validate_steady_state(
                t_block=m,
                deriv_var_set=ComponentSet([m.v_fail]),
                steady_state_derivative_tolerance=1,
                constraint_tolerance=1e-5,
            )

        assert "Nonzero derivative variables" in str(exc_info.value)
        assert "v_fail: 5.00e-01 (scaling factor=3.00e+00)" in str(exc_info.value)

    def test_steady_state_validation_scaled_constraint_success(self, pyomo_system):
        m = pyomo_system
        m.scaling_factor[m.con_fail] = 0.1  # Scaled residual = 0.05

        m.con_ok.deactivate()
        m.con_ineq.deactivate()

        _validate_steady_state(
            t_block=m,
            deriv_var_set=ComponentSet(),
            steady_state_derivative_tolerance=1e-5,
            constraint_tolerance=0.2,
        )

    def test_steady_state_validation_scaled_constraint_failure(self, pyomo_system):
        m = pyomo_system
        m.scaling_factor[m.con_fail] = 10.0  # Scaled residual = 5.0

        m.con_ok.deactivate()
        m.con_ineq.deactivate()

        with pytest.raises(ValueError) as exc_info:
            _validate_steady_state(
                t_block=m,
                deriv_var_set=ComponentSet(),
                steady_state_derivative_tolerance=1e-5,
                constraint_tolerance=0.8,
            )

        assert "Nonzero constraint residuals" in str(exc_info.value)
        assert "con_fail: 5.00e-01 (scaling factor=1.00e+01)" in str(exc_info.value)

    def test_steady_state_validation_burnt_toast_exception(self, pyomo_system):
        m = pyomo_system
        m.con_ok.deactivate()
        m.con_fail.deactivate()

        with pytest.raises(BurntToast) as exc_info:
            _validate_steady_state(
                t_block=m,
                deriv_var_set=ComponentSet(),
                steady_state_derivative_tolerance=1e-5,
                constraint_tolerance=1e-5,
            )

        assert "Encountered the inequality constraint con_ineq" in str(exc_info.value)


@pytest.mark.unit
class TestLinearizeSystemUnits(object):
    @pytest.fixture
    def scalar_dae_model(self):
        """Creates a minimal DAE model: dx/dt = -3*x + 2*u."""
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=2.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.u = Var(m.time, initialize=1.0)
        m.y = Var(m.time, initialize=0.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == -3.0 * b.x[t] + 2.0 * b.u[t]

        @m.Constraint(m.time)
        def output_eqn(b, t):
            return b.y[t] == 5.0 * b.x[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")
        return m

    def test_linearize_system_output_fields_and_types(self):
        """
        Verify that linearize_system returns a dictionary containing all
        of the specified output keys, with correct types and exact contents
        as defined by the updated function docstring.
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=0.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.u = Var(m.time, initialize=1.0)
        m.y = Var(m.time, initialize=0.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == -b.x[t] + b.u[t]

        @m.Constraint(m.time)
        def out_eqn(b, t):
            return b.y[t] == 2.0 * b.x[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        m.u.fix(1.0)
        for t in m.time:
            m.x[t].set_value(1.0)
            m.xdot[t].set_value(0.0)
            m.y[t].set_value(2.0)

        # Representative time t1 = 1
        t1 = 1
        out = linearize_system(
            m,
            m.time,
            representative_time=t1,
            input_vars=[m.u],
            output_vars=[m.y],
            scaled=False,
        )

        expected_keys = {
            "scaled_jac",
            "nlp",
            "diff_vars",
            "alg_vars",
            "input_vars",
            "disturbance_vars",
            "output_vars",
            "A",
            "B",
            "Bd",
            "C",
            "D",
            "Dd",
        }

        # 1. Assert all keys exist in the returned dictionary
        for key in expected_keys:
            assert (
                key in out
            ), f"Expected key '{key}' was not found in the output dictionary."

        # 2. Check types of the diagnostic keys

        assert isinstance(out["scaled_jac"], (coo_matrix, csc_matrix, csr_matrix))
        assert isinstance(out["nlp"], PyomoNLP)

        # 3. Check types and specific identity of the variable lists
        assert isinstance(out["diff_vars"], list)
        assert isinstance(out["alg_vars"], list)
        assert isinstance(out["input_vars"], list)
        assert isinstance(out["disturbance_vars"], list)
        assert isinstance(out["output_vars"], list)

        # Pyomo VarData objects must be compared using 'is' via our helper
        assert len(out["diff_vars"]) == 1
        assert _component_index(out["diff_vars"], m.x[t1]) == 0

        assert len(out["input_vars"]) == 1
        assert _component_index(out["input_vars"], m.u[t1]) == 0

        assert len(out["output_vars"]) == 1
        assert _component_index(out["output_vars"], m.y[t1]) == 0

        # m.y[t1] is the algebraic variable in this system
        assert len(out["alg_vars"]) == 1
        assert _component_index(out["alg_vars"], m.y[t1]) == 0

        assert len(out["disturbance_vars"]) == 0

        # 4. Check types and dimensions of the output matrices
        assert isinstance(out["A"], ndarray)
        assert isinstance(out["B"], ndarray)
        assert isinstance(out["Bd"], ndarray)
        assert isinstance(out["C"], ndarray)
        assert isinstance(out["D"], ndarray)
        assert isinstance(out["Dd"], ndarray)

        assert out["A"].shape == (1, 1)
        assert out["B"].shape == (1, 1)
        assert out["Bd"].shape == (1, 0)
        assert out["C"].shape == (1, 1)
        assert out["D"].shape == (1, 1)
        assert out["Dd"].shape == (1, 0)

    def test_linearize_system_standard(self, scalar_dae_model):
        """Test basic successful linearization, unscaled output, and matrix sizes."""
        m = scalar_dae_model
        # Ensure the system is solved / at steady state (xdot = -3*2 + 2*3 = 0 -> u = 3)
        m.u.fix(3.0)
        for t in m.time:
            m.x[t].set_value(2.0)
            m.xdot[t].set_value(0.0)
            m.y[t].set_value(10.0)

        out = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.y],
            scaled=False,
        )

        assert out["A"] == pytest.approx(array([[-3.0]]))
        assert out["B"] == pytest.approx(array([[2.0]]))
        assert out["C"] == pytest.approx(array([[5.0]]))
        assert out["D"] == pytest.approx(array([[0.0]]))
        assert len(out["diff_vars"]) == 1
        assert len(out["alg_vars"]) == 1  # m.y is algebraic

    def test_linearize_system_restores_fixedness(self, scalar_dae_model):
        """Verify that variable fixedness is perfectly restored even if exceptions occur."""
        m = scalar_dae_model
        # Ensure the system is solved / at steady state (xdot = -3*2 + 2*3 = 0 -> u = 3)
        m.u.fix(3.0)
        for t in m.time:
            m.x[t].set_value(2.0)
            m.xdot[t].set_value(0.0)
            m.y[t].set_value(10.0)

        m.x[1].unfix()

        # Check pre-fixedness
        assert m.u[1].fixed is True
        assert m.x[1].fixed is False

        # Run linearization
        linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.y],
        )

        # Check post-fixedness (must remain unchanged)
        assert m.u[1].fixed is True
        assert m.x[1].fixed is False

    def test_linearize_system_non_square_dof_failure(self, scalar_dae_model):
        """Verify RuntimeError is raised if the timestep problem is not square (dof != 0)."""
        m = scalar_dae_model
        m.u.unfix()  # Leaves an extra degree of freedom

        with pytest.raises(RuntimeError) as exc_info:
            linearize_system(
                m,
                m.time,
                representative_time=1,
                input_vars=[],  # Exclude u so it doesn't get factored into square calculation
                output_vars=[m.y],
            )
        assert "degrees of freedom" in str(exc_info.value)

    def test_linearize_system_invalid_variable_category(self, scalar_dae_model):
        """Verify ValueError is raised if a variable is missing from the total set."""
        m = scalar_dae_model
        # Create an external variable completely disconnected from the active constraints
        m.external_var = Var(m.time, initialize=1.0)

        with pytest.raises(ValueError) as exc_info:
            linearize_system(
                m,
                m.time,
                representative_time=1,
                input_vars=[m.external_var],
                output_vars=[m.y],
            )
        assert "is not present in the total variables set" in str(exc_info.value)

    def test_linearize_system_nonzero_D_matrix(self):
        """Test a system where there is a direct feedthrough from u to y.
        dx/dt = -x + u
        y = 3*x + 2*u
        Should yield: A = [[-1]], B = [[1]], C = [[3]], D = [[2]]
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=1.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.u = Var(m.time, initialize=1.0)
        m.y = Var(m.time, initialize=5.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == -b.x[t] + b.u[t]

        @m.Constraint(m.time)
        def output_eqn(b, t):
            return b.y[t] == 3.0 * b.x[t] + 2.0 * b.u[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        m.u.fix(1.0)
        for t in m.time:
            m.x[t].set_value(1.0)
            m.xdot[t].set_value(0.0)
            m.y[t].set_value(5.0)

        out = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.y],
            scaled=False,
        )

        assert out["A"] == pytest.approx(array([[-1.0]]))
        assert out["B"] == pytest.approx(array([[1.0]]))
        assert out["C"] == pytest.approx(array([[3.0]]))
        assert out["D"] == pytest.approx(array([[2.0]]))
        assert out["Bd"].shape == (1, 0)
        assert out["Dd"].shape == (1, 0)

    def test_linearize_system_with_disturbances(self):
        """Test a system with inputs, disturbances, and outputs.
        dx/dt = -2*x + u - 3*d
        y = x + d
        Should yield: A = [[-2]], B = [[1]], Bd = [[-3]], C = [[1]], D = [[0]], Dd = [[1]]
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=2.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.u = Var(m.time, initialize=7.0)
        m.d = Var(m.time, initialize=1.0)
        m.y = Var(m.time, initialize=3.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == -2.0 * b.x[t] + b.u[t] - 3.0 * b.d[t]

        @m.Constraint(m.time)
        def output_eqn(b, t):
            return b.y[t] == b.x[t] + b.d[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        m.u.fix(7.0)
        m.d.fix(1.0)
        for t in m.time:
            m.x[t].set_value(2.0)
            m.xdot[t].set_value(0.0)
            m.y[t].set_value(3.0)

        out = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            disturbance_vars=[m.d],
            output_vars=[m.y],
            scaled=False,
        )

        assert out["A"] == pytest.approx(array([[-2.0]]))
        assert out["B"] == pytest.approx(array([[1.0]]))
        assert out["Bd"] == pytest.approx(array([[-3.0]]))
        assert out["C"] == pytest.approx(array([[1.0]]))
        assert out["D"] == pytest.approx(array([[0.0]]))
        assert out["Dd"] == pytest.approx(array([[1.0]]))

    def test_linearize_system_nontrivial_mass_matrix(self):
        """Test a coupled system where the derivatives are implicit (mass matrix M != I).
        2*x1dot + x2dot = -x1 + u  ==> eq1
        x1dot + 2*x2dot = -x2      ==> eq2

        Rewritten explicitly (M_inv * rhs):
        [x1dot] = [ 2  1 ]^-1 * [-x1 + u] = [ 2/3  -1/3 ] * [-x1 + u] = -2/3*x1 + 1/3*x2 + 2/3*u
        [x2dot]   [ 1  2 ]      [-x2    ]   [-1/3   2/3 ]   [-x2    ]   1/3*x1 - 2/3*x2 - 1/3*u

        Our A matrix should be:
        A = [[ -2/3,  1/3 ],
             [  1/3, -2/3 ]]
        B = [[  2/3 ],
             [ -1/3 ]]
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x1 = Var(m.time, initialize=1.0)
        m.x2 = Var(m.time, initialize=-1.0)
        m.x1dot = DerivativeVar(m.x1, wrt=m.time, initialize=0.0)
        m.x2dot = DerivativeVar(m.x2, wrt=m.time, initialize=0.0)
        m.u = Var(m.time, initialize=1.0)

        @m.Constraint(m.time)
        def ode_eqn1(b, t):
            return 2.0 * b.x1dot[t] + b.x2dot[t] == -b.x1[t] + b.u[t]

        @m.Constraint(m.time)
        def ode_eqn2(b, t):
            return b.x1dot[t] + 2.0 * b.x2dot[t] == -b.x2[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        m.u.fix(1.0)
        for t in m.time:
            m.x1[t].set_value(1.0)
            m.x2[t].set_value(0)
            m.x1dot[t].set_value(0.0)
            m.x2dot[t].set_value(0.0)

        out = linearize_system(
            m, m.time, representative_time=1, input_vars=[m.u], scaled=False
        )

        expected_A = array([[-2.0 / 3.0, 1.0 / 3.0], [1.0 / 3.0, -2.0 / 3.0]])
        expected_B = array([[2.0 / 3.0], [-1.0 / 3.0]])

        # Align indices of diff vars list with expected matrix orientation
        assert len(out["diff_vars"]) == 2
        idx1 = _component_index(out["diff_vars"], m.x1[1])
        idx2 = _component_index(out["diff_vars"], m.x2[1])

        # map actual A & B back to match [x1, x2] ordering
        A_actual = out["A"][[[idx1], [idx2]], [idx1, idx2]]
        B_actual = out["B"][[idx1, idx2], :]

        assert A_actual == pytest.approx(expected_A)
        assert B_actual == pytest.approx(expected_B)

    def test_linearize_system_preserves_ordering(self):
        """
        Verify that the order of the inputs and disturbances provided to
        linearize_system is perfectly preserved in the shape and positioning
        of columns in B, Bd, D, and Dd.
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=0.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)

        # Multiple inputs, disturbances, and outputs to test ordering permutations
        m.u1 = Var(m.time, initialize=0)
        m.u2 = Var(m.time, initialize=0)
        m.d1 = Var(m.time, initialize=0)
        m.d2 = Var(m.time, initialize=0)

        m.y1 = Var(m.time, initialize=0.0)
        m.y2 = Var(m.time, initialize=0.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            # dx/dt = -x + 10*u1 + 20*u2 - 30*d1 - 40*d2
            return (
                b.xdot[t]
                == -b.x[t]
                + 10.0 * b.u1[t]
                + 20.0 * b.u2[t]
                - 30.0 * b.d1[t]
                - 40.0 * b.d2[t]
            )

        @m.Constraint(m.time)
        def out_eqn1(b, t):
            # y1 = 2*u1 + 3*u2
            return b.y1[t] == 2.0 * b.u1[t] + 3.0 * b.u2[t]

        @m.Constraint(m.time)
        def out_eqn2(b, t):
            # y2 = 4*d1 + 5*d2
            return b.y2[t] == 4.0 * b.d1[t] + 5.0 * b.d2[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        # -- Permutation A --
        out_a = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u1, m.u2],
            disturbance_vars=[m.d1, m.d2],
            output_vars=[m.y1, m.y2],
            scaled=False,
        )

        # Output vars ordering check
        assert _component_index(out_a["output_vars"], m.y1[1]) == 0
        assert _component_index(out_a["output_vars"], m.y2[1]) == 1

        # Matrix columns mapped to: [u1, u2] and [d1, d2]
        # B = [10, 20]
        # Bd = [-30, -40]
        # D (rows: y1, y2; cols: u1, u2) = [[2, 3], [0, 0]]
        # Dd (rows: y1, y2; cols: d1, d2) = [[0, 0], [4, 5]]
        assert array(out_a["B"]) == pytest.approx(array([[10.0, 20.0]]))
        assert array(out_a["Bd"]) == pytest.approx(array([[-30.0, -40.0]]))
        assert array(out_a["D"]) == pytest.approx(array([[2.0, 3.0], [0.0, 0.0]]))
        assert array(out_a["Dd"]) == pytest.approx(array([[0.0, 0.0], [4.0, 5.0]]))

        # -- Permutation B (Reversed Lists) --
        out_b = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u2, m.u1],
            disturbance_vars=[m.d2, m.d1],
            output_vars=[m.y2, m.y1],
            scaled=False,
        )

        # Output vars ordering check
        assert _component_index(out_b["output_vars"], m.y2[1]) == 0
        assert _component_index(out_b["output_vars"], m.y1[1]) == 1

        # Matrix columns mapped to: [u2, u1] and [d2, d1]
        # B = [20, 10]
        # Bd = [-40, -30]
        # D (rows: y2, y1; cols: u2, u1) = [[0, 0], [3, 2]]
        # Dd (rows: y2, y1; cols: d2, d1) = [[5, 4], [0, 0]]
        assert array(out_b["B"]) == pytest.approx(array([[20.0, 10.0]]))
        assert array(out_b["Bd"]) == pytest.approx(array([[-40.0, -30.0]]))
        assert array(out_b["D"]) == pytest.approx(array([[0.0, 0.0], [3.0, 2.0]]))
        assert array(out_b["Dd"]) == pytest.approx(array([[5.0, 4.0], [0.0, 0.0]]))

    def test_linearize_system_representative_time(self):
        """
        Verify that the representative_time argument correctly dictates which
        time index is used to build the timestep problem, extract variables,
        and evaluate the linearized matrices.

        We set up a system with a time-varying parameter (p) that changes over time:
        dx/dt = -p[t]*x + u
        At t=1, p[1] = 1.0  => dx/dt = -1.0*x + u  => A = [[-1.0]]
        At t=2, p[2] = 5.0  => dx/dt = -5.0*x + u  => A = [[-5.0]]
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1, 2])
        m.x = Var(m.time, initialize=1.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.u = Var(m.time, initialize=1.0)
        m.p = Param(m.time, initialize={0: 1.0, 1: 1.0, 2: 5.0}, mutable=True)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == -b.p[t] * b.x[t] + b.u[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=2, scheme="BACKWARD")

        m.u.fix(1.0)

        # Establish operating steady-state for t=1 (p=1.0, x=1.0, u=1.0)
        m.x[1].set_value(1.0)
        m.xdot[1].set_value(0.0)

        # Establish operating steady-state for t=2 (p=5.0, x=0.2, u=1.0)
        m.x[2].set_value(0.2)
        m.xdot[2].set_value(0.0)

        # --- Case A: representative_time = 1 ---
        out_t1 = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.x],
            scaled=False,
        )

        # Assert correct variable slice is extracted (t = 1)
        assert len(out_t1["diff_vars"]) == 1
        assert _component_index(out_t1["diff_vars"], m.x[1]) == 0
        assert _component_index(out_t1["input_vars"], m.u[1]) == 0
        assert _component_index(out_t1["output_vars"], m.x[1]) == 0

        # Assert A matrix reflects p[1] = 1.0
        assert array(out_t1["A"]) == pytest.approx(array([[-1.0]]))

        # --- Case B: representative_time = 2 ---
        out_t2 = linearize_system(
            m,
            m.time,
            representative_time=2,
            input_vars=[m.u],
            output_vars=[m.x],
            scaled=False,
        )

        # Assert correct variable slice is extracted (t = 2)
        assert len(out_t2["diff_vars"]) == 1
        assert _component_index(out_t2["diff_vars"], m.x[2]) == 0
        assert _component_index(out_t2["input_vars"], m.u[2]) == 0
        assert _component_index(out_t2["output_vars"], m.x[2]) == 0

        # Assert A matrix reflects p[2] = 5.0
        assert array(out_t2["A"]) == pytest.approx(array([[-5.0]]))

    def test_linearize_system_implicit_representative_time(self):
        """
        Verify that omitting representative_time uses the second index of the
        Time ContinuousSet (time.at(2)) automatically.
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0.0, 10.0, 20.0])  # time.at(2) will be 10.0
        m.x = Var(m.time, initialize=0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0)
        m.p = Param(m.time, initialize={0: -1, 10: -2, 20: -3})

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == b.p[t] * b.x[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        # Omit representative_time
        out = linearize_system(
            m, m.time, input_vars=[], output_vars=[m.x], scaled=False
        )

        # Confirm variables and equations evaluated at t = 10.0
        assert _component_index(out["diff_vars"], m.x[10.0]) == 0
        assert _component_index(out["output_vars"], m.x[10.0]) == 0
        assert array(out["A"]) == pytest.approx(array([[-2]]))

    def test_linearize_system_scaled_behavior(self):
        """
        Verify that when scaled=True, variables are adjusted by their scaling factors,
        and when scaled=False, the returned matrices are correctly unscaled.

        Let:
        dx/dt = -3*x + 5*u
        y = 4*x
        Unscaled continuous system matrices (scaled=False):
        A = [[-3.0]], B = [[5.0]], C = [[4.0]], D = [[0.0]]

        Let non-unity scaling factors be:
        s_x = 10.0       (scaling factor for state x)
        s_xdot = 2.0     (scaling factor for derivative xdot)
        s_u = 0.5        (scaling factor for input u)
        s_y = 8.0        (scaling factor for output y)

        Expected Scaled matrices (scaled=True):
        A_scaled = diag(s_xdot) @ A @ diag(1/s_x)    = [[ 2.0 * (-3.0) * (1 / 10.0) ]] = [[-0.6]]
        B_scaled = diag(s_xdot) @ B @ diag(1/s_u)    = [[ 2.0 * (5.0)  * (1 / 0.5)  ]] = [[20.0]]
        C_scaled = diag(s_y)    @ C @ diag(1/s_x)    = [[ 8.0 * (4.0)  * (1 / 10.0) ]] = [[ 3.2]]
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=1.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.u = Var(m.time, initialize=1.0)
        m.y = Var(m.time, initialize=4.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == -3.0 * b.x[t] + 5.0 * b.u[t]

        @m.Constraint(m.time)
        def output_eqn(b, t):
            return b.y[t] == 4.0 * b.x[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        m.u.fix(3.0)
        for t in m.time:
            m.x[t].set_value(5.0)
            m.xdot[t].set_value(0.0)
            m.y[t].set_value(20.0)

        # Declare scaling factors in the model's Export suffix
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        m.scaling_factor[m.x[1]] = 10.0
        m.scaling_factor[m.xdot[1]] = 2.0
        m.scaling_factor[m.u[1]] = 0.5
        m.scaling_factor[m.y[1]] = 8.0

        # --- Case A: scaled = False (returns physical matrices) ---
        out_unscaled = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.y],
            scaled=False,
        )

        assert array(out_unscaled["A"]) == pytest.approx(array([[-3.0]]))
        assert array(out_unscaled["B"]) == pytest.approx(array([[5.0]]))
        assert array(out_unscaled["C"]) == pytest.approx(array([[4.0]]))
        assert array(out_unscaled["D"]) == pytest.approx(array([[0.0]]))

        # --- Case B: scaled = True (returns mathematically scaled matrices) ---
        out_scaled = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.y],
            scaled=True,
        )

        assert array(out_scaled["A"]) == pytest.approx(array([[-0.6]]))
        assert array(out_scaled["B"]) == pytest.approx(array([[20.0]]))
        assert array(out_scaled["C"]) == pytest.approx(array([[3.2]]))
        assert array(out_scaled["D"]) == pytest.approx(array([[0.0]]))

    def test_linearize_system_constraint_scaling_invariance(self):
        """
        Verify that constraint scaling factors DO NOT change the resulting unscaled
        state-space matrices (A, B, C, D).

        System definition:
        dx/dt = -3*x + 5*u    ==> Equation 1 (ode_eqn)
        y = 2*x               ==> Equation 2 (output_eqn)

        Unscaled continuous system matrices (A = [[-3.0]], B = [[5.0]], C = [[2.0]], D = [[0.0]])
        should remain exactly the same even if Equation 1 is scaled by 100.0 and Equation 2
        is scaled by 0.01.
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=1.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.u = Var(m.time, initialize=1.0)
        m.y = Var(m.time, initialize=2.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == -3.0 * b.x[t] + 5.0 * b.u[t]

        @m.Constraint(m.time)
        def output_eqn(b, t):
            return b.y[t] == 2.0 * b.x[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        m.u.fix(3.0)
        for t in m.time:
            m.x[t].set_value(5.0)
            m.xdot[t].set_value(0.0)
            m.y[t].set_value(10.0)

        # 1. Base run without any constraint scaling factors
        out_base = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.y],
            scaled=False,
        )

        assert array(out_base["A"]) == pytest.approx(array([[-3.0]]))
        assert array(out_base["B"]) == pytest.approx(array([[5.0]]))
        assert array(out_base["C"]) == pytest.approx(array([[2.0]]))
        assert array(out_base["D"]) == pytest.approx(array([[0.0]]))

        # 2. Add highly non-unity scaling factors to the constraints
        m.scaling_factor = Suffix(direction=Suffix.EXPORT)
        # Apply scaling factors to the specific constraint datas evaluated at t=1
        m.scaling_factor[m.ode_eqn[1]] = 100.0
        m.scaling_factor[m.output_eqn[1]] = 0.01

        out_scaled_constraints = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.y],
            scaled=False,
        )

        # The unscaled matrices must remain exactly identical to the base run
        assert array(out_scaled_constraints["A"]) == pytest.approx(array(out_base["A"]))
        assert array(out_scaled_constraints["B"]) == pytest.approx(array(out_base["B"]))
        assert array(out_scaled_constraints["C"]) == pytest.approx(array(out_base["C"]))
        assert array(out_scaled_constraints["D"]) == pytest.approx(array(out_base["D"]))

        # 3. Finally, the results should be unaffected when scaled=True because constraint
        # scaling factors are cancelled when inverting the Jacobian
        out_scaled_constraints = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            output_vars=[m.y],
            scaled=True,
        )

        # The unscaled matrices must remain exactly identical to the base run
        assert array(out_scaled_constraints["A"]) == pytest.approx(array(out_base["A"]))
        assert array(out_scaled_constraints["B"]) == pytest.approx(array(out_base["B"]))
        assert array(out_scaled_constraints["C"]) == pytest.approx(array(out_base["C"]))
        assert array(out_scaled_constraints["D"]) == pytest.approx(array(out_base["D"]))

    def test_linearize_system_output_categories(self):
        """
        Verify that output matrices (C, D, Dd) are correctly constructed for
        each of the five possible mathematical variable categories:
        1) Differential (x)
        2) Algebraic (z)
        3) Derivative (xdot)
        4) Input (u)
        5) Disturbance (d)

        System equations:
        dx/dt = -3*x + 5*u - 2*d   (ode_eqn)
        z = 4*x + 6*u + 7*d        (alg_eqn)
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=1.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.z = Var(m.time, initialize=0.0)
        m.u = Var(m.time, initialize=1.0)
        m.d = Var(m.time, initialize=1.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == -3.0 * b.x[t] + 5.0 * b.u[t] - 2.0 * b.d[t]

        @m.Constraint(m.time)
        def alg_eqn(b, t):
            return b.z[t] == 4.0 * b.x[t] + 6.0 * b.u[t] + 7.0 * b.d[t]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        m.u.fix(1.0)
        m.d.fix(1.0)
        for t in m.time:
            m.x[t].set_value(1.0)
            m.xdot[t].set_value(0.0)
            m.z[t].set_value(17.0)

        # We pass one variable of each category as an output variable
        out = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u],
            disturbance_vars=[m.d],
            output_vars=[m.x, m.z, m.xdot, m.u, m.d],
            scaled=False,
        )

        # Output list ordering: [x, z, xdot, u, d]
        assert len(out["output_vars"]) == 5
        assert _component_index(out["output_vars"], m.x[1]) == 0
        assert _component_index(out["output_vars"], m.z[1]) == 1
        assert _component_index(out["output_vars"], m.xdot[1]) == 2
        assert _component_index(out["output_vars"], m.u[1]) == 3
        assert _component_index(out["output_vars"], m.d[1]) == 4

        # 1) Differential Output (index 0)
        # y = x  => C[0, :] = [[1]], D[0, :] = [[0]], Dd[0, :] = [[0]]
        assert array(out["C"])[0, :] == pytest.approx(array([1.0]))
        assert array(out["D"])[0, :] == pytest.approx(array([0.0]))
        assert array(out["Dd"])[0, :] == pytest.approx(array([0.0]))

        # 2) Algebraic Output (index 1)
        # y = 4*x + 6*u + 7*d  => C[1, :] = [[4]], D[1, :] = [[6]], Dd[1, :] = [[7]]
        assert array(out["C"])[1, :] == pytest.approx(array([4.0]))
        assert array(out["D"])[1, :] == pytest.approx(array([6.0]))
        assert array(out["Dd"])[1, :] == pytest.approx(array([7.0]))

        # 3) Derivative Output (index 2)
        # y = xdot = -3*x + 5*u - 2*d  => C[2, :] = [[-3.0]], D[2, :] = [[5.0]], Dd[2, :] = [[-2.0]]
        assert array(out["C"])[2, :] == pytest.approx(array([-3.0]))
        assert array(out["D"])[2, :] == pytest.approx(array([5.0]))
        assert array(out["Dd"])[2, :] == pytest.approx(array([-2.0]))

        # 4) Input Output (index 3)
        # y = u  => C[3, :] = [[0]], D[3, :] = [[1]], Dd[3, :] = [[0]]
        assert array(out["C"])[3, :] == pytest.approx(array([0.0]))
        assert array(out["D"])[3, :] == pytest.approx(array([1.0]))
        assert array(out["Dd"])[3, :] == pytest.approx(array([0.0]))

        # 5) Disturbance Output (index 4)
        # y = d  => C[4, :] = [[0]], D[4, :] = [[0]], Dd[4, :] = [[1]]
        assert array(out["C"])[4, :] == pytest.approx(array([0.0]))
        assert array(out["D"])[4, :] == pytest.approx(array([0.0]))
        assert array(out["Dd"])[4, :] == pytest.approx(array([1.0]))

    def test_linearize_system_output_categories_multivariable(self):
        """
        Verify that output matrices (C, D, Dd) are correctly constructed for
        each of the five possible mathematical variable categories when there are
        multiple input and disturbance variables to verify column mapping.

        System equations:
        dx/dt = -3*x + 5*u1 + 10*u2 - 2*d1 - 8*d2   (ode_eqn)
        z = 4*x + 6*u1 + 12*u2 + 7*d1 + 14*d2       (alg_eqn)

        To enforce steady state (dx/dt = 0):
        3*x = 5(1) + 10(1.5) - 2(1) - 8(2) = 2
        x = 2/3

        z = 4(2/3) + 6(1) + 12(1.5) + 7(1) + 14(2) = 8/3 + 59 = 185/3
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.x = Var(m.time, initialize=2.0 / 3.0)
        m.xdot = DerivativeVar(m.x, wrt=m.time, initialize=0.0)
        m.z = Var(m.time, initialize=185.0 / 3.0)
        m.u1 = Var(m.time, initialize=1.0)
        m.u2 = Var(m.time, initialize=1.5)
        m.d1 = Var(m.time, initialize=1.0)
        m.d2 = Var(m.time, initialize=2.0)

        @m.Constraint(m.time)
        def ode_eqn(b, t):
            return b.xdot[t] == (
                -3.0 * b.x[t]
                + 5.0 * b.u1[t]
                + 10.0 * b.u2[t]
                - 2.0 * b.d1[t]
                - 8.0 * b.d2[t]
            )

        @m.Constraint(m.time)
        def alg_eqn(b, t):
            return b.z[t] == (
                4.0 * b.x[t]
                + 6.0 * b.u1[t]
                + 12.0 * b.u2[t]
                + 7.0 * b.d1[t]
                + 14.0 * b.d2[t]
            )

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        m.u1.fix(1.0)
        m.u2.fix(1.5)
        m.d1.fix(1.0)
        m.d2.fix(2.0)
        for t in m.time:
            m.x[t].set_value(2.0 / 3.0)
            m.xdot[t].set_value(0.0)
            m.z[t].set_value(185.0 / 3.0)

        # We pass one variable of each category as an output variable
        # Inputs: [u1, u2]
        # Disturbances: [d1, d2]
        out = linearize_system(
            m,
            m.time,
            representative_time=1,
            input_vars=[m.u1, m.u2],
            disturbance_vars=[m.d1, m.d2],
            output_vars=[m.x, m.z, m.xdot, m.u2, m.d1],
            scaled=False,
        )

        # Output list ordering: [x, z, xdot, u2, d1]
        assert len(out["output_vars"]) == 5
        assert _component_index(out["output_vars"], m.x[1]) == 0
        assert _component_index(out["output_vars"], m.z[1]) == 1
        assert _component_index(out["output_vars"], m.xdot[1]) == 2
        assert _component_index(out["output_vars"], m.u2[1]) == 3
        assert _component_index(out["output_vars"], m.d1[1]) == 4

        # 1) Differential Output (index 0)
        # y = x  => C = [[1]], D = [[0, 0]], Dd = [[0, 0]]
        assert array(out["C"])[0, :] == pytest.approx(array([1.0]))
        assert array(out["D"])[0, :] == pytest.approx(array([0.0, 0.0]))
        assert array(out["Dd"])[0, :] == pytest.approx(array([0.0, 0.0]))

        # 2) Algebraic Output (index 1)
        # y = 4*x + 6*u1 + 12*u2 + 7*d1 + 14*d2
        # C = [[4]], D = [[6, 12]], Dd = [[7, 14]]
        assert array(out["C"])[1, :] == pytest.approx(array([4.0]))
        assert array(out["D"])[1, :] == pytest.approx(array([6.0, 12.0]))
        assert array(out["Dd"])[1, :] == pytest.approx(array([7.0, 14.0]))

        # 3) Derivative Output (index 2)
        # y = xdot = -3*x + 5*u1 + 10*u2 - 2*d1 - 8*d2
        # C = [[-3]], D = [[5, 10]], Dd = [[-2, -8]]
        assert array(out["C"])[2, :] == pytest.approx(array([-3.0]))
        assert array(out["D"])[2, :] == pytest.approx(array([5.0, 10.0]))
        assert array(out["Dd"])[2, :] == pytest.approx(array([-2.0, -8.0]))

        # 4) Input Output (index 3 checking u2)
        # y = u2  => C = [[0]], D = [[0, 1]] (since u2 is the 2nd input), Dd = [[0, 0]]
        assert array(out["C"])[3, :] == pytest.approx(array([0.0]))
        assert array(out["D"])[3, :] == pytest.approx(array([0.0, 1.0]))
        assert array(out["Dd"])[3, :] == pytest.approx(array([0.0, 0.0]))

        # 5) Disturbance Output (index 4 checking d1)
        # y = d1  => C = [[0]], D = [[0, 0]], Dd = [[1, 0]] (since d1 is the 1st disturbance)
        assert array(out["C"])[4, :] == pytest.approx(array([0.0]))
        assert array(out["D"])[4, :] == pytest.approx(array([0.0, 0.0]))
        assert array(out["Dd"])[4, :] == pytest.approx(array([1.0, 0.0]))

    def test_linearize_system_with_sliced_references(self):
        """
        Verify that linearize_system works transparently when passed Reference slices
        to represent variables indexed by sets in addition to the time continuous set.
        """
        m = ConcreteModel()
        m.time = ContinuousSet(initialize=[0, 1])
        m.space = Set(initialize=["inlet", "outlet"])

        # A state variable indexed by both time AND space
        m.x_spatial = Var(m.time, m.space, initialize=0)
        m.xdot_spatial = DerivativeVar(m.x_spatial, wrt=m.time, initialize=0)

        @m.Constraint(m.time, m.space)
        def ode_eqn(b, t, s):
            # dx/dt = -2*x
            return b.xdot_spatial[t, s] == -2.0 * b.x_spatial[t, s]

        fd_factory = TransformationFactory("dae.finite_difference")
        fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

        # Create Reference slices indexed only by time
        x_outlet_ref = Reference(m.x_spatial[:, "outlet"])

        out = linearize_system(
            m, m.time, representative_time=1, output_vars=[x_outlet_ref], scaled=False
        )

        assert len(out["diff_vars"]) == 2  # Spatial states at 'inlet' and 'outlet'
        assert len(out["output_vars"]) == 1
        assert _component_index(out["output_vars"], m.x_spatial[1, "outlet"]) == 0


@pytest.mark.unit
class TestC2dConversion(object):
    def test_c2d_conversion(self):
        """
        Verify the continuous-to-discrete LTI system conversion function 'c2d'
        with states, inputs, and disturbances.
        """
        continuous_sys = {
            "A": array([[-2.0]]),
            "B": array([[1.0]]),
            "Bd": array([[3.0]]),
        }

        dt = 0.5
        discrete_sys = c2d(continuous_sys, dt=dt)

        expected_A = array([[0.36787944]])
        expected_B = array([[0.31606028]])
        expected_Bd = array([[0.94818083]])

        assert discrete_sys["A"] == pytest.approx(expected_A, rel=1e-5)
        assert discrete_sys["B"] == pytest.approx(expected_B, rel=1e-5)
        assert discrete_sys["Bd"] == pytest.approx(expected_Bd, rel=1e-5)


@pytest.mark.component
def test_ogunnaike_ray_example_10_4():
    """
    Example 10.4 from Ogunnaike and Ray (1994) for level control
     conical tank. The cone's point faces downwards (like
    an ice cream cone), so the cone's "base" is at its top.
    """
    m = ConcreteModel()
    m.time = ContinuousSet(initialize=[0, 1])

    m.h = Var(m.time, initialize=1, doc="Tank level")
    m.alpha = Param(initialize=1, mutable=True, doc="Cone shape parameter")
    m.c = Param(initialize=1, mutable=True, doc="Cone outflow parameter")
    m.beta = Var(initialize=1, doc="Cone outflow parameter")
    m.hdot = DerivativeVar(m.h, wrt=m.time, initialize=0)

    m.Fi = Var(m.time, initialize=1, doc="Inlet flowrate")

    @m.Constraint()
    def beta_eqn(b):
        # Why do we need both beta and c? Ask the authors.
        return b.beta == b.alpha * b.c

    @m.Constraint(m.time)
    def hdot_eqn(b, t):
        return b.hdot[t] == b.alpha * b.Fi[t] / b.h[t] ** 2 - b.beta * b.h[t] ** -1.5

    fd_factory = TransformationFactory("dae.finite_difference")
    fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

    for alpha in logspace(-1, 1, 3):
        for c in logspace(-1, 1, 3):
            for Fi in logspace(-1, 1, 3):
                m.alpha.set_value(alpha)
                m.c.set_value(c)
                m.Fi.fix(Fi)
                m.beta.set_value(c * alpha)
                for t in m.time:
                    m.h[t].set_value((Fi / c) ** 2)

                out = linearize_system(
                    m,
                    m.time,
                    representative_time=1,
                    disturbance_vars=[m.Fi],
                    output_vars=[m.h],
                )

                K = value(2 / c * m.h[1] ** 0.5)
                tau = value(2 / m.beta * m.h[1] ** 2.5)

                assert out["A"] == pytest.approx(-1 / tau)
                assert out["B"].shape == (1, 0)
                assert out["Bd"] == pytest.approx(K / tau)
                assert out["C"] == pytest.approx(1)
                assert out["D"].shape == (1, 0)
                assert out["Dd"] == pytest.approx(0)


@pytest.mark.component
def test_rawlings_mayne_example():
    """
    Example 1.11 from Rawlings, Mayne, and Diehl (2024) involving a
    nonisothermal CSTR with varying level and a cooling jacket.

    Rawlings, J. B., Mayne, D. Q., Diehl, M. (2024).
    Model Predictive Control: Theory,Computation, and Design.
    United States: Nob Hill Publishing, LLC.
    """
    m = ConcreteModel()
    m.time = ContinuousSet(initialize=[0, 1])

    # State variables
    m.c = Var(
        m.time,
        initialize=0.878,
        units=pyunits.kmol / pyunits.m**3,
        doc="Tank concentration",
    )
    m.T = Var(m.time, initialize=324.5, units=pyunits.K, doc="Tank temperature")
    m.h = Var(m.time, initialize=0.659, units=pyunits.m, doc="Fluid height in tank")

    # Derivative variables
    m.cdot = DerivativeVar(
        m.c, wrt=m.time, initialize=0, units=pyunits.kmol / pyunits.m**3 / pyunits.min
    )
    m.Tdot = DerivativeVar(m.T, wrt=m.time, initialize=0, units=pyunits.K / pyunits.min)
    m.hdot = DerivativeVar(m.h, wrt=m.time, initialize=0, units=pyunits.m / pyunits.m)

    # Control variables
    m.Tc = Var(
        m.time, initialize=300, units=pyunits.K, doc="Cooling jacket temperature"
    )
    m.F = Var(
        m.time, initialize=0.1, units=pyunits.m**3 / pyunits.min, doc="Outlet flowrate"
    )

    # Disturbance variables
    m.F0 = Var(
        m.time, initialize=0.1, units=pyunits.m**3 / pyunits.min, doc="Inlet flowrate"
    )
    m.T0 = Var(m.time, initialize=350, units=pyunits.K, doc="Inlet temperature")
    m.c0 = Var(
        m.time,
        initialize=1,
        units=pyunits.kmol / pyunits.m**3,
        doc="Inlet concentration",
    )

    # Fixed parameters
    m.r = Param(initialize=0.219, mutable=True, units=pyunits.m, doc="Tank_radius")
    m.k0 = Param(
        initialize=7.2e10,
        mutable=True,
        units=1 / pyunits.min,
        doc="Rate constant preexponential factor",
    )
    m.E_over_R = Param(
        initialize=8750,
        mutable=True,
        units=pyunits.K,
        doc="Rate constant exponential term",
    )
    m.U = Param(
        initialize=54.94,
        mutable=True,
        units=pyunits.kJ / pyunits.min / pyunits.m**2 / pyunits.K,
        doc="Overall heat transfer coefficient",
    )
    m.rho = Param(
        initialize=1000,
        mutable=True,
        units=pyunits.kg / pyunits.m**3,
        doc="Liquid Density",
    )
    m.Cp = Param(
        initialize=0.239,
        mutable=True,
        units=pyunits.kJ / pyunits.kg / pyunits.K,
        doc="Liquid heat capacity",
    )
    m.dH_rxn = Param(
        initialize=-5e4,
        mutable=True,
        units=pyunits.kJ / pyunits.kmol,
        doc="Enthalpy of reaction",
    )

    # Intermediate quantities
    m.A = Var(
        initialize=Constants.pi * 0.219**2, units=pyunits.m**2, doc="Tank base area"
    )
    m.rxn_rate = Var(
        m.time, initialize=1, units=pyunits.kmol / pyunits.m**3 / pyunits.min
    )

    # Equations
    @m.Constraint()
    def A_eqn(b):
        return b.A == Constants.pi * m.r**2

    @m.Constraint(m.time)
    def rxn_rate_eqn(b, t):
        return b.rxn_rate[t] == b.k0 * exp(-m.E_over_R / b.T[t]) * b.c[t]

    @m.Constraint(m.time)
    def c_eqn(b, t):
        return b.cdot[t] == (
            # F0 multiplying both c0 and c is the correct
            # ODE, it is not a typo. It is the result of
            # converting a material balance equation into
            # an equation for concentration
            b.F0[t] * (b.c0[t] - b.c[t]) / (b.A * b.h[t])
            - b.rxn_rate[t]
        )

    @m.Constraint(m.time)
    def T_eqn(b, t):
        return b.Tdot[t] == (
            # Same deal for temperature as with concentration
            b.F0[t] * (b.T0[t] - b.T[t]) / (b.A * b.h[t])
            - b.dH_rxn / (b.rho * b.Cp) * b.rxn_rate[t]
            + 2 * b.U / b.r / b.rho / b.Cp * (b.Tc[t] - b.T[t])
        )

    @m.Constraint(m.time)
    def h_eqn(b, t):
        return b.hdot[t] == (b.F0[t] - b.F[t]) / b.A

    fd_factory = TransformationFactory("dae.finite_difference")
    fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

    m.c0.fix()
    m.T0.fix()

    m.F0.fix()
    m.Tc.fix()
    m.F.fix()

    m.c[0].fix()
    m.T[0].fix()
    m.h[0].fix()

    for t in m.time:
        m.c[t].set_value(0.877773845919673)
        m.T[t].set_value(324.50842862814466)
        m.h[t].set_value(0.6586682605839007)
        calculate_variable_from_constraint(m.rxn_rate[t], m.rxn_rate_eqn[t])

    cont_sys = linearize_system(
        m,
        m.time,
        representative_time=1,
        input_vars=[m.Tc, m.F],
        disturbance_vars=[m.F0],
        output_vars=[m.c, m.T, m.h],
    )

    assert cont_sys["C"] == pytest.approx(eye(3))
    assert cont_sys["D"] == pytest.approx(zeros((3, 2)))
    assert cont_sys["Dd"] == pytest.approx(zeros((3, 1)))

    disc_sys = c2d(cont_sys, dt=1)
    assert disc_sys["A"] == pytest.approx(
        array([[0.2681, -0.00338, -0.00728], [9.703, 0.3279, -25.44], [0, 0, 1]]),
        rel=1e-2,
    )
    assert disc_sys["B"] == pytest.approx(
        array([[-0.00537, 0.1655], [1.297, 97.91], [0, -6.637]]), rel=1e-2
    )
    assert disc_sys["Bd"] == pytest.approx(
        array([[-0.1175], [69.74], [6.637]]), rel=1e-2
    )


@pytest.mark.component
def test_pendulum_example():
    """
    Testing an example of a pendulum system taken from a set of slides
    posted online by Jovana Andrejevic and Catherine Ding.
    https://people.math.wisc.edu/~chr/am205/g_act/DAE_slides.pdf

    This test covers cases where discretization equations are deactivated
    as part of an index reduction scheme, as well as creating LTI models
    of individual blocks (instead of the entire model).
    """
    m = ConcreteModel()
    m.time = ContinuousSet(initialize=[0, 1])
    m.f1 = Block()

    init = {(0, "x"): 0, (0, "y"): -1, (1, "x"): 0, (1, "y"): -1}
    m.f1.s = Var(
        m.time,
        ["x", "y"],
        initialize=init,
        units=pyunits.m,
        doc="displacement of pendulum",
    )
    m.f1.v = Var(
        m.time,
        ["x", "y"],
        initialize=0,
        units=pyunits.m / pyunits.s,
        doc="velocity of pendulum",
    )
    # The example called T "tension", but it does not have the units of tension
    m.f1.T = Var(
        m.time,
        initialize=1,
        units=1 / pyunits.s**2,
        doc="centripedal tension of pendulum",
    )

    m.f1.sdot = DerivativeVar(m.f1.s, wrt=m.time, initialize=0)
    m.f1.vdot = DerivativeVar(m.f1.v, wrt=m.time, initialize=0)

    @m.f1.Constraint(m.time, ["x", "y"])
    def sdot_eqn(b, t, c):
        return b.sdot[t, c] == b.v[t, c]

    @m.f1.Constraint(m.time, ["x", "y"])
    def vdot_eqn(b, t, c):
        if c == "y":
            term = -Constants.acceleration_gravity
        else:
            term = 0 * pyunits.m / pyunits.s**2
        return b.vdot[t, c] == -b.T[t] * b.s[t, c] + term

    @m.f1.Constraint(m.time)
    def pendulum_rod_eqn(b, t):
        return b.s[t, "x"] ** 2 + b.s[t, "y"] ** 2 == 1

    # There is also a two variable encoding using polar coordinates
    m.f2 = Block()
    m.f2.theta = Var(m.time, initialize=0)
    m.f2.thetadot = DerivativeVar(m.f2.theta, wrt=m.time, initialize=0)
    m.f2.omega = Var(m.time, initialize=0)
    m.f2.omegadot = DerivativeVar(m.f2.omega, wrt=m.time, initialize=0)

    @m.f2.Constraint(m.time)
    def thetadot_eqn(b, t):
        return b.thetadot[t] == b.omega[t]

    @m.f2.Constraint(m.time)
    def omegadot_eqn(b, t):
        return b.omegadot[t] == -Constants.acceleration_gravity * sin(b.theta[t])

    fd_factory = TransformationFactory("dae.finite_difference")
    fd_factory.apply_to(m, wrt=m.time, nfe=1, scheme="BACKWARD")

    # We want to reduce the index of the first formulation by getting rid of the y
    # equations and replacing them with derivatives of the rod equation

    m.f1.sdot_disc_eq[:, "y"].deactivate()
    m.f1.vdot_disc_eq[:, "y"].deactivate()

    @m.f1.Constraint(m.time)
    def first_deriv_eqn(b, t):
        return 0 == b.s[t, "x"] * b.v[t, "x"] + b.s[t, "y"] * b.v[t, "y"]

    @m.f1.Constraint(m.time)
    def second_deriv_eqn(b, t):
        return (
            0
            == b.v[t, "x"] ** 2
            + b.v[t, "y"] ** 2
            - b.T[t]
            - Constants.acceleration_gravity * b.s[t, "y"]
        )

    # Now calculate the tension
    for t in m.time:
        calculate_variable_from_constraint(m.f1.T[t], m.f1.second_deriv_eqn[t])

    # Both of these encodings should produce the same linearized approximation
    # for small angles because x = L * sin(theta), and the pendulum length L
    # equals 1 here.
    for blk in [m.f1, m.f2]:
        out = linearize_system(
            blk,
            m.time,
            representative_time=1,
            input_vars=[],
            disturbance_vars=[],
            output_vars=[],
        )

        assert out["A"] == pytest.approx(
            array([[0, 1], [-value(Constants.acceleration_gravity), 0]])
        )
        assert out["B"].shape == (2, 0)
        assert out["Bd"].shape == (2, 0)
        assert out["C"].shape == (0, 2)
        assert out["D"].shape == (0, 0)
        assert out["Dd"].shape == (0, 0)
