import pytest
from itertools import combinations

from numpy import array, eye, logspace, zeros

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    exp,
    Param,
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

from idaes.core.solvers import get_solver
from idaes.core.solvers.petsc_object import PETScIntegrator
from idaes.core.util.constants import Constants
from idaes.core.util.dynamics.linearization import (
    _validate_vardata_collections,
    _validate_steady_state,
    c2d,
    linearize_system,
)
from idaes.core.util.exceptions import BurntToast


@pytest.mark.unit
class TestValidateVardataCollections(object):
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
