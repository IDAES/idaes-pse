import pytest
from itertools import combinations

from numpy import array, eye, zeros

from pyomo.environ import (
    ConcreteModel,
    exp,
    Param,
    TransformationFactory,
    units as pyunits,
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
    c2d,
    linearize_system,
)


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
        assert output_partition["deriv"] == [base_setup[0].deriv_v]
        assert output_partition["alg"] == [base_setup[0].alg_v]
        assert output_partition["dist"] == [base_setup[0].dist_v]
        assert output_partition["diff"] == []
        assert output_partition["input"] == []

    def test_output_partition_assignment(self, base_setup):
        """Verifies that input and differential variables that show up as outputs
        are correctly assigned to their respective partitioned lists."""
        m, vardata_lists, vardata_sets = base_setup

        # Configure outputs to explicitly include a differential and an input variable
        # alongside an algebraic variable
        vardata_lists["output"] = [m.input_v, m.alg_v, m.diff_v]
        vardata_sets["output"] = ComponentSet(vardata_lists["output"])

        output_partition = _validate_vardata_collections(vardata_lists, vardata_sets)

        # Assert that they are correctly partitioned based on their true classification
        assert output_partition["input"] == [m.input_v]
        assert output_partition["alg"] == [m.alg_v]
        assert output_partition["diff"] == [m.diff_v]

        # Ensure other categories remain empty
        assert output_partition["deriv"] == []
        assert output_partition["dist"] == []

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


@pytest.mark.component
def test_ogunnaike_ray_example_10_4():
    """
    Example 10.4 from Ogunnaike and Ray (1994) for level control
     conical tank. The cone's point faces downwards (like
    an ice cream cone), so the cone's "base" is at its top.
    """
    m = ConcreteModel()
    m.time = ContinuousSet(initialize=[0, 1])

    # Tank geometry
    m.h = Var(m.time, initialize=1, doc="Tank level")
    m.A = Var(m.time, initialize=1, doc="Cross sectional area of the tank at height h")
    m.alpha = Param(initialize=2, doc="Cone shape parameter")
    m.c = Param(initialize=0.5, doc="Cone outflow parameter")
    m.beta = Var(initialize=1, doc="Cone outflow parameter")
    # m.R = Param(initialize=1, doc="Cone base radius")
    # m.H = Param(initialize=1, doc="Cone's overall height")
    m.hdot = DerivativeVar(m.h, wrt=m.time)

    m.Fi = Var(m.time, initialize=1, doc="Inlet flowrate")

    @m.Constraint(m.time)
    def A_eqn(b, t):
        return b.A[t] == b.h[t] ** 2 / b.alpha

    @m.Constraint()
    def beta_eqn(b):
        # Why do we need both beta and c? Ask the authors.
        return b.beta == b.alpha * b.c

    @m.Constraint(m.time)
    def hdot_eqn(b, t):
        b.hdot[t] == b.alpha * b.Fi[t] / b.h[t] ** 2 - b.beta * b.h[t] ** -1.5


@pytest.mark.component
def test_red_book_example():
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

    # Solve for steady state---there is rounding error in the
    # values transcribed from the red book.
    m.cdot[0].fix(0)
    m.Tdot[0].fix(0)
    m.hdot[0].fix(0)

    solver_obj = get_solver("ipopt_v2", solver_options={"constr_viol_tol": 1e-8})
    solver_obj.solve(m)

    for t in m.time:
        m.cdot[t].unfix()
        m.Tdot[t].unfix()
        m.hdot[t].unfix()

    m.c[0].fix()
    m.T[0].fix()
    m.h[0].fix()

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
