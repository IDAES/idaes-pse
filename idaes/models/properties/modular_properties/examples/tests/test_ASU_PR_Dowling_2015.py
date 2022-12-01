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
"""
Author: Andrew Lee, Alejandro Garciadiego
"""
import pytest
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.unittest import assertStructuredAlmostEqual

from idaes.core import Component
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE

from idaes.models.properties.modular_properties.examples.ASU_PR import (
    configuration_Dowling_2015,
)

from idaes.models.properties.tests.test_harness import PropertyTestHarness


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def _as_quantity(x):
    unit = pyunits.get_units(x)
    if unit is None:
        unit = pyunits.dimensionless
    return value(x) * unit._get_pint_unit()


class TestASUPR(PropertyTestHarness):
    def configure(self):
        self.prop_pack = GenericParameterBlock
        self.param_args = configuration_Dowling_2015
        self.prop_args = {}
        self.has_density_terms = False


# Test for configuration dictionaries with parameters from Properties of Gases
# and liquids 3rd edition and Dowling 2015
# Stream values extracted from:
# A framework for efficient large scale equation-oriented flowsheet optimization (2015)
#     Computers and Chemical Engineering - Alexander W. Dowling
class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration_Dowling_2015)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 3
        for i in model.params.component_list:
            assert i in ["nitrogen", "argon", "oxygen"]
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 6
        for i in model.params._phase_component_set:
            assert i in [
                ("Liq", "nitrogen"),
                ("Liq", "argon"),
                ("Liq", "oxygen"),
                ("Vap", "nitrogen"),
                ("Vap", "argon"),
                ("Vap", "oxygen"),
            ]

        assert model.params.config.state_definition == FTPx

        assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (10, 300, 350, pyunits.K),
                "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 3
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"nitrogen": ("Vap", "Liq")},
            "PE2": {"argon": ("Vap", "Liq")},
            "PE3": {"oxygen": ("Vap", "Liq")},
        }

        assert model.params.pressure_ref.value == 101325
        assert model.params.temperature_ref.value == 298.15

        assert model.params.nitrogen.mw.value == 28.0135e-3
        assert model.params.nitrogen.pressure_crit.value == 33.943875e5
        assert model.params.nitrogen.temperature_crit.value == 126.2

        assert model.params.argon.mw.value == 39.948e-3
        assert model.params.argon.pressure_crit.value == 48.737325e5
        assert model.params.argon.temperature_crit.value == 150.86

        assert model.params.oxygen.mw.value == 31.999e-3
        assert model.params.oxygen.pressure_crit.value == 50.45985e5
        assert model.params.oxygen.temperature_crit.value == 154.58

        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration_Dowling_2015)

        model.props = model.params.build_state_block([1], defined_state=True)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 100
        assert model.props[1].flow_mol.ub == 1000
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 1e5
        assert model.props[1].pressure.ub == 1e7
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 300
        assert model.props[1].temperature.ub == 350
        assert model.props[1].temperature.lb == 10

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 3
        for i in model.props[1].mole_frac_comp:
            assert value(model.props[1].mole_frac_comp[i]) == 1 / 3

        assert_units_consistent(model)

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "Total Molar Flowrate",
                "Total Mole Fraction",
                "Temperature",
                "Pressure",
            ]

    @pytest.mark.component
    def test_initialize(self, model):
        # Fix state
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(77.5295)
        model.props[1].pressure.fix(101000)
        model.props[1].mole_frac_comp["nitrogen"].fix(0.9823)
        model.props[1].mole_frac_comp["argon"].fix(0.00803)
        model.props[1].mole_frac_comp["oxygen"].fix(0.00967)

        assert degrees_of_freedom(model.props[1]) == 0

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={"tol": 1e-6})

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.props[1].mole_frac_phase_comp[
            "Liq", "nitrogen"
        ].value == pytest.approx(0.9663, abs=1e-3)
        assert model.props[1].mole_frac_phase_comp[
            "Vap", "nitrogen"
        ].value == pytest.approx(0.9890, abs=1e-3)
        assert model.props[1].phase_frac["Vap"].value == pytest.approx(0.703, abs=1e-3)
        assert value(model.props[1].enth_mol_phase["Vap"]) == pytest.approx(
            -6567, abs=1e1
        )

    @pytest.mark.unit
    def test_report(self, model):
        model.props[1].report()

    @pytest.mark.component
    def test_SF0(self, model):
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(300.00)
        model.props[1].pressure.fix(353140)
        model.props[1].mole_frac_comp["nitrogen"].fix(0.7800)
        model.props[1].mole_frac_comp["argon"].fix(0.0100)
        model.props[1].mole_frac_comp["oxygen"].fix(0.2100)

        assert degrees_of_freedom(model.props[1]) == 0

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={"tol": 1e-6})

        assert degrees_of_freedom(model) == 0

        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

        assert model.props[1].mole_frac_phase_comp[
            "Vap", "nitrogen"
        ].value == pytest.approx(0.7800, abs=1e-3)
        assert model.props[1].phase_frac["Vap"].value == pytest.approx(1.0, abs=1e-3)
        assert value(model.props[1].enth_mol_phase["Vap"]) == pytest.approx(
            25.7, abs=1e1
        )

    @pytest.mark.component
    def test_SFIL3(self, model):
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(92.88)
        model.props[1].pressure.fix(353140)
        model.props[1].mole_frac_comp["nitrogen"].fix(0.6653)
        model.props[1].mole_frac_comp["argon"].fix(0.0140)
        model.props[1].mole_frac_comp["oxygen"].fix(0.3207)

        assert degrees_of_freedom(model.props[1]) == 0

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={"tol": 1e-6})

        assert degrees_of_freedom(model) == 0

        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

        assert model.props[1].mole_frac_phase_comp[
            "Liq", "nitrogen"
        ].value == pytest.approx(0.6653, abs=1e-3)
        assert model.props[1].phase_frac["Vap"].value == pytest.approx(0.0, abs=1e-3)
        assert value(model.props[1].enth_mol_phase["Liq"]) == pytest.approx(
            -11662.4, abs=1e1
        )
