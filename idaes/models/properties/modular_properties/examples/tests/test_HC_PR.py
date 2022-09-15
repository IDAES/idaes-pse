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
    assert_optimal_termination,
    ConcreteModel,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
import pyomo.common.unittest as unittest

from idaes.core import Component
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver
from idaes.core.util.performance import PerformanceBaseClass

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE

from idaes.models.properties.modular_properties.examples.HC_PR import configuration


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def _as_quantity(x):
    unit = pyunits.get_units(x)
    if unit is None:
        unit = pyunits.dimensionless
    return value(x) * unit._get_pint_unit()


# Test for configuration dictionaries with parameters from Properties of Gases
# and liquids 4th edition
def build_model():
    model = ConcreteModel()
    model.params = GenericParameterBlock(**configuration)

    model.props = model.params.build_state_block([1], defined_state=True)

    # Fix state
    model.props[1].flow_mol.fix(1)
    model.props[1].temperature.fix(295.00)
    model.props[1].pressure.fix(1e5)
    model.props[1].mole_frac_comp["hydrogen"].fix(0.077)
    model.props[1].mole_frac_comp["methane"].fix(0.077)
    model.props[1].mole_frac_comp["ethane"].fix(0.077)
    model.props[1].mole_frac_comp["propane"].fix(0.077)
    model.props[1].mole_frac_comp["nbutane"].fix(0.077)
    model.props[1].mole_frac_comp["ibutane"].fix(0.077)
    model.props[1].mole_frac_comp["ethylene"].fix(0.077)
    model.props[1].mole_frac_comp["propene"].fix(0.077)
    model.props[1].mole_frac_comp["butene"].fix(0.077)
    model.props[1].mole_frac_comp["pentene"].fix(0.077)
    model.props[1].mole_frac_comp["hexene"].fix(0.077)
    model.props[1].mole_frac_comp["heptene"].fix(0.077)
    model.props[1].mole_frac_comp["octene"].fix(0.076)

    return model


def initialize_model(model):
    model.props.initialize(optarg={"tol": 1e-6})


@pytest.mark.performance
class Test_HC_PR_Performance(PerformanceBaseClass, unittest.TestCase):
    def build_model(self):
        return build_model()

    def initialize_model(self, model):
        initialize_model(model)


class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 13
        for i in model.params.component_list:
            assert i in [
                "hydrogen",
                "methane",
                "ethane",
                "propane",
                "nbutane",
                "ibutane",
                "ethylene",
                "propene",
                "butene",
                "pentene",
                "hexene",
                "heptene",
                "octene",
            ]

            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 24
        for i in model.params._phase_component_set:
            assert i in [
                ("Liq", "ethane"),
                ("Vap", "hydrogen"),
                ("Vap", "methane"),
                ("Vap", "ethane"),
                ("Liq", "propane"),
                ("Liq", "nbutane"),
                ("Liq", "ibutane"),
                ("Vap", "propane"),
                ("Vap", "nbutane"),
                ("Vap", "ibutane"),
                ("Liq", "ethylene"),
                ("Liq", "propene"),
                ("Liq", "butene"),
                ("Vap", "ethylene"),
                ("Vap", "propene"),
                ("Vap", "butene"),
                ("Liq", "pentene"),
                ("Liq", "hexene"),
                ("Liq", "heptene"),
                ("Vap", "pentene"),
                ("Vap", "hexene"),
                ("Vap", "heptene"),
                ("Liq", "octene"),
                ("Vap", "octene"),
            ]

        assert model.params.config.state_definition == FTPx

        unittest.assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 1500, pyunits.K),
                "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 11
        for i in model.params.phase_equilibrium_idx:
            assert i in [
                "PE1",
                "PE2",
                "PE3",
                "PE4",
                "PE5",
                "PE6",
                "PE7",
                "PE8",
                "PE9",
                "PE10",
                "PE11",
            ]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"ethane": ("Vap", "Liq")},
            "PE2": {"propane": ("Vap", "Liq")},
            "PE3": {"nbutane": ("Vap", "Liq")},
            "PE4": {"ibutane": ("Vap", "Liq")},
            "PE5": {"ethylene": ("Vap", "Liq")},
            "PE6": {"propene": ("Vap", "Liq")},
            "PE7": {"butene": ("Vap", "Liq")},
            "PE8": {"pentene": ("Vap", "Liq")},
            "PE9": {"hexene": ("Vap", "Liq")},
            "PE10": {"heptene": ("Vap", "Liq")},
            "PE11": {"octene": ("Vap", "Liq")},
        }

        assert model.params.pressure_ref.value == 101325
        assert model.params.temperature_ref.value == 298.15

        assert model.params.hydrogen.mw.value == 2.016e-3
        assert model.params.hydrogen.pressure_crit.value == 12.9e5
        assert model.params.hydrogen.temperature_crit.value == 33.2

        assert model.params.methane.mw.value == 16.043e-3
        assert model.params.methane.pressure_crit.value == 46e5
        assert model.params.methane.temperature_crit.value == 190.4

        assert model.params.ethane.mw.value == 30.070e-3
        assert model.params.ethane.pressure_crit.value == 48.8e5
        assert model.params.ethane.temperature_crit.value == 305.4

        assert model.params.propane.mw.value == 44.094e-3
        assert model.params.propane.pressure_crit.value == 42.5e5
        assert model.params.propane.temperature_crit.value == 369.8

        assert model.params.nbutane.mw.value == 58.124e-3
        assert model.params.nbutane.pressure_crit.value == 38.0e5
        assert model.params.nbutane.temperature_crit.value == 425.2

        assert model.params.ibutane.mw.value == 58.124e-3
        assert model.params.ibutane.pressure_crit.value == 36.5e5
        assert model.params.ibutane.temperature_crit.value == 408.2

        assert model.params.ethylene.mw.value == 28.054e-3
        assert model.params.ethylene.pressure_crit.value == 50.5e5
        assert model.params.ethylene.temperature_crit.value == 282.4

        assert model.params.propene.mw.value == 42.081e-3
        assert model.params.propene.pressure_crit.value == 46.2e5
        assert model.params.propene.temperature_crit.value == 365.0

        assert model.params.butene.mw.value == 56.104e-3
        assert model.params.butene.pressure_crit.value == 40.2e5
        assert model.params.butene.temperature_crit.value == 419.3

        assert model.params.pentene.mw.value == 70.135e-3
        assert model.params.pentene.pressure_crit.value == 40.5e5
        assert model.params.pentene.temperature_crit.value == 464.7

        assert model.params.hexene.mw.value == 84.162e-3
        assert model.params.hexene.pressure_crit.value == 31.7e5
        assert model.params.hexene.temperature_crit.value == 504.0

        assert model.params.heptene.mw.value == 98.189e-3
        assert model.params.heptene.pressure_crit.value == 25.4e5
        assert model.params.heptene.temperature_crit.value == 537.2

        assert model.params.octene.mw.value == 112.216e-3
        assert model.params.octene.pressure_crit.value == 26.2e5
        assert model.params.octene.temperature_crit.value == 566.6
        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        return build_model()

    @pytest.mark.integration
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 1
        assert model.props[1].flow_mol.ub == 1000
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 1e5
        assert model.props[1].pressure.ub == 1e7
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 295
        assert model.props[1].temperature.ub == 1500
        assert model.props[1].temperature.lb == 273.15

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 13
        for i in model.props[1].mole_frac_comp:
            assert value(model.props[1].mole_frac_comp[i]) == pytest.approx(
                0.077, abs=1e-2
            )

    @pytest.mark.integration
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.integration
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.integration
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.integration
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

    @pytest.mark.integration
    def test_initialize(self, model):
        assert degrees_of_freedom(model.props[1]) == 0

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        initialize_model(model)

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.integration
    def test_solve(self, model):
        results = solver.solve(model)

        assert_optimal_termination(results)

    @pytest.mark.integration
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.props[1].mole_frac_phase_comp[
            "Vap", "hydrogen"
        ].value == pytest.approx(0.09996, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp[
            "Liq", "propene"
        ].value == pytest.approx(0.01056, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp[
            "Vap", "propene"
        ].value == pytest.approx(0.09681, abs=1e-4)
        assert model.props[1].phase_frac["Vap"].value == pytest.approx(
            0.77026, abs=1e-4
        )

    @pytest.mark.integration
    def test_report(self, model):
        model.props[1].report()
