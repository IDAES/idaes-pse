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

from idaes.core import Component, FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver


from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.unit_models import Flash
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE

from idaes.models.properties.modular_properties.examples.CO2_bmimPF6_PR import (
    configuration,
)


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
class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()

        model.param = GenericParameterBlock(**configuration)

        assert isinstance(model.param.phase_list, Set)
        assert len(model.param.phase_list) == 2
        for i in model.param.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.param.Liq.is_liquid_phase()
        assert model.param.Vap.is_vapor_phase()

        assert isinstance(model.param.component_list, Set)
        assert len(model.param.component_list) == 2
        for i in model.param.component_list:
            assert i in ["bmimPF6", "carbon_dioxide"]
            assert isinstance(model.param.get_component(i), Component)

        assert isinstance(model.param._phase_component_set, Set)
        assert len(model.param._phase_component_set) == 3
        for i in model.param._phase_component_set:
            assert i in [
                ("Liq", "bmimPF6"),
                ("Liq", "carbon_dioxide"),
                ("Vap", "carbon_dioxide"),
            ]

        assert model.param.config.state_definition == FTPx

        assertStructuredAlmostEqual(
            model.param.config.state_bounds,
            {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (10, 300, 500, pyunits.K),
                "pressure": (5e-4, 1e5, 1e10, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.param.config.phase_equilibrium_state == {("Vap", "Liq"): SmoothVLE}

        assert isinstance(model.param.phase_equilibrium_idx, Set)
        assert len(model.param.phase_equilibrium_idx) == 1
        for i in model.param.phase_equilibrium_idx:
            assert i in ["PE1"]

        assert model.param.phase_equilibrium_list == {
            "PE1": {"carbon_dioxide": ("Vap", "Liq")}
        }

        assert model.param.pressure_ref.value == 101325
        assert model.param.temperature_ref.value == 298.15

        assert model.param.bmimPF6.mw.value == 284.18e-3
        assert model.param.bmimPF6.pressure_crit.value == 24e5
        assert model.param.bmimPF6.temperature_crit.value == 860

        assert model.param.carbon_dioxide.mw.value == 44.010e-3
        assert model.param.carbon_dioxide.pressure_crit.value == 71.8e5
        assert model.param.carbon_dioxide.temperature_crit.value == 304.1

        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()

        model.fs = FlowsheetBlock(dynamic=False)

        model.fs.param = GenericParameterBlock(**configuration)

        model.fs.props = model.fs.param.build_state_block([1], defined_state=True)

        # Fix state
        model.fs.props[1].flow_mol.fix(1)
        model.fs.props[1].temperature.fix(298.15)
        model.fs.props[1].pressure.fix(1214713.75)
        model.fs.props[1].mole_frac_comp["carbon_dioxide"].fix(0.2)
        model.fs.props[1].mole_frac_comp["bmimPF6"].fix(0.8)

        assert degrees_of_freedom(model.fs.props[1]) == 0

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.fs.props[1].flow_mol, Var)
        assert value(model.fs.props[1].flow_mol) == 1
        assert model.fs.props[1].flow_mol.ub == 1000
        assert model.fs.props[1].flow_mol.lb == 0

        assert isinstance(model.fs.props[1].pressure, Var)
        assert value(model.fs.props[1].pressure) == 1214713.75
        assert model.fs.props[1].pressure.ub == 1e10
        assert model.fs.props[1].pressure.lb == 5e-4

        assert isinstance(model.fs.props[1].temperature, Var)
        assert value(model.fs.props[1].temperature) == 298.15
        assert model.fs.props[1].temperature.ub == 500
        assert model.fs.props[1].temperature.lb == 10

        assert isinstance(model.fs.props[1].mole_frac_comp, Var)
        assert len(model.fs.props[1].mole_frac_comp) == 2
        assert value(model.fs.props[1].mole_frac_comp["carbon_dioxide"]) == 0.2
        assert value(model.fs.props[1].mole_frac_comp["bmimPF6"]) == 0.8

        assert_units_consistent(model)

    @pytest.mark.unit
    def test_basic_scaling(self, model):
        model.fs.props[1].calculate_scaling_factors()

        assert len(model.fs.props[1].scaling_factor) == 18
        assert model.fs.props[1].scaling_factor[model.fs.props[1].flow_mol] == 1e-2
        assert (
            model.fs.props[1].scaling_factor[model.fs.props[1].flow_mol_phase["Liq"]]
            == 1e-2
        )
        assert (
            model.fs.props[1].scaling_factor[model.fs.props[1].flow_mol_phase["Vap"]]
            == 1e-2
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].flow_mol_phase_comp["Liq", "bmimPF6"]
            ]
            == 1e-2
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].flow_mol_phase_comp["Liq", "carbon_dioxide"]
            ]
            == 1e-2
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].flow_mol_phase_comp["Vap", "carbon_dioxide"]
            ]
            == 1e-2
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].mole_frac_comp["bmimPF6"]
            ]
            == 1000
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].mole_frac_comp["carbon_dioxide"]
            ]
            == 1000
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].mole_frac_phase_comp["Liq", "bmimPF6"]
            ]
            == 1000
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].mole_frac_phase_comp["Liq", "carbon_dioxide"]
            ]
            == 1000
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].mole_frac_phase_comp["Vap", "carbon_dioxide"]
            ]
            == 1000
        )
        assert model.fs.props[1].scaling_factor[model.fs.props[1].pressure] == 1e-5
        assert model.fs.props[1].scaling_factor[model.fs.props[1].temperature] == 1e-2
        assert (
            model.fs.props[1].scaling_factor[model.fs.props[1]._teq["Vap", "Liq"]]
            == 1e-2
        )
        assert model.fs.props[1].scaling_factor[model.fs.props[1]._t1_Vap_Liq] == 1e-2

        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1]._mole_frac_tbub["Vap", "Liq", "bmimPF6"]
            ]
            == 1
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1]._mole_frac_tbub["Vap", "Liq", "carbon_dioxide"]
            ]
            == 1000
        )
        assert (
            model.fs.props[1].scaling_factor[
                model.fs.props[1].temperature_bubble["Vap", "Liq"]
            ]
            == 1e-2
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.fs.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.fs.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.fs.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "Total Molar Flowrate",
                "Total Mole Fraction",
                "Temperature",
                "Pressure",
            ]

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.fs.props[1].report()


class TestFlashIntegration(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()

        model.fs = FlowsheetBlock(dynamic=False)

        model.fs.param = GenericParameterBlock(**configuration)

        model.fs.unit = Flash(
            property_package=model.fs.param,
            has_heat_transfer=False,
            has_pressure_change=False,
        )
        # Fix state
        model.fs.unit.inlet.flow_mol.fix(1)
        model.fs.unit.inlet.temperature.fix(200.00)
        model.fs.unit.inlet.pressure.fix(101325)
        model.fs.unit.inlet.mole_frac_comp[0, "carbon_dioxide"].fix(1 / 2)
        model.fs.unit.inlet.mole_frac_comp[0, "bmimPF6"].fix(1 / 2)

        assert degrees_of_freedom(model.fs) == 0

        # Apply scaling - model will not solver without this
        model.fs.unit.control_volume.properties_in[0].calculate_scaling_factors()
        model.fs.unit.control_volume.properties_out[0].calculate_scaling_factors()

        return model

    @pytest.mark.component
    def test_initialize(self, model):

        initialization_tester(model, dof=0)

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):
        # Check phase equilibrium results
        assert value(
            model.fs.unit.liq_outlet.mole_frac_comp[0, "carbon_dioxide"]
        ) == pytest.approx(0.3119, abs=1e-4)
        assert value(
            model.fs.unit.vap_outlet.mole_frac_comp[0, "carbon_dioxide"]
        ) == pytest.approx(1.0000, abs=1e-4)
        assert value(
            model.fs.unit.vap_outlet.flow_mol[0] / model.fs.unit.liq_outlet.flow_mol[0]
        ) == pytest.approx(0.37619, abs=1e-4)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.fs.unit.report()
