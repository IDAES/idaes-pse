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
Author: Brandon Paul
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

from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.properties.modular_properties.state_definitions import FPhx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import \
        IdealBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import \
    fugacity

from idaes.models.properties.examples.methanol_ideal_VLE import config_dict


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
        model.params = GenericParameterBlock(default=config_dict)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 4
        for i in model.params.component_list:
            assert i in [
                "CH4",
                "CO",
                "H2",
                "CH3OH",
            ]

            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 5
        for i in model.params._phase_component_set:
            assert i in [
                ("Vap", "CH4"),
                ("Vap", "CO"),
                ("Vap", "H2"),
                ("Vap", "CH3OH"),
                ("Liq", "CH3OH"),
            ]

        # NOTE: would be nice to reference block instead, the attribute
        # model.params.CH3OH.phase_equilibrium_form does not exist
        assert config_dict["components"]["CH3OH"]["phase_equilibrium_form"] \
            == {
                ("Vap", "Liq"): fugacity
                }

        assert model.params.config.phases == {
            'Liq': {"type": LiquidPhase,
                    "equation_of_state": Ideal},
            'Vap': {"type": VaporPhase,
                    "equation_of_state": Ideal}
            }

        assert model.params.config.state_definition == FPhx

        assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            {
                "flow_mol": (1e-10, 100, 1e10, pyunits.mol/pyunits.s),
                "enth_mol": (-1e10, 100, 1e10, pyunits.J/pyunits.mol),
                "temperature": (198.15, 298.15, 512.75, pyunits.K),
                "pressure": (1e-10, 1e5, 1e10, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert model.params.config.bubble_dew_method == IdealBubbleDew

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 1
        for i in model.params.phase_equilibrium_idx:
            assert i in [
                "PE1",
            ]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"CH3OH": ("Vap", "Liq")},
        }

        assert model.params.pressure_ref.value == 1E5
        assert model.params.temperature_ref.value == 298.15

        assert model.params.CH4.mw.value == 16.043E-3
        assert model.params.CH4.pressure_crit.value == 46.0e5
        assert model.params.CH4.temperature_crit.value == 190.4

        assert model.params.CO.mw.value == 28.010E-3
        assert model.params.CO.pressure_crit.value == 35.0e5
        assert model.params.CO.temperature_crit.value == 132.9

        assert model.params.H2.mw.value == 2.016E-3
        assert model.params.H2.pressure_crit.value == 12.9e5
        assert model.params.H2.temperature_crit.value == 33.0

        assert model.params.CH3OH.mw.value == 32.042E-3
        assert model.params.CH3OH.pressure_crit.value == 80.9e5
        assert model.params.CH3OH.temperature_crit.value == 512.6
        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=config_dict)

        model.props = model.params.build_state_block(
            [1], default={"defined_state": True}
        )

        # Fix state
        model.props[1].flow_mol.fix(415.44)
        model.props[1].enth_mol.fix(-1.3643e5)
        model.props[1].pressure.fix(5.1e6)
        model.props[1].mole_frac_comp["CH4"].fix(2.2964e-6)
        model.props[1].mole_frac_comp["CO"].fix(0.11438)
        model.props[1].mole_frac_comp["H2"].fix(0.23743)
        model.props[1].mole_frac_comp["CH3OH"].fix(0.64818)

        assert degrees_of_freedom(model.props[1]) == 0

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 415.44
        assert model.props[1].flow_mol.ub == 1e10
        assert model.props[1].flow_mol.lb == 1e-10

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 5.1e6
        assert model.props[1].pressure.ub == 1e10
        assert model.props[1].pressure.lb == 1e-10

        assert isinstance(model.props[1].enth_mol, Var)
        assert value(model.props[1].enth_mol) == -1.3643e5
        assert model.props[1].enth_mol.ub == 1e10
        assert model.props[1].enth_mol.lb == -1e10

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 4
        assert model.props[1].mole_frac_comp["CH4"].value == 2.2964e-6
        assert model.props[1].mole_frac_comp["CO"].value == 0.11438
        assert model.props[1].mole_frac_comp["H2"].value == 0.23743
        assert model.props[1].mole_frac_comp["CH3OH"].value == 0.64818

    @pytest.mark.integration
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "enth_mol", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "Total Molar Flowrate",
                "Total Mole Fraction",
                "Molar Enthalpy",
                "Pressure",
            ]

    @pytest.mark.integration
    def test_initialize(self, model):
        # Fix state
        model.props[1].flow_mol.fix(415.44)
        model.props[1].enth_mol.fix(-1.3643e5)
        model.props[1].pressure.fix(5.1e6)
        model.props[1].mole_frac_comp["CH4"].fix(2.2964e-6)
        model.props[1].mole_frac_comp["CO"].fix(0.11438)
        model.props[1].mole_frac_comp["H2"].fix(0.23743)
        model.props[1].mole_frac_comp["CH3OH"].fix(0.64818)

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

    @pytest.mark.integration
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.integration
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.props[1].mole_frac_phase_comp[
            "Vap", "CH3OH"
        ].value == pytest.approx(0.64818, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp[
            "Liq", "CH3OH"
        ].value == pytest.approx(1.00000, abs=1e-4)

    @pytest.mark.unit
    def test_report(self, model):
        model.props[1].report()
