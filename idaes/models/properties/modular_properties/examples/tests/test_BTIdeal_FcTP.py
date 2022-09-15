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
Author: Andrew Lee
"""
import pytest
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Expression,
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

from idaes.core import LiquidPhase, VaporPhase

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.properties.modular_properties.state_definitions import FcTP
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity

import idaes.models.properties.modular_properties.pure.Perrys as Perrys
import idaes.models.properties.modular_properties.pure.RPP4 as RPP4
import idaes.models.properties.modular_properties.pure.NIST as NIST

from idaes.models.properties.tests.test_harness import PropertyTestHarness


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def _as_quantity(x):
    unit = pyunits.get_units(x)
    if unit is None:
        unit = pyunits.dimensionless
    return value(x) * unit._get_pint_unit()


config_dict = {
    "components": {
        "benzene": {
            "type": Component,
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": 78.1136e-3,  # [1]
                "pressure_crit": 48.9e5,  # [1]
                "temperature_crit": 562.2,  # [1]
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": 1.0162,  # [2] pg. 2-98
                    "2": 0.2655,
                    "3": 562.16,
                    "4": 0.28212,
                },
                "cp_mol_ig_comp_coeff": {
                    "A": -3.392e1,  # [1]
                    "B": 4.739e-1,
                    "C": -3.017e-4,
                    "D": 7.130e-8,
                },
                "cp_mol_liq_comp_coeff": {
                    "1": 1.29e2,  # [2]
                    "2": -1.7e-1,
                    "3": 6.48e-4,
                    "4": 0,
                    "5": 0,
                },
                "enth_mol_form_liq_comp_ref": 49.0e3,  # [3]
                "enth_mol_form_vap_comp_ref": 82.9e3,  # [3]
                "pressure_sat_comp_coeff": {
                    "A": 4.72583,  # [NIST]
                    "B": 1660.652,
                    "C": -1.461,
                },
            },
        },
        "toluene": {
            "type": Component,
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": 92.1405e-3,  # [1]
                "pressure_crit": 41e5,  # [1]
                "temperature_crit": 591.8,  # [1]
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": 0.8488,  # [2] pg. 2-98
                    "2": 0.26655,
                    "3": 591.8,
                    "4": 0.2878,
                },
                "cp_mol_ig_comp_coeff": {
                    "A": -2.435e1,
                    "B": 5.125e-1,
                    "C": -2.765e-4,
                    "D": 4.911e-8,
                },
                "cp_mol_liq_comp_coeff": {
                    "1": 1.40e2,  # [2]
                    "2": -1.52e-1,
                    "3": 6.95e-4,
                    "4": 0,
                    "5": 0,
                },
                "enth_mol_form_liq_comp_ref": 12.0e3,  # [3]
                "enth_mol_form_vap_comp_ref": 50.1e3,  # [3]
                "pressure_sat_comp_coeff": {
                    "A": 4.07827,  # [NIST]
                    "B": 1343.943,
                    "C": -53.773,
                },
            },
        },
    },
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
    },
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "state_definition": FcTP,
    "state_bounds": {
        "flow_mol_comp": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 450, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,
}


class TestBTIdeal_FcTP(PropertyTestHarness):
    def configure(self):
        self.prop_pack = GenericParameterBlock
        self.param_args = config_dict
        self.prop_args = {}
        self.has_density_terms = False


class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**config_dict)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ["benzene", "toluene"]
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 4
        for i in model.params._phase_component_set:
            assert i in [
                ("Liq", "benzene"),
                ("Liq", "toluene"),
                ("Vap", "benzene"),
                ("Vap", "toluene"),
            ]

        assert model.params.config.state_definition == FcTP

        assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            {
                "flow_mol_comp": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 450, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 2
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"benzene": ("Vap", "Liq")},
            "PE2": {"toluene": ("Vap", "Liq")},
        }

        assert model.params.pressure_ref.value == 1e5
        assert model.params.temperature_ref.value == 300

        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**config_dict)

        model.props = model.params.state_block_class(
            [1], parameters=model.params, defined_state=True
        )

        model.props[1].calculate_scaling_factors()

        # Fix state
        model.props[1].flow_mol_comp.fix(0.5)
        model.props[1].temperature.fix(368)
        model.props[1].pressure.fix(101325)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Expression)
        assert isinstance(model.props[1].flow_mol_comp, Var)
        for j in model.params.component_list:
            assert value(model.props[1].flow_mol_comp[j]) == 0.5
            assert model.props[1].flow_mol_comp[j].ub == 1000
            assert model.props[1].flow_mol_comp[j].lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 101325
        assert model.props[1].pressure.ub == 1e6
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 368
        assert model.props[1].temperature.ub == 450
        assert model.props[1].temperature.lb == 273.15

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 2
        for i in model.props[1].mole_frac_comp:
            assert value(model.props[1].mole_frac_comp[i]) == 0.5

        assert_units_consistent(model)

    @pytest.mark.unit
    def test_basic_scaling(self, model):

        assert len(model.props[1].scaling_factor) == 25
        assert model.props[1].scaling_factor[model.props[1].flow_mol] == 1e-2
        assert (
            model.props[1].scaling_factor[model.props[1].flow_mol_phase["Liq"]] == 1e-2
        )
        assert (
            model.props[1].scaling_factor[model.props[1].flow_mol_phase["Vap"]] == 1e-2
        )
        assert (
            model.props[1].scaling_factor[model.props[1].flow_mol_comp["benzene"]]
            == 1e-2
        )
        assert (
            model.props[1].scaling_factor[model.props[1].flow_mol_comp["toluene"]]
            == 1e-2
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].flow_mol_phase_comp["Liq", "benzene"]
            ]
            == 1e-2
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].flow_mol_phase_comp["Liq", "toluene"]
            ]
            == 1e-2
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].flow_mol_phase_comp["Vap", "benzene"]
            ]
            == 1e-2
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].flow_mol_phase_comp["Vap", "toluene"]
            ]
            == 1e-2
        )
        assert (
            model.props[1].scaling_factor[model.props[1].mole_frac_comp["benzene"]]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[model.props[1].mole_frac_comp["toluene"]]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].mole_frac_phase_comp["Liq", "benzene"]
            ]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].mole_frac_phase_comp["Liq", "toluene"]
            ]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].mole_frac_phase_comp["Vap", "benzene"]
            ]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].mole_frac_phase_comp["Vap", "toluene"]
            ]
            == 1000
        )
        assert model.props[1].scaling_factor[model.props[1].pressure] == 1e-5
        assert model.props[1].scaling_factor[model.props[1].temperature] == 1e-2
        assert model.props[1].scaling_factor[model.props[1]._teq["Vap", "Liq"]] == 1e-2
        assert model.props[1].scaling_factor[model.props[1]._t1_Vap_Liq] == 1e-2

        assert (
            model.props[1].scaling_factor[
                model.props[1]._mole_frac_tbub["Vap", "Liq", "benzene"]
            ]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1]._mole_frac_tbub["Vap", "Liq", "toluene"]
            ]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1]._mole_frac_tdew["Vap", "Liq", "benzene"]
            ]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1]._mole_frac_tdew["Vap", "Liq", "toluene"]
            ]
            == 1000
        )
        assert (
            model.props[1].scaling_factor[
                model.props[1].temperature_bubble["Vap", "Liq"]
            ]
            == 1e-2
        )
        assert (
            model.props[1].scaling_factor[model.props[1].temperature_dew["Vap", "Liq"]]
            == 1e-2
        )

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 3
        for i in sv:
            assert i in ["flow_mol_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 3
        for i in sv:
            assert i in ["Molar Flowrate", "Temperature", "Pressure"]

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.props[1]) == 0

    @pytest.mark.component
    def test_initialize(self, model):
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
            "Liq", "benzene"
        ].value == pytest.approx(0.4070, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp[
            "Vap", "benzene"
        ].value == pytest.approx(0.6296, abs=1e-4)
        assert model.props[1].phase_frac["Vap"].value == pytest.approx(0.4177, abs=1e-4)

        assert value(model.props[1].enth_mol) == pytest.approx(48220.5, 1e-5)

    @pytest.mark.unit
    def test_report(self, model):
        model.props[1].report()
