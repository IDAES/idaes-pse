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
General tests for generic proeprties with Henry components present

Author: Andrew Lee
"""
# Import Python libraries
import pytest

# Import Pyomo components
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Set,
    value,
    units as pyunits,
)

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.phase_equil.henry import ConstantH
import idaes.models.properties.modular_properties.pure.Perrys as Perrys
import idaes.models.properties.modular_properties.pure.RPP4 as RPP4
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.henry import HenryType

import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)
solver = get_solver()


class TestNoHenryComps(object):
    configuration = {
        # Specifying components
        "components": {
            "A": {
                "type": Component,
                "enth_mol_liq_comp": Perrys,
                "enth_mol_ig_comp": RPP4,
                "pressure_sat_comp": RPP4,
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                "parameter_data": {
                    "mw": 78.1136e-3,
                    "pressure_crit": 48.9e5,
                    "temperature_crit": 562.2,
                    "cp_mol_ig_comp_coeff": {
                        "A": -3.392e1,
                        "B": 4.739e-1,
                        "C": -3.017e-4,
                        "D": 7.130e-8,
                    },
                    "cp_mol_liq_comp_coeff": {
                        "1": 1.29e2,
                        "2": -1.7e-1,
                        "3": 6.48e-4,
                        "4": 0,
                        "5": 0,
                    },
                    "enth_mol_form_liq_comp_ref": 49.0e3,
                    "enth_mol_form_vap_comp_ref": 82.9e3,
                    "pressure_sat_comp_coeff": {
                        "A": -6.98273,
                        "B": 1.33213,
                        "C": -2.62863,
                        "D": -3.33399,
                    },
                },
            },
            "B": {
                "type": Component,
                "enth_mol_liq_comp": Perrys,
                "enth_mol_ig_comp": RPP4,
                "pressure_sat_comp": RPP4,
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                "parameter_data": {
                    "mw": 92.1405e-3,
                    "pressure_crit": 41e5,
                    "temperature_crit": 591.8,
                    "cp_mol_ig_comp_coeff": {
                        "A": -2.435e1,
                        "B": 5.125e-1,
                        "C": -2.765e-4,
                        "D": 4.911e-8,
                    },
                    "cp_mol_liq_comp_coeff": {
                        "1": 1.40e2,
                        "2": -1.52e-1,
                        "3": 6.95e-4,
                        "4": 0,
                        "5": 0,
                    },
                    "enth_mol_form_liq_comp_ref": 12.0e3,
                    "enth_mol_form_vap_comp_ref": 50.1e3,
                    "pressure_sat_comp_coeff": {
                        "A": -7.28607,
                        "B": 1.38091,
                        "C": -2.83433,
                        "D": -2.79168,
                    },
                },
            },
        },
        # Specifying phases
        "phases": {
            "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
            "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
        },
        # Declare a base units dict to save code later
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 100, 1000),
            "temperature": (273.15, 300, 450),
            "pressure": (5e4, 1e5, 1e6),
        },
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        # Defining phase equilibria
        "phases_in_equilibrium": [("Vap", "Liq")],
        "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
        "bubble_dew_method": IdealBubbleDew,
    }

    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**self.configuration)

        model.props = model.params.build_state_block([1], defined_state=True)

        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(368)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["A"].fix(0.5)
        model.props[1].mole_frac_comp["B"].fix(0.5)

        # Trigger construction of some things we will test later
        model.props[1].pressure_bubble
        model.props[1].pressure_dew

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ["A", "B"]
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 4
        for i in model.params._phase_component_set:
            assert i in [("Liq", "A"), ("Liq", "B"), ("Vap", "A"), ("Vap", "B")]

        assert model.params.config.state_definition == FTPx

        assert model.params.config.state_bounds == {
            "flow_mol": (0, 100, 1000),
            "temperature": (273.15, 300, 450),
            "pressure": (5e4, 1e5, 1e6),
        }

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 2
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"A": ("Vap", "Liq")},
            "PE2": {"B": ("Vap", "Liq")},
        }

        assert model.params.pressure_ref.value == 1e5
        assert model.params.temperature_ref.value == 300

    @pytest.mark.unit
    def test_init_bubble_temperature(self, model):
        model.props._init_Tbub(model.props[1], pyunits.K)

        assert pytest.approx(365.35, abs=0.01) == value(
            model.props[1].temperature_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.7137, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2863, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "B")]
        )

    @pytest.mark.unit
    def test_init_dew_temperature(self, model):
        model.props._init_Tdew(model.props[1], pyunits.K)

        assert pytest.approx(372.02, abs=0.01) == value(
            model.props[1].temperature_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2909, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.7091, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "B")]
        )

    @pytest.mark.unit
    def test_init_bubble_pressure(self, model):
        model.props._init_Pbub(model.props[1], pyunits.K)

        assert pytest.approx(109479, abs=1) == value(
            model.props[1].pressure_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.7119, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2881, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "B")]
        )

    @pytest.mark.unit
    def test_init_dew_pressure(self, model):
        model.props._init_Pdew(model.props[1], pyunits.K)

        assert pytest.approx(89820, abs=1) == value(
            model.props[1].pressure_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2881, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.7119, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "B")]
        )

    @pytest.mark.component
    def test_solve_vle(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

        assert pytest.approx(365.35, abs=0.01) == value(
            model.props[1].temperature_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.7137, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2863, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "B")]
        )

        assert pytest.approx(372.02, abs=0.01) == value(
            model.props[1].temperature_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2909, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.7091, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "B")]
        )

        assert pytest.approx(109479, abs=1) == value(
            model.props[1].pressure_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.7119, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2881, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "B")]
        )

        assert pytest.approx(89820, abs=1) == value(
            model.props[1].pressure_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2881, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.7119, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "B")]
        )


configuration = {
    # Specifying components
    "components": {
        "A": {
            "type": Component,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": 78.1136e-3,
                "pressure_crit": 48.9e5,
                "temperature_crit": 562.2,
                "cp_mol_ig_comp_coeff": {
                    "A": -3.392e1,
                    "B": 4.739e-1,
                    "C": -3.017e-4,
                    "D": 7.130e-8,
                },
                "cp_mol_liq_comp_coeff": {
                    "1": 1.29e2,
                    "2": -1.7e-1,
                    "3": 6.48e-4,
                    "4": 0,
                    "5": 0,
                },
                "enth_mol_form_liq_comp_ref": 49.0e3,
                "enth_mol_form_vap_comp_ref": 82.9e3,
                "pressure_sat_comp_coeff": {
                    "A": -6.98273,
                    "B": 1.33213,
                    "C": -2.62863,
                    "D": -3.33399,
                },
            },
        },
        "B": {
            "type": Component,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": 92.1405e-3,
                "pressure_crit": 41e5,
                "temperature_crit": 591.8,
                "cp_mol_ig_comp_coeff": {
                    "A": -2.435e1,
                    "B": 5.125e-1,
                    "C": -2.765e-4,
                    "D": 4.911e-8,
                },
                "cp_mol_liq_comp_coeff": {
                    "1": 1.40e2,
                    "2": -1.52e-1,
                    "3": 6.95e-4,
                    "4": 0,
                    "5": 0,
                },
                "enth_mol_form_liq_comp_ref": 12.0e3,
                "enth_mol_form_vap_comp_ref": 50.1e3,
                "pressure_sat_comp_coeff": {
                    "A": -7.28607,
                    "B": 1.38091,
                    "C": -2.83433,
                    "D": -2.79168,
                },
            },
        },
        "C": {
            "type": Component,
            "henry_component": {"Liq": {"method": ConstantH, "type": HenryType.Kpx}},
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": 92.1405e-3,
                "pressure_crit": 41e5,
                "temperature_crit": 591.8,
                "henry_ref": {"Liq": 2e5},
                "cp_mol_ig_comp_coeff": {
                    "A": -2.435e1,
                    "B": 5.125e-1,
                    "C": -2.765e-4,
                    "D": 4.911e-8,
                },
                "cp_mol_liq_comp_coeff": {
                    "1": 1.40e2,
                    "2": -1.52e-1,
                    "3": 6.95e-4,
                    "4": 0,
                    "5": 0,
                },
                "enth_mol_form_liq_comp_ref": 12.0e3,
                "enth_mol_form_vap_comp_ref": 50.1e3,
            },
        },
    },
    # Specifying phases
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
    },
    # Declare a base units dict to save code later
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 100, 1000),
        "temperature": (273.15, 300, 450),
        "pressure": (5e4, 1e5, 1e6),
    },
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,
}


class TestHenryComps0(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration)

        model.props = model.params.build_state_block([1], defined_state=True)

        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(368)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["A"].fix(0.5)
        model.props[1].mole_frac_comp["B"].fix(0.5)
        model.props[1].mole_frac_comp["C"].fix(1e-10)

        # Trigger construction of some things we will test later
        model.props[1].pressure_bubble
        model.props[1].pressure_dew

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 3
        for i in model.params.component_list:
            assert i in ["A", "B", "C"]
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 6
        for i in model.params._phase_component_set:
            assert i in [
                ("Liq", "A"),
                ("Liq", "B"),
                ("Liq", "C"),
                ("Vap", "A"),
                ("Vap", "B"),
                ("Vap", "C"),
            ]

        assert model.params.config.state_definition == FTPx

        assert model.params.config.state_bounds == {
            "flow_mol": (0, 100, 1000),
            "temperature": (273.15, 300, 450),
            "pressure": (5e4, 1e5, 1e6),
        }

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 3
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"A": ("Vap", "Liq")},
            "PE2": {"B": ("Vap", "Liq")},
            "PE3": {"C": ("Vap", "Liq")},
        }

        assert model.params.pressure_ref.value == 1e5
        assert model.params.temperature_ref.value == 300
        assert model.params.C.henry_ref_Liq.value == 2e5

    @pytest.mark.unit
    def test_init_bubble_temperature(self, model):
        model.props._init_Tbub(model.props[1], pyunits.K)

        assert pytest.approx(365.35, abs=0.01) == value(
            model.props[1].temperature_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.7137, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2863, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "B")]
        )

    @pytest.mark.unit
    def test_init_dew_temperature(self, model):
        model.props._init_Tdew(model.props[1], pyunits.K)

        assert pytest.approx(372.02, abs=0.01) == value(
            model.props[1].temperature_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2909, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.7091, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "B")]
        )

    @pytest.mark.unit
    def test_init_bubble_pressure(self, model):
        model.props._init_Pbub(model.props[1], pyunits.K)

        assert pytest.approx(109479, abs=1) == value(
            model.props[1].pressure_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.7119, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2881, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "B")]
        )

    @pytest.mark.unit
    def test_init_dew_pressure(self, model):
        model.props._init_Pdew(model.props[1], pyunits.K)

        assert pytest.approx(89820, abs=1) == value(
            model.props[1].pressure_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2881, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.7119, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "B")]
        )

    @pytest.mark.component
    def test_solve_vle(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

        assert pytest.approx(365.35, abs=0.01) == value(
            model.props[1].temperature_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.7137, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2863, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "B")]
        )

        assert pytest.approx(372.02, abs=0.01) == value(
            model.props[1].temperature_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2909, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.7091, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "B")]
        )

        assert pytest.approx(109479, abs=1) == value(
            model.props[1].pressure_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.7119, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2881, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "B")]
        )

        assert pytest.approx(89820, abs=1) == value(
            model.props[1].pressure_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2881, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.7119, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "B")]
        )


class TestHenryComps(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration)

        model.props = model.params.build_state_block([1], defined_state=True)

        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(368)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["A"].fix(0.45)
        model.props[1].mole_frac_comp["B"].fix(0.45)
        model.props[1].mole_frac_comp["C"].fix(0.1)

        # Trigger construction of some things we will test later
        model.props[1].pressure_bubble
        model.props[1].pressure_dew

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 3
        for i in model.params.component_list:
            assert i in ["A", "B", "C"]
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 6
        for i in model.params._phase_component_set:
            assert i in [
                ("Liq", "A"),
                ("Liq", "B"),
                ("Liq", "C"),
                ("Vap", "A"),
                ("Vap", "B"),
                ("Vap", "C"),
            ]

        assert model.params.config.state_definition == FTPx

        assert model.params.config.state_bounds == {
            "flow_mol": (0, 100, 1000),
            "temperature": (273.15, 300, 450),
            "pressure": (5e4, 1e5, 1e6),
        }

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 3
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"A": ("Vap", "Liq")},
            "PE2": {"B": ("Vap", "Liq")},
            "PE3": {"C": ("Vap", "Liq")},
        }

        assert model.params.pressure_ref.value == 1e5
        assert model.params.temperature_ref.value == 300
        assert model.params.C.henry_ref_Liq.value == 2e5

    @pytest.mark.unit
    def test_init_bubble_temperature(self, model):
        model.props._init_Tbub(model.props[1], pyunits.K)

        assert pytest.approx(361.50, abs=0.01) == value(
            model.props[1].temperature_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.5750, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2276, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "B")]
        )
        assert pytest.approx(0.1974, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "C")]
        )

    @pytest.mark.unit
    def test_init_dew_temperature(self, model):
        model.props._init_Tdew(model.props[1], pyunits.K)

        assert pytest.approx(370.23, abs=0.01) == value(
            model.props[1].temperature_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2750, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.6744, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "B")]
        )
        assert pytest.approx(0.0506, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "C")]
        )

    @pytest.mark.unit
    def test_init_bubble_pressure(self, model):
        model.props._init_Pbub(model.props[1], pyunits.K)

        assert pytest.approx(118531, abs=1) == value(
            model.props[1].pressure_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.5918, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2395, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "B")]
        )
        assert pytest.approx(0.1687, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "C")]
        )

    @pytest.mark.unit
    def test_init_dew_pressure(self, model):
        model.props._init_Pdew(model.props[1], pyunits.K)

        assert pytest.approx(95056, abs=1) == value(
            model.props[1].pressure_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2744, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.6780, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "B")]
        )
        assert pytest.approx(0.0476, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "C")]
        )

    @pytest.mark.component
    def test_solve_vle(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

        assert pytest.approx(361.50, abs=0.01) == value(
            model.props[1].temperature_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.5750, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2276, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "B")]
        )
        assert pytest.approx(0.1974, abs=1e-4) == value(
            model.props[1]._mole_frac_tbub[("Vap", "Liq", "C")]
        )

        assert pytest.approx(370.23, abs=0.01) == value(
            model.props[1].temperature_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2750, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.6744, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "B")]
        )
        assert pytest.approx(0.0506, abs=1e-4) == value(
            model.props[1]._mole_frac_tdew[("Vap", "Liq", "C")]
        )

        assert pytest.approx(118531, abs=1) == value(
            model.props[1].pressure_bubble[("Vap", "Liq")]
        )
        assert pytest.approx(0.5918, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.2395, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "B")]
        )
        assert pytest.approx(0.1687, abs=1e-4) == value(
            model.props[1]._mole_frac_pbub[("Vap", "Liq", "C")]
        )

        assert pytest.approx(95056, abs=1) == value(
            model.props[1].pressure_dew[("Vap", "Liq")]
        )
        assert pytest.approx(0.2744, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "A")]
        )
        assert pytest.approx(0.6780, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "B")]
        )
        assert pytest.approx(0.0476, abs=1e-4) == value(
            model.props[1]._mole_frac_pdew[("Vap", "Liq", "C")]
        )
