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
from pyomo.environ import (ConcreteModel,
                           Set,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           units as pyunits)

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from idaes.core.phases import PhaseType as PT

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)
from idaes.core.util import get_solver

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.phase_equil import SmoothVLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import \
        IdealBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import fugacity
from idaes.generic_models.properties.core.phase_equil.henry import ConstantH
import idaes.generic_models.properties.core.pure.Perrys as Perrys
import idaes.generic_models.properties.core.pure.RPP4 as RPP4

import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)
solver = get_solver()


# -----------------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system with N2 and CH4
configuration = {
    # Specifying components
    "components": {
        "A": {"type": Component,
              "enth_mol_liq_comp": Perrys,
              "enth_mol_ig_comp": RPP4,
              "pressure_sat_comp": RPP4,
              "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
              "parameter_data": {
                  "mw": 78.1136E-3,
                  "pressure_crit": 48.9e5,
                  "temperature_crit": 562.2,
                  "cp_mol_ig_comp_coeff": {'A': -3.392E1,
                                           'B': 4.739E-1,
                                           'C': -3.017E-4,
                                           'D': 7.130E-8},
                  "cp_mol_liq_comp_coeff": {'1': 1.29E2,
                                            '2': -1.7E-1,
                                            '3': 6.48E-4,
                                            '4': 0,
                                            '5': 0},
                  "enth_mol_form_liq_comp_ref": 49.0e3,
                  "enth_mol_form_vap_comp_ref": 82.9e3,
                  "pressure_sat_comp_coeff": {'A': -6.98273,
                                              'B': 1.33213,
                                              'C': -2.62863,
                                              'D': -3.33399}}},
        "B": {"type": Component,
              "enth_mol_liq_comp": Perrys,
              "enth_mol_ig_comp": RPP4,
              "pressure_sat_comp": RPP4,
              "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
              "parameter_data": {
                  "mw": 92.1405E-3,
                  "pressure_crit": 41e5,
                  "temperature_crit": 591.8,
                  "cp_mol_ig_comp_coeff": {'A': -2.435E1,
                                           'B': 5.125E-1,
                                           'C': -2.765E-4,
                                           'D': 4.911E-8},
                  "cp_mol_liq_comp_coeff": {'1': 1.40E2,
                                            '2': -1.52E-1,
                                            '3': 6.95E-4,
                                            '4': 0,
                                            '5': 0},
                  "enth_mol_form_liq_comp_ref": 12.0e3,
                  "enth_mol_form_vap_comp_ref": 50.1e3,
                  "pressure_sat_comp_coeff": {'A': -7.28607,
                                              'B': 1.38091,
                                              'C': -2.83433,
                                              'D': -2.79168}}},
        "C": {"type": Component,
              "henry_component": {"Liq": ConstantH},
              "enth_mol_liq_comp": Perrys,
              "enth_mol_ig_comp": RPP4,
              "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
              "parameter_data": {
                  "mw": 92.1405E-3,
                  "pressure_crit": 41e5,
                  "temperature_crit": 591.8,
                  "henry_ref": {"Liq": 2e5},
                  "cp_mol_ig_comp_coeff": {'A': -2.435E1,
                                           'B': 5.125E-1,
                                           'C': -2.765E-4,
                                           'D': 4.911E-8},
                  "cp_mol_liq_comp_coeff": {'1': 1.40E2,
                                            '2': -1.52E-1,
                                            '3': 6.95E-4,
                                            '4': 0,
                                            '5': 0},
                  "enth_mol_form_liq_comp_ref": 12.0e3,
                  "enth_mol_form_vap_comp_ref": 50.1e3}}},

    # Specifying phases
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Ideal},
                'Vap': {"type": VaporPhase,
                        "equation_of_state": Ideal}},

    # Declare a base units dict to save code later
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 100, 1000),
                     "temperature": (273.15, 300, 450),
                     "pressure": (5e4, 1e5, 1e6)},
    "pressure_ref": 1e5,
    "temperature_ref": 300,

    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE}}


class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration)

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
            assert i in [("Liq", "A"), ("Liq", "B"), ("Liq", "C"),
                         ("Vap", "A"), ("Vap", "B"), ("Vap", "C")]

        assert model.params.config.state_definition == FTPx

        assert model.params.config.state_bounds == {
                "flow_mol": (0, 100, 1000),
                "temperature": (273.15, 300, 450),
                "pressure": (5e4, 1e5, 1e6)}

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE}

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 3
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"A": ("Vap", "Liq")},
            "PE2": {"B": ("Vap", "Liq")},
            "PE3": {"C": ("Vap", "Liq")}}

        assert model.params.pressure_ref.value == 1e5
        assert model.params.temperature_ref.value == 300
        assert model.params.C.henry_ref_Liq.value == 2e5


class TestBubbleTemperature(object):
    @pytest.fixture
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration)

        model.props = model.params.build_state_block(
                [1],
                default={"defined_state": True})

        return model

    @pytest.mark.component
    def test_bubble_temperature(self, model):
        # Fix state such that bubble point is known
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(360)
        model.props[1].pressure.fix(76727.7)
        model.props[1].mole_frac_comp["A"].fix(0.5360)
        model.props[1].mole_frac_comp["B"].fix(0.2034)
        model.props[1].mole_frac_comp["C"].fix(0.2607)

        # Trigger build of temperature bubble (if not already present)
        model.props[1].temperature_bubble

        model.props[1]._init_Tbub()

        model.props[1].temperature_bubble.display()
        model.props[1]._mole_frac_tbub.display()

        assert False
