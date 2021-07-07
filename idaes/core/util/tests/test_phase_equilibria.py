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
This module contains phase equilibria utility functions for use in IDAES models.
"""

__author__ = "Alejandro Garciadiego"

import pytest

import idaes.logger as idaeslog

from pyomo.environ import ConcreteModel
# Import Pyomo units
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_equivalent

from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.core.phases import PhaseType as PT

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import \
        IdealBubbleDew, LogBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import \
    fugacity, log_fugacity

import idaes.generic_models.properties.core.pure.NIST as NIST
import idaes.generic_models.properties.core.pure.RPP5 as RPP5

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from idaes.core.util.phase_equilibria import (TXYDataClass, Txy_data)


@pytest.mark.unit
def test_Txy_dataclass():

    TD = TXYDataClass('benzene', 'toluene',
                      pyunits.kg / pyunits.m / pyunits.s ** 2, pyunits.K,
                      101325)
    TD.TBubb = [353.3205, 365.3478, 383.8817]
    TD.TDew = [353.3237, 372.0203, 383.8845]
    TD.x = [0.9999, 0.5, 0.01]

    assert TD.Component_1 == 'benzene'
    assert TD.Component_2 == 'toluene'
    assert_units_equivalent(TD.Punits, pyunits.kg / pyunits.m / pyunits.s ** 2)
    assert_units_equivalent(TD.Tunits, pyunits.K)
    assert TD.P == 101325
    assert TD.TBubb == [pytest.approx(353.3205, abs=1e-4),
                        pytest.approx(365.3478, abs=1e-4),
                        pytest.approx(383.8817, abs=1e-4)]
    assert TD.TDew == [pytest.approx(353.3237, abs=1e-4),
                       pytest.approx(372.0203, abs=1e-4),
                       pytest.approx(383.8845, abs=1e-4)]
    assert TD.x == [pytest.approx(0.9999, abs=1e-4),
                    pytest.approx(0.5, abs=1e-4),
                    pytest.approx(0.01, abs=1e-4)]


# Author: Alejandro Garciadiego
@pytest.mark.component
def test_Txy_data():
    configuration = {
        # Specifying components
        "components": {
            'benzene': {"type": Component,
                        "pressure_sat_comp": NIST,
                        "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                        "parameter_data": {
                            "mw": (78.1136E-3, pyunits.kg/pyunits.mol),  # [1]
                            "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                            "temperature_crit": (562.2, pyunits.K),  # [1]
                            "pressure_sat_comp_coeff":{"A": (4.72583, None),  # [2]
                                                        "B": (1660.652, pyunits.K),
                                                        "C": (-1.461, pyunits.K)}}},
            'toluene': {"type": Component,
                        "pressure_sat_comp": NIST,
                        "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                        "parameter_data": {
                            "mw": (92.1405E-3, pyunits.kg/pyunits.mol),  # [1]
                            "pressure_crit": (41e5, pyunits.Pa),  # [1]
                            "temperature_crit": (591.8, pyunits.K),  # [1]
                            "pressure_sat_comp_coeff":{"A": (4.07827, None),  # [2]
                                                        "B": (1343.943, pyunits.K),
                                                        "C": (-53.773, pyunits.K)}}}},

        # Specifying phases
        "phases":  {'Liq': {"type": LiquidPhase,
                            "equation_of_state": Ideal},
                    'Vap': {"type": VaporPhase,
                            "equation_of_state": Ideal}},

        # Set base units of measurement
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},

        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                         "temperature": (273.15, 300, 450, pyunits.K),
                         "pressure": (5e4, 1e5, 1e6, pyunits.Pa)},
        "pressure_ref": (1e5, pyunits.Pa),
        "temperature_ref": (300, pyunits.K),

        # Defining phase equilibria
        "phases_in_equilibrium": [("Vap", "Liq")],
        "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
        "bubble_dew_method": IdealBubbleDew}

    model = ConcreteModel()

    model.params = GenericParameterBlock(default=configuration)

    TD = Txy_data(model, 'benzene', 'toluene', 101325,
                  num_points=3, temperature=298.15,
                  print_level=idaeslog.CRITICAL, solver_op={'tol': 1e-6})

    assert TD.Component_1 == 'benzene'
    assert TD.Component_2 == 'toluene'
    assert_units_equivalent(TD.Punits, pyunits.kg / pyunits.m / pyunits.s ** 2)
    assert_units_equivalent(TD.Tunits, pyunits.K)
    assert TD.P == 101325
    assert TD.TBubb == [pytest.approx(353.4853, abs=1e-4),
                        pytest.approx(365.2127, abs=1e-4),
                        pytest.approx(383.2909, abs=1e-4)]
    assert TD.TDew == [pytest.approx(353.7978, abs=1e-2),
                       pytest.approx(371.8702, abs=1e-4),
                       pytest.approx(383.5685, abs=1e-4)]
    assert TD.x == [pytest.approx(0.99, abs=1e-4),
                    pytest.approx(0.5, abs=1e-4),
                    pytest.approx(0.01, abs=1e-4)]


# Author: Alejandro Garciadiego
@pytest.mark.component
def test_Txy_data_no_dew():
    configuration = {
        # Specifying components
    "components": {
        "bmimPF6": {"type": Component,
                    "valid_phase_types": PT.liquidPhase,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (284.18E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (24e5, pyunits.Pa),  # [1]
                        "temperature_crit": (860, pyunits.K),  # [1]
                        "omega": 0.7917}},

        "carbon_dioxide": {"type": Component,
                  "pressure_sat_comp": NIST,
                  "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                  "parameter_data": {
                      "mw": (44.010E-3, pyunits.kg/pyunits.mol),
                      "pressure_crit": (71.8e5, pyunits.Pa),
                      "temperature_crit": (304.1, pyunits.K),
                      "omega": 0.239,
                      "pressure_sat_comp_coeff": {"A": (6.81228, None),
                                                  "B": (1301.679, pyunits.K),
                                                  "C": (-3.494, pyunits.K)}}}},

        # Specifying phases
        "phases":  {'Liq': {"type": LiquidPhase,
                            "equation_of_state": Cubic,
                            "equation_of_state_options": {
                                "type": CubicType.PR}},
                    'Vap': {"type": VaporPhase,
                            "equation_of_state": Cubic,
                            "equation_of_state_options": {
                                "type": CubicType.PR}}},

        # Set base units of measurement
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},

        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                         "temperature": (10, 300, 500, pyunits.K),
                         "pressure": (5e-4, 1e5, 1e10, pyunits.Pa)},
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),

        # Defining phase equilibria
        "phases_in_equilibrium": [("Vap", "Liq")],
        "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
        "bubble_dew_method": LogBubbleDew,
        "parameter_data": {"PR_kappa": {("bmimPF6", "bmimPF6"): 0.000,
                                        ("bmimPF6", "carbon_dioxide"): -0.05718,
                                        ("carbon_dioxide", "carbon_dioxide"): 0.000,
                                        ("carbon_dioxide", "bmimPF6"): 0.0293}}}

    model = ConcreteModel()

    model.params = GenericParameterBlock(default=configuration)

    TD = Txy_data(model, "carbon_dioxide", "bmimPF6", 101325,
                  num_points=3, temperature=150.15,
                  print_level=idaeslog.CRITICAL, solver_op={'tol': 1e-6})

    assert TD.Component_1 == 'carbon_dioxide'
    assert TD.Component_2 == 'bmimPF6'
    assert_units_equivalent(TD.Punits, pyunits.kg / pyunits.m / pyunits.s ** 2)
    assert_units_equivalent(TD.Tunits, pyunits.K)
    assert TD.P == 101325
    assert TD.TBubb == [pytest.approx(185.6468, abs=1e-4),
                        pytest.approx(191.4697, abs=1e-4),
                        pytest.approx(331.1625, abs=1e-4)]
    assert TD.TDew == []
    assert TD.x == [pytest.approx(0.99, abs=1e-4),
                    pytest.approx(0.5, abs=1e-4),
                    pytest.approx(0.01, abs=1e-4)]


# Author: Alejandro Garciadiego
@pytest.mark.component
def test_Txy_data_no_liq():
    configuration = {
        # Specifying components
        "components": {
        "methane": {"type": Component,
                    "valid_phase_types": PT.vaporPhase,
                    "parameter_data": {
                        "mw": (16.043E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (46e5, pyunits.Pa),  # [1]
                        "temperature_crit": (190.4, pyunits.K),  # [1]
                        "omega": 0.011}},

        "ethane": {"type": Component,
                    "pressure_sat_comp": RPP5,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (30.070E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (48.8e5, pyunits.Pa),  # [1]
                        "temperature_crit": (305.4, pyunits.K),  # [1]
                        "omega": 0.099,
                        "pressure_sat_comp_coeff": {"A": (3.95405, None),  # [4]
                                                    "B": (663.720, pyunits.K),
                                                    "C": (256.681, pyunits.K)}}}},

    # Specifying phases
    "phases":  {"Liq": {"type": LiquidPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.PR}},
                "Vap": {"type": VaporPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.PR}}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                     "temperature": (10, 300, 1500, pyunits.K),
                     "pressure": (5e-4, 1e5, 1e10, pyunits.Pa)},
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),

    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {"PR_kappa": {("methane", "methane"): 0.000,
                                    ("methane", "ethane"): 0.000,
                                    ("ethane", "ethane"): 0.000,
                                    ("ethane", "methane"): 0.000}}}

    model = ConcreteModel()

    model.params = GenericParameterBlock(default=configuration)

    TD = Txy_data(model, 'methane', 'ethane', 101325,
                  num_points=3, temperature=298,
                  print_level=idaeslog.CRITICAL, solver_op={'tol': 1e-6})

    assert TD.Component_1 == 'methane'
    assert TD.Component_2 == 'ethane'
    assert_units_equivalent(TD.Punits, pyunits.kg / pyunits.m / pyunits.s ** 2)
    assert_units_equivalent(TD.Tunits, pyunits.K)
    assert TD.P == 101325
    assert TD.TBubb == []
    assert TD.TDew == [pytest.approx(126.9025, abs=1e-2),
                       pytest.approx(172.1489, abs=1e-4),
                       pytest.approx(184.2534, abs=1e-4)]
    assert TD.x == [pytest.approx(0.99, abs=1e-4),
                    pytest.approx(0.5, abs=1e-4),
                    pytest.approx(0.01, abs=1e-4)]
