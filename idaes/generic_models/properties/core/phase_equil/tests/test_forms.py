##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Library of common forms for phase equilibrium constraints
"""
from pyomo.environ import ConcreteModel, Var, units as pyunits

from idaes.generic_models.properties.core.generic.generic_property import \
    GenericParameterBlock
from idaes.generic_models.properties.core.state_definitions import FTPx

from idaes.generic_models.properties.core.phase_equil.forms import *
import pytest


# Dummy EoS to use for fugacity calls
class DummyEoS(object):
    def common(self):
        pass

    def build_parameters(b):
        pass

    def fug_phase_comp_eq(b, p, j, pp):
        return b.x[p, j]

    def log_fug_phase_comp_eq(b, p, j, pp):
        return 42*b.x[p, j]


@pytest.mark.unit
def test_fugacity():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    m.params = GenericParameterBlock(default={
        "components": {"H2O": {"parameter_data": {"temperature_crit": 647.3},
                               "phase_equilibrium_form": {
                                   ("Vap", "Liq"): fugacity}}},
        "phases": {"Liq": {"equation_of_state": DummyEoS},
                   "Vap": {"equation_of_state": DummyEoS}},
        "state_definition": FTPx,
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K}})

    assert str(fugacity(m, "Vap", "Liq", "H2O")) == str(
        m.x["Vap", "H2O"] == m.x["Liq", "H2O"])


@pytest.mark.unit
def test_log_fugacity():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    m.params = GenericParameterBlock(default={
        "components": {"H2O": {"parameter_data": {"temperature_crit": 647.3},
                               "phase_equilibrium_form": {
                                   ("Vap", "Liq"): log_fugacity}}},
        "phases": {"Liq": {"equation_of_state": DummyEoS},
                   "Vap": {"equation_of_state": DummyEoS}},
        "state_definition": FTPx,
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K}})

    assert str(log_fugacity(m, "Vap", "Liq", "H2O")) == str(
        42*m.x["Vap", "H2O"] == 42*m.x["Liq", "H2O"])
