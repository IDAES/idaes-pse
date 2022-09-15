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
Library of common forms for phase equilibrium constraints
"""
from pyomo.environ import ConcreteModel, Var, units as pyunits

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx

from idaes.models.properties.modular_properties.phase_equil.forms import *
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
        return 42 * b.x[p, j]


@pytest.mark.unit
def test_fugacity():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    m.params = GenericParameterBlock(
        components={
            "H2O": {
                "parameter_data": {"temperature_crit": 647.3},
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            }
        },
        phases={
            "Liq": {"equation_of_state": DummyEoS},
            "Vap": {"equation_of_state": DummyEoS},
        },
        state_definition=FTPx,
        pressure_ref=100000.0,
        temperature_ref=300,
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
    )

    assert str(fugacity.return_expression(m, "Vap", "Liq", "H2O")) == str(
        m.x["Vap", "H2O"] == m.x["Liq", "H2O"]
    )


@pytest.mark.unit
def test_log_fugacity():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    m.params = GenericParameterBlock(
        components={
            "H2O": {
                "parameter_data": {"temperature_crit": 647.3},
                "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            }
        },
        phases={
            "Liq": {"equation_of_state": DummyEoS},
            "Vap": {"equation_of_state": DummyEoS},
        },
        state_definition=FTPx,
        pressure_ref=100000.0,
        temperature_ref=300,
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
    )

    assert str(log_fugacity.return_expression(m, "Vap", "Liq", "H2O")) == str(
        42 * m.x["Vap", "H2O"] == 42 * m.x["Liq", "H2O"]
    )
