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
Tests for smooth VLE formulation

Authors: Andrew Lee
"""

import pytest

from pyomo.environ import (ConcreteModel,
                           Constraint,
                           Expression,
                           Param,
                           value,
                           Var,
                           units as pyunits)

from idaes.generic_models.properties.core.generic.generic_property import \
    GenericParameterBlock
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.phase_equil import SmoothVLE
from idaes.generic_models.properties.core.phase_equil.forms import fugacity
from idaes.core.util.exceptions import ConfigurationError


# Dummy EoS to use for fugacity calls
class DummyEoS(object):
    def common(self, pobj):
        pass

    def build_parameters(b):
        pass

    def fug_phase_comp(b, p, j, pp):
        return b.fug_phase_comp[p, j]


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = GenericParameterBlock(default={
        "components": {"H2O": {"parameter_data": {"temperature_crit": 647.3},
                               "phase_equilibrium_form":
                                   {("Vap", "Liq"): fugacity}}},
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

    # Create a dummy state block
    m.props = m.params.state_block_class([1], default={"parameters": m.params})

    m.props[1].temperature_bubble = Var([("Liq", "Vap")], initialize=300)
    m.props[1].temperature_dew = Var([("Liq", "Vap")], initialize=300)
    m.props[1]._teq = Var([("Liq", "Vap")], initialize=300)

    m.props[1].fug_phase_comp = Var(m.params.phase_list,
                                    m.params.component_list,
                                    initialize=10)

    SmoothVLE.phase_equil(m.props[1], ("Liq", "Vap"))

    return m


@pytest.mark.unit
def test_build(frame):
    assert isinstance(frame.props[1].eps_1_Liq_Vap, Param)
    assert value(frame.props[1].eps_1_Liq_Vap) == 0.01

    assert isinstance(frame.props[1].eps_2_Liq_Vap, Param)
    assert value(frame.props[1].eps_2_Liq_Vap) == 0.0005

    assert isinstance(frame.props[1]._t1_Liq_Vap, Var)

    assert isinstance(frame.props[1]._t1_constraint_Liq_Vap, Constraint)
    assert isinstance(frame.props[1]._teq_constraint_Liq_Vap, Constraint)


@pytest.mark.unit
def test_t1(frame):
    # Test that T1 is the max(T, T_bubble)
    # Can't check directly, but see that residual of constraint is correct
    for t in [200, 300, 400, 500]:
        for tb in [200, 300, 400, 500]:
            frame.props[1].temperature.value = t
            frame.props[1].temperature_bubble[("Liq", "Vap")].value = tb
            frame.props[1]._t1_Liq_Vap.value = max(t, tb)
            assert value(frame.props[1]._t1_constraint_Liq_Vap.body) == \
                pytest.approx(0, abs=5e-3)


@pytest.mark.unit
def test_t_eq(frame):
    # Test that Teq is the min(T1, T_dew)
    # Can't check directly, but see that residual of constraint is correct
    for t1 in [200, 300, 400, 500]:
        for td in [200, 300, 400, 500]:
            frame.props[1]._t1_Liq_Vap.value = t1
            frame.props[1].temperature_dew[("Liq", "Vap")].value = td
            frame.props[1]._teq[("Liq", "Vap")].value = min(t1, td)
            assert value(frame.props[1]._teq_constraint_Liq_Vap.body) == \
                pytest.approx(0, abs=5e-3)


@pytest.mark.unit
def test_non_VLE_pair():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = GenericParameterBlock(default={
        "components": {"H2O": {"parameter_data": {"temperature_crit": 647.3},
                               "phase_equilibrium_form":
                                   {("Sol", "Liq"): fugacity}}},
        "phases": {"Sol": {"equation_of_state": DummyEoS},
                   "Liq": {"equation_of_state": DummyEoS}},
        "state_definition": FTPx,
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K}})

    # Create a dummy state block
    m.props = m.params.state_block_class([1], default={"parameters": m.params})

    with pytest.raises(ConfigurationError,
                       match="params Generic Property Package phase pair "
                       "Liq-Sol was set to use Smooth VLE formulation, "
                       "however this is not a vapor-liquid pair."):
        SmoothVLE.phase_equil(m.props[1], ("Liq", "Sol"))
