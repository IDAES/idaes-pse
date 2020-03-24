##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
                           Block,
                           Expression,
                           Param,
                           Set,
                           value,
                           Var)
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.core.util.misc import add_object_reference


# Dummy EoS to use for fugacity calls
class DummyEoS(object):
    def fug_phase_comp(b, p, j):
        return b.fug_phase_comp[p, j]


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()
    m.params.config = ConfigBlock()
    m.params.config.declare("equation_of_state", ConfigValue(
            default={"Liq": DummyEoS, "Vap": DummyEoS}))

    m.params.component_list = Set(initialize=["H2O"])
    m.params.phase_list = Set(initialize=["Liq", "Vap"])

    m.params.temperature_crit_comp = Var(["H2O"], initialize=647.3)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=300)
    m.props[1].temperature_bubble = Var(initialize=300)
    m.props[1].temperature_dew = Var(initialize=300)

    m.props[1].fug_phase_comp = Var(m.params.phase_list,
                                    m.params.component_list,
                                    initialize=10)

    smooth_VLE.phase_equil(m.props[1])

    return m


def test_build(frame):
    assert isinstance(frame.props[1].eps_1, Param)
    assert value(frame.props[1].eps_1) == 0.01

    assert isinstance(frame.props[1].eps_2, Param)
    assert value(frame.props[1].eps_2) == 0.0005

    assert isinstance(frame.props[1]._t1, Var)
    assert isinstance(frame.props[1]._teq, Var)

    assert isinstance(frame.props[1]._t1_constraint, Constraint)
    assert isinstance(frame.props[1]._teq_constraint, Constraint)
    assert isinstance(frame.props[1].equilibrium_constraint, Constraint)
    for k in frame.props[1].equilibrium_constraint:
        assert k in frame.params.component_list
        assert str(frame.props[1].equilibrium_constraint[k].body) == str(
                frame.props[1].fug_phase_comp["Vap", k] -
                frame.props[1].fug_phase_comp["Liq", k])

    assert isinstance(frame.props[1]._tr_eq, Expression)
    assert len(frame.props[1]._tr_eq) == 1
    assert value(frame.props[1]._tr_eq["H2O"]) == value(
            frame.props[1]._teq/frame.params.temperature_crit_comp["H2O"])


def test_t1(frame):
    # Test that T1 is the max(T, T_bubble)
    # Can't check directly, but see that residual of constraint is correct
    for t in [200, 300, 400, 500]:
        for tb in [200, 300, 400, 500]:
            frame.props[1].temperature.value = t
            frame.props[1].temperature_bubble.value = tb
            frame.props[1]._t1.value = max(t, tb)
            assert value(frame.props[1]._t1_constraint.body) == \
                pytest.approx(0, abs=5e-3)


def test_t_eq(frame):
    # Test that Teq is the min(T1, T_dew)
    # Can't check directly, but see that residual of constraint is correct
    for t1 in [200, 300, 400, 500]:
        for td in [200, 300, 400, 500]:
            frame.props[1]._t1.value = t1
            frame.props[1].temperature_dew.value = td
            frame.props[1]._teq.value = min(t1, td)
            assert value(frame.props[1]._teq_constraint.body) == \
                pytest.approx(0, abs=5e-3)
