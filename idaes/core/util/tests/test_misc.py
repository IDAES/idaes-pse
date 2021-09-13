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
This module contains miscalaneous utility functions for use in IDAES models.
"""

import pytest

from pyomo.environ import ConcreteModel, Expression, Set, Block, Var, units
from pyomo.network import Port, Arc
from pyomo.common.config import ConfigBlock
from pyomo.core.base.units_container import UnitsError

from idaes.core.util.misc import (add_object_reference, copy_port_values,
                                  TagReference, VarLikeExpression,
                                  set_param_from_config)
import idaes.logger as idaeslog


# Author: Andrew Lee
@pytest.mark.unit
def test_add_object_reference():
    m = ConcreteModel()

    m.s = Set(initialize=[1, 2, 3])

    add_object_reference(m, "test_ref", m.s)

    assert hasattr(m, "test_ref")
    assert m.test_ref == m.s


# Author: Andrew Lee
@pytest.mark.unit
def test_add_object_reference_fail():
    m = ConcreteModel()

    with pytest.raises(AttributeError):
        add_object_reference(m, "test_ref", m.s)


# Author: John Eslick
@pytest.mark.unit
def test_port_copy():
    m = ConcreteModel()
    m.b1 = Block()
    m.b2 = Block()
    m.b1.x = Var(initialize=3)
    m.b1.y = Var([0, 1], initialize={0: 4, 1: 5})
    m.b1.z = Var([0, 1], ["A", "B"], initialize={
        (0, "A"): 6, (0, "B"): 7, (1, "A"): 8, (1, "B"): 9})
    m.b2.x = Var(initialize=1)
    m.b2.y = Var([0, 1], initialize=1)
    m.b2.z = Var([0, 1], ["A", "B"], initialize=1)
    m.b1.port = Port()
    m.b2.port = Port()
    m.b1.port.add(m.b1.x, "x")
    m.b1.port.add(m.b1.y, "y")
    m.b1.port.add(m.b1.z, "z")
    m.b2.port.add(m.b2.x, "x")
    m.b2.port.add(m.b2.y, "y")
    m.b2.port.add(m.b2.z, "z")

    def assert_copied_right():
        assert(m.b2.x.value == 3)
        assert(m.b2.y[0].value == 4)
        assert(m.b2.y[1].value == 5)
        assert(m.b2.z[0, "A"].value == 6)
        assert(m.b2.z[0, "B"].value == 7)
        assert(m.b2.z[1, "A"].value == 8)
        assert(m.b2.z[1, "B"].value == 9)

    def reset():
        m.b2.x = 0
        m.b2.y[0] = 0
        m.b2.y[1] = 0
        m.b2.z[0, "A"] = 0
        m.b2.z[0, "B"] = 0
        m.b2.z[1, "A"] = 0
        m.b2.z[1, "B"] = 0
    m.arc = Arc(source=m.b1.port, dest=m.b2.port)

    copy_port_values(m.b2.port, m.b1.port)
    assert_copied_right()
    reset()

    copy_port_values(source=m.b1.port, destination=m.b2.port)
    assert_copied_right()
    reset()

    copy_port_values(m.arc)
    assert_copied_right()
    reset()

    copy_port_values(arc=m.arc)
    assert_copied_right()
    reset()

    with pytest.raises(AttributeError):
        copy_port_values(arc=m.b1.port)

    with pytest.raises(RuntimeError):
        copy_port_values(source=m.b1.port, destination=m.b2.port, arc=m.arc)

    with pytest.raises(RuntimeError):
        copy_port_values(source=m.b1.port, arc=m.arc)

    with pytest.raises(AttributeError):
        copy_port_values(source=m.b1.port, destination=m.arc)


# Author: John Eslick
@pytest.mark.unit
def test_tag_reference():
    m = ConcreteModel()
    m.z = Var([0, 1], ["A", "B"], initialize={
        (0, "A"): 6, (0, "B"): 7, (1, "A"): 8, (1, "B"): 9})
    test_tag = {}
    test_tag["MyTag34&@!e.5"] = TagReference(m.z[:, "A"], description="z tag")
    assert(len(test_tag["MyTag34&@!e.5"]) == 2)
    assert(test_tag["MyTag34&@!e.5"][0].value == 6)
    assert(test_tag["MyTag34&@!e.5"][1].value == 8)
    assert(test_tag["MyTag34&@!e.5"].description == "z tag")
    m.b = Block([0, 1])
    m.b[0].y = Var(initialize=1)
    m.b[1].y = Var(initialize=2)
    test_tag = TagReference(m.b[:].y, description="y tag")
    assert(test_tag[0].value == 1)
    assert(test_tag[1].value == 2)
    assert(test_tag.description == "y tag")


@pytest.mark.unit
def test_SimpleVarLikeExpression():
    m = ConcreteModel()

    # Need a Var to use in the Expression to avoid being able to set the value
    # of a float
    m.v = Var(initialize=42)

    m.e = VarLikeExpression(expr=m.v)

    assert m.e.type() is Expression
    assert not m.e.is_indexed()
    assert m.e.value == 42

    with pytest.raises(TypeError,
                       match="e is an Expression and does not have a value "
                       "which can be set."):
        m.e.value = 10

    with pytest.raises(TypeError,
                       match="e is an Expression and can not have bounds. "
                       "Use an inequality Constraint instead."):
        m.e.setub(10)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not have bounds. "
                       "Use an inequality Constraint instead."):
        m.e.setlb(0)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not be fixed. "
                       "Use an equality Constraint instead."):
        m.e.fix(8)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not be unfixed."):
        m.e.unfix()


@pytest.mark.unit
def test_IndexedVarLikeExpression():
    m = ConcreteModel()

    # Need a Var to use in the Expression to avoid being able to set the value
    # of a float
    m.v = Var(initialize=42)

    m.e = VarLikeExpression([1, 2, 3, 4], expr=m.v)

    assert m.e.type() is Expression
    assert m.e.is_indexed()

    with pytest.raises(TypeError,
                       match="e is an Expression and can not have bounds. "
                       "Use inequality Constraints instead."):
        m.e.setub(10)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not have bounds. "
                       "Use inequality Constraints instead."):
        m.e.setlb(0)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not be fixed. "
                       "Use equality Constraints instead."):
        m.e.fix(8)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not be unfixed."):
        m.e.unfix()

    for i in m.e:
        assert m.e[i].value == 42

        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and does not have "
                           "a value which can be set.".format(i)):
            m.e[i].value = 10

        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and can not have "
                           "bounds. Use an inequality Constraint instead."
                           .format(i)):
            m.e[i].setub(10)
        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and can not have "
                           "bounds. Use an inequality Constraint instead."
                           .format(i)):
            m.e[i].setlb(0)
        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and can not be "
                           "fixed. Use an equality Constraint instead."
                           .format(i)):
            m.e[i].fix(8)
        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and can not be "
                           "unfixed.".format(i)):
            m.e[i].unfix()


@pytest.mark.unit
class TestSetParamFromConfig():
    def test_default_config(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": 42}

        m.b.test_param = Var(initialize=1)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

        assert ("b no units provided for parameter test_param - assuming "
                "default units" in caplog.text)

    def test_specified_config(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config2 = ConfigBlock(implicit=True)
        m.b.config2.parameter_data = {"test_param": 42}

        m.b.test_param = Var(initialize=1)

        set_param_from_config(
            m.b, "test_param", config=m.b.config2)

        assert m.b.test_param.value == 42

    def test_no_config(self):
        m = ConcreteModel()
        m.b = Block()

        m.b.test_param = Var(initialize=1)

        with pytest.raises(AttributeError,
                           match="b - set_param_from_config method was "
                           "not provided with a config argument, but no "
                           "default Config block exists. Please specify the "
                           "Config block to use via the config argument."):
            set_param_from_config(m.b, "test_param")

    def test_invalid_config(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = "foo"

        m.b.test_param = Var(initialize=1)

        with pytest.raises(TypeError,
                           match="b - set_param_from_config - config argument "
                           "provided is not an instance of a Config Block."):
            set_param_from_config(m.b, "test_param")

    def test_no_param(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": 42}

        with pytest.raises(AttributeError,
                           match="b - set_param_from_config method was "
                           "provided with param argument test_param, but no "
                           "attribute of that name exists."):
            set_param_from_config(m.b, "test_param")

    def test_no_parameter_data(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {}

        m.b.test_param = Var(initialize=1)

        with pytest.raises(KeyError,
                           match="b - set_param_from_config method was "
                           "provided with param argument test_param, but the "
                           "config block does not contain a value for this "
                           "parameter."):
            set_param_from_config(m.b, "test_param")

    def test_indexed(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": {"1": 42}}

        m.b.test_param_1 = Var(initialize=1)

        set_param_from_config(m.b, "test_param", index="1")

        assert m.b.test_param_1.value == 42

    def test_no_param_indexed(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": {"1": 42}}

        m.b.test_param = Var(initialize=1)

        with pytest.raises(AttributeError,
                           match="b - set_param_from_config method was "
                           "provided with param and index arguments "
                           "test_param 1, but no attribute with that "
                           "combination \(test_param_1\) exists."):
            set_param_from_config(m.b, "test_param", index="1")

    def test_no_parameter_data_indexed(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": {"2": 42}}

        m.b.test_param_1 = Var(initialize=1)

        with pytest.raises(KeyError,
                           match="b - set_param_from_config method was "
                           "provided with param and index arguments "
                           "test_param 1, but the config block does not "
                           "contain a value for this parameter and index."):
            set_param_from_config(m.b, "test_param", index="1")

    def test_dimensionless_default(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": 42}

        m.b.test_param = Var(initialize=1, units=units.dimensionless)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

        assert ("b no units provided for parameter test_param - assuming "
                "default units" in caplog.text)

    def test_dimensionless_defined(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.dimensionless)}

        m.b.test_param = Var(initialize=1, units=units.dimensionless)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_dimensionless_defined_none(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, None)}

        m.b.test_param = Var(initialize=1, units=units.dimensionless)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_none_defined(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, None)}

        m.b.test_param = Var(initialize=1, units=None)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_none_defined_dimensionless(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.dimensionless)}

        m.b.test_param = Var(initialize=1, units=None)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_consistent_units(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.m)}

        m.b.test_param = Var(initialize=1, units=units.m)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_inconsistent_units(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.m)}

        m.b.test_param = Var(initialize=1, units=units.s)

        with pytest.raises(UnitsError,
                           match="Cannot convert m to s. Units are not "
                           "compatible."):
            set_param_from_config(m.b, "test_param")

    def test_inconsistent_units_dimensionless(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.dimensionless)}

        m.b.test_param = Var(initialize=1, units=units.s)

        with pytest.raises(UnitsError,
                           match="Cannot convert dimensionless to s. Units "
                           "are not compatible."):
            set_param_from_config(m.b, "test_param")

    def test_inconsistent_units_none(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, None)}

        m.b.test_param = Var(initialize=1, units=units.s)

        with pytest.raises(UnitsError,
                           match="Cannot convert None to s. Units "
                           "are not compatible."):
            set_param_from_config(m.b, "test_param")

    def test_inconsistent_units_dimensionless_2(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.s)}

        m.b.test_param = Var(initialize=1, units=units.dimensionless)

        with pytest.raises(UnitsError,
                           match="Cannot convert s to None. Units "
                           "are not compatible."):
            set_param_from_config(m.b, "test_param")

    def test_inconsistent_units_none_2(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.s)}

        m.b.test_param = Var(initialize=1, units=None)

        with pytest.raises(UnitsError,
                           match="Cannot convert s to None. Units "
                           "are not compatible."):
            set_param_from_config(m.b, "test_param")

    def test_unitted_default(self, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": 42}

        m.b.test_param = Var(initialize=1, units=units.m)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

        assert ("b no units provided for parameter test_param - assuming "
                "default units" in caplog.text)
