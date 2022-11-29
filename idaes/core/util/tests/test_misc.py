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

from pyomo.environ import ConcreteModel, Set, Block, Var, Param, units
from pyomo.network import Port, Arc
from pyomo.common.config import ConfigBlock
from pyomo.core.base.units_container import UnitsError

from idaes.core.util.misc import (
    add_object_reference,
    set_param_from_config,
    is_constant_up_to_units,
)
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


@pytest.mark.unit
class TestSetParamFromConfig:
    def test_default_config(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": 42}

        m.b.test_param = Var(initialize=1)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

        assert (
            "b no units provided for parameter test_param - assuming "
            "default units" in caplog.text
        )

    def test_specified_config(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config2 = ConfigBlock(implicit=True)
        m.b.config2.parameter_data = {"test_param": 42}

        m.b.test_param = Var(initialize=1)

        set_param_from_config(m.b, "test_param", config=m.b.config2)

        assert m.b.test_param.value == 42

    def test_no_config(self):
        m = ConcreteModel()
        m.b = Block()

        m.b.test_param = Var(initialize=1)

        with pytest.raises(
            AttributeError,
            match="b - set_param_from_config method was "
            "not provided with a config argument, but no "
            "default Config block exists. Please specify the "
            "Config block to use via the config argument.",
        ):
            set_param_from_config(m.b, "test_param")

    def test_invalid_config(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = "foo"

        m.b.test_param = Var(initialize=1)

        with pytest.raises(
            TypeError,
            match="b - set_param_from_config - config argument "
            "provided is not an instance of a Config Block.",
        ):
            set_param_from_config(m.b, "test_param")

    def test_no_param(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": 42}

        with pytest.raises(
            AttributeError,
            match="b - set_param_from_config method was "
            "provided with param argument test_param, but no "
            "attribute of that name exists.",
        ):
            set_param_from_config(m.b, "test_param")

    def test_no_parameter_data(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {}

        m.b.test_param = Var(initialize=1)

        with pytest.raises(
            KeyError,
            match="b - set_param_from_config method was "
            "provided with param argument test_param, but the "
            "config block does not contain a value for this "
            "parameter.",
        ):
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

        with pytest.raises(
            AttributeError,
            match="b - set_param_from_config method was "
            "provided with param and index arguments "
            "test_param 1, but no attribute with that "
            "combination \(test_param_1\) exists.",
        ):
            set_param_from_config(m.b, "test_param", index="1")

    def test_no_parameter_data_indexed(self):
        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": {"2": 42}}

        m.b.test_param_1 = Var(initialize=1)

        with pytest.raises(
            KeyError,
            match="b - set_param_from_config method was "
            "provided with param and index arguments "
            "test_param 1, but the config block does not "
            "contain a value for this parameter and index.",
        ):
            set_param_from_config(m.b, "test_param", index="1")

    def test_dimensionless_default(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": 42}

        m.b.test_param = Var(initialize=1, units=units.dimensionless)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

        assert (
            "b no units provided for parameter test_param - assuming "
            "default units" in caplog.text
        )

    def test_dimensionless_defined(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.dimensionless)}

        m.b.test_param = Var(initialize=1, units=units.dimensionless)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_dimensionless_defined_none(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, None)}

        m.b.test_param = Var(initialize=1, units=units.dimensionless)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_none_defined(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, None)}

        m.b.test_param = Var(initialize=1, units=None)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_none_defined_dimensionless(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.dimensionless)}

        m.b.test_param = Var(initialize=1, units=None)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_consistent_units(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.m)}

        m.b.test_param = Var(initialize=1, units=units.m)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

    def test_inconsistent_units(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.m)}

        m.b.test_param = Var(initialize=1, units=units.s)

        with pytest.raises(
            UnitsError, match="Cannot convert m to s. Units are not " "compatible."
        ):
            set_param_from_config(m.b, "test_param")

    def test_inconsistent_units_dimensionless(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.dimensionless)}

        m.b.test_param = Var(initialize=1, units=units.s)

        with pytest.raises(
            UnitsError,
            match="Cannot convert dimensionless to s. Units " "are not compatible.",
        ):
            set_param_from_config(m.b, "test_param")

    def test_inconsistent_units_none(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, None)}

        m.b.test_param = Var(initialize=1, units=units.s)

        with pytest.raises(
            UnitsError, match="Cannot convert None to s. Units " "are not compatible."
        ):
            set_param_from_config(m.b, "test_param")

    def test_inconsistent_units_dimensionless_2(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.s)}

        m.b.test_param = Var(initialize=1, units=units.dimensionless)

        with pytest.raises(
            UnitsError, match="Cannot convert s to None. Units " "are not compatible."
        ):
            set_param_from_config(m.b, "test_param")

    def test_inconsistent_units_none_2(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": (42, units.s)}

        m.b.test_param = Var(initialize=1, units=None)

        with pytest.raises(
            UnitsError, match="Cannot convert s to None. Units " "are not compatible."
        ):
            set_param_from_config(m.b, "test_param")

    def test_unitted_default(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger=("idaes.core.util.misc"))

        m = ConcreteModel()
        m.b = Block()
        m.b.config = ConfigBlock(implicit=True)
        m.b.config.parameter_data = {"test_param": 42}

        m.b.test_param = Var(initialize=1, units=units.m)

        set_param_from_config(m.b, "test_param")

        assert m.b.test_param.value == 42

        assert (
            "b no units provided for parameter test_param - assuming "
            "default units" in caplog.text
        )


@pytest.mark.unit
def test_is_constant_up_to_units():

    # test float
    assert is_constant_up_to_units(42.0)

    # test int
    assert is_constant_up_to_units(42)

    # test unit-only expression
    assert is_constant_up_to_units(42.0 * units.m**2 / units.s)

    m = ConcreteModel()
    m.fixed_param = Param(initialize=42.0)

    # test non-mutable param
    assert is_constant_up_to_units(m.fixed_param)

    # test non-mutable param with unit-only expression
    assert is_constant_up_to_units(m.fixed_param * units.m**2 / units.s)

    m.fixed_param_2 = Param(initialize=6.28)

    # test expression of fixed params with units
    assert is_constant_up_to_units(
        (m.fixed_param**2) * m.fixed_param_2 * units.m**2
    )

    m.mutable_param = Param(initialize=42, mutable=True)

    # test mutable param
    assert not is_constant_up_to_units(m.mutable_param)

    m.mutable_param_units = Param(initialize=42, units=units.m**2 / units.s)
    # test mutable param with units specified
    assert not is_constant_up_to_units(m.mutable_param_units)

    m.fixed_variable = Var(initialize=42)
    m.fixed_variable.fix()
    m.fixed_variable_units = Var(initialize=42, units=units.m**2 / units.s)
    m.fixed_variable_units.fix()

    m.unfixed_variable = Var(initialize=6.28)
    m.unfixed_variable_units = Var(initialize=6.28, units=units.m**2 / units.s)

    # test variables
    assert not is_constant_up_to_units(m.fixed_variable)
    assert not is_constant_up_to_units(m.fixed_variable_units)
    assert not is_constant_up_to_units(m.unfixed_variable)
    assert not is_constant_up_to_units(m.unfixed_variable_units)

    # test combinations
    assert not is_constant_up_to_units(m.fixed_variable * m.fixed_param * units.m**2)
    assert not is_constant_up_to_units(m.fixed_variable * 42.0 * units.m**2)

    assert not is_constant_up_to_units(
        m.unfixed_variable * m.fixed_param * units.m**2
    )
    assert not is_constant_up_to_units(m.unfixed_variable * 42.0 * units.m**2)

    assert not is_constant_up_to_units(
        m.fixed_variable_units * m.fixed_param * units.m**2
    )
    assert not is_constant_up_to_units(m.fixed_variable_units * 42.0 * units.m**2)

    assert not is_constant_up_to_units(
        m.unfixed_variable_units * m.fixed_param * units.m**2
    )
    assert not is_constant_up_to_units(m.unfixed_variable_units * 42.0 * units.m**2)

    assert not is_constant_up_to_units(
        (1 / m.mutable_param_units) * m.fixed_param * units.m**2
    )
    assert not is_constant_up_to_units(
        (1 / m.mutable_param_units) * m.fixed_variable_units * 42.0 * units.m**2
    )
