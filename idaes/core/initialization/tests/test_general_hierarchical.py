#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for general hierarchical initialization routines
"""
import pytest
import types

from pyomo.environ import Block, ConcreteModel, Constraint, value, Var

from idaes.core.initialization.general_hierarchical import (
    SingleControlVolumeUnitInitializer,
)

from idaes.core.util.exceptions import InitializationError


__author__ = "Andrew Lee"


# TODO: Unit testing of methods
@pytest.mark.unit
def test_nonstandard_model():
    m = ConcreteModel()
    m.initialization_order = [m]

    initializer = SingleControlVolumeUnitInitializer()

    with pytest.raises(
        TypeError,
        match="Model unknown does not appear to be a standard form unit model. "
        "Please use an Initializer specific to the model being initialized.",
    ):
        initializer.initialize(m)


class DummyInit:
    def plugin_prepare(self, plugin, output="default"):
        plugin._test = output

    def plugin_initialize(self, model, **kwargs):
        model._plugin_initialized = True

    def plugin_finalize(self, plugin, output="default"):
        plugin._test2 = output

    def initialize(self, model, **kwargs):
        model._initialized = True


@pytest.mark.unit
def test_init_props_0D_no_copy():
    m = ConcreteModel()
    m.control_volume = Block()
    m.control_volume.properties_in = Block()
    m.control_volume.properties_in.v1 = Var()

    m.control_volume.properties_out = Block()
    m.control_volume.properties_out.v1 = Var()
    m.control_volume.properties_out.v2 = Var()

    def estimate_outlet_state(block, *args, **kwargs):
        block._estimated = True

    m.control_volume.estimate_outlet_state = types.MethodType(
        estimate_outlet_state, m.control_volume
    )

    initializer = SingleControlVolumeUnitInitializer()

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer._init_props_0D(m.control_volume, False)

    assert m.control_volume.properties_in._initialized
    assert m.control_volume.properties_out._initialized
    assert m.control_volume._estimated


@pytest.mark.unit
def test_init_props_0D_copy():
    m = ConcreteModel()
    m.control_volume = Block()
    m.control_volume.properties_in = Block()
    m.control_volume.properties_in.v1 = Var(initialize=12)

    m.control_volume.properties_out = Block()
    m.control_volume.properties_out.v1 = Var(initialize=6)
    m.control_volume.properties_out.v2 = Var(initialize=7)

    initializer = SingleControlVolumeUnitInitializer()

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer._init_props_0D(m.control_volume, True)

    assert m.control_volume.properties_in._initialized
    assert not hasattr(m.control_volume.properties_out, "_initialized")
    assert value(m.control_volume.properties_out.v1) == 12
    assert value(m.control_volume.properties_out.v2) == 7


@pytest.mark.unit
def test_init_props_1D():
    m = ConcreteModel()
    m.control_volume = Block()
    m.control_volume.properties = Block()

    initializer = SingleControlVolumeUnitInitializer()

    def estimate_states(block, *args, **kwargs):
        block._estimated = True

    m.control_volume.estimate_states = types.MethodType(
        estimate_states, m.control_volume
    )

    initializer.add_submodel_initializer(m.control_volume.properties, DummyInit())
    initializer._init_props_1D(m.control_volume)

    assert m.control_volume.properties._initialized
    assert m.control_volume._estimated


@pytest.mark.unit
def test_init_rxns():
    m = ConcreteModel()
    m.control_volume = Block()
    m.control_volume.reactions = Block()

    # Add a constraint that to the reactions that used a var from outside the block
    m.control_volume.var = Var()
    m.control_volume.reactions.cons = Constraint(expr=m.control_volume.var == 10)

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())
    initializer._init_rxns(m.control_volume)

    assert m.control_volume.reactions._initialized
    assert not m.control_volume.var.fixed


@pytest.mark.unit
def test_initialize_control_volume_no_copy():
    m = ConcreteModel()
    m.control_volume = Block()
    m.control_volume.properties_in = Block()
    m.control_volume.properties_in.v1 = Var(initialize=12)

    m.control_volume.properties_out = Block()
    m.control_volume.properties_out.v1 = Var(initialize=6)
    m.control_volume.properties_out.v2 = Var(initialize=7)

    m.control_volume.reactions = Block()
    m.control_volume.reactions.v1 = Var()

    m.control_volume.c1 = Constraint(
        expr=m.control_volume.properties_in.v1 == m.control_volume.properties_out.v1
    )
    m.control_volume.c2 = Constraint(
        expr=m.control_volume.reactions.v1 == 0.5 * m.control_volume.properties_out.v1
    )

    m.control_volume.properties_in.v1.fix()

    def estimate_outlet_state(block, *args, **kwargs):
        block._estimated = True

    m.control_volume.estimate_outlet_state = types.MethodType(
        estimate_outlet_state, m.control_volume
    )

    initializer = SingleControlVolumeUnitInitializer()

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())

    initializer.initialize_control_volume(m.control_volume, False)

    assert value(m.control_volume.properties_in.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v1) == pytest.approx(6, rel=1e-5)
    assert value(m.control_volume.properties_out.v2) == pytest.approx(7, rel=1e-5)

    assert m.control_volume.properties_in._initialized
    assert m.control_volume.properties_out._initialized
    assert m.control_volume.reactions._initialized
    assert m.control_volume._estimated


@pytest.mark.unit
def test_initialize_control_volume_copy():
    m = ConcreteModel()
    m.control_volume = Block()
    m.control_volume.properties_in = Block()
    m.control_volume.properties_in.v1 = Var(initialize=12)

    m.control_volume.properties_out = Block()
    m.control_volume.properties_out.v1 = Var(initialize=6)
    m.control_volume.properties_out.v2 = Var(initialize=7)

    m.control_volume.reactions = Block()
    m.control_volume.reactions.v1 = Var()

    m.control_volume.c1 = Constraint(
        expr=m.control_volume.properties_in.v1 == m.control_volume.properties_out.v1
    )
    m.control_volume.c2 = Constraint(
        expr=m.control_volume.reactions.v1 == 0.5 * m.control_volume.properties_out.v1
    )

    m.control_volume.properties_in.v1.fix()

    initializer = SingleControlVolumeUnitInitializer()

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())

    initializer.initialize_control_volume(m.control_volume, True)

    assert value(m.control_volume.properties_in.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v2) == pytest.approx(7, rel=1e-5)

    assert m.control_volume.properties_in._initialized
    assert not hasattr(m.control_volume.properties_out, "_initialized")
    assert m.control_volume.reactions._initialized


@pytest.mark.unit
def test_initialize_submodels_no_order():
    m = ConcreteModel()
    m.initialization_order = []

    initializer = SingleControlVolumeUnitInitializer()

    with pytest.raises(
        InitializationError,
        match="Main model \(unknown\) was not initialized \(no results returned\). "
        "This is likely due to an error in the model.initialization_order.",
    ):
        initializer.initialize_submodels(m, {}, {}, copy_inlet_state=False)


@pytest.mark.unit
def test_initialize_submodels_no_plugins_no_copy():
    m = ConcreteModel()
    m.initialization_order = [m]

    m.control_volume = Block()
    m.control_volume.properties_in = Block()
    m.control_volume.properties_in.v1 = Var(initialize=12)

    m.control_volume.properties_out = Block()
    m.control_volume.properties_out.v1 = Var(initialize=6)
    m.control_volume.properties_out.v2 = Var(initialize=7)

    m.control_volume.reactions = Block()
    m.control_volume.reactions.v1 = Var()

    # plug-in block, but not in initialization order so should be skipped
    m.plugin = Block()

    m.control_volume.c1 = Constraint(
        expr=m.control_volume.properties_in.v1 == m.control_volume.properties_out.v1
    )
    m.control_volume.c2 = Constraint(
        expr=m.control_volume.reactions.v1 == 0.5 * m.control_volume.properties_out.v1
    )

    m.control_volume.properties_in.v1.fix()

    def estimate_outlet_state(block, *args, **kwargs):
        block._estimated = True

    m.control_volume.estimate_outlet_state = types.MethodType(
        estimate_outlet_state, m.control_volume
    )

    initializer = SingleControlVolumeUnitInitializer()

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())
    initializer.add_submodel_initializer(m.plugin, DummyInit())

    initializer.initialize_submodels(
        m, {m.plugin: {}}, {m.plugin: DummyInit()}, copy_inlet_state=False
    )

    assert value(m.control_volume.properties_in.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v2) == pytest.approx(7, rel=1e-5)
    assert value(m.control_volume.reactions.v1) == pytest.approx(6, rel=1e-5)

    assert m.control_volume.properties_in._initialized
    assert m.control_volume.properties_out._initialized
    assert m.control_volume.reactions._initialized
    assert m.control_volume._estimated
    assert not hasattr(m.plugin, "_initialized")
    assert not hasattr(m.plugin, "_plugin_initialized")


@pytest.mark.unit
def test_initialize_submodels_no_plugins_copy():
    m = ConcreteModel()
    m.initialization_order = [m]

    m.control_volume = Block()
    m.control_volume.properties_in = Block()
    m.control_volume.properties_in.v1 = Var(initialize=12)

    m.control_volume.properties_out = Block()
    m.control_volume.properties_out.v1 = Var(initialize=6)
    m.control_volume.properties_out.v2 = Var(initialize=7)

    m.control_volume.reactions = Block()
    m.control_volume.reactions.v1 = Var()

    # plug-in block, but not in initialization order so should be skipped
    m.plugin = Block()

    m.control_volume.c1 = Constraint(
        expr=m.control_volume.properties_in.v1 == m.control_volume.properties_out.v1
    )
    m.control_volume.c2 = Constraint(
        expr=m.control_volume.reactions.v1 == 0.5 * m.control_volume.properties_out.v1
    )

    m.control_volume.properties_in.v1.fix()

    initializer = SingleControlVolumeUnitInitializer()

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())
    initializer.add_submodel_initializer(m.plugin, DummyInit())

    initializer.initialize_submodels(
        m, {m.plugin: {}}, {m.plugin: DummyInit()}, copy_inlet_state=True
    )

    assert value(m.control_volume.properties_in.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v2) == pytest.approx(7, rel=1e-5)
    assert value(m.control_volume.reactions.v1) == pytest.approx(6, rel=1e-5)

    assert m.control_volume.properties_in._initialized
    assert not hasattr(m.control_volume.properties_out, "_initialized")
    assert m.control_volume.reactions._initialized
    assert not hasattr(m.plugin, "_initialized")
    assert not hasattr(m.plugin, "_plugin_initialized")


@pytest.mark.unit
def test_initialize_submodels_w_plugins_copy():
    m = ConcreteModel()
    m.initialization_order = [m]

    m.control_volume = Block()
    m.control_volume.properties_in = Block()
    m.control_volume.properties_in.v1 = Var(initialize=12)

    m.control_volume.properties_out = Block()
    m.control_volume.properties_out.v1 = Var(initialize=6)
    m.control_volume.properties_out.v2 = Var(initialize=7)

    m.control_volume.reactions = Block()
    m.control_volume.reactions.v1 = Var()

    m.plugin = Block()
    m.initialization_order.append(m.plugin)

    m.control_volume.c1 = Constraint(
        expr=m.control_volume.properties_in.v1 == m.control_volume.properties_out.v1
    )
    m.control_volume.c2 = Constraint(
        expr=m.control_volume.reactions.v1 == 0.5 * m.control_volume.properties_out.v1
    )

    m.control_volume.properties_in.v1.fix()

    initializer = SingleControlVolumeUnitInitializer()

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())

    initializer.initialize_submodels(
        m, {m.plugin: {}}, {m.plugin: DummyInit()}, copy_inlet_state=True
    )

    assert value(m.control_volume.properties_in.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v2) == pytest.approx(7, rel=1e-5)
    assert value(m.control_volume.reactions.v1) == pytest.approx(6, rel=1e-5)

    assert m.control_volume.properties_in._initialized
    assert not hasattr(m.control_volume.properties_out, "_initialized")
    assert m.control_volume.reactions._initialized
    assert not hasattr(m.plugin, "_initialized")
    assert m.plugin._plugin_initialized
