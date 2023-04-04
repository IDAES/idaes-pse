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
def test_prepare_plugins_none():
    m = ConcreteModel()
    m.initialization_order = [m]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    args = {"foo": "bar"}
    subinit, plugin_args = initializer._prepare_plugins(m, args, log)

    assert subinit == {}
    assert plugin_args == {"foo": "bar"}
    # Make sure we didn't change the original args
    assert args == {"foo": "bar"}


@pytest.mark.unit
def test_prepare_plugins_no_args():
    m = ConcreteModel()
    m.b = Block()
    m.initialization_order = [m, m.b]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)
    initializer.add_submodel_initializer(m.b, DummyInit)

    args = {"foo": "bar"}
    subinit, plugin_args = initializer._prepare_plugins(m, args, log)

    assert len(subinit) == 1
    assert isinstance(subinit[m.b], DummyInit)
    assert plugin_args == {"foo": "bar", m.b: {}}
    assert m.b._test == "default"
    # Make sure we didn't change the original args
    assert args == {"foo": "bar"}


@pytest.mark.unit
def test_prepare_plugins_w_args():
    m = ConcreteModel()
    m.b = Block()
    m.initialization_order = [m, m.b]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)
    initializer.add_submodel_initializer(m.b, DummyInit)

    args = {"foo": "bar", m.b: {"output": "checkval"}}
    subinit, plugin_args = initializer._prepare_plugins(m, args, log)

    assert len(subinit) == 1
    assert isinstance(subinit[m.b], DummyInit)
    assert plugin_args == {"foo": "bar", m.b: {"output": "checkval"}}
    assert m.b._test == "checkval"
    # Make sure we didn't change the original args
    assert args == {"foo": "bar", m.b: {"output": "checkval"}}


@pytest.mark.unit
def test_solve_full_model_no_plugins():
    m = ConcreteModel()
    m.initialization_order = [m]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    results = initializer._solve_full_model(m, log, "test_results")

    assert results == "test_results"


@pytest.mark.unit
def test_solve_full_model_w_plugins():
    m = ConcreteModel()
    m.b = Block()
    m.b.v = Var()
    m.b.c = Constraint(expr=m.b.v == 12)

    m.initialization_order = [m, m.b]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    results = initializer._solve_full_model(m, log, "test_results")

    assert value(m.b.v) == pytest.approx(12, rel=1e-8)
    assert isinstance(results, dict)
    assert "Solution" in results


@pytest.mark.unit
def test_cleanup_none():
    m = ConcreteModel()
    m.initialization_order = [m]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    args = {"foo": "bar"}
    subinit = {}
    initializer._cleanup(m, args, subinit, log)

    # No-op - what to assert?


@pytest.mark.unit
def test_cleanup_no_args():
    m = ConcreteModel()
    m.b = Block()
    m.initialization_order = [m, m.b]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    args = {"foo": "bar", m.b: {}}
    subinit = {m.b: DummyInit()}
    initializer._cleanup(m, args, subinit, log)

    assert m.b._test2 == "default"


@pytest.mark.unit
def test_cleanup_w_args():
    m = ConcreteModel()
    m.b = Block()
    m.initialization_order = [m, m.b]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    args = {"foo": "bar", m.b: {"output": "forty-two"}}
    subinit = {m.b: DummyInit()}
    initializer._cleanup(m, args, subinit, log)

    assert m.b._test2 == "forty-two"


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
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer._init_props_0D(m, False)

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
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer._init_props_0D(m, True)

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
    log = initializer.get_logger(m)

    def estimate_states(block, *args, **kwargs):
        block._estimated = True

    m.control_volume.estimate_states = types.MethodType(
        estimate_states, m.control_volume
    )

    initializer.add_submodel_initializer(m.control_volume.properties, DummyInit())
    initializer._init_props_1D(m)

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
    initializer._init_rxns(m)

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
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())

    results = initializer._initialize_control_volume(m, False, log)

    assert isinstance(results, dict)
    assert "Solution" in results

    assert value(m.control_volume.properties_in.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v2) == pytest.approx(7, rel=1e-5)
    assert value(m.control_volume.reactions.v1) == pytest.approx(6, rel=1e-5)

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
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())

    results = initializer._initialize_control_volume(m, True, log)

    assert isinstance(results, dict)
    assert "Solution" in results

    assert value(m.control_volume.properties_in.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v2) == pytest.approx(7, rel=1e-5)
    assert value(m.control_volume.reactions.v1) == pytest.approx(6, rel=1e-5)

    assert m.control_volume.properties_in._initialized
    assert not hasattr(m.control_volume.properties_out, "_initialized")
    assert m.control_volume.reactions._initialized


# TODO: init props 1D


@pytest.mark.unit
def test_initialize_submodels_no_order():
    m = ConcreteModel()
    m.initialization_order = []

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    with pytest.raises(
        InitializationError,
        match="Main model \(unknown\) was not initialized \(no results returned\). "
        "This is likely due to an error in the model.initialization_order.",
    ):
        initializer._initialize_submodels(m, {}, False, {}, log)


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
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())
    initializer.add_submodel_initializer(m.plugin, DummyInit())

    initializer._initialize_submodels(
        m, {m.plugin: {}}, False, {m.plugin: DummyInit()}, log
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
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())
    initializer.add_submodel_initializer(m.plugin, DummyInit())

    initializer._initialize_submodels(
        m, {m.plugin: {}}, True, {m.plugin: DummyInit()}, log
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
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())

    initializer._initialize_submodels(
        m, {m.plugin: {}}, True, {m.plugin: DummyInit()}, log
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
