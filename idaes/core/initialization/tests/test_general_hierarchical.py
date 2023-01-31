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
Tests for general hierarchical initialization routines
"""
import pytest
import types

from pyomo.environ import Block, ConcreteModel, Constraint, value, Var

from idaes.core.initialization.general_hierarchical import (
    SingleControlVolumeUnitInitializer,
)
from idaes.core.initialization.initializer_base import InitializationStatus
import idaes.logger as idaeslog
from idaes.core.util.exceptions import InitializationError

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.cstr import CSTR
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)

__author__ = "Andrew Lee"


@pytest.mark.integration
def test_workflow():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = SaponificationParameterBlock()
    m.fs.reactions = SaponificationReactionParameterBlock(
        property_package=m.fs.properties
    )

    m.fs.unit = CSTR(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        has_equilibrium_reactions=False,
        has_heat_transfer=True,
        has_heat_of_reaction=True,
        has_pressure_change=True,
    )

    m.fs.unit.inlet.flow_vol.fix(1.0e-03)
    m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.unit.inlet.temperature.fix(303.15)
    m.fs.unit.inlet.pressure.fix(101325.0)

    m.fs.unit.volume.fix(1.5e-03)
    m.fs.unit.heat_duty.fix(0)
    m.fs.unit.deltaP.fix(0)

    initializer = SingleControlVolumeUnitInitializer()
    initializer.initialize(m.fs.unit)

    assert initializer.summary[m.fs.unit]["status"] == InitializationStatus.Ok

    assert value(m.fs.unit.outlet.flow_vol[0]) == pytest.approx(1e-3, rel=1e-6)
    assert value(m.fs.unit.outlet.conc_mol_comp[0, "H2O"]) == pytest.approx(
        55388, rel=1e-6
    )
    assert value(m.fs.unit.outlet.conc_mol_comp[0, "NaOH"]) == pytest.approx(
        20.316235, rel=1e-6
    )
    assert value(m.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]) == pytest.approx(
        20.316235, rel=1e-6
    )
    assert value(m.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]) == pytest.approx(
        79.684995, rel=1e-6
    )
    assert value(m.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]) == pytest.approx(
        79.684995, rel=1e-6
    )
    assert value(m.fs.unit.outlet.temperature[0]) == pytest.approx(304.0856, rel=1e-6)
    assert value(m.fs.unit.outlet.pressure[0]) == pytest.approx(101325, rel=1e-6)


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

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.properties_in, DummyInit())
    initializer.add_submodel_initializer(m.control_volume.properties_out, DummyInit())
    initializer._init_props_0D(m, False)

    assert m.control_volume.properties_in._initialized
    assert m.control_volume.properties_out._initialized


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

    initializer.add_submodel_initializer(m.control_volume.properties, DummyInit())
    initializer._init_props_1D(m)

    assert m.control_volume.properties._initialized


@pytest.mark.unit
def test_init_rxns():
    m = ConcreteModel()
    m.control_volume = Block()
    m.control_volume.reactions = Block()

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    initializer.add_submodel_initializer(m.control_volume.reactions, DummyInit())
    initializer._init_rxns(m)

    assert m.control_volume.reactions._initialized


@pytest.mark.unit
def test_initialize_main_model_no_copy():
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

    results = initializer._initialize_main_model(m, False, log)

    assert isinstance(results, dict)
    assert "Solution" in results

    assert value(m.control_volume.properties_in.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v1) == pytest.approx(12, rel=1e-5)
    assert value(m.control_volume.properties_out.v2) == pytest.approx(7, rel=1e-5)
    assert value(m.control_volume.reactions.v1) == pytest.approx(6, rel=1e-5)

    assert m.control_volume.properties_in._initialized
    assert m.control_volume.properties_out._initialized
    assert m.control_volume.reactions._initialized


@pytest.mark.unit
def test_initialize_main_model_copy():
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

    results = initializer._initialize_main_model(m, True, log)

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
