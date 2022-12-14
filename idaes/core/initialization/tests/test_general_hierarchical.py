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
        20.31609, rel=1e-6
    )
    assert value(m.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]) == pytest.approx(
        20.31609, rel=1e-6
    )
    assert value(m.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]) == pytest.approx(
        79.683910, rel=1e-6
    )
    assert value(m.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]) == pytest.approx(
        79.683910, rel=1e-6
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
    def addon_prepare(self, addon, output="default"):
        addon._test = output

    def addon_finalize(self, addon, output="default"):
        addon._test2 = output


@pytest.mark.unit
def test_prepare_addons_none():
    m = ConcreteModel()
    m.initialization_order = [m]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    args = {"foo": "bar"}
    subinit, addon_args = initializer._prepare_addons(m, args, log)

    assert subinit == {}
    assert addon_args == {"foo": "bar"}
    # Make sure we didn't change the original args
    assert args == {"foo": "bar"}


@pytest.mark.unit
def test_prepare_addons_no_args():
    m = ConcreteModel()
    m.b = Block()
    m.initialization_order = [m, m.b]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)
    initializer.add_submodel_initializer(m.b, DummyInit)

    args = {"foo": "bar"}
    subinit, addon_args = initializer._prepare_addons(m, args, log)

    assert len(subinit) == 1
    assert isinstance(subinit[m.b], DummyInit)
    assert addon_args == {"foo": "bar", m.b: {}}
    assert m.b._test == "default"
    # Make sure we didn't change the original args
    assert args == {"foo": "bar"}


@pytest.mark.unit
def test_prepare_addons_w_args():
    m = ConcreteModel()
    m.b = Block()
    m.initialization_order = [m, m.b]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)
    initializer.add_submodel_initializer(m.b, DummyInit)

    args = {"foo": "bar", m.b: {"output": "checkval"}}
    subinit, addon_args = initializer._prepare_addons(m, args, log)

    assert len(subinit) == 1
    assert isinstance(subinit[m.b], DummyInit)
    assert addon_args == {"foo": "bar", m.b: {"output": "checkval"}}
    assert m.b._test == "checkval"
    # Make sure we didn't change the original args
    assert args == {"foo": "bar", m.b: {"output": "checkval"}}


@pytest.mark.unit
def test_solve_full_model_no_addons():
    m = ConcreteModel()
    m.initialization_order = [m]

    initializer = SingleControlVolumeUnitInitializer()
    log = initializer.get_logger(m)

    results = initializer._solve_full_model(m, log, "test_results")

    assert results == "test_results"


@pytest.mark.unit
def test_solve_full_model_w_addons():
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
