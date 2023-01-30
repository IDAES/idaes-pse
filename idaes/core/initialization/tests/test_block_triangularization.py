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
Tests for Block Triangularization initialization
"""
import pytest
import types

from pyomo.environ import ConcreteModel, Constraint, value, Var

from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.initialization.initializer_base import InitializationStatus

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

    initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
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


class TestBTSubMethods:
    @pytest.mark.unit
    def test_config(self):
        initializer = BlockTriangularizationInitializer()

        assert hasattr(initializer, "config")

        assert "constraint_tolerance" in initializer.config

        assert "block_solver" in initializer.config
        assert "block_solver_options" in initializer.config
        assert "calculate_variable_options" in initializer.config

    # TODO: Tests for prechecks and initialization_routine stand alone

    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.v1 = Var()
        m.v2 = Var()
        m.v3 = Var()
        m.v4 = Var()

        m.c1 = Constraint(expr=m.v1 == m.v2)
        m.c2 = Constraint(expr=2 * m.v2 == m.v3 + m.v4)
        m.c3 = Constraint(expr=m.v3 - m.v4 == 0)

        # Add a dummy method for fixing initialization states
        def fix_initialization_states(blk):
            blk.v1.fix(4)

        m.fix_initialization_states = types.MethodType(fix_initialization_states, m)

        return m

    @pytest.mark.component
    def test_workflow(self, model):
        initializer = BlockTriangularizationInitializer()

        status = initializer.initialize(model)

        assert model.v1.value == 4
        assert model.v2.value == 4
        assert model.v3.value == 4
        assert model.v4.value == 4

        assert not model.v1.fixed

        assert status == InitializationStatus.Ok
