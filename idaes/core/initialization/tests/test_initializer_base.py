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
Tests for InitializerBase class
"""
import pytest

from pyomo.environ import ConcreteModel

from idaes.core.initialization.initializer_base import InitializerBase

from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import InitializationError
from idaes.models.unit_models.cstr import CSTR
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)

__author__ = "Andrew Lee"


class TestBTSubMethods:
    @pytest.mark.unit
    def test_init(self):
        initializer = InitializerBase()

        assert initializer.postcheck_summary == {}

        assert hasattr(initializer, "config")

        assert "constraint_tolerance" in initializer.config

    @pytest.fixture
    def model(self):
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

        return m

    @pytest.mark.unit
    def test_postcheck(self, model):
        # TODO: Use a simpler test model
        initializer = InitializerBase()

        with pytest.raises(
            InitializationError,
            match=f"fs.unit failed to initialize successfully: uninitialized variables or "
            "unconverged equality constraints detected. Please check postcheck summary for "
            "more information.",
        ):
            initializer.postcheck(model.fs.unit)

        assert len(initializer.postcheck_summary["uninitialized_vars"]) == 0
        assert len(initializer.postcheck_summary["unconverged_constraints"]) == 0
