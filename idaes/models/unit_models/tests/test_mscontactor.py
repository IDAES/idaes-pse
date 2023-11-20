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
Tests for multi-state contactor unit model.
Authors: Andrew Lee
"""

import pytest
from types import MethodType

from pyomo.environ import (
    assert_optimal_termination,
    Block,
    ConcreteModel,
    Constraint,
    Expression,
    log,
    RangeSet,
    Set,
    TransformationFactory,
    units,
    value,
    Var,
)
from pyomo.network import Arc, Port
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from pyomo.dae import DerivativeVar

from idaes.core import (
    FlowsheetBlock,
    FlowDirection,
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    Component,
    Phase,
    MaterialFlowBasis,
    MaterialBalanceType,
)
from idaes.models.unit_models import Mixer, MixingType, MomentumMixingType
from idaes.models.unit_models.mscontactor import (
    MSContactor,
    MSContactorData,
    _get_state_blocks,
    MSContactorInitializer,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.misc import add_object_reference
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, PropertyNotSupportedError
from idaes.core.util.initialization import (
    propagate_state,
    fix_state_vars,
    revert_state_vars,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
    ReactionBlock,
)
from idaes.core.initialization import InitializationStatus
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)


# -----------------------------------------------------------------------------
# Property packages for testing
@declare_process_block_class("Parameters1")
class Parameter1Data(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent1 = Component()
        self.solute1 = Component()
        self.solute2 = Component()
        self.solute3 = Component()

        self._state_block_class = StateBlock1

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class SBlock1Base(StateBlock):
    def initialize(blk, **kwargs):
        pass

    def release_state(blk, **kwargs):
        pass


@declare_process_block_class("StateBlock1", block_class=SBlock1Base)
class StateBlock1Data(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super().build()

        self.flow_mol_phase_comp = Var(
            self.phase_component_set,
            units=units.mol / units.s,
        )
        self.enth_flow = Var(
            units=units.J / units.s,
        )
        self.pressure = Var(units=units.Pa)

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        return self.enth_flow

    def get_material_density_terms(self, p, j):
        return 42

    def get_energy_density_terms(self, p):
        return 43

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mol_phase_comp": self.flow_mol_phase_comp,
            "enth_flow": self.enth_flow,
            "pressure": self.pressure,
        }


@declare_process_block_class("Parameters2")
class Parameter2Data(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent2 = Component()
        self.solute1 = Component()
        self.solute2 = Component()

        self._state_block_class = StateBlock2

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class SBlock2Base(StateBlock):
    def initialize(blk, **kwargs):
        pass

    def release_state(blk, **kwargs):
        pass


@declare_process_block_class("StateBlock2", block_class=SBlock1Base)
class StateBlock2Data(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super().build()

        self.flow_mol_phase_comp = Var(
            self.phase_component_set,
            units=units.mol / units.s,
        )
        self.enth_flow = Var(
            units=units.J / units.s,
        )
        self.pressure = Var(units=units.Pa)

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_phase_comp[p, j]

    def get_enthalpy_flow_terms(self, p):
        return self.enth_flow

    def get_material_density_terms(self, p, j):
        return 52

    def get_energy_density_terms(self, p):
        return 53

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_mol_phase_comp": self.flow_mol_phase_comp,
            "enth_flow": self.enth_flow,
            "pressure": self.pressure,
        }


@declare_process_block_class("Parameters3")
class Parameter3Data(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent4 = Component()

        self._state_block_class = StateBlock3

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class SBlock3Base(StateBlock):
    def initialize(blk, **kwargs):
        pass

    def release_state(blk, **kwargs):
        pass


@declare_process_block_class("StateBlock3", block_class=SBlock1Base)
class StateBlock3Data(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super().build()

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass


@declare_process_block_class("Parameters4")
class Parameter3Data(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()
        self.phase2 = Phase()

        self.solvent1 = Component()
        self.solute1 = Component()
        self.solute2 = Component()
        self.solute3 = Component()

        self._state_block_class = StateBlock4

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


@declare_process_block_class("StateBlock4", block_class=SBlock1Base)
class StateBlock4Data(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super().build()

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar


# -----------------------------------------------------------------------------
# Frame class for unit testing
@declare_process_block_class("ECFrame")
class ECFrameData(MSContactorData):
    def build(self):
        super(MSContactorData, self).build()


# -----------------------------------------------------------------------------
class TestBuild:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()
        m.fs.properties2 = Parameters2()

        m.fs.unit = ECFrame(
            number_of_finite_elements=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {
                    "property_package": m.fs.properties2,
                    "flow_direction": FlowDirection.backward,
                },
            },
        )

        return m

    @pytest.fixture
    def dynamic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(
            dynamic=True,
            time_set=[0, 1],
            time_units=units.s,
        )

        m.fs.properties1 = Parameters1()
        m.fs.properties2 = Parameters2()

        m.fs.unit = ECFrame(
            number_of_finite_elements=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {
                    "property_package": m.fs.properties2,
                    "flow_direction": FlowDirection.backward,
                },
            },
        )

        return m

    @pytest.mark.unit
    def test_config(self, model):
        assert not model.fs.unit.config.dynamic
        assert not model.fs.unit.config.has_holdup
        assert model.fs.unit.config.number_of_finite_elements == 2
        assert "stream1" in model.fs.unit.config.streams
        assert "stream2" in model.fs.unit.config.streams
        assert model.fs.unit.config.interacting_streams is None

        assert (
            model.fs.unit.config.streams["stream1"].property_package
            is model.fs.properties1
        )
        assert model.fs.unit.config.streams["stream1"].property_package_args == {}
        assert model.fs.unit.config.streams["stream1"].reaction_package is None
        assert model.fs.unit.config.streams["stream1"].reaction_package_args == {}
        assert (
            model.fs.unit.config.streams["stream1"].flow_direction
            is FlowDirection.forward
        )
        assert model.fs.unit.config.streams["stream1"].has_feed
        assert not model.fs.unit.config.streams["stream1"].has_rate_reactions
        assert not model.fs.unit.config.streams["stream1"].has_equilibrium_reactions
        assert model.fs.unit.config.streams["stream1"].side_streams is None
        assert model.fs.unit.config.streams["stream1"].has_energy_balance
        assert model.fs.unit.config.streams["stream1"].has_pressure_balance
        assert not model.fs.unit.config.streams["stream1"].has_pressure_change

        assert (
            model.fs.unit.config.streams["stream2"].property_package
            is model.fs.properties2
        )
        assert model.fs.unit.config.streams["stream2"].property_package_args == {}
        assert model.fs.unit.config.streams["stream2"].reaction_package is None
        assert model.fs.unit.config.streams["stream2"].reaction_package_args == {}
        assert (
            model.fs.unit.config.streams["stream2"].flow_direction
            is FlowDirection.backward
        )
        assert model.fs.unit.config.streams["stream2"].has_feed
        assert not model.fs.unit.config.streams["stream2"].has_rate_reactions
        assert not model.fs.unit.config.streams["stream2"].has_equilibrium_reactions
        assert model.fs.unit.config.streams["stream2"].side_streams is None
        assert model.fs.unit.config.streams["stream2"].has_energy_balance
        assert model.fs.unit.config.streams["stream2"].has_pressure_balance
        assert not model.fs.unit.config.streams["stream2"].has_pressure_change

    @pytest.mark.unit
    def test_verify_inputs_too_few_streams(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()

        m.fs.unit = ECFrame(
            number_of_finite_elements=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
            },
        )

        with pytest.raises(
            ConfigurationError,
            match="MSContactor models must define at least two streams; received "
            "\['stream1'\]",
        ):
            m.fs.unit._verify_inputs()

    @pytest.mark.unit
    def test_verify_inputs_no_common_components(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()
        m.fs.properties3 = Parameters3()

        m.fs.unit = ECFrame(
            number_of_finite_elements=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {"property_package": m.fs.properties3},
            },
        )

        with pytest.raises(
            ConfigurationError,
            match="No common components found in property packages and no heterogeneous reactions "
            "specified. The MSContactor model assumes that mass transfer occurs between "
            "components with the same name in different streams or due to heterogeneous reactions.",
        ):
            m.fs.unit._verify_inputs()

        # Should pass if heterogeneous reaction argument provided
        m.fs.unit.config.heterogeneous_reactions = True
        m.fs.unit._verify_inputs()

    @pytest.mark.unit
    def test_verify_inputs_reactions_with_no_package(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()

        m.fs.unit = ECFrame(
            number_of_finite_elements=2,
            streams={
                "stream1": {
                    "property_package": m.fs.properties1,
                    "has_rate_reactions": True,
                },
                "stream2": {"property_package": m.fs.properties1},
            },
        )

        with pytest.raises(
            ConfigurationError,
            match="Stream stream1 was set to include reactions, "
            "but no reaction package was provided.",
        ):
            m.fs.unit._verify_inputs()

    @pytest.mark.unit
    def test_verify_inputs_construct_components(self, model):
        model.fs.unit._verify_inputs()

        assert isinstance(model.fs.unit.elements, RangeSet)
        assert len(model.fs.unit.elements) == 2

        assert isinstance(model.fs.unit.stream_interactions, Set)
        assert len(model.fs.unit.stream_interactions) == 1
        assert model.fs.unit.stream_interactions == [("stream1", "stream2")]

        assert isinstance(model.fs.unit.stream_component_interactions, Set)
        # One stream pair with two common components
        assert len(model.fs.unit.stream_component_interactions) == 2
        for k in model.fs.unit.stream_component_interactions:
            assert k in [
                ("stream1", "stream2", "solute1"),
                ("stream1", "stream2", "solute2"),
            ]

    @pytest.mark.unit
    def test_build_state_blocks(self, model):
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()

        assert flow_basis == MaterialFlowBasis.molar
        assert uom == model.fs.properties1.get_metadata().derived_units

        assert isinstance(model.fs.unit.stream1, StateBlock1)
        assert len(model.fs.unit.stream1) == 2
        assert not model.fs.unit.stream1[0, 1].config.defined_state
        assert not model.fs.unit.stream1[0, 2].config.defined_state

        assert isinstance(model.fs.unit.stream1_inlet_state, StateBlock1)
        assert len(model.fs.unit.stream1_inlet_state) == 1
        assert model.fs.unit.stream1_inlet_state[0].config.defined_state

        assert not hasattr(model.fs.unit, "stream1_side_stream_set")
        assert not hasattr(model.fs.unit, "stream1_side_stream_state")

        assert isinstance(model.fs.unit.stream2, StateBlock2)
        assert len(model.fs.unit.stream2) == 2
        assert not model.fs.unit.stream2[0, 1].config.defined_state
        assert not model.fs.unit.stream2[0, 2].config.defined_state

        assert isinstance(model.fs.unit.stream2_inlet_state, StateBlock2)
        assert len(model.fs.unit.stream2_inlet_state) == 1
        assert model.fs.unit.stream2_inlet_state[0].config.defined_state

        assert not hasattr(model.fs.unit, "stream2_side_stream_set")
        assert not hasattr(model.fs.unit, "stream2_side_stream_state")

    @pytest.mark.unit
    def test_build_state_blocks_no_feed(self, model):
        model.fs.unit.config.streams["stream2"].has_feed = False
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()

        assert flow_basis == MaterialFlowBasis.molar
        assert uom == model.fs.properties1.get_metadata().derived_units

        assert isinstance(model.fs.unit.stream1, StateBlock1)
        assert len(model.fs.unit.stream1) == 2
        assert not model.fs.unit.stream1[0, 1].config.defined_state
        assert not model.fs.unit.stream1[0, 2].config.defined_state

        assert isinstance(model.fs.unit.stream1_inlet_state, StateBlock1)
        assert len(model.fs.unit.stream1_inlet_state) == 1
        assert model.fs.unit.stream1_inlet_state[0].config.defined_state

        assert not hasattr(model.fs.unit, "stream1_side_stream_set")
        assert not hasattr(model.fs.unit, "stream1_side_stream_state")

        assert isinstance(model.fs.unit.stream2, StateBlock2)
        assert len(model.fs.unit.stream2) == 2
        assert not model.fs.unit.stream2[0, 1].config.defined_state
        assert not model.fs.unit.stream2[0, 2].config.defined_state

        assert not hasattr(model.fs.unit, "stream2_inlet_state")
        assert not hasattr(model.fs.unit, "stream2_side_stream_set")
        assert not hasattr(model.fs.unit, "stream2_side_stream_state")

    @pytest.mark.unit
    def test_build_state_blocks_side_stream(self, model):
        model.fs.unit.config.streams["stream2"].side_streams = [1]
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()

        assert flow_basis == MaterialFlowBasis.molar
        assert uom == model.fs.properties1.get_metadata().derived_units

        assert isinstance(model.fs.unit.stream1, StateBlock1)
        assert len(model.fs.unit.stream1) == 2
        assert not model.fs.unit.stream1[0, 1].config.defined_state
        assert not model.fs.unit.stream1[0, 2].config.defined_state

        assert isinstance(model.fs.unit.stream1_inlet_state, StateBlock1)
        assert len(model.fs.unit.stream1_inlet_state) == 1
        assert model.fs.unit.stream1_inlet_state[0].config.defined_state

        assert not hasattr(model.fs.unit, "stream1_side_stream_set")
        assert not hasattr(model.fs.unit, "stream1_side_stream_state")

        assert isinstance(model.fs.unit.stream2, StateBlock2)
        assert len(model.fs.unit.stream2) == 2
        assert not model.fs.unit.stream2[0, 1].config.defined_state
        assert not model.fs.unit.stream2[0, 2].config.defined_state

        assert isinstance(model.fs.unit.stream2_inlet_state, StateBlock2)
        assert len(model.fs.unit.stream2_inlet_state) == 1
        assert model.fs.unit.stream2_inlet_state[0].config.defined_state

        assert isinstance(model.fs.unit.stream2_side_stream_set, Set)
        assert len(model.fs.unit.stream2_side_stream_set) == 1
        assert isinstance(model.fs.unit.stream2_side_stream_state, StateBlock2)
        assert len(model.fs.unit.stream2_side_stream_state) == 1
        assert not model.fs.unit.stream2_side_stream_state[0, 1].config.defined_state

    @pytest.mark.unit
    def test_build_state_blocks_side_stream_invalid(self, model):
        model.fs.unit.config.streams["stream2"].side_streams = [10]
        model.fs.unit._verify_inputs()

        with pytest.raises(
            ConfigurationError,
            match="side_streams must be a sub-set of the set of elements. "
            "Found 10 in side_streams which is not in elements.",
        ):
            model.fs.unit._build_state_blocks()

    @pytest.mark.unit
    def test_build_state_blocks_different_flow_basis(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()
        m.fs.properties3 = Parameters3()

        m.fs.unit = ECFrame(
            number_of_finite_elements=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                # Properties3 has a different flow basis
                "stream2": {"property_package": m.fs.properties3},
            },
        )

        m.fs.unit.elements = Set(initialize=[1, 2])

        with pytest.raises(
            ConfigurationError,
            match="Property packages use different flow bases: ExtractionCascade "
            "requires all property packages to use the same basis. stream2 uses "
            "MaterialFlowBasis.mass, whilst first stream uses "
            "MaterialFlowBasis.molar.",
        ):
            m.fs.unit._build_state_blocks()

    @pytest.mark.unit
    def test_get_state_blocks(self, model):
        model.fs.unit._verify_inputs()
        model.fs.unit._build_state_blocks()

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 1, "stream1"
        )
        assert in_state is model.fs.unit.stream1_inlet_state[0]
        assert out_state is model.fs.unit.stream1[0, 1]
        assert side_state is None

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 2, "stream1"
        )
        assert in_state is model.fs.unit.stream1[0, 1]
        assert out_state is model.fs.unit.stream1[0, 2]
        assert side_state is None

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 1, "stream2"
        )
        assert in_state is model.fs.unit.stream2[0, 2]
        assert out_state is model.fs.unit.stream2[0, 1]
        assert side_state is None

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 2, "stream2"
        )
        assert in_state is model.fs.unit.stream2_inlet_state[0]
        assert out_state is model.fs.unit.stream2[0, 2]
        assert side_state is None

    @pytest.mark.unit
    def test_get_state_blocks_no_feed(self, model):
        model.fs.unit.config.streams["stream1"].has_feed = False
        model.fs.unit.config.streams["stream2"].has_feed = False
        model.fs.unit._verify_inputs()
        model.fs.unit._build_state_blocks()

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 1, "stream1"
        )
        assert in_state is None
        assert out_state is model.fs.unit.stream1[0, 1]
        assert side_state is None

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 2, "stream1"
        )
        assert in_state is model.fs.unit.stream1[0, 1]
        assert out_state is model.fs.unit.stream1[0, 2]
        assert side_state is None

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 1, "stream2"
        )
        assert in_state is model.fs.unit.stream2[0, 2]
        assert out_state is model.fs.unit.stream2[0, 1]
        assert side_state is None

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 2, "stream2"
        )
        assert in_state is None
        assert out_state is model.fs.unit.stream2[0, 2]
        assert side_state is None

    @pytest.mark.unit
    def test_get_state_blocks_side_streams(self, model):
        model.fs.unit.config.streams["stream1"].side_streams = [1]
        model.fs.unit.config.streams["stream2"].side_streams = [2]
        model.fs.unit._verify_inputs()
        model.fs.unit._build_state_blocks()

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 1, "stream1"
        )
        assert in_state is model.fs.unit.stream1_inlet_state[0]
        assert out_state is model.fs.unit.stream1[0, 1]
        assert side_state is model.fs.unit.stream1_side_stream_state[0, 1]

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 2, "stream1"
        )
        assert in_state is model.fs.unit.stream1[0, 1]
        assert out_state is model.fs.unit.stream1[0, 2]
        assert side_state is None

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 1, "stream2"
        )
        assert in_state is model.fs.unit.stream2[0, 2]
        assert out_state is model.fs.unit.stream2[0, 1]
        assert side_state is None

        in_state, out_state, side_state = _get_state_blocks(
            model.fs.unit, 0, 2, "stream2"
        )
        assert in_state is model.fs.unit.stream2_inlet_state[0]
        assert out_state is model.fs.unit.stream2[0, 2]
        assert side_state is model.fs.unit.stream2_side_stream_state[0, 2]

    @pytest.mark.unit
    def test_add_geometry_no_holdup(self, model):
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._add_geometry(uom)

        assert not hasattr(model.fs.unit, "volume")
        assert not hasattr(model.fs.unit, "volume_frac_stream")
        assert not hasattr(model.fs.unit, "sum_volume_frac")

        assert not hasattr(model.fs.unit, "stream1_phase_fraction")
        assert not hasattr(model.fs.unit, "stream1_sum_phase_fractions")

        assert not hasattr(model.fs.unit, "stream2_phase_fraction")
        assert not hasattr(model.fs.unit, "stream2_sum_phase_fractions")

    @pytest.mark.unit
    def test_add_geometry_holdup_single_phase(self, dynamic):
        dynamic.fs.unit._verify_inputs()
        flow_basis, uom = dynamic.fs.unit._build_state_blocks()
        dynamic.fs.unit._add_geometry(uom)

        assert isinstance(dynamic.fs.unit.volume, Var)
        assert len(dynamic.fs.unit.volume) == 2
        assert isinstance(dynamic.fs.unit.volume_frac_stream, Var)
        assert len(dynamic.fs.unit.volume_frac_stream) == 2 * 2 * 2
        assert isinstance(dynamic.fs.unit.sum_volume_frac, Constraint)
        assert len(dynamic.fs.unit.sum_volume_frac) == 2 * 2 * 1

        assert isinstance(dynamic.fs.unit.stream1_phase_fraction, Expression)
        assert isinstance(dynamic.fs.unit.stream2_phase_fraction, Expression)
        assert not hasattr(dynamic.fs.unit, "stream1_sum_phase_fractions")
        assert not hasattr(dynamic.fs.unit, "stream2_sum_phase_fractions")

        for i in dynamic.fs.unit.stream1_phase_fraction.values():
            assert i.expr == 1
        for i in dynamic.fs.unit.stream2_phase_fraction.values():
            assert i.expr == 1

    @pytest.mark.unit
    def test_add_geometry_holdup_multi_phase(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(
            dynamic=True,
            time_set=[0, 1],
            time_units=units.s,
        )

        m.fs.properties1 = Parameters1()
        m.fs.properties2 = Parameters4()

        m.fs.unit = ECFrame(
            number_of_finite_elements=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {
                    "property_package": m.fs.properties2,
                    "flow_direction": FlowDirection.backward,
                },
            },
        )

        m.fs.unit._verify_inputs()
        flow_basis, uom = m.fs.unit._build_state_blocks()
        m.fs.unit._add_geometry(uom)

        assert isinstance(m.fs.unit.volume, Var)
        assert len(m.fs.unit.volume) == 2
        assert isinstance(m.fs.unit.volume_frac_stream, Var)
        assert len(m.fs.unit.volume_frac_stream) == 2 * 2 * 2
        assert isinstance(m.fs.unit.sum_volume_frac, Constraint)
        assert len(m.fs.unit.sum_volume_frac) == 2 * 2 * 1

        assert isinstance(m.fs.unit.stream1_phase_fraction, Expression)
        assert isinstance(m.fs.unit.stream2_phase_fraction, Var)
        assert not hasattr(m.fs.unit, "stream1_sum_phase_fractions")
        assert isinstance(m.fs.unit.stream2_sum_phase_fractions, Constraint)

        for i in m.fs.unit.stream1_phase_fraction.values():
            assert i.expr == 1

        for (t, e), con in m.fs.unit.stream2_sum_phase_fractions.items():
            assert str(con.expr) == str(
                1
                == sum(
                    m.fs.unit.stream2_phase_fraction[t, e, p]
                    for p in ["phase1", "phase2"]
                )
            )

    @pytest.mark.unit
    def test_material_balances(self, model):
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert isinstance(model.fs.unit.material_transfer_term, Var)
        # One stream pair with two common components over two elements and 1 time point
        assert len(model.fs.unit.material_transfer_term) == 4
        assert_units_equivalent(
            model.fs.unit.material_transfer_term._units, units.mol / units.s
        )

        assert isinstance(model.fs.unit.stream1_material_balance, Constraint)
        # 1 time point, 2 elements, 4 components
        assert len(model.fs.unit.stream1_material_balance) == 8

        for j in ["solvent1", "solute3"]:  # no mass transfer, forward flow
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 2].flow_mol_phase_comp["phase1", j]
            )
        for j in ["solute1", "solute2"]:  # has +ve mass transfer, forward flow
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 2].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )

        assert isinstance(model.fs.unit.stream2_material_balance, Constraint)
        # 1 time point, 2 elements, 3 components
        assert len(model.fs.unit.stream2_material_balance) == 6
        for j in ["solvent2"]:  # no mass transfer, reverse flow
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 1].flow_mol_phase_comp["phase1", j]
            )
        for j in ["solute1", "solute2"]:  # has -ve mass transfer, reverse flow
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )

    @pytest.mark.unit
    def test_material_balances_dynamic(self, dynamic):
        dynamic.fs.unit._verify_inputs()
        flow_basis, uom = dynamic.fs.unit._build_state_blocks()
        dynamic.fs.unit._add_geometry(uom)
        dynamic.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert isinstance(dynamic.fs.unit.material_transfer_term, Var)
        # One stream pair with two common components over two elements and 2 time point
        assert len(dynamic.fs.unit.material_transfer_term) == 8
        assert_units_equivalent(
            dynamic.fs.unit.material_transfer_term._units, units.mol / units.s
        )

        assert isinstance(dynamic.fs.unit.stream1_material_holdup, Var)
        assert len(dynamic.fs.unit.stream1_material_holdup) == 16
        assert isinstance(dynamic.fs.unit.stream1_material_accumulation, DerivativeVar)
        assert len(dynamic.fs.unit.stream1_material_accumulation) == 16
        assert isinstance(
            dynamic.fs.unit.stream1_material_holdup_constraint, Constraint
        )
        assert len(dynamic.fs.unit.stream1_material_holdup_constraint) == 16
        for (
            t,
            x,
            p,
            j,
        ), con in dynamic.fs.unit.stream1_material_holdup_constraint.items():
            assert str(con.expr) == str(
                dynamic.fs.unit.stream1_material_holdup[t, x, p, j]
                == dynamic.fs.unit.volume[x]
                * dynamic.fs.unit.volume_frac_stream[t, x, "stream1"]
                * dynamic.fs.unit.stream1_phase_fraction[t, x, p]
                * 42
            )

        assert isinstance(dynamic.fs.unit.stream2_material_holdup, Var)
        assert len(dynamic.fs.unit.stream2_material_holdup) == 12
        assert isinstance(dynamic.fs.unit.stream2_material_accumulation, DerivativeVar)
        assert len(dynamic.fs.unit.stream2_material_accumulation) == 12
        assert isinstance(
            dynamic.fs.unit.stream2_material_holdup_constraint, Constraint
        )
        assert len(dynamic.fs.unit.stream2_material_holdup_constraint) == 12
        for (
            t,
            x,
            p,
            j,
        ), con in dynamic.fs.unit.stream2_material_holdup_constraint.items():
            assert str(con.expr) == str(
                dynamic.fs.unit.stream2_material_holdup[t, x, p, j]
                == dynamic.fs.unit.volume[x]
                * dynamic.fs.unit.volume_frac_stream[t, x, "stream2"]
                * dynamic.fs.unit.stream2_phase_fraction[t, x, p]
                * 52
            )

        assert isinstance(dynamic.fs.unit.stream1_material_balance, Constraint)
        # 2 time point, 2 elements, 4 components
        assert len(dynamic.fs.unit.stream1_material_balance) == 16

        for t in dynamic.fs.time:
            for j in ["solvent1", "solute3"]:  # no mass transfer, forward flow
                assert str(
                    dynamic.fs.unit.stream1_material_balance[t, 1, j].expr
                ) == str(
                    dynamic.fs.unit.stream1_material_accumulation[t, 1, "phase1", j]
                    == dynamic.fs.unit.stream1_inlet_state[t].flow_mol_phase_comp[
                        "phase1", j
                    ]
                    - dynamic.fs.unit.stream1[t, 1].flow_mol_phase_comp["phase1", j]
                )
                assert str(
                    dynamic.fs.unit.stream1_material_balance[t, 2, j].expr
                ) == str(
                    dynamic.fs.unit.stream1_material_accumulation[t, 2, "phase1", j]
                    == dynamic.fs.unit.stream1[t, 1].flow_mol_phase_comp["phase1", j]
                    - dynamic.fs.unit.stream1[t, 2].flow_mol_phase_comp["phase1", j]
                )
            for j in ["solute1", "solute2"]:  # has +ve mass transfer, forward flow
                assert str(
                    dynamic.fs.unit.stream1_material_balance[t, 1, j].expr
                ) == str(
                    dynamic.fs.unit.stream1_material_accumulation[t, 1, "phase1", j]
                    == dynamic.fs.unit.stream1_inlet_state[t].flow_mol_phase_comp[
                        "phase1", j
                    ]
                    - dynamic.fs.unit.stream1[t, 1].flow_mol_phase_comp["phase1", j]
                    + dynamic.fs.unit.material_transfer_term[
                        t, 1, "stream1", "stream2", j
                    ]
                )
                assert str(
                    dynamic.fs.unit.stream1_material_balance[t, 2, j].expr
                ) == str(
                    dynamic.fs.unit.stream1_material_accumulation[t, 2, "phase1", j]
                    == dynamic.fs.unit.stream1[t, 1].flow_mol_phase_comp["phase1", j]
                    - dynamic.fs.unit.stream1[t, 2].flow_mol_phase_comp["phase1", j]
                    + dynamic.fs.unit.material_transfer_term[
                        t, 2, "stream1", "stream2", j
                    ]
                )

        assert isinstance(dynamic.fs.unit.stream2_material_balance, Constraint)
        # 2 time point, 2 elements, 3 components
        assert len(dynamic.fs.unit.stream2_material_balance) == 12

        for t in dynamic.fs.time:
            for j in ["solvent2"]:  # no mass transfer, reverse flow
                assert str(
                    dynamic.fs.unit.stream2_material_balance[t, 2, j].expr
                ) == str(
                    dynamic.fs.unit.stream2_material_accumulation[t, 2, "phase1", j]
                    == dynamic.fs.unit.stream2_inlet_state[t].flow_mol_phase_comp[
                        "phase1", j
                    ]
                    - dynamic.fs.unit.stream2[t, 2].flow_mol_phase_comp["phase1", j]
                )
                assert str(
                    dynamic.fs.unit.stream2_material_balance[t, 1, j].expr
                ) == str(
                    dynamic.fs.unit.stream2_material_accumulation[t, 1, "phase1", j]
                    == dynamic.fs.unit.stream2[t, 2].flow_mol_phase_comp["phase1", j]
                    - dynamic.fs.unit.stream2[t, 1].flow_mol_phase_comp["phase1", j]
                )
            for j in ["solute1", "solute2"]:  # has -ve mass transfer, reverse flow
                assert str(
                    dynamic.fs.unit.stream2_material_balance[t, 2, j].expr
                ) == str(
                    dynamic.fs.unit.stream2_material_accumulation[t, 2, "phase1", j]
                    == dynamic.fs.unit.stream2_inlet_state[t].flow_mol_phase_comp[
                        "phase1", j
                    ]
                    - dynamic.fs.unit.stream2[t, 2].flow_mol_phase_comp["phase1", j]
                    - dynamic.fs.unit.material_transfer_term[
                        t, 2, "stream1", "stream2", j
                    ]
                )
                assert str(
                    dynamic.fs.unit.stream2_material_balance[t, 1, j].expr
                ) == str(
                    dynamic.fs.unit.stream2_material_accumulation[t, 1, "phase1", j]
                    == dynamic.fs.unit.stream2[t, 2].flow_mol_phase_comp["phase1", j]
                    - dynamic.fs.unit.stream2[t, 1].flow_mol_phase_comp["phase1", j]
                    - dynamic.fs.unit.material_transfer_term[
                        t, 1, "stream1", "stream2", j
                    ]
                )

    @pytest.mark.unit
    def test_build_material_balances_no_feed(self, model):
        model.fs.unit.config.streams["stream2"].has_feed = False
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert isinstance(model.fs.unit.material_transfer_term, Var)
        # One stream pair with two common components over two elements and 1 time point
        assert len(model.fs.unit.material_transfer_term) == 4
        assert_units_equivalent(
            model.fs.unit.material_transfer_term._units, units.mol / units.s
        )

        assert isinstance(model.fs.unit.stream1_material_balance, Constraint)
        # 1 time point, 2 elements, 4 components
        assert len(model.fs.unit.stream1_material_balance) == 8

        for j in ["solvent1", "solute3"]:  # no mass transfer, forward flow
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 2].flow_mol_phase_comp["phase1", j]
            )
        for j in ["solute1", "solute2"]:  # has +ve mass transfer, forward flow
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 2].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )

        assert isinstance(model.fs.unit.stream2_material_balance, Constraint)
        # 1 time point, 2 elements, 3 components
        assert len(model.fs.unit.stream2_material_balance) == 6
        for j in ["solvent2"]:  # no mass transfer, reverse flow
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == -model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 1].flow_mol_phase_comp["phase1", j]
            )
        for j in ["solute1", "solute2"]:  # has -ve mass transfer, reverse flow
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == -model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )

    @pytest.mark.unit
    def test_material_balances_side_stream(self, model):
        model.fs.unit.config.streams["stream2"].side_streams = [1]
        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert isinstance(model.fs.unit.material_transfer_term, Var)
        # One stream pair with two common components over two elements and 1 time point
        assert len(model.fs.unit.material_transfer_term) == 4
        assert_units_equivalent(
            model.fs.unit.material_transfer_term._units, units.mol / units.s
        )

        assert isinstance(model.fs.unit.stream1_material_balance, Constraint)
        # 1 time point, 2 elements, 4 components
        assert len(model.fs.unit.stream1_material_balance) == 8

        for j in ["solvent1", "solute3"]:  # no mass transfer, forward flow
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 2].flow_mol_phase_comp["phase1", j]
            )
        for j in ["solute1", "solute2"]:  # has +ve mass transfer, forward flow
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream1[0, 1].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream1[0, 2].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )

        assert isinstance(model.fs.unit.stream2_material_balance, Constraint)
        # 1 time point, 2 elements, 3 components
        assert len(model.fs.unit.stream2_material_balance) == 6
        for j in ["solvent2"]:  # no mass transfer, reverse flow
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 1].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.stream2_side_stream_state[0, 1].flow_mol_phase_comp[
                    "phase1", j
                ]
            )
        for j in ["solute1", "solute2"]:  # has -ve mass transfer, reverse flow
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2_inlet_state[0].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == model.fs.unit.stream2[0, 2].flow_mol_phase_comp["phase1", j]
                - model.fs.unit.stream2[0, 1].flow_mol_phase_comp["phase1", j]
                + model.fs.unit.stream2_side_stream_state[0, 1].flow_mol_phase_comp[
                    "phase1", j
                ]
                - model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )

    @pytest.mark.unit
    def test_energy_balances(self, model):
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_energy_balance_constraints(uom)

        assert isinstance(model.fs.unit.energy_transfer_term, Var)
        # 1 stream interaction, 2 elements
        assert len(model.fs.unit.energy_transfer_term) == 2
        for k in model.fs.unit.energy_transfer_term:
            assert k in [
                (0, 1, "stream1", "stream2"),
                (0, 2, "stream1", "stream2"),
            ]

        assert isinstance(model.fs.unit.stream1_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_energy_balance) == 2

        assert str(model.fs.unit.stream1_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].enth_flow
                - model.fs.unit.stream1[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream1_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1[0, 1].enth_flow
                - model.fs.unit.stream1[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )

        assert isinstance(model.fs.unit.stream2_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream2_energy_balance) == 2

        assert str(model.fs.unit.stream2_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream2_inlet_state[0].enth_flow
                - model.fs.unit.stream2[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream2_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream2[0, 2].enth_flow
                - model.fs.unit.stream2[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )

    @pytest.mark.unit
    def test_energy_balances_dynamic(self, dynamic):
        dynamic.fs.unit._verify_inputs()
        _, uom = dynamic.fs.unit._build_state_blocks()
        dynamic.fs.unit._add_geometry(uom)
        dynamic.fs.unit._build_energy_balance_constraints(uom)

        assert isinstance(dynamic.fs.unit.energy_transfer_term, Var)
        # 1 stream interaction, 2 elements
        assert len(dynamic.fs.unit.energy_transfer_term) == 4
        for k in dynamic.fs.unit.energy_transfer_term:
            assert k in [
                (0, 1, "stream1", "stream2"),
                (0, 2, "stream1", "stream2"),
                (1, 1, "stream1", "stream2"),
                (1, 2, "stream1", "stream2"),
            ]

        assert isinstance(dynamic.fs.unit.stream1_energy_holdup, Var)
        assert len(dynamic.fs.unit.stream1_energy_holdup) == 4
        assert isinstance(dynamic.fs.unit.stream1_energy_accumulation, DerivativeVar)
        assert len(dynamic.fs.unit.stream1_energy_accumulation) == 4
        assert isinstance(dynamic.fs.unit.stream1_energy_holdup_constraint, Constraint)
        assert len(dynamic.fs.unit.stream1_energy_holdup_constraint) == 4
        for (
            t,
            x,
            p,
        ), con in dynamic.fs.unit.stream1_energy_holdup_constraint.items():
            assert str(con.expr) == str(
                dynamic.fs.unit.stream1_energy_holdup[t, x, p]
                == dynamic.fs.unit.volume[x]
                * dynamic.fs.unit.volume_frac_stream[t, x, "stream1"]
                * dynamic.fs.unit.stream1_phase_fraction[t, x, p]
                * 43
            )

        assert isinstance(dynamic.fs.unit.stream2_energy_holdup, Var)
        assert len(dynamic.fs.unit.stream2_energy_holdup) == 4
        assert isinstance(dynamic.fs.unit.stream2_energy_accumulation, DerivativeVar)
        assert len(dynamic.fs.unit.stream2_energy_accumulation) == 4
        assert isinstance(dynamic.fs.unit.stream2_energy_holdup_constraint, Constraint)
        assert len(dynamic.fs.unit.stream2_energy_holdup_constraint) == 4
        for (
            t,
            x,
            p,
        ), con in dynamic.fs.unit.stream2_energy_holdup_constraint.items():
            assert str(con.expr) == str(
                dynamic.fs.unit.stream2_energy_holdup[t, x, p]
                == dynamic.fs.unit.volume[x]
                * dynamic.fs.unit.volume_frac_stream[t, x, "stream2"]
                * dynamic.fs.unit.stream2_phase_fraction[t, x, p]
                * 53
            )

        assert isinstance(dynamic.fs.unit.stream1_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(dynamic.fs.unit.stream1_energy_balance) == 4

        for t in dynamic.fs.time:
            assert str(dynamic.fs.unit.stream1_energy_balance[t, 1].expr) == str(
                dynamic.fs.unit.stream1_energy_accumulation[t, 1, "phase1"]
                == units.convert(
                    dynamic.fs.unit.stream1_inlet_state[t].enth_flow
                    - dynamic.fs.unit.stream1[t, 1].enth_flow,
                    units.kg * units.m**2 / units.s**3,
                )
                + dynamic.fs.unit.energy_transfer_term[t, 1, "stream1", "stream2"]
            )
            assert str(dynamic.fs.unit.stream1_energy_balance[t, 2].expr) == str(
                dynamic.fs.unit.stream1_energy_accumulation[t, 2, "phase1"]
                == units.convert(
                    dynamic.fs.unit.stream1[t, 1].enth_flow
                    - dynamic.fs.unit.stream1[t, 2].enth_flow,
                    units.kg * units.m**2 / units.s**3,
                )
                + dynamic.fs.unit.energy_transfer_term[t, 2, "stream1", "stream2"]
            )

        assert isinstance(dynamic.fs.unit.stream2_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(dynamic.fs.unit.stream2_energy_balance) == 4

        for t in dynamic.fs.time:
            assert str(dynamic.fs.unit.stream2_energy_balance[t, 2].expr) == str(
                dynamic.fs.unit.stream2_energy_accumulation[t, 2, "phase1"]
                == units.convert(
                    dynamic.fs.unit.stream2_inlet_state[t].enth_flow
                    - dynamic.fs.unit.stream2[t, 2].enth_flow,
                    units.kg * units.m**2 / units.s**3,
                )
                - dynamic.fs.unit.energy_transfer_term[t, 2, "stream1", "stream2"]
            )
            assert str(dynamic.fs.unit.stream2_energy_balance[t, 1].expr) == str(
                dynamic.fs.unit.stream2_energy_accumulation[t, 1, "phase1"]
                == units.convert(
                    dynamic.fs.unit.stream2[t, 2].enth_flow
                    - dynamic.fs.unit.stream2[t, 1].enth_flow,
                    units.kg * units.m**2 / units.s**3,
                )
                - dynamic.fs.unit.energy_transfer_term[t, 1, "stream1", "stream2"]
            )

    @pytest.mark.unit
    def test_energy_balances_has_heat_transfer(self, model):
        model.fs.unit.config.streams["stream2"].has_heat_transfer = True
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_energy_balance_constraints(uom)

        assert not hasattr(model.fs.unit, "stream1_heat")
        assert isinstance(model.fs.unit.stream1_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_energy_balance) == 2

        assert str(model.fs.unit.stream1_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].enth_flow
                - model.fs.unit.stream1[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream1_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1[0, 1].enth_flow
                - model.fs.unit.stream1[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )

        assert isinstance(model.fs.unit.stream2_heat, Var)
        assert len(model.fs.unit.stream2_heat) == 2
        assert_units_equivalent(model.fs.unit.stream2_heat, units.watt)

        assert isinstance(model.fs.unit.stream2_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream2_energy_balance) == 2

        assert str(model.fs.unit.stream2_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream2_inlet_state[0].enth_flow
                - model.fs.unit.stream2[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
            + model.fs.unit.stream2_heat[0, 2]
        )
        assert str(model.fs.unit.stream2_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream2[0, 2].enth_flow
                - model.fs.unit.stream2[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
            + model.fs.unit.stream2_heat[0, 1]
        )

    @pytest.mark.unit
    def test_energy_balances_no_feed(self, model):
        model.fs.unit.config.streams["stream2"].has_feed = False
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_energy_balance_constraints(uom)

        assert isinstance(model.fs.unit.stream1_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_energy_balance) == 2

        assert str(model.fs.unit.stream1_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].enth_flow
                - model.fs.unit.stream1[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream1_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1[0, 1].enth_flow
                - model.fs.unit.stream1[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )

        assert isinstance(model.fs.unit.stream2_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream2_energy_balance) == 2

        assert str(model.fs.unit.stream2_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                -model.fs.unit.stream2[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream2_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream2[0, 2].enth_flow
                - model.fs.unit.stream2[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )

    @pytest.mark.unit
    def test_energy_balances_has_energy_balance_false(self, model):
        model.fs.unit.config.streams["stream2"].has_energy_balance = False
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_energy_balance_constraints(uom)

        assert isinstance(model.fs.unit.stream1_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_energy_balance) == 2

        assert str(model.fs.unit.stream1_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].enth_flow
                - model.fs.unit.stream1[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream1_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1[0, 1].enth_flow
                - model.fs.unit.stream1[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )

        assert not hasattr(model.fs.unit, "stream2_energy_balance")

    @pytest.mark.unit
    def test_energy_balances_side_stream(self, model):
        model.fs.unit.config.streams["stream2"].side_streams = [1]
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_energy_balance_constraints(uom)

        assert isinstance(model.fs.unit.stream1_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_energy_balance) == 2

        assert str(model.fs.unit.stream1_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].enth_flow
                - model.fs.unit.stream1[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream1_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream1[0, 1].enth_flow
                - model.fs.unit.stream1[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )

        assert isinstance(model.fs.unit.stream2_energy_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream2_energy_balance) == 2

        assert str(model.fs.unit.stream2_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream2_inlet_state[0].enth_flow
                - model.fs.unit.stream2[0, 2].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream2_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                model.fs.unit.stream2[0, 2].enth_flow
                - model.fs.unit.stream2[0, 1].enth_flow
                + model.fs.unit.stream2_side_stream_state[0, 1].enth_flow,
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )

    @pytest.mark.unit
    def test_pressure_balances(self, model):
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_pressure_balance_constraints(uom)

        assert isinstance(model.fs.unit.stream1_pressure_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_pressure_balance) == 2

        assert str(model.fs.unit.stream1_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].pressure
                - model.fs.unit.stream1[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream1_pressure_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1[0, 1].pressure
                - model.fs.unit.stream1[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert isinstance(model.fs.unit.stream2_pressure_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream2_pressure_balance) == 2

        assert str(model.fs.unit.stream2_pressure_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream2_inlet_state[0].pressure
                - model.fs.unit.stream2[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream2_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream2[0, 2].pressure
                - model.fs.unit.stream2[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert not hasattr(model.fs.unit, "stream1_side_stream_pressure_balance")
        assert not hasattr(model.fs.unit, "stream2_side_stream_pressure_balance")

    @pytest.mark.unit
    def test_pressure_balances_deltaP(self, model):
        model.fs.unit.config.streams["stream2"].has_pressure_change = True
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_pressure_balance_constraints(uom)

        assert not hasattr(model.fs.unit, "stream1_deltaP")
        assert isinstance(model.fs.unit.stream1_pressure_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_pressure_balance) == 2

        assert str(model.fs.unit.stream1_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].pressure
                - model.fs.unit.stream1[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream1_pressure_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1[0, 1].pressure
                - model.fs.unit.stream1[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert isinstance(model.fs.unit.stream2_deltaP, Var)
        assert len(model.fs.unit.stream2_deltaP) == 2
        assert_units_equivalent(model.fs.unit.stream2_deltaP, units.Pa)

        assert isinstance(model.fs.unit.stream2_pressure_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream2_pressure_balance) == 2

        assert str(model.fs.unit.stream2_pressure_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream2_inlet_state[0].pressure
                - model.fs.unit.stream2[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
            + model.fs.unit.stream2_deltaP[0, 2]
        )
        assert str(model.fs.unit.stream2_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream2[0, 2].pressure
                - model.fs.unit.stream2[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
            + model.fs.unit.stream2_deltaP[0, 1]
        )

        assert not hasattr(model.fs.unit, "stream1_side_stream_pressure_balance")
        assert not hasattr(model.fs.unit, "stream2_side_stream_pressure_balance")

    @pytest.mark.unit
    def test_pressure_balances_no_feed(self, model):
        model.fs.unit.config.streams["stream2"].has_feed = False
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_pressure_balance_constraints(uom)

        assert isinstance(model.fs.unit.stream1_pressure_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_pressure_balance) == 2

        assert str(model.fs.unit.stream1_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].pressure
                - model.fs.unit.stream1[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream1_pressure_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1[0, 1].pressure
                - model.fs.unit.stream1[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert isinstance(model.fs.unit.stream2_pressure_balance, Constraint)
        # 1 time point, 1 elements; No balance at feed end
        assert len(model.fs.unit.stream2_pressure_balance) == 1

        assert str(model.fs.unit.stream2_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream2[0, 2].pressure
                - model.fs.unit.stream2[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert not hasattr(model.fs.unit, "stream1_side_stream_pressure_balance")
        assert not hasattr(model.fs.unit, "stream2_side_stream_pressure_balance")

    @pytest.mark.unit
    def test_pressure_balances_has_pressure_balance_false(self, model):
        model.fs.unit.config.streams["stream2"].has_pressure_balance = False
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_pressure_balance_constraints(uom)

        assert isinstance(model.fs.unit.stream1_pressure_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_pressure_balance) == 2

        assert str(model.fs.unit.stream1_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].pressure
                - model.fs.unit.stream1[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream1_pressure_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1[0, 1].pressure
                - model.fs.unit.stream1[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert not hasattr(model.fs.unit, "stream2_pressure_balance")

        assert not hasattr(model.fs.unit, "stream1_side_stream_pressure_balance")
        assert not hasattr(model.fs.unit, "stream2_side_stream_pressure_balance")

    @pytest.mark.unit
    def test_pressure_balances_side_stream(self, model):
        model.fs.unit.config.streams["stream2"].side_streams = [1]
        model.fs.unit._verify_inputs()
        _, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_pressure_balance_constraints(uom)

        assert isinstance(model.fs.unit.stream1_pressure_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream1_pressure_balance) == 2

        assert str(model.fs.unit.stream1_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1_inlet_state[0].pressure
                - model.fs.unit.stream1[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream1_pressure_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream1[0, 1].pressure
                - model.fs.unit.stream1[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert isinstance(model.fs.unit.stream2_pressure_balance, Constraint)
        # 1 time point, 2 elements
        assert len(model.fs.unit.stream2_pressure_balance) == 2

        assert str(model.fs.unit.stream2_pressure_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream2_inlet_state[0].pressure
                - model.fs.unit.stream2[0, 2].pressure,
                units.kg / units.m / units.s**2,
            )
        )
        assert str(model.fs.unit.stream2_pressure_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**-1 * units.s**-2)
            == units.convert(
                model.fs.unit.stream2[0, 2].pressure
                - model.fs.unit.stream2[0, 1].pressure,
                units.kg / units.m / units.s**2,
            )
        )

        assert not hasattr(model.fs.unit, "stream1_side_stream_pressure_balance")
        assert isinstance(
            model.fs.unit.stream2_side_stream_pressure_balance, Constraint
        )
        assert len(model.fs.unit.stream2_side_stream_pressure_balance) == 1
        assert str(
            model.fs.unit.stream2_side_stream_pressure_balance[0, 1].expr
        ) == str(
            model.fs.unit.stream2[0, 1].pressure
            == model.fs.unit.stream2_side_stream_state[0, 1].pressure
        )

    @pytest.mark.unit
    def test_ports(self, model):
        model.fs.unit._verify_inputs()
        model.fs.unit._build_state_blocks()
        model.fs.unit._build_ports()

        assert isinstance(model.fs.unit.stream1_inlet, Port)
        for p, j in model.fs.unit.stream1.phase_component_set:
            assert (
                model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream1_inlet.enth_flow[0]
            is model.fs.unit.stream1_inlet_state[0].enth_flow
        )
        assert (
            model.fs.unit.stream1_inlet.pressure[0]
            is model.fs.unit.stream1_inlet_state[0].pressure
        )

        assert isinstance(model.fs.unit.stream1_outlet, Port)
        for p, j in model.fs.unit.stream1.phase_component_set:
            assert (
                model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream1[0, 2].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream1_outlet.enth_flow[0]
            is model.fs.unit.stream1[0, 2].enth_flow
        )
        assert (
            model.fs.unit.stream1_outlet.pressure[0]
            is model.fs.unit.stream1[0, 2].pressure
        )

        assert isinstance(model.fs.unit.stream2_inlet, Port)
        for p, j in model.fs.unit.stream2.phase_component_set:
            assert (
                model.fs.unit.stream2_inlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream2_inlet_state[0].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream2_inlet.enth_flow[0]
            is model.fs.unit.stream2_inlet_state[0].enth_flow
        )
        assert (
            model.fs.unit.stream2_inlet.pressure[0]
            is model.fs.unit.stream2_inlet_state[0].pressure
        )

        assert isinstance(model.fs.unit.stream2_outlet, Port)
        for p, j in model.fs.unit.stream2.phase_component_set:
            assert (
                model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream2[0, 1].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream2_outlet.enth_flow[0]
            is model.fs.unit.stream2[0, 1].enth_flow
        )
        assert (
            model.fs.unit.stream2_outlet.pressure[0]
            is model.fs.unit.stream2[0, 1].pressure
        )

    @pytest.mark.unit
    def test_ports_no_feed(self, model):
        model.fs.unit.config.streams["stream2"].has_feed = False
        model.fs.unit._verify_inputs()
        model.fs.unit._build_state_blocks()
        model.fs.unit._build_ports()

        assert isinstance(model.fs.unit.stream1_inlet, Port)
        for p, j in model.fs.unit.stream1.phase_component_set:
            assert (
                model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream1_inlet_state[0].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream1_inlet.enth_flow[0]
            is model.fs.unit.stream1_inlet_state[0].enth_flow
        )
        assert (
            model.fs.unit.stream1_inlet.pressure[0]
            is model.fs.unit.stream1_inlet_state[0].pressure
        )

        assert isinstance(model.fs.unit.stream1_outlet, Port)
        for p, j in model.fs.unit.stream1.phase_component_set:
            assert (
                model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream1[0, 2].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream1_outlet.enth_flow[0]
            is model.fs.unit.stream1[0, 2].enth_flow
        )
        assert (
            model.fs.unit.stream1_outlet.pressure[0]
            is model.fs.unit.stream1[0, 2].pressure
        )

        assert not hasattr(model.fs.unit, "stream2_inlet")

        assert isinstance(model.fs.unit.stream2_outlet, Port)
        for p, j in model.fs.unit.stream2.phase_component_set:
            assert (
                model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, p, j]
                is model.fs.unit.stream2[0, 1].flow_mol_phase_comp[p, j]
            )
        assert (
            model.fs.unit.stream2_outlet.enth_flow[0]
            is model.fs.unit.stream2[0, 1].enth_flow
        )
        assert (
            model.fs.unit.stream2_outlet.pressure[0]
            is model.fs.unit.stream2[0, 1].pressure
        )


class TestReactions:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()
        m.fs.reactions = ReactionParameterTestBlock(property_package=m.fs.properties)

        m.fs.unit = ECFrame(
            number_of_finite_elements=2,
            streams={
                "stream1": {"property_package": m.fs.properties},
                "stream2": {
                    "property_package": m.fs.properties,
                    "flow_direction": FlowDirection.backward,
                },
            },
        )

        return m

    @pytest.mark.unit
    def test_inherent_reactions(self, model):
        # Activate inherent reactions for testing
        model.fs.properties._has_inherent_reactions = True

        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert isinstance(model.fs.unit.stream1_inherent_reaction_extent, Var)
        assert len(model.fs.unit.stream1_inherent_reaction_extent) == 4
        for k in model.fs.unit.stream1_inherent_reaction_extent:
            assert k in [(0, 1, "i1"), (0, 1, "i2"), (0, 2, "i1"), (0, 2, "i2")]

        assert isinstance(model.fs.unit.stream1_inherent_reaction_generation, Var)
        assert len(model.fs.unit.stream1_inherent_reaction_generation) == 8
        for k in model.fs.unit.stream1_inherent_reaction_generation:
            assert k in [
                (0, 1, "p1", "c1"),
                (0, 1, "p1", "c2"),
                (0, 1, "p2", "c1"),
                (0, 1, "p2", "c2"),
                (0, 2, "p1", "c1"),
                (0, 2, "p1", "c2"),
                (0, 2, "p2", "c1"),
                (0, 2, "p2", "c2"),
            ]

        assert isinstance(
            model.fs.unit.stream1_inherent_reaction_constraint, Constraint
        )
        assert len(model.fs.unit.stream1_inherent_reaction_constraint) == 8
        for k in model.fs.unit.stream1_inherent_reaction_constraint:
            assert k in [
                (0, 1, "p1", "c1"),
                (0, 1, "p1", "c2"),
                (0, 1, "p2", "c1"),
                (0, 1, "p2", "c2"),
                (0, 2, "p1", "c1"),
                (0, 2, "p1", "c2"),
                (0, 2, "p2", "c1"),
                (0, 2, "p2", "c2"),
            ]

        for j in [
            "c1",
            "c2",
        ]:  # has +ve mass transfer, forward flow, inherent reactions
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream1_inlet_state[0].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                + model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream1_inherent_reaction_generation[0, 1, p, j]
                    for p in ["p1", "p2"]
                )
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream1[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                + model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream1_inherent_reaction_generation[0, 2, p, j]
                    for p in ["p1", "p2"]
                )
            )

        assert isinstance(model.fs.unit.stream2_inherent_reaction_extent, Var)
        assert len(model.fs.unit.stream2_inherent_reaction_extent) == 4
        for k in model.fs.unit.stream2_inherent_reaction_extent:
            assert k in [(0, 1, "i1"), (0, 1, "i2"), (0, 2, "i1"), (0, 2, "i2")]

        assert isinstance(model.fs.unit.stream2_inherent_reaction_generation, Var)
        assert len(model.fs.unit.stream2_inherent_reaction_generation) == 8
        for k in model.fs.unit.stream2_inherent_reaction_generation:
            assert k in [
                (0, 1, "p1", "c1"),
                (0, 1, "p1", "c2"),
                (0, 1, "p2", "c1"),
                (0, 1, "p2", "c2"),
                (0, 2, "p1", "c1"),
                (0, 2, "p1", "c2"),
                (0, 2, "p2", "c1"),
                (0, 2, "p2", "c2"),
            ]

        assert isinstance(
            model.fs.unit.stream2_inherent_reaction_constraint, Constraint
        )
        assert len(model.fs.unit.stream2_inherent_reaction_constraint) == 8
        for k in model.fs.unit.stream2_inherent_reaction_constraint:
            assert k in [
                (0, 1, "p1", "c1"),
                (0, 1, "p1", "c2"),
                (0, 1, "p2", "c1"),
                (0, 1, "p2", "c2"),
                (0, 2, "p1", "c1"),
                (0, 2, "p1", "c2"),
                (0, 2, "p2", "c1"),
                (0, 2, "p2", "c2"),
            ]

        for j in [
            "c1",
            "c2",
        ]:  # has -ve mass transfer, forward flow, inherent reactions
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream2_inlet_state[0].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream2_inherent_reaction_generation[0, 2, p, j]
                    for p in ["p1", "p2"]
                )
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream2[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream2_inherent_reaction_generation[0, 1, p, j]
                    for p in ["p1", "p2"]
                )
            )

    @pytest.mark.unit
    def test_reaction_blocks(self, model):
        model.fs.unit.config.streams["stream1"].reaction_package = model.fs.reactions
        model.fs.unit.config.streams["stream2"].reaction_package = model.fs.reactions
        model.fs.unit.config.streams["stream2"].has_equilibrium_reactions = True

        model.fs.unit._verify_inputs()
        model.fs.unit._build_state_blocks()

        assert isinstance(model.fs.unit.stream1_reactions, ReactionBlock)
        assert len(model.fs.unit.stream1_reactions) == 2
        for k, b in model.fs.unit.stream1_reactions.items():
            assert k in [(0, 1), (0, 2)]
            assert not b.config.has_equilibrium

        assert isinstance(model.fs.unit.stream2_reactions, ReactionBlock)
        assert len(model.fs.unit.stream2_reactions) == 2
        for k, b in model.fs.unit.stream2_reactions.items():
            assert k in [(0, 1), (0, 2)]
            assert b.config.has_equilibrium

    @pytest.mark.unit
    def test_equilibrium_reactions(self, model):
        model.fs.unit.config.streams["stream2"].reaction_package = model.fs.reactions
        model.fs.unit.config.streams["stream2"].has_equilibrium_reactions = True

        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert not hasattr(model.fs.unit, "stream1_equilibrium_reaction_extent")
        assert not hasattr(model.fs.unit, "stream1_equilibrium_reaction_generation")
        assert not hasattr(model.fs.unit, "stream1_equilibrium_reaction_constraint")

        assert isinstance(model.fs.unit.stream2_equilibrium_reaction_extent, Var)
        assert len(model.fs.unit.stream2_equilibrium_reaction_extent) == 4
        for k in model.fs.unit.stream2_equilibrium_reaction_extent:
            assert k in [(0, 1, "e1"), (0, 1, "e2"), (0, 2, "e1"), (0, 2, "e2")]

        assert isinstance(model.fs.unit.stream2_equilibrium_reaction_generation, Var)
        assert len(model.fs.unit.stream2_equilibrium_reaction_generation) == 8
        for k in model.fs.unit.stream2_equilibrium_reaction_generation:
            assert k in [
                (0, 1, "p1", "c1"),
                (0, 1, "p1", "c2"),
                (0, 1, "p2", "c1"),
                (0, 1, "p2", "c2"),
                (0, 2, "p1", "c1"),
                (0, 2, "p1", "c2"),
                (0, 2, "p2", "c1"),
                (0, 2, "p2", "c2"),
            ]

        assert isinstance(
            model.fs.unit.stream2_equilibrium_reaction_constraint, Constraint
        )
        assert len(model.fs.unit.stream2_equilibrium_reaction_constraint) == 8
        for k in model.fs.unit.stream2_equilibrium_reaction_constraint:
            assert k in [
                (0, 1, "p1", "c1"),
                (0, 1, "p1", "c2"),
                (0, 1, "p2", "c1"),
                (0, 1, "p2", "c2"),
                (0, 2, "p1", "c1"),
                (0, 2, "p1", "c2"),
                (0, 2, "p2", "c1"),
                (0, 2, "p2", "c2"),
            ]

        for j in [
            "c1",
            "c2",
        ]:  # has +ve mass transfer, forward flow, no reactions
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream1_inlet_state[0].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                + model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream1[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                + model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )

        for j in [
            "c1",
            "c2",
        ]:  # has -ve mass transfer, forward flow, equilibrium reactions
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream2_inlet_state[0].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream2_equilibrium_reaction_generation[0, 2, p, j]
                    for p in ["p1", "p2"]
                )
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream2[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream2_equilibrium_reaction_generation[0, 1, p, j]
                    for p in ["p1", "p2"]
                )
            )

    @pytest.mark.unit
    def test_heterogeneous_reactions_no_build_method(self, model):
        model.fs.unit.config.heterogeneous_reactions = True

        model.fs.unit._verify_inputs()
        _, _ = model.fs.unit._build_state_blocks()

        with pytest.raises(
            ConfigurationError,
            match="Heterogeneous reaction package has not implemented a "
            "build_reaction_block method. Please ensure that your "
            "reaction block conforms to the required standards.",
        ):
            model.fs.unit._build_heterogeneous_reaction_blocks()

    @pytest.mark.unit
    def test_heterogeneous_reactions_no_rxn_index(self, model):
        model.hetero_dummy = Block()

        def build_reaction_block(*args, **kwargs):
            pass

        model.hetero_dummy.build_reaction_block = MethodType(
            build_reaction_block, model.hetero_dummy
        )

        model.fs.unit.config.heterogeneous_reactions = model.hetero_dummy

        model.fs.unit._verify_inputs()
        _, _ = model.fs.unit._build_state_blocks()

        with pytest.raises(
            PropertyNotSupportedError,
            match="Heterogeneous reaction package does not contain a list of "
            "reactions \(reaction_idx\).",
        ):
            model.fs.unit._build_heterogeneous_reaction_blocks()

    @pytest.mark.unit
    def test_heterogeneous_reactions(self, model):
        model.fs.hetero_dummy = Block()
        model.fs.hetero_dummy.reaction_idx = Set(initialize=["R1", "R2", "R3", "R4"])
        model.fs.hetero_dummy.reaction_stoichiometry = {
            ("R1", "p1", "c1"): 1,
            ("R2", "p1", "c2"): 1,
            ("R3", "p2", "c1"): 1,
            ("R4", "p2", "c2"): 1,
        }

        model.fs.unit.config.heterogeneous_reactions = model.fs.hetero_dummy

        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()

        model.fs.unit.heterogeneous_reactions = Block(
            model.fs.time,
            model.fs.unit.elements,
        )
        for e in model.fs.unit.elements:
            add_object_reference(
                model.fs.unit.heterogeneous_reactions[0, e],
                "params",
                model.fs.hetero_dummy,
            )

        model.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert isinstance(model.fs.unit.heterogeneous_reaction_extent, Var)
        for k in model.fs.unit.heterogeneous_reaction_extent.keys():
            assert k in [
                (0, 1, "R1"),
                (0, 1, "R2"),
                (0, 1, "R3"),
                (0, 1, "R4"),
                (0, 2, "R1"),
                (0, 2, "R2"),
                (0, 2, "R3"),
                (0, 2, "R4"),
            ]

        for s in ["stream1", "stream2"]:
            gen = getattr(model.fs.unit, s + "_heterogeneous_reactions_generation")
            assert isinstance(gen, Var)
            for k in gen:
                assert k in [
                    (0, 1, "p1", "c1"),
                    (0, 1, "p1", "c2"),
                    (0, 1, "p2", "c1"),
                    (0, 1, "p2", "c2"),
                    (0, 2, "p1", "c1"),
                    (0, 2, "p1", "c2"),
                    (0, 2, "p2", "c1"),
                    (0, 2, "p2", "c2"),
                ]

            con = getattr(model.fs.unit, s + "_heterogeneous_reaction_constraint")
            assert isinstance(con, Constraint)
            for k, c in con.items():
                assert k in [
                    (0, 1, "p1", "c1"),
                    (0, 1, "p1", "c2"),
                    (0, 1, "p2", "c1"),
                    (0, 1, "p2", "c2"),
                    (0, 2, "p1", "c1"),
                    (0, 2, "p1", "c2"),
                    (0, 2, "p2", "c1"),
                    (0, 2, "p2", "c2"),
                ]

                if k[2] == "p1" and k[3] == "c1":
                    r = "R1"
                elif k[2] == "p1" and k[3] == "c2":
                    r = "R2"
                elif k[2] == "p2" and k[3] == "c1":
                    r = "R3"
                else:
                    r = "R4"

                expr = str(
                    gen[k] - model.fs.unit.heterogeneous_reaction_extent[0, k[1], r]
                )
                assert str(c.body) == expr

        for j in [
            "c1",
            "c2",
        ]:  # has +ve mass transfer, forward flow, heterogeneous reactions
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream1_inlet_state[0].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                + model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream1_heterogeneous_reactions_generation[0, 1, p, j]
                    for p in ["p1", "p2"]
                )
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream1[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                + model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream1_heterogeneous_reactions_generation[0, 2, p, j]
                    for p in ["p1", "p2"]
                )
            )

        for j in [
            "c1",
            "c2",
        ]:  # has -ve mass transfer, forward flow, heterogeneous reactions
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream2_inlet_state[0].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream2_heterogeneous_reactions_generation[0, 2, p, j]
                    for p in ["p1", "p2"]
                )
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream2[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream2_heterogeneous_reactions_generation[0, 1, p, j]
                    for p in ["p1", "p2"]
                )
            )

    @pytest.mark.unit
    def test_rate_reactions(self, model):
        model.fs.unit.config.streams["stream2"].reaction_package = model.fs.reactions
        model.fs.unit.config.streams["stream2"].has_rate_reactions = True

        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)

        assert not hasattr(model.fs.unit, "stream1_rate_reaction_extent")
        assert not hasattr(model.fs.unit, "stream1_rate_reaction_generation")
        assert not hasattr(model.fs.unit, "stream1_rate_reaction_constraint")

        assert isinstance(model.fs.unit.stream2_rate_reaction_extent, Var)
        assert len(model.fs.unit.stream2_rate_reaction_extent) == 4
        for k in model.fs.unit.stream2_rate_reaction_extent:
            assert k in [(0, 1, "r1"), (0, 1, "r2"), (0, 2, "r1"), (0, 2, "r2")]

        assert isinstance(model.fs.unit.stream2_rate_reaction_generation, Var)
        assert len(model.fs.unit.stream2_rate_reaction_generation) == 8
        for k in model.fs.unit.stream2_rate_reaction_generation:
            assert k in [
                (0, 1, "p1", "c1"),
                (0, 1, "p1", "c2"),
                (0, 1, "p2", "c1"),
                (0, 1, "p2", "c2"),
                (0, 2, "p1", "c1"),
                (0, 2, "p1", "c2"),
                (0, 2, "p2", "c1"),
                (0, 2, "p2", "c2"),
            ]

        assert isinstance(model.fs.unit.stream2_rate_reaction_constraint, Constraint)
        assert len(model.fs.unit.stream2_rate_reaction_constraint) == 8
        for k in model.fs.unit.stream2_rate_reaction_constraint:
            assert k in [
                (0, 1, "p1", "c1"),
                (0, 1, "p1", "c2"),
                (0, 1, "p2", "c1"),
                (0, 1, "p2", "c2"),
                (0, 2, "p1", "c1"),
                (0, 2, "p1", "c2"),
                (0, 2, "p2", "c1"),
                (0, 2, "p2", "c2"),
            ]

        for j in [
            "c1",
            "c2",
        ]:  # has +ve mass transfer, forward flow, no reactions
            assert str(model.fs.unit.stream1_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream1_inlet_state[0].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                + model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
            )
            assert str(model.fs.unit.stream1_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream1[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                + model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
            )

        for j in [
            "c1",
            "c2",
        ]:  # has -ve mass transfer, forward flow, rate reactions
            assert str(model.fs.unit.stream2_material_balance[0, 2, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream2_inlet_state[0].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - model.fs.unit.material_transfer_term[0, 2, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream2_rate_reaction_generation[0, 2, p, j]
                    for p in ["p1", "p2"]
                )
            )
            assert str(model.fs.unit.stream2_material_balance[0, 1, j].expr) == str(
                0 * (units.mol * units.s**-1)
                == sum(
                    model.fs.unit.stream2[0, 2].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 1].get_material_flow_terms(p, j)
                    for p in ["p1", "p2"]
                )
                - model.fs.unit.material_transfer_term[0, 1, "stream1", "stream2", j]
                + sum(
                    model.fs.unit.stream2_rate_reaction_generation[0, 1, p, j]
                    for p in ["p1", "p2"]
                )
            )

    @pytest.mark.unit
    def test_heat_of_reaction_rate(self, model):
        model.fs.unit.config.streams["stream2"].reaction_package = model.fs.reactions
        model.fs.unit.config.streams["stream2"].has_rate_reactions = True
        model.fs.unit.config.streams["stream2"].has_heat_of_reaction = True

        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)
        model.fs.unit._build_energy_balance_constraints(uom)

        assert str(model.fs.unit.stream1_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                sum(
                    model.fs.unit.stream1_inlet_state[0].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 1].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                ),
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream1_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                sum(
                    model.fs.unit.stream1[0, 1].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 2].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                ),
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )

        assert str(model.fs.unit.stream2_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                sum(
                    model.fs.unit.stream2_inlet_state[0].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 2].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                ),
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
            - sum(
                model.fs.unit.stream2_rate_reaction_extent[0, 2, r]
                * model.fs.unit.stream2_reactions[0, 2].dh_rxn[r]
                for r in ["r1", "r2"]
            )
        )
        assert str(model.fs.unit.stream2_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                sum(
                    model.fs.unit.stream2[0, 2].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 1].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                ),
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
            - sum(
                model.fs.unit.stream2_rate_reaction_extent[0, 1, r]
                * model.fs.unit.stream2_reactions[0, 1].dh_rxn[r]
                for r in ["r1", "r2"]
            )
        )

    @pytest.mark.unit
    def test_heat_of_reaction_equilibrium(self, model):
        model.fs.unit.config.streams["stream2"].reaction_package = model.fs.reactions
        model.fs.unit.config.streams["stream2"].has_equilibrium_reactions = True
        model.fs.unit.config.streams["stream2"].has_heat_of_reaction = True

        model.fs.unit._verify_inputs()
        flow_basis, uom = model.fs.unit._build_state_blocks()
        model.fs.unit._build_material_balance_constraints(flow_basis, uom)
        model.fs.unit._build_energy_balance_constraints(uom)

        assert str(model.fs.unit.stream1_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                sum(
                    model.fs.unit.stream1_inlet_state[0].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 1].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                ),
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
        )
        assert str(model.fs.unit.stream1_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                sum(
                    model.fs.unit.stream1[0, 1].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream1[0, 2].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                ),
                units.kg * units.m**2 / units.s**3,
            )
            + model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
        )

        assert str(model.fs.unit.stream2_energy_balance[0, 2].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                sum(
                    model.fs.unit.stream2_inlet_state[0].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 2].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                ),
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 2, "stream1", "stream2"]
            - sum(
                model.fs.unit.stream2_equilibrium_reaction_extent[0, 2, e]
                * model.fs.unit.stream2_reactions[0, 2].dh_rxn[e]
                for e in ["e1", "e2"]
            )
        )
        assert str(model.fs.unit.stream2_energy_balance[0, 1].expr) == str(
            0 * (units.kg * units.m**2 * units.s**-3)
            == units.convert(
                sum(
                    model.fs.unit.stream2[0, 2].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                )
                - sum(
                    model.fs.unit.stream2[0, 1].get_enthalpy_flow_terms(p)
                    for p in ["p1", "p2"]
                ),
                units.kg * units.m**2 / units.s**3,
            )
            - model.fs.unit.energy_transfer_term[0, 1, "stream1", "stream2"]
            - sum(
                model.fs.unit.stream2_equilibrium_reaction_extent[0, 1, e]
                * model.fs.unit.stream2_reactions[0, 1].dh_rxn[e]
                for e in ["e1", "e2"]
            )
        )


class TestToyProblem:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties1 = Parameters1()
        m.fs.properties2 = Parameters2()

        m.fs.unit = MSContactor(
            number_of_finite_elements=2,
            streams={
                "stream1": {"property_package": m.fs.properties1},
                "stream2": {
                    "property_package": m.fs.properties2,
                    "flow_direction": FlowDirection.backward,
                },
            },
        )

        return m

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        # Expect 17 DoF:
        # 6 stream1 inlets, 5 stream2 inlets
        # 2x2 mass transfer terms and 2 energy transfer terms
        assert (degrees_of_freedom(model)) == 17

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_toy_problem(self, model):
        model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, "phase1", "solvent1"].fix(2)
        model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, "phase1", "solute1"].fix(3)
        model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, "phase1", "solute2"].fix(4)
        model.fs.unit.stream1_inlet.flow_mol_phase_comp[0, "phase1", "solute3"].fix(5)
        model.fs.unit.stream1_inlet.enth_flow[0].fix(5000)
        model.fs.unit.stream1_inlet.pressure[0].fix(1.2e5)

        model.fs.unit.stream2_inlet.flow_mol_phase_comp[0, "phase1", "solvent2"].fix(11)
        model.fs.unit.stream2_inlet.flow_mol_phase_comp[0, "phase1", "solute1"].fix(12)
        model.fs.unit.stream2_inlet.flow_mol_phase_comp[0, "phase1", "solute2"].fix(13)
        model.fs.unit.stream2_inlet.enth_flow[0].fix(7000)
        model.fs.unit.stream2_inlet.pressure[0].fix(2e5)

        model.fs.unit.material_transfer_term[0, :, "stream1", "stream2", "solute1"].fix(
            0.5
        )
        model.fs.unit.material_transfer_term[0, :, "stream1", "stream2", "solute2"].fix(
            -0.5
        )

        model.fs.unit.energy_transfer_term[0, :, "stream1", "stream2"].fix(100)

        assert (degrees_of_freedom(model)) == 0

        solver = get_solver()
        results = solver.solve(model)

        assert_optimal_termination(results)

        assert value(
            model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, "phase1", "solvent1"]
        ) == pytest.approx(2, rel=1e-5)
        assert value(
            model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, "phase1", "solute1"]
        ) == pytest.approx(
            4, rel=1e-5
        )  # 3 + 0.5 + 0.5
        assert value(
            model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, "phase1", "solute2"]
        ) == pytest.approx(
            3, rel=1e-5
        )  # 3 - 0.5 - 0.5
        assert value(
            model.fs.unit.stream1_outlet.flow_mol_phase_comp[0, "phase1", "solute3"]
        ) == pytest.approx(5, rel=1e-5)

        assert value(
            model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, "phase1", "solvent2"]
        ) == pytest.approx(11, rel=1e-5)
        assert value(
            model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, "phase1", "solute1"]
        ) == pytest.approx(
            11, rel=1e-5
        )  # 12 - 0.5 - 0.5
        assert value(
            model.fs.unit.stream2_outlet.flow_mol_phase_comp[0, "phase1", "solute2"]
        ) == pytest.approx(
            14, rel=1e-5
        )  # 13 + 0.5 + 0.5

        assert value(model.fs.unit.stream1_outlet.enth_flow[0]) == pytest.approx(
            5200, rel=1e-5
        )
        assert value(model.fs.unit.stream2_outlet.enth_flow[0]) == pytest.approx(
            6800, rel=1e-5
        )

        assert value(model.fs.unit.stream1_outlet.pressure[0]) == pytest.approx(
            1.2e5, rel=1e-5
        )
        assert value(model.fs.unit.stream2_outlet.pressure[0]) == pytest.approx(
            2e5, rel=1e-5
        )


# -----------------------------------------------------------------------------
# Li-Co Diafiltration example
@declare_process_block_class("LiCoParameters")
class LiCoParameterData(PhysicalParameterBlock):
    def build(self):
        super().build()

        self.phase1 = Phase()

        self.solvent = Component()
        self.Li = Component()
        self.Co = Component()

        self._state_block_class = LiCoStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.hour,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class LiCoSBlockBase(StateBlock):
    def initialize(blk, *args, hold_state=False, **kwargs):
        flags = fix_state_vars(blk, {})

        if hold_state is True:
            return flags
        else:
            blk.release_state(flags)

    def release_state(blk, flags, **kwargs):
        if flags is None:
            return
        # Unfix state variables
        revert_state_vars(blk, flags)


@declare_process_block_class("LiCoStateBlock", block_class=LiCoSBlockBase)
class LiCoStateBlock1Data(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=units.m**3 / units.hour,
            bounds=(1e-8, None),
        )
        self.conc_mass_solute = Var(
            ["Li", "Co"],
            units=units.kg / units.m**3,
            bounds=(1e-8, None),
        )

    def get_material_flow_terms(self, p, j):
        if j == "solvent":
            # Assume constant density of pure water
            return self.flow_vol * 1000 * units.kg / units.m**3
        else:
            return self.flow_vol * self.conc_mass_solute[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_solute": self.conc_mass_solute,
        }


class TestMSContactorInitializer:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        # Add separation stages
        m.fs.contactor = MSContactor(
            number_of_finite_elements=10,
            streams={
                "s1": {
                    "property_package": m.fs.properties,
                },
                "s2": {
                    "property_package": m.fs.properties,
                },
            },
        )

        m.fs.contactor.s1_inlet.flow_vol[0].set_value(1.0e-03)
        m.fs.contactor.s1_inlet.conc_mol_comp[0, "H2O"].set_value(55388.0)
        m.fs.contactor.s1_inlet.conc_mol_comp[0, "NaOH"].set_value(100.0)
        m.fs.contactor.s1_inlet.conc_mol_comp[0, "EthylAcetate"].set_value(100.0)
        m.fs.contactor.s1_inlet.conc_mol_comp[0, "SodiumAcetate"].set_value(0.0)
        m.fs.contactor.s1_inlet.conc_mol_comp[0, "Ethanol"].set_value(0.0)
        m.fs.contactor.s1_inlet.temperature[0].set_value(303.15)
        m.fs.contactor.s1_inlet.pressure[0].set_value(101325.0)

        m.fs.contactor.s2_inlet.flow_vol[0].set_value(2.0e-03)
        m.fs.contactor.s2_inlet.conc_mol_comp[0, "H2O"].set_value(55388.0)
        m.fs.contactor.s2_inlet.conc_mol_comp[0, "NaOH"].set_value(50.0)
        m.fs.contactor.s2_inlet.conc_mol_comp[0, "EthylAcetate"].set_value(50.0)
        m.fs.contactor.s2_inlet.conc_mol_comp[0, "SodiumAcetate"].set_value(50.0)
        m.fs.contactor.s2_inlet.conc_mol_comp[0, "Ethanol"].set_value(50.0)
        m.fs.contactor.s2_inlet.temperature[0].set_value(323.15)
        m.fs.contactor.s2_inlet.pressure[0].set_value(2e5)

        m.fs.contactor.material_transfer_term.fix(0)
        m.fs.contactor.energy_transfer_term.fix(0)

        return m

    @pytest.mark.unit
    def test_default_initializer(self, model):
        assert MSContactorData.default_initializer is MSContactorInitializer
        assert model.fs.contactor.default_initializer is MSContactorInitializer

    @pytest.mark.component
    def test_MSInitializer(self, model):
        initializer = MSContactorInitializer()
        initializer.initialize(model.fs.contactor)

        assert (
            initializer.summary[model.fs.contactor]["status"] == InitializationStatus.Ok
        )

        for x in model.fs.contactor.elements:
            assert value(model.fs.contactor.s1[0, x].flow_vol) == pytest.approx(
                1.0e-03, rel=1e-6
            )
            assert value(
                model.fs.contactor.s1[0, x].conc_mol_comp["H2O"]
            ) == pytest.approx(55388, rel=1e-6)
            assert value(
                model.fs.contactor.s1[0, x].conc_mol_comp["NaOH"]
            ) == pytest.approx(1e2, rel=1e-6)
            assert value(
                model.fs.contactor.s1[0, x].conc_mol_comp["EthylAcetate"]
            ) == pytest.approx(1e2, rel=1e-6)
            assert value(
                model.fs.contactor.s1[0, x].conc_mol_comp["SodiumAcetate"]
            ) == pytest.approx(0, abs=1e-4)
            assert value(
                model.fs.contactor.s1[0, x].conc_mol_comp["Ethanol"]
            ) == pytest.approx(0, abs=1e-4)
            assert value(model.fs.contactor.s1[0, x].temperature) == pytest.approx(
                303.15, rel=1e-6
            )
            assert value(model.fs.contactor.s1[0, x].pressure) == pytest.approx(
                101325, rel=1e-6
            )

            assert value(model.fs.contactor.s2[0, x].flow_vol) == pytest.approx(
                2.0e-03, rel=1e-6
            )
            assert value(
                model.fs.contactor.s2[0, x].conc_mol_comp["H2O"]
            ) == pytest.approx(55388, rel=1e-6)
            assert value(
                model.fs.contactor.s2[0, x].conc_mol_comp["NaOH"]
            ) == pytest.approx(50, rel=1e-6)
            assert value(
                model.fs.contactor.s2[0, x].conc_mol_comp["EthylAcetate"]
            ) == pytest.approx(50, rel=1e-6)
            assert value(
                model.fs.contactor.s2[0, x].conc_mol_comp["SodiumAcetate"]
            ) == pytest.approx(50, rel=1e-6)
            assert value(
                model.fs.contactor.s2[0, x].conc_mol_comp["Ethanol"]
            ) == pytest.approx(50, rel=1e-6)
            assert value(model.fs.contactor.s2[0, x].temperature) == pytest.approx(
                323.15, rel=1e-6
            )
            assert value(model.fs.contactor.s2[0, x].pressure) == pytest.approx(
                2e5, rel=1e-6
            )

        assert not model.fs.contactor.s1_inlet.flow_vol[0].fixed
        assert not model.fs.contactor.s1_inlet.conc_mol_comp[0, "H2O"].fixed
        assert not model.fs.contactor.s1_inlet.conc_mol_comp[0, "NaOH"].fixed
        assert not model.fs.contactor.s1_inlet.conc_mol_comp[0, "EthylAcetate"].fixed
        assert not model.fs.contactor.s1_inlet.conc_mol_comp[0, "SodiumAcetate"].fixed
        assert not model.fs.contactor.s1_inlet.conc_mol_comp[0, "Ethanol"].fixed
        assert not model.fs.contactor.s1_inlet.temperature[0].fixed
        assert not model.fs.contactor.s1_inlet.pressure[0].fixed

        assert not model.fs.contactor.s2_inlet.flow_vol[0].fixed
        assert not model.fs.contactor.s2_inlet.conc_mol_comp[0, "H2O"].fixed
        assert not model.fs.contactor.s2_inlet.conc_mol_comp[0, "NaOH"].fixed
        assert not model.fs.contactor.s2_inlet.conc_mol_comp[0, "EthylAcetate"].fixed
        assert not model.fs.contactor.s2_inlet.conc_mol_comp[0, "SodiumAcetate"].fixed
        assert not model.fs.contactor.s2_inlet.conc_mol_comp[0, "Ethanol"].fixed
        assert not model.fs.contactor.s2_inlet.temperature[0].fixed
        assert not model.fs.contactor.s2_inlet.pressure[0].fixed

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, model):
        perf_dict = model.fs.contactor._get_performance_contents()

        assert perf_dict == {}

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, model):
        stable = model.fs.contactor._get_stream_table_contents()

        expected = {
            "Units": {
                "Volumetric Flowrate": getattr(units.pint_registry, "m**3/second"),
                "Molar Concentration H2O": getattr(units.pint_registry, "mole/m**3"),
                "Molar Concentration NaOH": getattr(units.pint_registry, "mole/m**3"),
                "Molar Concentration EthylAcetate": getattr(
                    units.pint_registry, "mole/m**3"
                ),
                "Molar Concentration SodiumAcetate": getattr(
                    units.pint_registry, "mole/m**3"
                ),
                "Molar Concentration Ethanol": getattr(
                    units.pint_registry, "mole/m**3"
                ),
                "Temperature": getattr(units.pint_registry, "K"),
                "Pressure": getattr(units.pint_registry, "Pa"),
            },
            "s1 Inlet": {
                "Volumetric Flowrate": pytest.approx(0.001, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(5.5388e4, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(0, abs=1e-6),
                "Molar Concentration Ethanol": pytest.approx(0, abs=1e-6),
                "Temperature": pytest.approx(303.15, rel=1e-4),
                "Pressure": pytest.approx(101325, rel=1e-4),
            },
            "s1 Outlet": {
                "Volumetric Flowrate": pytest.approx(1, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(100, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(100, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(100, rel=1e-4),
                "Temperature": pytest.approx(298.15, rel=1e-4),
                "Pressure": pytest.approx(101325, rel=1e-4),
            },
            "s2 Inlet": {
                "Volumetric Flowrate": pytest.approx(0.002, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(5.5388e4, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(50, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(50, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(50, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(50, rel=1e-4),
                "Temperature": pytest.approx(323.15, rel=1e-4),
                "Pressure": pytest.approx(2e5, rel=1e-4),
            },
            "s2 Outlet": {
                "Volumetric Flowrate": pytest.approx(1, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(100, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(100, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(100, rel=1e-4),
                "Temperature": pytest.approx(298.15, rel=1e-4),
                "Pressure": pytest.approx(101325, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected


class TestLiCODiafiltration:
    """
    Test case based on:

    Wamble, N.P., Eugene, E.A., Phillip, W.A., Dowling, A.W.,
    'Optimal Diafiltration Membrane Cascades Enable Green Recycling
    of Spent Lithium-Ion Batteries',
    ACS Sustainable Chem. Eng. 2022, 10, 12207−12225

    Configuration and results based on Figure 2, Case III
    """

    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = LiCoParameters()

        # Add separation stages
        m.fs.stage1 = MSContactor(
            number_of_finite_elements=10,
            streams={
                "retentate": {
                    "property_package": m.fs.properties,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
                "permeate": {
                    "property_package": m.fs.properties,
                    "has_feed": False,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
            },
        )

        m.fs.stage2 = MSContactor(
            number_of_finite_elements=10,
            streams={
                "retentate": {
                    "property_package": m.fs.properties,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
                "permeate": {
                    "property_package": m.fs.properties,
                    "has_feed": False,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
            },
        )

        m.fs.stage3 = MSContactor(
            number_of_finite_elements=10,
            streams={
                "retentate": {
                    "property_package": m.fs.properties,
                    "side_streams": [10],
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
                "permeate": {
                    "property_package": m.fs.properties,
                    "has_feed": False,
                    "has_energy_balance": False,
                    "has_pressure_balance": False,
                },
            },
        )

        # Add mixers
        m.fs.mix1 = Mixer(
            num_inlets=2,
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            energy_mixing_type=MixingType.none,
            momentum_mixing_type=MomentumMixingType.none,
        )
        m.fs.mix2 = Mixer(
            num_inlets=2,
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            energy_mixing_type=MixingType.none,
            momentum_mixing_type=MomentumMixingType.none,
        )

        # Connect units
        m.fs.stream1 = Arc(
            source=m.fs.stage1.permeate_outlet,
            destination=m.fs.mix1.inlet_2,
        )
        m.fs.stream2 = Arc(
            source=m.fs.mix1.outlet,
            destination=m.fs.stage2.retentate_inlet,
        )
        m.fs.stream3 = Arc(
            source=m.fs.stage2.permeate_outlet,
            destination=m.fs.mix2.inlet_2,
        )
        m.fs.stream4 = Arc(
            source=m.fs.mix2.outlet,
            destination=m.fs.stage3.retentate_inlet,
        )
        m.fs.stream5 = Arc(
            source=m.fs.stage2.retentate_outlet,
            destination=m.fs.stage1.retentate_inlet,
        )
        m.fs.stream6 = Arc(
            source=m.fs.stage3.retentate_outlet,
            destination=m.fs.mix1.inlet_1,
        )

        TransformationFactory("network.expand_arcs").apply_to(m)

        # Global constants
        J = 0.1 * units.m / units.hour
        w = 1.5 * units.m
        rho = 1000 * units.kg / units.m**3

        # Add mass transfer variables and constraints
        m.fs.solutes = Set(initialize=["Li", "Co"])

        m.fs.sieving_coefficient = Var(
            m.fs.solutes,
            units=units.dimensionless,
        )
        m.fs.sieving_coefficient["Li"].fix(1.3)
        m.fs.sieving_coefficient["Co"].fix(0.5)

        m.fs.stage1.length = Var(units=units.m)
        m.fs.stage2.length = Var(units=units.m)
        m.fs.stage3.length = Var(units=units.m)

        # Start by initializing with a short length
        # Too long and we deplete solvent due to lack of recycles
        m.fs.stage1.length.fix(10)
        m.fs.stage2.length.fix(10)
        m.fs.stage3.length.fix(10)

        def solvent_rule(b, s):
            return (
                b.material_transfer_term[0, s, "permeate", "retentate", "solvent"]
                == J * b.length * w * rho / 10
            )

        def solute_rule(b, s, j):
            if s == 1:
                in_state = b.retentate_inlet_state[0]
            else:
                sp = b.elements.prev(s)
                in_state = b.retentate[0, sp]

            return log(b.retentate[0, s].conc_mass_solute[j]) + (
                m.fs.sieving_coefficient[j] - 1
            ) * log(in_state.flow_vol) == log(in_state.conc_mass_solute[j]) + (
                m.fs.sieving_coefficient[j] - 1
            ) * log(
                b.retentate[0, s].flow_vol
            )

        m.fs.stage1.solvent_flux = Constraint(
            m.fs.stage1.elements,
            rule=solvent_rule,
        )
        m.fs.stage1.solute_sieving = Constraint(
            m.fs.stage1.elements,
            m.fs.solutes,
            rule=solute_rule,
        )

        m.fs.stage2.solvent_flux = Constraint(
            m.fs.stage2.elements,
            rule=solvent_rule,
        )
        m.fs.stage2.solute_sieving = Constraint(
            m.fs.stage2.elements,
            m.fs.solutes,
            rule=solute_rule,
        )

        # For stage 3, we need to account for the side feed at element 10
        def stage3_solute_rule(b, s, j):
            if s == 1:
                q_in = b.retentate_inlet_state[0].flow_vol
                c_in = b.retentate_inlet_state[0].conc_mass_solute[j]
            elif s == 10:
                sp = b.elements.prev(s)
                q_in = (
                    b.retentate[0, sp].flow_vol
                    + b.retentate_side_stream_state[0, 10].flow_vol
                )
                c_in = (
                    b.retentate[0, sp].conc_mass_solute[j] * b.retentate[0, sp].flow_vol
                    + b.retentate_side_stream_state[0, 10].conc_mass_solute[j]
                    * b.retentate_side_stream_state[0, 10].flow_vol
                ) / q_in
            else:
                sp = b.elements.prev(s)
                q_in = b.retentate[0, sp].flow_vol
                c_in = b.retentate[0, sp].conc_mass_solute[j]

            return log(b.retentate[0, s].conc_mass_solute[j]) + (
                m.fs.sieving_coefficient[j] - 1
            ) * log(q_in) == log(c_in) + (m.fs.sieving_coefficient[j] - 1) * log(
                b.retentate[0, s].flow_vol
            )

        m.fs.stage3.solvent_flux = Constraint(
            m.fs.stage3.elements,
            rule=solvent_rule,
        )
        m.fs.stage3.solute_sieving = Constraint(
            m.fs.stage3.elements,
            m.fs.solutes,
            rule=stage3_solute_rule,
        )

        return m

    @pytest.mark.component
    def test_diafiltration_build(self, model):
        # TODO: More checks here
        assert isinstance(model.fs.stage3.retentate_inlet, Port)
        assert isinstance(model.fs.stage3.retentate_outlet, Port)
        assert not hasattr(model.fs.stage3, "permeate_inlet")
        assert isinstance(model.fs.stage3.permeate_outlet, Port)

    @pytest.mark.integration
    def test_initialize_and_solve(self, model):
        # Start with stage 3
        # Initial feed guess is pure diafiltrate (no recycle)
        model.fs.stage3.retentate_inlet.flow_vol[0].fix(30)
        model.fs.stage3.retentate_inlet.conc_mass_solute[0, "Li"].fix(0.1)
        model.fs.stage3.retentate_inlet.conc_mass_solute[0, "Co"].fix(0.2)

        # Fresh feed stream is side feed at element 10
        model.fs.stage3.retentate_side_stream_state[0, 10].flow_vol.fix(100)
        model.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Li"].fix(
            1.7
        )
        model.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Co"].fix(
            17
        )

        # Initialize flow and conc values to avoid log(0)
        for s in model.fs.stage3.retentate.values():
            s.flow_vol.set_value(30)
            s.conc_mass_solute["Li"].set_value(0.1)
            s.conc_mass_solute["Co"].set_value(0.2)

        assert degrees_of_freedom(model.fs.stage3) == 0

        initializer = MSContactorInitializer()
        initializer.initialize(model.fs.stage3)

        # Unfix feed guesses
        model.fs.stage3.retentate_inlet.flow_vol[0].unfix()
        model.fs.stage3.retentate_inlet.conc_mass_solute[0, "Li"].unfix()
        model.fs.stage3.retentate_inlet.conc_mass_solute[0, "Co"].unfix()

        # Stage 2 next - feed is retentate of stage 3
        Q = value(model.fs.stage3.retentate_outlet.flow_vol[0])
        C_Li = value(model.fs.stage3.retentate_outlet.conc_mass_solute[0, "Li"])
        C_Co = value(model.fs.stage3.retentate_outlet.conc_mass_solute[0, "Co"])

        model.fs.stage2.retentate_inlet.flow_vol[0].fix(Q)
        model.fs.stage2.retentate_inlet.conc_mass_solute[0, "Li"].fix(C_Li)
        model.fs.stage2.retentate_inlet.conc_mass_solute[0, "Co"].fix(C_Co)

        # Initialize flow and conc values to avoid log(0)
        for s in model.fs.stage2.retentate.values():
            s.flow_vol.set_value(Q)
            s.conc_mass_solute["Li"].set_value(C_Li)
            s.conc_mass_solute["Co"].set_value(C_Co)

        assert degrees_of_freedom(model.fs.stage2) == 0

        initializer.initialize(model.fs.stage2)

        # Unfix feed guesses
        model.fs.stage2.retentate_inlet.flow_vol[0].unfix()
        model.fs.stage2.retentate_inlet.conc_mass_solute[0, "Li"].unfix()
        model.fs.stage2.retentate_inlet.conc_mass_solute[0, "Co"].unfix()

        # Initialize Mixer 2
        # Inlet 1 is fresh diafilatrate
        model.fs.mix2.inlet_1.flow_vol[0].fix(30)
        model.fs.mix2.inlet_1.conc_mass_solute[0, "Li"].fix(0.1)
        model.fs.mix2.inlet_1.conc_mass_solute[0, "Co"].fix(0.2)

        propagate_state(
            destination=model.fs.mix2.inlet_2,
            source=model.fs.stage2.permeate_outlet,
        )

        model.fs.mix2.initialize()  # TODO: Update to Initializer object

        # Initialize first stage - feed is retentate of stage 2
        Q = value(model.fs.stage2.retentate_outlet.flow_vol[0])
        C_Li = value(model.fs.stage2.retentate_outlet.conc_mass_solute[0, "Li"])
        C_Co = value(model.fs.stage2.retentate_outlet.conc_mass_solute[0, "Co"])

        model.fs.stage1.retentate_inlet.flow_vol[0].fix(Q)
        model.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].fix(C_Li)
        model.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].fix(C_Co)

        # Initialize flow and conc values to avoid log(0)
        for s in model.fs.stage1.retentate.values():
            s.flow_vol.set_value(Q)
            s.conc_mass_solute["Li"].set_value(C_Li)
            s.conc_mass_solute["Co"].set_value(C_Co)

        assert degrees_of_freedom(model.fs.stage1) == 0

        initializer.initialize(model.fs.stage1)

        # Unfix feed guesses
        model.fs.stage1.retentate_inlet.flow_vol[0].unfix()
        model.fs.stage1.retentate_inlet.conc_mass_solute[0, "Li"].unfix()
        model.fs.stage1.retentate_inlet.conc_mass_solute[0, "Co"].unfix()

        # Initialize Mixer 1
        propagate_state(
            destination=model.fs.mix1.inlet_1,
            source=model.fs.stage3.retentate_outlet,
        )

        propagate_state(
            destination=model.fs.mix1.inlet_2,
            source=model.fs.stage1.permeate_outlet,
        )

        model.fs.mix1.initialize()  # TODO: Update to Initializer object

        # Solve the full model
        assert degrees_of_freedom(model) == 0

        solver = get_solver()
        res = solver.solve(model, tee=True)
        assert_optimal_termination(res)

        # Increase stage length and re-solve
        L = 756.4  # isotropic stages
        model.fs.stage1.length.fix(L)
        model.fs.stage2.length.fix(L)
        model.fs.stage3.length.fix(L)

        res = solver.solve(model, tee=True)
        assert_optimal_termination(res)

        # Check conservation
        # Solvent
        assert value(
            model.fs.mix2.inlet_1.flow_vol[0]
            + model.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
        ) == pytest.approx(
            value(
                model.fs.stage3.permeate_outlet.flow_vol[0]
                + model.fs.stage1.retentate_outlet.flow_vol[0]
            ),
            rel=1e-5,
        )
        # Lithium
        assert value(
            model.fs.mix2.inlet_1.flow_vol[0]
            * model.fs.mix2.inlet_1.conc_mass_solute[0, "Li"]
            + model.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * model.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Li"]
        ) == pytest.approx(
            value(
                model.fs.stage3.permeate_outlet.flow_vol[0]
                * model.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"]
                + model.fs.stage1.retentate_outlet.flow_vol[0]
                * model.fs.stage1.retentate_outlet.conc_mass_solute[0, "Li"]
            ),
            rel=1e-5,
        )
        # Cobalt
        assert value(
            model.fs.mix2.inlet_1.flow_vol[0]
            * model.fs.mix2.inlet_1.conc_mass_solute[0, "Co"]
            + model.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
            * model.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute["Co"]
        ) == pytest.approx(
            value(
                model.fs.stage3.permeate_outlet.flow_vol[0]
                * model.fs.stage3.permeate_outlet.conc_mass_solute[0, "Co"]
                + model.fs.stage1.retentate_outlet.flow_vol[0]
                * model.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"]
            ),
            rel=1e-5,
        )

        # Calculate recovery
        R_Li = value(
            model.fs.stage3.permeate_outlet.flow_vol[0]
            * model.fs.stage3.permeate_outlet.conc_mass_solute[0, "Li"]
            / (
                model.fs.mix2.inlet_1.flow_vol[0]
                * model.fs.mix2.inlet_1.conc_mass_solute[0, "Li"]
                + model.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
                * model.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute[
                    "Li"
                ]
            )
        )
        R_Co = value(
            model.fs.stage1.retentate_outlet.flow_vol[0]
            * model.fs.stage1.retentate_outlet.conc_mass_solute[0, "Co"]
            / (
                model.fs.mix2.inlet_1.flow_vol[0]
                * model.fs.mix2.inlet_1.conc_mass_solute[0, "Co"]
                + model.fs.stage3.retentate_side_stream_state[0, 10].flow_vol
                * model.fs.stage3.retentate_side_stream_state[0, 10].conc_mass_solute[
                    "Co"
                ]
            )
        )
        assert R_Li == pytest.approx(0.9451, rel=1e-4)
        assert R_Co == pytest.approx(0.6378, rel=1e-4)
