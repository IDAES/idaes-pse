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
Tests for ControlVolumeBlockData.

Author: Andrew Lee
"""
import pytest
import pandas

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    RangeSet,
    Set,
    Var,
    value,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from pyomo.network import Port
from pyomo.common.config import ConfigBlock

from idaes.core import (
    FlowsheetBlock,
    declare_process_block_class,
    StateBlockData,
    StateBlock,
    PhysicalParameterBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    Phase,
    Component,
)
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties import iapws95
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)

from idaes.models.unit_models.mixer import (
    Mixer,
    MixerData,
    MixingType,
    MomentumMixingType,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    InitializationError,
    PropertyNotSupportedError,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    TestStateBlock,
    initialization_tester,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

# TODO: Should have a test for this that does not requrie models_extra
from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
# Unit Tests for Mixer
class TestMixer(object):
    @declare_process_block_class("MixerFrame")
    class MixerFrameData(MixerData):
        def build(self):
            super(MixerData, self).build()

    @pytest.fixture(scope="function")
    def mixer_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.mix = MixerFrame(property_package=m.fs.pp)

        return m

    @pytest.mark.unit
    def test_mixer_config(self, mixer_frame):
        assert len(mixer_frame.fs.mix.config) == 12
        assert mixer_frame.fs.mix.config.dynamic is False
        assert mixer_frame.fs.mix.config.has_holdup is False
        assert mixer_frame.fs.mix.config.property_package == mixer_frame.fs.pp
        assert isinstance(mixer_frame.fs.mix.config.property_package_args, ConfigBlock)
        assert len(mixer_frame.fs.mix.config.property_package_args) == 0
        assert mixer_frame.fs.mix.config.inlet_list is None
        assert mixer_frame.fs.mix.config.num_inlets is None
        assert mixer_frame.fs.mix.config.has_phase_equilibrium is False
        assert mixer_frame.fs.mix.config.energy_mixing_type == MixingType.extensive
        assert (
            mixer_frame.fs.mix.config.momentum_mixing_type
            == MomentumMixingType.minimize
        )
        assert mixer_frame.fs.mix.config.mixed_state_block is None
        assert mixer_frame.fs.mix.config.construct_ports is True
        assert (
            mixer_frame.fs.mix.config.material_balance_type
            == MaterialBalanceType.useDefault
        )

    @pytest.mark.unit
    def test_inherited_methods(self, mixer_frame):
        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        assert hasattr(mixer_frame.fs.mix.config.property_package, "phase_list")

    @pytest.mark.unit
    def test_create_inlet_list_default(self, mixer_frame):
        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()

        for i in inlet_list:
            assert i in ["inlet_1", "inlet_2"]

    @pytest.mark.unit
    def test_create_inlet_list_inlet_list(self, mixer_frame):
        mixer_frame.fs.mix.config.inlet_list = ["foo", "bar"]

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()

        for i in inlet_list:
            assert i in ["foo", "bar"]

    @pytest.mark.unit
    def test_create_inlet_list_num_inlets(self, mixer_frame):
        mixer_frame.fs.mix.config.num_inlets = 3

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()

        for i in inlet_list:
            assert i in ["inlet_1", "inlet_2", "inlet_3"]

    @pytest.mark.unit
    def test_create_inlet_list_both_args_consistent(self, mixer_frame):
        mixer_frame.fs.mix.config.inlet_list = ["foo", "bar"]
        mixer_frame.fs.mix.config.num_inlets = 2

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()

        for i in inlet_list:
            assert i in ["foo", "bar"]

    @pytest.mark.unit
    def test_create_inlet_list_both_args_inconsistent(self, mixer_frame):
        mixer_frame.fs.mix.config.inlet_list = ["foo", "bar"]
        mixer_frame.fs.mix.config.num_inlets = 3

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        with pytest.raises(ConfigurationError):
            mixer_frame.fs.mix.create_inlet_list()

    @pytest.mark.unit
    def test_add_inlet_state_blocks(self, mixer_frame):
        mixer_frame.fs.mix.config.inlet_list = ["foo", "bar"]

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)

        assert isinstance(mixer_frame.fs.mix.foo_state, StateBlock)
        assert isinstance(mixer_frame.fs.mix.bar_state, StateBlock)

        assert len(inlet_blocks) == 2
        for i in inlet_blocks:
            assert isinstance(i, StateBlock)
            assert i.local_name in ["foo_state", "bar_state"]
            assert i[0].config.has_phase_equilibrium is False
            assert i[0].config.defined_state is True
            assert len(i[0].config) == 3

    @pytest.mark.unit
    def test_add_inlet_state_blocks_prop_pack_args(self, mixer_frame):
        mixer_frame.fs.mix.config.property_package_args = {"test": 1}
        mixer_frame.fs.mix.config.inlet_list = ["foo", "bar"]

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)

        assert isinstance(mixer_frame.fs.mix.foo_state, StateBlock)
        assert isinstance(mixer_frame.fs.mix.bar_state, StateBlock)

        assert len(inlet_blocks) == 2
        for i in inlet_blocks:
            assert isinstance(i, StateBlock)
            assert i.local_name in ["foo_state", "bar_state"]
            assert i[0].config.has_phase_equilibrium is False
            assert i[0].config.defined_state is True
            assert len(i[0].config) == 4
            assert i[0].config.test == 1

    @pytest.mark.unit
    def test_add_mixed_state_block(self, mixer_frame):
        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        assert isinstance(mixed_block, StateBlock)
        assert hasattr(mixer_frame.fs.mix, "mixed_state")
        assert not mixer_frame.fs.mix.mixed_state[0].config.has_phase_equilibrium
        assert mixer_frame.fs.mix.mixed_state[0].config.defined_state is False
        assert len(mixer_frame.fs.mix.mixed_state[0].config) == 3

    @pytest.mark.unit
    def test_add_mixed_state_block_prop_pack_args(self, mixer_frame):
        mixer_frame.fs.mix.config.property_package_args = {"test": 1}

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        assert isinstance(mixed_block, StateBlock)
        assert hasattr(mixer_frame.fs.mix, "mixed_state")
        assert not mixer_frame.fs.mix.mixed_state[0].config.has_phase_equilibrium
        assert mixer_frame.fs.mix.mixed_state[0].config.defined_state is False
        assert len(mixer_frame.fs.mix.mixed_state[0].config) == 4
        assert mixer_frame.fs.mix.mixed_state[0].config.test == 1

    @pytest.mark.unit
    def test_get_mixed_state_block(self, mixer_frame):
        mixer_frame.fs.sb = TestStateBlock(
            mixer_frame.fs.time, parameters=mixer_frame.fs.pp
        )

        mixer_frame.fs.mix.config.mixed_state_block = mixer_frame.fs.sb

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        mixed_block = mixer_frame.fs.mix.get_mixed_state_block()

        assert mixed_block == mixer_frame.fs.sb

    @pytest.mark.unit
    def test_get_mixed_state_block_none(self, mixer_frame):
        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        with pytest.raises(BurntToast):
            mixer_frame.fs.mix.get_mixed_state_block()

    @pytest.mark.unit
    def test_get_mixed_state_block_mismatch(self, mixer_frame):
        mixer_frame.fs.sb = TestStateBlock(
            mixer_frame.fs.time, parameters=mixer_frame.fs.pp
        )

        # Change parameters arg to create mismatch
        mixer_frame.fs.sb[0].config.parameters = None

        mixer_frame.fs.mix.config.mixed_state_block = mixer_frame.fs.sb

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        with pytest.raises(ConfigurationError):
            mixer_frame.fs.mix.get_mixed_state_block()

    # Test mixing equation methods
    @pytest.mark.unit
    def test_add_material_mixing_equations_pc(self, mixer_frame):
        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_material_mixing_equations(
            inlet_blocks, mixed_block, MaterialBalanceType.componentPhase
        )

        assert isinstance(mixer_frame.fs.mix.material_mixing_equations, Constraint)
        assert len(mixer_frame.fs.mix.material_mixing_equations) == 4

    @pytest.mark.unit
    def test_add_material_mixing_equations_pc_equilibrium(self, mixer_frame):
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_material_mixing_equations(
            inlet_blocks, mixed_block, MaterialBalanceType.componentPhase
        )

        assert isinstance(mixer_frame.fs.mix.phase_equilibrium_generation, Var)
        assert isinstance(mixer_frame.fs.mix.material_mixing_equations, Constraint)
        assert len(mixer_frame.fs.mix.material_mixing_equations) == 4

    @pytest.mark.unit
    def test_add_material_mixing_equations_pc_equilibrium_not_supported(
        self, mixer_frame
    ):
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        # Remove phase equilibrium list to trigger error
        mixer_frame.fs.pp.del_component(mixer_frame.fs.pp.phase_equilibrium_idx)

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        with pytest.raises(PropertyNotSupportedError):
            mixer_frame.fs.mix.add_material_mixing_equations(
                inlet_blocks, mixed_block, MaterialBalanceType.componentPhase
            )

    @pytest.mark.unit
    def test_add_material_mixing_equations_tc(self, mixer_frame):
        mixer_frame.fs.mix.config.material_balance_type = (
            MaterialBalanceType.componentTotal
        )

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_material_mixing_equations(
            inlet_blocks, mixed_block, MaterialBalanceType.componentTotal
        )

        assert isinstance(mixer_frame.fs.mix.material_mixing_equations, Constraint)
        assert len(mixer_frame.fs.mix.material_mixing_equations) == 2

    @pytest.mark.unit
    def test_add_material_mixing_equations_tc_equilibrium(self, mixer_frame):
        mixer_frame.fs.mix.config.material_balance_type = (
            MaterialBalanceType.componentTotal
        )
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_material_mixing_equations(
            inlet_blocks, mixed_block, MaterialBalanceType.componentTotal
        )

        assert isinstance(mixer_frame.fs.mix.material_mixing_equations, Constraint)
        assert len(mixer_frame.fs.mix.material_mixing_equations) == 2

    @pytest.mark.unit
    def test_add_material_mixing_equations_t(self, mixer_frame):
        mixer_frame.fs.mix.config.material_balance_type = MaterialBalanceType.total

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_material_mixing_equations(
            inlet_blocks, mixed_block, MaterialBalanceType.total
        )

        assert isinstance(mixer_frame.fs.mix.material_mixing_equations, Constraint)
        assert len(mixer_frame.fs.mix.material_mixing_equations) == 1

    @pytest.mark.unit
    def test_add_material_mixing_equations_t_equilibrium(self, mixer_frame):
        mixer_frame.fs.mix.config.material_balance_type = MaterialBalanceType.total
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_material_mixing_equations(
            inlet_blocks, mixed_block, MaterialBalanceType.total
        )

        assert isinstance(mixer_frame.fs.mix.material_mixing_equations, Constraint)
        assert len(mixer_frame.fs.mix.material_mixing_equations) == 1

    @pytest.mark.unit
    def test_add_material_mixing_equations_e(self, mixer_frame):
        mixer_frame.fs.mix.config.material_balance_type = (
            MaterialBalanceType.elementTotal
        )
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        with pytest.raises(ConfigurationError):
            mixer_frame.fs.mix.add_material_mixing_equations(
                inlet_blocks, mixed_block, MaterialBalanceType.elementTotal
            )

    @pytest.mark.unit
    def test_add_material_mixing_equations_none(self, mixer_frame):
        mixer_frame.fs.mix.config.material_balance_type = MaterialBalanceType.none
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_material_mixing_equations(
            inlet_blocks, mixed_block, MaterialBalanceType.none
        )

        assert not hasattr(mixer_frame.fs.mix, "material_mixing_equations")

    @pytest.mark.unit
    def test_add_energy_mixing_equations(self, mixer_frame):
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_energy_mixing_equations(inlet_blocks, mixed_block)

        assert isinstance(mixer_frame.fs.mix.enthalpy_mixing_equations, Constraint)
        assert len(mixer_frame.fs.mix.enthalpy_mixing_equations) == 1

    @pytest.mark.unit
    def test_add_pressure_minimization_equations(self, mixer_frame):
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_pressure_minimization_equations(
            inlet_blocks, mixed_block
        )

        assert mixer_frame.fs.mix.inlet_idx.ctype is RangeSet
        assert isinstance(mixer_frame.fs.mix.minimum_pressure, Var)
        assert len(mixer_frame.fs.mix.minimum_pressure) == 2
        assert isinstance(mixer_frame.fs.mix.eps_pressure, Param)
        assert isinstance(mixer_frame.fs.mix.minimum_pressure_constraint, Constraint)
        assert len(mixer_frame.fs.mix.minimum_pressure) == 2
        assert isinstance(mixer_frame.fs.mix.mixture_pressure, Constraint)

    @pytest.mark.unit
    def test_add_pressure_equality_equations(self, mixer_frame):
        mixer_frame.fs.mix.config.has_phase_equilibrium = True

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_pressure_equality_equations(inlet_blocks, mixed_block)

        assert isinstance(mixer_frame.fs.mix.pressure_equality_constraints, Constraint)
        assert len(mixer_frame.fs.mix.pressure_equality_constraints) == 2

    @pytest.mark.unit
    def test_add_port_objects(self, mixer_frame):
        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_port_objects(inlet_list, inlet_blocks, mixed_block)

        assert isinstance(mixer_frame.fs.mix.inlet_1, Port)
        assert isinstance(mixer_frame.fs.mix.inlet_2, Port)
        assert isinstance(mixer_frame.fs.mix.outlet, Port)

    @pytest.mark.unit
    def test_add_port_objects_construct_ports_False(self, mixer_frame):
        mixer_frame.fs.mix.config.construct_ports = False

        mixer_frame.fs.mix._get_property_package()
        mixer_frame.fs.mix._get_indexing_sets()

        inlet_list = mixer_frame.fs.mix.create_inlet_list()
        inlet_blocks = mixer_frame.fs.mix.add_inlet_state_blocks(inlet_list)
        mixed_block = mixer_frame.fs.mix.add_mixed_state_block()

        mixer_frame.fs.mix.add_port_objects(inlet_list, inlet_blocks, mixed_block)

        assert hasattr(mixer_frame.fs.mix, "inlet_1") is False
        assert hasattr(mixer_frame.fs.mix, "inlet_2") is False
        assert hasattr(mixer_frame.fs.mix, "outlet") is False

    # -------------------------------------------------------------------------
    # Test build method
    @pytest.mark.build
    @pytest.mark.unit
    def test_build_default(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.mix = Mixer(property_package=m.fs.pp)

        assert isinstance(m.fs.mix.material_mixing_equations, Constraint)
        assert len(m.fs.mix.material_mixing_equations) == 4
        assert hasattr(m.fs.mix, "phase_equilibrium_idx_ref") is False

        assert isinstance(m.fs.mix.enthalpy_mixing_equations, Constraint)
        assert len(m.fs.mix.enthalpy_mixing_equations) == 1

        assert m.fs.mix.inlet_idx.ctype is RangeSet
        assert isinstance(m.fs.mix.minimum_pressure, Var)
        assert len(m.fs.mix.minimum_pressure) == 2
        assert isinstance(m.fs.mix.eps_pressure, Param)
        assert isinstance(m.fs.mix.minimum_pressure_constraint, Constraint)
        assert len(m.fs.mix.minimum_pressure) == 2
        assert isinstance(m.fs.mix.mixture_pressure, Constraint)

        assert isinstance(m.fs.mix.inlet_1, Port)
        assert isinstance(m.fs.mix.inlet_2, Port)
        assert isinstance(m.fs.mix.outlet, Port)

    @pytest.mark.build
    @pytest.mark.unit
    def test_build_phase_equilibrium(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.mix = Mixer(property_package=m.fs.pp, has_phase_equilibrium=True)

        assert isinstance(m.fs.mix.material_mixing_equations, Constraint)
        assert len(m.fs.mix.material_mixing_equations) == 4
        assert isinstance(m.fs.mix.phase_equilibrium_generation, Var)

        assert isinstance(m.fs.mix.enthalpy_mixing_equations, Constraint)
        assert len(m.fs.mix.enthalpy_mixing_equations) == 1

        assert m.fs.mix.inlet_idx.ctype is RangeSet
        assert isinstance(m.fs.mix.minimum_pressure, Var)
        assert len(m.fs.mix.minimum_pressure) == 2
        assert isinstance(m.fs.mix.eps_pressure, Param)
        assert isinstance(m.fs.mix.minimum_pressure_constraint, Constraint)
        assert len(m.fs.mix.minimum_pressure) == 2
        assert isinstance(m.fs.mix.mixture_pressure, Constraint)

        assert isinstance(m.fs.mix.inlet_1, Port)
        assert isinstance(m.fs.mix.inlet_2, Port)
        assert isinstance(m.fs.mix.outlet, Port)

    @pytest.mark.build
    @pytest.mark.unit
    def test_build_phase_pressure_equality(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.mix = Mixer(
            property_package=m.fs.pp, momentum_mixing_type=MomentumMixingType.equality
        )

        assert isinstance(m.fs.mix.material_mixing_equations, Constraint)
        assert len(m.fs.mix.material_mixing_equations) == 4

        assert isinstance(m.fs.mix.enthalpy_mixing_equations, Constraint)
        assert len(m.fs.mix.enthalpy_mixing_equations) == 1

        assert isinstance(m.fs.mix.pressure_equality_constraints, Constraint)
        assert len(m.fs.mix.pressure_equality_constraints) == 2

        assert isinstance(m.fs.mix.inlet_1, Port)
        assert isinstance(m.fs.mix.inlet_2, Port)
        assert isinstance(m.fs.mix.outlet, Port)

    # -------------------------------------------------------------------------
    # Test models checks, initialize and release state methods
    @pytest.mark.unit
    def test_model_checks(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.mix = Mixer(
            property_package=m.fs.pp, momentum_mixing_type=MomentumMixingType.equality
        )

        m.fs.mix.model_check()

        assert m.fs.mix.inlet_1_state[0].check is True
        assert m.fs.mix.inlet_2_state[0].check is True
        assert m.fs.mix.mixed_state[0].check is True

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()
        m.fs.sb = TestStateBlock(m.fs.time, parameters=m.fs.pp)

        m.fs.mix = Mixer(property_package=m.fs.pp)

        stable = m.fs.mix._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "component_flow_phase ('p1', 'c1')": getattr(
                        pyunits.pint_registry, "mole/second"
                    ),
                    "component_flow_phase ('p1', 'c2')": getattr(
                        pyunits.pint_registry, "mole/second"
                    ),
                    "component_flow_phase ('p2', 'c1')": getattr(
                        pyunits.pint_registry, "mole/second"
                    ),
                    "component_flow_phase ('p2', 'c2')": getattr(
                        pyunits.pint_registry, "mole/second"
                    ),
                    "temperature": getattr(pyunits.pint_registry, "K"),
                    "pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "inlet_1": {
                    "component_flow_phase ('p1', 'c1')": 2.00,
                    "component_flow_phase ('p1', 'c2')": 2.00,
                    "component_flow_phase ('p2', 'c1')": 2.00,
                    "component_flow_phase ('p2', 'c2')": 2.00,
                    "temperature": 300,
                    "pressure": 1e5,
                },
                "inlet_2": {
                    "component_flow_phase ('p1', 'c1')": 2.00,
                    "component_flow_phase ('p1', 'c2')": 2.00,
                    "component_flow_phase ('p2', 'c1')": 2.00,
                    "component_flow_phase ('p2', 'c2')": 2.00,
                    "temperature": 300,
                    "pressure": 1e5,
                },
                "Outlet": {
                    "component_flow_phase ('p1', 'c1')": 2.00,
                    "component_flow_phase ('p1', 'c2')": 2.00,
                    "component_flow_phase ('p2', 'c1')": 2.00,
                    "component_flow_phase ('p2', 'c2')": 2.00,
                    "temperature": 300,
                    "pressure": 1e5,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.initialization
    @pytest.mark.component
    def test_initialize(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()
        m.fs.sb = TestStateBlock(m.fs.time, parameters=m.fs.pp)

        m.fs.mix = Mixer(property_package=m.fs.pp, mixed_state_block=m.fs.sb)

        # Change one inlet pressure to check initialization calculations
        m.fs.mix.inlet_1_state[0].pressure = 8e4

        f = m.fs.mix.initialize(hold_state=True)

        assert m.fs.mix.inlet_1_state[0].init_test is True
        assert m.fs.mix.inlet_2_state[0].init_test is True
        assert m.fs.sb[0].init_test is True
        assert m.fs.mix.inlet_1_state[0].hold_state is True
        assert m.fs.mix.inlet_2_state[0].hold_state is True
        assert m.fs.sb[0].hold_state is False

        assert m.fs.sb[0].flow_mol_phase_comp["p1", "c1"].value == 4
        assert m.fs.sb[0].flow_mol_phase_comp["p1", "c2"].value == 4
        assert m.fs.sb[0].flow_mol_phase_comp["p2", "c1"].value == 4
        assert m.fs.sb[0].flow_mol_phase_comp["p2", "c2"].value == 4

        assert m.fs.sb[0].temperature.value == 300

        assert m.fs.sb[0].pressure.value == 8e4

        m.fs.mix.release_state(flags=f)

        assert m.fs.mix.inlet_1_state[0].hold_state is False
        assert m.fs.mix.inlet_2_state[0].hold_state is False
        assert m.fs.sb[0].hold_state is False


# -----------------------------------------------------------------------------
class TestBTX(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(valid_phase="Liq")

        m.fs.unit = Mixer(property_package=m.fs.properties)

        m.fs.unit.inlet_1.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.inlet_1.temperature[0].fix(365)  # K
        m.fs.unit.inlet_1.pressure[0].fix(2e5)  # Pa
        m.fs.unit.inlet_1.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet_1.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.inlet_2.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.inlet_2.temperature[0].fix(300)  # K
        m.fs.unit.inlet_2.pressure[0].fix(101325)  # Pa
        m.fs.unit.inlet_2.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet_2.mole_frac_comp[0, "toluene"].fix(0.5)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "inlet_1")
        assert len(btx.fs.unit.inlet_1.vars) == 4
        assert hasattr(btx.fs.unit.inlet_1, "flow_mol")
        assert hasattr(btx.fs.unit.inlet_1, "mole_frac_comp")
        assert hasattr(btx.fs.unit.inlet_1, "temperature")
        assert hasattr(btx.fs.unit.inlet_1, "pressure")

        assert hasattr(btx.fs.unit, "inlet_2")
        assert len(btx.fs.unit.inlet_2.vars) == 4
        assert hasattr(btx.fs.unit.inlet_2, "flow_mol")
        assert hasattr(btx.fs.unit.inlet_2, "mole_frac_comp")
        assert hasattr(btx.fs.unit.inlet_2, "temperature")
        assert hasattr(btx.fs.unit.inlet_2, "pressure")

        assert hasattr(btx.fs.unit, "outlet")
        assert len(btx.fs.unit.outlet.vars) == 4
        assert hasattr(btx.fs.unit.outlet, "flow_mol")
        assert hasattr(btx.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.outlet, "temperature")
        assert hasattr(btx.fs.unit.outlet, "pressure")

        assert number_variables(btx) == 35
        assert number_total_constraints(btx) == 25
        assert number_unused_variables(btx) == 0

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "flow_mol": getattr(pyunits.pint_registry, "mole/second"),
                    "mole_frac_comp benzene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "mole_frac_comp toluene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "temperature": getattr(pyunits.pint_registry, "kelvin"),
                    "pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "inlet_1": {
                    "flow_mol": 5.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 365,
                    "pressure": 2e5,
                },
                "inlet_2": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 300,
                    "pressure": 101325.0,
                },
                "Outlet": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 298.15,
                    "pressure": 101325.0,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert pytest.approx(6, abs=1e-3) == value(btx.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(354.7, abs=1e-1) == value(
            btx.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e3) == value(btx.fs.unit.outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.inlet_1.flow_mol[0]
                    + btx.fs.unit.inlet_2.flow_mol[0]
                    - btx.fs.unit.outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        assert 1e-6 >= abs(
            value(
                btx.fs.unit.inlet_1.flow_mol[0]
                * btx.fs.unit.inlet_1_state[0].enth_mol_phase["Liq"]
                + btx.fs.unit.inlet_2.flow_mol[0]
                * btx.fs.unit.inlet_2_state[0].enth_mol_phase["Liq"]
                - btx.fs.unit.outlet.flow_mol[0]
                * btx.fs.unit.mixed_state[0].enth_mol_phase["Liq"]
            )
        )


# -----------------------------------------------------------------------------
# Tests for Mixer in cases where proeprties do not support pressure
@declare_process_block_class("NoPressureTestBlock")
class _NoPressureParameterBlock(PhysicalParameterBlock):
    def build(self):
        super(_NoPressureParameterBlock, self).build()

        self.p1 = Phase()
        self.p2 = Phase()
        self.c1 = Component()
        self.c2 = Component()

        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])

        self.phase_equilibrium_list = {
            "e1": ["c1", ("p1", "p2")],
            "e2": ["c2", ("p1", "p2")],
        }

        self._state_block_class = NoPressureStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.g,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


@declare_process_block_class("NoPressureStateBlock", block_class=StateBlock)
class NoPressureStateBlockData(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(NoPressureStateBlockData, self).build()

        self.flow_vol = Var(initialize=20)
        self.flow_mol_phase_comp = Var(
            self._params.phase_list, self._params.component_list, initialize=2
        )
        self.temperature = Var(initialize=300)
        self.test_var = Var(initialize=1)

    def get_material_flow_terms(b, p, j):
        return b.test_var

    def get_enthalpy_flow_terms(b, p):
        return b.test_var

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal


class TestMixer_NoPressure(object):
    @declare_process_block_class("MixerFrame")
    class MixerFrameData(MixerData):
        def build(self):
            super(MixerData, self).build()

    @pytest.mark.build
    @pytest.mark.unit
    def test_pressure_minimization_unsupported(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = NoPressureTestBlock()

        with pytest.raises(
            PropertyNotSupportedError,
            match="fs.mix The property package supplied for this unit "
            "does not appear to support pressure, which is required for "
            "momentum mixing. Please set momentum_mixing_type to "
            "MomentumMixingType.none or provide a property package which "
            "supports pressure.",
        ):
            m.fs.mix = Mixer(
                property_package=m.fs.pp,
                momentum_mixing_type=MomentumMixingType.minimize,
                construct_ports=False,
            )

    @pytest.mark.build
    @pytest.mark.unit
    def test_pressure_equal_unsupported(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = NoPressureTestBlock()

        with pytest.raises(
            PropertyNotSupportedError,
            match="fs.mix The property package supplied for this unit "
            "does not appear to support pressure, which is required for "
            "momentum mixing. Please set momentum_mixing_type to "
            "MomentumMixingType.none or provide a property package which "
            "supports pressure.",
        ):
            m.fs.mix = Mixer(
                property_package=m.fs.pp,
                momentum_mixing_type=MomentumMixingType.equality,
                construct_ports=False,
            )

    @pytest.mark.build
    @pytest.mark.unit
    def test_pressure_equal_and_min_unsupported(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = NoPressureTestBlock()

        with pytest.raises(
            PropertyNotSupportedError,
            match="fs.mix The property package supplied for this unit "
            "does not appear to support pressure, which is required for "
            "momentum mixing. Please set momentum_mixing_type to "
            "MomentumMixingType.none or provide a property package which "
            "supports pressure.",
        ):
            m.fs.mix = Mixer(
                property_package=m.fs.pp,
                momentum_mixing_type=MomentumMixingType.minimize_and_equality,
                construct_ports=False,
            )

    @pytest.mark.build
    @pytest.mark.unit
    def test_pressure_none_unsupported(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = NoPressureTestBlock()

        m.fs.mix = Mixer(
            property_package=m.fs.pp,
            momentum_mixing_type=MomentumMixingType.none,
            construct_ports=False,
        )


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = Mixer(
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentTotal,
            momentum_mixing_type=MomentumMixingType.equality,
        )

        m.fs.unit.inlet_1.flow_mol[0].fix(100)
        m.fs.unit.inlet_1.enth_mol[0].fix(5500)
        m.fs.unit.inlet_1.pressure[0].fix(101325)

        m.fs.unit.inlet_2.flow_mol[0].fix(100)
        m.fs.unit.inlet_2.enth_mol[0].fix(5000)
        m.fs.unit.inlet_2.pressure[0].value = 1e5

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iapws):
        assert len(iapws.fs.unit.inlet_1.vars) == 3
        assert hasattr(iapws.fs.unit.inlet_1, "flow_mol")
        assert hasattr(iapws.fs.unit.inlet_1, "enth_mol")
        assert hasattr(iapws.fs.unit.inlet_1, "pressure")

        assert len(iapws.fs.unit.inlet_2.vars) == 3
        assert hasattr(iapws.fs.unit.inlet_2, "flow_mol")
        assert hasattr(iapws.fs.unit.inlet_2, "enth_mol")
        assert hasattr(iapws.fs.unit.inlet_2, "pressure")

        assert hasattr(iapws.fs.unit, "outlet")
        assert len(iapws.fs.unit.outlet.vars) == 3
        assert hasattr(iapws.fs.unit.outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.outlet, "pressure")

    @pytest.mark.component
    def test_units(self, iapws):
        assert_units_consistent(iapws)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, iapws):
        stable = iapws.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "Molar Flow": getattr(pyunits.pint_registry, "mole/second"),
                    "Mass Flow": getattr(pyunits.pint_registry, "kg/second"),
                    "T": getattr(pyunits.pint_registry, "K"),
                    "P": getattr(pyunits.pint_registry, "Pa"),
                    "Vapor Fraction": getattr(pyunits.pint_registry, "dimensionless"),
                    "Molar Enthalpy": getattr(pyunits.pint_registry, "J/mole"),
                },
                "inlet_1": {
                    "Molar Flow": 100,
                    "Mass Flow": 1.8015,
                    "T": 346.05,
                    "P": 101325,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 5500.0,
                },
                "inlet_2": {
                    "Molar Flow": 100,
                    "Mass Flow": 1.8015,
                    "T": 339.43,
                    "P": 1e5,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 5000,
                },
                "Outlet": {
                    "Molar Flow": 1,
                    "Mass Flow": 1.8015e-2,
                    "T": 270.4877112932641,
                    "P": 11032305.8275,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 0.01102138712926277,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(200, abs=1e-5) == value(iapws.fs.unit.outlet.flow_mol[0])

        assert pytest.approx(5250, abs=1e0) == value(iapws.fs.unit.outlet.enth_mol[0])

        assert pytest.approx(101325, abs=1e2) == value(iapws.fs.unit.outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert (
            abs(
                value(
                    iapws.fs.unit.inlet_1.flow_mol[0]
                    + iapws.fs.unit.inlet_2.flow_mol[0]
                    - iapws.fs.unit.outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    iapws.fs.unit.inlet_1.flow_mol[0]
                    * iapws.fs.unit.inlet_1.enth_mol[0]
                    + iapws.fs.unit.inlet_2.flow_mol[0]
                    * iapws.fs.unit.inlet_2.enth_mol[0]
                    - iapws.fs.unit.outlet.flow_mol[0]
                    * iapws.fs.unit.outlet.enth_mol[0]
                )
            )
            <= 1e-6
        )


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = Mixer(property_package=m.fs.properties)

        m.fs.unit.inlet_1.flow_vol[0].fix(1e-3)
        m.fs.unit.inlet_1.temperature[0].fix(320)
        m.fs.unit.inlet_1.pressure[0].fix(101325)
        m.fs.unit.inlet_1.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet_1.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet_1.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet_1.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet_1.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.inlet_2.flow_vol[0].fix(1e-3)
        m.fs.unit.inlet_2.temperature[0].fix(300)
        m.fs.unit.inlet_2.pressure[0].fix(101325)
        m.fs.unit.inlet_2.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet_2.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet_2.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet_2.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet_2.conc_mol_comp[0, "Ethanol"].fix(0.0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):
        assert len(sapon.fs.unit.inlet_1.vars) == 4
        assert hasattr(sapon.fs.unit.inlet_1, "flow_vol")
        assert hasattr(sapon.fs.unit.inlet_1, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.inlet_1, "temperature")
        assert hasattr(sapon.fs.unit.inlet_1, "pressure")

        assert len(sapon.fs.unit.inlet_2.vars) == 4
        assert hasattr(sapon.fs.unit.inlet_2, "flow_vol")
        assert hasattr(sapon.fs.unit.inlet_2, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.inlet_2, "temperature")
        assert hasattr(sapon.fs.unit.inlet_2, "pressure")

        assert len(sapon.fs.unit.outlet.vars) == 4
        assert hasattr(sapon.fs.unit.outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.outlet, "temperature")
        assert hasattr(sapon.fs.unit.outlet, "pressure")

        assert number_variables(sapon) == 26
        assert number_total_constraints(sapon) == 10
        assert number_unused_variables(sapon) == 0

    @pytest.mark.component
    def test_units(self, sapon):
        assert_units_consistent(sapon)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, sapon):
        stable = sapon.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "Volumetric Flowrate": getattr(
                        pyunits.pint_registry, "m**3/second"
                    ),
                    "Molar Concentration H2O": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Molar Concentration NaOH": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Molar Concentration EthylAcetate": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Molar Concentration SodiumAcetate": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Molar Concentration Ethanol": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Temperature": getattr(pyunits.pint_registry, "K"),
                    "Pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "inlet_1": {
                    "Volumetric Flowrate": 1e-3,
                    "Molar Concentration H2O": 55388,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 0,
                    "Molar Concentration Ethanol": 0,
                    "Temperature": 320,
                    "Pressure": 1.0132e05,
                },
                "inlet_2": {
                    "Volumetric Flowrate": 1e-3,
                    "Molar Concentration H2O": 55388,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 0,
                    "Molar Concentration Ethanol": 0,
                    "Temperature": 300,
                    "Pressure": 1.0132e05,
                },
                "Outlet": {
                    "Volumetric Flowrate": 1.00,
                    "Molar Concentration H2O": 100.00,
                    "Molar Concentration NaOH": 100.0,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 100.00,
                    "Molar Concentration Ethanol": 100.00,
                    "Temperature": 298.15,
                    "Pressure": 1.0132e05,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.component
    def test_scaling(self, sapon):
        sapon.fs.unit.calculate_scaling_factors()

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(2e-3, abs=1e-6) == value(sapon.fs.unit.outlet.flow_vol[0])

        assert pytest.approx(55388.0, abs=1e0) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, abs=1e-3) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, abs=1e-3) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]
        )

        assert pytest.approx(310.0, abs=1e-1) == value(
            sapon.fs.unit.outlet.temperature[0]
        )

        assert pytest.approx(101325, abs=1e2) == value(sapon.fs.unit.outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet_1.flow_vol[0]
                    + sapon.fs.unit.inlet_2.flow_vol[0]
                    - sapon.fs.unit.outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    sapon.fs.unit.inlet_1.flow_vol[0]
                    * sapon.fs.properties.dens_mol
                    * sapon.fs.properties.cp_mol
                    * (
                        sapon.fs.unit.inlet_1.temperature[0]
                        - sapon.fs.properties.temperature_ref
                    )
                    + sapon.fs.unit.inlet_2.flow_vol[0]
                    * sapon.fs.properties.dens_mol
                    * sapon.fs.properties.cp_mol
                    * (
                        sapon.fs.unit.inlet_2.temperature[0]
                        - sapon.fs.properties.temperature_ref
                    )
                    - sapon.fs.unit.outlet.flow_vol[0]
                    * sapon.fs.properties.dens_mol
                    * sapon.fs.properties.cp_mol
                    * (
                        sapon.fs.unit.outlet.temperature[0]
                        - sapon.fs.properties.temperature_ref
                    )
                )
            )
            <= 1e-3
        )


@pytest.mark.component
def test_construction_component_not_in_phase():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()
    m.fs.prop_params = GenericParameterBlock(**get_prop(["H2O", "H2"], ["Liq", "Vap"]))
    m.fs.inject1 = Mixer(
        property_package=m.fs.prop_params,
        inlet_list=["in1", "in2"],
        momentum_mixing_type=MomentumMixingType.none,
    )
    iscale.calculate_scaling_factors(m)


@pytest.mark.unit
def test_initialization_error():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = Mixer(property_package=m.fs.pp)

    m.fs.mix.inlet_1_state[0].material_flow_mol.fix(10)
    m.fs.mix.inlet_2_state[0].material_flow_mol.fix(10)
    m.fs.mix.mixed_state[0].material_flow_mol.fix(100)

    with pytest.raises(InitializationError):
        m.fs.mix.initialize()
