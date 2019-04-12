##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for ControlVolumeBlockData.

Author: Andrew Lee
"""
import pytest

from pyomo.environ import ConcreteModel, Constraint, Param, Set, Var
from pyomo.network import Port
from pyomo.common.config import ConfigBlock

from idaes.core import (FlowsheetBlockData,
                        declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlock,
                        StateBlockData)
from idaes.unit_models.mixer import (Mixer,
                                     MixerData,
                                     MixingType,
                                     MomentumMixingType)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError,
                                        PropertyNotSupportedError)


# -----------------------------------------------------------------------------
# Mockup classes for testing
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("PhysicalParameterTestBlock")
class _PhysicalParameterBlock(PhysicalParameterBlock):
    def build(self):
        super(_PhysicalParameterBlock, self).build()

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])
        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])

        self.state_block_class = TestStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class SBlockBase(StateBlock):
    def initialize(blk, outlvl=0, optarg=None, solver=None,
                   hold_state=False, **state_args):
        for k in blk.keys():
            blk[k].init_test = True
            blk[k].hold_state = hold_state

    def release_state(blk, flags=None, outlvl=0):
        for k in blk.keys():
            blk[k].hold_state = not blk[k].hold_state


@declare_process_block_class("TestStateBlock", block_class=SBlockBase)
class StateTestBlockData(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(StateTestBlockData, self).build()

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])
        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])
        self.phase_equilibrium_list = \
            {"e1": ["c1", ("p1", "p2")],
             "e2": ["c2", ("p1", "p2")]}

        self.pressure = Var(initialize=1e5)
        self.flow_mol_phase_comp = Var(self.phase_list,
                                       self.component_list,
                                       initialize=1)
        self.enth_mol_phase = Var(self.phase_list, initialize=2)

    def get_material_flow_terms(b, p, j):
        return b.flow_mol_phase_comp[p, j]

    def get_enthalpy_flow_terms(b, p):
        return b.enth_mol_phase[p]

    def define_state_vars(self):
        return {"component_flow": self.flow_mol_phase_comp,
                "enthalpy": self.enth_mol_phase,
                "pressure": self.pressure}

    def model_check(self):
        self.check = True


@declare_process_block_class("MixerFrame")
class MixerFrameData(MixerData):
    def build(self):
        super(MixerData, self).build()


# -----------------------------------------------------------------------------
# Test Mixer unit model
def test_mixer_config():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp})

    assert len(m.fs.mix.config) == 12
    assert m.fs.mix.config.dynamic is False
    assert m.fs.mix.config.has_holdup is False
    assert m.fs.mix.config.property_package == m.fs.pp
    assert isinstance(m.fs.mix.config.property_package_args, ConfigBlock)
    assert len(m.fs.mix.config.property_package_args) == 0
    assert m.fs.mix.config.inlet_list is None
    assert m.fs.mix.config.num_inlets is None
    assert m.fs.mix.config.has_phase_equilibrium is False
    assert m.fs.mix.config.material_mixing_type == MixingType.extensive
    assert m.fs.mix.config.energy_mixing_type == MixingType.extensive
    assert m.fs.mix.config.momentum_mixing_type == MomentumMixingType.minimize
    assert m.fs.mix.config.mixed_state_block is None
    assert m.fs.mix.config.construct_ports is True


def test_inherited_methods():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    assert hasattr(m.fs.mix, "phase_list_ref")
    assert hasattr(m.fs.mix, "component_list_ref")


def test_create_inlet_list_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()

    for i in inlet_list:
        assert i in ["inlet_1", "inlet_2"]


def test_create_inlet_list_inlet_list():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "inlet_list": ["foo", "bar"]})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()

    for i in inlet_list:
        assert i in ["foo", "bar"]


def test_create_inlet_list_num_inlets():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "num_inlets": 3})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()

    for i in inlet_list:
        assert i in ["inlet_1", "inlet_2", "inlet_3"]


def test_create_inlet_list_both_args_consistent():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "inlet_list": ["foo", "bar"],
                                   "num_inlets": 2})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()

    for i in inlet_list:
        assert i in ["foo", "bar"]


def test_create_inlet_list_both_args_inconsistent():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "inlet_list": ["foo", "bar"],
                                   "num_inlets": 3})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    with pytest.raises(ConfigurationError):
        m.fs.mix.create_inlet_list()


def test_add_inlet_state_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "inlet_list": ["foo", "bar"]})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)

    assert isinstance(m.fs.mix.foo_state, StateBlock)
    assert isinstance(m.fs.mix.bar_state, StateBlock)

    assert len(inlet_blocks) == 2
    for i in inlet_blocks:
        assert isinstance(i, StateBlock)
        assert i.local_name in ["foo_state", "bar_state"]
        assert i[0].config.has_phase_equilibrium is False
        assert i[0].config.defined_state is True
        assert len(i[0].config) == 3


def test_add_inlet_state_blocks_prop_pack_args():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "property_package_args": {"test": 1},
                                   "inlet_list": ["foo", "bar"]})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)

    assert isinstance(m.fs.mix.foo_state, StateBlock)
    assert isinstance(m.fs.mix.bar_state, StateBlock)

    assert len(inlet_blocks) == 2
    for i in inlet_blocks:
        assert isinstance(i, StateBlock)
        assert i.local_name in ["foo_state", "bar_state"]
        assert i[0].config.has_phase_equilibrium is False
        assert i[0].config.defined_state is True
        assert len(i[0].config) == 4
        assert i[0].config.test == 1


def test_add_mixed_state_block():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    mixed_block = m.fs.mix.add_mixed_state_block()

    assert isinstance(mixed_block, StateBlock)
    assert hasattr(m.fs.mix, "mixed_state")
    assert m.fs.mix.mixed_state[0].config.has_phase_equilibrium is False
    assert m.fs.mix.mixed_state[0].config.defined_state is False
    assert len(m.fs.mix.mixed_state[0].config) == 3


def test_add_mixed_state_block_prop_pack_args():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "property_package_args": {"test": 1}})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    mixed_block = m.fs.mix.add_mixed_state_block()

    assert isinstance(mixed_block, StateBlock)
    assert hasattr(m.fs.mix, "mixed_state")
    assert m.fs.mix.mixed_state[0].config.has_phase_equilibrium is False
    assert m.fs.mix.mixed_state[0].config.defined_state is False
    assert len(m.fs.mix.mixed_state[0].config) == 4
    assert m.fs.mix.mixed_state[0].config.test == 1


def test_get_mixed_state_block():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    mixed_block = m.fs.mix.get_mixed_state_block()

    assert mixed_block == m.fs.sb


def test_get_mixed_state_block_none():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    with pytest.raises(BurntToast):
        m.fs.mix.get_mixed_state_block()


def test_get_mixed_state_block_mismatch():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    # Change parameters arg to create mismatch
    m.fs.sb[0].config.parameters = None

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    with pytest.raises(ConfigurationError):
        m.fs.mix.get_mixed_state_block()


# Test mixing equation methods
def test_add_material_mixing_equations():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)
    mixed_block = m.fs.mix.get_mixed_state_block()

    m.fs.mix.add_material_mixing_equations(inlet_blocks, mixed_block)

    assert isinstance(m.fs.mix.material_mixing_equations, Constraint)
    assert len(m.fs.mix.material_mixing_equations) == 4


def test_add_material_mixing_equations_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb,
                                   "has_phase_equilibrium": True})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)
    mixed_block = m.fs.mix.get_mixed_state_block()

    m.fs.mix.add_material_mixing_equations(inlet_blocks, mixed_block)

    assert hasattr(m.fs.mix, "phase_equilibrium_idx_ref")
    assert isinstance(m.fs.mix.phase_equilibrium_generation, Var)
    assert isinstance(m.fs.mix.material_mixing_equations, Constraint)
    assert len(m.fs.mix.material_mixing_equations) == 4


def test_add_material_mixing_equations_equilibrium_not_supported():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb,
                                   "has_phase_equilibrium": True})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)
    mixed_block = m.fs.mix.get_mixed_state_block()

    with pytest.raises(PropertyNotSupportedError):
        m.fs.mix.add_material_mixing_equations(inlet_blocks, mixed_block)


def test_add_energy_mixing_equations():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb,
                                   "has_phase_equilibrium": True})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)
    mixed_block = m.fs.mix.get_mixed_state_block()

    m.fs.mix.add_energy_mixing_equations(inlet_blocks, mixed_block)

    assert isinstance(m.fs.mix.enthalpy_mixing_equations, Constraint)
    assert len(m.fs.mix.enthalpy_mixing_equations) == 1


def test_add_pressure_minimization_equations():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb,
                                   "has_phase_equilibrium": True})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)
    mixed_block = m.fs.mix.get_mixed_state_block()

    m.fs.mix.add_pressure_minimization_equations(inlet_blocks, mixed_block)

    assert isinstance(m.fs.mix.inlet_idx, Set)
    assert isinstance(m.fs.mix.minimum_pressure, Var)
    assert len(m.fs.mix.minimum_pressure) == 2
    assert isinstance(m.fs.mix.eps_pressure, Param)
    assert isinstance(m.fs.mix.minimum_pressure_constraint, Constraint)
    assert len(m.fs.mix.minimum_pressure) == 2
    assert isinstance(m.fs.mix.mixture_pressure, Constraint)


def test_add_pressure_equality_equations():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb,
                                   "has_phase_equilibrium": True})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)
    mixed_block = m.fs.mix.get_mixed_state_block()

    m.fs.mix.add_pressure_equality_equations(inlet_blocks, mixed_block)

    assert isinstance(m.fs.mix.pressure_equality_constraints, Constraint)
    assert len(m.fs.mix.pressure_equality_constraints) == 2


def test_add_port_objects():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)
    mixed_block = m.fs.mix.add_mixed_state_block()

    m.fs.mix.add_port_objects(inlet_list, inlet_blocks, mixed_block)

    assert isinstance(m.fs.mix.inlet_1, Port)
    assert isinstance(m.fs.mix.inlet_2, Port)
    assert isinstance(m.fs.mix.outlet, Port)


def test_add_port_objects_construct_ports_False():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = MixerFrame(default={"property_package": m.fs.pp,
                                   "construct_ports": False})

    m.fs.mix._get_property_package()
    m.fs.mix._get_indexing_sets()

    inlet_list = m.fs.mix.create_inlet_list()
    inlet_blocks = m.fs.mix.add_inlet_state_blocks(inlet_list)
    mixed_block = m.fs.mix.add_mixed_state_block()

    m.fs.mix.add_port_objects(inlet_list, inlet_blocks, mixed_block)

    assert hasattr(m.fs.mix, "inlet_1") is False
    assert hasattr(m.fs.mix, "inlet_2") is False
    assert hasattr(m.fs.mix, "outlet") is False


# -----------------------------------------------------------------------------
# Test build method
def test_build_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = Mixer(default={"property_package": m.fs.pp})

    assert isinstance(m.fs.mix.material_mixing_equations, Constraint)
    assert len(m.fs.mix.material_mixing_equations) == 4
    assert hasattr(m.fs.mix, "phase_equilibrium_idx_ref") is False

    assert isinstance(m.fs.mix.enthalpy_mixing_equations, Constraint)
    assert len(m.fs.mix.enthalpy_mixing_equations) == 1

    assert isinstance(m.fs.mix.inlet_idx, Set)
    assert isinstance(m.fs.mix.minimum_pressure, Var)
    assert len(m.fs.mix.minimum_pressure) == 2
    assert isinstance(m.fs.mix.eps_pressure, Param)
    assert isinstance(m.fs.mix.minimum_pressure_constraint, Constraint)
    assert len(m.fs.mix.minimum_pressure) == 2
    assert isinstance(m.fs.mix.mixture_pressure, Constraint)

    assert isinstance(m.fs.mix.inlet_1, Port)
    assert isinstance(m.fs.mix.inlet_2, Port)
    assert isinstance(m.fs.mix.outlet, Port)


def test_build_phase_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = Mixer(default={"property_package": m.fs.pp,
                              "has_phase_equilibrium": True})

    assert isinstance(m.fs.mix.material_mixing_equations, Constraint)
    assert len(m.fs.mix.material_mixing_equations) == 4
    assert hasattr(m.fs.mix, "phase_equilibrium_idx_ref")
    assert isinstance(m.fs.mix.phase_equilibrium_generation, Var)

    assert isinstance(m.fs.mix.enthalpy_mixing_equations, Constraint)
    assert len(m.fs.mix.enthalpy_mixing_equations) == 1

    assert isinstance(m.fs.mix.inlet_idx, Set)
    assert isinstance(m.fs.mix.minimum_pressure, Var)
    assert len(m.fs.mix.minimum_pressure) == 2
    assert isinstance(m.fs.mix.eps_pressure, Param)
    assert isinstance(m.fs.mix.minimum_pressure_constraint, Constraint)
    assert len(m.fs.mix.minimum_pressure) == 2
    assert isinstance(m.fs.mix.mixture_pressure, Constraint)

    assert isinstance(m.fs.mix.inlet_1, Port)
    assert isinstance(m.fs.mix.inlet_2, Port)
    assert isinstance(m.fs.mix.outlet, Port)


def test_build_phase_pressure_equality():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = Mixer(default={
            "property_package": m.fs.pp,
            "momentum_mixing_type": MomentumMixingType.equality})

    assert isinstance(m.fs.mix.material_mixing_equations, Constraint)
    assert len(m.fs.mix.material_mixing_equations) == 4

    assert isinstance(m.fs.mix.enthalpy_mixing_equations, Constraint)
    assert len(m.fs.mix.enthalpy_mixing_equations) == 1

    assert isinstance(m.fs.mix.pressure_equality_constraints, Constraint)
    assert len(m.fs.mix.pressure_equality_constraints) == 2

    assert isinstance(m.fs.mix.inlet_1, Port)
    assert isinstance(m.fs.mix.inlet_2, Port)
    assert isinstance(m.fs.mix.outlet, Port)


# -----------------------------------------------------------------------------
# Test models checks, initialize and release state methods
def test_model_checks():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.mix = Mixer(default={
            "property_package": m.fs.pp,
            "momentum_mixing_type": MomentumMixingType.equality})

    m.fs.mix.model_check()

    assert m.fs.mix.inlet_1_state[0].check is True
    assert m.fs.mix.inlet_2_state[0].check is True
    assert m.fs.mix.mixed_state[0].check is True


def test_initialize():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.mix = Mixer(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb})

    # Change one inlet pressure to check initialization calculations
    m.fs.mix.inlet_1_state[0].pressure = 8e4

    f = m.fs.mix.initialize()

    assert m.fs.mix.inlet_1_state[0].init_test is True
    assert m.fs.mix.inlet_2_state[0].init_test is True
    assert m.fs.sb[0].init_test is True
    assert m.fs.mix.inlet_1_state[0].hold_state is True
    assert m.fs.mix.inlet_2_state[0].hold_state is True
    assert m.fs.sb[0].hold_state is False

    assert m.fs.sb[0].flow_mol_phase_comp["p1", "c1"].value == 2
    assert m.fs.sb[0].flow_mol_phase_comp["p1", "c2"].value == 2
    assert m.fs.sb[0].flow_mol_phase_comp["p2", "c1"].value == 2
    assert m.fs.sb[0].flow_mol_phase_comp["p2", "c2"].value == 2

    assert m.fs.sb[0].enth_mol_phase["p1"].value == 2
    assert m.fs.sb[0].enth_mol_phase["p2"].value == 2

    assert m.fs.sb[0].pressure.value == 8e4

    m.fs.mix.release_state(flags=f)

    assert m.fs.mix.inlet_1_state[0].hold_state is False
    assert m.fs.mix.inlet_2_state[0].hold_state is False
    assert m.fs.sb[0].hold_state is False
