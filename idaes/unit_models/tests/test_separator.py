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

from pyomo.environ import (ConcreteModel, Constraint, Set, SolverFactory, Var,
                           value)
from pyomo.network import Port
from pyomo.common.config import ConfigBlock

from idaes.core import (FlowsheetBlock,
                        FlowsheetBlockData,
                        declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlock,
                        StateBlockData,
                        MaterialBalanceType)
from idaes.unit_models.separator import (Separator,
                                         SeparatorData,
                                         SplittingType,
                                         EnergySplittingType)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.property_models.examples.saponification_thermo import (
    SaponificationParameterBlock)

from idaes.property_models import iapws95


# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


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
        self.enth_mol = Var(initialize=2) # total molar enthalpy (both phases)
        self.temperature = Var(initialize=5)

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


@declare_process_block_class("SeparatorFrame")
class SeparatorFrameData(SeparatorData):
    def build(self):
        super(SeparatorData, self).build()


# -----------------------------------------------------------------------------
# Basic tests of Separator unit model
def test_separator_config():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

    assert len(m.fs.sep.config) == 14
    assert m.fs.sep.config.dynamic is False
    assert m.fs.sep.config.has_holdup is False
    assert m.fs.sep.config.property_package == m.fs.pp
    assert isinstance(m.fs.sep.config.property_package_args, ConfigBlock)
    assert len(m.fs.sep.config.property_package_args) == 0
    assert m.fs.sep.config.outlet_list is None
    assert m.fs.sep.config.num_outlets is None
    assert m.fs.sep.config.split_basis == SplittingType.totalFlow
    assert m.fs.sep.config.ideal_separation is True
    assert m.fs.sep.config.ideal_split_map is None
    assert m.fs.sep.config.mixed_state_block is None
    assert m.fs.sep.config.construct_ports is True
    assert m.fs.sep.config.material_balance_type == \
        MaterialBalanceType.componentPhase
    assert m.fs.sep.config.has_phase_equilibrium is False


def test_validate_config_arguments():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "has_phase_equilibrium": True,
                                       "ideal_separation": True})

    with pytest.raises(ConfigurationError):
        m.fs.sep._validate_config_arguments()


def test_inherited_methods():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    assert hasattr(m.fs.sep.config.property_package, "phase_list")

def test_create_outlet_list_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()

    for o in outlet_list:
        assert o in ["outlet_1", "outlet_2"]


def test_create_outlet_list_outlet_list():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "outlet_list": ["foo", "bar"]})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()

    for o in outlet_list:
        assert o in ["foo", "bar"]


def test_create_outlet_list_num_outlets():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "num_outlets": 3})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()

    for o in outlet_list:
        assert o in ["outlet_1", "outlet_2", "outlet_3"]


def test_create_outlet_list_both_args_consistent():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "outlet_list": ["foo", "bar"],
                                       "num_outlets": 2})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()

    for o in outlet_list:
        assert o in ["foo", "bar"]


def test_create_outlet_list_both_args_inconsistent():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "outlet_list": ["foo", "bar"],
                                       "num_outlets": 3})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    with pytest.raises(ConfigurationError):
        m.fs.sep.create_outlet_list()


def test_add_outlet_state_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "outlet_list": ["foo", "bar"]})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)

    assert isinstance(m.fs.sep.foo_state, StateBlock)
    assert isinstance(m.fs.sep.bar_state, StateBlock)

    assert len(outlet_blocks) == 2
    for o in outlet_blocks:
        assert isinstance(o, StateBlock)
        assert o.local_name in ["foo_state", "bar_state"]
        assert o[0].config.has_phase_equilibrium is False
        assert o[0].config.defined_state is False
        assert len(o[0].config) == 3


def test_add_outlet_state_blocks_prop_pack_args():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "property_package_args": {"test": 1},
                                       "outlet_list": ["foo", "bar"]})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)

    assert isinstance(m.fs.sep.foo_state, StateBlock)
    assert isinstance(m.fs.sep.bar_state, StateBlock)

    assert len(outlet_blocks) == 2
    for o in outlet_blocks:
        assert isinstance(o, StateBlock)
        assert o.local_name in ["foo_state", "bar_state"]
        assert o[0].config.has_phase_equilibrium is False
        assert o[0].config.defined_state is False
        assert len(o[0].config) == 4
        assert o[0].config.test == 1


def test_add_mixed_state_block():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    mixed_block = m.fs.sep.add_mixed_state_block()

    assert isinstance(mixed_block, StateBlock)
    assert hasattr(m.fs.sep, "mixed_state")
    assert m.fs.sep.mixed_state[0].config.has_phase_equilibrium is False
    assert m.fs.sep.mixed_state[0].config.defined_state is True
    assert len(m.fs.sep.mixed_state[0].config) == 3


def test_add_mixed_state_block_prop_pack_args():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "property_package_args": {"test": 1}})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    mixed_block = m.fs.sep.add_mixed_state_block()

    assert isinstance(mixed_block, StateBlock)
    assert hasattr(m.fs.sep, "mixed_state")
    assert m.fs.sep.mixed_state[0].config.has_phase_equilibrium is False
    assert m.fs.sep.mixed_state[0].config.defined_state is True
    assert len(m.fs.sep.mixed_state[0].config) == 4
    assert m.fs.sep.mixed_state[0].config.test == 1


def test_get_mixed_state_block():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    mixed_block = m.fs.sep.get_mixed_state_block()

    assert mixed_block == m.fs.sb


def test_get_mixed_state_block_none():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    with pytest.raises(BurntToast):
        m.fs.sep.get_mixed_state_block()


def test_get_mixed_state_block_mismatch():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    # Change parameters arg to create mismatch
    m.fs.sb[0].config.parameters = None

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    with pytest.raises(ConfigurationError):
        m.fs.sep.get_mixed_state_block()


# -----------------------------------------------------------------------------
# Test non-ideal splitting equation methods
def test_add_split_fractions_total():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()

    m.fs.sep.add_split_fractions(outlet_list)

    assert isinstance(m.fs.sep.outlet_idx, Set)
    assert len(m.fs.sep.outlet_idx) == len(outlet_list)
    assert isinstance(m.fs.sep.split_fraction, Var)
    assert len(m.fs.sep.split_fraction) == 2
    assert isinstance(m.fs.sep.sum_split_frac, Constraint)
    assert len(m.fs.sep.sum_split_frac) == 1


def test_add_split_fractions_phase():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis": SplittingType.phaseFlow})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()

    m.fs.sep.add_split_fractions(outlet_list)

    assert isinstance(m.fs.sep.outlet_idx, Set)
    assert len(m.fs.sep.outlet_idx) == len(outlet_list)
    assert isinstance(m.fs.sep.split_fraction, Var)
    assert len(m.fs.sep.split_fraction) == 4

    for t in m.fs.time:
        for o in m.fs.sep.outlet_idx:
            for p in m.fs.sep.config.property_package.phase_list:
                assert m.fs.sep.split_fraction[t, o, p] == 0.5

    assert isinstance(m.fs.sep.sum_split_frac, Constraint)
    assert len(m.fs.sep.sum_split_frac) == 2


def test_add_split_fractions_component():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis":
                                           SplittingType.componentFlow})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()

    m.fs.sep.add_split_fractions(outlet_list)

    assert isinstance(m.fs.sep.outlet_idx, Set)
    assert len(m.fs.sep.outlet_idx) == len(outlet_list)
    assert isinstance(m.fs.sep.split_fraction, Var)
    assert len(m.fs.sep.split_fraction) == 4

    for t in m.fs.time:
        for o in m.fs.sep.outlet_idx:
            for j in m.fs.sep.config.property_package.component_list:
                assert m.fs.sep.split_fraction[t, o, j] == 0.5

    assert isinstance(m.fs.sep.sum_split_frac, Constraint)
    assert len(m.fs.sep.sum_split_frac) == 2


def test_add_split_fractions_phase_component():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis":
                                           SplittingType.phaseComponentFlow})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()

    m.fs.sep.add_split_fractions(outlet_list)

    assert isinstance(m.fs.sep.outlet_idx, Set)
    assert len(m.fs.sep.outlet_idx) == len(outlet_list)
    assert isinstance(m.fs.sep.split_fraction, Var)
    assert len(m.fs.sep.split_fraction) == 8

    for t in m.fs.time:
        for o in m.fs.sep.outlet_idx:
            for p in m.fs.sep.config.property_package.phase_list:
                for j in m.fs.sep.config.property_package.component_list:
                    assert m.fs.sep.split_fraction[t, o, p, j] == 0.5

    assert isinstance(m.fs.sep.sum_split_frac, Constraint)
    assert len(m.fs.sep.sum_split_frac) == 4


def test_add_material_splitting_constraints_pc_total_no_equil():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert not hasattr(m.fs.sep, "phase_equilibrium_generation")


def test_add_material_splitting_constraints_pc_phase_no_equil():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis": SplittingType.phaseFlow})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert not hasattr(m.fs.sep, "phase_equilibrium_generation")


def test_add_material_splitting_constraints_pc_component_no_equil():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis":
                                           SplittingType.componentFlow})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert not hasattr(m.fs.sep, "phase_equilibrium_generation")


def test_add_material_splitting_constraints_pc_phase_component_no_equil():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis":
                                           SplittingType.phaseComponentFlow})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert not hasattr(m.fs.sep, "phase_equilibrium_generation")


def test_add_material_splitting_constraints_pc_total_equil():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "has_phase_equilibrium": True})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert isinstance(m.fs.sep.phase_equilibrium_generation, Var)
    assert len(m.fs.sep.phase_equilibrium_generation) == 2


def test_add_material_splitting_constraints_pc_phase_equil():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis": SplittingType.phaseFlow,
                                       "has_phase_equilibrium": True})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert isinstance(m.fs.sep.phase_equilibrium_generation, Var)
    assert len(m.fs.sep.phase_equilibrium_generation) == 2


def test_add_material_splitting_constraints_pc_component_equil():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis":
                                           SplittingType.componentFlow,
                                       "has_phase_equilibrium": True})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert isinstance(m.fs.sep.phase_equilibrium_generation, Var)
    assert len(m.fs.sep.phase_equilibrium_generation) == 2


def test_add_material_splitting_constraints_pc_phase_component_equil():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "mixed_state_block": m.fs.sb,
                                       "split_basis":
                                           SplittingType.phaseComponentFlow,
                                       "has_phase_equilibrium": True})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert isinstance(m.fs.sep.phase_equilibrium_generation, Var)
    assert len(m.fs.sep.phase_equilibrium_generation) == 2


def test_add_material_splitting_constraints_tc_total():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "material_balance_type": MaterialBalanceType.componentTotal})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 4


def test_add_material_splitting_constraints_tc_phase():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.phaseFlow,
            "material_balance_type": MaterialBalanceType.componentTotal})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 4


def test_add_material_splitting_constraints_tc_component():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.componentFlow,
            "material_balance_type": MaterialBalanceType.componentTotal})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 4


def test_add_material_splitting_constraints_tc_phase_component():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.phaseComponentFlow,
            "material_balance_type": MaterialBalanceType.componentTotal})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 4


def test_add_material_splitting_constraints_t_total():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "material_balance_type": MaterialBalanceType.total})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 2


def test_add_material_splitting_constraints_t_phase():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.phaseFlow,
            "material_balance_type": MaterialBalanceType.total})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 2


def test_add_material_splitting_constraints_t_component():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.componentFlow,
            "material_balance_type": MaterialBalanceType.total})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 2


def test_add_material_splitting_constraints_t_phase_component():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.phaseComponentFlow,
            "material_balance_type": MaterialBalanceType.total})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 2


def test_add_material_splitting_constraints_te_total():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "material_balance_type": MaterialBalanceType.elementTotal})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    with pytest.raises(ConfigurationError):
        m.fs.sep.add_material_splitting_constraints(mixed_block)


def test_add_material_splitting_constraints_none_total():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "material_balance_type": MaterialBalanceType.none})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert not hasattr(m.fs.sep, "material_splitting_eqn")


def test_add_material_splitting_constraints_none_phase():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.phaseFlow,
            "material_balance_type": MaterialBalanceType.none})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert not hasattr(m.fs.sep, "material_splitting_eqn")


def test_add_material_splitting_constraints_none_component():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.componentFlow,
            "material_balance_type": MaterialBalanceType.none})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert not hasattr(m.fs.sep, "material_splitting_eqn")


def test_add_material_splitting_constraints_none_phase_component():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "split_basis": SplittingType.phaseComponentFlow,
            "material_balance_type": MaterialBalanceType.none})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)

    m.fs.sep.add_material_splitting_constraints(mixed_block)

    assert not hasattr(m.fs.sep, "material_splitting_eqn")


def test_add_energy_splitting_constraints():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb})
    assert(m.fs.sep.config.energy_split_basis ==
           EnergySplittingType.equal_temperature)
    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)
    m.fs.sep.add_energy_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.temperature_equality_eqn, Constraint)
    assert len(m.fs.sep.temperature_equality_eqn) == 2


def test_add_energy_splitting_constraints_enthalpy():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={
        "property_package": m.fs.pp,
        "mixed_state_block": m.fs.sb,
        "energy_split_basis": EnergySplittingType.equal_molar_enthalpy})
    assert(m.fs.sep.config.energy_split_basis ==
           EnergySplittingType.equal_molar_enthalpy)
    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)
    m.fs.sep.add_energy_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.molar_enthalpy_equality_eqn, Constraint)
    assert len(m.fs.sep.molar_enthalpy_equality_eqn) == 2


def test_add_momentum_splitting_constraints():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                   "mixed_state_block": m.fs.sb})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.get_mixed_state_block()
    m.fs.sep.add_split_fractions(outlet_list)
    m.fs.sep.add_momentum_splitting_constraints(mixed_block)

    assert isinstance(m.fs.sep.pressure_equality_eqn, Constraint)
    assert len(m.fs.sep.pressure_equality_eqn) == 2


def test_add_inlet_port_objects():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.add_mixed_state_block()

    m.fs.sep.add_inlet_port_objects(mixed_block)

    assert isinstance(m.fs.sep.inlet, Port)


def test_add_inlet_port_objects_construct_ports_False():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "construct_ports": False})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.add_mixed_state_block()

    m.fs.sep.add_inlet_port_objects(mixed_block)

    assert hasattr(m.fs.sep, "inlet") is False


def test_add_outlet_port_objects():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.add_mixed_state_block()

    m.fs.sep.add_outlet_port_objects(outlet_list, outlet_blocks)

    assert isinstance(m.fs.sep.outlet_1, Port)
    assert isinstance(m.fs.sep.outlet_2, Port)


def test_add_outlet_port_objects_construct_ports_False():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp,
                                       "construct_ports": False})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.add_mixed_state_block()

    m.fs.sep.add_outlet_port_objects(outlet_list, outlet_blocks)

    assert hasattr(m.fs.sep, "outlet_1") is False
    assert hasattr(m.fs.sep, "outlet_2") is False


def test_build_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = Separator(default={"property_package": m.fs.pp,
                                  "ideal_separation": False})

    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
    assert len(m.fs.sep.material_splitting_eqn) == 8
    assert isinstance(m.fs.sep.temperature_equality_eqn, Constraint)
    assert len(m.fs.sep.temperature_equality_eqn) == 2
    assert isinstance(m.fs.sep.pressure_equality_eqn, Constraint)
    assert len(m.fs.sep.pressure_equality_eqn) == 2

    assert isinstance(m.fs.sep.outlet_1, Port)
    assert isinstance(m.fs.sep.outlet_2, Port)
    assert isinstance(m.fs.sep.inlet, Port)


def test_model_checks():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = Separator(default={
            "property_package": m.fs.pp,
            "ideal_separation": False})

    m.fs.sep.model_check()

    assert m.fs.sep.outlet_1_state[0].check is True
    assert m.fs.sep.outlet_2_state[0].check is True
    assert m.fs.sep.mixed_state[0].check is True


def test_initialize():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = Separator(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "ideal_separation": False,
            "split_basis": SplittingType.phaseComponentFlow})

    # Change one outlet pressure to check initialization calculations
    m.fs.sep.outlet_1_state[0].pressure = 8e4

    f = m.fs.sep.initialize(hold_state=True)

    assert m.fs.sep.outlet_1_state[0].init_test is True
    assert m.fs.sep.outlet_2_state[0].init_test is True
    assert m.fs.sb[0].init_test is True
    assert m.fs.sep.outlet_1_state[0].hold_state is False
    assert m.fs.sep.outlet_2_state[0].hold_state is False
    assert m.fs.sb[0].hold_state is True

    assert m.fs.sep.outlet_1.component_flow[0, "p1", "c1"].value == 0.5
    assert m.fs.sep.outlet_1.component_flow[0, "p1", "c2"].value == 0.5
    assert m.fs.sep.outlet_1.component_flow[0, "p2", "c1"].value == 0.5
    assert m.fs.sep.outlet_1.component_flow[0, "p2", "c2"].value == 0.5
    assert m.fs.sep.outlet_1.enthalpy[0, "p1"].value == 2
    assert m.fs.sep.outlet_1.enthalpy[0, "p2"].value == 2
    assert m.fs.sep.outlet_1.pressure[0].value == 1e5

    assert m.fs.sep.outlet_2.component_flow[0, "p1", "c1"].value == 0.5
    assert m.fs.sep.outlet_2.component_flow[0, "p1", "c2"].value == 0.5
    assert m.fs.sep.outlet_2.component_flow[0, "p2", "c1"].value == 0.5
    assert m.fs.sep.outlet_2.component_flow[0, "p2", "c2"].value == 0.5
    assert m.fs.sep.outlet_2.enthalpy[0, "p1"].value == 2
    assert m.fs.sep.outlet_2.enthalpy[0 ,"p2"].value == 2
    assert m.fs.sep.outlet_2.pressure[0].value == 1e5

    m.fs.sep.release_state(flags=f)

    assert m.fs.sep.outlet_1_state[0].hold_state is False
    assert m.fs.sep.outlet_2_state[0].hold_state is False
    assert m.fs.sb[0].hold_state is False


def test_initialize_inconsistent_keys():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})

    m.fs.sep = Separator(default={
            "property_package": m.fs.pp,
            "mixed_state_block": m.fs.sb,
            "ideal_separation": False,
            "split_basis": SplittingType.phaseFlow})

    # Change one outlet pressure to check initialization calculations
    m.fs.sep.outlet_1_state[0].pressure = 8e4

    with pytest.raises(KeyError):
        m.fs.sep.initialize()


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize_total_flow():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()

    m.fs.sb = Separator(default={
            "property_package": m.fs.properties,
            "ideal_separation": False,
            "split_basis": SplittingType.totalFlow})

    m.fs.sb.inlet.flow_vol.fix(1.0e-03)
    m.fs.sb.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.sb.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.sb.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.sb.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.sb.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.sb.inlet.temperature.fix(303.15)
    m.fs.sb.inlet.pressure.fix(101325.0)

    m.fs.sb.split_fraction[0, "outlet_1"].fix(0.2)

    assert degrees_of_freedom(m) == 0

    m.fs.sb.initialize(outlvl=5, optarg={'tol': 1e-6})

    assert (pytest.approx(0.2, abs=1e-3) ==
             m.fs.sb.split_fraction[0, "outlet_1"].value)
    assert (pytest.approx(0.8, abs=1e-3) ==
             m.fs.sb.split_fraction[0, "outlet_2"].value)

    assert (pytest.approx(101325.0, abs=1e-2) ==
            m.fs.sb.outlet_1.pressure[0].value)
    assert (pytest.approx(303.15, abs=1e-2) ==
            m.fs.sb.outlet_1.temperature[0].value)
    assert (pytest.approx(2e-4, abs=1e-6) ==
            m.fs.sb.outlet_1.flow_vol[0].value)
    assert (pytest.approx(55388.0, abs=1e-2) ==
            m.fs.sb.outlet_1.conc_mol_comp[0, "H2O"].value)
    assert (pytest.approx(100.0, abs=1e-2) ==
            m.fs.sb.outlet_1.conc_mol_comp[0, "NaOH"].value)
    assert (pytest.approx(100.0, abs=1e-2) ==
            m.fs.sb.outlet_1.conc_mol_comp[0, "EthylAcetate"].value)
    assert (pytest.approx(0.0, abs=1e-2) ==
            m.fs.sb.outlet_1.conc_mol_comp[0, "SodiumAcetate"].value)
    assert (pytest.approx(0.0, abs=1e-2) ==
            m.fs.sb.outlet_1.conc_mol_comp[0, "Ethanol"].value)

    assert (pytest.approx(101325.0, abs=1e-2) ==
            m.fs.sb.outlet_2.pressure[0].value)
    assert (pytest.approx(303.15, abs=1e-2) ==
            m.fs.sb.outlet_2.temperature[0].value)
    assert (pytest.approx(8e-4, abs=1e-6) ==
            m.fs.sb.outlet_2.flow_vol[0].value)
    assert (pytest.approx(55388.0, abs=1e-2) ==
            m.fs.sb.outlet_2.conc_mol_comp[0, "H2O"].value)
    assert (pytest.approx(100.0, abs=1e-2) ==
            m.fs.sb.outlet_2.conc_mol_comp[0, "NaOH"].value)
    assert (pytest.approx(100.0, abs=1e-2) ==
            m.fs.sb.outlet_2.conc_mol_comp[0, "EthylAcetate"].value)
    assert (pytest.approx(0.0, abs=1e-2) ==
            m.fs.sb.outlet_2.conc_mol_comp[0, "SodiumAcetate"].value)
    assert (pytest.approx(0.0, abs=1e-2) ==
            m.fs.sb.outlet_2.conc_mol_comp[0, "Ethanol"].value)


def test_report():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()

    m.fs.sb = Separator(default={
            "property_package": m.fs.properties,
            "ideal_separation": False,
            "split_basis": SplittingType.totalFlow})

    m.fs.sb.report()


# -----------------------------------------------------------------------------
# Testing of ideal splitting methods
@declare_process_block_class("PhysicalParameterBlock2")
class _PhysicalParameterBlock2(PhysicalParameterBlock):
    def build(self):
        super(_PhysicalParameterBlock2, self).build()

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])
        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])

        self.state_block_class = StateBlock2

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


@declare_process_block_class("StateBlock2", block_class=SBlockBase)
class StateBlockData2(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(StateBlockData2, self).build()

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
        self.temperature = Var(initialize=5)

    def get_material_flow_terms(b, p, j):
        return b.flow_mol_phase_comp[p, j]

    def get_enthalpy_flow_terms(b, p):
        return b.enth_mol_phase[p]

    def define_state_vars(self):
        return {"component_flow": self.flow_mol_phase_comp,
                "temperature": self.temperature,
                "pressure": self.pressure}

    def model_check(self):
        self.check = True


#def test_ideal_splitting_methods():
#    m = ConcreteModel()
#    m.fs = Flowsheet(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterBlock2()
#
#    m.fs.sep = SeparatorFrame(default={
#            "property_package": m.fs.pp,
#            "split_basis": SplittingType.phaseComponentFlow,
#            "ideal_split_map": {("p1", "c1"): "outlet_1",
#                                ("p1", "c2"): "outlet_2",
#                                ("p2", "c1"): "outlet_2",
#                                ("p2", "c2"): "outlet_2"}})
#
#    m.fs.sep._get_property_package()
#    m.fs.sep._get_indexing_sets()
#
#    outlet_list = m.fs.sep.create_outlet_list()
#    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
#    mixed_block = m.fs.sep.add_mixed_state_block()
#
#    m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)
#
#    assert isinstance(m.fs.sep.outlet_1, Port)
#    assert isinstance(m.fs.sep.outlet_2, Port)
#
#    assert m.fs.sep.outlet_1[0].component_flow["p1", "c1"].value == 1.0
#    assert m.fs.sep.outlet_1[0].component_flow["p1", "c2"].value == 0.0
#    assert m.fs.sep.outlet_1[0].component_flow["p2", "c1"].value == 0.0
#    assert m.fs.sep.outlet_1[0].component_flow["p2", "c2"].value == 0.0
#    assert m.fs.sep.outlet_1[0].temperature.value == 5
#    assert m.fs.sep.outlet_1[0].pressure.value == 1e5
#
#    assert m.fs.sep.outlet_2[0].component_flow["p1", "c1"].value == 0.0
#    assert m.fs.sep.outlet_2[0].component_flow["p1", "c2"].value == 1.0
#    assert m.fs.sep.outlet_2[0].component_flow["p2", "c1"].value == 1.0
#    assert m.fs.sep.outlet_2[0].component_flow["p2", "c2"].value == 1.0
#    assert m.fs.sep.outlet_2[0].temperature.value == 5
#    assert m.fs.sep.outlet_2[0].pressure.value == 1e5
#
#
#def test_ideal_splitting_methods_key_error():
#    m = ConcreteModel()
#    m.fs = Flowsheet(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterBlock2()
#
#    m.fs.sep = SeparatorFrame(default={
#            "property_package": m.fs.pp,
#            "split_basis": SplittingType.componentFlow,
#            "ideal_split_map": {("p1", "c1"): "outlet_1",
#                                ("p1", "c2"): "outlet_2",
#                                ("p2", "c1"): "outlet_2",
#                                ("p2", "c2"): "outlet_2"}})
#
#    m.fs.sep._get_property_package()
#    m.fs.sep._get_indexing_sets()
#
#    outlet_list = m.fs.sep.create_outlet_list()
#    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
#    mixed_block = m.fs.sep.add_mixed_state_block()
#
#    with pytest.raises(KeyError):
#        m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)


def test_ideal_splitting_methods_no_map():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock2()

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "split_basis": SplittingType.phaseComponentFlow})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.add_mixed_state_block()

    with pytest.raises(ConfigurationError):
        m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)


def test_ideal_splitting_methods_totalFlow():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock2()

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "split_basis": SplittingType.totalFlow,
            "ideal_split_map": {("p1", "c1"): "outlet_1",
                                ("p1", "c2"): "outlet_2",
                                ("p2", "c1"): "outlet_2",
                                ("p2", "c2"): "outlet_2"}})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.add_mixed_state_block()

    with pytest.raises(ConfigurationError):
        m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)


def test_ideal_splitting_methods_no_ports():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock2()

    m.fs.sep = SeparatorFrame(default={
            "property_package": m.fs.pp,
            "split_basis": SplittingType.phaseComponentFlow,
            "ideal_split_map": {("p1", "c1"): "outlet_1",
                                ("p1", "c2"): "outlet_2",
                                ("p2", "c1"): "outlet_2",
                                ("p2", "c2"): "outlet_2"},
                                "construct_ports": False})

    m.fs.sep._get_property_package()
    m.fs.sep._get_indexing_sets()

    outlet_list = m.fs.sep.create_outlet_list()
    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
    mixed_block = m.fs.sep.add_mixed_state_block()

    with pytest.raises(ConfigurationError):
        m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)

@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS unavailable")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_sep_ph_mixed_default():
    """Test mixed phase form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(default={})
    m.fs.sep = Separator(default={
        "property_package":m.fs.properties,
        "split_basis":SplittingType.totalFlow,
        "ideal_separation":False,
        "num_outlets":2,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})

    m.fs.sep.inlet.enth_mol.fix(24000)
    m.fs.sep.inlet.flow_mol.fix(100)
    m.fs.sep.inlet.pressure.fix(101325)
    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
    m.fs.sep.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.sep.mixed_state[0]
    prop_out = m.fs.sep.outlet_1_state[0]
    prop_out2 =  m.fs.sep.outlet_2_state[0]
    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
    assert abs(value(prop_out.flow_mol - 90.0) <= 1e-6)
    assert abs(value(prop_out2.flow_mol - 10.0) <= 1e-6)

@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS unavailable")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_sep_ph_lg_default():
    """Test mixed phase form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
        "phase_presentation":iapws95.PhaseType.LG})
    m.fs.sep = Separator(default={
        "property_package":m.fs.properties,
        "split_basis":SplittingType.totalFlow,
        "ideal_separation":False,
        "num_outlets":2,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})

    m.fs.sep.inlet.enth_mol.fix(24000)
    m.fs.sep.inlet.flow_mol.fix(100)
    m.fs.sep.inlet.pressure.fix(101325)
    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
    m.fs.sep.initialize()
    # by default the separator writes phase mass balances so -2 degress of
    # freedom is right.  for IAPWS need total mass balance.
    assert degrees_of_freedom(m) == -2

@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS unavailable")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_sep_ph_lg_total_mb():
    """Test mixed phase form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
        "phase_presentation":iapws95.PhaseType.LG})
    m.fs.sep = Separator(default={
        "property_package":m.fs.properties,
        "split_basis":SplittingType.totalFlow,
        "ideal_separation":False,
        "num_outlets":2,
        "material_balance_type":MaterialBalanceType.total,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})

    m.fs.sep.inlet.enth_mol.fix(24000)
    m.fs.sep.inlet.flow_mol.fix(100)
    m.fs.sep.inlet.pressure.fix(101325)
    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
    m.fs.sep.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.sep.mixed_state[0]
    prop_out = m.fs.sep.outlet_1_state[0]
    prop_out2 =  m.fs.sep.outlet_2_state[0]
    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
    assert abs(value(prop_out.flow_mol - 90.0) <= 1e-6)
    assert abs(value(prop_out2.flow_mol - 10.0) <= 1e-6)

@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS unavailable")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_sep_ph_lg_has_phase_eq():
    """Test mixed phase form with P-H state vars and phase mass balances
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
        "phase_presentation":iapws95.PhaseType.LG})
    m.fs.sep = Separator(default={
        "has_phase_equilibrium":True,
        "property_package":m.fs.properties,
        "split_basis":SplittingType.totalFlow,
        "ideal_separation":False,
        "num_outlets":2,
        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})

    m.fs.sep.inlet.enth_mol.fix(24000)
    m.fs.sep.inlet.flow_mol.fix(100)
    m.fs.sep.inlet.pressure.fix(101325)
    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
    m.fs.sep.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.sep.mixed_state[0]
    prop_out = m.fs.sep.outlet_1_state[0]
    prop_out2 =  m.fs.sep.outlet_2_state[0]
    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
    assert abs(value(prop_out.flow_mol - 90.0) <= 1e-6)
    assert abs(value(prop_out2.flow_mol - 10.0) <= 1e-6)
