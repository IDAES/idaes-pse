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
Tests for Separator unit model.

Author: Andrew Lee
"""
import pytest

from pyomo.environ import (ConcreteModel,
                           Constraint,
                           TerminationCondition,
                           Set,
                           SolverStatus,
                           value,
                           Var)

from pyomo.network import Port
from pyomo.common.config import ConfigBlock

from idaes.core import (FlowsheetBlock,
                        declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlock,
                        MaterialBalanceType)
from idaes.unit_models.separator import (Separator,
                                         SeparatorData,
                                         SplittingType,
                                         EnergySplittingType)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError)

from idaes.property_models.examples.saponification_thermo import (
    SaponificationParameterBlock)
from idaes.property_models.ideal.BTX_ideal_VLE import BTXParameterBlock
from idaes.property_models import iapws95

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from idaes.core.util.testing import (get_default_solver,
                                     PhysicalParameterTestBlock,
                                     TestStateBlock)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = Separator(default={"property_package": m.fs.properties})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 14

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == \
        MaterialBalanceType.componentPhase
    assert m.fs.unit.config.energy_split_basis == \
        EnergySplittingType.equal_temperature
    assert m.fs.unit.config.outlet_list is None
    assert m.fs.unit.config.num_outlets == 2
    assert m.fs.unit.config.split_basis == SplittingType.totalFlow
    assert not m.fs.unit.config.ideal_separation
    assert m.fs.unit.config.ideal_split_map is None
    assert m.fs.unit.config.mixed_state_block is None
    assert m.fs.unit.config.construct_ports
    assert not m.fs.unit.config.has_phase_equilibrium
    assert m.fs.unit.config.property_package is m.fs.properties


# -----------------------------------------------------------------------------
# Mockup classes for testing
@declare_process_block_class("SeparatorFrame")
class SeparatorFrameData(SeparatorData):
    def build(self):
        super(SeparatorData, self).build()


# -----------------------------------------------------------------------------
# Tests of Separator unit model construction methods
@pytest.mark.build
class TestBaseConstruction(object):
    @pytest.fixture(scope="function")
    def build(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

        return m

    def test_separator_config(self, build):
        assert len(build.fs.sep.config) == 14
        assert build.fs.sep.config.dynamic is False
        assert build.fs.sep.config.has_holdup is False
        assert build.fs.sep.config.property_package == build.fs.pp
        assert isinstance(build.fs.sep.config.property_package_args,
                          ConfigBlock)
        assert len(build.fs.sep.config.property_package_args) == 0
        assert build.fs.sep.config.outlet_list is None
        assert build.fs.sep.config.num_outlets is None
        assert build.fs.sep.config.split_basis == SplittingType.totalFlow
        assert build.fs.sep.config.ideal_separation is False
        assert build.fs.sep.config.ideal_split_map is None
        assert build.fs.sep.config.mixed_state_block is None
        assert build.fs.sep.config.construct_ports is True
        assert build.fs.sep.config.material_balance_type == \
            MaterialBalanceType.componentPhase
        assert build.fs.sep.config.has_phase_equilibrium is False

    def test_validate_config_arguments(self, build):
        build.fs.sep.config.has_phase_equilibrium = True
        build.fs.sep.config.ideal_separation = True

        with pytest.raises(ConfigurationError):
            build.fs.sep._validate_config_arguments()

    def test_create_outlet_list_default(self, build):
        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()

        for o in outlet_list:
            assert o in ["outlet_1", "outlet_2"]

    def test_create_outlet_list_outlet_list(self, build):
        build.fs.sep.config.outlet_list = ["foo", "bar"]

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()

        for o in outlet_list:
            assert o in ["foo", "bar"]

    def test_create_outlet_list_num_outlets(self, build):
        build.fs.sep.config.num_outlets = 3

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()

        for o in outlet_list:
            assert o in ["outlet_1", "outlet_2", "outlet_3"]

    def test_create_outlet_list_both_args_consistent(self, build):
        build.fs.sep.config.outlet_list = ["foo", "bar"]
        build.fs.sep.config.num_outlets = 2

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()

        for o in outlet_list:
            assert o in ["foo", "bar"]

    def test_create_outlet_list_both_args_inconsistent(self, build):
        build.fs.sep.config.outlet_list = ["foo", "bar"]
        build.fs.sep.config.num_outlets = 3

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        with pytest.raises(ConfigurationError):
            build.fs.sep.create_outlet_list()

    def test_add_outlet_state_blocks(self, build):
        build.fs.sep.config.outlet_list = ["foo", "bar"]

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()
        outlet_blocks = build.fs.sep.add_outlet_state_blocks(outlet_list)

        assert isinstance(build.fs.sep.foo_state, StateBlock)
        assert isinstance(build.fs.sep.bar_state, StateBlock)

        assert len(outlet_blocks) == 2
        for o in outlet_blocks:
            assert isinstance(o, StateBlock)
            assert o.local_name in ["foo_state", "bar_state"]
            assert o[0].config.has_phase_equilibrium is False
            assert o[0].config.defined_state is False
            assert len(o[0].config) == 3

    def test_add_outlet_state_blocks_prop_pack_args(self, build):
        build.fs.sep.config.property_package_args = {"test": 1}
        build.fs.sep.config.outlet_list = ["foo", "bar"]

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()
        outlet_blocks = build.fs.sep.add_outlet_state_blocks(outlet_list)

        assert isinstance(build.fs.sep.foo_state, StateBlock)
        assert isinstance(build.fs.sep.bar_state, StateBlock)

        assert len(outlet_blocks) == 2
        for o in outlet_blocks:
            assert isinstance(o, StateBlock)
            assert o.local_name in ["foo_state", "bar_state"]
            assert o[0].config.has_phase_equilibrium is False
            assert o[0].config.defined_state is False
            assert len(o[0].config) == 4
            assert o[0].config.test == 1

    def test_add_mixed_state_block(self, build):
        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        mixed_block = build.fs.sep.add_mixed_state_block()

        assert isinstance(mixed_block, StateBlock)
        assert hasattr(build.fs.sep, "mixed_state")
        assert not build.fs.sep.mixed_state[0].config.has_phase_equilibrium
        assert build.fs.sep.mixed_state[0].config.defined_state
        assert len(build.fs.sep.mixed_state[0].config) == 3

    def test_add_mixed_state_block_prop_pack_args(self, build):
        build.fs.sep.config.property_package_args = {"test": 1}

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        mixed_block = build.fs.sep.add_mixed_state_block()

        assert isinstance(mixed_block, StateBlock)
        assert hasattr(build.fs.sep, "mixed_state")
        assert not build.fs.sep.mixed_state[0].config.has_phase_equilibrium
        assert build.fs.sep.mixed_state[0].config.defined_state
        assert len(build.fs.sep.mixed_state[0].config) == 4
        assert build.fs.sep.mixed_state[0].config.test == 1

    def test_get_mixed_state_block(self, build):
        build.fs.sb = TestStateBlock(build.fs.time,
                                     default={"parameters": build.fs.pp})

        build.fs.sep.config.mixed_state_block = build.fs.sb

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        mixed_block = build.fs.sep.get_mixed_state_block()

        assert mixed_block == build.fs.sb

    def test_get_mixed_state_block_none(self, build):
        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        with pytest.raises(BurntToast):
            build.fs.sep.get_mixed_state_block()

    def test_get_mixed_state_block_mismatch(self, build):
        build.fs.sb = TestStateBlock(build.fs.time,
                                     default={"parameters": build.fs.pp})

        # Change parameters arg to create mismatch
        build.fs.sb[0].config.parameters = None

        build.fs.sep.config.mixed_state_block = build.fs.sb

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        with pytest.raises(ConfigurationError):
            build.fs.sep.get_mixed_state_block()


# -----------------------------------------------------------------------------
# Tests of Separator unit model non-ideal construction methods
@pytest.mark.build
class TestSplitConstruction(object):
    @pytest.fixture(scope="function")
    def build(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(default={"property_package": m.fs.pp})

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.outlet_blocks = m.fs.sep.add_outlet_state_blocks(m.outlet_list)
        m.fs.sep.add_mixed_state_block()

        return m

    def test_add_split_fractions_total(self, build):
        build.fs.sep.add_split_fractions(build.outlet_list)

        assert isinstance(build.fs.sep.outlet_idx, Set)
        assert len(build.fs.sep.outlet_idx) == len(build.outlet_list)
        assert isinstance(build.fs.sep.split_fraction, Var)
        assert len(build.fs.sep.split_fraction) == 2
        assert isinstance(build.fs.sep.sum_split_frac, Constraint)
        assert len(build.fs.sep.sum_split_frac) == 1

    def test_add_split_fractions_phase(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow

        build.fs.sep.add_split_fractions(build.outlet_list)

        assert isinstance(build.fs.sep.outlet_idx, Set)
        assert len(build.fs.sep.outlet_idx) == len(build.outlet_list)
        assert isinstance(build.fs.sep.split_fraction, Var)
        assert len(build.fs.sep.split_fraction) == 4

        for t in build.fs.time:
            for o in build.fs.sep.outlet_idx:
                for p in build.fs.sep.config.property_package.phase_list:
                    assert build.fs.sep.split_fraction[t, o, p] == 0.5

        assert isinstance(build.fs.sep.sum_split_frac, Constraint)
        assert len(build.fs.sep.sum_split_frac) == 2

    def test_add_split_fractions_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow

        build.fs.sep.add_split_fractions(build.outlet_list)

        assert isinstance(build.fs.sep.outlet_idx, Set)
        assert len(build.fs.sep.outlet_idx) == len(build.outlet_list)
        assert isinstance(build.fs.sep.split_fraction, Var)
        assert len(build.fs.sep.split_fraction) == 4

        for t in build.fs.time:
            for o in build.fs.sep.outlet_idx:
                for j in build.fs.sep.config.property_package.component_list:
                    assert build.fs.sep.split_fraction[t, o, j] == 0.5

        assert isinstance(build.fs.sep.sum_split_frac, Constraint)
        assert len(build.fs.sep.sum_split_frac) == 2

    def test_add_split_fractions_phase_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow

        build.fs.sep.add_split_fractions(build.outlet_list)

        assert isinstance(build.fs.sep.outlet_idx, Set)
        assert len(build.fs.sep.outlet_idx) == len(build.outlet_list)
        assert isinstance(build.fs.sep.split_fraction, Var)
        assert len(build.fs.sep.split_fraction) == 8

        for t in build.fs.time:
            for o in build.fs.sep.outlet_idx:
                for p in build.fs.sep.config.property_package.phase_list:
                    for j in build.fs.sep.config.property_package.component_list:
                        assert build.fs.sep.split_fraction[t, o, p, j] == 0.5

        assert isinstance(build.fs.sep.sum_split_frac, Constraint)
        assert len(build.fs.sep.sum_split_frac) == 4

    def test_add_material_splitting_constraints_pc_total_no_equil(self, build):
        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert not hasattr(build.fs.sep, "phase_equilibrium_generation")

    def test_add_material_splitting_constraints_pc_phase_no_equil(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert not hasattr(build.fs.sep, "phase_equilibrium_generation")

    def test_add_material_splitting_constraints_pc_component_no_equil(self,
                                                                      build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert not hasattr(build.fs.sep, "phase_equilibrium_generation")

    def test_add_material_splitting_constraints_pc_phase_component_no_equil(
            self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert not hasattr(build.fs.sep, "phase_equilibrium_generation")

    def test_add_material_splitting_constraints_pc_total_equil(self, build):
        build.fs.sep.config.has_phase_equilibrium = True

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert isinstance(build.fs.sep.phase_equilibrium_generation, Var)
        assert len(build.fs.sep.phase_equilibrium_generation) == 4

    def test_add_material_splitting_constraints_pc_phase_equil(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow
        build.fs.sep.config.has_phase_equilibrium = True

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert isinstance(build.fs.sep.phase_equilibrium_generation, Var)
        assert len(build.fs.sep.phase_equilibrium_generation) == 4

    def test_add_material_splitting_constraints_pc_component_equil(self,
                                                                   build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow
        build.fs.sep.config.has_phase_equilibrium = True

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert isinstance(build.fs.sep.phase_equilibrium_generation, Var)
        assert len(build.fs.sep.phase_equilibrium_generation) == 4

    def test_add_material_splitting_constraints_pc_phase_component_equil(
            self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow
        build.fs.sep.config.has_phase_equilibrium = True

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert isinstance(build.fs.sep.phase_equilibrium_generation, Var)
        assert len(build.fs.sep.phase_equilibrium_generation) == 4

    def test_add_material_splitting_constraints_tc_total(self, build):
        build.fs.sep.config.material_balance_type = \
            MaterialBalanceType.componentTotal

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 4

    def test_add_material_splitting_constraints_tc_phase(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow
        build.fs.sep.config.material_balance_type = \
            MaterialBalanceType.componentTotal

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 4

    def test_add_material_splitting_constraints_tc_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow
        build.fs.sep.config.material_balance_type = \
            MaterialBalanceType.componentTotal

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 4

    def test_add_material_splitting_constraints_tc_phase_component(self,
                                                                   build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow
        build.fs.sep.config.material_balance_type = \
            MaterialBalanceType.componentTotal

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 4

    def test_add_material_splitting_constraints_t_total(self, build):
        build.fs.sep.config.material_balance_type = MaterialBalanceType.total

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 2

    def test_add_material_splitting_constraints_t_phase(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.total

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 2

    def test_add_material_splitting_constraints_t_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.total

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 2

    def test_add_material_splitting_constraints_t_phase_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.total

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 2

    def test_add_material_splitting_constraints_te_total(self, build):
        build.fs.sep.config.material_balance_type = \
            MaterialBalanceType.elementTotal

        build.fs.sep.add_split_fractions(build.outlet_list)

        with pytest.raises(ConfigurationError):
            build.fs.sep.add_material_splitting_constraints(
                    build.fs.sep.mixed_state)

    def test_add_material_splitting_constraints_none_total(self, build):
        build.fs.sep.config.material_balance_type = MaterialBalanceType.none

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert not hasattr(build.fs.sep, "material_splitting_eqn")

    def test_add_material_splitting_constraints_none_phase(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.none

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert not hasattr(build.fs.sep, "material_splitting_eqn")

    def test_add_material_splitting_constraints_none_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.none

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert not hasattr(build.fs.sep, "material_splitting_eqn")

    def test_add_material_splitting_constraints_none_phase_component(self,
                                                                     build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.none

        build.fs.sep.add_split_fractions(build.outlet_list)

        build.fs.sep.add_material_splitting_constraints(
                build.fs.sep.mixed_state)

        assert not hasattr(build.fs.sep, "material_splitting_eqn")

    def test_add_energy_splitting_constraints(self, build):
        assert(build.fs.sep.config.energy_split_basis ==
               EnergySplittingType.equal_temperature)

        build.fs.sep.add_split_fractions(build.outlet_list)
        build.fs.sep.add_energy_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.temperature_equality_eqn, Constraint)
        assert len(build.fs.sep.temperature_equality_eqn) == 2

    def test_add_energy_splitting_constraints_enthalpy(self, build):
        build.fs.sep.config.energy_split_basis = \
            EnergySplittingType.equal_molar_enthalpy
        assert(build.fs.sep.config.energy_split_basis ==
               EnergySplittingType.equal_molar_enthalpy)

        build.fs.sep.add_split_fractions(build.outlet_list)
        build.fs.sep.add_energy_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.molar_enthalpy_equality_eqn, Constraint)
        assert len(build.fs.sep.molar_enthalpy_equality_eqn) == 2

    def test_add_momentum_splitting_constraints(self, build):
        build.fs.sep.add_split_fractions(build.outlet_list)
        build.fs.sep.add_momentum_splitting_constraints(
                build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.pressure_equality_eqn, Constraint)
        assert len(build.fs.sep.pressure_equality_eqn) == 2

    def test_add_inlet_port_objects(self, build):
        build.fs.sep.add_inlet_port_objects(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.inlet, Port)

    def test_add_inlet_port_objects_construct_ports_False(self, build):
        build.fs.sep.config.construct_ports = False

        build.fs.sep.add_inlet_port_objects(build.fs.sep.mixed_state)

        assert hasattr(build.fs.sep, "inlet") is False

    def test_add_outlet_port_objects(self, build):
        build.fs.sep.add_outlet_port_objects(build.outlet_list,
                                             build.outlet_blocks)

        assert isinstance(build.fs.sep.outlet_1, Port)
        assert isinstance(build.fs.sep.outlet_2, Port)

    def test_add_outlet_port_objects_construct_ports_False(self, build):
        build.fs.sep.config.construct_ports = False

        build.fs.sep.add_outlet_port_objects(build.outlet_list,
                                             build.outlet_blocks)

        assert hasattr(build.fs.sep, "outlet_1") is False
        assert hasattr(build.fs.sep, "outlet_2") is False


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = Separator(default={
                "property_package": m.fs.properties,
                "material_balance_type": MaterialBalanceType.componentPhase,
                "split_basis": SplittingType.totalFlow,
                "outlet_list": ["a", "B", "c"],
                "ideal_separation": False,
                "has_phase_equilibrium": False})

        return m

    @pytest.mark.build
    def test_build(self, sapon):

        assert hasattr(sapon.fs.unit, "inlet")
        assert len(sapon.fs.unit.inlet.vars) == 4
        assert hasattr(sapon.fs.unit.inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.inlet, "temperature")
        assert hasattr(sapon.fs.unit.inlet, "pressure")

        assert hasattr(sapon.fs.unit, "a")
        assert len(sapon.fs.unit.a.vars) == 4
        assert hasattr(sapon.fs.unit.a, "flow_vol")
        assert hasattr(sapon.fs.unit.a, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.a, "temperature")
        assert hasattr(sapon.fs.unit.a, "pressure")

        assert hasattr(sapon.fs.unit, "B")
        assert len(sapon.fs.unit.B.vars) == 4
        assert hasattr(sapon.fs.unit.B, "flow_vol")
        assert hasattr(sapon.fs.unit.B, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.B, "temperature")
        assert hasattr(sapon.fs.unit.B, "pressure")

        assert hasattr(sapon.fs.unit, "c")
        assert len(sapon.fs.unit.c.vars) == 4
        assert hasattr(sapon.fs.unit.c, "flow_vol")
        assert hasattr(sapon.fs.unit.c, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.c, "temperature")
        assert hasattr(sapon.fs.unit.c, "pressure")

        assert isinstance(sapon.fs.unit.split_fraction, Var)

        assert number_variables(sapon) == 35
        assert number_total_constraints(sapon) == 25
        assert number_unused_variables(sapon) == 0

    def test_dof(self, sapon):
        sapon.fs.unit.inlet.flow_vol.fix(1)
        sapon.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        sapon.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        sapon.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        sapon.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        sapon.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        sapon.fs.unit.inlet.temperature.fix(303.15)
        sapon.fs.unit.inlet.pressure.fix(101325.0)

        sapon.fs.unit.split_fraction[0, "a"].fix(0.3)
        sapon.fs.unit.split_fraction[0, "B"].fix(0.5)

        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, sapon):
        orig_fixed_vars = fixed_variables_set(sapon)
        orig_act_consts = activated_constraints_set(sapon)

        sapon.fs.unit.initialize(optarg={'tol': 1e-6})

        assert degrees_of_freedom(sapon) == 0

        fin_fixed_vars = fixed_variables_set(sapon)
        fin_act_consts = activated_constraints_set(sapon)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, sapon):
        assert (pytest.approx(0.3, abs=1e-5) ==
                value(sapon.fs.unit.a.flow_vol[0]))
        assert (pytest.approx(101325.0, abs=1e-2) ==
                value(sapon.fs.unit.a.pressure[0]))
        assert (pytest.approx(303.15, abs=1e-2) ==
                value(sapon.fs.unit.a.temperature[0]))
        assert (pytest.approx(55388, abs=1e0) ==
                value(sapon.fs.unit.a.conc_mol_comp[0, "H2O"]))
        assert (pytest.approx(100.0, abs=1e-3) ==
                value(sapon.fs.unit.a.conc_mol_comp[0, "NaOH"]))
        assert (pytest.approx(100.0, abs=1e-3) ==
                value(sapon.fs.unit.a.conc_mol_comp[0, "EthylAcetate"]))
        assert (pytest.approx(0.0, abs=1e-3) ==
                value(sapon.fs.unit.a.conc_mol_comp[0, "SodiumAcetate"]))
        assert (pytest.approx(0.0, abs=1e-3) ==
                value(sapon.fs.unit.a.conc_mol_comp[0, "Ethanol"]))

        assert (pytest.approx(0.5, abs=1e-5) ==
                value(sapon.fs.unit.B.flow_vol[0]))
        assert (pytest.approx(101325.0, abs=1e-2) ==
                value(sapon.fs.unit.B.pressure[0]))
        assert (pytest.approx(303.15, abs=1e-2) ==
                value(sapon.fs.unit.B.temperature[0]))
        assert (pytest.approx(55388, abs=1e0) ==
                value(sapon.fs.unit.B.conc_mol_comp[0, "H2O"]))
        assert (pytest.approx(100.0, abs=1e-3) ==
                value(sapon.fs.unit.B.conc_mol_comp[0, "NaOH"]))
        assert (pytest.approx(100.0, abs=1e-3) ==
                value(sapon.fs.unit.B.conc_mol_comp[0, "EthylAcetate"]))
        assert (pytest.approx(0.0, abs=1e-3) ==
                value(sapon.fs.unit.B.conc_mol_comp[0, "SodiumAcetate"]))
        assert (pytest.approx(0.0, abs=1e-3) ==
                value(sapon.fs.unit.B.conc_mol_comp[0, "Ethanol"]))

        assert (pytest.approx(0.2, abs=1e-5) ==
                value(sapon.fs.unit.c.flow_vol[0]))
        assert (pytest.approx(101325.0, abs=1e-2) ==
                value(sapon.fs.unit.c.pressure[0]))
        assert (pytest.approx(303.15, abs=1e-2) ==
                value(sapon.fs.unit.c.temperature[0]))
        assert (pytest.approx(55388, abs=1e0) ==
                value(sapon.fs.unit.c.conc_mol_comp[0, "H2O"]))
        assert (pytest.approx(100.0, abs=1e-3) ==
                value(sapon.fs.unit.c.conc_mol_comp[0, "NaOH"]))
        assert (pytest.approx(100.0, abs=1e-3) ==
                value(sapon.fs.unit.c.conc_mol_comp[0, "EthylAcetate"]))
        assert (pytest.approx(0.0, abs=1e-3) ==
                value(sapon.fs.unit.c.conc_mol_comp[0, "SodiumAcetate"]))
        assert (pytest.approx(0.0, abs=1e-3) ==
                value(sapon.fs.unit.c.conc_mol_comp[0, "Ethanol"]))

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_conservation(self, sapon):
        assert abs(value(sapon.fs.unit.inlet.flow_vol[0] -
                         sapon.fs.unit.a.flow_vol[0] -
                         sapon.fs.unit.B.flow_vol[0] -
                         sapon.fs.unit.c.flow_vol[0])) <= 1e-6

        assert (abs(value(sapon.fs.unit.inlet.flow_vol[0] *
                          sum(sapon.fs.unit.inlet.conc_mol_comp[0, j]
                              for j in sapon.fs.properties.component_list) -
                          sapon.fs.unit.a.flow_vol[0] *
                          sum(sapon.fs.unit.a.conc_mol_comp[0, j]
                              for j in sapon.fs.properties.component_list) -
                          sapon.fs.unit.B.flow_vol[0] *
                          sum(sapon.fs.unit.B.conc_mol_comp[0, j]
                              for j in sapon.fs.properties.component_list) -
                          sapon.fs.unit.c.flow_vol[0] *
                          sum(sapon.fs.unit.c.conc_mol_comp[0, j]
                              for j in sapon.fs.properties.component_list)))
                <= 1e-5)

        assert abs(value(
                (sapon.fs.unit.inlet.flow_vol[0] *
                 sapon.fs.properties.dens_mol *
                 sapon.fs.properties.cp_mol *
                 (sapon.fs.unit.inlet.temperature[0] -
                    sapon.fs.properties.temperature_ref)) -
                (sapon.fs.unit.a.flow_vol[0] *
                 sapon.fs.properties.dens_mol *
                 sapon.fs.properties.cp_mol *
                 (sapon.fs.unit.a.temperature[0] -
                  sapon.fs.properties.temperature_ref)) -
                (sapon.fs.unit.B.flow_vol[0] *
                 sapon.fs.properties.dens_mol *
                 sapon.fs.properties.cp_mol *
                 (sapon.fs.unit.B.temperature[0] -
                  sapon.fs.properties.temperature_ref)) -
                (sapon.fs.unit.c.flow_vol[0] *
                 sapon.fs.properties.dens_mol *
                 sapon.fs.properties.cp_mol *
                 (sapon.fs.unit.c.temperature[0] -
                  sapon.fs.properties.temperature_ref)))) <= 1e-3

    @pytest.mark.ui
    def test_report(self, sapon):
        sapon.fs.unit.report()


# -----------------------------------------------------------------------------
class TestBTX(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = BTXParameterBlock()

        m.fs.unit = Separator(default={
                "property_package": m.fs.properties,
                "material_balance_type": MaterialBalanceType.componentPhase,
                "split_basis": SplittingType.phaseFlow,
                "ideal_separation": False,
                "has_phase_equilibrium": True})

        return m

    @pytest.mark.build
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "inlet")
        assert len(btx.fs.unit.inlet.vars) == 4
        assert hasattr(btx.fs.unit.inlet, "flow_mol")
        assert hasattr(btx.fs.unit.inlet, "mole_frac")
        assert hasattr(btx.fs.unit.inlet, "temperature")
        assert hasattr(btx.fs.unit.inlet, "pressure")

        assert hasattr(btx.fs.unit, "outlet_1")
        assert len(btx.fs.unit.outlet_1.vars) == 4
        assert hasattr(btx.fs.unit.outlet_1, "flow_mol")
        assert hasattr(btx.fs.unit.outlet_1, "mole_frac")
        assert hasattr(btx.fs.unit.outlet_1, "temperature")
        assert hasattr(btx.fs.unit.outlet_1, "pressure")

        assert hasattr(btx.fs.unit, "outlet_2")
        assert len(btx.fs.unit.outlet_2.vars) == 4
        assert hasattr(btx.fs.unit.outlet_2, "flow_mol")
        assert hasattr(btx.fs.unit.outlet_2, "mole_frac")
        assert hasattr(btx.fs.unit.outlet_2, "temperature")
        assert hasattr(btx.fs.unit.outlet_2, "pressure")

        assert isinstance(btx.fs.unit.split_fraction, Var)

        assert number_variables(btx) == 37
        assert number_total_constraints(btx) == 24
        assert number_unused_variables(btx) == 0

    def test_dof(self, btx):
        btx.fs.unit.inlet.flow_mol[0].fix(5)  # mol/s
        btx.fs.unit.inlet.temperature[0].fix(365)  # K
        btx.fs.unit.inlet.pressure[0].fix(101325)  # Pa
        btx.fs.unit.inlet.mole_frac[0, "benzene"].fix(0.5)
        btx.fs.unit.inlet.mole_frac[0, "toluene"].fix(0.5)

        btx.fs.unit.split_fraction[0, "outlet_1", "Vap"].fix(0.8)
        btx.fs.unit.split_fraction[0, "outlet_2", "Liq"].fix(0.8)

        btx.fs.unit.pprint()
        assert degrees_of_freedom(btx) == 0

#    @pytest.mark.initialize
#    @pytest.mark.solver
#    @pytest.mark.skipif(solver is None, reason="Solver not available")
#    def test_initialize(self, btx):
#        orig_fixed_vars = fixed_variables_set(btx)
#        orig_act_consts = activated_constraints_set(btx)
#
#        btx.fs.unit.initialize(optarg={'tol': 1e-6})
#
#        assert degrees_of_freedom(btx) == 0
#
#        fin_fixed_vars = fixed_variables_set(btx)
#        fin_act_consts = activated_constraints_set(btx)
#
#        assert len(fin_act_consts) == len(orig_act_consts)
#        assert len(fin_fixed_vars) == len(orig_fixed_vars)
#
#        for c in fin_act_consts:
#            assert c in orig_act_consts
#        for v in fin_fixed_vars:
#            assert v in orig_fixed_vars
#
#    @pytest.mark.solver
#    @pytest.mark.skipif(solver is None, reason="Solver not available")
#    def test_solve(self, btx):
#        results = solver.solve(btx)
#
#        # Check for optimal solution
#        assert results.solver.termination_condition == \
#            TerminationCondition.optimal
#        assert results.solver.status == SolverStatus.ok
#
#    @pytest.mark.initialize
#    @pytest.mark.solver
#    @pytest.mark.skipif(solver is None, reason="Solver not available")
#    def test_solution(self, btx):
#        assert (pytest.approx(5, abs=1e-3) ==
#                value(btx.fs.unit.outlet_1.flow_mol[0]))
#        assert (pytest.approx(359.5, abs=1e-1) ==
#                value(btx.fs.unit.outlet_1.temperature[0]))
#        assert (pytest.approx(101325, abs=1e-3) ==
#                value(btx.fs.unit.outlet_1.pressure[0]))
#
#        assert (pytest.approx(1, abs=1e-3) ==
#                value(btx.fs.unit.outlet_2.flow_mol[0]))
#        assert (pytest.approx(329.9, abs=1e-1) ==
#                value(btx.fs.unit.outlet_2.temperature[0]))
#        assert (pytest.approx(101325, abs=1e-3) ==
#                value(btx.fs.unit.outlet_2.pressure[0]))
#
#    @pytest.mark.initialize
#    @pytest.mark.solver
#    @pytest.mark.skipif(solver is None, reason="Solver not available")
#    def test_conservation(self, btx):
#        assert abs(value(btx.fs.unit.inlet_1.flow_mol[0] -
#                         btx.fs.unit.outlet_1.flow_mol[0])) <= 1e-6
#        assert abs(value(btx.fs.unit.inlet_2.flow_mol[0] -
#                         btx.fs.unit.outlet_2.flow_mol[0])) <= 1e-6
#
#        side_1 = value(
#                btx.fs.unit.outlet_1.flow_mol[0] *
#                (btx.fs.unit.side_1.properties_in[0].enth_mol_phase['Liq'] -
#                 btx.fs.unit.side_1.properties_out[0].enth_mol_phase['Liq']))
#        side_2 = value(
#                btx.fs.unit.outlet_2.flow_mol[0] *
#                (btx.fs.unit.side_2.properties_in[0].enth_mol_phase['Liq'] -
#                 btx.fs.unit.side_2.properties_out[0].enth_mol_phase['Liq']))
#        assert abs(side_1 + side_2) <= 1e-6

    @pytest.mark.ui
    def test_report(self, btx):
        btx.fs.unit.report()





#@pytest.mark.build
#def test_build_default():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterTestBlock()
#
#    m.fs.sep = Separator(default={"property_package": m.fs.pp,
#                                  "ideal_separation": False})
#
#    assert isinstance(m.fs.sep.material_splitting_eqn, Constraint)
#    assert len(m.fs.sep.material_splitting_eqn) == 8
#    assert isinstance(m.fs.sep.temperature_equality_eqn, Constraint)
#    assert len(m.fs.sep.temperature_equality_eqn) == 2
#    assert isinstance(m.fs.sep.pressure_equality_eqn, Constraint)
#    assert len(m.fs.sep.pressure_equality_eqn) == 2
#
#    assert isinstance(m.fs.sep.outlet_1, Port)
#    assert isinstance(m.fs.sep.outlet_2, Port)
#    assert isinstance(m.fs.sep.inlet, Port)
#
#
#def test_model_checks():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterTestBlock()
#
#    m.fs.sep = Separator(default={
#            "property_package": m.fs.pp,
#            "ideal_separation": False})
#
#    m.fs.sep.model_check()
#
#    assert m.fs.sep.outlet_1_state[0].check is True
#    assert m.fs.sep.outlet_2_state[0].check is True
#    assert m.fs.sep.mixed_state[0].check is True
#
#
#@pytest.mark.initialize
#def test_initialize():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterTestBlock()
#    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})
#
#    m.fs.sep = Separator(default={
#            "property_package": m.fs.pp,
#            "mixed_state_block": m.fs.sb,
#            "ideal_separation": False,
#            "split_basis": SplittingType.phaseComponentFlow})
#
#    # Change one outlet pressure to check initialization calculations
#    m.fs.sep.outlet_1_state[0].pressure = 8e4
#
#    f = m.fs.sep.initialize(hold_state=True)
#
#    assert m.fs.sep.outlet_1_state[0].init_test is True
#    assert m.fs.sep.outlet_2_state[0].init_test is True
#    assert m.fs.sb[0].init_test is True
#    assert m.fs.sep.outlet_1_state[0].hold_state is False
#    assert m.fs.sep.outlet_2_state[0].hold_state is False
#    assert m.fs.sb[0].hold_state is True
#
#    assert m.fs.sep.outlet_1.component_flow[0, "p1", "c1"].value == 0.5
#    assert m.fs.sep.outlet_1.component_flow[0, "p1", "c2"].value == 0.5
#    assert m.fs.sep.outlet_1.component_flow[0, "p2", "c1"].value == 0.5
#    assert m.fs.sep.outlet_1.component_flow[0, "p2", "c2"].value == 0.5
#    assert m.fs.sep.outlet_1.enthalpy[0, "p1"].value == 2
#    assert m.fs.sep.outlet_1.enthalpy[0, "p2"].value == 2
#    assert m.fs.sep.outlet_1.pressure[0].value == 1e5
#
#    assert m.fs.sep.outlet_2.component_flow[0, "p1", "c1"].value == 0.5
#    assert m.fs.sep.outlet_2.component_flow[0, "p1", "c2"].value == 0.5
#    assert m.fs.sep.outlet_2.component_flow[0, "p2", "c1"].value == 0.5
#    assert m.fs.sep.outlet_2.component_flow[0, "p2", "c2"].value == 0.5
#    assert m.fs.sep.outlet_2.enthalpy[0, "p1"].value == 2
#    assert m.fs.sep.outlet_2.enthalpy[0 ,"p2"].value == 2
#    assert m.fs.sep.outlet_2.pressure[0].value == 1e5
#
#    m.fs.sep.release_state(flags=f)
#
#    assert m.fs.sep.outlet_1_state[0].hold_state is False
#    assert m.fs.sep.outlet_2_state[0].hold_state is False
#    assert m.fs.sb[0].hold_state is False
#
#
#@pytest.mark.initialization
#def test_initialize_inconsistent_keys():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterTestBlock()
#    m.fs.sb = TestStateBlock(m.fs.time, default={"parameters": m.fs.pp})
#
#    m.fs.sep = Separator(default={
#            "property_package": m.fs.pp,
#            "mixed_state_block": m.fs.sb,
#            "ideal_separation": False,
#            "split_basis": SplittingType.phaseFlow})
#
#    # Change one outlet pressure to check initialization calculations
#    m.fs.sep.outlet_1_state[0].pressure = 8e4
#
#    with pytest.raises(KeyError):
#        m.fs.sep.initialize()
#
#
#@pytest.mark.initialization
#@pytest.mark.solver
#@pytest.mark.skipif(solver is None, reason="Solver not available")
#def test_initialize_total_flow():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#
#    m.fs.properties = SaponificationParameterBlock()
#
#    m.fs.sb = Separator(default={
#            "property_package": m.fs.properties,
#            "ideal_separation": False,
#            "split_basis": SplittingType.totalFlow})
#
#    m.fs.sb.inlet.flow_vol.fix(1.0e-03)
#    m.fs.sb.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
#    m.fs.sb.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
#    m.fs.sb.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
#    m.fs.sb.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
#    m.fs.sb.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)
#
#    m.fs.sb.inlet.temperature.fix(303.15)
#    m.fs.sb.inlet.pressure.fix(101325.0)
#
#    m.fs.sb.split_fraction[0, "outlet_1"].fix(0.2)
#
#    assert degrees_of_freedom(m) == 0
#
#    m.fs.sb.initialize(outlvl=5, optarg={'tol': 1e-6})
#
#    assert (pytest.approx(0.2, abs=1e-3) ==
#             m.fs.sb.split_fraction[0, "outlet_1"].value)
#    assert (pytest.approx(0.8, abs=1e-3) ==
#             m.fs.sb.split_fraction[0, "outlet_2"].value)
#
#    assert (pytest.approx(101325.0, abs=1e-2) ==
#            m.fs.sb.outlet_1.pressure[0].value)
#    assert (pytest.approx(303.15, abs=1e-2) ==
#            m.fs.sb.outlet_1.temperature[0].value)
#    assert (pytest.approx(2e-4, abs=1e-6) ==
#            m.fs.sb.outlet_1.flow_vol[0].value)
#    assert (pytest.approx(55388.0, abs=1e-2) ==
#            m.fs.sb.outlet_1.conc_mol_comp[0, "H2O"].value)
#    assert (pytest.approx(100.0, abs=1e-2) ==
#            m.fs.sb.outlet_1.conc_mol_comp[0, "NaOH"].value)
#    assert (pytest.approx(100.0, abs=1e-2) ==
#            m.fs.sb.outlet_1.conc_mol_comp[0, "EthylAcetate"].value)
#    assert (pytest.approx(0.0, abs=1e-2) ==
#            m.fs.sb.outlet_1.conc_mol_comp[0, "SodiumAcetate"].value)
#    assert (pytest.approx(0.0, abs=1e-2) ==
#            m.fs.sb.outlet_1.conc_mol_comp[0, "Ethanol"].value)
#
#    assert (pytest.approx(101325.0, abs=1e-2) ==
#            m.fs.sb.outlet_2.pressure[0].value)
#    assert (pytest.approx(303.15, abs=1e-2) ==
#            m.fs.sb.outlet_2.temperature[0].value)
#    assert (pytest.approx(8e-4, abs=1e-6) ==
#            m.fs.sb.outlet_2.flow_vol[0].value)
#    assert (pytest.approx(55388.0, abs=1e-2) ==
#            m.fs.sb.outlet_2.conc_mol_comp[0, "H2O"].value)
#    assert (pytest.approx(100.0, abs=1e-2) ==
#            m.fs.sb.outlet_2.conc_mol_comp[0, "NaOH"].value)
#    assert (pytest.approx(100.0, abs=1e-2) ==
#            m.fs.sb.outlet_2.conc_mol_comp[0, "EthylAcetate"].value)
#    assert (pytest.approx(0.0, abs=1e-2) ==
#            m.fs.sb.outlet_2.conc_mol_comp[0, "SodiumAcetate"].value)
#    assert (pytest.approx(0.0, abs=1e-2) ==
#            m.fs.sb.outlet_2.conc_mol_comp[0, "Ethanol"].value)
#
#
#@pytest.mark.ui
#def test_report():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#
#    m.fs.properties = SaponificationParameterBlock()
#
#    m.fs.sb = Separator(default={
#            "property_package": m.fs.properties,
#            "ideal_separation": False,
#            "split_basis": SplittingType.totalFlow})
#
#    m.fs.sb.report()
#

## -----------------------------------------------------------------------------
## Testing of ideal splitting methods
#@declare_process_block_class("PhysicalParameterBlock2")
#class _PhysicalParameterBlock2(PhysicalParameterBlock):
#    def build(self):
#        super(_PhysicalParameterBlock2, self).build()
#
#        self.phase_list = Set(initialize=["p1", "p2"])
#        self.component_list = Set(initialize=["c1", "c2"])
#        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])
#
#        self.state_block_class = StateBlock2
#
#    @classmethod
#    def define_metadata(cls, obj):
#        obj.add_default_units({'time': 's',
#                               'length': 'm',
#                               'mass': 'g',
#                               'amount': 'mol',
#                               'temperature': 'K',
#                               'energy': 'J',
#                               'holdup': 'mol'})
#
#
#@declare_process_block_class("StateBlock2", block_class=SBlockBase)
#class StateBlockData2(StateBlockData):
#    CONFIG = ConfigBlock(implicit=True)
#
#    def build(self):
#        super(StateBlockData2, self).build()
#
#        self.phase_list = Set(initialize=["p1", "p2"])
#        self.component_list = Set(initialize=["c1", "c2"])
#        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])
#        self.phase_equilibrium_list = \
#            {"e1": ["c1", ("p1", "p2")],
#             "e2": ["c2", ("p1", "p2")]}
#
#        self.pressure = Var(initialize=1e5)
#        self.flow_mol_phase_comp = Var(self.phase_list,
#                                       self.component_list,
#                                       initialize=1)
#        self.enth_mol_phase = Var(self.phase_list, initialize=2)
#        self.temperature = Var(initialize=5)
#
#    def get_material_flow_terms(b, p, j):
#        return b.flow_mol_phase_comp[p, j]
#
#    def get_enthalpy_flow_terms(b, p):
#        return b.enth_mol_phase[p]
#
#    def define_state_vars(self):
#        return {"component_flow": self.flow_mol_phase_comp,
#                "temperature": self.temperature,
#                "pressure": self.pressure}
#
#    def model_check(self):
#        self.check = True
#
#
##def test_ideal_splitting_methods():
##    m = ConcreteModel()
##    m.fs = FlowsheetBlock(default={"dynamic": False})
##    m.fs.pp = PhysicalParameterBlock2()
##
##    m.fs.sep = SeparatorFrame(default={
##            "property_package": m.fs.pp,
##            "split_basis": SplittingType.phaseComponentFlow,
##            "ideal_split_map": {("p1", "c1"): "outlet_1",
##                                ("p1", "c2"): "outlet_2",
##                                ("p2", "c1"): "outlet_2",
##                                ("p2", "c2"): "outlet_2"}})
##
##    m.fs.sep._get_property_package()
##    m.fs.sep._get_indexing_sets()
##
##    outlet_list = m.fs.sep.create_outlet_list()
##    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
##    mixed_block = m.fs.sep.add_mixed_state_block()
##
##    m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)
##
##    assert isinstance(m.fs.sep.outlet_1, Port)
##    assert isinstance(m.fs.sep.outlet_2, Port)
##
##    assert m.fs.sep.outlet_1[0].component_flow["p1", "c1"].value == 1.0
##    assert m.fs.sep.outlet_1[0].component_flow["p1", "c2"].value == 0.0
##    assert m.fs.sep.outlet_1[0].component_flow["p2", "c1"].value == 0.0
##    assert m.fs.sep.outlet_1[0].component_flow["p2", "c2"].value == 0.0
##    assert m.fs.sep.outlet_1[0].temperature.value == 5
##    assert m.fs.sep.outlet_1[0].pressure.value == 1e5
##
##    assert m.fs.sep.outlet_2[0].component_flow["p1", "c1"].value == 0.0
##    assert m.fs.sep.outlet_2[0].component_flow["p1", "c2"].value == 1.0
##    assert m.fs.sep.outlet_2[0].component_flow["p2", "c1"].value == 1.0
##    assert m.fs.sep.outlet_2[0].component_flow["p2", "c2"].value == 1.0
##    assert m.fs.sep.outlet_2[0].temperature.value == 5
##    assert m.fs.sep.outlet_2[0].pressure.value == 1e5
##
##
##def test_ideal_splitting_methods_key_error():
##    m = ConcreteModel()
##    m.fs = FlowsheetBlock(default={"dynamic": False})
##    m.fs.pp = PhysicalParameterBlock2()
##
##    m.fs.sep = SeparatorFrame(default={
##            "property_package": m.fs.pp,
##            "split_basis": SplittingType.componentFlow,
##            "ideal_split_map": {("p1", "c1"): "outlet_1",
##                                ("p1", "c2"): "outlet_2",
##                                ("p2", "c1"): "outlet_2",
##                                ("p2", "c2"): "outlet_2"}})
##
##    m.fs.sep._get_property_package()
##    m.fs.sep._get_indexing_sets()
##
##    outlet_list = m.fs.sep.create_outlet_list()
##    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
##    mixed_block = m.fs.sep.add_mixed_state_block()
##
##    with pytest.raises(KeyError):
##        m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)
#
#
#def test_ideal_splitting_methods_no_map():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterBlock2()
#
#    m.fs.sep = SeparatorFrame(default={
#            "property_package": m.fs.pp,
#            "split_basis": SplittingType.phaseComponentFlow})
#
#    m.fs.sep._get_property_package()
#    m.fs.sep._get_indexing_sets()
#
#    outlet_list = m.fs.sep.create_outlet_list()
#    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
#    mixed_block = m.fs.sep.add_mixed_state_block()
#
#    with pytest.raises(ConfigurationError):
#        m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)
#
#
#def test_ideal_splitting_methods_totalFlow():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterBlock2()
#
#    m.fs.sep = SeparatorFrame(default={
#            "property_package": m.fs.pp,
#            "split_basis": SplittingType.totalFlow,
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
#    with pytest.raises(ConfigurationError):
#        m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)
#
#
#def test_ideal_splitting_methods_no_ports():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.pp = PhysicalParameterBlock2()
#
#    m.fs.sep = SeparatorFrame(default={
#            "property_package": m.fs.pp,
#            "split_basis": SplittingType.phaseComponentFlow,
#            "ideal_split_map": {("p1", "c1"): "outlet_1",
#                                ("p1", "c2"): "outlet_2",
#                                ("p2", "c1"): "outlet_2",
#                                ("p2", "c2"): "outlet_2"},
#                                "construct_ports": False})
#
#    m.fs.sep._get_property_package()
#    m.fs.sep._get_indexing_sets()
#
#    outlet_list = m.fs.sep.create_outlet_list()
#    outlet_blocks = m.fs.sep.add_outlet_state_blocks(outlet_list)
#    mixed_block = m.fs.sep.add_mixed_state_block()
#
#    with pytest.raises(ConfigurationError):
#        m.fs.sep.partition_outlet_flows(mixed_block, outlet_list)
#
#@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS unavailable")
#@pytest.mark.skipif(solver is None, reason="Solver not available")
#def test_sep_ph_mixed_default():
#    """Test mixed phase form with P-H state vars and phase mass balances
#    """
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.properties = iapws95.Iapws95ParameterBlock(default={})
#    m.fs.sep = Separator(default={
#        "property_package":m.fs.properties,
#        "split_basis":SplittingType.totalFlow,
#        "ideal_separation":False,
#        "num_outlets":2,
#        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})
#
#    m.fs.sep.inlet.enth_mol.fix(24000)
#    m.fs.sep.inlet.flow_mol.fix(100)
#    m.fs.sep.inlet.pressure.fix(101325)
#    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
#    m.fs.sep.initialize()
#    assert degrees_of_freedom(m) == 0
#    solver.solve(m)
#    prop_in = m.fs.sep.mixed_state[0]
#    prop_out = m.fs.sep.outlet_1_state[0]
#    prop_out2 =  m.fs.sep.outlet_2_state[0]
#    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
#    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
#    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
#    assert abs(value(prop_out.flow_mol - 90.0) <= 1e-6)
#    assert abs(value(prop_out2.flow_mol - 10.0) <= 1e-6)
#
#@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS unavailable")
#@pytest.mark.skipif(solver is None, reason="Solver not available")
#def test_sep_ph_lg_default():
#    """Test mixed phase form with P-H state vars and phase mass balances
#    """
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
#        "phase_presentation":iapws95.PhaseType.LG})
#    m.fs.sep = Separator(default={
#        "property_package":m.fs.properties,
#        "split_basis":SplittingType.totalFlow,
#        "ideal_separation":False,
#        "num_outlets":2,
#        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})
#
#    m.fs.sep.inlet.enth_mol.fix(24000)
#    m.fs.sep.inlet.flow_mol.fix(100)
#    m.fs.sep.inlet.pressure.fix(101325)
#    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
#    m.fs.sep.initialize()
#    # by default the separator writes phase mass balances so -2 degress of
#    # freedom is right.  for IAPWS need total mass balance.
#    assert degrees_of_freedom(m) == -2
#
#@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS unavailable")
#@pytest.mark.skipif(solver is None, reason="Solver not available")
#def test_sep_ph_lg_total_mb():
#    """Test mixed phase form with P-H state vars and phase mass balances
#    """
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
#        "phase_presentation":iapws95.PhaseType.LG})
#    m.fs.sep = Separator(default={
#        "property_package":m.fs.properties,
#        "split_basis":SplittingType.totalFlow,
#        "ideal_separation":False,
#        "num_outlets":2,
#        "material_balance_type":MaterialBalanceType.total,
#        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})
#
#    m.fs.sep.inlet.enth_mol.fix(24000)
#    m.fs.sep.inlet.flow_mol.fix(100)
#    m.fs.sep.inlet.pressure.fix(101325)
#    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
#    m.fs.sep.initialize()
#    assert degrees_of_freedom(m) == 0
#    solver.solve(m)
#    prop_in = m.fs.sep.mixed_state[0]
#    prop_out = m.fs.sep.outlet_1_state[0]
#    prop_out2 =  m.fs.sep.outlet_2_state[0]
#    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
#    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
#    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
#    assert abs(value(prop_out.flow_mol - 90.0) <= 1e-6)
#    assert abs(value(prop_out2.flow_mol - 10.0) <= 1e-6)
#
#@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS unavailable")
#@pytest.mark.skipif(solver is None, reason="Solver not available")
#def test_sep_ph_lg_has_phase_eq():
#    """Test mixed phase form with P-H state vars and phase mass balances
#    """
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#    m.fs.properties = iapws95.Iapws95ParameterBlock(default={
#        "phase_presentation":iapws95.PhaseType.LG})
#    m.fs.sep = Separator(default={
#        "has_phase_equilibrium":True,
#        "property_package":m.fs.properties,
#        "split_basis":SplittingType.totalFlow,
#        "ideal_separation":False,
#        "num_outlets":2,
#        "energy_split_basis":EnergySplittingType.equal_molar_enthalpy})
#
#    m.fs.sep.inlet.enth_mol.fix(24000)
#    m.fs.sep.inlet.flow_mol.fix(100)
#    m.fs.sep.inlet.pressure.fix(101325)
#    m.fs.sep.split_fraction[0,"outlet_2"].fix(0.10)
#    m.fs.sep.initialize()
#    assert degrees_of_freedom(m) == 0
#    solver.solve(m)
#    prop_in = m.fs.sep.mixed_state[0]
#    prop_out = m.fs.sep.outlet_1_state[0]
#    prop_out2 =  m.fs.sep.outlet_2_state[0]
#    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
#    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
#    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
#    assert abs(value(prop_out.flow_mol - 90.0) <= 1e-6)
#    assert abs(value(prop_out2.flow_mol - 10.0) <= 1e-6)
