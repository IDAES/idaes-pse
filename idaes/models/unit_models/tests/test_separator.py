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
Tests for Separator unit model.

Author: Andrew Lee
"""
import pytest
import pandas

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Set,
    value,
    Var,
    units as pyunits,
)

from pyomo.network import Port
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    declare_process_block_class,
    MaterialBalanceType,
    StateBlockData,
    StateBlock,
    PhysicalParameterBlock,
    Phase,
    Component,
)
from idaes.models.unit_models.separator import (
    Separator,
    SeparatorData,
    SplittingType,
    EnergySplittingType,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    InitializationError,
)

from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties import iapws95

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
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


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
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(property_package=m.fs.pp)

        return m

    @pytest.mark.unit
    def test_separator_config(self, build):
        assert len(build.fs.sep.config) == 14
        assert build.fs.sep.config.dynamic is False
        assert build.fs.sep.config.has_holdup is False
        assert build.fs.sep.config.property_package == build.fs.pp
        assert isinstance(build.fs.sep.config.property_package_args, ConfigBlock)
        assert len(build.fs.sep.config.property_package_args) == 0
        assert build.fs.sep.config.outlet_list is None
        assert build.fs.sep.config.num_outlets is None
        assert build.fs.sep.config.split_basis == SplittingType.totalFlow
        assert build.fs.sep.config.ideal_separation is False
        assert build.fs.sep.config.ideal_split_map is None
        assert build.fs.sep.config.mixed_state_block is None
        assert build.fs.sep.config.construct_ports is True
        assert (
            build.fs.sep.config.material_balance_type == MaterialBalanceType.useDefault
        )
        assert build.fs.sep.config.has_phase_equilibrium is False

    @pytest.mark.unit
    def test_validate_config_arguments(self, build):
        build.fs.sep.config.has_phase_equilibrium = True
        build.fs.sep.config.ideal_separation = True

        with pytest.raises(ConfigurationError):
            build.fs.sep._validate_config_arguments()

    @pytest.mark.unit
    def test_create_outlet_list_default(self, build):
        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()

        for o in outlet_list:
            assert o in ["outlet_1", "outlet_2"]

    @pytest.mark.unit
    def test_create_outlet_list_outlet_list(self, build):
        build.fs.sep.config.outlet_list = ["foo", "bar"]

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()

        for o in outlet_list:
            assert o in ["foo", "bar"]

    @pytest.mark.unit
    def test_create_outlet_list_num_outlets(self, build):
        build.fs.sep.config.num_outlets = 3

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()

        for o in outlet_list:
            assert o in ["outlet_1", "outlet_2", "outlet_3"]

    @pytest.mark.unit
    def test_create_outlet_list_both_args_consistent(self, build):
        build.fs.sep.config.outlet_list = ["foo", "bar"]
        build.fs.sep.config.num_outlets = 2

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        outlet_list = build.fs.sep.create_outlet_list()

        for o in outlet_list:
            assert o in ["foo", "bar"]

    @pytest.mark.unit
    def test_create_outlet_list_both_args_inconsistent(self, build):
        build.fs.sep.config.outlet_list = ["foo", "bar"]
        build.fs.sep.config.num_outlets = 3

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        with pytest.raises(ConfigurationError):
            build.fs.sep.create_outlet_list()

    @pytest.mark.unit
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

    @pytest.mark.unit
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

    @pytest.mark.unit
    def test_add_mixed_state_block(self, build):
        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        mixed_block = build.fs.sep.add_mixed_state_block()

        assert isinstance(mixed_block, StateBlock)
        assert hasattr(build.fs.sep, "mixed_state")
        assert not build.fs.sep.mixed_state[0].config.has_phase_equilibrium
        assert build.fs.sep.mixed_state[0].config.defined_state
        assert len(build.fs.sep.mixed_state[0].config) == 3

    @pytest.mark.unit
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

    @pytest.mark.unit
    def test_get_mixed_state_block(self, build):
        build.fs.sb = TestStateBlock(build.fs.time, parameters=build.fs.pp)

        build.fs.sep.config.mixed_state_block = build.fs.sb

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        mixed_block = build.fs.sep.get_mixed_state_block()

        assert mixed_block == build.fs.sb

    @pytest.mark.unit
    def test_get_mixed_state_block_none(self, build):
        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        with pytest.raises(BurntToast):
            build.fs.sep.get_mixed_state_block()

    @pytest.mark.unit
    def test_get_mixed_state_block_mismatch(self, build):
        build.fs.sb = TestStateBlock(build.fs.time, parameters=build.fs.pp)

        # Change parameters arg to create mismatch
        build.fs.sb[0].config.parameters = None

        build.fs.sep.config.mixed_state_block = build.fs.sb

        build.fs.sep._get_property_package()
        build.fs.sep._get_indexing_sets()

        with pytest.raises(ConfigurationError):
            build.fs.sep.get_mixed_state_block()


# -----------------------------------------------------------------------------
# Tests of Separator unit model scaling factors
@pytest.mark.unit
class TestBaseScaling(object):
    """Test scaling calculations.  For now they just make sure there are no
    exceptions. This can be expanded in the future.
    """

    @pytest.fixture(scope="function")
    def m(self):
        b = ConcreteModel()
        b.fs = FlowsheetBlock(dynamic=False)
        b.fs.pp = PhysicalParameterTestBlock()
        return b

    def test_no_exception_scaling_calc_external_mixed_state(self, m):
        m.fs.sb = TestStateBlock(m.fs.time, parameters=m.fs.pp)
        m.fs.sep1 = Separator(property_package=m.fs.pp, mixed_state_block=m.fs.sb)
        iscale.calculate_scaling_factors(m)

    def test_no_exception_scaling_calc_internal_mixed_state(self, m):
        m.fs.sep1 = Separator(property_package=m.fs.pp)
        iscale.calculate_scaling_factors(m)


# -----------------------------------------------------------------------------
# Tests of Separator unit model non-ideal construction methods
@pytest.mark.build
class TestSplitConstruction(object):
    @pytest.fixture(scope="function")
    def build(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(property_package=m.fs.pp)

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.outlet_blocks = m.fs.sep.add_outlet_state_blocks(m.outlet_list)
        m.fs.sep.add_mixed_state_block()

        return m

    @pytest.mark.unit
    def test_add_split_fractions_total(self, build):
        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.outlet_idx, Set)
        assert len(build.fs.sep.outlet_idx) == len(build.outlet_list)
        assert isinstance(build.fs.sep.split_fraction, Var)
        assert len(build.fs.sep.split_fraction) == 2
        assert isinstance(build.fs.sep.sum_split_frac, Constraint)
        assert len(build.fs.sep.sum_split_frac) == 1

    @pytest.mark.unit
    def test_add_split_fractions_phase(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.outlet_idx, Set)
        assert len(build.fs.sep.outlet_idx) == len(build.outlet_list)
        assert isinstance(build.fs.sep.split_fraction, Var)
        assert len(build.fs.sep.split_fraction) == 4

        for t in build.fs.time:
            for o in build.fs.sep.outlet_idx:
                for p in build.fs.sep.config.property_package.phase_list:
                    assert build.fs.sep.split_fraction[t, o, p].value == 0.5

        assert isinstance(build.fs.sep.sum_split_frac, Constraint)
        assert len(build.fs.sep.sum_split_frac) == 2

    @pytest.mark.unit
    def test_add_split_fractions_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.outlet_idx, Set)
        assert len(build.fs.sep.outlet_idx) == len(build.outlet_list)
        assert isinstance(build.fs.sep.split_fraction, Var)
        assert len(build.fs.sep.split_fraction) == 4

        for t in build.fs.time:
            for o in build.fs.sep.outlet_idx:
                for j in build.fs.sep.config.property_package.component_list:
                    assert build.fs.sep.split_fraction[t, o, j].value == 0.5

        assert isinstance(build.fs.sep.sum_split_frac, Constraint)
        assert len(build.fs.sep.sum_split_frac) == 2

    @pytest.mark.unit
    def test_add_split_fractions_phase_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.outlet_idx, Set)
        assert len(build.fs.sep.outlet_idx) == len(build.outlet_list)
        assert isinstance(build.fs.sep.split_fraction, Var)
        assert len(build.fs.sep.split_fraction) == 8

        for t in build.fs.time:
            for o in build.fs.sep.outlet_idx:
                for p in build.fs.sep.config.property_package.phase_list:
                    for j in build.fs.sep.config.property_package.component_list:
                        assert 0.5 == build.fs.sep.split_fraction[t, o, p, j].value

        assert isinstance(build.fs.sep.sum_split_frac, Constraint)
        assert len(build.fs.sep.sum_split_frac) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_pc_total_no_equil(self, build):
        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert not hasattr(build.fs.sep, "phase_equilibrium_generation")

    @pytest.mark.unit
    def test_add_material_splitting_constraints_pc_phase_no_equil(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert not hasattr(build.fs.sep, "phase_equilibrium_generation")

    @pytest.mark.unit
    def test_add_material_splitting_constraints_pc_component_no_equil(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert not hasattr(build.fs.sep, "phase_equilibrium_generation")

    @pytest.mark.unit
    def test_add_material_splitting_constraints_pc_phase_component_no_equil(
        self, build
    ):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert not hasattr(build.fs.sep, "phase_equilibrium_generation")

    @pytest.mark.unit
    def test_add_material_splitting_constraints_pc_total_equil(self, build):
        build.fs.sep.config.has_phase_equilibrium = True

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert isinstance(build.fs.sep.phase_equilibrium_generation, Var)
        assert len(build.fs.sep.phase_equilibrium_generation) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_pc_phase_equil(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow
        build.fs.sep.config.has_phase_equilibrium = True

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert isinstance(build.fs.sep.phase_equilibrium_generation, Var)
        assert len(build.fs.sep.phase_equilibrium_generation) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_pc_component_equil(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow
        build.fs.sep.config.has_phase_equilibrium = True

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert isinstance(build.fs.sep.phase_equilibrium_generation, Var)
        assert len(build.fs.sep.phase_equilibrium_generation) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_pc_phase_component_equil(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow
        build.fs.sep.config.has_phase_equilibrium = True

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 8
        assert isinstance(build.fs.sep.phase_equilibrium_generation, Var)
        assert len(build.fs.sep.phase_equilibrium_generation) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_tc_total(self, build):
        build.fs.sep.config.material_balance_type = MaterialBalanceType.componentTotal

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_tc_phase(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.componentTotal

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_tc_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.componentTotal

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_tc_phase_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.componentTotal

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 4

    @pytest.mark.unit
    def test_add_material_splitting_constraints_t_total(self, build):
        build.fs.sep.config.material_balance_type = MaterialBalanceType.total

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 2

    @pytest.mark.unit
    def test_add_material_splitting_constraints_t_phase(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.total

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 2

    @pytest.mark.unit
    def test_add_material_splitting_constraints_t_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.total

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 2

    @pytest.mark.unit
    def test_add_material_splitting_constraints_t_phase_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.total

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.material_splitting_eqn, Constraint)
        assert len(build.fs.sep.material_splitting_eqn) == 2

    @pytest.mark.unit
    def test_add_material_splitting_constraints_te_total(self, build):
        build.fs.sep.config.material_balance_type = MaterialBalanceType.elementTotal

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        with pytest.raises(ConfigurationError):
            build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

    @pytest.mark.unit
    def test_add_material_splitting_constraints_none_total(self, build):
        build.fs.sep.config.material_balance_type = MaterialBalanceType.none

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert not hasattr(build.fs.sep, "material_splitting_eqn")

    @pytest.mark.unit
    def test_add_material_splitting_constraints_none_phase(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.none

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert not hasattr(build.fs.sep, "material_splitting_eqn")

    @pytest.mark.unit
    def test_add_material_splitting_constraints_none_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.componentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.none

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert not hasattr(build.fs.sep, "material_splitting_eqn")

    @pytest.mark.unit
    def test_add_material_splitting_constraints_none_phase_component(self, build):
        build.fs.sep.config.split_basis = SplittingType.phaseComponentFlow
        build.fs.sep.config.material_balance_type = MaterialBalanceType.none

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)

        build.fs.sep.add_material_splitting_constraints(build.fs.sep.mixed_state)

        assert not hasattr(build.fs.sep, "material_splitting_eqn")

    @pytest.mark.unit
    def test_add_energy_splitting_constraints(self, build):
        assert (
            build.fs.sep.config.energy_split_basis
            == EnergySplittingType.equal_temperature
        )

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)
        build.fs.sep.add_energy_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.temperature_equality_eqn, Constraint)
        assert len(build.fs.sep.temperature_equality_eqn) == 2

    @pytest.mark.unit
    def test_add_energy_splitting_constraints_enthalpy(self, build):
        build.fs.sep.config.energy_split_basis = (
            EnergySplittingType.equal_molar_enthalpy
        )
        assert (
            build.fs.sep.config.energy_split_basis
            == EnergySplittingType.equal_molar_enthalpy
        )

        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)
        build.fs.sep.add_energy_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.molar_enthalpy_equality_eqn, Constraint)
        assert len(build.fs.sep.molar_enthalpy_equality_eqn) == 2

    @pytest.mark.unit
    def test_add_momentum_splitting_constraints(self, build):
        build.fs.sep.add_split_fractions(build.outlet_list, build.fs.sep.mixed_state)
        build.fs.sep.add_momentum_splitting_constraints(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.pressure_equality_eqn, Constraint)
        assert len(build.fs.sep.pressure_equality_eqn) == 2

    @pytest.mark.unit
    def test_add_inlet_port_objects(self, build):
        build.fs.sep.add_inlet_port_objects(build.fs.sep.mixed_state)

        assert isinstance(build.fs.sep.inlet, Port)

    @pytest.mark.unit
    def test_add_inlet_port_objects_construct_ports_False(self, build):
        build.fs.sep.config.construct_ports = False

        build.fs.sep.add_inlet_port_objects(build.fs.sep.mixed_state)

        assert hasattr(build.fs.sep, "inlet") is False

    @pytest.mark.unit
    def test_add_outlet_port_objects(self, build):
        build.fs.sep.add_outlet_port_objects(build.outlet_list, build.outlet_blocks)

        assert isinstance(build.fs.sep.outlet_1, Port)
        assert isinstance(build.fs.sep.outlet_2, Port)

    @pytest.mark.unit
    def test_add_outlet_port_objects_construct_ports_False(self, build):
        build.fs.sep.config.construct_ports = False

        build.fs.sep.add_outlet_port_objects(build.outlet_list, build.outlet_blocks)

        assert hasattr(build.fs.sep, "outlet_1") is False
        assert hasattr(build.fs.sep, "outlet_2") is False


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = Separator(
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentPhase,
            split_basis=SplittingType.totalFlow,
            outlet_list=["a", "B", "c"],
            ideal_separation=False,
            has_phase_equilibrium=False,
        )

        m.fs.unit.inlet.flow_vol.fix(1)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(101325.0)

        m.fs.unit.split_fraction[0, "a"].fix(0.3)
        m.fs.unit.split_fraction[0, "B"].fix(0.5)

        return m

    @pytest.mark.build
    @pytest.mark.unit
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

    @pytest.mark.component
    def test_units(self, sapon):
        assert_units_consistent(sapon)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Split Fraction [('a',)]": sapon.fs.unit.split_fraction[0, "a"],
                "Split Fraction [('B',)]": sapon.fs.unit.split_fraction[0, "B"],
                "Split Fraction [('c',)]": sapon.fs.unit.split_fraction[0, "c"],
            }
        }

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
                "Inlet": {
                    "Volumetric Flowrate": 1.00,
                    "Molar Concentration H2O": 55388,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 0,
                    "Molar Concentration Ethanol": 0,
                    "Temperature": 303.15,
                    "Pressure": 1.0132e05,
                },
                "a": {
                    "Volumetric Flowrate": 1.00,
                    "Molar Concentration H2O": 100.00,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 100.00,
                    "Molar Concentration Ethanol": 100.00,
                    "Temperature": 298.15,
                    "Pressure": 1.0132e05,
                },
                "B": {
                    "Volumetric Flowrate": 1.00,
                    "Molar Concentration H2O": 100.00,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 100.00,
                    "Molar Concentration Ethanol": 100.00,
                    "Temperature": 298.15,
                    "Pressure": 1.0132e05,
                },
                "c": {
                    "Volumetric Flowrate": 1.00,
                    "Molar Concentration H2O": 100.00,
                    "Molar Concentration NaOH": 100.00,
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

    @pytest.mark.solver
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
        assert pytest.approx(0.3, abs=1e-5) == value(sapon.fs.unit.a.flow_vol[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(sapon.fs.unit.a.pressure[0])
        assert pytest.approx(303.15, abs=1e-2) == value(sapon.fs.unit.a.temperature[0])
        assert pytest.approx(55388, abs=1e0) == value(
            sapon.fs.unit.a.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, abs=1e-3) == value(
            sapon.fs.unit.a.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, abs=1e-3) == value(
            sapon.fs.unit.a.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.a.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.a.conc_mol_comp[0, "Ethanol"]
        )

        assert pytest.approx(0.5, abs=1e-5) == value(sapon.fs.unit.B.flow_vol[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(sapon.fs.unit.B.pressure[0])
        assert pytest.approx(303.15, abs=1e-2) == value(sapon.fs.unit.B.temperature[0])
        assert pytest.approx(55388, abs=1e0) == value(
            sapon.fs.unit.B.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, abs=1e-3) == value(
            sapon.fs.unit.B.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, abs=1e-3) == value(
            sapon.fs.unit.B.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.B.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.B.conc_mol_comp[0, "Ethanol"]
        )

        assert pytest.approx(0.2, abs=1e-5) == value(sapon.fs.unit.c.flow_vol[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(sapon.fs.unit.c.pressure[0])
        assert pytest.approx(303.15, abs=1e-2) == value(sapon.fs.unit.c.temperature[0])
        assert pytest.approx(55388, abs=1e0) == value(
            sapon.fs.unit.c.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, abs=1e-3) == value(
            sapon.fs.unit.c.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, abs=1e-3) == value(
            sapon.fs.unit.c.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.c.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.c.conc_mol_comp[0, "Ethanol"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0]
                    - sapon.fs.unit.a.flow_vol[0]
                    - sapon.fs.unit.B.flow_vol[0]
                    - sapon.fs.unit.c.flow_vol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.inlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                    - sapon.fs.unit.a.flow_vol[0]
                    * sum(
                        sapon.fs.unit.a.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                    - sapon.fs.unit.B.flow_vol[0]
                    * sum(
                        sapon.fs.unit.B.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                    - sapon.fs.unit.c.flow_vol[0]
                    * sum(
                        sapon.fs.unit.c.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                )
            )
            <= 1e-5
        )

        assert (
            abs(
                value(
                    (
                        sapon.fs.unit.inlet.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.inlet.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    - (
                        sapon.fs.unit.a.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.a.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    - (
                        sapon.fs.unit.B.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.B.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    - (
                        sapon.fs.unit.c.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.c.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                )
            )
            <= 1e-3
        )


# -----------------------------------------------------------------------------
class TestBTXIdeal(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal"
        )

        m.fs.unit = Separator(
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentPhase,
            split_basis=SplittingType.phaseFlow,
            ideal_separation=False,
            has_phase_equilibrium=True,
        )

        m.fs.unit.inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.inlet.temperature[0].fix(368)  # K
        m.fs.unit.inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.split_fraction[0, "outlet_1", "Vap"].fix(0.8)
        m.fs.unit.split_fraction[0, "outlet_2", "Liq"].fix(0.8)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "inlet")
        assert len(btx.fs.unit.inlet.vars) == 4
        assert hasattr(btx.fs.unit.inlet, "flow_mol")
        assert hasattr(btx.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.inlet, "temperature")
        assert hasattr(btx.fs.unit.inlet, "pressure")

        assert hasattr(btx.fs.unit, "outlet_1")
        assert len(btx.fs.unit.outlet_1.vars) == 4
        assert hasattr(btx.fs.unit.outlet_1, "flow_mol")
        assert hasattr(btx.fs.unit.outlet_1, "mole_frac_comp")
        assert hasattr(btx.fs.unit.outlet_1, "temperature")
        assert hasattr(btx.fs.unit.outlet_1, "pressure")

        assert hasattr(btx.fs.unit, "outlet_2")
        assert len(btx.fs.unit.outlet_2.vars) == 4
        assert hasattr(btx.fs.unit.outlet_2, "flow_mol")
        assert hasattr(btx.fs.unit.outlet_2, "mole_frac_comp")
        assert hasattr(btx.fs.unit.outlet_2, "temperature")
        assert hasattr(btx.fs.unit.outlet_2, "pressure")

        assert isinstance(btx.fs.unit.split_fraction, Var)

        assert number_variables(btx) == 59
        assert number_total_constraints(btx) == 52
        assert number_unused_variables(btx) == 0

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Split Fraction [('outlet_1', 'Liq')]": btx.fs.unit.split_fraction[
                    0, "outlet_1", "Liq"
                ],
                "Split Fraction [('outlet_1', 'Vap')]": btx.fs.unit.split_fraction[
                    0, "outlet_1", "Vap"
                ],
                "Split Fraction [('outlet_2', 'Liq')]": btx.fs.unit.split_fraction[
                    0, "outlet_2", "Liq"
                ],
                "Split Fraction [('outlet_2', 'Vap')]": btx.fs.unit.split_fraction[
                    0, "outlet_2", "Vap"
                ],
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "flow_mol": getattr(pyunits.pint_registry, "mole/s"),
                    "mole_frac_comp benzene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "mole_frac_comp toluene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "temperature": getattr(pyunits.pint_registry, "K"),
                    "pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "Inlet": {
                    "flow_mol": 1.00,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 368,
                    "pressure": 101325,
                },
                "outlet_1": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 298.15,
                    "pressure": 101325.0,
                },
                "outlet_2": {
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
    def test_initialiszation(self, btx):
        btx.fs.unit.initialize()

        assert pytest.approx(1, abs=1e-4) == value(btx.fs.unit.mixed_state[0].flow_mol)
        assert pytest.approx(0.604, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].flow_mol_phase["Liq"]
        )
        assert pytest.approx(0.396, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].flow_mol_phase["Vap"]
        )
        assert pytest.approx(368.0, abs=1e-1) == value(
            btx.fs.unit.mixed_state[0].temperature
        )
        assert pytest.approx(101325, abs=1e3) == value(
            btx.fs.unit.mixed_state[0].pressure
        )
        assert pytest.approx(0.412, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].mole_frac_phase_comp["Liq", "benzene"]
        )
        assert pytest.approx(0.588, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].mole_frac_phase_comp["Liq", "toluene"]
        )
        assert pytest.approx(0.634, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].mole_frac_phase_comp["Vap", "benzene"]
        )
        assert pytest.approx(0.366, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].mole_frac_phase_comp["Vap", "toluene"]
        )

        # Also trigger build of phase enthalpy vars.
        btx.fs.unit.mixed_state[0].enth_mol_phase["Vap"] = 0.5
        btx.fs.unit.outlet_1_state[0].enth_mol_phase["Vap"] = 0.5
        btx.fs.unit.outlet_2_state[0].enth_mol_phase["Vap"] = 0.5

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
        assert pytest.approx(0.438, abs=1e-3) == value(btx.fs.unit.outlet_1.flow_mol[0])
        assert pytest.approx(368.0, abs=1e-1) == value(
            btx.fs.unit.outlet_1.temperature[0]
        )
        assert pytest.approx(101325, abs=1e3) == value(btx.fs.unit.outlet_1.pressure[0])
        assert pytest.approx(0.573, abs=1e-3) == value(
            btx.fs.unit.outlet_1.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.427, abs=1e-3) == value(
            btx.fs.unit.outlet_1.mole_frac_comp[0, "toluene"]
        )

        assert pytest.approx(0.562, abs=1e-3) == value(btx.fs.unit.outlet_2.flow_mol[0])
        assert pytest.approx(368.0, abs=1e-1) == value(
            btx.fs.unit.outlet_2.temperature[0]
        )
        assert pytest.approx(101325, abs=1e3) == value(btx.fs.unit.outlet_2.pressure[0])
        assert pytest.approx(0.443, abs=1e-3) == value(
            btx.fs.unit.outlet_2.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.557, abs=1e-3) == value(
            btx.fs.unit.outlet_2.mole_frac_comp[0, "toluene"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.inlet.flow_mol[0]
                    - btx.fs.unit.outlet_1.flow_mol[0]
                    - btx.fs.unit.outlet_2.flow_mol[0]
                )
            )
            <= 1e-5
        )

        assert (
            abs(
                value(
                    btx.fs.unit.inlet.flow_mol[0]
                    * btx.fs.unit.inlet.mole_frac_comp[0, "benzene"]
                    - btx.fs.unit.outlet_1.flow_mol[0]
                    * btx.fs.unit.outlet_1.mole_frac_comp[0, "benzene"]
                    - btx.fs.unit.outlet_2.flow_mol[0]
                    * btx.fs.unit.outlet_2.mole_frac_comp[0, "benzene"]
                )
            )
            <= 1e-5
        )

        assert (
            abs(
                value(
                    btx.fs.unit.inlet.flow_mol[0]
                    * btx.fs.unit.inlet.mole_frac_comp[0, "toluene"]
                    - btx.fs.unit.outlet_1.flow_mol[0]
                    * btx.fs.unit.outlet_1.mole_frac_comp[0, "toluene"]
                    - btx.fs.unit.outlet_2.flow_mol[0]
                    * btx.fs.unit.outlet_2.mole_frac_comp[0, "toluene"]
                )
            )
            <= 1e-5
        )

        assert (
            abs(
                value(
                    btx.fs.unit.mixed_state[0].flow_mol_phase["Vap"]
                    * btx.fs.unit.mixed_state[0].enth_mol_phase["Vap"]
                    + btx.fs.unit.mixed_state[0].flow_mol_phase["Liq"]
                    * btx.fs.unit.mixed_state[0].enth_mol_phase["Liq"]
                    - btx.fs.unit.outlet_1_state[0].flow_mol_phase["Vap"]
                    * btx.fs.unit.outlet_1_state[0].enth_mol_phase["Vap"]
                    - btx.fs.unit.outlet_1_state[0].flow_mol_phase["Liq"]
                    * btx.fs.unit.outlet_1_state[0].enth_mol_phase["Liq"]
                    - btx.fs.unit.outlet_2_state[0].flow_mol_phase["Vap"]
                    * btx.fs.unit.outlet_2_state[0].enth_mol_phase["Vap"]
                    - btx.fs.unit.outlet_2_state[0].flow_mol_phase["Liq"]
                    * btx.fs.unit.outlet_2_state[0].enth_mol_phase["Liq"]
                )
            )
            <= 1e-1
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

        m.fs.unit = Separator(
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentPhase,
            split_basis=SplittingType.componentFlow,
            num_outlets=3,
            ideal_separation=False,
            has_phase_equilibrium=False,
        )

        m.fs.unit.inlet.flow_mol[0].fix(100)
        m.fs.unit.inlet.enth_mol[0].fix(5000)
        m.fs.unit.inlet.pressure[0].fix(101325)

        m.fs.unit.split_fraction[0, "outlet_1", "H2O"].fix(0.4)
        m.fs.unit.split_fraction[0, "outlet_2", "H2O"].fix(0.5)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iapws):
        assert len(iapws.fs.unit.inlet.vars) == 3
        assert hasattr(iapws.fs.unit.inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.inlet, "pressure")

        assert hasattr(iapws.fs.unit, "outlet_1")
        assert len(iapws.fs.unit.outlet_1.vars) == 3
        assert hasattr(iapws.fs.unit.outlet_1, "flow_mol")
        assert hasattr(iapws.fs.unit.outlet_1, "enth_mol")
        assert hasattr(iapws.fs.unit.outlet_1, "pressure")

        assert hasattr(iapws.fs.unit, "outlet_2")
        assert len(iapws.fs.unit.outlet_2.vars) == 3
        assert hasattr(iapws.fs.unit.outlet_2, "flow_mol")
        assert hasattr(iapws.fs.unit.outlet_2, "enth_mol")
        assert hasattr(iapws.fs.unit.outlet_2, "pressure")

        assert hasattr(iapws.fs.unit, "outlet_3")
        assert len(iapws.fs.unit.outlet_3.vars) == 3
        assert hasattr(iapws.fs.unit.outlet_3, "flow_mol")
        assert hasattr(iapws.fs.unit.outlet_3, "enth_mol")
        assert hasattr(iapws.fs.unit.outlet_3, "pressure")

        assert isinstance(iapws.fs.unit.split_fraction, Var)

        assert number_variables(iapws) == 15
        assert number_total_constraints(iapws) == 10
        assert number_unused_variables(iapws) == 0

    @pytest.mark.component
    def test_units(self, iapws):
        assert_units_consistent(iapws)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Split Fraction [('outlet_1', 'H2O')]": iapws.fs.unit.split_fraction[
                    0, "outlet_1", "H2O"
                ],
                "Split Fraction [('outlet_2', 'H2O')]": iapws.fs.unit.split_fraction[
                    0, "outlet_2", "H2O"
                ],
                "Split Fraction [('outlet_3', 'H2O')]": iapws.fs.unit.split_fraction[
                    0, "outlet_3", "H2O"
                ],
            }
        }

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
                "Inlet": {
                    "Molar Flow": 100,
                    "Mass Flow": 1.8015,
                    "T": 339.43,
                    "P": 101325,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 5000,
                },
                "outlet_1": {
                    "Molar Flow": 1,
                    "Mass Flow": 1.8015e-2,
                    "T": 270.4877112932641,
                    "P": 11032305.8275,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 0.01102138712926277,
                },
                "outlet_2": {
                    "Molar Flow": 1,
                    "Mass Flow": 1.8015e-2,
                    "T": 270.4877112932641,
                    "P": 11032305.8275,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 0.01102138712926277,
                },
                "outlet_3": {
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
    def test_initialiszation(self, iapws):
        iapws.fs.unit.initialize()

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
        assert pytest.approx(40, abs=1e-3) == value(iapws.fs.unit.outlet_1.flow_mol[0])
        assert pytest.approx(50, abs=1e-3) == value(iapws.fs.unit.outlet_2.flow_mol[0])
        assert pytest.approx(10, abs=1e-3) == value(iapws.fs.unit.outlet_3.flow_mol[0])

        assert pytest.approx(5000, abs=1e0) == value(iapws.fs.unit.outlet_1.enth_mol[0])
        assert pytest.approx(5000, abs=1e0) == value(iapws.fs.unit.outlet_2.enth_mol[0])
        assert pytest.approx(5000, abs=1e0) == value(iapws.fs.unit.outlet_3.enth_mol[0])

        assert pytest.approx(101325, abs=1e2) == value(
            iapws.fs.unit.outlet_1.pressure[0]
        )
        assert pytest.approx(101325, abs=1e2) == value(
            iapws.fs.unit.outlet_2.pressure[0]
        )
        assert pytest.approx(101325, abs=1e2) == value(
            iapws.fs.unit.outlet_3.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert (
            abs(
                value(
                    iapws.fs.unit.inlet.flow_mol[0]
                    - iapws.fs.unit.outlet_1.flow_mol[0]
                    - iapws.fs.unit.outlet_2.flow_mol[0]
                    - iapws.fs.unit.outlet_3.flow_mol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    iapws.fs.unit.inlet.flow_mol[0] * iapws.fs.unit.inlet.enth_mol[0]
                    - iapws.fs.unit.outlet_1.flow_mol[0]
                    * iapws.fs.unit.outlet_1.enth_mol[0]
                    - iapws.fs.unit.outlet_2.flow_mol[0]
                    * iapws.fs.unit.outlet_2.enth_mol[0]
                    - iapws.fs.unit.outlet_3.flow_mol[0]
                    * iapws.fs.unit.outlet_3.enth_mol[0]
                )
            )
            <= 1e-2
        )


# -----------------------------------------------------------------------------
# Define some generic Property Block classes for testing ideal separations
@declare_process_block_class("IdealTestBlock")
class _IdealParameterBlock(PhysicalParameterBlock):
    def build(self):
        super(_IdealParameterBlock, self).build()

        self.p1 = Phase()
        self.p2 = Phase()
        self.c1 = Component()
        self.c2 = Component()

        self._phase_component_set = Set(
            initialize=[("p1", "c1"), ("p1", "c2"), ("p2", "c1"), ("p2", "c2")]
        )

        self._state_block_class = IdealStateBlock

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


@declare_process_block_class("IdealStateBlock", block_class=StateBlock)
class IdealTestBlockData(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(IdealTestBlockData, self).build()

        # Add an attribute to allow us to change the state variable definition
        self._state_var_switch = 1

        self.flow_mol_phase_comp = Var(
            self.params.phase_list, self.params.component_list, initialize=2
        )
        self.flow_mol_phase = Var(self.params.phase_list, initialize=2)
        self.flow_mol_comp = Var(self.params.component_list, initialize=2)
        self.flow_mol = Var(initialize=2)

        self.pressure = Var(initialize=1e5)
        self.temperature = Var(initialize=300)

        self.mole_frac_comp = Var(self.params.component_list, initialize=0.5)
        self.mole_frac_phase_comp = Var(
            self.params.phase_list, self.params.component_list, initialize=0.5
        )

        self.test_var = Var(initialize=1)
        self.test_var_comp = Var(self.params.component_list, initialize=1)
        self.test_var_phase = Var(self.params.phase_list, initialize=1)
        self.test_var_phase_comp = Var(
            self.params.phase_list, self.params.component_list, initialize=1
        )

        # Set some values to make sure partitioning is correct
        self.flow_mol_phase_comp["p1", "c1"] = 1
        self.flow_mol_phase_comp["p1", "c2"] = 2
        self.flow_mol_phase_comp["p2", "c1"] = 3
        self.flow_mol_phase_comp["p2", "c2"] = 4

        self.flow_mol_phase["p1"] = 5
        self.flow_mol_phase["p2"] = 6

        self.flow_mol_comp["c1"] = 7
        self.flow_mol_comp["c2"] = 8

        self.flow_mol = 9

        self.mole_frac_phase_comp["p1", "c1"] = 0.9
        self.mole_frac_phase_comp["p1", "c2"] = 0.7
        self.mole_frac_phase_comp["p2", "c1"] = 0.5
        self.mole_frac_phase_comp["p2", "c2"] = 0.3

        self.test_var_comp["c1"] = 2000
        self.test_var_comp["c2"] = 3000

        self.test_var_phase["p1"] = 4000
        self.test_var_phase["p2"] = 5000

        self.test_var_phase_comp["p1", "c1"] = 6000
        self.test_var_phase_comp["p1", "c2"] = 7000
        self.test_var_phase_comp["p2", "c1"] = 8000
        self.test_var_phase_comp["p2", "c2"] = 9000

    def define_state_vars(self):
        if self._state_var_switch == 1:
            return {"mole_frac_comp": self.mole_frac_comp}
        elif self._state_var_switch == 2:
            return {"mole_frac_phase_comp": self.mole_frac_phase_comp}
        elif self._state_var_switch == 3:
            return {"flow_mol_phase_comp": self.flow_mol_phase_comp}
        elif self._state_var_switch == 4:
            return {"flow_mol_phase": self.flow_mol_phase}
        elif self._state_var_switch == 5:
            return {"flow_mol_comp": self.flow_mol_comp}
        elif self._state_var_switch == 6:
            return {"temperature": self.temperature, "pressure": self.pressure}
        elif self._state_var_switch == 7:
            return {"test_var": self.test_var}


# -----------------------------------------------------------------------------
# Tests of Separator unit model ideal construction methods
@pytest.mark.build
class TestIdealConstruction(object):
    @pytest.mark.unit
    def test_phase_component(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert isinstance(m.fs.sep.outlet_1, Port)
        assert isinstance(m.fs.sep.outlet_2, Port)
        assert isinstance(m.fs.sep.outlet_3, Port)
        assert isinstance(m.fs.sep.outlet_4, Port)

        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p1", "c1"]) == 2.0
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p2", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.temperature[0]) == 300
        assert value(m.fs.sep.outlet_1.pressure[0]) == 1e5

        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p1", "c2"]) == 2.0
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p2", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_2.temperature[0]) == 300
        assert value(m.fs.sep.outlet_2.pressure[0]) == 1e5

        assert value(m.fs.sep.outlet_3.component_flow_phase[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_3.component_flow_phase[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_3.component_flow_phase[0, "p2", "c1"]) == 2.0
        assert value(m.fs.sep.outlet_3.component_flow_phase[0, "p2", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_3.temperature[0]) == 300
        assert value(m.fs.sep.outlet_3.pressure[0]) == 1e5

        assert value(m.fs.sep.outlet_4.component_flow_phase[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.component_flow_phase[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_4.component_flow_phase[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.component_flow_phase[0, "p2", "c2"]) == 2.0
        assert value(m.fs.sep.outlet_4.temperature[0]) == 300
        assert value(m.fs.sep.outlet_4.pressure[0]) == 1e5

    @pytest.mark.unit
    def test_phase(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert isinstance(m.fs.sep.outlet_1, Port)
        assert isinstance(m.fs.sep.outlet_2, Port)

        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p1", "c1"]) == 2.0
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p1", "c2"]) == 2.0
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p2", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.temperature[0]) == 300
        assert value(m.fs.sep.outlet_1.pressure[0]) == 1e5

        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p2", "c1"]) == 2.0
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p2", "c2"]) == 2.0
        assert value(m.fs.sep.outlet_2.temperature[0]) == 300
        assert value(m.fs.sep.outlet_2.pressure[0]) == 1e5

    @pytest.mark.unit
    def test_component(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert isinstance(m.fs.sep.outlet_1, Port)
        assert isinstance(m.fs.sep.outlet_2, Port)

        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p1", "c1"]) == 2.0
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p2", "c1"]) == 2.0
        assert value(m.fs.sep.outlet_1.component_flow_phase[0, "p2", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.temperature[0]) == 300
        assert value(m.fs.sep.outlet_1.pressure[0]) == 1e5

        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p1", "c2"]) == 2.0
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.component_flow_phase[0, "p2", "c2"]) == 2.0
        assert value(m.fs.sep.outlet_2.temperature[0]) == 300
        assert value(m.fs.sep.outlet_2.pressure[0]) == 1e5

    @pytest.mark.unit
    def test_ideal_w_no_ports(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
            construct_ports=False,
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        with pytest.raises(ConfigurationError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_ideal_w_total_flow(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.totalFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        with pytest.raises(ConfigurationError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_ideal_w_no_split_map(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.totalFlow,
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        with pytest.raises(ConfigurationError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_phase_component_mismatch(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={("p1", "c1"): "outlet_1", ("p1", "c2"): "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        with pytest.raises(ConfigurationError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_component_mismatch(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        with pytest.raises(ConfigurationError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_phase_mismatch(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        with pytest.raises(ConfigurationError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_split_map_mismatch(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = PhysicalParameterTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=1,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        with pytest.raises(ConfigurationError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_mole_frac_w_component_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.mole_frac_comp[0, "c1"]) == 1
        assert value(m.fs.sep.outlet_1.mole_frac_comp[0, "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_2.mole_frac_comp[0, "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.mole_frac_comp[0, "c2"]) == 1

    @pytest.mark.unit
    def test_mole_frac_w_phase_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.mole_frac_comp[0, "c1"]) == 0.9
        assert value(m.fs.sep.outlet_1.mole_frac_comp[0, "c2"]) == 0.7
        assert value(m.fs.sep.outlet_2.mole_frac_comp[0, "c1"]) == 0.5
        assert value(m.fs.sep.outlet_2.mole_frac_comp[0, "c2"]) == 0.3

    @pytest.mark.unit
    def test_mole_frac_w_phase_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.mole_frac_comp[0, "c1"]) == 1
        assert value(m.fs.sep.outlet_1.mole_frac_comp[0, "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_2.mole_frac_comp[0, "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.mole_frac_comp[0, "c2"]) == 1
        assert value(m.fs.sep.outlet_3.mole_frac_comp[0, "c1"]) == 1
        assert value(m.fs.sep.outlet_3.mole_frac_comp[0, "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_4.mole_frac_comp[0, "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.mole_frac_comp[0, "c2"]) == 1

    @pytest.mark.unit
    def test_mole_frac_w_phase_split_no_fallback(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        # Delete mole_frac_phase_comp so that the fallback should fail
        m.fs.sep.mixed_state[0].del_component(
            m.fs.sep.mixed_state[0].mole_frac_phase_comp
        )

        with pytest.raises(AttributeError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_mole_frac_phase_w_component_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 2

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p1", "c1"]) == 1
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p2", "c1"]) == 1
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p2", "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p1", "c2"]) == 1
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p2", "c2"]) == 1

    @pytest.mark.unit
    def test_mole_frac_phase_w_phase_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 2

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p1", "c1"]) == 0.9
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p1", "c2"]) == 0.7
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p2", "c1"]) == 0.5
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p2", "c2"]) == 0.3

        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p1", "c1"]) == 0.9
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p1", "c2"]) == 0.7
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p2", "c1"]) == 0.5
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p2", "c2"]) == 0.3

    @pytest.mark.unit
    def test_mole_frac_phase_w_phase_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 2

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p1", "c1"]) == 1
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p2", "c1"]) == 1
        assert value(m.fs.sep.outlet_1.mole_frac_phase_comp[0, "p2", "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p1", "c2"]) == 1
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.mole_frac_phase_comp[0, "p2", "c2"]) == 1

        assert value(m.fs.sep.outlet_3.mole_frac_phase_comp[0, "p1", "c1"]) == 1
        assert value(m.fs.sep.outlet_3.mole_frac_phase_comp[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_3.mole_frac_phase_comp[0, "p2", "c1"]) == 1
        assert value(m.fs.sep.outlet_3.mole_frac_phase_comp[0, "p2", "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_4.mole_frac_phase_comp[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.mole_frac_phase_comp[0, "p1", "c2"]) == 1
        assert value(m.fs.sep.outlet_4.mole_frac_phase_comp[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.mole_frac_phase_comp[0, "p2", "c2"]) == 1

    @pytest.mark.unit
    def test_flow_phase_comp_w_phase_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 3

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p1", "c1"]) == 1
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p2", "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p1", "c2"]) == 2
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p2", "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_3.flow_mol_phase_comp[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_3.flow_mol_phase_comp[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_3.flow_mol_phase_comp[0, "p2", "c1"]) == 3
        assert value(m.fs.sep.outlet_3.flow_mol_phase_comp[0, "p2", "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_4.flow_mol_phase_comp[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.flow_mol_phase_comp[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_4.flow_mol_phase_comp[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.flow_mol_phase_comp[0, "p2", "c2"]) == 4

    @pytest.mark.unit
    def test_flow_phase_comp_w_phase_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 3

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p1", "c1"]) == 1
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p1", "c2"]) == 2
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p2", "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p2", "c1"]) == 3
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p2", "c2"]) == 4

    @pytest.mark.unit
    def test_flow_phase_comp_w_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 3

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p1", "c1"]) == 1
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p1", "c2"]) == 1e-8
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p2", "c1"]) == 3
        assert value(m.fs.sep.outlet_1.flow_mol_phase_comp[0, "p2", "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p1", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p1", "c2"]) == 2
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p2", "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_phase_comp[0, "p2", "c2"]) == 4

    @pytest.mark.unit
    def test_flow_phase_w_phase_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 4

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_phase[0, "p1"]) == 1
        assert value(m.fs.sep.outlet_1.flow_mol_phase[0, "p2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.flow_mol_phase[0, "p1"]) == 2
        assert value(m.fs.sep.outlet_2.flow_mol_phase[0, "p2"]) == 1e-8

        assert value(m.fs.sep.outlet_3.flow_mol_phase[0, "p1"]) == 1e-8
        assert value(m.fs.sep.outlet_3.flow_mol_phase[0, "p2"]) == 3

        assert value(m.fs.sep.outlet_4.flow_mol_phase[0, "p1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.flow_mol_phase[0, "p2"]) == 4

    @pytest.mark.unit
    def test_flow_phase_w_phase_comp_split_no_fallback(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 4
        m.fs.sep.mixed_state[0].del_component(
            m.fs.sep.mixed_state[0].flow_mol_phase_comp
        )

        with pytest.raises(AttributeError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_flow_phase_w_phase_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 4

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_phase[0, "p1"]) == 5
        assert value(m.fs.sep.outlet_1.flow_mol_phase[0, "p2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.flow_mol_phase[0, "p1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_phase[0, "p2"]) == 6

    @pytest.mark.unit
    def test_flow_phase_w_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 4

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_phase[0, "p1"]) == 1
        assert value(m.fs.sep.outlet_1.flow_mol_phase[0, "p2"]) == 3

        assert value(m.fs.sep.outlet_2.flow_mol_phase[0, "p1"]) == 2
        assert value(m.fs.sep.outlet_2.flow_mol_phase[0, "p2"]) == 4

    @pytest.mark.unit
    def test_flow_phase_w_comp_split_no_fallback(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 4
        m.fs.sep.mixed_state[0].del_component(
            m.fs.sep.mixed_state[0].flow_mol_phase_comp
        )

        with pytest.raises(AttributeError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_flow_comp_w_phase_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 5

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_comp[0, "c1"]) == 1
        assert value(m.fs.sep.outlet_1.flow_mol_comp[0, "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.flow_mol_comp[0, "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_comp[0, "c2"]) == 2

        assert value(m.fs.sep.outlet_3.flow_mol_comp[0, "c1"]) == 3
        assert value(m.fs.sep.outlet_3.flow_mol_comp[0, "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_4.flow_mol_comp[0, "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_4.flow_mol_comp[0, "c2"]) == 4

    @pytest.mark.unit
    def test_flow_comp_w_phase_comp_split_no_fallback(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 5
        m.fs.sep.mixed_state[0].del_component(
            m.fs.sep.mixed_state[0].flow_mol_phase_comp
        )

        with pytest.raises(AttributeError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_flow_comp_w_phase_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 5

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_comp[0, "c1"]) == 1
        assert value(m.fs.sep.outlet_1.flow_mol_comp[0, "c2"]) == 2

        assert value(m.fs.sep.outlet_2.flow_mol_comp[0, "c1"]) == 3
        assert value(m.fs.sep.outlet_2.flow_mol_comp[0, "c2"]) == 4

    @pytest.mark.unit
    def test_flow_comp_w_phase_split_no_fallback(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 5
        m.fs.sep.mixed_state[0].del_component(
            m.fs.sep.mixed_state[0].flow_mol_phase_comp
        )

        with pytest.raises(AttributeError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_flow_comp_w_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 5

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.flow_mol_comp[0, "c1"]) == 7
        assert value(m.fs.sep.outlet_1.flow_mol_comp[0, "c2"]) == 1e-8

        assert value(m.fs.sep.outlet_2.flow_mol_comp[0, "c1"]) == 1e-8
        assert value(m.fs.sep.outlet_2.flow_mol_comp[0, "c2"]) == 8

    @pytest.mark.unit
    def test_t_p(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 6

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.temperature[0]) == 300
        assert value(m.fs.sep.outlet_1.pressure[0]) == 1e5

        assert value(m.fs.sep.outlet_2.temperature[0]) == 300
        assert value(m.fs.sep.outlet_2.pressure[0]) == 1e5

    @pytest.mark.unit
    def test_general_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 7

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.test_var[0]) == 2000

        assert value(m.fs.sep.outlet_2.test_var[0]) == 3000

    @pytest.mark.unit
    def test_general_comp_split_fallback(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 7
        m.fs.sep.mixed_state[0].del_component(m.fs.sep.mixed_state[0].test_var_comp)

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.test_var[0]) == 14000

        assert value(m.fs.sep.outlet_2.test_var[0]) == 16000

    @pytest.mark.unit
    def test_general_comp_split_fallback_fail(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.componentFlow,
            ideal_split_map={"c1": "outlet_1", "c2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 7
        m.fs.sep.mixed_state[0].del_component(m.fs.sep.mixed_state[0].test_var_comp)
        m.fs.sep.mixed_state[0].del_component(
            m.fs.sep.mixed_state[0].test_var_phase_comp
        )

        with pytest.raises(AttributeError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_general_phase_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 7

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.test_var[0]) == 4000

        assert value(m.fs.sep.outlet_2.test_var[0]) == 5000

    @pytest.mark.unit
    def test_general_phase_split_fallback(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 7
        m.fs.sep.mixed_state[0].del_component(m.fs.sep.mixed_state[0].test_var_phase)

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.test_var[0]) == 13000

        assert value(m.fs.sep.outlet_2.test_var[0]) == 17000

    @pytest.mark.unit
    def test_general_phase_split_fallback_fail(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=2,
            ideal_separation=True,
            split_basis=SplittingType.phaseFlow,
            ideal_split_map={"p1": "outlet_1", "p2": "outlet_2"},
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 7
        m.fs.sep.mixed_state[0].del_component(m.fs.sep.mixed_state[0].test_var_phase)
        m.fs.sep.mixed_state[0].del_component(
            m.fs.sep.mixed_state[0].test_var_phase_comp
        )

        with pytest.raises(AttributeError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

    @pytest.mark.unit
    def test_general_phase_comp_split(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 7

        m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)

        assert value(m.fs.sep.outlet_1.test_var[0]) == 6000
        assert value(m.fs.sep.outlet_2.test_var[0]) == 7000
        assert value(m.fs.sep.outlet_3.test_var[0]) == 8000
        assert value(m.fs.sep.outlet_4.test_var[0]) == 9000

    @pytest.mark.unit
    def test_general_phase_comp_split_fallback_fail(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.pp = IdealTestBlock()

        m.fs.sep = SeparatorFrame(
            property_package=m.fs.pp,
            num_outlets=4,
            ideal_separation=True,
            split_basis=SplittingType.phaseComponentFlow,
            ideal_split_map={
                ("p1", "c1"): "outlet_1",
                ("p1", "c2"): "outlet_2",
                ("p2", "c1"): "outlet_3",
                ("p2", "c2"): "outlet_4",
            },
        )

        m.fs.sep._get_property_package()
        m.fs.sep._get_indexing_sets()

        m.outlet_list = m.fs.sep.create_outlet_list()
        m.fs.sep.add_mixed_state_block()

        m.fs.sep.mixed_state[0]._state_var_switch = 7
        m.fs.sep.mixed_state[0].del_component(
            m.fs.sep.mixed_state[0].test_var_phase_comp
        )

        with pytest.raises(AttributeError):
            m.fs.sep.partition_outlet_flows(m.fs.sep.mixed_state, m.outlet_list)


# -----------------------------------------------------------------------------
class TestBTX_Ideal(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock()

        m.fs.unit = Separator(
            property_package=m.fs.properties,
            material_balance_type=MaterialBalanceType.componentPhase,
            split_basis=SplittingType.phaseFlow,
            ideal_separation=True,
            ideal_split_map={"Vap": "outlet_1", "Liq": "outlet_2"},
            has_phase_equilibrium=False,
        )

        m.fs.unit.inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.inlet.temperature[0].fix(368)  # K
        m.fs.unit.inlet.pressure[0].fix(101325)  # Pa

        m.fs.unit.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "inlet")
        assert len(btx.fs.unit.inlet.vars) == 4
        assert hasattr(btx.fs.unit.inlet, "flow_mol")
        assert hasattr(btx.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.inlet, "temperature")
        assert hasattr(btx.fs.unit.inlet, "pressure")

        assert hasattr(btx.fs.unit, "outlet_1")
        assert len(btx.fs.unit.outlet_1.vars) == 4
        assert hasattr(btx.fs.unit.outlet_1, "flow_mol")
        assert hasattr(btx.fs.unit.outlet_1, "mole_frac_comp")
        assert hasattr(btx.fs.unit.outlet_1, "temperature")
        assert hasattr(btx.fs.unit.outlet_1, "pressure")

        assert hasattr(btx.fs.unit, "outlet_2")
        assert len(btx.fs.unit.outlet_2.vars) == 4
        assert hasattr(btx.fs.unit.outlet_2, "flow_mol")
        assert hasattr(btx.fs.unit.outlet_2, "mole_frac_comp")
        assert hasattr(btx.fs.unit.outlet_2, "temperature")
        assert hasattr(btx.fs.unit.outlet_2, "pressure")

        assert number_variables(btx) == 17
        assert number_total_constraints(btx) == 12
        assert number_unused_variables(btx) == 0

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict is None

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "flow_mol": getattr(pyunits.pint_registry, "mole/s"),
                    "mole_frac_comp benzene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "mole_frac_comp toluene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "temperature": getattr(pyunits.pint_registry, "K"),
                    "pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "inlet": {
                    "flow_mol": 1.00,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 368,
                    "pressure": 101325,
                },
                "outlet_1": {
                    "flow_mol": 0.5,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 368,
                    "pressure": 101325.0,
                },
                "outlet_2": {
                    "flow_mol": 0.5,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 368,
                    "pressure": 101325.0,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialiszation(self, btx):
        btx.fs.unit.initialize()

        assert pytest.approx(1, abs=1e-4) == value(btx.fs.unit.mixed_state[0].flow_mol)
        assert pytest.approx(0.604, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].flow_mol_phase["Liq"]
        )
        assert pytest.approx(0.396, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].flow_mol_phase["Vap"]
        )
        assert pytest.approx(368.0, abs=1e-1) == value(
            btx.fs.unit.mixed_state[0].temperature
        )
        assert pytest.approx(101325, abs=1e3) == value(
            btx.fs.unit.mixed_state[0].pressure
        )
        assert pytest.approx(0.412, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].mole_frac_phase_comp["Liq", "benzene"]
        )
        assert pytest.approx(0.588, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].mole_frac_phase_comp["Liq", "toluene"]
        )
        assert pytest.approx(0.634, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].mole_frac_phase_comp["Vap", "benzene"]
        )
        assert pytest.approx(0.366, abs=1e-3) == value(
            btx.fs.unit.mixed_state[0].mole_frac_phase_comp["Vap", "toluene"]
        )

        assert degrees_of_freedom(btx) == 0

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
        assert pytest.approx(0.396, abs=1e-3) == value(btx.fs.unit.outlet_1.flow_mol[0])
        assert pytest.approx(368.0, abs=1e-1) == value(
            btx.fs.unit.outlet_1.temperature[0]
        )
        assert pytest.approx(101325, abs=1e3) == value(btx.fs.unit.outlet_1.pressure[0])
        assert pytest.approx(0.634, abs=1e-3) == value(
            btx.fs.unit.outlet_1.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.366, abs=1e-3) == value(
            btx.fs.unit.outlet_1.mole_frac_comp[0, "toluene"]
        )

        assert pytest.approx(0.604, abs=1e-3) == value(btx.fs.unit.outlet_2.flow_mol[0])
        assert pytest.approx(368.0, abs=1e-1) == value(
            btx.fs.unit.outlet_2.temperature[0]
        )
        assert pytest.approx(101325, abs=1e3) == value(btx.fs.unit.outlet_2.pressure[0])
        assert pytest.approx(0.412, abs=1e-3) == value(
            btx.fs.unit.outlet_2.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.588, abs=1e-3) == value(
            btx.fs.unit.outlet_2.mole_frac_comp[0, "toluene"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.inlet.flow_mol[0]
                    - btx.fs.unit.outlet_1.flow_mol[0]
                    - btx.fs.unit.outlet_2.flow_mol[0]
                )
            )
            <= 1e-5
        )

        assert (
            abs(
                value(
                    btx.fs.unit.inlet.flow_mol[0]
                    * btx.fs.unit.inlet.mole_frac_comp[0, "benzene"]
                    - btx.fs.unit.outlet_1.flow_mol[0]
                    * btx.fs.unit.outlet_1.mole_frac_comp[0, "benzene"]
                    - btx.fs.unit.outlet_2.flow_mol[0]
                    * btx.fs.unit.outlet_2.mole_frac_comp[0, "benzene"]
                )
            )
            <= 1e-5
        )

        assert (
            abs(
                value(
                    btx.fs.unit.inlet.flow_mol[0]
                    * btx.fs.unit.inlet.mole_frac_comp[0, "toluene"]
                    - btx.fs.unit.outlet_1.flow_mol[0]
                    * btx.fs.unit.outlet_1.mole_frac_comp[0, "toluene"]
                    - btx.fs.unit.outlet_2.flow_mol[0]
                    * btx.fs.unit.outlet_2.mole_frac_comp[0, "toluene"]
                )
            )
            <= 1e-5
        )

        # Assume energy conservation is covered by control volume tests


@pytest.mark.unit
def test_initialization_error():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.sep = Separator(property_package=m.fs.pp)

    m.fs.sep.outlet_1_state[0].material_flow_mol.fix(10)
    m.fs.sep.outlet_2_state[0].material_flow_mol.fix(10)
    m.fs.sep.mixed_state[0].material_flow_mol.fix(100)

    m.fs.sep.split_fraction.fix()

    with pytest.raises(InitializationError):
        m.fs.sep.initialize()
