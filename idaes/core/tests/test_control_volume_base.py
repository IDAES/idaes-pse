##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
# 
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for ControlVolumeBase.

Author: Andrew Lee
"""
import inspect
import pytest
from pyomo.environ import ConcreteModel, Block, Set
from pyomo.common.config import ConfigBlock, ConfigValue
from idaes.core import (ControlVolumeBase, CONFIG_Base,
                        MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType, FlowDirection,
                        declare_process_block_class,
                        FlowsheetBlockData, UnitBlockData, useDefault,
                        PhysicalParameterBase, ReactionParameterBase)
from idaes.core.util.exceptions import (ConfigurationError, DynamicError,
                                        PropertyPackageError)


# -----------------------------------------------------------------------------
# Test Enumerators for balance type options
def test_material_balance_type():
    assert len(MaterialBalanceType) == 5

    # Test that error is raised when given non-member
    with pytest.raises(AttributeError):
        MaterialBalanceType.foo


def test_energy_balance_type():
    assert len(EnergyBalanceType) == 5

    # Test that error is raised when given non-member
    with pytest.raises(AttributeError):
        EnergyBalanceType.foo


def test_momentum_balance_type():
    assert len(MomentumBalanceType) == 5

    # Test that error is raised when given non-member
    with pytest.raises(AttributeError):
        MomentumBalanceType.foo


def testflow_direction():
    assert len(FlowDirection) == 2

    # Test that error is raised when given non-member
    with pytest.raises(AttributeError):
        FlowDirection.foo


# -----------------------------------------------------------------------------
# Test CONFIG_Base
def test_CONFIG_Base():
    c = CONFIG_Base()

    assert len(c) == 16

    for i in c:
        if i == "dynamic":
            assert c[i] == useDefault
        elif i == "material_balance_type":
            assert c[i] == MaterialBalanceType.componentPhase
        elif i == "energy_balance_type":
            assert c[i] == EnergyBalanceType.enthalpyPhase
        elif i == "momentum_balance_type":
            assert c[i] == MomentumBalanceType.pressureTotal
        elif i == "property_package":
            assert c[i] == useDefault
        elif i == "reaction_package":
            assert c[i] is None
        elif i in ["property_package_args", "reaction_package_args"]:
            assert isinstance(c[i], ConfigBlock)
            assert len(c[i]) == 0
        else:
            assert c[i] is False


def test_CONFIG_Base_validation_general():
    # No config argument takes a string, float/int or list
    c = CONFIG_Base()

    for i in c:
        with pytest.raises(ValueError):
            c[i] = "foo"
        with pytest.raises(ValueError):
            c[i] = 10.0
        with pytest.raises(ValueError):
            c[i] = [1, 2]


def test_CONFIG_Base_true_false():
    # Check arguments that accept True/False as values
    c = CONFIG_Base()

    for i in c:
        if i not in ["material_balance_type", "energy_balance_type",
                     "momentum_balance_type", "property_package",
                     "reaction_package", "property_package_args",
                     "reaction_package_args"]:
            c[i] = True
            c[i] = False


def test_CONFIG_Base_material_balance_type():
    c = CONFIG_Base()

    for i in MaterialBalanceType:
        c["material_balance_type"] = i


def test_CONFIG_Base_energy_balance_type():
    c = CONFIG_Base()

    for i in EnergyBalanceType:
        c["energy_balance_type"] = i


def test_CONFIG_Base_momentum_balance_type():
    c = CONFIG_Base()

    for i in MomentumBalanceType:
        c["momentum_balance_type"] = i


# -----------------------------------------------------------------------------
# Mockup classes for testing
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("Unit")
class _UnitData(UnitBlockData):
    CONFIG = UnitBlockData.CONFIG()
    CONFIG.declare("property_package",
                   ConfigValue(default=None))
    CONFIG.declare("property_package_args",
                   ConfigValue(default={}))

    def build(self):
        super(_UnitData, self).build()


# -----------------------------------------------------------------------------
# Testing ControlVolumeBase
@declare_process_block_class("CVFrame")
class CVFrameData(ControlVolumeBase):
    def build(self):
        super(ControlVolumeBase, self).build()


def test_config_block():
    cv = CVFrame(concrete=True)

    assert len(cv.config) == 6
    assert cv.config.dynamic == useDefault
    assert cv.config.property_package == useDefault
    assert isinstance(cv.config.property_package_args, ConfigBlock)
    assert len(cv.config.property_package_args) == 0
    assert cv.config.reaction_package is None
    assert isinstance(cv.config.reaction_package_args, ConfigBlock)
    assert len(cv.config.reaction_package_args) == 0
    assert cv.config.auto_construct is False


# -----------------------------------------------------------------------------
# Test _setup_dynamics
def test_setup_dynamics_use_parent_value():
    # Test that dynamic = None works correctly
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.u = Unit(default={"dynamic": False})
    m.fs.u.cv = CVFrame()

    m.fs.u.cv._setup_dynamics()

    assert m.fs.u.cv.config.dynamic is False
    assert m.fs.u.cv.time == [0]


def test_setup_dynamics_use_parent_value_fail_no_dynamic():
    # Test that dynamic = None works correctly
    fs = Flowsheet(default={"dynamic": False}, concrete=True)

    # Create a Block (with no dynamic attribute)
    fs.b = Block()
    # Add a time attribute to make sure the correct failure triggers
    fs.b.time = Set(initialize=[0])

    fs.b.cv = CVFrame()

    # _setup_dynamics should return DynamicError
    with pytest.raises(DynamicError):
        fs.b.cv._setup_dynamics()


def test_setup_dynamics_use_parent_value_fail_no_time():
    # Test that methods fails when no time domain in parent
    fs = Flowsheet(default={"dynamic": False}, concrete=True)

    # Mock-up an object with a dynamic flag but no time domain (using CVFrame)
    fs.u = CVFrame()
    # Set the config flag
    fs.u.config.dynamic = False

    # Add a second holdup block
    fs.u.cv = CVFrame()

    # _setup_dynamics should return DynamicError, as fs.u has no time
    with pytest.raises(DynamicError):
        fs.u.cv._setup_dynamics()


# -----------------------------------------------------------------------------
# Test _get_property_package
@declare_process_block_class("PropertyParameterBlock")
class _PropertyParameterBlock(PhysicalParameterBase):
    def build(self):
        super(_PropertyParameterBlock, self).build()

        frm = inspect.stack()[1]
        self._package_module = inspect.getmodule(frm[0])

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])


def test_get_property_package_set():
    m = ConcreteModel()
    m.pp = PropertyParameterBlock()
    m.cv = CVFrame(default={"property_package": m.pp})
    m.cv._get_property_package()

    assert m.cv._property_module == m.pp._package_module


def test_get_property_package_default_args():
    m = ConcreteModel()
    m.pp = PropertyParameterBlock(
                default={"default_arguments": {"test": "foo"}})
    m.cv = CVFrame(default={"property_package": m.pp})
    m.cv._get_property_package()

    assert m.cv.config.property_package_args["test"] == "foo"


def test_get_reaction_package_module_combine_args():
    # Test that local and default args combine correctly
    m = ConcreteModel()
    m.pp = PropertyParameterBlock(
            default={"default_arguments": {"test1": "foo",
                                           "test2": "bar"}})
    m.cv = CVFrame(default={"property_package": m.pp,
                            "property_package_args": {"test2": "baz",
                                                      "test3": "bar"}})
    m.cv._get_property_package()

    assert m.cv.config.property_package_args["test1"] == "foo"
    assert m.cv.config.property_package_args["test2"] == "baz"
    assert m.cv.config.property_package_args["test3"] == "bar"


# -----------------------------------------------------------------------------
# Test _get_default_prop_pack
def test_get_default_prop_pack_works():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PropertyParameterBlock()
    m.fs.config.default_property_package = m.fs.pp

    m.fs.cv = CVFrame()
    assert m.fs.cv._get_default_prop_pack() == m.fs.pp


# TODO : should test more failure modes
def test_get_default_prop_pack_no_default():
    m = ConcreteModel()
    m.fs = Flowsheet()

    m.fs.cv = CVFrame()
    with pytest.raises(ConfigurationError):
        m.fs.cv._get_default_prop_pack()


def test_get_property_package_call_to_get_default_prop_pack():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PropertyParameterBlock()
    m.fs.config.default_property_package = m.fs.pp

    m.fs.cv = CVFrame()
    m.fs.cv._get_property_package()
    assert m.fs.cv._property_module == m.fs.pp._package_module


# -----------------------------------------------------------------------------
# Test _get_indexing_sets
def test_get_indexing_sets():
    m = ConcreteModel()
    m.pp = PropertyParameterBlock()
    m.cv = CVFrame(default={"property_package": m.pp})
    m.cv._get_property_package()
    m.cv._get_indexing_sets()

    assert hasattr(m.cv, "phase_list")
    assert hasattr(m.cv, "component_list")


def test_get_indexing_sets_missing_phase_list():
    m = ConcreteModel()
    m.pp = PropertyParameterBlock()
    m.pp.del_component(m.pp.phase_list)
    m.cv = CVFrame(default={"property_package": m.pp})
    m.cv._get_property_package()

    with pytest.raises(PropertyPackageError):
        m.cv._get_indexing_sets()


def test_get_indexing_sets_missing_component_list():
    m = ConcreteModel()
    m.pp = PropertyParameterBlock()
    m.pp.del_component(m.pp.component_list)
    m.cv = CVFrame(default={"property_package": m.pp})
    m.cv._get_property_package()

    with pytest.raises(PropertyPackageError):
        m.cv._get_indexing_sets()
    assert hasattr(m.cv, "phase_list")


# -----------------------------------------------------------------------------
# Test _get_reaction_package
def test_get_reaction_package_none():
    m = ConcreteModel()
    m.r = CVFrame()

    m.r._get_reaction_package()

    assert hasattr(m.r, "reaction_module") is False


@declare_process_block_class("ReactionParameterBlock")
class _ReactionParameterBlock(ReactionParameterBase):
    def build(self):
        super(ReactionParameterBase, self).build()

        frm = inspect.stack()[1]
        self._package_module = inspect.getmodule(frm[0])


def test_get_reaction_package_module():
    m = ConcreteModel()
    m.rp = ReactionParameterBlock(
                default={"default_arguments": {"test": "foo"}})
    m.cv = CVFrame(default={"reaction_package": m.rp})

    m.cv._get_reaction_package()

    assert m.cv._reaction_module == m.rp._package_module
    assert m.cv.config.reaction_package_args["test"] == "foo"


def test_get_reaction_package_module_default_args():
    # Test that local and default args combine correctly
    m = ConcreteModel()
    m.rp = ReactionParameterBlock(
            default={"default_arguments": {"test1": "foo",
                                           "test2": "bar"}})
    m.cv = CVFrame(default={"reaction_package": m.rp,
                            "reaction_package_args": {"test2": "baz",
                                                      "test3": "bar"}})

    m.cv._get_reaction_package()

    assert m.cv.config.reaction_package_args["test1"] == "foo"
    assert m.cv.config.reaction_package_args["test2"] == "baz"
    assert m.cv.config.reaction_package_args["test3"] == "bar"


# -----------------------------------------------------------------------------
# Test build and auto_construct methods
def test_build():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PropertyParameterBlock()
    m.fs.cv = CVFrame(default={"property_package": m.fs.pp})

    super(CVFrameData, m.fs.cv).build()


def test_add_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.cv = CVFrame()

    with pytest.raises(NotImplementedError):
        m.fs.cv.add_geometry()


def test_auto_construct():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PropertyParameterBlock()
    m.fs.cv = CVFrame(default={"property_package": m.fs.pp,
                               "auto_construct": True})

    with pytest.raises(NotImplementedError):
        super(CVFrameData, m.fs.cv).build()

# -----------------------------------------------------------------------------
# Test _validate_add_balance_arguments
def test_validate_add_balance_arguments_both_false():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PropertyParameterBlock()
    m.fs.cv = CVFrame(default={"property_package": m.fs.pp})
    super(CVFrameData, m.fs.cv).build()

    d, h = m.fs.cv._validate_add_balance_arguments(dynamic=useDefault,
                                                   has_holdup=False)

    assert d is False
    assert h is False


def test_validate_add_balance_arguments_both_true():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PropertyParameterBlock()
    m.fs.cv = CVFrame(default={"property_package": m.fs.pp})
    super(CVFrameData, m.fs.cv).build()

    d, h = m.fs.cv._validate_add_balance_arguments(dynamic=True,
                                                   has_holdup=True)

    assert d is True
    assert h is True


def test_validate_add_balance_arguments_ss_with_holdup():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PropertyParameterBlock()
    m.fs.cv = CVFrame(default={"property_package": m.fs.pp})
    super(CVFrameData, m.fs.cv).build()

    d, h = m.fs.cv._validate_add_balance_arguments(dynamic=False,
                                                   has_holdup=True)

    assert d is False
    assert h is True


def test_validate_add_balance_arguments_invalid():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PropertyParameterBlock()
    m.fs.cv = CVFrame(default={"property_package": m.fs.pp})
    super(CVFrameData, m.fs.cv).build()

    with pytest.raises(ConfigurationError):
        d, h = m.fs.cv._validate_add_balance_arguments(dynamic=useDefault,
                                                       has_holdup=False)


def test_validate_add_balance_arguments_dynamic_error():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PropertyParameterBlock()
    m.fs.cv = CVFrame(default={"property_package": m.fs.pp})
    super(CVFrameData, m.fs.cv).build()

    with pytest.raises(DynamicError):
        d, h = m.fs.cv._validate_add_balance_arguments(dynamic=True,
                                                       has_holdup=True)

# -----------------------------------------------------------------------------
# Test NotImplementedErrors for all property and balance type methods
def test_add_state_blocks():
    m = ConcreteModel()
    m.cv = CVFrame()

    with pytest.raises(NotImplementedError):
        m.cv.add_state_blocks()


def test_add_reaction_blocks():
    m = ConcreteModel()
    m.cv = CVFrame()

    with pytest.raises(NotImplementedError):
        m.cv.add_reaction_blocks()


def test_add_material_balances():
    m = ConcreteModel()
    m.cv = CVFrame()

    for t in MaterialBalanceType:
        if t == MaterialBalanceType.none:
            assert m.cv.add_material_balances(t) is None
        else:
            with pytest.raises(NotImplementedError):
                m.cv.add_material_balances(t)


def test_add_energy_balances():
    m = ConcreteModel()
    m.cv = CVFrame()

    for t in EnergyBalanceType:
        if t == EnergyBalanceType.none:
            assert m.cv.add_energy_balances(t) is None
        else:
            with pytest.raises(NotImplementedError):
                m.cv.add_energy_balances(t)


def test_add_momentum_balances():
    m = ConcreteModel()
    m.cv = CVFrame()

    for t in MomentumBalanceType:
        if t == MomentumBalanceType.none:
            assert m.cv.add_momentum_balances(t) is None
        else:
            with pytest.raises(NotImplementedError):
                m.cv.add_momentum_balances(t)
