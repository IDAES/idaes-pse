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
Tests for unit_model.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel, Expression, Set, Var
from idaes.core import (FlowsheetBlockData, declare_process_block_class,
                        UnitBlockData)
#, PropertyParameterBase, PropertyBlockDataBase)
from pyomo.common.config import ConfigValue


@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("Unit")
class UnitData(UnitBlockData):
    def build(self):
        pass


def test_config_block():
    m = ConcreteModel()

    m.u = Unit()

    assert m.u.config.dynamic == 'use_parent_value'


def test_config_args():
    m = ConcreteModel()

    m.u = Unit(dynamic=True)

    assert m.u.config.dynamic is True


def test_config_args_invalid():
    # Test validation of config arguments
    m = ConcreteModel()

    m.u = Unit()

    m.u.config.dynamic = True
    m.u.config.dynamic = False
    m.u.config.dynamic = 'use_parent_value'

    # Test that Value error raised when given invalid config arguments
    with pytest.raises(ValueError):
        m.u.config.dynamic = "foo"  # invalid str
    with pytest.raises(ValueError):
        m.u.config.dynamic = 5  # invalid int
    with pytest.raises(ValueError):
        m.u.config.dynamic = 2.0  # invalid float
    with pytest.raises(ValueError):
        m.u.config.dynamic = [2.0]  # invalid list
    with pytest.raises(ValueError):
        m.u.config.dynamic = {'a': 2.0}  # invalid dict


def test_setup_dynamics1():
    # Test that _setup_dynamics gets argument from parent
    m = ConcreteModel()

    m.fs = Flowsheet(dynamic=False)

    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    assert m.fs.u.config.dynamic is False


def test_setup_dynamics2():
    # Test that _setup_dynamics returns an Attribute error when parent has no
    # dynamic config argument

    m = ConcreteModel()
    m.u = Unit()

    with pytest.raises(AttributeError):
        m.u._setup_dynamics()


def test_setup_dynamics_dynamic_in_steady_state():
    # Test that a Value Error is raised when a dynamic models is placed in a
    # steady-state parent
    m = ConcreteModel()

    m.fs = Flowsheet(dynamic=False)

    m.fs.u = Unit(dynamic=True)
    with pytest.raises(ValueError):
        m.fs.u._setup_dynamics()


def test_setup_dynamics_get_time():
    # Test that time domain is collected correctly
    m = ConcreteModel()

    m.fs = Flowsheet(dynamic=False)

    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    assert m.fs.u.time == m.fs.time


def test_setup_dynamics_get_time_fails():
    # Test that AttributeError is raised when parent does not have time domain
    m = ConcreteModel()

    m.u = Unit()
    with pytest.raises(AttributeError):
        m.fs.u._setup_dynamics()


def test_setup_dynamics_include_holdup():
    # Test that include holdup argument is True when dynamic is True
    m = ConcreteModel()

    m.fs = Flowsheet(dynamic=True)

    m.fs.u = Unit()
    m.fs.u.config.declare("include_holdup", ConfigValue(default=False))
    m.fs.u._setup_dynamics()

    assert m.fs.u.config.include_holdup is True


## -----------------------------------------------------------------------------
## Setup a basic unit model for testing
#@declare_process_block_class("PropertyParameterBlock")
#class PropertyParameterData(PropertyParameterBase):
#    def build(self):
#        super(PropertyParameterData, self).build()
#
#    def get_supported_properties(self):
#        return {'pressure': {'method': None, 'units': 'Pa'},
#                'temperature': {'method': None, 'units': 'K'}}
#
#    def get_package_units(self):
#        return {'time': 's',
#                'length': 'm',
#                'mass': 'g',
#                'amount': 'mol',
#                'temperature': 'K',
#                'energy': 'J',
#                'holdup': 'mol'}
#
#@declare_process_block_class("PropertyBlock")
#class PropertyBlockData(PropertyBlockDataBase):
#    def build(self):
#        super(PropertyBlockData, self).build()
#        self.phase_list = Set(initialize=["phase1", "phase2"])
#        self.component_list = Set(initialize=["comp1", "comp2"])
#        self.pressure = Var(initialize=1.0)
#        self.temperature = Var(initialize=1.0)
#
#        def material_terms(b, p, j):
#            return 1
#        self.material_balance_term = Expression(self.phase_list,
#                                                self.component_list,
#                                                rule=material_terms)
#
#        def energy_terms(b, p):
#            return 2
#        self.energy_balance_term = Expression(self.phase_list,
#                                              rule=energy_terms)
#
#        def material_density_terms(b, p, j):
#            return 3
#        self.material_density_term = Expression(self.phase_list,
#                                                self.component_list,
#                                                rule=material_density_terms)
#
#        def energy_density_terms(b, p):
#            return 4
#        self.energy_density_term = Expression(self.phase_list,
#                                              rule=energy_density_terms)
#
#    def declare_port_members(b):
#        members = {"1": 1,
#                   "2": 2,
#                   "3": 3,
#                   "pressure": b.pressure}
#        return members
#
#    def model_check(b):
#        b.test_check = True
#
#
#@pytest.fixture()
#def m():
#    """
#    Build a basic Unit model with one Holdup for testing.
#    """
#    m = ConcreteModel()
#    m.fs = Flowsheet(dynamic=False)
#    m.fs.prop = PropertyParameterBlock()
#    m.fs.config.default_property_package = m.fs.prop
#    m.fs.u = Unit()
#    m.fs.u.config.declare("property_package", ConfigValue(default=m.fs.prop))
#    m.fs.u.config.declare("property_package_args", ConfigValue(default={}))
#    m.fs.u._setup_dynamics()
#    m.fs.u.holdup = Holdup0D()
#    return m
## -----------------------------------------------------------------------------
#
#
#def test_build_inlets_no_holdup():
#    # Check that build_inlets returns error when holdup does not exist
#    m = ConcreteModel()
#    m.fs = Flowsheet(dynamic=False)
#    m.fs.u = Unit()
#    m.fs.u._setup_dynamics()
#
#    with pytest.raises(AttributeError):
#        m.fs.u.build_inlets()
#
#
#def test_build_inlets_no_args(m):
#    m.fs.u.build_inlets()
#
#    assert m.fs.u.inlet[0].vars["1"] == 1
#    assert m.fs.u.inlet[0].vars["2"] == 2
#    assert m.fs.u.inlet[0].vars["3"] == 3
#    assert hasattr(m.fs.u.holdup, "inlet_mixer") is False
#
#
#def test_build_inlets_holdup_arg(m):
#    m.fs.u.build_inlets(holdup="holdup")
#
#    assert m.fs.u.holdup_inlet[0].vars["1"] == 1
#    assert m.fs.u.holdup_inlet[0].vars["2"] == 2
#    assert m.fs.u.holdup_inlet[0].vars["3"] == 3
#    assert hasattr(m.fs.u.holdup, "inlet_mixer") is False
#
#
#def test_build_inlets_invalid_holdup(m):
#    with pytest.raises(AttributeError):
#        m.fs.u.build_inlets(holdup="foo")
#
#
#def test_build_inlets_multiple_inlets(m):
#    m.fs.u.build_inlets(inlets=["1", "2"])
#
#    assert m.fs.u.inlet[0, "1"].vars["1"] == 1
#    assert m.fs.u.inlet[0, "1"].vars["2"] == 2
#    assert m.fs.u.inlet[0, "1"].vars["3"] == 3
#    assert m.fs.u.inlet[0, "2"].vars["1"] == 1
#    assert m.fs.u.inlet[0, "2"].vars["2"] == 2
#    assert m.fs.u.inlet[0, "2"].vars["3"] == 3
#    assert hasattr(m.fs.u.holdup, "inlet_mixer")
#
#
#def test_build_inlets_invalid_inlet_arg1(m):
#    with pytest.raises(TypeError):
#        m.fs.u.build_inlets(inlets=1)
#
#
#def test_build_inlets_invalid_inlet_arg2(m):
#    with pytest.raises(TypeError):
#        m.fs.u.build_inlets(inlets="foo")
#
#
#def test_build_inlets_invalid_inlet_arg3(m):
#    with pytest.raises(TypeError):
#        m.fs.u.build_inlets(inlets={"1": 1})
#
#
#def test_build_inlets_duplicate_inlet_names(m):
#    with pytest.raises(ValueError):
#        m.fs.u.build_inlets(inlets=["1", "2", "1"])
#
#
#def test_build_inlets_num_inlets(m):
#    m.fs.u.build_inlets(num_inlets=2)
#
#    assert m.fs.u.inlet[0, "1"].vars["1"] == 1
#    assert m.fs.u.inlet[0, "1"].vars["2"] == 2
#    assert m.fs.u.inlet[0, "1"].vars["3"] == 3
#    assert m.fs.u.inlet[0, "2"].vars["1"] == 1
#    assert m.fs.u.inlet[0, "2"].vars["2"] == 2
#    assert m.fs.u.inlet[0, "2"].vars["3"] == 3
#    assert hasattr(m.fs.u.holdup, "inlet_mixer")
#
#
#def test_build_inlets_inconsistent_args(m):
#    with pytest.raises(ValueError):
#        m.fs.u.build_inlets(inlets=["1", "2", "3"], num_inlets=2)
#
#
#def test_build_inlets_both_args_consistent(m):
#    m.fs.u.build_inlets(inlets=["1", "2", "3"], num_inlets=3)
#
#
## Test that no errors arise when working with Holdup1D
#def test_build_inlets_1d_forward():
#    m = ConcreteModel()
#    m.fs = Flowsheet(dynamic=False)
#    m.fs.prop = PropertyParameterBlock()
#    m.fs.config.default_property_package = m.fs.prop
#    m.fs.u = Unit()
#    m.fs.u.config.declare("property_package", ConfigValue(default=m.fs.prop))
#    m.fs.u.config.declare("property_package_args", ConfigValue(default={}))
#    m.fs.u._setup_dynamics()
#    m.fs.u.holdup = Holdup1D(flow_direction="forward")
#
#    m.fs.u.build_inlets()
#    m.fs.u.holdup.properties[0, 0].pressure = 10
#    assert m.fs.u.inlet[0].vars["pressure"].value == 10
#
#
#def test_build_inlets_1d_backward():
#    m = ConcreteModel()
#    m.fs = Flowsheet(dynamic=False)
#    m.fs.prop = PropertyParameterBlock()
#    m.fs.config.default_property_package = m.fs.prop
#    m.fs.u = Unit()
#    m.fs.u.config.declare("property_package", ConfigValue(default=m.fs.prop))
#    m.fs.u.config.declare("property_package_args", ConfigValue(default={}))
#    m.fs.u._setup_dynamics()
#    m.fs.u.holdup = Holdup1D(flow_direction="backward")
#
#    m.fs.u.build_inlets()
#    m.fs.u.holdup.properties[0, 1].pressure = 10
#    assert m.fs.u.inlet[0].vars["pressure"].value == 10
#
## -----------------------------------------------------------------------------
#
#
#def test_build_outlets_no_holdup():
#    # Check that build_outlets returns error when holdup does not exist
#    m = ConcreteModel()
#    m.fs = Flowsheet(dynamic=False)
#    m.fs.u = Unit()
#    m.fs.u._setup_dynamics()
#
#    with pytest.raises(AttributeError):
#        m.fs.u.build_outlets()
#
#
#def test_build_outlets_no_args(m):
#    m.fs.u.build_outlets()
#
#    assert m.fs.u.outlet[0].vars["1"] == 1
#    assert m.fs.u.outlet[0].vars["2"] == 2
#    assert m.fs.u.outlet[0].vars["3"] == 3
#    assert hasattr(m.fs.u.holdup, "outlet_splitter") is False
#
#
#def test_build_outlets_holdup_arg(m):
#    m.fs.u.build_outlets(holdup="holdup")
#
#    assert m.fs.u.holdup_outlet[0].vars["1"] == 1
#    assert m.fs.u.holdup_outlet[0].vars["2"] == 2
#    assert m.fs.u.holdup_outlet[0].vars["3"] == 3
#    assert hasattr(m.fs.u.holdup, "outlet_splitter") is False
#
#
#def test_build_outlets_invalid_holdup(m):
#    with pytest.raises(AttributeError):
#        m.fs.u.build_outlets(holdup="foo")
#
#
#def test_build_outlets_multiple_outlets(m):
#    m.fs.u.build_outlets(outlets=["1", "2"])
#
#    assert m.fs.u.outlet[0, "1"].vars["1"] == 1
#    assert m.fs.u.outlet[0, "1"].vars["2"] == 2
#    assert m.fs.u.outlet[0, "1"].vars["3"] == 3
#    assert m.fs.u.outlet[0, "2"].vars["1"] == 1
#    assert m.fs.u.outlet[0, "2"].vars["2"] == 2
#    assert m.fs.u.outlet[0, "2"].vars["3"] == 3
#    assert hasattr(m.fs.u.holdup, "outlet_splitter")
#
#
#def test_build_outlets_invalid_outlet_arg1(m):
#    with pytest.raises(TypeError):
#        m.fs.u.build_outlets(outlets=1)
#
#
#def test_build_outlets_invalid_outlet_arg2(m):
#    with pytest.raises(TypeError):
#        m.fs.u.build_outlets(outlets="foo")
#
#
#def test_build_outlets_invalid_outlet_arg3(m):
#    with pytest.raises(TypeError):
#        m.fs.u.build_outlets(outlets={"1": 1})
#
#
#def test_build_outlets_duplicate_outlet_names(m):
#    with pytest.raises(ValueError):
#        m.fs.u.build_outlets(outlets=["1", "2", "1"])
#
#
#def test_build_outlets_num_outlets(m):
#    m.fs.u.build_outlets(num_outlets=2)
#
#    assert m.fs.u.outlet[0, "1"].vars["1"] == 1
#    assert m.fs.u.outlet[0, "1"].vars["2"] == 2
#    assert m.fs.u.outlet[0, "1"].vars["3"] == 3
#    assert m.fs.u.outlet[0, "2"].vars["1"] == 1
#    assert m.fs.u.outlet[0, "2"].vars["2"] == 2
#    assert m.fs.u.outlet[0, "2"].vars["3"] == 3
#    assert hasattr(m.fs.u.holdup, "outlet_splitter")
#
#
#def test_build_outlets_inconsistent_args(m):
#    with pytest.raises(ValueError):
#        m.fs.u.build_outlets(outlets=["1", "2", "3"], num_outlets=2)
#
#
#def test_build_outlets_both_args_consistent(m):
#    m.fs.u.build_outlets(outlets=["1", "2", "3"], num_outlets=3)
#
#
## Test that no errors arise when working with Holdup1D
#def test_build_outlets_1d_forward():
#    m = ConcreteModel()
#    m.fs = Flowsheet(dynamic=False)
#    m.fs.prop = PropertyParameterBlock()
#    m.fs.config.default_property_package = m.fs.prop
#    m.fs.u = Unit()
#    m.fs.u.config.declare("property_package", ConfigValue(default=m.fs.prop))
#    m.fs.u.config.declare("property_package_args", ConfigValue(default={}))
#    m.fs.u._setup_dynamics()
#    m.fs.u.holdup = Holdup1D(flow_direction="forward")
#
#    m.fs.u.build_outlets()
#    m.fs.u.holdup.properties[0, 1].pressure = 10
#    assert m.fs.u.outlet[0].vars["pressure"].value == 10
#
#
#def test_build_outlets_1d_backward():
#    m = ConcreteModel()
#    m.fs = Flowsheet(dynamic=False)
#    m.fs.prop = PropertyParameterBlock()
#    m.fs.config.default_property_package = m.fs.prop
#    m.fs.u = Unit()
#    m.fs.u.config.declare("property_package", ConfigValue(default=m.fs.prop))
#    m.fs.u.config.declare("property_package_args", ConfigValue(default={}))
#    m.fs.u._setup_dynamics()
#    m.fs.u.holdup = Holdup1D(flow_direction="backward")
#
#    m.fs.u.build_outlets()
#    m.fs.u.holdup.properties[0, 0].pressure = 10
#    assert m.fs.u.outlet[0].vars["pressure"].value == 10
#
#
## -----------------------------------------------------------------------------
## Test model_check
#def test_model_check(m):
#    m.fs.u.model_check()
#
#    assert m.fs.u.holdup.properties_in[0].test_check
#    assert m.fs.u.holdup.properties_out[0].test_check
