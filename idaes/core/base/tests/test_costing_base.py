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
Tests for costing base classes
"""
import pytest

from pyomo.environ import ConcreteModel, Constraint, Set, units as pyunits, Var, Param
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core import (
    FlowsheetCostingBlock,
    FlowsheetCostingBlockData,
    UnitModelCostingBlock,
    register_idaes_currency_units,
)

# TODO : Tests for cases with multiple costing packages
pyunits.load_definitions_from_strings(["USD_test = [test_currency]"])


@pytest.mark.unit
def test_register_idaes_currency_units():
    # Register standard currency units
    register_idaes_currency_units()

    assert "USD_CE500" in pyunits.pint_registry

    CEI = {
        "USD_1990": 357.6,
        "USD_1991": 361.3,
        "USD_1992": 358.2,
        "USD_1993": 359.2,
        "USD_1994": 368.1,
        "USD_1995": 381.1,
        "USD_1996": 381.7,
        "USD_1997": 386.5,
        "USD_1998": 389.5,
        "USD_1999": 390.6,
        "USD_2000": 394.1,
        "USD_2001": 394.3,
        "USD_2002": 395.6,
        "USD_2003": 402.0,
        "USD_2004": 444.2,
        "USD_2005": 468.2,
        "USD_2006": 499.6,
        "USD_2007": 525.4,
        "USD_2008": 575.4,
        "USD_2009": 521.9,
        "USD_2010": 550.8,
        "USD_2011": 585.7,
        "USD_2012": 584.6,
        "USD_2013": 567.3,
        "USD_2014": 576.1,
        "USD_2015": 556.8,
        "USD_2016": 541.7,
        "USD_2017": 567.5,
        "USD_2018": 603.1,
        "USD_2019": 607.5,
        "USD_2020": 596.2,
    }

    for c, conv in CEI.items():
        assert c in pyunits.pint_registry

        assert pytest.approx(conv / 500, rel=1e-10) == pyunits.convert_value(
            1, pyunits.USD_CE500, getattr(pyunits, c)
        )


class TestFlowsheetCostingBlock:
    @pytest.mark.unit
    def test_basic_attributes(self):
        m = ConcreteModel()

        with pytest.raises(
            ValueError,
            match="costing - costing package has not specified the base "
            "currency units to use for costing.",
        ):
            m.costing = FlowsheetCostingBlock()

        assert m.costing.base_currency is None
        assert m.costing.base_period is pyunits.year
        assert m.costing.defined_flows == {}
        assert FlowsheetCostingBlockData.unit_mapping == {}

    @pytest.mark.unit
    def test_build_global_params(self):
        with pytest.raises(
            NotImplementedError,
            match="Derived class has not defined a " "build_global_params method.",
        ):
            FlowsheetCostingBlockData.build_global_params(self)

    @pytest.mark.unit
    def test_build_process_costs(self):
        with pytest.raises(
            NotImplementedError,
            match="Derived class has not defined a " "build_process_costs method.",
        ):
            FlowsheetCostingBlockData.build_process_costs(self)

    @pytest.mark.unit
    def test_initialize_build(self):
        with pytest.raises(
            NotImplementedError,
            match="Derived class has not defined an " "initialize_build method.",
        ):
            FlowsheetCostingBlockData.initialize(self)


# Create some dummy classes to represent inherited unit models
@declare_process_block_class("TypeA")
class TypeAData(UnitModelBlockData):
    def build(self):
        super().build()
        self.class_type = "A"


@declare_process_block_class("TypeB")
class TypeBData(TypeAData):
    def build(self):
        super().build()
        self.class_type = "B"


@declare_process_block_class("TypeC")
class TypeCData(TypeBData):
    def build(self):
        super().build()
        self.class_type = "C"


@declare_process_block_class("TypeD")
class TypeDData(TypeAData):
    def build(self):
        super().build()
        self.class_type = "D"


@declare_process_block_class("TypeE")
class TypeEData(UnitModelBlockData):
    def build(self):
        super().build()
        self.class_type = "E"


@declare_process_block_class("TestCostingPackage")
class TestCostingPackageData(FlowsheetCostingBlockData):
    def build_global_params(self):
        self.base_currency = pyunits.USD_test
        self.base_period = pyunits.year

        self.test_flow_2_cost = Param(initialize=0.07, units=pyunits.kW)

        self.defined_flows = {
            "test_flow_1": 0.2 * pyunits.J,
            "test_flow_2": self.test_flow_2_cost,
        }

        self._bgp = True

    def build_process_costs(self):
        self._bpc = True

    def initialize_build(self):
        self._init = True

    def method_1(blk):
        blk.cost_method = 1

    def method_2(blk):
        blk.cost_method = 2

    def method_3(blk):
        blk.cost_method = 3

    def method_4(blk):
        blk.cost_method = 4

    unit_mapping = {TypeA: method_1, TypeB: method_2, TypeC: method_3}


class TestFlowsheetCostingBlock:
    @pytest.mark.unit
    def test_costing_package_no_base_currency(self):
        @declare_process_block_class("TestCostingPackage2")
        class TestCostingPackage2Data(TestCostingPackageData):
            def build_global_params(self):
                super().build_global_params()
                self.base_currency = None

        m = ConcreteModel()

        with pytest.raises(
            ValueError,
            match="costing - costing package has not specified the base "
            "currency units to use for costing.",
        ):
            m.costing = TestCostingPackage2()

    @pytest.fixture(scope="class")
    def costing(self):
        m = ConcreteModel()
        m.costing = TestCostingPackage()

        return m

    @pytest.mark.unit
    def test_basic_attributes(self, costing):
        assert costing.costing._registered_unit_costing == []
        assert isinstance(costing.costing.flow_types, Set)
        assert len(costing.costing.flow_types) == 2
        assert "test_flow_1" in costing.costing.flow_types
        assert costing.costing._registered_flows == {
            "test_flow_1": [],
            "test_flow_2": [],
        }

        assert costing.costing._costing_methods_map == {
            TypeAData: TestCostingPackageData.method_1,
            TypeBData: TestCostingPackageData.method_2,
            TypeCData: TestCostingPackageData.method_3,
        }

        # Check that test_flow_1 was properly defined
        assert isinstance(costing.costing.test_flow_1_cost, Var)
        assert costing.costing.test_flow_1_cost.value == 0.2
        assert_units_equivalent(costing.costing.test_flow_1_cost.get_units(), pyunits.J)

        assert isinstance(costing.costing.test_flow_2_cost, Param)

        # Test that build_global_parameters was called successfully
        assert costing.costing._bgp

    @pytest.mark.unit
    def test_register_flow_type(self, costing):
        costing.costing.register_flow_type(
            "test_flow", 42 * pyunits.USD_test / pyunits.mol
        )

        assert isinstance(costing.costing.test_flow_cost, Var)
        assert costing.costing.test_flow_cost.value == 42
        assert_units_equivalent(
            costing.costing.test_flow_cost.get_units(), pyunits.USD_test / pyunits.mol
        )
        assert "test_flow" in costing.costing.flow_types

        assert costing.costing._registered_flows == {
            "test_flow_1": [],
            "test_flow_2": [],
            "test_flow": [],
        }

    @pytest.mark.unit
    def test_register_flow_component_exists(self, costing):
        costing.costing.test_flow_3_cost = Var()
        with pytest.raises(
            RuntimeError,
            match="Component test_flow_3_cost already exists on costing but is not 42.",
        ):
            costing.costing.register_flow_type("test_flow_3", 42)
        # cleanup for next test
        costing.costing.flow_types.remove("test_flow_3")
        costing.costing.del_component(costing.costing.test_flow_3_cost)

    @pytest.mark.unit
    def test_cost_flow_invalid_type(self, costing):
        with pytest.raises(
            ValueError,
            match="foo is not a recognized flow type. Please "
            "check your spelling and that the flow type has "
            "been registered with the FlowsheetCostingBlock.",
        ):
            costing.costing.cost_flow(42, "foo")

    @pytest.mark.unit
    def test_cost_flow_indexed_var(self, costing):
        costing.indexed_var = Var(
            [1, 2, 3], initialize=1, units=pyunits.mol / pyunits.s
        )
        with pytest.raises(
            TypeError,
            match="indexed_var is an indexed component. Flow "
            "costing only supports unindexed components.",
        ):
            costing.costing.cost_flow(costing.indexed_var, "test_flow")

    @pytest.mark.unit
    def test_cost_flow_unbounded_var(self, costing, caplog):
        costing.costing.cost_flow(costing.indexed_var[1], "test_flow")

        warn_str = (
            "indexed_var[1] has a lower bound of less "
            "than zero. Costing requires that all flows have a "
            "lower bound equal to or greater than zero to "
            "avoid negative costs."
        )

        assert warn_str in caplog.text

    @pytest.mark.unit
    def test_cost_flow_var(self, costing):
        costing.indexed_var[1].setlb(0)

        costing.costing.cost_flow(costing.indexed_var[1], "test_flow")

        assert costing.indexed_var[1] in costing.costing._registered_flows["test_flow"]

    @pytest.mark.unit
    def test_cost_flow_unbounded_expr(self, costing, caplog):
        costing.costing.cost_flow(-costing.indexed_var[2], "test_flow")

        warn_str = (
            "flow_expr is an expression with a lower "
            "bound of less than zero. Costing requires that "
            "all flows have a lower bound equal to or greater "
            "than zero to avoid negative costs."
        )

        assert warn_str in caplog.text

    @pytest.mark.unit
    def test_cost_flow_expr(self, costing):
        costing.indexed_var[2].setub(0)

        costing.costing.cost_flow(-costing.indexed_var[2], "test_flow")

        assert str(-costing.indexed_var[2]) == str(
            costing.costing._registered_flows["test_flow"][-1]
        )

    @pytest.mark.unit
    def test_get_costing_method_for(self, costing):
        costing.unit_a = TypeA()
        costing.unit_b = TypeB()
        costing.unit_c = TypeC()
        costing.unit_d = TypeD()
        costing.unit_e = TypeE()

        assert isinstance(costing.unit_a, TypeAData)

        assert (
            costing.costing._get_costing_method_for(costing.unit_a)
            is TestCostingPackageData.method_1
        )

        assert (
            costing.costing._get_costing_method_for(costing.unit_b)
            is TestCostingPackageData.method_2
        )

        assert (
            costing.costing._get_costing_method_for(costing.unit_c)
            is TestCostingPackageData.method_3
        )

        # TypeD not registered with property package, but inherits from TypeA
        # Should get method_1
        assert (
            costing.costing._get_costing_method_for(costing.unit_d)
            is TestCostingPackageData.method_1
        )

        # TypeE not registered with property package and no inheritance
        # Should return RuntimeError
        with pytest.raises(
            RuntimeError,
            match="Could not identify default costing method "
            "for unit_e. This implies the unit model's class "
            "and parent classes do not exist in the default "
            "mapping provided by the costing package. Please "
            "provide a specific costing method for this unit.",
        ):
            costing.costing._get_costing_method_for(costing.unit_e)

    @pytest.mark.unit
    def test_cost_unit_first(self, costing):
        costing.unit_a.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing.costing
        )

        assert costing.unit_a.costing.cost_method == 1

        assert costing.unit_a.costing in costing.unit_a._initialization_order
        assert costing.unit_a.costing in costing.costing._registered_unit_costing

    @pytest.mark.unit
    def test_del_unit_costing(self, costing):
        costing.unit_a.del_component(costing.unit_a.costing)

        assert not hasattr(costing.unit_a, "costing")
        assert costing.unit_a._initialization_order == []
        assert costing.costing._registered_unit_costing == []

    @pytest.mark.unit
    def test_cost_unit_duplicate(self, costing):
        costing.unit_a.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing.costing
        )

        # First, check implicit replacement
        # This should work
        costing.unit_a.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing.costing,
            costing_method=TestCostingPackageData.method_2,
        )

        assert costing.unit_a.costing.cost_method == 2

        assert costing.unit_a.costing in costing.unit_a._initialization_order
        assert costing.unit_a.costing in costing.costing._registered_unit_costing

        # Then check double costing
        with pytest.raises(
            RuntimeError,
            match="Unit model unit_a already has a costing block "
            "registered: unit_a.costing. Each unit may only have a single "
            "UnitModelCostingBlock associated with it.",
        ):
            costing.unit_a.costing2 = UnitModelCostingBlock(
                flowsheet_costing_block=costing.costing
            )

        # Clean everything up at the end
        costing.unit_a.del_component(costing.unit_a.costing)
        costing.unit_a.del_component(costing.unit_a.costing2)

        # Make sure we cleaned up
        assert not hasattr(costing.unit_a, "costing")
        assert not hasattr(costing.unit_a, "costing2")

    @pytest.mark.unit
    def test_cost_unit_custom_method(self, costing):
        def custom_method(blk):
            blk.capital_cost = Var(
                initialize=1, bounds=(0, 1e10), units=pyunits.USD_test
            )
            blk.fixed_operating_cost = Var(
                initialize=1, bounds=(0, 1e10), units=pyunits.USD_test / pyunits.year
            )
            blk.variable_operating_cost = Var(
                initialize=1, bounds=(0, 1e10), units=pyunits.USD_test / pyunits.year
            )

            blk.capital_cost_constraint = Constraint(
                expr=blk.capital_cost == 4.2e6 * pyunits.USD_test
            )
            blk.fixed_operating_cost_constraint = Constraint(
                expr=blk.fixed_operating_cost == 1e2 * pyunits.USD_test / pyunits.year
            )
            blk.variable_operating_cost_constraint = Constraint(
                expr=blk.variable_operating_cost
                == 7e4 * pyunits.USD_test / pyunits.year
            )

            blk._checkvar = True

        costing.unit_a.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing.costing, costing_method=custom_method
        )

        assert isinstance(costing.unit_a.costing, UnitModelCostingBlock)
        assert costing.unit_a.costing in costing.unit_a._initialization_order
        assert costing.unit_a.costing in costing.costing._registered_unit_costing

        assert isinstance(costing.unit_a.costing.capital_cost, Var)
        assert isinstance(costing.unit_a.costing.variable_operating_cost, Var)
        assert isinstance(costing.unit_a.costing.fixed_operating_cost, Var)
        assert costing.unit_a.costing._checkvar

    @pytest.mark.unit
    def test_cost_unit_capital_cost_not_var(self, costing):
        def dummy_method(blk):
            blk.capital_cost = "foo"

        with pytest.raises(
            TypeError,
            match="unit_b capital_cost component must be a "
            "Var. Please check the costing package you are "
            "using to ensure that all costing components are "
            "declared as variables.",
        ):
            costing.unit_b.costing = UnitModelCostingBlock(
                flowsheet_costing_block=costing.costing, costing_method=dummy_method
            )

        # Clean up for next test
        costing.unit_b.del_component(costing.unit_b.costing)

    @pytest.mark.unit
    def test_cost_unit_capital_cost_lb(self, costing, caplog):
        def dummy_method(blk):
            blk.capital_cost = Var()

        costing.unit_b.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing.costing, costing_method=dummy_method
        )

        assert (
            "unit_b capital_cost component has a lower bound less than "
            "zero. Be aware that this may result in negative costs during "
            "optimization." in caplog.text
        )

        # Clean up for next test
        costing.unit_b.del_component(costing.unit_b.costing)

    @pytest.mark.unit
    def test_cost_unit_fixed_operating_cost_not_var(self, costing):
        def dummy_method(blk):
            blk.fixed_operating_cost = "foo"

        with pytest.raises(
            TypeError,
            match="unit_b fixed_operating_cost component must "
            "be a Var. Please check the costing package you "
            "are using to ensure that all costing components "
            "are declared as variables.",
        ):
            costing.unit_b.costing = UnitModelCostingBlock(
                flowsheet_costing_block=costing.costing, costing_method=dummy_method
            )

        # Clean up for next test
        costing.unit_b.del_component(costing.unit_b.costing)

    @pytest.mark.unit
    def test_cost_unit_fixed_operating_cost_lb(self, costing, caplog):
        def dummy_method(blk):
            blk.fixed_operating_cost = Var()

        costing.unit_b.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing.costing, costing_method=dummy_method
        )

        assert (
            "unit_b fixed_operating_cost component has a lower bound less "
            "than zero. Be aware that this may result in negative costs "
            "during optimization." in caplog.text
        )

        # Clean up for next test
        costing.unit_b.del_component(costing.unit_b.costing)

    @pytest.mark.unit
    def test_cost_unit_variable_operating_cost_not_var(self, costing):
        def dummy_method(blk):
            blk.variable_operating_cost = "foo"

        with pytest.raises(
            TypeError,
            match="unit_b variable_operating_cost component "
            "must be a Var. Please check the costing package "
            "you are using to ensure that all costing "
            "components are declared as variables.",
        ):
            costing.unit_b.costing = UnitModelCostingBlock(
                flowsheet_costing_block=costing.costing, costing_method=dummy_method
            )

        # Clean up for next test
        costing.unit_b.del_component(costing.unit_b.costing)

    @pytest.mark.unit
    def test_cost_unit_variable_operating_cost_lb(self, costing, caplog):
        def dummy_method(blk):
            blk.variable_operating_cost = Var()

        costing.unit_b.costing = UnitModelCostingBlock(
            flowsheet_costing_block=costing.costing, costing_method=dummy_method
        )

        assert (
            "unit_b variable_operating_cost component has a lower bound "
            "less than zero. Be aware that this may result in negative "
            "costs during optimization." in caplog.text
        )

        # Clean up for next test
        costing.unit_b.del_component(costing.unit_b.costing)

    @pytest.mark.unit
    def test_cost_process(self, costing):
        costing.costing.cost_process()

        # Check that build_process_costs was called from costing package
        assert costing.costing._bpc

        # Then check aggregation
        assert isinstance(costing.costing.aggregate_capital_cost, Var)
        assert str(costing.costing.aggregate_capital_cost.get_units()) == str(
            pyunits.USD_test
        )
        assert isinstance(costing.costing.aggregate_capital_cost_constraint, Constraint)

        assert isinstance(costing.costing.aggregate_fixed_operating_cost, Var)
        assert str(pyunits.USD_test / pyunits.year) == str(
            costing.costing.aggregate_fixed_operating_cost.get_units()
        )
        assert isinstance(
            costing.costing.aggregate_fixed_operating_cost_constraint, Constraint
        )

        assert isinstance(costing.costing.aggregate_variable_operating_cost, Var)
        assert str(pyunits.USD_test / pyunits.year) == str(
            costing.costing.aggregate_variable_operating_cost.get_units()
        )
        assert isinstance(
            costing.costing.aggregate_variable_operating_cost_constraint, Constraint
        )

        assert isinstance(costing.costing.aggregate_flow_test_flow, Var)
        assert str(pyunits.mol / pyunits.s) == str(
            costing.costing.aggregate_flow_test_flow.get_units()
        )
        assert isinstance(
            costing.costing.aggregate_flow_test_flow_constraint, Constraint
        )

        # We also have a test_flow_1 type registered, but no flows costed
        # This should have been skipped
        assert not hasattr(costing.costing, "aggregate_flow_test_flow_1")
        assert not hasattr(costing.costing, "aggregate_flow_test_flow_1_constraint")
        # unused flows do not get added to aggregate_flow_costs
        assert "test_flow_1" not in costing.costing.aggregate_flow_costs

        assert isinstance(costing.costing.aggregate_flow_costs, Var)
        assert str(pyunits.USD_test / pyunits.year) == str(
            costing.costing.aggregate_flow_costs.get_units()
        )
        assert len(costing.costing.aggregate_flow_costs) == 1
        assert isinstance(costing.costing.aggregate_flow_costs_constraint, Constraint)
        assert len(costing.costing.aggregate_flow_costs_constraint) == 1

    @pytest.mark.unit
    def test_unit_consistency(self, costing):
        assert_units_consistent(costing)

    @pytest.mark.unit
    def test_degrees_of_freedom(self, costing):
        costing.indexed_var[1].fix(2)
        costing.indexed_var[2].fix(-3)
        costing.indexed_var[3].fix(0)

        assert degrees_of_freedom(costing) == 0

    @pytest.mark.unit
    def test_initialize(self, costing):
        costing.costing.initialize()

        # Check that initialize was called from costing package
        assert costing.costing._init

        # Check that unit-level vars were initialized
        assert costing.unit_a.costing.capital_cost.value == 4.2e6
        assert costing.unit_a.costing.fixed_operating_cost.value == 100
        assert costing.unit_a.costing.variable_operating_cost.value == 7e4

        # Check that aggregate vars were initialized
        # Capital and operating costs should equal the unit level ones
        assert costing.costing.aggregate_capital_cost.value == 4.2e6
        assert costing.costing.aggregate_fixed_operating_cost.value == 100
        assert costing.costing.aggregate_variable_operating_cost.value == 7e4

        assert costing.costing.aggregate_flow_test_flow.value == 10

        assert pytest.approx(
            costing.costing.aggregate_flow_costs["test_flow"].value, rel=1e-12
        ) == (
            pyunits.convert_value(
                10 * 42, from_units=1 / pyunits.s, to_units=1 / pyunits.year
            )
        )
