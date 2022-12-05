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
Tests for flowsheet_model.

Author: Andrew Lee
"""
import pytest

from pyomo.environ import ConcreteModel, Set, TransformationFactory, units
from pyomo.dae import ContinuousSet
from pyomo.network import Arc

from idaes.core import (
    FlowsheetBlockData,
    FlowsheetBlock,
    UnitModelBlock,
    declare_process_block_class,
    useDefault,
)
from idaes.core.util.exceptions import DynamicError
from idaes.core.util.testing import PhysicalParameterTestBlock

from idaes.models.unit_models import Heater


@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(FlowsheetBlockData, self).build()


class TestConfig(object):
    @pytest.fixture()
    def model(self):
        m = ConcreteModel()
        m.fs = Flowsheet()

        return m

    @pytest.mark.unit
    def test_config(self, model):
        assert len(model.fs.config) == 5
        assert model.fs.config.dynamic is useDefault
        assert model.fs.config.time is None
        assert model.fs.config.time_set == [0]
        assert model.fs.config.default_property_package is None
        assert model.fs.config.time_units is None

    @pytest.mark.unit
    def test_config_validation_dynamic(self, model):
        # Test validation of dynamic argument
        model.fs.config.dynamic = False
        assert model.fs.config.dynamic is False

        model.fs.config.dynamic = True
        assert model.fs.config.dynamic is True

        with pytest.raises(ValueError):
            model.fs.config.dynamic = "foo"
        with pytest.raises(ValueError):
            model.fs.config.dynamic = 2
        with pytest.raises(ValueError):
            model.fs.config.dynamic = 2.0
        with pytest.raises(ValueError):
            model.fs.config.dynamic = [1]
        with pytest.raises(ValueError):
            model.fs.config.dynamic = {"foo": 1}

    @pytest.mark.unit
    def test_config_validation_time(self, model):
        # Test validation of time argument
        model.test_set = Set(initialize=[0, 1, 2])
        model.test_contset = ContinuousSet(bounds=[0, 1])

        model.fs.config.time = model.test_set
        assert model.fs.config.time == model.test_set

        model.fs.config.time = model.test_contset
        assert model.fs.config.time == model.test_contset

        with pytest.raises(ValueError):
            model.fs.config.time = "foo"
        with pytest.raises(ValueError):
            model.fs.config.time = 1
        with pytest.raises(ValueError):
            model.fs.config.time = 2.0
        with pytest.raises(ValueError):
            model.fs.config.time = [1]
        with pytest.raises(ValueError):
            model.fs.config.time = {"foo": 1}

    @pytest.mark.unit
    def test_config_validation_time_set(self, model):
        # Test validation of time_set argument
        model.fs.config.time_set = [1, 2, 3]
        assert model.fs.config.time_set == [1, 2, 3]
        model.fs.config.time_set = 5
        assert model.fs.config.time_set == [5]
        model.fs.config.time_set = 2.0
        assert model.fs.config.time_set == [2.0]

        with pytest.raises(ValueError):
            model.fs.config.time_set = "foo"  # invalid str
        with pytest.raises(ValueError):
            model.fs.config.time_set = {"a": 2.0}  # invalid dict

    @pytest.mark.unit
    def test_config_validation_default_property_package(self, model):
        # Test default_property_package attribute
        model.fs.p = PhysicalParameterTestBlock()

        model.fs.config.default_property_package = model.fs.p

        # Test default_property_package - invalid values
        with pytest.raises(ValueError):
            model.fs.config.default_property_package = "foo"  # invalid str
        with pytest.raises(ValueError):
            model.fs.config.default_property_package = 5  # invalid int
        with pytest.raises(ValueError):
            model.fs.config.default_property_package = 2.0  # invalid float
        with pytest.raises(ValueError):
            model.fs.config.default_property_package = [2.0]  # invalid list
        with pytest.raises(ValueError):
            model.fs.config.default_property_package = {"a": 2.0}  # invalid dict


class TestBuild(object):
    # Test that build method works for all combinations of config arguments
    @pytest.mark.unit
    def test_is_flowsheet(self):
        # Test that flowsheet has is_flowsheet method and that it returns True
        m = ConcreteModel()
        m.fs = FlowsheetBlock()

        assert hasattr(m.fs, "is_flowsheet")
        assert m.fs.is_flowsheet()

    @pytest.mark.unit
    def test_flowsheet(self):
        # Should return None
        m = ConcreteModel()
        m.fs = FlowsheetBlock()

        assert m.fs.flowsheet() is None

    @pytest.mark.unit
    def test_default(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock()

        assert m.fs.config.dynamic is False
        assert isinstance(m.fs.time, Set)
        assert m.fs.time == [0]
        assert m.fs.config.time is m.fs.time
        assert m.fs.time_units is None

    @pytest.mark.unit
    def test_ss_default(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        assert m.fs.config.dynamic is False
        assert isinstance(m.fs.time, Set)
        assert m.fs.time == [0]
        assert m.fs.config.time is m.fs.time
        assert m.fs.time_units is None

    @pytest.mark.unit
    def test_ss_time_set(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False, time_set=[1, 2, 3])

        assert m.fs.config.dynamic is False
        assert isinstance(m.fs.time, Set)
        for t in m.fs.time:
            assert t in [1, 2, 3]
        assert len(m.fs.time) == 3
        assert m.fs.config.time is m.fs.time
        assert m.fs.time_units is None

    @pytest.mark.unit
    def test_dynamic_default(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)

        assert m.fs.config.dynamic is True
        assert isinstance(m.fs.time, ContinuousSet)
        for t in m.fs.time:
            assert t in [0, 1]
        assert m.fs.config.time is m.fs.time
        assert m.fs.time_units == units.s

    @pytest.mark.unit
    def test_dynamic_time_set(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_set=[1, 2], time_units=units.s)

        assert m.fs.config.dynamic is True
        assert isinstance(m.fs.time, ContinuousSet)
        for t in m.fs.time:
            assert t in [1, 2]
        assert m.fs.config.time is m.fs.time
        assert m.fs.time_units == units.s

    @pytest.mark.unit
    def test_dynamic_time_set_invalid(self):
        m = ConcreteModel()

        with pytest.raises(
            DynamicError,
            match="Flowsheet provided with invalid "
            "time_set attribute - must have at "
            "least two values \(start and end\).",
        ):
            m.fs = FlowsheetBlock(dynamic=True, time_set=1, time_units=units.s)

    @pytest.mark.unit
    def test_ss_external_time(self):
        m = ConcreteModel()
        m.s = Set(initialize=[4, 5])
        m.fs = FlowsheetBlock(dynamic=False, time=m.s)

        assert m.fs.config.dynamic is False
        assert m.fs.config.time is m.s
        assert m.fs.time is m.s
        assert m.fs.time_units is None

    @pytest.mark.unit
    def test_dynamic_external_time_continuous(self):
        m = ConcreteModel()
        m.s = ContinuousSet(initialize=[4, 5])
        m.fs = FlowsheetBlock(dynamic=False, time=m.s)

        assert m.fs.config.dynamic is False
        assert m.fs.config.time is m.s
        assert m.fs.time is m.s
        assert m.fs.time_units is None

    @pytest.mark.unit
    def test_dynamic_external_time(self):
        m = ConcreteModel()
        m.s = ContinuousSet(initialize=[4, 5])
        m.fs = FlowsheetBlock(dynamic=True, time=m.s, time_units=units.s)

        assert m.fs.config.dynamic is True
        assert m.fs.config.time is m.s
        assert m.fs.time is m.s
        assert m.fs.time_units == units.s

    @pytest.mark.unit
    def test_dynamic_external_time_invalid(self):
        m = ConcreteModel()
        m.s = Set(initialize=[4, 5])

        with pytest.raises(DynamicError):
            m.fs = FlowsheetBlock(dynamic=True, time=m.s, time_units=units.s)

    @pytest.mark.unit
    def test_ss_external_time_and_time_set(self):
        # Should ignore time set
        m = ConcreteModel()
        m.s = Set(initialize=[4, 5])
        m.fs = FlowsheetBlock(dynamic=False, time=m.s, time_set=[1, 2])

        assert m.fs.config.dynamic is False
        assert m.fs.config.time is m.s
        assert m.fs.time is m.s
        assert m.fs.time_units is None

    @pytest.mark.unit
    def test_dynamic_external_time_and_time_set(self):
        # Should ignore time set
        m = ConcreteModel()
        m.s = ContinuousSet(initialize=[4, 5])
        m.fs = FlowsheetBlock(
            dynamic=True, time=m.s, time_set=[1, 2], time_units=units.s
        )

        assert m.fs.config.dynamic is True
        assert m.fs.config.time is m.s
        assert m.fs.time is m.s
        assert m.fs.time_units == units.s

    @pytest.mark.unit
    def test_time_units_ss(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False, time_units=units.s)

        assert m.fs.time_units is units.s

    @pytest.mark.unit
    def test_time_units_dynamic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)

        assert m.fs.time_units is units.s

    @pytest.mark.unit
    def test_time_units_external(self):
        # Should ignore time set
        m = ConcreteModel()
        m.s = ContinuousSet(initialize=[4, 5])
        m.fs = FlowsheetBlock(dynamic=True, time=m.s, time_units=units.s)

        assert m.fs.time_units is units.s


class TestSubFlowsheetBuild(object):
    # Test that build method works with nested flowsheets
    @pytest.mark.unit
    def test_flowsheet(self):
        # With nested flowsheet, flowsheet should return parent
        m = ConcreteModel()
        m.fs = FlowsheetBlock()
        m.fs.sub = FlowsheetBlock()

        assert m.fs.sub.flowsheet() is m.fs

    @pytest.mark.unit
    def test_default(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock()
        m.fs.sub = FlowsheetBlock()

        assert m.fs.sub.config.dynamic is False
        assert m.fs.sub.config.time is m.fs.config.time
        assert m.fs.sub.time_units is None

    @pytest.mark.unit
    def test_parent_dynamic_inherit(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)
        m.fs.sub = FlowsheetBlock()

        assert m.fs.sub.config.dynamic is True
        assert m.fs.sub.config.time is m.fs.config.time
        assert m.fs.sub.time_units == units.s

    @pytest.mark.unit
    def test_both_dynamic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)
        m.fs.sub = FlowsheetBlock(dynamic=True)

        assert m.fs.sub.config.dynamic is True
        assert m.fs.sub.config.time is m.fs.config.time
        assert m.fs.sub.time_units == units.s

    @pytest.mark.unit
    def test_ss_in_dynamic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)
        m.fs.sub = FlowsheetBlock(dynamic=False)

        assert m.fs.sub.config.dynamic is False
        assert m.fs.sub.config.time is m.fs.config.time
        assert m.fs.sub.time_units == units.s

    @pytest.mark.unit
    def test_dynamic_in_ss(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(DynamicError):
            m.fs.sub = FlowsheetBlock(dynamic=True)

    @pytest.mark.unit
    def test_ss_external_time(self):
        m = ConcreteModel()
        m.s = Set(initialize=[4, 5])
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)
        m.fs.sub = FlowsheetBlock(dynamic=False, time=m.s)

        assert m.fs.sub.config.dynamic is False
        assert m.fs.sub.config.time is m.s
        assert m.fs.sub.time_units is None

    @pytest.mark.unit
    def test_dynamic_external_time(self):
        m = ConcreteModel()
        m.s = ContinuousSet(initialize=[4, 5])
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)
        m.fs.sub = FlowsheetBlock(dynamic=True, time=m.s, time_units=units.minute)

        assert m.fs.sub.config.dynamic is True
        assert m.fs.sub.config.time is m.s
        assert m.fs.sub.time_units == units.minute

    @pytest.mark.unit
    def test_dynamic_external_time_invalid(self):
        m = ConcreteModel()
        m.s = Set(initialize=[4, 5])
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)
        with pytest.raises(DynamicError):
            m.fs.sub = FlowsheetBlock(dynamic=True, time=m.s, time_units=units.minute)

    @pytest.mark.unit
    def testtime_units_inherit(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)
        # Set differnt time units here to make sure they are ignored
        m.fs.sub = FlowsheetBlock(time_units=units.min)

        # Time should come from parent, not local settings
        assert m.fs.sub.time_units is units.s


class TestOther(object):
    @pytest.mark.unit
    def test_model_checks(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = PhysicalParameterTestBlock()
        m.fs.config.default_property_package = m.fs.props

        m.fs.unit1 = UnitModelBlock()

        m.fs.model_check()


class TestVisualisation(object):
    @pytest.mark.unit
    def test_report_empty(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.report()

    @pytest.mark.unit
    def test_get_stream_table_contents(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = PhysicalParameterTestBlock()
        m.fs.config.default_property_package = m.fs.props

        m.fs.unit1 = Heater()
        m.fs.unit2 = Heater()

        m.fs.stream = Arc(source=m.fs.unit1.outlet, destination=m.fs.unit2.inlet)
        TransformationFactory("network.expand_arcs").apply_to(m)

        df = m.fs._get_stream_table_contents()

        assert df.loc["pressure"]["stream"] == 1e5
        assert df.loc["temperature"]["stream"] == 300
        assert df.loc["component_flow_phase ('p1', 'c1')"]["stream"] == 2.0
        assert df.loc["component_flow_phase ('p1', 'c2')"]["stream"] == 2.0
        assert df.loc["component_flow_phase ('p2', 'c1')"]["stream"] == 2.0
        assert df.loc["component_flow_phase ('p2', 'c2')"]["stream"] == 2.0

        m.fs.report()
