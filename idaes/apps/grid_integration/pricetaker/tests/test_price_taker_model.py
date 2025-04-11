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

import pandas as pd
import pyomo.environ as pyo
import pytest

import idaes.logger as idaeslog
from idaes.apps.grid_integration import DesignModel, OperationModel
from idaes.apps.grid_integration.pricetaker.price_taker_model import PriceTakerModel
from idaes.apps.grid_integration.pricetaker.unit_commitment import UnitCommitmentData

# pylint: disable = unused-import
from idaes.apps.grid_integration.pricetaker.tests.test_clustering import (
    dummy_data_fixture,
    sklearn_avail,
)
from idaes.core.util.exceptions import ConfigurationError

pytest.importorskip("sklearn", reason="sklearn not available")


def simple_flowsheet_func(m):
    """Dummy flowsheet function for testing"""
    m.x = pyo.Var()
    m.y = pyo.Var()
    m.con1 = pyo.Constraint(expr=m.x + m.y == 1)
    m.LMP = pyo.Param(initialize=0.0, mutable=True)

    m.blk = pyo.Block()
    m.blk.power = pyo.Var()
    m.blk.op_mode = pyo.Var(within=pyo.Binary)
    m.blk.LMP = pyo.Param(initialize=0.0, mutable=True)


def foo_design_model(m, max_power, min_power):
    """Dummy design model"""
    m.PMAX = pyo.Param(initialize=max_power, mutable=False)
    m.PMIN = pyo.Param(initialize=min_power, mutable=False)

    m.capacity = pyo.Var()

    # DesignModel automatically declares install_unit variables
    m.low_capacity_limit = pyo.Constraint(expr=m.PMIN * m.install_unit <= m.capacity)
    m.up_capacity_limit = pyo.Constraint(expr=m.capacity <= m.PMAX * m.install_unit)

    # Capital investment cost ($)
    m.capex = 1000.0
    # Fixed operating and investment ($ / year)
    m.fom = 30.0


def foo_operation_model(m, design_blk):
    """Dummy operation model"""
    # Operation Variables
    m.power = pyo.Var(domain=pyo.NonNegativeReals, bounds=(0, design_blk.PMAX))

    # NOTE: OperationModel automatically declares op_mode,
    # startup and shutdown binary variables, and LMP Parameter
    # Surrogate cost constraint (combined fuel_cost and non_fuel_vom)
    m.fuel_cost = pyo.Expression(expr=23 * m.power + 50 * m.op_mode)
    m.elec_revenue = pyo.Expression(expr=m.power * m.LMP)
    m.su_sd_costs = pyo.Expression(expr=5 * m.startup + 4 * m.shutdown)


def build_foo_flowsheet(m, design_blk):
    """Dummy flowsheet model for testing"""
    m.op_blk = OperationModel(
        model_func=foo_operation_model,
        model_args={"design_blk": design_blk},
    )

    m.fixed_rev = pyo.Expression(expr=100)
    m.fixed_cost = pyo.Expression(expr=50)


@pytest.mark.unit
def test_num_representative_days():
    """Tests the num_representative_days property"""
    m = PriceTakerModel()

    # By default, num_representative_days should be None
    assert m.num_representative_days is None

    # Assign number of representative days
    m.num_representative_days = 20
    assert m.num_representative_days == 20

    # Test overwrite error
    with pytest.raises(
        ConfigurationError,
        match=(
            "num_representative_days is already defined as 20 "
            "and it cannot be overwritten."
            "\n\tInstantiate a new PriceTakerModel object."
        ),
    ):
        m.num_representative_days = 25

    # Ensure that the value remained the same
    assert m.num_representative_days == 20

    # Test error if an invalid value is specified
    m = PriceTakerModel()
    with pytest.raises(ValueError):
        m.num_representative_days = -0.5


@pytest.mark.unit
def test_horizon_length(caplog):
    """Tests the horizon_length property"""
    m = PriceTakerModel()

    # Test default value warning is not printed, it it is assigned
    with caplog.at_level(idaeslog.WARNING):
        m.horizon_length = 48
        assert m.horizon_length == 48
        assert (
            "Attribute horizon_length is not specified. Using 24 hours "
            "as the horizon length."
        ) not in caplog.text

    # Test overwrite error
    with pytest.raises(
        ConfigurationError,
        match=(
            "horizon_length is already defined as 48 and it cannot be "
            "overwritten.\n\tInstantiate a new PriceTakerModel object."
        ),
    ):
        m.horizon_length = 72

    # Ensure that the value remained the same
    assert m.horizon_length == 48

    # Test the default value warning
    m = PriceTakerModel()
    with caplog.at_level(idaeslog.WARNING):
        assert m.horizon_length == 24
        assert (
            "Attribute horizon_length is not specified. Using 24 hours "
            "as the horizon length."
        ) in caplog.text

    # Test error message associated with an infeasible value
    with pytest.raises(ValueError):
        m.horizon_length = -4


@pytest.mark.unit
def test_append_lmp_data():
    "Tests the append_lmp_data method"
    m = PriceTakerModel()
    data = [-4, 5.0, 2, 50, -3.2]

    m.append_lmp_data(
        lmp_data=pd.Series(data),
        num_representative_days=3,
        horizon_length=2,
        seed=30,
    )

    # pylint: disable = protected-access
    assert m._config.lmp_data == data
    assert m.num_representative_days == 3
    assert m.horizon_length == 2
    assert m._config.seed == 30

    # Test the error associated with overwriting LMP data
    with pytest.raises(
        ConfigurationError,
        match=(
            "Attempted to overwrite the LMP data. Instantiate a "
            "new PriceTakerModel object to change the LMP data."
        ),
    ):
        m.append_lmp_data(lmp_data=data + [4, 7, 8])


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_get_optimal_representative_days(dummy_data):
    """Tests the get_optimal_representative_days method"""
    m = PriceTakerModel()

    # Test lmp data not available error
    with pytest.raises(
        ConfigurationError,
        match=(
            "LMP data is missing. Please append the LMP data using the "
            "`append_lmp_data` method."
        ),
    ):
        m.get_optimal_representative_days()

    # Append LMP data
    m.append_lmp_data(lmp_data=dummy_data, horizon_length=2)

    # Test invalid kmin and kmax error
    with pytest.raises(ValueError):
        m.get_optimal_representative_days(kmin=-5, kmax=-3.5)

    # Test if the method works
    assert 3 == m.get_optimal_representative_days(
        kmin=2, kmax=6, method="silhouette", generate_elbow_plot=False
    )


@pytest.mark.unit
def test_mp_model_full_year(dummy_data):
    """
    Tests the build_multiperiod_model using the full year price signal
    """
    m = PriceTakerModel()

    # Test lmp data not available error
    with pytest.raises(
        ConfigurationError,
        match=(
            "LMP data is missing. Please append the LMP data using the "
            "`append_lmp_data` method."
        ),
    ):
        m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    # Build the multiperiod model
    m.append_lmp_data(lmp_data=dummy_data)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    # Test if the model construction is successful
    assert m.horizon_length == len(dummy_data)
    assert m.num_representative_days == 1
    assert len(m.period) == len(dummy_data)
    assert isinstance(m.set_days, pyo.RangeSet)
    assert isinstance(m.set_time, pyo.RangeSet)
    assert len(m.set_days) == 1
    assert len(m.set_time) == len(dummy_data)
    assert m.rep_days_weights == {1: 1}
    assert m.rep_days_lmp == {1: {t + 1: val for t, val in enumerate(dummy_data)}}
    for d, t in m.period:
        assert isinstance(m.period[d, t].x, pyo.Var)
        assert isinstance(m.period[d, t].y, pyo.Var)
        assert isinstance(m.period[d, t].con1, pyo.Constraint)

        # Test if the LMP data is updated at flowsheet level
        assert m.period[d, t].LMP.value == pytest.approx(dummy_data[t - 1])

        # Test if the LMP data is updated for each operational block
        assert m.period[d, t].blk.LMP.value == pytest.approx(dummy_data[t - 1])

    # Test model overwrite error
    with pytest.raises(
        ConfigurationError,
        match=(
            "A multiperiod model might already exist, as the object has "
            "`period` attribute."
        ),
    ):
        m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_mp_model_rep_days(dummy_data):
    """
    Tests the build_multiperiod_model using representative days
    """
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data, horizon_length=2, num_representative_days=3)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    # Test if the model construction is successful
    assert m.horizon_length == 2
    assert m.num_representative_days == 3
    assert len(m.period) == 2 * 3
    assert len(m.set_days) == 3
    assert len(m.set_time) == 2
    assert m.rep_days_weights == {1: 4, 2: 4, 3: 4}

    rep_lmp_data = {
        1: {1: 10.5, 2: 0.5},
        2: {1: 5.5, 2: 10.5},
        3: {1: 0.5, 2: 0.5},
    }
    assert m.rep_days_lmp == rep_lmp_data
    for d, t in m.period:
        assert isinstance(m.period[d, t].x, pyo.Var)
        assert isinstance(m.period[d, t].y, pyo.Var)
        assert isinstance(m.period[d, t].con1, pyo.Constraint)

        # Test if the LMP data is updated at flowsheet level
        assert m.period[d, t].LMP.value == pytest.approx(rep_lmp_data[d][t])

        # Test if the LMP data is updated for each operational block
        assert m.period[d, t].blk.LMP.value == pytest.approx(rep_lmp_data[d][t])


@pytest.mark.unit
def test_get_operation_blocks(dummy_data):
    """Tests the _get_operation_blocks method"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)

    # Test multiperiod model does not exist error
    with pytest.raises(
        ConfigurationError,
        match=(
            "Unable to find the multiperiod model. Please use the "
            "build_multiperiod_model method to construct one."
        ),
    ):
        # pylint: disable = protected-access
        m._get_operation_blocks(blk_name="blk", attribute_list=["op_mode"])

    # Build multiperiod model
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    # pylint: disable = protected-access
    op_blocks = m._get_operation_blocks(blk_name="blk", attribute_list=["op_mode"])

    assert len(op_blocks) == 1
    assert len(op_blocks[1]) == 24

    for d, t in m.period:
        assert op_blocks[d][t] is m.period[d, t].blk

    # Test operational block not found error
    with pytest.raises(
        AttributeError,
        match="Operational block foo_blk does not exist.",
    ):
        m._get_operation_blocks(blk_name="foo_blk", attribute_list=["op_mode"])

    # Test missing attribute error
    with pytest.raises(
        AttributeError,
        match=(
            "Required attribute startup is not found in " "the operational block blk."
        ),
    ):
        m._get_operation_blocks(blk_name="blk", attribute_list=["op_mode", "startup"])


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_get_operation_blocks_rep_days(dummy_data):
    """Tests the get_operation_blocks method with representative days"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data, horizon_length=2, num_representative_days=3)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    # pylint: disable = protected-access
    op_blocks = m._get_operation_blocks(blk_name="blk", attribute_list=["op_mode"])

    assert len(op_blocks) == 3  # 3 representative days
    for d in m.rep_days_weights:
        assert len(op_blocks[d]) == 2  # horizon length = 2

    for d, t in m.period:
        assert op_blocks[d][t] is m.period[d, t].blk


@pytest.mark.unit
def test_get_operation_vars(dummy_data):
    """Tests the _get_operation_vars method"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)

    # Test multiperiod model does not exist error
    with pytest.raises(
        ConfigurationError,
        match=(
            "Unable to find the multiperiod model. Please use the "
            "build_multiperiod_model method to construct one."
        ),
    ):
        # pylint: disable = protected-access
        m._get_operation_vars(var_name="blk.power")

    # Build multiperiod model
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    # pylint: disable = protected-access
    op_vars = m._get_operation_vars(var_name="blk.power")

    assert len(op_vars) == 1
    assert len(op_vars[1]) == 24

    for d, t in m.period:
        assert op_vars[d][t] is m.period[d, t].blk.power

    # Test operational block not found error
    with pytest.raises(
        AttributeError,
        match="Variable foo_blk.foo does not exist in the multiperiod model.",
    ):
        m._get_operation_vars(var_name="foo_blk.foo")


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_get_operation_vars_rep_days(dummy_data):
    """Tests the get_operation_vars method with representative days"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data, horizon_length=2, num_representative_days=3)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    # pylint: disable = protected-access
    op_vars = m._get_operation_vars(var_name="blk.power")

    assert len(op_vars) == 3  # 3 representative days
    for d in m.rep_days_weights:
        assert len(op_vars[d]) == 2  # horizon length = 2

    for d, t in m.period:
        assert op_vars[d][t] is m.period[d, t].blk.power


@pytest.mark.unit
def test_add_linking_constraints(dummy_data):
    """Tests the add_linking_constraints method"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)
    m.add_linking_constraints(
        previous_time_var="blk.power", current_time_var="blk.op_mode"
    )

    assert len(m.variable_linking_constraints_1) == 24 - 1
    assert str(m.variable_linking_constraints_1[1, 10].expr) == (
        "period[1,9].blk.power  ==  period[1,10].blk.op_mode"
    )
    assert str(m.variable_linking_constraints_1[1, 24].expr) == (
        "period[1,23].blk.power  ==  period[1,24].blk.op_mode"
    )

    # Add another set of linking constraints
    m.add_linking_constraints(previous_time_var="x", current_time_var="y")
    assert len(m.variable_linking_constraints_2) == 24 - 1
    assert str(m.variable_linking_constraints_2[1, 10].expr) == (
        "period[1,9].x  ==  period[1,10].y"
    )
    assert str(m.variable_linking_constraints_2[1, 24].expr) == (
        "period[1,23].x  ==  period[1,24].y"
    )


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_add_linking_constraints_rep_days(dummy_data):
    """Tests the add_linking_constraints method with representative days"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data, horizon_length=12, num_representative_days=2)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)
    m.add_linking_constraints(
        previous_time_var="blk.power", current_time_var="blk.op_mode"
    )

    assert len(m.variable_linking_constraints_1) == 2 * (12 - 1)
    assert str(m.variable_linking_constraints_1[1, 10].expr) == (
        "period[1,9].blk.power  ==  period[1,10].blk.op_mode"
    )
    assert str(m.variable_linking_constraints_1[2, 12].expr) == (
        "period[2,11].blk.power  ==  period[2,12].blk.op_mode"
    )


@pytest.mark.unit
def test_retrieve_uc_data():
    """Tests the _retrieve_uc_data method"""
    m = PriceTakerModel()

    # pylint: disable = protected-access
    assert len(m._op_blk_uc_data) == 0
    uc_data = m._retrieve_uc_data(op_blk="blk_1", commodity="power", op_range_lb=0.2)

    assert len(m._op_blk_uc_data) == 1
    assert uc_data is m._op_blk_uc_data["blk_1", "power"]
    assert isinstance(uc_data, UnitCommitmentData)

    # Test updating the data
    uc_data_2 = m._retrieve_uc_data(op_blk="blk_1", commodity="power", startup_rate=0.2)
    assert len(m._op_blk_uc_data) == 1
    # It should return the existing object.
    assert uc_data_2 is uc_data
    assert uc_data_2.config.startup_rate == pytest.approx(0.2)

    # Test invalid data error
    with pytest.raises(
        ConfigurationError,
        match=(
            "For commidity power in operational "
            "block blk_1, \n\tthe shutdown rate is less than "
            "the minimum stable operation value."
        ),
    ):
        m._retrieve_uc_data(op_blk="blk_1", commodity="power", shutdown_rate=0.1)

    # Test addition of different commodity to the same operational block
    uc_data_3 = m._retrieve_uc_data(
        op_blk="blk_1", commodity="ng_flow", startup_rate=0.2
    )
    assert len(m._op_blk_uc_data) == 2
    assert ("blk_1", "ng_flow") in m._op_blk_uc_data
    assert uc_data_3 is not uc_data
    assert uc_data_3 is not uc_data_2

    # Test addition of different operational block
    uc_data_4 = m._retrieve_uc_data(
        op_blk="blk_new", commodity="ng_flow", startup_rate=0.2
    )
    assert len(m._op_blk_uc_data) == 3
    assert ("blk_new", "ng_flow") in m._op_blk_uc_data
    assert uc_data_4 is not uc_data
    assert uc_data_4 is not uc_data_2
    assert uc_data_4 is not uc_data_3

    assert list(m._op_blk_uc_data.keys()) == [
        ("blk_1", "power"),
        ("blk_1", "ng_flow"),
        ("blk_new", "ng_flow"),
    ]


@pytest.mark.unit
def test_add_capacity_limits(dummy_data):
    """Tests the add_capacity_limits method"""

    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    method_args = {
        "op_block_name": "blk",
        "commodity": "power",
        "capacity": 100,
        "op_range_lb": 0.3,
    }

    # Build the multiperiod model
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    # Test operation block does not exist error
    with pytest.raises(
        AttributeError,
        match="Operational block foo_blk does not exist.",
    ):
        method_args_2 = method_args.copy()
        method_args_2["op_block_name"] = "foo_blk"
        m.add_capacity_limits(**method_args_2)

    # Test if capacity limits are added correctly
    m.add_capacity_limits(**method_args)

    assert len(m.blk_power_limits) == 1
    assert len(m.blk_power_limits[1].capacity_low_limit_con) == 24
    assert len(m.blk_power_limits[1].capacity_high_limit_con) == 24

    # Test overwriting error
    with pytest.raises(
        ConfigurationError,
        match="Attempting to overwrite capacity limits for power in blk.",
    ):
        m.add_capacity_limits(**method_args)


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_add_capacity_limits_rep_days(dummy_data):
    """Tests the add_capacity_limits method with representative days"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data, horizon_length=2, num_representative_days=3)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)
    m.add_capacity_limits(
        op_block_name="blk", commodity="power", capacity=100, op_range_lb=0.3
    )

    # Test if capacity limits are added correctly
    assert len(m.blk_power_limits) == 3  # 3 representative days
    # Following must be valid because the horizon length is 2
    for d in [1, 2, 3]:
        assert len(m.blk_power_limits[d].capacity_low_limit_con) == 2
        assert len(m.blk_power_limits[d].capacity_high_limit_con) == 2


@pytest.mark.unit
def test_add_ramping_limits(dummy_data):
    """Tests the add_ramping_limits method"""

    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.design_blk = DesignModel(
        model_func=foo_design_model, model_args={"max_power": 400, "min_power": 300}
    )
    method_args = {
        "op_block_name": "op_blk",
        "commodity": "power",
        "capacity": m.design_blk.capacity,
        "startup_rate": 0.3,
        "shutdown_rate": 0.4,
        "rampup_rate": 0.2,
        "rampdown_rate": 0.2,
    }

    # Build the multiperiod model
    m.build_multiperiod_model(
        flowsheet_func=build_foo_flowsheet,
        flowsheet_options={"design_blk": m.design_blk},
    )

    # Test necessary attributes does not exist error
    with pytest.raises(
        AttributeError,
        match="Operational block foo_blk does not exist.",
    ):
        method_args_2 = method_args.copy()
        method_args_2["op_block_name"] = "foo_blk"
        m.add_ramping_limits(**method_args_2)

    # Test if ramping limits are added correctly
    m.add_ramping_limits(**method_args)

    # pylint: disable = protected-access
    uc_data = m._op_blk_uc_data["op_blk", "power"].config
    assert uc_data.startup_rate == pytest.approx(method_args["startup_rate"])
    assert uc_data.shutdown_rate == pytest.approx(method_args["shutdown_rate"])
    assert uc_data.rampup_rate == pytest.approx(method_args["rampup_rate"])
    assert uc_data.rampdown_rate == pytest.approx(method_args["rampdown_rate"])
    assert uc_data.capacity is m.design_blk.capacity

    assert len(m.op_blk_power_ramping) == 1
    assert len(m.op_blk_power_ramping[1].ramp_up_con) == 24 - 1
    assert len(m.op_blk_power_ramping[1].ramp_down_con) == 24 - 1

    # Test overwriting error
    with pytest.raises(
        ConfigurationError,
        match="Attempting to overwrite ramping limits for power in op_blk.",
    ):
        m.add_ramping_limits(**method_args)


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_add_ramping_limits_rep_days(dummy_data):
    """Tests the add_ramping_limits method with representative days"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data, horizon_length=12, num_representative_days=2)
    m.design_blk = DesignModel(
        model_func=foo_design_model, model_args={"max_power": 400, "min_power": 300}
    )
    m.build_multiperiod_model(
        flowsheet_func=build_foo_flowsheet,
        flowsheet_options={"design_blk": m.design_blk},
    )
    m.add_ramping_limits(
        op_block_name="op_blk",
        commodity="power",
        capacity=m.design_blk.capacity,
        startup_rate=0.3,
        shutdown_rate=0.4,
        rampup_rate=0.2,
        rampdown_rate=0.2,
    )
    assert len(m.op_blk_power_ramping) == 2  # 2 representative days
    # Following must be valid because the horizon length is 12
    for d in [1, 2]:
        assert len(m.op_blk_power_ramping[d].ramp_up_con) == 12 - 1
        assert len(m.op_blk_power_ramping[d].ramp_down_con) == 12 - 1


@pytest.mark.unit
def test_start_shut_no_install(dummy_data):
    """Tests the add_startup_shutdown method without install_unit variable."""

    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.design_blk = DesignModel(
        model_func=foo_design_model, model_args={"max_power": 400, "min_power": 300}
    )
    # Build the multiperiod model
    m.build_multiperiod_model(
        flowsheet_func=build_foo_flowsheet,
        flowsheet_options={"design_blk": m.design_blk},
    )

    # Test if startup and shutdown constraints are added correctly
    m.add_startup_shutdown(
        op_block_name="op_blk", minimum_up_time=3, minimum_down_time=4
    )

    assert hasattr(m, "op_blk_startup_shutdown")
    assert hasattr(m.config, "minimum_up_time")
    assert hasattr(m.config, "minimum_down_time")
    assert len(m.op_blk_startup_shutdown) == 1
    assert len(m.op_blk_startup_shutdown[1].binary_relationship_con) == 24 - 1
    assert len(m.op_blk_startup_shutdown[1].minimum_up_time_con) == 24 - (3 - 1)
    assert len(m.op_blk_startup_shutdown[1].minimum_down_time_con) == 24 - (4 - 1)

    con3_1_10 = (
        "period[1,7].op_blk.shutdown + period[1,8].op_blk.shutdown + "
        "period[1,9].op_blk.shutdown + period[1,10].op_blk.shutdown  <=  "
        "1 - period[1,10].op_blk.op_mode"
    )
    assert con3_1_10 == str(m.op_blk_startup_shutdown[1].minimum_down_time_con[10].expr)

    # Test overwriting error
    with pytest.raises(
        ConfigurationError,
        match=(
            "Attempting to overwrite startup/shutdown constraints "
            "for operation block op_blk."
        ),
    ):
        m.add_startup_shutdown(
            op_block_name="op_blk", minimum_up_time=5, minimum_down_time=3
        )


@pytest.mark.unit
def test_start_shut_with_install(dummy_data):
    """Tests the add_startup_shutdown method with install_unit variable."""

    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.design_blk = DesignModel(
        model_func=foo_design_model, model_args={"max_power": 400, "min_power": 300}
    )
    # Build the multiperiod model
    m.build_multiperiod_model(
        flowsheet_func=build_foo_flowsheet,
        flowsheet_options={"design_blk": m.design_blk},
    )

    # Test the error associated with the absence of install_unit variable
    method_args = {
        "op_block_name": "op_blk",
        "des_block_name": "design_blk",
        "minimum_up_time": 3,
        "minimum_down_time": 4,
    }
    m.design_blk.del_component(m.design_blk.install_unit)
    assert m.find_component("design_blk.install_unit") is None

    with pytest.raises(
        AttributeError,
        match=(
            "Binary variable associated with unit installation is not found "
            "in design_blk. \n\tDo not specify des_block_name argument if "
            "installation of the unit is not a decision variable."
        ),
    ):
        m.add_startup_shutdown(**method_args)

    # Construct the constraints with install_unit variable
    m.design_blk.install_unit = pyo.Var(within=pyo.Binary)
    assert m.find_component("design_blk.install_unit") is not None
    m.add_startup_shutdown(**method_args)

    con3_1_10 = (
        "period[1,7].op_blk.shutdown + period[1,8].op_blk.shutdown + "
        "period[1,9].op_blk.shutdown + period[1,10].op_blk.shutdown  <=  "
        "design_blk.install_unit - period[1,10].op_blk.op_mode"
    )
    assert con3_1_10 == str(m.op_blk_startup_shutdown[1].minimum_down_time_con[10].expr)


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_start_shut_with_rep_days(dummy_data):
    """
    Tests the add_startup_shutdown method with install_unit variable,
    and with representative days
    """

    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data, horizon_length=12, num_representative_days=2)
    m.design_blk = DesignModel(
        model_func=foo_design_model, model_args={"max_power": 400, "min_power": 300}
    )
    m.build_multiperiod_model(
        flowsheet_func=build_foo_flowsheet,
        flowsheet_options={"design_blk": m.design_blk},
    )
    m.add_startup_shutdown(
        op_block_name="op_blk",
        des_block_name="design_blk",
        minimum_up_time=3,
        minimum_down_time=4,
    )

    assert len(m.op_blk_startup_shutdown) == 2  # 2 representative days
    # Following must be valid because the horizon length is 12
    for d in [1, 2]:
        assert len(m.op_blk_startup_shutdown[d].binary_relationship_con) == 12 - 1
        assert len(m.op_blk_startup_shutdown[d].minimum_up_time_con) == 12 - (3 - 1)
        assert len(m.op_blk_startup_shutdown[d].minimum_down_time_con) == 12 - (4 - 1)

    con3_1_10 = (
        "period[1,7].op_blk.shutdown + period[1,8].op_blk.shutdown + "
        "period[1,9].op_blk.shutdown + period[1,10].op_blk.shutdown  <=  "
        "design_blk.install_unit - period[1,10].op_blk.op_mode"
    )
    con3_2_10 = (
        "period[2,7].op_blk.shutdown + period[2,8].op_blk.shutdown + "
        "period[2,9].op_blk.shutdown + period[2,10].op_blk.shutdown  <=  "
        "design_blk.install_unit - period[2,10].op_blk.op_mode"
    )
    assert con3_1_10 == str(m.op_blk_startup_shutdown[1].minimum_down_time_con[10].expr)
    assert con3_2_10 == str(m.op_blk_startup_shutdown[2].minimum_down_time_con[10].expr)


@pytest.mark.unit
def test_add_hourly_cashflows_warnings(dummy_data, caplog):
    """Tests the add_hourly_cashflows method with empty args"""

    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.design_blk = DesignModel(
        model_func=foo_design_model, model_args={"max_power": 400, "min_power": 300}
    )
    # Build the multiperiod model
    m.build_multiperiod_model(
        flowsheet_func=build_foo_flowsheet,
        flowsheet_options={"design_blk": m.design_blk},
    )

    with caplog.at_level(idaeslog.WARNING):
        m.add_hourly_cashflows()

        assert (
            "Argument operational_costs is not specified, so the total "
            "operational cost will be set to 0."
        ) in caplog.text

        assert (
            "Argument revenue_streams is not specified, so the total "
            "revenue will be set to 0."
        ) in caplog.text

    for d, t in m.period:
        assert str(m.period[d, t].total_hourly_cost.expr) == "0.0"
        assert str(m.period[d, t].total_hourly_revenue.expr) == "0.0"
        assert str(m.period[d, t].net_hourly_cash_inflow.expr) == "0 - 0"


@pytest.mark.unit
def test_add_hourly_cashflows(dummy_data, caplog):
    """Tests the add_hourly_cashflows method"""

    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.design_blk = DesignModel(
        model_func=foo_design_model, model_args={"max_power": 400, "min_power": 300}
    )
    # Build the multiperiod model
    m.build_multiperiod_model(
        flowsheet_func=build_foo_flowsheet,
        flowsheet_options={"design_blk": m.design_blk},
    )

    with caplog.at_level(idaeslog.WARNING):
        m.add_hourly_cashflows(
            operational_costs=["fuel_cost", "su_sd_costs", "fixed_cost"],
            revenue_streams=["elec_revenue", "fixed_rev"],
        )

        assert (
            "Argument operational_costs is not specified, so the total "
            "operational cost will be set to 0."
        ) not in caplog.text

        assert (
            "Argument revenue_streams is not specified, so the total "
            "revenue will be set to 0."
        ) not in caplog.text

    for d, t in m.period:
        assert str(m.period[d, t].total_hourly_cost.expr) == (
            f"23*period[{d},{t}].op_blk.power + 50*period[{d},{t}].op_blk.op_mode"
            f" + (5*period[{d},{t}].op_blk.startup + 4*period[{d},{t}].op_blk.shutdown)"
            f" + 50"
        )
        assert str(m.period[d, t].total_hourly_revenue.expr) == (
            f"period[{d},{t}].op_blk.LMP*period[{d},{t}].op_blk.power + 100"
        )


@pytest.mark.unit
def test_add_overall_cashflows_warnings(dummy_data, caplog):
    """Tests the warning/error messages in add_overall_cashflows"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    with pytest.raises(
        ConfigurationError,
        match=(
            "Hourly cashflows are not added to the model. Please run "
            "add_hourly_cashflows method before calling the "
            "add_overall_cashflows method."
        ),
    ):
        m.add_overall_cashflows()

    with caplog.at_level(idaeslog.WARNING):
        m.add_hourly_cashflows()
        m.add_overall_cashflows()

        assert (
            "Argument operational_costs is not specified, so the total "
            "operational cost will be set to 0."
        ) in caplog.text

        assert (
            "Argument revenue_streams is not specified, so the total "
            "revenue will be set to 0."
        ) in caplog.text

        assert (
            "No design blocks were found, so the overall capital cost "
            "(capex) and the fixed O&M cost (fom) are set to 0."
        ) in caplog.text

    cf = m.cashflows
    assert str(cf.capex_calculation.expr) == "cashflows.capex  ==  0"
    assert str(cf.fom_calculation.expr) == "cashflows.fom  ==  0"
    assert str(cf.depreciation_calculation.expr) == (
        "cashflows.depreciation  ==  0.03333333333333333*cashflows.capex"
    )
    assert hasattr(cf, "net_cash_inflow")
    assert hasattr(cf, "net_cash_inflow_calculation")
    assert str(cf.corporate_tax_calculation.expr) == (
        "0.2*(cashflows.net_cash_inflow - cashflows.fom - cashflows.depreciation)"
        "  <=  cashflows.corporate_tax"
    )

    assert str(cf.net_profit_calculation.expr) == (
        "cashflows.net_profit  ==  cashflows.net_cash_inflow - cashflows.fom"
        " - cashflows.corporate_tax"
    )
    assert str(cf.npv.expr) == (
        "cashflows.net_profit - 0.08882743338727227*cashflows.capex"
    )
    assert str(cf.lifetime_npv.expr) == (
        "11.257783343127485*cashflows.net_profit - cashflows.capex"
    )


@pytest.mark.unit
def test_add_overall_cashflows(dummy_data):
    """Tests the add_overall_cashflows method"""

    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.design_blk = DesignModel(
        model_func=foo_design_model, model_args={"max_power": 400, "min_power": 300}
    )
    # Build the multiperiod model
    m.build_multiperiod_model(
        flowsheet_func=build_foo_flowsheet,
        flowsheet_options={"design_blk": m.design_blk},
    )
    m.add_hourly_cashflows(
        operational_costs=["fuel_cost", "su_sd_costs", "fixed_cost"],
        revenue_streams=["elec_revenue", "fixed_rev"],
    )
    m.add_overall_cashflows()
    cf = m.cashflows
    assert str(cf.capex_calculation.expr) == "cashflows.capex  ==  1000.0"
    assert str(cf.fom_calculation.expr) == "cashflows.fom  ==  30.0"


@pytest.mark.unit
def test_add_objective_function(dummy_data):
    """Tests the add_objective_function method"""
    m = PriceTakerModel()
    m.append_lmp_data(lmp_data=dummy_data)
    m.build_multiperiod_model(flowsheet_func=simple_flowsheet_func)

    with pytest.raises(
        ConfigurationError,
        match=(
            "Overall cashflows are not appended. Please run the "
            "add_overall_cashflows method."
        ),
    ):
        m.add_objective_function(objective_type="npv")

    # Test the error associated with unsupported objective function
    m.add_hourly_cashflows()
    m.add_overall_cashflows()

    with pytest.raises(
        ConfigurationError,
        match=(
            "foo is not a supported objective function."
            "Please specify either npv, or lifetime_npv, or net_profit "
            "as the objective_type."
        ),
    ):
        m.add_objective_function(objective_type="foo")

    # Test if the objective function is added
    m.add_objective_function(objective_type="net_profit")
    assert hasattr(m, "obj")
