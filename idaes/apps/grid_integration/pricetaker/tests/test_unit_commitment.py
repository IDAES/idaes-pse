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

import pyomo.environ as pyo
import pytest

import idaes.apps.grid_integration.pricetaker.unit_commitment as uc


@pytest.mark.unit
def test_unit_commitment_data():
    """Tests the UnitCommitmentData class"""
    m = uc.UnitCommitmentData(
        blk_name="ngcc",
        commodity_name="power",
        op_range_lb=0.3,
    )
    assert m.config.startup_rate is None
    assert m.config.shutdown_rate is None
    assert m.config.rampup_rate is None
    assert m.config.rampdown_rate is None
    assert m.config.op_range_lb == pytest.approx(0.3)
    assert m.blk_name == "ngcc"
    assert m.commodity_name == "power"

    # Test the update method
    mdl = pyo.ConcreteModel()
    mdl.capacity = pyo.Var()

    m.update(
        startup_rate=0.4,
        rampup_rate=0.4,
        capacity=mdl.capacity,
    )
    assert m.config.startup_rate == pytest.approx(0.4)
    assert m.config.rampup_rate == pytest.approx(0.4)
    assert isinstance(m.config.capacity, pyo.Var)
    assert m.config.capacity.value is None

    # Test assertion error associated with invalid startup_rate
    with pytest.raises(
        uc.ConfigurationError,
        match=(
            "For commidity power in operational "
            "block ngcc, \n\tthe startup rate is less than "
            "the minimum stable operation value."
        ),
    ):
        m.update(startup_rate=0.2)

    # Test assertion error associated with invalid shutdown_rate
    with pytest.raises(
        uc.ConfigurationError,
        match=(
            "For commidity power in operational "
            "block ngcc, \n\tthe shutdown rate is less than "
            "the minimum stable operation value."
        ),
    ):
        m.update(startup_rate=0.4, shutdown_rate=0.2)

    # Test assertion error associated with incomplete data
    with pytest.raises(
        uc.ConfigurationError,
        match=(
            "Necessary arguments needed for the ramping constraints " "are missing."
        ),
    ):
        # rampdown_rate value is missing, so it should throw an error.
        m.assert_ramping_args_present()


@pytest.mark.unit
def test_startup_shutdown_constraints():
    """Tests the startup_shutdown_constraints function"""
    m = pyo.ConcreteModel()
    m.set_time = pyo.RangeSet(10)

    @m.Block(m.set_time)
    def op_blk(b, _):
        b.op_mode = pyo.Var(within=pyo.Binary)
        b.startup = pyo.Var(within=pyo.Binary)
        b.shutdown = pyo.Var(within=pyo.Binary)

    m.startup_shutdown = pyo.Block()

    # Append startup and shutdown constraints
    uc.startup_shutdown_constraints(
        blk=m.startup_shutdown,
        op_blocks=m.op_blk,
        install_unit=1,
        minimum_up_time=3,
        minimum_down_time=4,
        set_time=m.set_time,
    )

    ss = m.startup_shutdown

    # Tests both the existence of the constraints and
    # the number of constraints
    # pylint: disable=no-member
    assert len(ss.binary_relationship_con) == 9
    assert len(ss.minimum_up_time_con) == 8
    assert len(ss.minimum_down_time_con) == 7

    # Check constraint is added for correct indices
    assert 1 not in ss.binary_relationship_con
    assert 1 not in ss.minimum_up_time_con
    assert 2 not in ss.minimum_up_time_con
    assert 1 not in ss.minimum_down_time_con
    assert 2 not in ss.minimum_down_time_con
    assert 3 not in ss.minimum_down_time_con

    # Check if the constraints are implemented correctly
    con1_2 = (
        "op_blk[2].op_mode - op_blk[1].op_mode  ==  "
        "op_blk[2].startup - op_blk[2].shutdown"
    )
    con1_10 = (
        "op_blk[10].op_mode - op_blk[9].op_mode  ==  "
        "op_blk[10].startup - op_blk[10].shutdown"
    )
    assert con1_2 == str(ss.binary_relationship_con[2].expr)
    assert con1_10 == str(ss.binary_relationship_con[10].expr)

    con2_3 = (
        "op_blk[1].startup + op_blk[2].startup + "
        "op_blk[3].startup  <=  op_blk[3].op_mode"
    )
    con2_10 = (
        "op_blk[8].startup + op_blk[9].startup + "
        "op_blk[10].startup  <=  op_blk[10].op_mode"
    )
    assert con2_3 == str(ss.minimum_up_time_con[3].expr)
    assert con2_10 == str(ss.minimum_up_time_con[10].expr)

    con3_4 = (
        "op_blk[1].shutdown + op_blk[2].shutdown + "
        "op_blk[3].shutdown + op_blk[4].shutdown  <=  "
        "1 - op_blk[4].op_mode"
    )
    con3_10 = (
        "op_blk[7].shutdown + op_blk[8].shutdown + "
        "op_blk[9].shutdown + op_blk[10].shutdown  <=  "
        "1 - op_blk[10].op_mode"
    )
    assert con3_4 == str(ss.minimum_down_time_con[4].expr)
    assert con3_10 == str(ss.minimum_down_time_con[10].expr)


@pytest.mark.unit
def test_capacity_limits_with_float():
    """
    Tests the capacity_limits function with a constant
    for capacity variable
    """
    m = pyo.ConcreteModel()
    m.set_time = pyo.RangeSet(10)

    @m.Block(m.set_time)
    def op_blk(b, _):
        b.power = pyo.Var(within=pyo.NonNegativeReals)
        b.op_mode = pyo.Var(within=pyo.Binary)

    m.cap_limits = pyo.Block()
    uc_data = uc.UnitCommitmentData(
        blk_name="op_blk", commodity_name="power", op_range_lb=0.3, capacity=650.0
    )
    uc.capacity_limits(
        blk=m.cap_limits,
        op_blocks=m.op_blk,
        uc_data=uc_data,
        set_time=m.set_time,
    )

    cl = m.cap_limits
    # Check the existence (and the number) of constraints
    # pylint: disable=no-member
    assert len(cl.capacity_low_limit_con) == 10
    assert len(cl.capacity_high_limit_con) == 10

    # Check if the constraint is implemented correctly
    con1_1 = "195.0*op_blk[1].op_mode  <=  op_blk[1].power"
    con1_9 = "195.0*op_blk[9].op_mode  <=  op_blk[9].power"
    assert con1_1 == str(cl.capacity_low_limit_con[1].expr)
    assert con1_9 == str(cl.capacity_low_limit_con[9].expr)

    con2_1 = "op_blk[1].power  <=  650.0*op_blk[1].op_mode"
    con2_9 = "op_blk[9].power  <=  650.0*op_blk[9].op_mode"
    assert con2_1 == str(cl.capacity_high_limit_con[1].expr)
    assert con2_9 == str(cl.capacity_high_limit_con[9].expr)


@pytest.mark.unit
def test_capacity_limits_with_var():
    """
    Tests the capacity_limits function with a Pyomo Var
    for capacity variable
    """
    m = pyo.ConcreteModel()
    m.set_time = pyo.RangeSet(10)
    m.capacity = pyo.Var(bounds=(150, 600))

    @m.Block(m.set_time)
    def op_blk(b, _):
        b.power = pyo.Var(within=pyo.NonNegativeReals)
        b.op_mode = pyo.Var(within=pyo.Binary)

    m.cap_limits = pyo.Block()
    uc_data = uc.UnitCommitmentData(
        blk_name="op_blk", commodity_name="power", op_range_lb=0.3, capacity=m.capacity
    )
    uc.capacity_limits(
        blk=m.cap_limits,
        op_blocks=m.op_blk,
        uc_data=uc_data,
        set_time=m.set_time,
    )

    cl = m.cap_limits
    # Check if the constraint is implemented correctly
    # pylint: disable=no-member
    con1_1 = "0.3*capacity*op_blk[1].op_mode  <=  op_blk[1].power"
    con1_9 = "0.3*capacity*op_blk[9].op_mode  <=  op_blk[9].power"
    assert con1_1 == str(cl.capacity_low_limit_con[1].expr)
    assert con1_9 == str(cl.capacity_low_limit_con[9].expr)

    con2_1 = "op_blk[1].power  <=  capacity*op_blk[1].op_mode"
    con2_9 = "op_blk[9].power  <=  capacity*op_blk[9].op_mode"
    assert con2_1 == str(cl.capacity_high_limit_con[1].expr)
    assert con2_9 == str(cl.capacity_high_limit_con[9].expr)


@pytest.mark.unit
def test_ramping_limits():
    """Tests the ramping_limits constraint"""
    m = pyo.ConcreteModel()
    m.set_time = pyo.RangeSet(10)
    m.capacity = pyo.Var(bounds=(150, 600))

    @m.Block(m.set_time)
    def op_blk(b, _):
        b.power = pyo.Var(within=pyo.NonNegativeReals)
        b.op_mode = pyo.Var(within=pyo.Binary)
        b.startup = pyo.Var(within=pyo.Binary)
        b.shutdown = pyo.Var(within=pyo.Binary)

    m.ramp_limits = pyo.Block()
    uc_data = uc.UnitCommitmentData(
        blk_name="op_blk",
        commodity_name="power",
        capacity=m.capacity,
        startup_rate=0.3,
        shutdown_rate=0.3,
        rampup_rate=0.2,
        rampdown_rate=0.15,
    )
    uc.ramping_limits(
        blk=m.ramp_limits,
        op_blocks=m.op_blk,
        uc_data=uc_data,
        set_time=m.set_time,
    )

    rl = m.ramp_limits
    # Check the existence (and the number) of constraints
    # pylint: disable=no-member
    assert 1 not in rl.ramp_up_con
    assert 1 not in rl.ramp_down_con
    assert len(rl.ramp_up_con) == 9
    assert len(rl.ramp_down_con) == 9

    # Check if the constraint is implemented correctly
    con1_2 = (
        "op_blk[2].power - op_blk[1].power  <=  "
        "0.3*capacity*op_blk[2].startup + "
        "0.2*capacity*op_blk[1].op_mode"
    )
    con1_10 = (
        "op_blk[10].power - op_blk[9].power  <=  "
        "0.3*capacity*op_blk[10].startup + "
        "0.2*capacity*op_blk[9].op_mode"
    )
    assert con1_2 == str(rl.ramp_up_con[2].expr)
    assert con1_10 == str(rl.ramp_up_con[10].expr)

    con2_2 = (
        "op_blk[1].power - op_blk[2].power  <=  "
        "0.3*capacity*op_blk[2].shutdown + "
        "0.15*capacity*op_blk[2].op_mode"
    )
    con2_10 = (
        "op_blk[9].power - op_blk[10].power  <=  "
        "0.3*capacity*op_blk[10].shutdown + "
        "0.15*capacity*op_blk[10].op_mode"
    )
    assert con2_2 == str(rl.ramp_down_con[2].expr)
    assert con2_10 == str(rl.ramp_down_con[10].expr)
