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
"""
Tests for the PIDController

Author: Douglas Allan
"""
import pytest
from io import StringIO

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    Param,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from pyomo.dae import DerivativeVar

from idaes.core import FlowsheetBlock, EnergyBalanceType, MomentumBalanceType
from idaes.models.control import (
    PIDController,
    ControllerType,
    ControllerMVBoundType,
    ControllerAntiwindupType,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    fixed_variables_set,
    activated_constraints_set,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError
from idaes.models_extra.power_generation.unit_models.soc_submodels.testing import (
    _build_test_utility,
)


def common_components(nt):
    return {
        Var: {
            "setpoint": nt,
            "process_var": nt,
            "manipulated_var": nt,
            "mv_ref": nt,
            "gain_p": nt,
        },
        Constraint: {"mv_eqn": nt},
        Expression: {"mv_unbounded": nt},
        Param: {},
    }


def add_integral_components(comp_dict, nt):
    comp_dict[Var] = {
        **comp_dict[Var],
        **{
            "gain_i": nt,
            "mv_integral_component": nt,
            "mv_integral_component_dot": nt,
        },
    }
    comp_dict[Constraint] = {
        **comp_dict[Constraint],
        **{
            "mv_integration_eqn": nt,
        },
    }


def add_derivative_components(comp_dict, nt):
    comp_dict[Var]["gain_d"] = nt
    comp_dict[Var]["derivative_term"] = nt


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_not_dynamic_exception():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_units=pyunits.dimensionless)

    m.fs.process_var = Var(m.fs.time, units=pyunits.dimensionless)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.dimensionless)

    with pytest.raises(
        ConfigurationError, match="PIDControllers work only with dynamic flowsheets."
    ):
        m.fs.unit = PIDController(
            process_var=m.fs.process_var, manipulated_var=m.fs.manipulated_var
        )


@pytest.mark.unit
def test_pi_control():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 1], time_units=pyunits.s)

    nt = 2  # Two timepoints

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var, manipulated_var=m.fs.manipulated_var
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    add_integral_components(comp_dict, nt)
    comp_dict[Expression]["error"] = nt
    comp_dict[Constraint]["initial_integral_error_eqn"] = 1

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )
    assert unit.mv_integration_eqn[0].active is False


@pytest.mark.unit
def test_p_control():
    m = ConcreteModel()
    nt = 3  # Three timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.P,
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    comp_dict[Expression]["error"] = nt

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )


@pytest.mark.unit
def test_pd_control():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PD,
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    add_derivative_components(comp_dict, nt)
    comp_dict[Expression]["error"] = nt
    comp_dict[Var]["negative_pv"] = nt
    comp_dict[Constraint]["negative_pv_eqn"] = nt

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )


@pytest.mark.unit
def test_pid_control():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PID,
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    add_integral_components(comp_dict, nt)
    add_derivative_components(comp_dict, nt)
    comp_dict[Expression]["error"] = nt
    comp_dict[Constraint]["initial_integral_error_eqn"] = 1
    comp_dict[Var]["negative_pv"] = nt
    comp_dict[Constraint]["negative_pv_eqn"] = nt

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )


@pytest.mark.component
def test_pid_units_consistent():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PID,
    )
    assert_units_consistent(m)
    assert_units_equivalent(m.fs.unit.setpoint, pyunits.m)
    assert_units_equivalent(m.fs.unit.mv_ref, pyunits.K)


@pytest.mark.unit
def test_derivative_on_error():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PID,
        derivative_on_error=True,
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    add_integral_components(comp_dict, nt)
    add_derivative_components(comp_dict, nt)
    comp_dict[Var]["error"] = nt
    comp_dict[Constraint]["initial_integral_error_eqn"] = 1
    comp_dict[Constraint]["error_eqn"] = nt

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )


@pytest.mark.component
def test_derivative_on_error_units_consistent():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PID,
        derivative_on_error=True,
    )

    assert_units_consistent(m)
    assert_units_equivalent(m.fs.unit.setpoint, pyunits.m)
    assert_units_equivalent(m.fs.unit.mv_ref, pyunits.K)


@pytest.mark.unit
def test_logistic_mv_bound():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PID,
        mv_bound_type=ControllerMVBoundType.LOGISTIC,
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    add_integral_components(comp_dict, nt)
    add_derivative_components(comp_dict, nt)
    comp_dict[Expression]["error"] = nt
    comp_dict[Constraint]["initial_integral_error_eqn"] = 1
    comp_dict[Var]["negative_pv"] = nt
    comp_dict[Constraint]["negative_pv_eqn"] = nt
    comp_dict[Param] = {"mv_lb": 1, "mv_ub": 1, "logistic_bound_k": 1}

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )


@pytest.mark.component
def test_logistic_mv_bound_units_consistent():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PID,
        mv_bound_type=ControllerMVBoundType.LOGISTIC,
    )
    assert_units_equivalent(m.fs.unit.setpoint, pyunits.m)
    assert_units_equivalent(m.fs.unit.mv_ref, pyunits.K)
    assert_units_consistent(m)


@pytest.mark.unit
def test_smooth_mv_bound():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PD,
        mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    add_derivative_components(comp_dict, nt)
    comp_dict[Expression]["error"] = nt
    comp_dict[Var]["negative_pv"] = nt
    comp_dict[Constraint]["negative_pv_eqn"] = nt
    comp_dict[Param] = {"mv_lb": 1, "mv_ub": 1, "smooth_eps": 1}

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )


@pytest.mark.component
def test_smooth_mv_bound_units_consistent():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PD,
        mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
    )
    assert_units_equivalent(m.fs.unit.setpoint, pyunits.m)
    assert_units_equivalent(m.fs.unit.mv_ref, pyunits.K)
    assert_units_consistent(m)


@pytest.mark.unit
def test_conditional_integration():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PI,
        mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
        antiwindup_type=ControllerAntiwindupType.CONDITIONAL_INTEGRATION,
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    add_integral_components(comp_dict, nt)
    comp_dict[Expression]["error"] = nt
    comp_dict[Constraint]["initial_integral_error_eqn"] = 1
    comp_dict[Param] = {
        "mv_lb": 1,
        "mv_ub": 1,
        "smooth_eps": 1,
        "conditional_integration_k": 1,
    }

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )


@pytest.mark.component
def test_conditional_integration_units_consistent():
    m = ConcreteModel()
    nt = 5  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PI,
        mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
        antiwindup_type=ControllerAntiwindupType.CONDITIONAL_INTEGRATION,
    )
    assert_units_equivalent(m.fs.unit.setpoint, pyunits.m)
    assert_units_equivalent(m.fs.unit.mv_ref, pyunits.K)
    assert_units_consistent(m)


@pytest.mark.unit
def test_back_calculation():
    m = ConcreteModel()
    nt = 11  # Eleven timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PID,
        mv_bound_type=ControllerMVBoundType.LOGISTIC,
        antiwindup_type=ControllerAntiwindupType.BACK_CALCULATION,
    )

    # Check unit config arguments
    assert len(unit.config) == 9

    assert unit.config.dynamic
    assert unit.config.has_holdup

    assert unit.process_var.referent is m.fs.process_var
    assert unit.manipulated_var.referent is m.fs.manipulated_var

    comp_dict = common_components(nt)
    add_integral_components(comp_dict, nt)
    add_derivative_components(comp_dict, nt)
    comp_dict[Expression]["error"] = nt
    comp_dict[Constraint]["initial_integral_error_eqn"] = 1
    comp_dict[Var]["negative_pv"] = nt
    comp_dict[Constraint]["negative_pv_eqn"] = nt
    comp_dict[Var]["gain_b"] = nt
    comp_dict[Param] = {"mv_lb": 1, "mv_ub": 1, "logistic_bound_k": 1}

    _build_test_utility(
        block=m.fs.unit,
        comp_dict=comp_dict,
        references=["process_var", "manipulated_var"],
    )


@pytest.mark.component
def test_back_calculation_units_consistent():
    m = ConcreteModel()
    nt = 11  # Eleven timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    unit = m.fs.unit = PIDController(
        process_var=m.fs.process_var,
        manipulated_var=m.fs.manipulated_var,
        controller_type=ControllerType.PID,
        mv_bound_type=ControllerMVBoundType.LOGISTIC,
        antiwindup_type=ControllerAntiwindupType.BACK_CALCULATION,
    )
    assert_units_equivalent(m.fs.unit.setpoint, pyunits.m)
    assert_units_equivalent(m.fs.unit.mv_ref, pyunits.K)
    assert_units_consistent(m)


@pytest.mark.unit
def test_antiwindup_no_bounds():
    m = ConcreteModel()
    nt = 11  # Eleven timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    with pytest.raises(
        ConfigurationError, match="User specified antiwindup method for unbounded MV."
    ):
        unit = m.fs.unit = PIDController(
            process_var=m.fs.process_var,
            manipulated_var=m.fs.manipulated_var,
            controller_type=ControllerType.PID,
            mv_bound_type=ControllerMVBoundType.NONE,
            antiwindup_type=ControllerAntiwindupType.BACK_CALCULATION,
        )


@pytest.mark.unit
def test_antiwindup_no_integration():
    m = ConcreteModel()
    nt = 7  # Five timepoints
    m.fs = FlowsheetBlock(dynamic=True, time_set=range(nt), time_units=pyunits.s)

    m.fs.process_var = Var(m.fs.time, units=pyunits.m)
    m.fs.manipulated_var = Var(m.fs.time, units=pyunits.K)

    with pytest.raises(
        ConfigurationError,
        match="User specified antiwindup method for controller without integral action.",
    ):
        unit = m.fs.unit = PIDController(
            process_var=m.fs.process_var,
            manipulated_var=m.fs.manipulated_var,
            controller_type=ControllerType.PD,
            mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
            antiwindup_type=ControllerAntiwindupType.BACK_CALCULATION,
        )


# if __name__ == "__main__":
#     test_config()
