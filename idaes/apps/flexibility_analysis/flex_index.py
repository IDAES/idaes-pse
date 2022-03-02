from pyomo.core.base.block import _BlockData
import pyomo.environ as pe
from typing import Sequence, Union, Mapping, MutableMapping, Optional
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.base.param import _ParamData
from .sampling import SamplingStrategy
from .flextest import (
    FlexTestConfig,
    solve_flextest,
    FlexTestMethod,
    FlexTestTermination,
    FlexTest,
)
import math
import logging


logger = logging.getLogger(__name__)


def _get_param_bounds(orig_param_bounds, nominal_values, flex_index):
    tmp_param_bounds = pe.ComponentMap()
    for p, (p_lb, p_ub) in orig_param_bounds.items():
        p_nom = nominal_values[p]
        tmp_param_bounds[p] = (
            p_nom - (p_nom - p_lb) * flex_index,
            p_nom + (p_ub - p_nom) * flex_index,
        )
    return tmp_param_bounds


def _add_table_row(
    log_level,
    outer_iter,
    flex_index_lower,
    flex_index_upper,
    fi_lb_max_viol,
    fi_ub_max_viol,
):
    if flex_index_lower is None:
        lb_str = str(flex_index_lower)
    else:
        lb_str = f"{flex_index_lower:12.3e}"
    if flex_index_upper is None:
        ub_str = str(flex_index_upper)
    else:
        ub_str = f"{flex_index_upper:12.3e}"
    if fi_lb_max_viol is None:
        lb_mv_str = str(fi_lb_max_viol)
    else:
        lb_mv_str = f"{fi_lb_max_viol:12.3e}"
    if fi_ub_max_viol is None:
        ub_mv_str = str(fi_ub_max_viol)
    else:
        ub_mv_str = f"{fi_ub_max_viol:12.3e}"
    logger.log(
        log_level,
        f"{outer_iter:<12}{lb_str:<12}{ub_str:<12}{lb_mv_str:<15}{ub_mv_str:<15}",
    )


def solve_flex_index(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]],
    controls: Sequence[_GeneralVarData],
    valid_var_bounds: MutableMapping[_GeneralVarData, Sequence[float]],
    in_place: bool = False,
    cap_index_at_1: bool = True,
    reconstruct_decision_rule: Optional[bool] = None,
    config: Optional[FlexTestConfig] = None,
    log_level: int = logging.INFO,
) -> float:
    original_uncertain_params = uncertain_params
    original_param_nominal_values = param_nominal_values
    original_param_bounds = param_bounds
    original_controls = controls
    original_valid_var_bounds = valid_var_bounds
    if not in_place:
        m = m.clone()
        uncertain_params = [m.find_component(i) for i in original_uncertain_params]
        param_nominal_values = pe.ComponentMap(
            (p, original_param_nominal_values[orig_p])
            for orig_p, p in zip(original_uncertain_params, uncertain_params)
        )
        param_bounds = pe.ComponentMap(
            (p, original_param_bounds[orig_p])
            for orig_p, p in zip(original_uncertain_params, uncertain_params)
        )
        controls = [m.find_component(i) for i in original_controls]
        valid_var_bounds = pe.ComponentMap(
            (m.find_component(v), bnds) for v, bnds in original_valid_var_bounds.items()
        )

    if (
        config.method == FlexTestMethod.linear_decision_rule
        and reconstruct_decision_rule is None
    ):
        reconstruct_decision_rule = True
    elif (
        config.method == FlexTestMethod.relu_decision_rule
        and reconstruct_decision_rule is None
    ):
        reconstruct_decision_rule = False
    elif reconstruct_decision_rule is None:
        reconstruct_decision_rule = False

    logger.log(
        log_level,
        f"{'Iter':<12}{'FI LB':<12}{'FI UB':<12}{'FI LB Max Viol':<15}{'FI UB Max Viol':<15}",
    )

    fi_lb_max_viol = -math.inf
    fi_ub_max_viol = math.inf

    outer_iter = 0
    _add_table_row(log_level, outer_iter, None, None, None, None)

    flex_index_lower = 0
    # make sure the nominal point is feasible
    nominal_bounds = pe.ComponentMap()
    for p, val in param_nominal_values.items():
        nominal_bounds[p] = (val, val)
    nominal_config: FlexTestConfig = config()
    nominal_config.method = FlexTestMethod.sampling
    nominal_config.terminate_early = True
    nominal_config.sampling_config.strategy = SamplingStrategy.grid
    nominal_config.sampling_config.num_points = 1
    nominal_config.sampling_config.enable_progress_bar = False
    nominal_res = solve_flextest(
        m=m,
        uncertain_params=uncertain_params,
        param_nominal_values=param_nominal_values,
        param_bounds=nominal_bounds,
        controls=controls,
        valid_var_bounds=valid_var_bounds,
        in_place=False,
        config=nominal_config,
    )
    if nominal_res.termination != FlexTestTermination.proven_feasible:
        raise RuntimeError("Nominal point is infeasible")

    outer_iter += 1
    fi_lb_max_viol = nominal_res.max_constraint_violation
    _add_table_row(
        log_level, outer_iter, flex_index_lower, None, fi_lb_max_viol, fi_ub_max_viol
    )

    flextest_config: FlexTestConfig = config()
    flextest_config.terminate_early = True
    flextest_config.sampling_config.enable_progress_bar = False
    flex_index_upper = 1

    if not cap_index_at_1:
        # Find an upper bound on the flexibility index (i.e., a point where the flextest fails)
        found_infeasible_point = False
        for _iter in range(10):
            tmp_param_bounds = _get_param_bounds(
                param_bounds, param_nominal_values, flex_index_upper
            )
            upper_res = solve_flextest(
                m=m,
                uncertain_params=uncertain_params,
                param_nominal_values=param_nominal_values,
                param_bounds=tmp_param_bounds,
                controls=controls,
                valid_var_bounds=valid_var_bounds,
                in_place=False,
                config=flextest_config,
            )
            outer_iter += 1
            if upper_res.termination == FlexTestTermination.found_infeasible_point:
                fi_ub_max_viol = upper_res.max_constraint_violation
                _add_table_row(
                    log_level,
                    outer_iter,
                    flex_index_lower,
                    flex_index_upper,
                    fi_lb_max_viol,
                    fi_ub_max_viol,
                )
                found_infeasible_point = True
                break
            elif upper_res.termination == FlexTestTermination.proven_feasible:
                flex_index_lower = flex_index_upper
                flex_index_upper *= 2
                fi_lb_max_viol = upper_res.max_constraint_violation
                _add_table_row(
                    log_level,
                    outer_iter,
                    flex_index_lower,
                    None,
                    fi_lb_max_viol,
                    fi_ub_max_viol,
                )
            else:
                raise RuntimeError("Unexpected termination from flexibility test")

        if not found_infeasible_point:
            raise RuntimeError("Could not find an upper bound on the flexibility index")

    max_param_bounds = _get_param_bounds(
        param_bounds, param_nominal_values, flex_index_upper
    )
    if cap_index_at_1:
        if reconstruct_decision_rule:
            res = solve_flextest(
                m=m,
                uncertain_params=uncertain_params,
                param_nominal_values=param_nominal_values,
                param_bounds=max_param_bounds,
                controls=controls,
                valid_var_bounds=valid_var_bounds,
                in_place=False,
                config=flextest_config,
            )
        else:
            ft = FlexTest(
                m=m,
                uncertain_params=uncertain_params,
                param_nominal_values=param_nominal_values,
                max_param_bounds=max_param_bounds,
                controls=controls,
                valid_var_bounds=valid_var_bounds,
                config=flextest_config,
            )
            res = ft.solve(max_param_bounds)
        outer_iter += 1
        if res.termination == FlexTestTermination.proven_feasible:
            flex_index_lower = 1
            fi_lb_max_viol = res.max_constraint_violation
            fi_ub_max_viol = res.max_constraint_violation
        elif res.termination == FlexTestTermination.found_infeasible_point:
            fi_ub_max_viol = res.max_constraint_violation
        else:
            raise RuntimeError("Unexpected termination from flexibility test")
        _add_table_row(
            log_level,
            outer_iter,
            flex_index_lower,
            flex_index_upper,
            fi_lb_max_viol,
            fi_ub_max_viol,
        )
    elif not reconstruct_decision_rule:
        ft = FlexTest(
            m=m,
            uncertain_params=uncertain_params,
            param_nominal_values=param_nominal_values,
            max_param_bounds=max_param_bounds,
            controls=controls,
            valid_var_bounds=valid_var_bounds,
            config=flextest_config,
        )

    while (flex_index_upper - flex_index_lower) > 1e-3:
        midpoint = 0.5 * (flex_index_lower + flex_index_upper)
        tmp_param_bounds = _get_param_bounds(
            param_bounds, param_nominal_values, midpoint
        )
        if reconstruct_decision_rule:
            res = solve_flextest(
                m=m,
                uncertain_params=uncertain_params,
                param_nominal_values=param_nominal_values,
                param_bounds=tmp_param_bounds,
                controls=controls,
                valid_var_bounds=valid_var_bounds,
                in_place=False,
                config=flextest_config,
            )
        else:
            res = ft.solve(param_bounds=tmp_param_bounds)
        outer_iter += 1
        if res.termination == FlexTestTermination.proven_feasible:
            flex_index_lower = midpoint
            fi_lb_max_viol = res.max_constraint_violation
        elif res.termination == FlexTestTermination.found_infeasible_point:
            flex_index_upper = midpoint
            fi_ub_max_viol = res.max_constraint_violation
        else:
            raise RuntimeError("Unexpected termination from flexibility test")
        _add_table_row(
            log_level,
            outer_iter,
            flex_index_lower,
            flex_index_upper,
            fi_lb_max_viol,
            fi_ub_max_viol,
        )

    return flex_index_lower
