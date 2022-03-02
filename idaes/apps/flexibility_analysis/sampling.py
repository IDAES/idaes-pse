from pyomo.core.base.block import _BlockData
import pyomo.environ as pe
import math
import numpy as np
import itertools
from typing import Sequence, Union, Mapping, Optional, MutableMapping, Tuple
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.base.param import _ParamData
from pyomo.contrib.appsi.base import PersistentSolver
from tqdm import tqdm
from .uncertain_params import _replace_uncertain_params
from .inner_problem import _build_inner_problem
import enum
from idaes.surrogate.pysmo.sampling import LatinHypercubeSampling
from .indices import _VarIndex
from pyomo.common.config import ConfigDict, ConfigValue, InEnum


class SamplingStrategy(enum.Enum):
    grid = enum.auto()
    lhs = enum.auto()


def _grid_sampling(uncertain_params: Sequence[_GeneralVarData], num_points: int, seed: int):
    uncertain_params_values = pe.ComponentMap()
    for p in uncertain_params:
        uncertain_params_values[p] = list(set([float(i) for i in np.linspace(p.lb, p.ub, num_points)]))
        uncertain_params_values[p].sort()

    sample_points = pe.ComponentMap()
    for p in uncertain_params:
        sample_points[p] = list()

    n_samples = 0
    for sample in itertools.product(*uncertain_params_values.values()):
        for p, v in zip(uncertain_params, sample):
            sample_points[p].append(float(v))
        n_samples += 1

    return n_samples, sample_points


def _lhs_sampling(uncertain_params: Sequence[_GeneralVarData], num_points: int, seed: int):
    lb_list = list()
    ub_list = list()
    for p in uncertain_params:
        lb_list.append(p.lb)
        ub_list.append(p.ub)
    sampler = LatinHypercubeSampling([lb_list, ub_list], number_of_samples=num_points, sampling_type='creation')

    np.random.seed(seed)
    sample_array = sampler.sample_points()

    sample_points = pe.ComponentMap()
    for ndx, p in enumerate(uncertain_params):
        sample_points[p] = [float(i) for i in sample_array[:, ndx]]

    return num_points, sample_points


_sample_strategy_map = dict()
_sample_strategy_map[SamplingStrategy.grid] = _grid_sampling
_sample_strategy_map[SamplingStrategy.lhs] = _lhs_sampling


class SamplingConfig(ConfigDict):
    def __init__(self, description=None, doc=None, implicit=False, implicit_domain=None, visibility=0):
        super().__init__(description=description, doc=doc, implicit=implicit, implicit_domain=implicit_domain,
                         visibility=visibility)
        self.strategy: SamplingStrategy = self.declare('strategy', ConfigValue(domain=InEnum(SamplingStrategy),
                                                                               default=SamplingStrategy.lhs))
        self.lhs_seed: int = self.declare('lhs_seed', ConfigValue(domain=int, default=0))
        self.solver = self.declare('solver', ConfigValue(default=pe.SolverFactory('appsi_ipopt')))
        self.num_points: int = self.declare('num_points', ConfigValue(domain=int, default=100))
        self.enable_progress_bar: bool = self.declare('enable_progress_bar', ConfigValue(domain=bool, default=True))


def _perform_sampling(m: _BlockData,
                      uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
                      controls: Optional[Sequence[_GeneralVarData]],
                      config: SamplingConfig) -> Tuple[MutableMapping[Union[_GeneralVarData, _ParamData],
                                                                      Sequence[float]],
                                                       Sequence[float],
                                                       MutableMapping[_GeneralVarData, Sequence[float]]]:
    if isinstance(config.solver, PersistentSolver):
        using_persistent = True
    else:
        using_persistent = False

    unc_param_vars = list()
    for p in uncertain_params:
        ndx = _VarIndex(p, None)
        p_var = m.unc_param_vars[ndx]
        unc_param_vars.append(p_var)
    n_samples, sample_points = _sample_strategy_map[config.strategy](unc_param_vars, config.num_points, config.lhs_seed)

    if using_persistent:
        config.solver.set_instance(m)
        original_update_config = config.solver.update_config()
        config.solver.update_config.check_for_new_or_removed_constraints = False
        config.solver.update_config.check_for_new_or_removed_vars = False
        config.solver.update_config.check_for_new_or_removed_params = False
        config.solver.update_config.update_constraints = False
        config.solver.update_config.update_vars = False
        config.solver.update_config.update_params = False
        config.solver.update_config.update_named_expressions = False

    max_violation_values = list()

    control_values = pe.ComponentMap()
    for v in controls:
        control_values[v] = list()

    for sample_ndx in tqdm(list(range(n_samples)), ncols=100, desc='Sampling', disable=not config.enable_progress_bar):
        for p, p_vals in sample_points.items():
            p.fix(p_vals[sample_ndx])

        if using_persistent:
            config.solver.update_variables(unc_param_vars)

        res = config.solver.solve(m)
        pe.assert_optimal_termination(res)
        max_violation_values.append(m.max_constraint_violation.value)

        for v in controls:
            control_values[v].append(v.value)

    if using_persistent:
        config.solver.update_config.set_value(original_update_config)

    unc_param_var_to_unc_param_map = pe.ComponentMap(zip(unc_param_vars, uncertain_params))
    sample_points = pe.ComponentMap((unc_param_var_to_unc_param_map[p], vals) for p, vals in sample_points.items())

    return sample_points, max_violation_values, control_values


def perform_sampling(m: _BlockData,
                     uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
                     param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
                     param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]],
                     controls: Optional[Sequence[_GeneralVarData]],
                     in_place: bool,
                     config: SamplingConfig) -> Tuple[MutableMapping[Union[_GeneralVarData, _ParamData],
                                                                     Sequence[float]],
                                                      Sequence[float],
                                                      MutableMapping[_GeneralVarData, Sequence[float]]]:
    original_model = m
    if not in_place:
        m = m.clone()
        uncertain_params = [m.find_component(p) for p in uncertain_params]
        param_nominal_values = pe.ComponentMap((m.find_component(p), val) for p, val in param_nominal_values.items())
        param_bounds = pe.ComponentMap((m.find_component(p), bnds) for p, bnds in param_bounds.items())
        controls = [m.find_component(v) for v in controls]

    _replace_uncertain_params(m, uncertain_params, param_nominal_values, param_bounds)
    _build_inner_problem(m=m, enforce_equalities=True, unique_constraint_violations=False, valid_var_bounds=None)
    sample_points, max_violation_values, control_values = _perform_sampling(m=m, uncertain_params=uncertain_params,
                                                                            controls=controls, config=config)

    sample_points = pe.ComponentMap((original_model.find_component(p), vals) for p, vals in sample_points.items())
    control_values = pe.ComponentMap((original_model.find_component(v), vals) for v, vals in control_values.items())

    return sample_points, max_violation_values, control_values
