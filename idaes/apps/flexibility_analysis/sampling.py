from __future__ import annotations
from pyomo.core.base.block import _BlockData
import pyomo.environ as pe
import numpy as np
from typing import Sequence, Union, Mapping, Optional, MutableMapping, Tuple, List
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.base.param import _ParamData
from .uncertain_params import _replace_uncertain_params
from .inner_problem import _build_inner_problem
import enum
from idaes.core.surrogate.pysmo.sampling import LatinHypercubeSampling
from .indices import _VarIndex
from pyomo.common.config import ConfigDict, ConfigValue, InEnum
from pyomo.contrib.solver.util import get_objective
from pyomo.common.errors import ApplicationError
from .check_optimal import assert_optimal_termination

try:
    from tqdm import tqdm
except ImportError:

    def tqdm(items, ncols, desc, disable):
        return items


class SamplingStrategy(enum.Enum):
    grid = "grid"
    lhs = "lhs"


SamplingStrategy.grid.__doc__ = r"Use evenly spaced samples"
SamplingStrategy.lhs.__doc__ = r"Use latin hypercube sampling"


class SamplingInitStrategy(enum.Enum):
    none = "none"
    square = "square"
    min_control_deviation = "min_control_deviation"
    all = "all"


SamplingInitStrategy.none.__doc__ = (
    r"Use the solution from the previous sample to initialize the inner problem"
)
SamplingInitStrategy.square.__doc__ = (
    r"Fix the controls and solve a square problem to initialize the inner problem"
)
SamplingInitStrategy.min_control_deviation.__doc__ = r"Fix the maximum constraint violation to 0 and minimized the square of the differences between the controls and their current values"
SamplingInitStrategy.all.__doc__ = r"Try both square and min_control_deviation"


class _GridSamplingState(enum.Enum):
    increment = "increment"
    decrement = "decrement"


class _ParamIterator(object):
    def __init__(
        self,
        param: _GeneralVarData,
        num_points: int,
        next_param: Optional[_ParamIterator],
    ):
        self.state = _GridSamplingState.increment
        self.ndx = 0
        self.pts = list(
            set([float(i) for i in np.linspace(param.lb, param.ub, num_points)])
        )
        self.pts.sort()
        self.next_param = next_param

    def reset(self):
        self.state = _GridSamplingState.increment
        self.ndx = 0

    def get_value(self):
        res = self.pts[self.ndx]
        return res

    def swap_state(self):
        if self.state == _GridSamplingState.increment:
            self.state = _GridSamplingState.decrement
        else:
            assert self.state == _GridSamplingState.decrement
            self.state = _GridSamplingState.increment

    def step(self) -> bool:
        if self.state == _GridSamplingState.increment:
            if self.ndx == len(self.pts) - 1:
                if self.next_param is None:
                    done = True
                else:
                    done = self.next_param.step()
                    self.swap_state()
            else:
                self.ndx += 1
                done = False
        else:
            assert self.state == _GridSamplingState.decrement
            if self.ndx == 0:
                assert self.next_param is not None
                done = self.next_param.step()
                self.swap_state()
            else:
                self.ndx -= 1
                done = False
        return done


class _GridSamplingIterator(object):
    def __init__(self, uncertain_params: Sequence[_GeneralVarData], num_points: int):
        self.params = list(uncertain_params)
        self.param_iterators: List[Optional[_ParamIterator]] = [None] * len(self.params)
        self.param_iterators[-1] = _ParamIterator(
            param=self.params[-1], num_points=num_points, next_param=None
        )
        for ndx in reversed(range(len(self.params) - 1)):
            self.param_iterators[ndx] = _ParamIterator(
                param=self.params[ndx],
                num_points=num_points,
                next_param=self.param_iterators[ndx + 1],
            )
        self.done = False

    def __next__(self):
        if self.done:
            raise StopIteration

        res = [i.get_value() for i in self.param_iterators]
        self.done = self.param_iterators[0].step()

        return res

    def __iter__(self):
        [i.reset() for i in self.param_iterators]
        self.done = False
        return self


def _grid_sampling(
    uncertain_params: Sequence[_GeneralVarData], num_points: int, seed: int
):
    it = _GridSamplingIterator(uncertain_params=uncertain_params, num_points=num_points)

    sample_points = pe.ComponentMap()
    for p in uncertain_params:
        sample_points[p] = list()

    n_samples = 0
    for sample in it:
        for p, v in zip(uncertain_params, sample):
            sample_points[p].append(float(v))
        n_samples += 1

    return n_samples, sample_points


def _lhs_sampling(
    uncertain_params: Sequence[_GeneralVarData], num_points: int, seed: int
):
    lb_list = list()
    ub_list = list()
    for p in uncertain_params:
        lb_list.append(p.lb)
        ub_list.append(p.ub)
    sampler = LatinHypercubeSampling(
        [lb_list, ub_list], number_of_samples=num_points, sampling_type="creation"
    )

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
    r"""
    A class for specifying options for sampling the uncertain parameter values
    and solving the inner problem of the flexibility test.

    Attributes
    ----------
    strategy: SamplingStrategy
        The method for sampling the uncertain parameters. (default: SamplingStrategy.lhs)
    lhs_seed: int
        The seed used for latin hypercube sampling (default: 0)
    solver: Union[Solver, OptSolver]
        The solver to use for the inner problem of the flexibility test problem
    num_points: int
        The number of samples of uncertain parameter values to use (default: 100)
    enable_progress_bar: bool
        If False, no progress bar will be shown (default: True)
    initialization_strategy: SamplingInitStrategy
        The initialization strategy to use for the inner problems of the
        flexibility test at each sample of the uncertain parameter values.
        (default: SamplingInitStrategy.none)
    total_violation: bool
        If True, the objective of the flexibility test will be the sum of
        the constraint violations instead of the maximum violation. (default: False)
    """

    def __init__(
        self,
        description=None,
        doc=None,
        implicit=False,
        implicit_domain=None,
        visibility=0,
    ):
        super().__init__(
            description=description,
            doc=doc,
            implicit=implicit,
            implicit_domain=implicit_domain,
            visibility=visibility,
        )

        self.strategy: SamplingStrategy = self.declare(
            "strategy",
            ConfigValue(domain=InEnum(SamplingStrategy), default=SamplingStrategy.lhs),
        )
        self.lhs_seed: int = self.declare(
            "lhs_seed", ConfigValue(domain=int, default=0)
        )
        self.solver = self.declare(
            "solver", ConfigValue(default=pe.SolverFactory("appsi_ipopt"))
        )
        self.num_points: int = self.declare(
            "num_points", ConfigValue(domain=int, default=100)
        )
        self.enable_progress_bar: bool = self.declare(
            "enable_progress_bar", ConfigValue(domain=bool, default=True)
        )
        self.initialization_strategy: SamplingInitStrategy = self.declare(
            "initialization_strategy",
            ConfigValue(
                domain=InEnum(SamplingInitStrategy), default=SamplingInitStrategy.none
            ),
        )
        self.total_violation: bool = self.declare(
            "total_violation", ConfigValue(domain=bool, default=False)
        )


def _deactivate_inequalities(m: _BlockData):
    deactivated_cons = list()
    for c in m.component_data_objects(pe.Constraint, descend_into=True, active=True):
        if not c.equality and c.lb != c.ub:
            deactivated_cons.append(c)
            c.deactivate()
    return deactivated_cons


def _init_with_square_problem(m: _BlockData, controls, solver):
    for v in controls:
        if v.value is None:
            raise RuntimeError(
                "Cannot initialize sampling problem with square problem because the "
                "control values are not initialized."
            )
        v.fix()
    if m.max_constraint_violation.value is None:
        m.max_constraint_violation.value = 0
    m.max_constraint_violation.fix()
    deactivated_cons = _deactivate_inequalities(m)
    try:
        res = solver.solve(m, tee=False, load_solutions=False)
    except ApplicationError:
        res = None
    for c in deactivated_cons:
        c.activate()
    for v in controls:
        v.unfix()
    if res is not None and pe.check_optimal_termination(res):
        m.solutions.load_from(res)
        max_viol = 0
        for c in deactivated_cons:
            assert c.lb is None
            assert c.ub == 0
            body_val = pe.value(c.body)
            if body_val > max_viol:
                max_viol = body_val
        m.max_constraint_violation.value = max_viol
    else:
        max_viol = None
    m.max_constraint_violation.unfix()
    return max_viol


def _solve_with_max_viol_fixed(m: _BlockData, controls, solver):
    orig_obj = get_objective(m)
    orig_obj.deactivate()
    orig_max_viol_value = m.max_constraint_violation.value
    m.max_constraint_violation.fix(0)

    obj_expr = 0
    for v in controls:
        if v.value is None:
            obj_expr += v**2
        else:
            obj_expr += (v - v.value) ** 2
    m.control_setpoints_obj = pe.Objective(expr=obj_expr)

    try:
        res = solver.solve(m, tee=False, load_solutions=False)
        if pe.check_optimal_termination(res):
            m.solutions.load_from(res)
            feasible = True
            control_vals = [v.value for v in controls]
        else:
            feasible = False
            control_vals = None
            m.max_constraint_violation.value = orig_max_viol_value
    except ApplicationError:
        feasible = False
        control_vals = None
        m.max_constraint_violation.value = orig_max_viol_value

    del m.control_setpoints_obj
    m.max_constraint_violation.unfix()
    orig_obj.activate()

    return feasible, control_vals


def _perform_sampling(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    controls: Optional[Sequence[_GeneralVarData]],
    config: SamplingConfig,
) -> Tuple[
    MutableMapping[Union[_GeneralVarData, _ParamData], Sequence[float]],
    Sequence[float],
    MutableMapping[_GeneralVarData, Sequence[float]],
]:
    unc_param_vars = list()
    for p in uncertain_params:
        ndx = _VarIndex(p, None)
        p_var = m.unc_param_vars[ndx]
        unc_param_vars.append(p_var)
    n_samples, sample_points = _sample_strategy_map[config.strategy](
        unc_param_vars, config.num_points, config.lhs_seed
    )

    obj_values = list()
    obj = get_objective(m)

    control_values = pe.ComponentMap()
    for v in controls:
        control_values[v] = list()

    if not config.total_violation:
        m.max_constraint_violation.value = 0
        orig_max_constraint_violation_ub = m.max_constraint_violation.ub

    for sample_ndx in tqdm(
        list(range(n_samples)),
        ncols=100,
        desc="Sampling",
        disable=not config.enable_progress_bar,
    ):
        for p, p_vals in sample_points.items():
            p.fix(p_vals[sample_ndx])

        if not config.total_violation and config.initialization_strategy in {
            SamplingInitStrategy.square,
            SamplingInitStrategy.all,
        }:
            max_viol_ub = _init_with_square_problem(m, controls, config.solver)
            if max_viol_ub is not None:
                m.max_constraint_violation.setub(max_viol_ub)
        if not config.total_violation and config.initialization_strategy in {
            SamplingInitStrategy.min_control_deviation,
            SamplingInitStrategy.all,
        }:
            feasible, control_vals = _solve_with_max_viol_fixed(
                m, controls, config.solver
            )
        else:
            feasible = False
            control_vals = None
        if feasible:
            obj_values.append(0)
            for v, val in zip(controls, control_vals):
                control_values[v].append(val)
        else:
            res = config.solver.solve(m)
            assert_optimal_termination(res)
            obj_values.append(pe.value(obj.expr))

            for v in controls:
                control_values[v].append(v.value)

        if not config.total_violation:
            m.max_constraint_violation.setub(orig_max_constraint_violation_ub)

    unc_param_var_to_unc_param_map = pe.ComponentMap(
        zip(unc_param_vars, uncertain_params)
    )
    sample_points = pe.ComponentMap(
        (unc_param_var_to_unc_param_map[p], vals) for p, vals in sample_points.items()
    )

    return sample_points, obj_values, control_values


def perform_sampling(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]],
    controls: Optional[Sequence[_GeneralVarData]],
    in_place: bool,
    config: SamplingConfig,
) -> Tuple[
    MutableMapping[Union[_GeneralVarData, _ParamData], Sequence[float]],
    Sequence[float],
    MutableMapping[_GeneralVarData, Sequence[float]],
]:
    original_model = m
    if not in_place:
        m = m.clone()
        uncertain_params = [m.find_component(p) for p in uncertain_params]
        param_nominal_values = pe.ComponentMap(
            (m.find_component(p), val) for p, val in param_nominal_values.items()
        )
        param_bounds = pe.ComponentMap(
            (m.find_component(p), bnds) for p, bnds in param_bounds.items()
        )
        controls = [m.find_component(v) for v in controls]

    _replace_uncertain_params(m, uncertain_params, param_nominal_values, param_bounds)
    _build_inner_problem(
        m=m,
        enforce_equalities=True,
        unique_constraint_violations=False,
        valid_var_bounds=None,
        total_violation=config.total_violation,
        total_violation_disjunctions=False,
    )
    sample_points, obj_values, control_values = _perform_sampling(
        m=m, uncertain_params=uncertain_params, controls=controls, config=config
    )

    sample_points = pe.ComponentMap(
        (original_model.find_component(p), vals) for p, vals in sample_points.items()
    )
    control_values = pe.ComponentMap(
        (original_model.find_component(v), vals) for v, vals in control_values.items()
    )

    return sample_points, obj_values, control_values
