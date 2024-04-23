import numpy as np
from .kkt import add_kkt_with_milp_complementarity_conditions
from pyomo.core.base.block import _BlockData
from pyomo.common.dependencies import attempt_import

coramin, coramin_available = attempt_import(
    "coramin", "coramin is required for flexibility analysis"
)
import pyomo.environ as pe
from .var_utils import (
    get_used_unfixed_variables,
    BoundsManager,
    _remove_var_bounds,
    _apply_var_bounds,
)
from .indices import _VarIndex, _ConIndex
from .uncertain_params import _replace_uncertain_params
from .inner_problem import _build_inner_problem
from pyomo.util.report_scaling import report_scaling
import logging
from typing import Sequence, Union, Mapping, MutableMapping, Optional, Tuple
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.base.param import _ParamData
from .decision_rules.linear_dr import construct_linear_decision_rule
from .decision_rules.dr_config import DRConfig
from pyomo.common.dependencies import attempt_import
from .sampling import (
    SamplingStrategy,
    perform_sampling,
    SamplingConfig,
    _perform_sampling,
)
import enum
from pyomo.common.config import (
    ConfigDict,
    ConfigValue,
    PositiveFloat,
    InEnum,
    MarkImmutable,
    NonNegativeFloat,
)
from pyomo.contrib.appsi.base import Solver
from pyomo.opt.base import OptSolver

relu_dr, relu_dr_available = attempt_import(
    "idaes.apps.flexibility_analysis.decision_rules.relu_dr",
    "The ReLU decision rule requires Tensorflow and OMLT",
)


logger = logging.getLogger(__name__)


def _get_longest_name(comps):
    longest_name = 0

    for i in comps:
        i_len = len(str(i))
        if i_len > longest_name:
            longest_name = i_len

    if longest_name > 195:
        longest_name = 195
    if longest_name < 12:
        longest_name = 12

    return longest_name


class FlexTestMethod(enum.Enum):
    active_constraint = enum.auto()
    linear_decision_rule = enum.auto()
    relu_decision_rule = enum.auto()
    vertex_enumeration = enum.auto()
    sampling = enum.auto()


FlexTestMethod.active_constraint.__doc__ = r"Solve the flexibility test using the active constraint method described in [Grossmann1987]_."
FlexTestMethod.linear_decision_rule.__doc__ = r"Solve the flexibility test by converting the inner minimization problem to a square problem by removing all degrees of freedom by creating a linear decision rule of the form :math:`z = A \theta + b`"
FlexTestMethod.relu_decision_rule.__doc__ = r"Solve the flexibility test by converting the inner minimization problem to a square problem by removing all degrees of freedom by creating a decision rule of the form :math:`z = f(\theta)` where :math:`f(\theta)` is a nueral network with ReLU activation functions."
FlexTestMethod.vertex_enumeration.__doc__ = r"Solve the flexibility test by solving the inner minimization problem at every vertex of the hyperrectangle defined by :math:`(\underline{\theta}, \overline{\theta})`."
FlexTestMethod.sampling.__doc__ = r"Solve the flexibility test by solving the inner minimization problem at random samples of :math:`\theta \in [\underline{\theta}, \overline{\theta}]`."


class ActiveConstraintConfig(ConfigDict):
    r"""
    A class for specifying options for the active constraint method for the
    flexibility test problem.

    Attributes
    ----------
    use_haar_conditions: bool
        If False, no constraint will be added to constraint the number of
        active inequalities. (default: True)
    default_BigM: float
        Default value for the bigM parameter used to reformulate the
        complimentarity conditions in the KKT system. (default: None)
    enforce_equalities: bool
        If False, :math:`h(x, z, \theta) = 0` is treated as two inequalities
        (performance constraints) that can be violated (:math:`h_{i}(x, z, \theta) \leq u`
        and :math:`-h_{i}(x, z, \theta) \leq u`) (default: True)
    skip_scaling_check: bool
        If True, the model scaling will not be checked. (default: False)
    total_violation: bool
        If True, the objective of the flexibility test will be the sum of
        the constraint violations instead of the maximum violation. (default: False)
    """

    def __init__(self):
        super().__init__(
            description=None,
            doc=None,
            implicit=False,
            implicit_domain=None,
            visibility=0,
        )
        self.use_haar_conditions: bool = self.declare(
            "use_haar_conditions", ConfigValue(domain=bool, default=True)
        )
        self.default_BigM: Optional[float] = self.declare(
            "default_BigM", ConfigValue(domain=NonNegativeFloat, default=None)
        )
        self.enforce_equalities: bool = self.declare(
            "enforce_equalities", ConfigValue(domain=bool, default=True)
        )
        self.skip_scaling_check: bool = self.declare(
            "skip_scaling_check", ConfigValue(domain=bool, default=False)
        )
        self.total_violation: bool = self.declare(
            "total_violation", ConfigValue(domain=bool, default=False)
        )


class FlexTestConfig(ConfigDict):
    r"""
    A class for specifying options for solving the flexibility test.

    Attributes
    ----------
    feasibility_tol: float
        Tolerance for considering constraints to be satisfied. In particular, if the
        maximum constraint violation is less than or equal to :py:attr:`feasibility_tol<FlexTestConfig.feasibility_tol>`, then
        the flexibility test passes. (default: 1e-6)
    terminate_early: bool
        If True, the specified algorithm should terminate as soon as a point
        (:math:`\theta`) is found that confirms the flexibility test fails. If
        False, the specified algorithm will continue until the :math:`\theta`
        that maximizes the constraint violation is found. (default: False)
    method: FlexTestMethod
        The method that should be used to solve the flexibility test. (default: :py:attr:`active_constraint<FlexTestMethod.active_constraint>`)
    minlp_solver: Union[Solver, OptSolver]
        A Pyomo solver interface appropriate for solving MINLPs
    sampling_config: SamplingConfig
        A config object for specifying how sampling should be performed when either
        generating data to create a decision rule or using sampling to solve the
        flexibility test.
    decision_rule_config: DRConfig
        Only used if method is one of the decision rules. Should be either a LinearDRConfig
        or a ReluDRConfig.
    active_constraint_config: ActiveConstraintConfig
        Only used if :py:attr:`method<FlexTestConfig.method>` is :py:attr:`active_constraint<FlexTestMethod.active_constraint>`
    total_violation: bool
        If False, the maximum constraint violation is considered. If True, the sum
        of the violations of all constraints is considered. Should normally be False. (default: False)
    """

    def __init__(self):
        super().__init__(
            description=None,
            doc=None,
            implicit=False,
            implicit_domain=None,
            visibility=0,
        )
        self.feasibility_tol: float = self.declare(
            "feasibility_tol", ConfigValue(domain=PositiveFloat, default=1e-6)
        )
        self.terminate_early: bool = self.declare(
            "terminate_early", ConfigValue(domain=bool, default=False)
        )
        self.method: FlexTestMethod = self.declare(
            "method",
            ConfigValue(
                domain=InEnum(FlexTestMethod), default=FlexTestMethod.active_constraint
            ),
        )
        self.minlp_solver = self.declare("minlp_solver", ConfigValue())
        self.sampling_config: SamplingConfig = self.declare(
            "sampling_config", SamplingConfig()
        )
        self.decision_rule_config = self.declare(
            "decision_rule_config", ConfigValue(default=None)
        )
        self.active_constraint_config: ActiveConstraintConfig = self.declare(
            "active_constraint_config", ActiveConstraintConfig()
        )
        self.total_violation: bool = self.declare(
            "total_violation", ConfigValue(domain=bool, default=False)
        )


class FlexTestTermination(enum.Enum):
    found_infeasible_point = enum.auto()
    proven_feasible = enum.auto()
    uncertain = enum.auto()


FlexTestTermination.found_infeasible_point.__doc__ = r"The meaning of this member depends on the method used to solve the flexibility/feasibility test, but it generally means that the flexibility test failed. If the solution method is not conservative (:py:attr:`FlexTestMethod.vertex_enumeration<idaes.apps.flexibility_analysis.FlexTestMethod.vertex_enumeration>`, :py:attr:`FlexTestMethod.sampling<idaes.apps.flexibility_analysis.FlexTestMethod.sampling>`), then :py:attr:`FlexTestTermination.found_infeasible_point<idaes.apps.flexibility_analysis.FlexTestTermination.found_infeasible_point>` indicates that a value of :math:`\theta` was found where at least one performance constraint (:math:`g_{j}(x, z, \theta) \leq 0`) is violated. Otherwise, :py:attr:`FlexTestTermination.found_infeasible_point<idaes.apps.flexibility_analysis.FlexTestTermination.found_infeasible_point>` indicates that a point was found where the performance constraints might be violated."
FlexTestTermination.proven_feasible.__doc__ = r"The meaning of this member depends on the method used to solve the flexibility/feasibility test, but it generally means that the flexibility test passed. If the solution method is conservative (:py:attr:`FlexTestMethod.active_constraint<idaes.apps.flexibility_analysis.FlexTestMethod.active_constraint>`, :py:attr:`FlexTestMethod.linear_decision_rule<idaes.apps.flexibility_analysis.FlexTestMethod.linear_decision_rule>`, :py:attr:`FlexTestMethod.relu_decision_rule<idaes.apps.flexibility_analysis.FlexTestMethod.relu_decision_rule>`), then :py:attr:`FlexTestTermination.proven_feasible<idaes.apps.flexibility_analysis.FlexTestTermination.proven_feasible>` indicates that, for any :math:`\theta \in [\underline{\theta}, \overline{\theta}]`, there exists a :math:`z` such that all of the performance constraints (:math:`g_{j}(x, z, \theta) \leq 0`) are satisfied. Otherwise, :py:attr:`FlexTestTermination.proven_feasible<idaes.apps.flexibility_analysis.FlexTestTermination.proven_feasible>` just indicates that no :math:`\theta \in [\underline{\theta}, \overline{\theta}]` was found that violates the performance constraints for all :math:`z`."
FlexTestTermination.uncertain.__doc__ = r"Cannot definitively say whether the flexibility test passes or fails. This usually indicates an error was encountered."


class FlexTestResults(object):
    r"""
    Results for the flexibility test problem.

    Attributes
    ----------
    termination: FlexTestTermination
    max_constraint_violation: float
        The largest constraint violation found (:math:`u`)
    unc_param_values_at_max_violation: Optional[MutableMapping[Union[_GeneralVarData, _ParamData], float]]
        The values of the uncertain parameters that generated the maximum constraint violation
    """

    def __init__(self):
        self.termination = FlexTestTermination.uncertain
        self.max_constraint_violation: Optional[float] = None
        self.unc_param_values_at_max_violation: Optional[
            MutableMapping[Union[_GeneralVarData, _ParamData], float]
        ] = None

    def __str__(self):
        s = f"Termination: {self.termination}\n"
        s += f"Maximum constraint violation: {self.max_constraint_violation}\n"
        if self.unc_param_values_at_max_violation is not None:
            s += f"Uncertain parameter values at maximum constraint violation: \n"
            longest_param_name = _get_longest_name(
                self.unc_param_values_at_max_violation.keys()
            )
            s += f'{"Param":<{longest_param_name + 5}}{"Value":>12}\n'
            for k, v in self.unc_param_values_at_max_violation.items():
                s += f"{str(k):<{longest_param_name + 5}}{v:>12.2e}\n"
        return s


def _get_dof(m: _BlockData):
    n_cons = len(
        set(
            i
            for i in m.component_data_objects(
                pe.Constraint, active=True, descend_into=True
            )
            if i.equality
        )
    )
    n_vars = len(get_used_unfixed_variables(m))
    return n_vars - n_cons


dr_construction_map = dict()
dr_construction_map[FlexTestMethod.linear_decision_rule] = (
    construct_linear_decision_rule
)


def build_flextest_with_dr(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Tuple[float, float]],
    controls: Sequence[_GeneralVarData],
    valid_var_bounds: MutableMapping[_GeneralVarData, Tuple[float, float]],
    config: FlexTestConfig,
):
    config.sampling_config.total_violation = config.total_violation

    # this has to be here in case tensorflow or omlt are not installed
    dr_construction_map[FlexTestMethod.relu_decision_rule] = (
        relu_dr.construct_relu_decision_rule
    )

    # enforce_equalities must be true for this method, or the resulting
    # problem will be unbounded; the key is degrees of freedom

    # perform sampling
    tmp = perform_sampling(
        m=m,
        uncertain_params=uncertain_params,
        param_nominal_values=param_nominal_values,
        param_bounds=param_bounds,
        controls=controls,
        in_place=False,
        config=config.sampling_config,
    )

    _sample_points, max_violation_values, control_values = tmp

    # replace uncertain parameters with variables
    _replace_uncertain_params(m, uncertain_params, param_nominal_values, param_bounds)
    for v in m.unc_param_vars.values():
        valid_var_bounds[v] = (v.lb, v.ub)
        v.fix()  # these should be fixed before we check the degrees of freedom
    for p, p_bnds in param_bounds.items():
        if p.is_variable_type():
            valid_var_bounds[p] = p_bnds

    if _get_dof(m) != len(controls):
        raise ValueError(
            "The number of controls must match the number of degrees of freedom"
        )

    # check the scaling of the model
    # this has to be done with valid_var_bounds (original bounds removed) to ensure we have
    # an entry in valid_var_bounds for every variable
    for v in m.unc_param_vars.values():
        v.unfix()
    bounds_manager = BoundsManager(m)
    bounds_manager.save_bounds()
    _remove_var_bounds(m)
    _apply_var_bounds(valid_var_bounds)
    passed = report_scaling(m)
    if not passed:
        raise ValueError(
            "Please scale the model. If a scaling report was not shown, "
            "set the logging level to INFO."
        )
    bounds_manager.pop_bounds()

    # construct the decision rule
    # the keys of sample_points need to be the new variables
    # that replaced the uncertain parameters
    sample_points: MutableMapping[_GeneralVarData, Sequence[float]] = pe.ComponentMap()
    for p in uncertain_params:
        ndx = _VarIndex(p, None)
        p_var = m.unc_param_vars[ndx]
        sample_points[p_var] = _sample_points[p]

    dr = dr_construction_map[config.method](
        input_vals=sample_points,
        output_vals=control_values,
        config=config.decision_rule_config,
    )

    if config.total_violation:
        total_violation_disjunctions = True
    else:
        total_violation_disjunctions = False
    _build_inner_problem(
        m=m,
        enforce_equalities=True,
        unique_constraint_violations=True,
        valid_var_bounds=valid_var_bounds,
        total_violation=config.total_violation,
        total_violation_disjunctions=total_violation_disjunctions,
    )
    _apply_var_bounds(valid_var_bounds)
    m.decision_rule = dr

    obj = coramin.utils.get_objective(m)
    obj.deactivate()

    if config.total_violation:
        m.max_total_violation_obj = pe.Objective(
            expr=sum(m.constraint_violation.values()), sense=pe.maximize
        )
    else:
        m.max_constraint_violation_obj = pe.Objective(
            expr=m.max_constraint_violation, sense=pe.maximize
        )


def build_active_constraint_flextest(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Tuple[float, float]],
    valid_var_bounds: MutableMapping[_GeneralVarData, Tuple[float, float]],
    config: Optional[ActiveConstraintConfig] = None,
):
    if config is None:
        config = ActiveConstraintConfig()

    _replace_uncertain_params(m, uncertain_params, param_nominal_values, param_bounds)
    for v in m.unc_param_vars.values():
        valid_var_bounds[v] = v.bounds
    for p, p_bnds in param_bounds.items():
        if p.is_variable_type():
            valid_var_bounds[p] = p_bnds

    # TODO: make this a context manager or try-finally
    if not config.skip_scaling_check:
        bounds_manager = BoundsManager(m)
        bounds_manager.save_bounds()
        _remove_var_bounds(m)
        _apply_var_bounds(valid_var_bounds)
        passed = report_scaling(m)
        if not passed:
            raise ValueError(
                "Please scale the model. If a scaling report was not "
                "shown, set the logging level to INFO."
            )
        bounds_manager.pop_bounds()

    # TODO: constraint.equality does not check for range constraints with equal bounds
    orig_equality_cons = [
        c
        for c in m.component_data_objects(pe.Constraint, descend_into=True, active=True)
        if c.equality
    ]

    _build_inner_problem(
        m=m,
        enforce_equalities=config.enforce_equalities,
        unique_constraint_violations=False,
        valid_var_bounds=valid_var_bounds,
        total_violation=config.total_violation,
        total_violation_disjunctions=False,
    )

    for v in m.unc_param_vars.values():
        v.fix()
    n_dof = _get_dof(m)
    for v in m.unc_param_vars.values():
        v.unfix()

    add_kkt_with_milp_complementarity_conditions(
        m=m,
        uncertain_params=list(m.unc_param_vars.values()),
        valid_var_bounds=valid_var_bounds,
        default_M=config.default_BigM,
    )

    # TODO: to control the namespace and reduce cloning:
    #  take the users model and stick it on a new block as a sub-block

    if not config.enforce_equalities and not config.total_violation:
        m.equality_cuts = pe.ConstraintList()
        max_viol_lb, max_viol_ub = valid_var_bounds[m.max_constraint_violation]
        for c in orig_equality_cons:
            key1 = _ConIndex(c, "lb")
            key2 = _ConIndex(m.ineq_violation_cons[key1], "ub")
            y1 = m.active_indicator[key2]
            key1 = _ConIndex(c, "ub")
            key2 = _ConIndex(m.ineq_violation_cons[key1], "ub")
            y2 = m.active_indicator[key2]
            m.equality_cuts.add(
                m.max_constraint_violation <= (1 - y1 * y2) * max_viol_ub
            )
            m.equality_cuts.add(
                m.max_constraint_violation >= (1 - y1 * y2) * max_viol_lb
            )

    if config.use_haar_conditions and not config.total_violation:
        m.n_active_ineqs = pe.Constraint(expr=sum(m.active_indicator.values()) == n_dof)

    if config.total_violation:
        m.max_total_violation_obj = pe.Objective(
            expr=sum(m.constraint_violation.values()), sense=pe.maximize
        )
    else:
        m.max_constraint_violation_obj = pe.Objective(
            expr=m.max_constraint_violation, sense=pe.maximize
        )


def _solve_flextest_active_constraint(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Tuple[float, float]],
    controls: Sequence[_GeneralVarData],
    valid_var_bounds: MutableMapping[_GeneralVarData, Tuple[float, float]],
    config: Optional[FlexTestConfig] = None,
) -> FlexTestResults:
    build_active_constraint_flextest(
        m=m,
        uncertain_params=uncertain_params,
        param_nominal_values=param_nominal_values,
        param_bounds=param_bounds,
        valid_var_bounds=valid_var_bounds,
        config=config.active_constraint_config,
    )
    opt = config.minlp_solver
    res = opt.solve(m)
    pe.assert_optimal_termination(res)

    results = FlexTestResults()
    results.max_constraint_violation = m.max_constraint_violation.value
    if results.max_constraint_violation > config.feasibility_tol:
        results.termination = FlexTestTermination.found_infeasible_point
    else:
        results.termination = FlexTestTermination.proven_feasible
    results.unc_param_values_at_max_violation = pe.ComponentMap()
    for key, v in m.unc_param_vars.items():
        results.unc_param_values_at_max_violation[key.var] = v.value
    return results


def _solve_flextest_decision_rule(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Tuple[float, float]],
    controls: Sequence[_GeneralVarData],
    valid_var_bounds: MutableMapping[_GeneralVarData, Tuple[float, float]],
    config: Optional[FlexTestConfig] = None,
) -> FlexTestResults:
    build_flextest_with_dr(
        m=m,
        uncertain_params=uncertain_params,
        param_nominal_values=param_nominal_values,
        param_bounds=param_bounds,
        controls=controls,
        valid_var_bounds=valid_var_bounds,
        config=config,
    )
    opt = config.minlp_solver
    res = opt.solve(m, tee=True)
    pe.assert_optimal_termination(res)

    results = FlexTestResults()
    results.max_constraint_violation = m.max_constraint_violation.value
    if results.max_constraint_violation > config.feasibility_tol:
        results.termination = FlexTestTermination.found_infeasible_point
    else:
        results.termination = FlexTestTermination.proven_feasible
    results.unc_param_values_at_max_violation = pe.ComponentMap()
    for key, v in m.unc_param_vars.items():
        results.unc_param_values_at_max_violation[key.var] = v.value
    return results


def _solve_flextest_sampling(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Tuple[float, float]],
    controls: Sequence[_GeneralVarData],
    valid_var_bounds: MutableMapping[_GeneralVarData, Tuple[float, float]],
    config: Optional[FlexTestConfig] = None,
) -> FlexTestResults:
    config.sampling_config.total_violation = config.total_violation
    tmp = perform_sampling(
        m=m,
        uncertain_params=uncertain_params,
        param_nominal_values=param_nominal_values,
        param_bounds=param_bounds,
        controls=controls,
        in_place=True,
        config=config.sampling_config,
    )
    sample_points, max_violation_values, control_values = tmp
    max_viol_ndx = int(np.argmax(max_violation_values))

    results = FlexTestResults()
    results.max_constraint_violation = max_violation_values[max_viol_ndx]
    if results.max_constraint_violation > config.feasibility_tol:
        results.termination = FlexTestTermination.found_infeasible_point
    else:
        results.termination = FlexTestTermination.proven_feasible
    results.unc_param_values_at_max_violation = pe.ComponentMap()
    for key, vals in sample_points.items():
        results.unc_param_values_at_max_violation[key] = vals[max_viol_ndx]
    return results


def _solve_flextest_vertex_enumeration(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Tuple[float, float]],
    controls: Sequence[_GeneralVarData],
    valid_var_bounds: MutableMapping[_GeneralVarData, Tuple[float, float]],
    config: Optional[FlexTestConfig] = None,
) -> FlexTestResults:
    config: FlexTestConfig = config()
    config.sampling_config.num_points = 2
    config.sampling_config.strategy = SamplingStrategy.grid
    config.sampling_config.total_violation = config.total_violation
    tmp = perform_sampling(
        m=m,
        uncertain_params=uncertain_params,
        param_nominal_values=param_nominal_values,
        param_bounds=param_bounds,
        controls=controls,
        in_place=True,
        config=config.sampling_config,
    )
    sample_points, max_violation_values, control_values = tmp
    max_viol_ndx = int(np.argmax(max_violation_values))

    results = FlexTestResults()
    results.max_constraint_violation = max_violation_values[max_viol_ndx]
    if results.max_constraint_violation > config.feasibility_tol:
        results.termination = FlexTestTermination.found_infeasible_point
    else:
        results.termination = FlexTestTermination.proven_feasible
    results.unc_param_values_at_max_violation = pe.ComponentMap()
    for key, vals in sample_points.items():
        results.unc_param_values_at_max_violation[key] = vals[max_viol_ndx]
    return results


_flextest_map = dict()
_flextest_map[FlexTestMethod.active_constraint] = _solve_flextest_active_constraint
_flextest_map[FlexTestMethod.sampling] = _solve_flextest_sampling
_flextest_map[FlexTestMethod.vertex_enumeration] = _solve_flextest_vertex_enumeration
_flextest_map[FlexTestMethod.linear_decision_rule] = _solve_flextest_decision_rule
_flextest_map[FlexTestMethod.relu_decision_rule] = _solve_flextest_decision_rule


def solve_flextest(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Tuple[float, float]],
    controls: Sequence[_GeneralVarData],
    valid_var_bounds: MutableMapping[_GeneralVarData, Tuple[float, float]],
    in_place: bool = False,
    config: Optional[FlexTestConfig] = None,
) -> FlexTestResults:
    r"""
    Parameters
    ----------
    m: _BlockData
        The pyomo model to be used for the feasibility/flexibility test.
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]]
        A sequence (e.g., list) defining the set of uncertain parameters (:math:`\theta`).
        These can be pyomo variables (Var) or parameters (param). However, if parameters are used,
        they must be mutable.
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float]
        A mapping (e.g., ComponentMap) from the uncertain parameters (:math:`\theta`) to their
        nominal values (:math:`\theta^{N}`).
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Tuple[float, float]]
        A mapping (e.g., ComponentMap) from the uncertain parameters (:math:`\theta`) to their
        bounds (:math:`\underline{\theta}`, :math:`\overline{\theta}`).
    controls: Sequence[_GeneralVarData]
        A sequence (e.g., list) defining the set of control variables (:math:`z`).
    valid_var_bounds: MutableMapping[_GeneralVarData, Tuple[float, float]]
        A mapping (e.g., ComponentMap) defining bounds for all variables (:math:`x` and :math:`z`) that
        should be valid for any :math:`\theta` between :math:`\underline{\theta}` and
        :math:`\overline{\theta}`. These are only used to make the resulting flexibility test problem
        more computationally tractable. All variable bounds in the model `m` are treated as performance
        constraints and relaxed (:math:`g_{j}(x, z, \theta) \leq u`). The bounds in `valid_var_bounds`
        are applied to the single-level problem generated from the active constraint method or one of
        the decision rules. This argument is not necessary for vertex enumeration or sampling.
    in_place: bool
        If True, m is modified in place to generate the model for solving the flexibility test. If False,
        the model is cloned first.
    config: Optional[FlexTestConfig]
        An object defining options for how the flexibility test should be solved.
    """
    if config is None:
        config = FlexTestConfig()

    original_model = m
    original_uncertain_params = uncertain_params
    original_param_nominal_values = param_nominal_values
    original_param_bounds = param_bounds
    original_controls = controls
    original_valid_var_bounds = valid_var_bounds
    if not in_place:
        # TODO:
        #  tmp_name = pyomo.common.modeling.unique_component_name(m, 'tmp_data')
        #  setattr(m, tmp_name, uncertain_params)
        #  new_m = m.clone()
        #  old_to_new_params = ComponentMap(zip(getattr(m, tmp_name),
        #                                       getattr(new_m, tmp_name)))
        #  delattr(m, tmp_name)
        #  m = new_m
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
    results = _flextest_map[config.method](
        m=m,
        uncertain_params=uncertain_params,
        param_nominal_values=param_nominal_values,
        param_bounds=param_bounds,
        controls=controls,
        valid_var_bounds=valid_var_bounds,
        config=config,
    )
    if not in_place:
        unc_param_values = pe.ComponentMap()
        for v, val in results.unc_param_values_at_max_violation.items():
            unc_param_values[original_model.find_component(v)] = val
        results.unc_param_values_at_max_violation = unc_param_values
    return results


class FlexTest(object):
    def __init__(
        self,
        m: _BlockData,
        uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
        param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
        max_param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]],
        controls: Sequence[_GeneralVarData],
        valid_var_bounds: MutableMapping[_GeneralVarData, Sequence[float]],
        config: Optional[FlexTestConfig] = None,
    ):
        if config is None:
            self.config: FlexTestConfig = FlexTestConfig()
        else:
            self.config: FlexTestConfig = config()
        MarkImmutable(self.config.get("method"))
        if self.config.method == FlexTestMethod.vertex_enumeration:
            self.config.sampling_config.strategy = SamplingStrategy.grid
            self.config.sampling_config.num_points = 2
            MarkImmutable(self.config.sampling_config.get("strategy"))
            MarkImmutable(self.config.sampling_config.get("num_points"))

        self._original_model = m
        self._model = m.clone()
        m = self._model
        self._uncertain_params = [m.find_component(i) for i in uncertain_params]
        self._param_nominal_values = pe.ComponentMap(
            (p, param_nominal_values[orig_p])
            for orig_p, p in zip(uncertain_params, self._uncertain_params)
        )
        self._max_param_bounds = pe.ComponentMap(
            (p, max_param_bounds[orig_p])
            for orig_p, p in zip(uncertain_params, self._uncertain_params)
        )
        self._controls = [m.find_component(i) for i in controls]
        self._valid_var_bounds = pe.ComponentMap(
            (m.find_component(v), bnds) for v, bnds in valid_var_bounds.items()
        )

        self._orig_param_clone_param_map = pe.ComponentMap(
            (i, j) for i, j in zip(uncertain_params, self._uncertain_params)
        )
        self._clone_param_orig_param_map = pe.ComponentMap(
            (i, j) for i, j in zip(self._uncertain_params, uncertain_params)
        )

        assert self.config.method in FlexTestMethod
        if self.config.method == FlexTestMethod.active_constraint:
            self._build_active_constraint_model()
        elif self.config.method == FlexTestMethod.linear_decision_rule:
            self._build_flextest_with_dr()
        elif self.config.method == FlexTestMethod.relu_decision_rule:
            self._build_flextest_with_dr()
        elif self.config.method in {
            FlexTestMethod.sampling,
            FlexTestMethod.vertex_enumeration,
        }:
            self._build_sampling()
        else:
            raise ValueError(f"Unrecognized method: {self.config.method}")

    def _set_param_bounds(
        self, param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]]
    ):
        for orig_p, clone_p in self._orig_param_clone_param_map.items():
            p_lb, p_ub = param_bounds[orig_p]
            ndx = _VarIndex(clone_p, None)
            p_var = self._model.unc_param_vars[ndx]
            p_var.setlb(p_lb)
            p_var.setub(p_ub)

    def _build_active_constraint_model(self):
        build_active_constraint_flextest(
            m=self._model,
            uncertain_params=self._uncertain_params,
            param_nominal_values=self._param_nominal_values,
            param_bounds=self._max_param_bounds,
            valid_var_bounds=self._valid_var_bounds,
        )

    def _build_flextest_with_dr(self):
        build_flextest_with_dr(
            m=self._model,
            uncertain_params=self._uncertain_params,
            param_nominal_values=self._param_nominal_values,
            param_bounds=self._max_param_bounds,
            controls=self._controls,
            valid_var_bounds=self._valid_var_bounds,
            config=self.config,
        )

    def _build_sampling(self):
        _replace_uncertain_params(
            m=self._model,
            uncertain_params=self._uncertain_params,
            param_nominal_values=self._param_nominal_values,
            param_bounds=self._max_param_bounds,
        )
        _build_inner_problem(
            m=self._model,
            enforce_equalities=True,
            unique_constraint_violations=False,
            valid_var_bounds=None,
        )

    def _solve_maximization(
        self, param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]]
    ) -> FlexTestResults:
        self._set_param_bounds(param_bounds=param_bounds)

        opt = self.config.minlp_solver
        res = opt.solve(self._model)
        pe.assert_optimal_termination(res)

        results = FlexTestResults()
        results.max_constraint_violation = self._model.max_constraint_violation.value
        if results.max_constraint_violation > self.config.feasibility_tol:
            results.termination = FlexTestTermination.found_infeasible_point
        else:
            results.termination = FlexTestTermination.proven_feasible
        results.unc_param_values_at_max_violation = pe.ComponentMap()
        for key, v in self._model.unc_param_vars.items():
            results.unc_param_values_at_max_violation[
                self._clone_param_orig_param_map[key.var]
            ] = v.value
        return results

    def _solve_sampling(
        self, param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]]
    ) -> FlexTestResults:
        self._set_param_bounds(param_bounds=param_bounds)
        tmp = _perform_sampling(
            m=self._model,
            uncertain_params=self._uncertain_params,
            controls=self._controls,
            config=self.config.sampling_config,
        )
        sample_points, max_violation_values, control_values = tmp
        sample_points = pe.ComponentMap(
            (self._clone_param_orig_param_map[p], vals)
            for p, vals in sample_points.items()
        )

        results = FlexTestResults()
        max_viol_ndx = int(np.argmax(max_violation_values))
        results.max_constraint_violation = max_violation_values[max_viol_ndx]
        if results.max_constraint_violation > self.config.feasibility_tol:
            results.termination = FlexTestTermination.found_infeasible_point
        else:
            results.termination = FlexTestTermination.proven_feasible
        results.unc_param_values_at_max_violation = pe.ComponentMap()
        for key, vals in sample_points.items():
            results.unc_param_values_at_max_violation[key] = vals[max_viol_ndx]
        return results

    def solve(
        self, param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]]
    ) -> FlexTestResults:
        if self.config.method in {
            FlexTestMethod.active_constraint,
            FlexTestMethod.linear_decision_rule,
            FlexTestMethod.relu_decision_rule,
        }:
            return self._solve_maximization(param_bounds)
        elif self.config.method in {
            FlexTestMethod.sampling,
            FlexTestMethod.vertex_enumeration,
        }:
            return self._solve_sampling(param_bounds=param_bounds)
        else:
            raise ValueError(f"Unrecognized method: {self.config.method}")
