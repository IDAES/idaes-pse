import pyomo.environ as pe
from typing import MutableMapping, Sequence
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.base.block import _BlockData
from pyomo.core.expr.numeric_expr import LinearExpression
from .dr_config import DRConfig
from pyomo.common.config import ConfigValue
from ..check_optimal import assert_optimal_termination


class LinearDRConfig(DRConfig):
    r"""
    A class for specifying options for constructing linear decision rules
    for use in the flexibility test problem.

    Attributes
    ----------
    solver: Union[Solver, OptSolver]
        The solver to use for building the linear decision rule (an LP solver).
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
        self.solver = self.declare(
            "solver", ConfigValue(default=pe.SolverFactory("ipopt"))
        )


def construct_linear_decision_rule(
    input_vals: MutableMapping[_GeneralVarData, Sequence[float]],
    output_vals: MutableMapping[_GeneralVarData, Sequence[float]],
    config: LinearDRConfig,
) -> _BlockData:
    n_inputs = len(input_vals)
    n_outputs = len(output_vals)

    res = pe.Block(concrete=True)
    res.output_set = pe.Set(initialize=list(range(n_outputs)))
    res.decision_rule = pe.Constraint(res.output_set)

    for out_ndx, (output_var, out_samples) in enumerate(output_vals.items()):
        trainer = pe.ConcreteModel()

        n_samples = len(out_samples)
        trainer.input_set = pe.Set(initialize=list(range(n_inputs)))
        trainer.sample_set = pe.Set(initialize=list(range(n_samples)))

        trainer.const = pe.Var()
        trainer.coefs = pe.Var(trainer.input_set)
        trainer.out_est = pe.Var(trainer.sample_set)

        obj_expr = sum(
            (trainer.out_est[i] - out_samples[i]) ** 2 for i in trainer.sample_set
        )
        trainer.objective = pe.Objective(expr=obj_expr)

        trainer.est_cons = pe.Constraint(trainer.sample_set)
        for ndx in trainer.sample_set:
            lin_coefs = [v[ndx] for k, v in input_vals.items()]
            lin_vars = list(trainer.coefs.values())
            lin_coefs.append(1)
            lin_vars.append(trainer.const)
            lin_coefs.append(-1)
            lin_vars.append(trainer.out_est[ndx])
            expr = LinearExpression(linear_coefs=lin_coefs, linear_vars=lin_vars)
            trainer.est_cons[ndx] = (expr, 0)

        results = config.solver.solve(trainer)
        assert_optimal_termination(results)

        lin_coefs = [v.value for v in trainer.coefs.values()]
        lin_vars = [v for v in input_vals.keys()]
        lin_coefs.append(-1)
        lin_vars.append(output_var)
        dr_expr = LinearExpression(
            constant=trainer.const.value, linear_coefs=lin_coefs, linear_vars=lin_vars
        )
        res.decision_rule[out_ndx] = (dr_expr, 0)

    return res
