import pyomo.environ as pe
from pyomo.contrib import appsi


def assert_optimal_termination(results):
    if hasattr(results, "termination_condition"):
        assert results.termination_condition == appsi.base.TerminationCondition.optimal
    else:
        pe.assert_optimal_termination(results)
