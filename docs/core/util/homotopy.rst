Homotopy Meta-Solver
====================

The IDAES homotopy meta-solver is useful for cases where a user has a feasible solution to a well-defined (i.e. square) problem at one set of conditions (i.e. value of fixed variables), and wishes to find a feasible solution to the same problem at a different set of conditions. In many situations this can be achieved by directly changing the values of the fixed variables to their new values and solving the problem, but cases exist where this is challenging. Homotopy solvers try to find a feasible path to the new solution by taking smaller steps in the value of the fixed variables to progressively find a solution at the new point.

.. note:: A homotopy solver should not be considered a fix to a poorly posed or
       ill-conditioned problem, and users should first consider whether their
       problem can be reformulated for better performance.

Homotopy Routine
----------------

The IDAES homotopy routine starts from a feasible solution to the problem at the initial values for the fixed variables (:math:`v_0`) and a set of target values for these (:math:`t`). The routine then calculates a set of new values for the fixed variables during the first homotopy evaluation based on an initial step size :math:`s_0` such that:

.. math:: v_1 = t \times s_0 + v_0 \times (1-s_0)

The problem is then passed to Ipopt to try to find a solution at the current values for the fixed variables. Based on the success or failure of the solver step, the following occurs:

1. If the solver returns an optimal solution, the step is accepted and the solution to the current state of the model is saved (to provide a feasible point to revert to in case a future step fails). If the current meta-solver progress is 1 (i.e. it has converged to the target values), the meta-solver terminates otherwise the meta-solver progress (:math:`p_i`) is then updated, :math:`p_i = p_{i-1} + s_i`, and the size of the next homotopy step is then calculated based on an adaptive step size method such that:

.. math:: s_{i+1} = s_i \times \left(1 + a \times \left[\frac{I_t}{I_a}-1\right]\right)

where :math:`I_a` is the number of solver iterations required in the current homotopy step, :math:`I_t` is the desired number of solver iterations per homotopy step (an input parameter to the homotopy routine) and :math:`a` is a step size acceleration factor (another input parameter). As such, the size of the homotopy step is adjusted to try to achieve a desired number of solver iterations per step as a proxy for difficulty in solving each step. If new step would overshoot the target values, then the step size is cut back to match the target values. The user can also specify a maximum and/or minimum size for the homotopy which can be used to limit the homotopy step.

A new set of values for the fixed variables is calculated using :math:`v_{i+1} = t \times (p_i+s_{i+1}) + v_0 \times (1-(p_i+s_{i+1}))` and the process repeated.

2. If the solver fails to find an optimal solution (for any reason), the current step is rejected and solution to the previous successful step is reloaded. If the last homotopy step was equal to the minimum homotopy step size, the meta-solver terminates, otherwise, a reduced homotopy step is calculated using:

.. math:: s_{i+1} = s_i \times c

where :math:`c` is a step cut factor (an input parameter between 0.1 and 0.9). If the new step homotopy step is less than the minimum homotopy step size, the minimum step is used instead.

A new set of fixed variable values are then calculated and another attempt to solve the problem is made.

Possible Termination Conditions
-------------------------------

The homotopy meta-solver has the following possible termination conditions (using the Pyomo `TerminationCondition` Enum):

* `TerminationCondition.optimal` - meta-solver successfully converged at the target values for the fixed variables.
* `TerminationCondition.other` - the meta-solver successfully converged at the target values for the fixed variables, but with regularization of during final step. Users are recommended to discard this solution.
* `TerminationCondition.minStepLength` - the meta-solver was unable to find a feasible path to the target values, as the solver failed to find a solution using the minimum homotopy step size.
* `TerminationCondition.maxEvaluations` - the meta-solver terminated due to reaching the maximum allowed number of attempted homotopy steps
* `TerminationCondition.infeasible` - could not find feasible solution to the problem at the initial values for the fixed variables.

Available Methods
^^^^^^^^^^^^^^^^^

.. automodule:: idaes.core.util.homotopy
    :members:

