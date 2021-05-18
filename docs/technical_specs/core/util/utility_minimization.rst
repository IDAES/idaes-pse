Utility minimization
====================

The IDAES Process Modeling Framework includes support for incorporating utility
minimization into a flowsheet to allow for calculation and optimization of utilities.
The formulation is implemented from M. A. Duran, I. E. Grossmann, Simultaneous
optimization and heat integration of chemical processes, AIChE Journal (1986).

This model is effective because the combinatorial search for the pinch candidate
is cast as the inequality constraint, which is efficiently handled
in equation oriented process optimization.

To create the constraints and utility variables ``Qs`` and ``Qw`` we utilize the function:

.. autofunction:: idaes.core.util.utility_minimization.min_utility

To effectively utilize the formulation, we need to add an objective function to minimize
the new variables, or a cost associated with them.

To generate composite curves, we must utilize the function ``heat_ex_data`` show bellow
to extract the heat exchanger and utility data.

.. autofunction:: idaes.core.util.utility_minimization.heat_ex_data

We can just run the function utilizing the data class previously obtained to print
composite curves with ``heat_ex_data``:

.. autofunction:: idaes.core.util.utility_minimization.generate_curves
