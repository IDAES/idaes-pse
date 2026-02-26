Diagnosing Scaling Issues
=========================

As mentioned in the previous section, the :ref:`IDAES Diagnostics Toolbox<explanations/model_diagnostics/index:Model Diagnostics Workflow>` contains a number of methods that can be used to help identify potential scaling issues. Scaling issues depend on the numerical state of the model, an thus the `report_numerical_issues` method is the place to start when looking for scaling issues.

Some common signs of poor scaling which can be seen in the numerical issues report include:

* large Jacobian Condition Numbers (>1e10),
* variables with very large or very small values,
* variables with values close to zero,
* variables with values close to a bound (conditionally),
* constraints with mismatched terms,
* constraints with potential cancellations,
* variables and constraints with extreme Jacobian norms,
* extreme entries in Jacobian matrix (conditionally).

If you see any of these warnings in the model diagnostics output, it is a sign that you have potential scaling issues which should be investigated in order to improve the performance, robustness and accuracy of your model.

.. Note::
  Not all scaling issues can be resolved through the application of scaling factors. In some cases, such as constraints with mismatched terms or possible cancellations, the constraint itself may be inherently poorly posed and thus may need to be refactored to resolve the scaling issue.
