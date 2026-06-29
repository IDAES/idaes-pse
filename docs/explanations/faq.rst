Troubleshooting and FAQ
=======================

Please visit the `IDAES discussions board on Github <https://github.com/IDAES/idaes-pse/discussions>`_ for more questions and answers.

.. contents::
    :depth: 2
    :local:

Frequently Asked Questions
--------------------------

Help! My model does not solve!
``````````````````````````````

The IDAES Team has limited resources and bandwidth for providing direct support for assisting users in debugging solver issues.  Before contacting the IDAES team for help, make sure you have tried all the steps below.

  1. First, check the output logs for any warnings or cautions as these may indicate that there are problems with your model or code that need to be addressed before you try to solve it. Also make sure that the code is making it to the solver call and that the issue does not lie with the model construction.
  2. Second, check the final solver status to understand how the solver is failing as this can give you hints as to what is wrong. If the solver status is infeasible, it generally indicates that the model is actually infeasible and that you need to fix something in your model. It is also helpful to look at the solver output (`solver.solve(model, tee=True)`) as there is often useful information in the solver traces.
  3. Next, try using the :ref:`IDAES Diagnostics workflow and toolbox <explanations/model_diagnostics/index:Model Diagnostics Workflow>` and ensure there are no structural or numerical issues in your model.
  4. If the model still fails to solve, try using the :ref:`Degeneracy Hunter <explanations/model_diagnostics/degeneracy_hunter:Degeneracy Hunter>` or :ref:`SVD <explanations/model_diagnostics/svd_analysis:SVD Analysis>` toolboxes to look for more complex issues.
  5. In some cases, it can help to try a different solver (if you have access to one). Different solvers use different approaches to solving the problem, and sometimes one will succeed where others fail (especially on mildly ill posed problems).
  6. If all else fails, see the “How do I get more help?” section below. When reaching out for help, please include the full output of the IDAES and solver logs, as well as any issues identified by the diagnostics toolbox.


Help! My model raises Exceptions during construction!
`````````````````````````````````````````````````````

Before contacting the IDAES team for assistance or filing a bug report, please do the following checks:

  1. First, check the error message you received and try to understand what it is telling you and where it is coming from (IDAES, Pyomo, or Python). Often, the trick to understanding the error message is not understanding what went wrong, but how it happened. If the error appears to be coming from Python or Pyomo, then it is most likely a bug somewhere in the code.
  2. Second, check the full Python traceback of the error to try to identify where in the code it is originating. This will help narrow down the possible issues and allow you to isolate the offending code for testing. It is always easier to debug a smaller model than a large flowsheet (especially for an IDAES developer who is not familiar with the flowsheet you are trying to solve).
  3. If you cannot find the source of the error yourself, try to make a minimum failing example which replicates the problem that you can share with the IDAES team, as this will greatly increase the chances of us having the time to test your example and isolate the issue.
  4. To get help, see the “How do I get more help?” section below. When reaching out for help, please include the full Python traceback and any warnings from the logger, as well as the minimum failing example.


How do I run some examples?
```````````````````````````

See :doc:`../tutorials/tutorials_examples`.

How do I get more help?
```````````````````````

Use the website to `register <https://idaes.org/register/>`_ for the IDAES support mailing list.
Then you can send questions to idaes-support@idaes.org. For more specific technical questions, we recommend
the `IDAES discussions board on Github <https://github.com/IDAES/idaes-pse/discussions>`_.

How can I cite IDAES?
`````````````````````

If you use this software in your work, please cite the version of the IDAES Toolset you used along with the following paper:

Lee, A., Ghouse, J. H., Eslick, J. C., Laird, C. D., Siirola, J. D., Zamarripa, M. A., Gunter, D., Shinn, J. H., Dowling, A. W., Bhattacharyya, D., Biegler, L. T., Burgard, A. P., & Miller, D. C. (2021). The IDAES process modeling framework and model library—Flexibility for process simulation and optimization. Journal of Advanced Manufacturing and Processing, 3. https://doi.org/doi/10.1002/amp2.10095

Each release of the IDAES Toolset is published on the `U.S. Department of Energy Office of Scientific and Technical Information webpage <https://www.osti.gov/search/semantic:idaes-pse>`_ and has its own DOI.

Troubleshooting
---------------

Missing win32api DLL
````````````````````

For Python 3.8 and maybe others, you can get an error when running Jupyter on Windows 10 about
missing the win32api DLL. There is a relatively easy fix::

  pip uninstall pywin32
  pip install pywin32==225

.. IMPORTANT::
   Python 3.8 is no longer officially supported.
