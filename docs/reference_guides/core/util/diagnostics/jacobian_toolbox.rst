Jacobian Analysis Toolbox
=========================

The IDAES Jacobian Analysis Toolbox is an advanced diagnostics tool for helping to identify scaling, degeneracy and ill-conditioning issues in model Jacobians and is accessible through the ``DiagnosticsToolbox``.

Further discussion of how to use and interpret the results of the Jacobian Analysis tools can be found :ref:`here<explanations/model_diagnostics/jacobian_analysis:Jacobian Analysis>`.

.. autoclass:: idaes.core.util.model_diagnostics.JacobianAnalysisToolbox
    :members:
    
SVD Callbacks
-------------

The Jacobian Analysis Toolbox supports callbacks to select the SVD analysis tool to use. Two callbacks are provided to make use of methods available in Scipy.

.. automodule:: idaes.core.util.model_diagnostics
    :noindex:
    :members: svd_dense, svd_sparse
