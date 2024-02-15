SVD Toolbox
===========

The IDAES SVD Toolbox is an advanced diagnostics tool for helping to identify scaling and degeneracy issues in models using Singular Value Decomposition and is accessible through the ``DiagnosticsToolbox``. Further discussion of how to use and interpret SVD results can be found :ref:`here<explanations/model_diagnostics/svd_analysis:SVD Analysis>`.

.. autoclass:: idaes.core.util.model_diagnostics.SVDToolbox
    :members:
    
SVD Callbacks
-------------

The SVD Toolbox supports callbacks to select the SVD analysis tool to use. Two callbacks are provided to make use of methods available in Scipy.

.. automodule:: idaes.core.util.model_diagnostics
    :noindex:
    :members: svd_dense, svd_sparse
