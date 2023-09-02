Full-Factorial Sampling
===============================

The ``pysmo.sampling.UniformSampling`` method carries out Uniform (full-factorial) sampling. This can be done in two modes:

* The samples can be selected from a user-provided dataset, or
* The samples can be generated from a set of provided bounds.

Available Methods
------------------

.. autoclass:: idaes.core.surrogate.pysmo.sampling.UniformSampling
    :members: __init__, sample_points

References
------------
[1] Loeven et al paper titled "A Probabilistic Radial Basis Function Approach for Uncertainty Quantification"
https://pdfs.semanticscholar.org/48a0/d3797e482e37f73e077893594e01e1c667a2.pdf