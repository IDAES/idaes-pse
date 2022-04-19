Halton Sampling
==================
Halton sampling is a low-discrepancy sampling method. It is a deterministic sampling method based on the Halton sequence, a sequence constructed by a set of co-prime bases. The Halton
sequence is an n-dimensional extension of the Van der Corput sequence; each individual Halton sequence is based on a radix inverse function defined on a prime number.

The ``pysmo.sampling.HaltonSampling`` method carries out Halton sampling. This can be done in two modes:

* The samples can be selected from a user-provided dataset, or
* The samples can be generated from a set of provided bounds.

The Halton sampling method is only available for low-dimensional problems :math:`n \leq 10`. At higher dimensions, the performance of the sampling method has been shown to degrade.

Available Methods
------------------

.. autoclass:: idaes.core.surrogate.pysmo.sampling.HaltonSampling
    :members: __init__, sample_points

References
------------
[1] Loeven et al paper titled "A Probabilistic Radial Basis Function Approach for Uncertainty Quantification"
https://pdfs.semanticscholar.org/48a0/d3797e482e37f73e077893594e01e1c667a2.pdf

[2] Webpage on low discrepancy sampling methods: 
http://planning.cs.uiuc.edu/node210.html
