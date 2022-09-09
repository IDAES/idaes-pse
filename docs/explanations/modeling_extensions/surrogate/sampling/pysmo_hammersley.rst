Hammersley Sampling
======================
Hammersley sampling is a low-discrepancy sampling method based on the Hammersley sequence. The Hammersley sequence is the same as the Halton sequence 
except in the first dimension where points are located equidistant from each other.

The ``pysmo.sampling.HammersleySampling`` method carries out Hammersley sampling. This can be done in two modes:

* The samples can be selected from a user-provided dataset, or
* The samples can be generated from a set of provided bounds.

The Hammersley sampling method is only available for low-dimensional problems :math:`n \leq 10`. At higher dimensions, the performance of the sampling method has been shown to degrade.

Available Methods
------------------

.. autoclass:: idaes.core.surrogate.pysmo.sampling.HammersleySampling
    :members: __init__, sample_points

References
------------
[1] Loeven et al paper titled "A Probabilistic Radial Basis Function Approach for Uncertainty Quantification"
https://pdfs.semanticscholar.org/48a0/d3797e482e37f73e077893594e01e1c667a2.pdf

[2] Webpage on low discrepancy sampling methods: 
http://planning.cs.uiuc.edu/node210.html

[3] Holger Dammertz's webpage titled "Hammersley Points on the Hemisphere" which discusses  Hammersley point set generation in two dimensional spaces, 
http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html
