Latin Hypercube Sampling (LHS)
===========================================
LHS is a stratified random sampling method originally developed for efficient uncertainty assessment. LHS partitions the parameter space 
into bins of equal probability with the goal of attaining a more even distribution of sample points in the parameter space that would be possible with pure random sampling.

 
The ``pysmo.sampling.LatinHypercubeSampling`` method carries out Latin Hypercube sampling. This can be done in two modes:

* The samples can be selected from a user-provided dataset, or
* The samples can be generated from a set of provided bounds.

Available Methods
------------------

.. autoclass:: idaes.core.surrogate.pysmo.sampling.LatinHypercubeSampling
    :members: __init__, sample_points

References
------------

[1] Loeven et al paper titled "A Probabilistic Radial Basis Function Approach for Uncertainty Quantification"
https://pdfs.semanticscholar.org/48a0/d3797e482e37f73e077893594e01e1c667a2.pdf

[2] Webpage on low discrepancy sampling methods:
http://planning.cs.uiuc.edu/node210.html

[3] Swiler, Laura and Slepoy, Raisa and Giunta, Anthony: "Evaluation of sampling methods in constructing response surface approximations"
https://arc.aiaa.org/doi/abs/10.2514/6.2006-1827
