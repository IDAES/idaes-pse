Custom Sampling
===========================================
With this method, users can explicitly define the distribution for the sampling of each input variable explicitly.
 
The ``pysmo.sampling.CustomSampling`` method carries out the user-defined sampling strategy. This can be done in two modes:

* The samples can be selected from a user-provided dataset, or
* The samples can be generated from a set of provided bounds.

We currently support three distributions options for sampling:

* "random", for sampling from a random distribution
* "uniform", for sampling from a uniform distribution
* "normal", for sampling from a Gaussian distribution

Available Methods
------------------

.. autoclass:: idaes.core.surrogate.pysmo.sampling.CustomSampling
    :members: __init__, sample_points


