Custom Sampling
===========================================
With this method, users can explicitly define the distribution for the sampling of each input variable explicitly.
 
The ``pysmo.sampling.CustomSampling`` method carries out the user-defined sampling strategy. This can be done in two modes:

* The samples can be selected from a user-provided dataset, or
* The samples can be generated from a set of provided bounds.

We currently support three distributions options for sampling:

* "random", for sampling from a random distribution.
* "uniform", for sampling from a uniform distribution.
* "normal", for sampling from a normal (i.e. Gaussian) distribution.

.. warning::
    **A note on Gaussian-based sampling**

    To remain consistent with the other sampling methods and distributions, bounds are required for specifying normal distributions, rather than the mean (:math:`\bar{x}`) and standard deviation (:math:`\sigma`). For a normal distribution, 99.7% of the points/sample fall within three standard deviations of the mean. Thus, the bounds of the distribution ay be computed as:

    .. math::
        \begin{equation}
        LB = \bar{x} - 3\sigma
        \end{equation}

    .. math::
        \begin{equation}
        UB = \bar{x} + 3\sigma
        \end{equation}

    While almost all of the points generated will typically fall between LB and UB, a few points may be generated outside the bounds (as should be expected from a normal distribution). However, users can choose to enforce the bounds as hard constraints by setting the boolean option  **strictly_enforce_gaussian_bounds**  to True during initialization. In that case, values exceeding the bounds are replaced by new values generated from the distributions. However, this may affect  the underlying distribution.


Available Methods
------------------

.. autoclass:: idaes.core.surrogate.pysmo.sampling.CustomSampling
    :members: __init__, sample_points


