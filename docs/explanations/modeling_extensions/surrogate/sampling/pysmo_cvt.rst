Centroidal voronoi tessellation (CVT) sampling
===================================================
In CVT, the generating point of each Voronoi cell coincides with its center of mass; CVT sampling locates the design samples at the centroids of each Voronoi cell in
the input space. CVT sampling is a geometric, space-filling sampling method which is similar to k-means clustering in its simplest form.

The ``pysmo.sampling.CVTSampling`` method carries out CVT sampling. This can be done in two modes:

* The samples can be selected from a user-provided dataset, or
* The samples can be generated from a set of provided bounds.

The CVT sampling algorithm implemented here is based on McQueen's method which involves a series of random sampling and averaging steps, 
see http://kmh-lanl.hansonhub.com/uncertainty/meetings/gunz03vgr.pdf.

Available Methods
------------------

.. autoclass:: idaes.core.surrogate.pysmo.sampling.CVTSampling
    :members: __init__, sample_points

References
------------
[1] Loeven et al paper titled "A Probabilistic Radial Basis Function Approach for Uncertainty Quantification"
https://pdfs.semanticscholar.org/48a0/d3797e482e37f73e077893594e01e1c667a2.pdf

[2] Centroidal Voronoi Tessellations: Applications and Algorithms by Qiang Du, Vance Faber, and Max Gunzburger 
https://doi.org/10.1137/S0036144599352836

[3] D. G. Loyola, M. Pedergnana, S. G. García, "Smart sampling and incremental function learning for very large high dimensional data"
https://www.sciencedirect.com/science/article/pii/S0893608015001768?via%3Dihub

