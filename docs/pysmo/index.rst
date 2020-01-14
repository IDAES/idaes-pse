.. PySMO documentation master file, created by
   sphinx-quickstart on Mon Jan 13 09:00:21 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PySMO: Pyomo-based Surrogate Modelling Objects
==============================================

The PySMO toolbox provides tools for generating different types of reduced order models. It provides IDAES users with
a set of sampling surrogate modeling tools which supports flowsheeting and direct integration into an equation-oriented
modeling framework. It allows users to directly integrate reduced order models with algebraic high-fidelity process
models within an single IDAES flowsheet.


Tools
-----
PySMO provides two sets of tools necessary for:

* Sampling
* Surrogate model generation

Sampling
*********
The PySMO package offers five common sampling methods for one-shot design:

* Latin hypercube sampling
* Uniform grid (full-factorial) sampling
* Halton sampling
* Hammersley sampling
* Centroidal voronoi tessellation (CVT) sampling

The sampling methods are able to generate samples based from variable bounds or select samples from
a user-provided dataset.

The following code snippet shows basic usage of the package for generating samples from a set of bounds:

.. doctest::

   Required imports
   >>> from idaes.surrogates.pysmo import sampling as sp

   Declaration of lower and upper bounds of 3D space to be sampled
   >>> bounds = [[0, 0, 0], [1.2, 0.1, 1]]

   Initialize the Halton sampling method and generate 10 samples
   >>> space_init = sp.HaltonSampling(bounds_list, sampling_type='creation', number_of_samples=10)
   >>> samples = space_init.sample_points()


The following code snippet shows basic usage of the package for selecting sample points from an existing dataset:

.. doctest::

   Required imports
   >>> from idaes.surrogates.pysmo import sampling as sp
   >>> import pandas as pd

   Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   Initialize the CVT sampling method and generate 25 samples
   >>> space_init = sp.CVTSampling(xy_data, sampling_type='selection', number_of_samples=25)
   >>> samples = space_init.sample_points()


.. note::
   The results of the sampling process will be a Numpy array or Pandas dataframe, depending on the
   format of the input data.

Further information about the sampling tools and their input options may be found by accessing the individual
sampling methods.

Surrogate Generation
*********************

PySMO offers tools for generating three types of surrogates:

.. toctree::
   :maxdepth: 1

   pysmo_polyregression
   pysmo_radialbasisfunctions
   pysmo_kriging


.. toctree::
   :maxdepth: 2
   :caption: Contents:




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
