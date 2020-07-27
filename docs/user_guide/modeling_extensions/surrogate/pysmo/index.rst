.. PySMO documentation master file, created by
   sphinx-quickstart on Mon Jan 13 09:00:21 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PySMO: Python-based Surrogate Modelling Objects
=====================================================

The PySMO toolbox provides tools for generating different types of reduced order models. It provides IDAES users with
a set of surrogate modeling tools which supports flowsheeting and direct integration into an equation-oriented
modeling framework. It allows users to directly integrate reduced order models with algebraic high-fidelity process
models within an single IDAES flowsheet.

PySMO provides two sets of tools necessary for sampling and surrogate model generation.

Surrogate Generation
---------------------------

PySMO offers tools for generating three types of surrogates:

.. toctree::
   :maxdepth: 1

   pysmo_polyregression
   pysmo_radialbasisfunctions
   pysmo_kriging

Sampling
----------
The PySMO package offers five common sampling methods for one-shot design:

.. toctree::
   :maxdepth: 1

   pysmo_lhs
   pysmo_uniform
   pysmo_halton
   pysmo_hammersley
   pysmo_cvt
   pysmo_sampling_properties


Further information about the sampling tools and their input options may be found by accessing the individual
sampling methods. Examples and details of the characteristics of the sampling methods may be found at
:ref:`sampling_details`.