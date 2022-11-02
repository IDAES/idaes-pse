.. _sampling_details:

More Information about PySMO's Sampling Methods
===================================================
The sampling methods are able to generate samples based from variable bounds or select samples from a user-provided dataset. To use any of the method, the class is first initialized with the required parameters, and then the ``sample_points`` method is called.

.. note::
   The results of the sampling process will be a Numpy array or Pandas dataframe, depending on the format of the input data.

Examples
=========

Sample generation from scratch
-------------------------------


The following code snippet shows basic usage of the package for generating samples from a set of bounds:

.. code:: python

   # Required imports
   >>> from idaes.core.surrogate.pysmo import sampling as sp

   # Declaration of lower and upper bounds of 3D space to be sampled
   >>> bounds_list = [[0, 0, 0], [1.2, 0.1, 1]]

   # Initialize the Halton sampling method and generate 10 samples
   >>> space_init = sp.HaltonSampling(bounds_list, sampling_type='creation', number_of_samples=10)
   >>> samples = space_init.sample_points()


Sample selection from a dataset
-------------------------------
The following code snippet shows basic usage of the package for selecting sample points from an existing dataset:

.. code:: python

   # Required imports
   >>> from idaes.core.surrogate.pysmo import sampling as sp
   >>> import pandas as pd

   # Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   # Initialize the CVT sampling method and generate 25 samples
   >>> space_init = sp.CVTSampling(xy_data, sampling_type='selection', number_of_samples=25)
   >>> samples = space_init.sample_points()

.. note::
  In the above case, it is assumed that there is only one output variable AND it is in the last column of the dataset. When there are multiple output columns, the xlabels and ylabels options must be specified (see below example).


The following code snippet shows basic usage of the package for sub-sampling a dataset with multiple output variables:

.. code:: python

   # Required imports
   >>> from idaes.core.surrogate.pysmo import sampling as sp
   >>> import pandas as pd

   # Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   # Initialize the CVT sampling method and generate 25 samples
   >>> space_init = sp.CVTSampling(xy_data, sampling_type='selection', number_of_samples=25, xlabels=['x1', 'x2'], ylabels=['y1', 'y2'])
   >>> samples = space_init.sample_points()

.. warning::
   The user must take care to be explicit when dealing with the subsampling of multi-output data:

   * When both ``xlabels`` and ``ylabels`` are specified, any columns not present in them will be dropped from the dataset.

   * If only one of ``xlabels`` or ``ylabels`` is specified, all the columns not present in the provided list are automatically assigned to the unspecified option (``ylabels`` or ``xlabels``). 

.. note::
  If the input data to be sub-sampled is in a numpy array, ``xlabels`` and ``ylabels`` should be a list of the column numbers.


Characteristics of sampling methods available in PySMO
---------------------------------------------------------

.. list-table:: Characteristics of the different sampling methods
   :widths: 15 15 15 20 20 15
   :header-rows: 1

   * -
     - Deterministic
     - Stochastic
     - Low-discrepancy
     - Space-filling
     - Geometric
   * - LHS
     -
     - :math:`\checkmark`
     -
     -
     - :math:`\checkmark`
   * - Full-factorial
     - :math:`\checkmark`
     -
     -
     -
     - :math:`\checkmark`
   * - Halton
     - :math:`\checkmark`
     -
     - :math:`\checkmark`
     -
     -
   * - Hammersley
     - :math:`\checkmark`
     -
     - :math:`\checkmark`
     -
     -
   * - CVT
     - :math:`\checkmark`
     -
     -
     - :math:`\checkmark`
     - :math:`\checkmark`
