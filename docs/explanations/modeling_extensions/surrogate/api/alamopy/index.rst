.. index::
    pair: alamo;alamopy


ALAMOPY: ALAMO Python
======================

.. toctree::
    :maxdepth: 1

    alamopy-options

The purpose of ALAMOPY (Automatic Learning of Algebraic MOdels PYthon wrapper) is to provide a wrapper for the software ALAMO which generates algebraic surrogate models of black-box systems for which a simulator or experimental setup is available. Consider a system for which the outputs **z** are an unknown function **f** of the system inputs **x**.  The software identifies a function **f**, i.e., a relationship between the inputs and outputs of the system, that best matches data (pairs of **x** and corresponding **z** values) that are collected via simulation or experimentation.

Basic Usage
-----------

ALAMOPY's main functions are **alamopy.AlamoTrainer**, which calls ALAMO to train surrogates from passed data, and **alamopy.AlamoSurrogate**, which populates an IDAES `SurrogateObject` with the ALAMO model. This object may then be passed directly to other IDAES methods for visualization or flowsheet integration (see the sections for Visualization and Examples below).

Data can be read in or simulated using available Python packages. The main arguments of the `alamopy.AlamoTrainer`` Python function are inputs and outputs, which are 2D arrays of data with associated variable labels. Once a trained surrogate object exists, `alamopy.AlamoSurrogate` takes the model expressions, variable labels and input bounds as arguments. For example,

.. code-block:: python

  # after reading or generating a DataFrame object called `data_training`
  trainer = AlamoTrainer(input_labels=['x1', 'x2'], output_labels=['z1', 'z2'], training_dataframe=data_training)
  trainer.config.[Alamo Option] = [Valid Option Choice]  # see below for more details
  success, alm_surr, msg = trainer.train_surrogate()

  surrogate_expressions = trainer._results['Model']
  input_labels = trainer._input_labels
  output_labels = trainer._output_labels
  xmin, xmax = [0.1, 0.8], [0.8, 1.2]
  input_bounds = {input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))}

  alm_surr = AlamoSurrogate(surrogate_expressions, input_labels, output_labels, input_bounds)

where [Alamo Option] is a valid keyword argument that can be passed to the ALAMO Python function subject a range of available [Valid Option Choices] values to customize the basis function set, names of output files, and other options available in ALAMO.

User may save their trained surrogate objects by serializing to JSON, and load into a different script, notebook or environment. For example,

.. code-block:: python

  # to save a model
  model = alm_surr.save_to_file('alamo_surrogate.json', overwrite=True)

  # to load a model
  surrogate = AlamoSurrogate.load_from_file('alamo_surrogate.json')

Options for *alamopy.AlamoTrainer*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: ALAMOPY.ALAMO Options

Below are some common arguments users may pass to ALAMO through `alamopy.AlamoTrainer` that govern the behavior of the regression procedure, solver, file characteristics and other options. See :ref: `ALAMOPY.ALAMO Options<explanations.modeling_extensions/surrogate/alamopy/alamopy-options:ALAMOPY.ALAMO Options>` for more details on valid options and values.

* input_labels - list of strings to label the input variables
* output_labels - list of strings to label the output variables
* logfcns, expfcns, cosfcns, sinfcns, linfcns, constant - these are '0-1' options to activate these functions
* monomialpower, multi2power, multi3power, ratiopower - list of terms to be used in the respective basis functions
* maxterms - maximum number of basis terms allowed in the "best fit" model
* filename - path to ALAMO .alm input file
* overwrite_files - True/False flag whether to overwrite files
* modeler - integer 1-8 determines the choice of fitness metric
* screener - integer 0-2 determines if and how the size of the model will be minimized during training
* solvemip - '0-1' option that will force the solving of the .gms file

Visualization
-------------

.. rubric:: Visualizing Surrogate Model Results

For visualizing ALAMO-trained surrogates via parity and residual plots, see :ref:`Visualizing Surrogate Model Results<explanations/modeling_extensions/surrogate/plotting/index:Visualizing Surrogate Model Results>`.

Custom Basis Functions
----------------------

Similar to ALAMO, there are advanced capabilities for customization and constrained regression facilitated by methods in ALAMOPY including custom basis functions. These methods interact with the regression using the ALAMO option file.

Custom basis functions can be added to the built-in functions to expand the functional forms available. To use this advanced capability in ALAMOPY, the following function is called. Note it is necessary to use the xlabels assigned to the input parameters.

.. code-block:: python
  
  trainer.config.custom_basis_functions = ["x1^2 * x2^2", "...", "..." ...]

.. rubric:: ALAMOPY: ALAMO Python

The prior form of ALAMOPY contains several advanced capabilities that are not yet supported by the new framework. See :ref: `ALAMOPY.ALAMO Options<explanations.modeling_extensions/surrogate/alamopy/alamopy-depr:ALAMOPY: ALAMO Python>` for more details on these methods.

ALAMOPY Examples
----------------

For an example of optimizing a flowsheet containing an ALAMO-trained surrogate model, see [Autothermal Reformer Flowsheet Optimization with ALAMO Surrogate Object](https://github.com/IDAES/examples-pse/blob/main/src/Examples/SurrMod/FlowsheetOptimization/ALAMO_flowsheet_optimization.ipynb).