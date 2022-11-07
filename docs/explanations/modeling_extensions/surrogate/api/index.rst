IDAES Surrogates API
====================

.. contents:: Contents
    :depth: 4

The IDAES Surrogates API provides a common interface for training surrogate modelings via supported tools and populating blocks with those models.

IDAES currently supports the following surrogate trainers:

.. toctree::
   :maxdepth: 1

   alamopy/index
   pysmo/index
   omlt-keras/index

API Structure and Basic Usage
-----------------------------

Users may call a Python wrapper for an available surrogate tool by passing variable data, labels and bounds, along with tool-specific settings. The API supports the workflow in the following sections.

Preparing Data For Training
^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, users must load or generate a dataset to train against, typically as a Pandas DataFrame with column labels. In best practice, the input variable columns are on the left and the output variable columns are on the right. For example, to read data from a CSV file:

.. code-block:: python

    # there exists a data file called data.csv
    import pd
    data = pd.read_csv(r'data.csv')  # load training data
    input_data, output_data = data.iloc[:, :2], data.iloc[:, 2:]  # 2 outputs
    input_labels, output_labels = input_data.columns, output_data.columns]

.. rubric:: Sampling

If desired, users may sample additional points or perform a training/validation split via :ref:`IDAES Sampling tools<explanations/modeling_extensions/surrogate/sampling/index:Sampling>` within PySMO. Note that sampling tools within PySMO are independent of surrogate training within PySMO, and may be called regardless of trainer selected. For example, to define training and validation datasets:

.. code-block:: python

    from idaes.core.surrogate.sampling.data_utils import split_training_validation
    n_data = data[input_labels[0]].size  # get the number of data points
    data_training, data_validation = split_training_validation(data, 0.8, seed=n_data)  # in this case, randomly separate 20% of the data for training and 80% of the data for validation

Training Surrogates
^^^^^^^^^^^^^^^^^^^

Once the training data is defined, users call the desired trainer method (`AlamoTrainer`, `PysmoPolyTrainer`, `PysmoRBFTrainer`, `PysmoKrigingTrainer`) to return a `SurrogateTrainer` object that is ready to train. For example, an `AlamoTrainer` object is created as below:

.. code-block:: python

    trainer = AlamoTrainer(input_labels=input_labels, output_labels=output_labels, training_dataframe=data_training)

Similarly, a PySMO RBF `PysmoRBFTrainer` object would be created as below:

.. code-block:: python

    trainer = PysmoRBFTrainer(input_labels=input_labels, output_labels=output_labels, training_dataframe=data_training)


Trainers allow their own configuration options, which alter the regression or file behavior and belong to the trainer's `CONFIG` block. For example, ALAMO supports a configuration option for allowing linear basis functions:

.. code-block:: python

    trainer.config.linfcns = True

Once the `SurrogateTrainer` is fully defined, calling `trainer.train_surrogate` performs regression and returns the model results as a dictionary object named `trainer._results`. This dictionary contains generated Pyomo model expressions in `trainer._results['Model']`.

Alternatively, users may train a TensorFlow Keras neural network model. Users should refer to the specific documentation in the table of contents at the top of this document for details on implementing each of the five supported training tools.

Building an IDAES Surrogate Object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add the model to an IDAES flowsheet, the `SurrogateTrainer` object needs to be transformed into an IDAES `SurrogateObject` object. This is done by calling an appropriate method (`AlamoSurrogate`, `PySMOSurrogate`, `KerasSurrogate`) and passing the generated surrogate expressions, along with variable labels and bounds. For example, an `AlamoSurrogate` object is created as below:

.. code-block:: python

    surr = AlamoSurrogate(expressions, input_labels, output_labels, input_bounds)

A similar object can be created for PySMO or Keras by replacing `AlamoSurrogate` with `PySMOSurrogate` or `KerasSurrogate`.

Flowsheet Integration
---------------------

After generating an IDAES `SurrogateObject`, the models may be combined with other IDAES method to assess and run the models.

Visualizing Surrogate Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: Visualizing Surrogate Model Results

Once the models are trained, users can leverage :ref:`IDAES Visualization tools<explanations/modeling_extensions/surrogate/plotting/index:Visualizing Surrogate Model Results>` to assess the model fits. For example, the following code generates 2D scatter, parity and residual plots for all output variables in the surrogate model and prints the plots to PDF:

.. code-block:: python

    from idaes.core.surrogate.plotting.sm_plotter import *
    surrogate_scatter2D(surr, data_validation, filename='scatter2D.pdf')
    surrogate_parity(surr, data_validation, filename='parity.pdf')
    surrogate_residual(surr, data_validation, filename='residual.pdf')

Saving/Loading Surrogates and Populating Surrogate Blocks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The IDAES Surrogates API offers simple methods to export/import trained surrogates and add them to IDAES flowsheets. After creating a `SurrogateObject`, calling `save_to_file()` on ALAMO and PySMO surrogates or `save_to_folder()` on Keras surrogates serializes the models to JSON. Then, users may load the surrogates into any script or flowsheet using the `load_from_file()` call for ALAMO and PySMO surrogates or `load_from_folder` for Keras surrogates. Keras neural networks are complex with model weights, training parameters and metadata, and as such require a folder of files rather than a single file when serialized. For example, to save and load an ALAMO `SurrogateObject`:

.. code-block:: python

    surr.save_to_file('alamo_surrogate.json', overwrite=True)
    alamo_surrogate = AlamoSurrogate.load_from_file('alamo_surrogate.json')

Then, the surrogate may be added to an IDAES `SurrogateBlock`:

.. code-block:: python
    
    m.fs.surrogate = SurrogateBlock()
    # define inputs and outputs as lists of Pyomo Var objects
    m.fs.surrogate.build_model(surrogate, input_vars=inputs, output_vars=outputs)

The surrogate must be able to find a list of Pyomo variable objects corresponding to the model input and output variables, and will automatically link the flowsheet variables to the model expressions. For a detailed demonstration of training, visualizing, optimizing and comparing available surrogate trainers, see the example [ML/AI Best Practices: "Selecting Surrogate Model Form/Size for Optimization"](https://github.com/IDAES/examples-pse/blob/main/src/Examples/SurrMod/FlowsheetOptimization/Best_practices_optimization.ipynb).