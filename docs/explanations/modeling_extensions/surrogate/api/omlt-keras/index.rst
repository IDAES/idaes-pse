.. _omlt:

OMLT: Optimization and Machine Learning Toolkit
===============================================

.. toctree::
    :maxdepth: 1

    omlt-keras-options

Keras is a deep learning framework that integrates with TensorFlow's structure for building and training artificial neural networks, and minimizes the number of user actions required to construct accurate networks. OMLT (Optimization and Machine Learning Toolkit) provides an interface to formulate machine learning models and import Keras or ONNX models as Pyomo blocks. The provided tools include an interface for accessing the Keras module (via the publicly available Python package) within IDAES flowsheets.

Basic Usage
-----------

Keras-OMLT main function is **keras_surrogate.KerasSurrogate**, which populates an IDAES `SurrogateObject` with the OMLT model. This object may then be passed directly to other IDAES methods for visualization or flowsheet integration (see the sections for Visualization and Examples below).

Data can be read in or simulated using available Python packages. The main arguments of the `keras_surrogate.KerasSurrogate` Python function are inputs and outputs. For example,

.. code-block:: python
  
  # import tensorflow
  import tensorflow as tf
  
  # selected settings for regression
  activation, optimizer, n_hidden_layers, n_nodes_per_layer = "tanh", "Adam", 2, 40
  loss, metrics = "mse", ["mae", "mse"]  

  # Create data objects for training using scalar normalization
  model = tf.keras.Sequential()
  model.add(
      tf.keras.layers.Dense(units=n_nodes_per_layer, input_dim=len(input_labels), activation=activation))
  # Create n hidden layers
  for i in range(1, n_hidden_layers):
      model.add(tf.keras.layers.Dense(units=n_nodes_per_layer, activation=activation))
  model.add(tf.keras.layers.Dense(units=len(output_labels)))

  # Train surrogate (calls optimizer on neural network and solves for weights)
  model.compile(loss=loss, optimizer=optimizer, metrics=metrics)
  # Loading the weights file
  mcp_save = tf.keras.callbacks.ModelCheckpoint(
      ".mdl_wts.hdf5", save_best_only=True, monitor="val_loss", mode="min")
  history = model.fit(
      x=x, y=y, validation_split=0.2, verbose=1, epochs=1000, callbacks=[mcp_save])

  keras_surrogate = KerasSurrogate(
      model,
      input_labels=list(input_labels),
      output_labels=list(output_labels),
      input_bounds=input_bounds,
      input_scaler=input_scaler,
      output_scaler=output_scaler,
  )

For an example on inputs, outputs, bounds and scalers see the `Autothermal Reformer Flowsheet Optimization with OMLT (TensorFlow Keras) Surrogate Object <https://github.com/IDAES/examples/blob/main/idaes_examples/notebooks/docs/surrogates/omlt/keras_flowsheet_optimization_src.ipynb>`_.

Saving and Loading OMLT-keras models
------------------------------------

The user may save their neural network objects by serializing to JSON, and load into a different script, notebook or environment. For example,

.. code-block:: python

  # To save a model
  model = keras_surrogate.save_to_folder("keras_surrogate")

  # To load a model
  keras_surrogate = KerasSurrogate.load_from_folder("keras_surrogate")

Visualizing Surrogate Model Results
-----------------------------------

For visualizing TensorFlow Keras neural networks via parity and residual plots, see :ref:`Visualizing Surrogate Model Results<explanations/modeling_extensions/surrogate/plotting/index:Visualizing Surrogate Model Results>`.


OMLT Example
------------

For an example of optimizing a flowsheet containing TensorFlow Keras neural networks utilizing the OMLT package, see the `Autothermal Reformer Flowsheet Optimization with OMLT (TensorFlow Keras) Surrogate Object <https://github.com/IDAES/examples/blob/main/idaes_examples/notebooks/docs/surrogates/omlt/keras_flowsheet_optimization_src.ipynb>`_.