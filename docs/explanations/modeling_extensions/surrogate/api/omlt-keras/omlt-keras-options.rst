keras_surrogate Options
=======================

This page lists in more detail the kerras OMLT options.

.. contents::
    :depth: 2

Installing OMLT
----------------

OMLT (Optimization and Machine Learning Toolkit) is an optional dependency: and specific examples through a user manual and installation guide. Alternatively, users may directly access the user guide here: https://omlt.readthedocs.io/en/latest/index.html.

More details on OMLT options may be found in the user guide documentation linked above. If users encounter specific error codes while running the OMLT tool in IDAES, the user guide contains detailed descriptions of each termination condition and error message.

Basic keras surrogate options
-----------------------------

Data Arguments
^^^^^^^^^^^^^^

The following arguments are required by the `KerasSurrogate` method:

* **keras_model**: this is the Keras Sequential model that will be loaded. Note that specialized layers may not be supported at this time
* **input_labels**: user-specified labels given to the inputs
* **output_labels**: user-specified labels given to the outputs
* **input_bounds**: minimum/maximum bounds for each input variable to constraint 
* **input_scaler**: the scaler to be used for the inputs
* **output_scaler**: the scaler to be used for the outputs

.. code-block:: python

  .. keras_model = model
  .. input_labels = input_data.columns
  .. output_labels = output_data.columns
  .. input_bounds = {input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))}
  .. input_scaler = OffsetScaler.create_normalizing_scaler(input_data)
  .. output_scaler = OffsetScaler.create_normalizing_scaler(output_data)

Provided Formulations
^^^^^^^^^^^^^^^^^^^^^

OMLT can formulate what we call full-space and reduced-space neural network representations using the **FULL_SPACE** object (for full-space) and **REDUCED_SPACE** object (for reduced-space). The full space formulation supports non-smooth **RELU_BIGM** activation functions and OMLT uses **RELU_COMPLEMENTARITY** to specify that ReLU activation functions should be formulated using complementarity conditions. For example,

.. code-block:: python

  m.fs.surrogate.build_model(
    keras_surrogate,
    formulation=KerasSurrogate.Formulation.FULL_SPACE,
    input_vars=inputs,
    output_vars=outputs)


OMLT Layers
^^^^^^^^^^^

* **ConvLayer2D**: two-dimensional convolutional layer
* **DenseLayer**: dense layer implementing output = activation(dot(input, weights) + biases)
* **IndexMapper**: map indexes from one layer to the other
* **InputLayer**: the first layer in any network
* **Layer2D**: abstract two-dimensional layer that downsamples values in a kernel to a single value
* **PoolingLayer2D**: two-dimensional pooling layer


Layer Functions
^^^^^^^^^^^^^^^

* **full_space.full_space_conv2d_layer**
* **full_space.full_space_dense_layer**
* **full_space.full_space_maxpool2d_layer** 
* **reduced_space.reduced_space_dense_layer**
* **partition_based.default_partition_split_func**
* **partition_based.partition_based_dense_relu_layer**


Activation Functions
^^^^^^^^^^^^^^^^^^^^

* **linear**: applies the linear activation function
* **sigmoid**: applies the sigmoid function
* **softplus**: applies the softplus function
* **tanh**: applies the tanh function
