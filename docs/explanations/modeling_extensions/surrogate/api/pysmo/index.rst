PySMO: Python-based Surrogate Modeling Objects
==============================================

The PySMO toolbox provides tools for generating different conventional types of reduced order models with minimal tuning. It provides IDAES users with a set of surrogate modeling tools which supports flowsheeting and direct integration into an equation-oriented modeling framework.

.. rubric:: Sampling

PySMO provides tools necessary for sampling and surrogate model generation. For more information about its sampling tools, see :ref:`Sampling<explanations/modeling_extensions/surrogate/sampling/index:Sampling>`.

Surrogate Generation Tools
----------------------------

PySMO offers tools for generating three types of surrogates:

.. toctree::
   :maxdepth: 1

   pysmo_polyregression
   pysmo_radialbasisfunctions
   pysmo_kriging


Basic Usage
-----------

.. note::
   The IDAES surrogate API is a wrapper around the original PySMO surrogate interface.

PySMO's main functions are a set of trainers (**pysmo_surrogate.PysmoPolyTrainer**, **pysmo_surrogate.PysmoKrigingTrainer**, **pysmo_surrogate.PysmoRBFTrainer**) , which call PySMO to train surrogates from passed data, and **pysmo_surrogate.PysmoSurrogate**, which populates an IDAES `SurrogateObject` with the trained PySMO model. This object may then be passed directly to other IDAES methods for visualization or flowsheet integration (see the sections for Visualization and Examples below).

Data can be read in or simulated using available Python packages. The main arguments of the three Python functions are inputs and outputs, which are 2D arrays of data with associated variable labels. Once a trained surrogate object exists, `pysmo_surrogate.PysmoSurrogate` takes the model expressions, variable labels and input bounds as arguments. For example, for a Kriging model:

.. code-block:: python

  # After reading or generating a DataFrame object called `data_training`
  trainer = PysmoKrigingTrainer(input_labels=['x1', 'x2'], output_labels=['z1', 'z2'], training_dataframe=data_training) # Could be any of the trainers
  trainer.config.[PySMO Kriging Option] = [Valid Option Choice]  # see below for more details
  pysmo_surr_expr = trainer.train_surrogate()

  input_labels = trainer._input_labels
  output_labels = trainer._output_labels
  xmin, xmax = [0.1, 0.8], [0.8, 1.2]
  input_bounds = {input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))}

  pysmo_surr = PysmoSurrogate(pysmo_surr_expr, input_labels, output_labels, input_bounds)

where [PySMO Kriging Option] is a valid keyword argument that can be passed to the PySMO Kriging Python function to customize the model. Each PySMO model type requires takes different optional arguments; a list of arguments for each PySMO model type (Polynomial Regression, Radial Basis Functions or Kriging) may be found on its personalized page.

Saving and Loading PySMO models
--------------------------------
The user may save their trained surrogate objects by serializing to JSON, and load into a different script, notebook or environment. For example,

.. code-block:: python

  # To save a model
  model = pysmo_surr.save_to_file('pysmo_surrogate.json', overwrite=True)

  # To load a model
  surrogate = PysmoSurrogate.load_from_file('pysmo_surrogate.json')

Visualizing Surrogate Model Results
------------------------------------
For visualizing PySMO-trained surrogates via parity and residual plots, see :ref:`Visualizing Surrogate Model Results<explanations/modeling_extensions/surrogate/plotting/index:Visualizing Surrogate Model Results>`.


PySMO Examples
----------------

For an example of optimizing a flowsheet containing an PySMO-trained surrogate model, see the `Autothermal reformer flowsheet optimization example <https://github.com/IDAES/examples-pse/blob/main/src/Examples/SurrMod/FlowsheetOptimization/PySMO_flowsheet_optimization.ipynb>`_.
