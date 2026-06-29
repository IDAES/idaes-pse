Applying Scaling Tools
======================

.. contents:: :local:
    :depth: 2

Basic Usage
-----------

.. Note::

  Modelers should endeavor to provide as much scaling information as possible before calling Scalers in order to provide as much information on your particular case as possible.

All Scalers in the Scaling Toolbox support a ``scale_model`` method that can be called to apply scaling to a model of the appropriate type. The ``scale_model`` method can be called as shown below:

.. code-block:: python

  # Import Scaler class from library
  from ideas.core.scaling import AutoScaler, set_scaling_factor
  
  # Set some scaling factors
  set_scaling_factor(my_model.my_var, 1e-3)

  # Create instance of Scaler object
  my_scaler = AutoScaler()

  # Apply Scaler to model
  my_scaler.scale_model(my_model)

Scaler Options
''''''''''''''

Many Scalers will support additional optional arguments which can be used to provide finer control over the scaling routine. See the documentation for the Scaler you are using for more details. Documentation for the core Scalers can be found here:

* :ref:`AutoScaler Class<reference_guides/scaling/autoscaler:AutoScaler Class>`
* :ref:`CustomScalerBase Class<reference_guides/scaling/custom_scaler:CustomScalerBase Class>`

Advanced Usage
--------------

In most cases, Scaler classes will have individual methods for scaling variables and constraints that can be called separately. Advanced modelers may wish to make use of these to gain even finer control over the scaling process, and may wish to experiment with mixing-and-matching routines from different Scalers (e.g., a modeler may wish to use the AutoScaler for variables but combine it with a more customized routine for constraint scaling from a Custom Scaler).

The CustomScalerBase class also contains a number of methods for common scaling approaches that that can be applied to individual variables and constraints, allowing advanced modelers to construct their own custom scaling routines.

For example, the simple example above can also be implemented by calling methods to scale the variables and constraints as shown below:

.. code-block:: python

  # Import Scaler class from library
  from ideas.core.scaling import AutoScaler, set_scaling_factor
  
  # Set some scaling factors
  set_scaling_factor(my_model.my_var, 1e-3)

  # Create instance of Scaler object
  my_scaler = AutoScaler()

  # Apply Scaler to model
  my_scaler.scale_variables_by_magnitude(my_model)
  my_scaler.scale_constraints_by_jacobian_norm(my_model)

