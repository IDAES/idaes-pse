Defining Property Packages
==========================

.. contents:: Contents 
    :depth: 2


Introduction
------------

In order to define a property package using the IDAES Generic Property Package Framework, users need to define a :ref:`Property Parameter Block<core/state_block:Physical Parameter Blocks>` in order to describe their material and its properties. The class should inherit from the IDAES `GenericParameterData` class and contain two methods;

1. `configure`, which defines the users selection of sub-models, and
2. `build`, which defines the parameters necessary for the selected property methods.

A basic outline of a user defined Property Parameter Block is shown below.

.. code-block:: python

    @declare_process_block_class("UserParameterBlock")
    class UserParameterData(GenericParameterData):
        def configure(self):
            # Set configuration options
        def build(self):
            # Define parameters

Users should populate the `configure` and `build` methods as discussed below.

Configure
---------

The 'configure` method is used to select the sub-models and methods to be used when constructing the `StateBlocks` associated with this property package within a process model. These options define the behavior of the property package, and allow users to customize the property package to their needs.

Configuration options are set by assigning an appropriate method to `self.configure.option_name` within the `configure` method. A full list of the available property options is given :ref:`here<property_models/general/generic_options:Configuration Options>`.

Build
-----

The `build` method is used to define and specify values for all the parameters associated with the property calculations. All property calculations depend on a set of empirically derived parameters to describe the behavior of the material. The list of parameters which need to be defined will depend upon the configuration options chosen, and the documentation for each method lists the expected parameters which need to be defined in this section. Users need only define those parameters required by the options they have chosen.

Property parameters can be defined as either Pyomo `Params` or `Vars` depending upon the application. Whilst `Params` would seem to be the logical choice, be aware that for parameter estimation problems, the parameters being estimated need to be defined as `Vars` so that the solver is free to vary them. 

.. note::

   If using `Params`, users should consider whether these should be `mutable` or not - `Params` that are not mutable have their value defined upon creation and this cannot be changed later.

   If using `Vars`, remember that you will need to fix the value unless you are trying to estimate the value of that parameter.

Property parameters need to have the correct set of indices and follow the naming convention laid out in the :ref:`standard naming conventions<standards:Standard Variable Names>` and described in the documentation for each property method. Property parameters are created using Pyomo code as shown below:

Param
^^^^^

.. code-block:: python

    self.parameter = Param([indices], initialize=value(s), mutable=True/False)

Var
^^^

.. code-block:: python

    self.parameter = Var([indices], initialize=value(s))
    self.parameter.fix()

Examples
--------

Examples of using the IDAES Generic Property Package Framework can be found in the `idaes/property_models/core/examples` folder.
