Scaling Workflow
================

.. contents:: :local:
    :depth: 1

General Approach
------------------

When scaling models, it is often difficult to know where to start. ``AutoScalers`` may appear to be attractive for this as they can scale an entire model (of any size) in one go and need minimal user input, but their narrow focus on the current model state is often insufficient for optimization purposes. Rather, it is generally necessary (and strongly encouraged) that modelers try to provide as much information about scaling as they can through their understanding of the system being modeled, and where possible the model structure.

As a starting point, the following workflow is recommended when scaling any model:

1. Understand as much about the model as possible, including expected input and output values and where possible the formulation of the constraints.
2. Make use of the Diagnostics Toolbox to ensure there are no structural issues and to identify potential scaling issues that must be resolved. This also provides a reference point for checking to see that your scaling factors are improving the state of the model. Modelers are encouraged to use these tools throughout the process to monitor their progress, however note that a partially scaled model will often have more issues than a completely unscaled model (this is often expected, and not necessarily a sign that you are going the wrong way). Of particular note are the ``display_variables_with_extreme_jacobians`` and ``display_constraints_with_extreme_jacobians`` methods (as well as the ``SVDToolbox`` for advanced users).
3. Start by scaling those variables you have the most information about – these will generally be variables process inputs, design and operating conditions, etc.
4. Working from what you already know, try to project expected scaling for other variables, repeating as necessary.
5. Once you have established scaling for all the variables (or as many as you can), start scaling constraints in a similar fashion (start with what you understand best). Make use of the scaling methods provided by the ``CustomScalerBase`` class to assist you in this.

Scaling Hierarchical Models
---------------------------

The hierarchical nature of IDAES models adds an additional challenge to scaling workflows, but also provides opportunities for modularization and the creation of Scalers dedicated to scaling specific sub-sections of a model (e.g., unit and property models).

Flowsheet Scaling Workflow
''''''''''''''''''''''''''

When scaling flowsheets, an approach similar to initialization can be used, where a modeler starts at the inlet (feed) to their process and scales each unit in a sequential fashion, propagating scaling factors along the connecting ``Arcs`` as they go (see ``propagate_state_scaling`` method in the ``CustomScalerBase`` class). Each unit model in the process can then be scaled in isolation applying Scalers suited to that unit model. Recycle loops bring an additional challenge to this, however careful consideration of the expected state of the recycle can help guide this process, or traditional iterative techniques can be applied.

Scaling Unit and Property Models
''''''''''''''''''''''''''''''''

Unit models in turn are hierarchical constructs, which depend on sub-models such as the ``StateBlocks`` used to calculate physical properties. Each of these sub-models can have Scalers of their own, and thus a hierarchical approach to scaling can be applied where the unit model first calls a Scaler for the inlet state, then propagates the scaling to the outlet state and calls a Scaler for that StateBlock, and then finally uses this information to inform the scaling of the unit level variables and constraints (including those in any control volume(s)). The ``CustomScalerBase`` class contains a number of methods and tools to assist in this process, or the experienced modeler may wish to perform these steps manually.


