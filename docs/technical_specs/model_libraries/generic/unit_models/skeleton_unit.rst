Skeleton Unit Model
===============================

The IDAES Skeleton Unit model represents a bare bones unit (hence, the name skeleton) that lets the user define a unit model with a custom set
of variables and constraints that best suits their needs but that can be connected to other models from the unit model library. This unit can be used
to either define a surrogate for a complex unit model or to define a first-principles based model that does not use IDAES
ControlVolumes or property packages. To facilitate connection to other models in the unit model library, the user needs to add ports and populate with variables that match
the port members of the models this will be connected to. To enable this, a special method is provided to enable adding ports with a specified set of variables.

The two methods available for this unit are:

1. ``add_ports``: Allows the user to define the inlet and outlet ports for the model and populate the port with a set of variables.
2. ``initialize``: A default initialize method that will attempt to solve the model as defined or that can be substituted with a custom callback function provided by the user.

Degrees of Freedom
------------------

There are no variables or constraints as this is a bare minimum skeleton that will be populated with variables and constraints
as defined by the user. Therefore, the number of degrees of freedom will be determined based on the model the user implements.

Model Structure
---------------

Given that this unit is only a skeleton, the model structure is whatever the user intends to be. For example, the
model has no restrictions on number of inlets and outlets.

Variables
---------

Variables will be defined by the user. Note that when defining variables, they will need to be indexed by the flowsheet time index irrespective if it is steady-state
or dynamic. This is to ensure that port members can be connected as there is still a 0 time index even when the flowsheet is steady-state.

Custom Callback
---------------

Note that the default value for the ``initializer`` config argument is set to trigger the ``_default_initialize`` method. To use a custom callback instead, the user needs to first define a function
and set the ``obj.config.initializer`` to the custom callback function they defined. In the ``initialize`` method, after checking that degrees of freedom is zero, the following is triggered:

.. code-block:: Python

   self.config.initializer(self, opt=opt, init_log=init_log, solve_log=solve_log)

Note that the function signature in the above line expects an instance, ``opt`` which is the solver (default is set to ipopt), ``init_log`` and ``solve_log`` which are the idaes logger objects for initialization
and solve. The user can either pass these arguments in their custom function or simply do the following:

.. code-block:: Python

   def custom_func(unit, **kwargs):

       # custom initialization sequence
       # unit - instance of SkeletonUnitModel

The `unit` refers to the instance and `**kwargs` will handle the additional function arguments that will be passed in the ``initialize`` method.

.. module:: idaes.generic_models.unit_models.skeleton_model

SkeletonUnitModel Class
-----------------------

.. autoclass:: SkeletonUnitModel
  :members:

SkeletonUnitModelData Class
---------------------------

.. autoclass:: SkeletonUnitModelData
  :members:
