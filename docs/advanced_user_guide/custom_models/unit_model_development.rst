Developing Custom Unit Models
=============================

.. contents:: :local:

UnitModelBlockData
------------------

The starting point for all unit models within the IDAES Process Modeling Framework is the `UnitModelBlockData` base class. This class contains a number of methods to assist users with creating new unit models including:

* checking and validating the “dynamic” and “has_holdup” configuration arguments to ensure consistency,
* adding `Port` objects to the model,
* create simple material balances between states (equal flow of each component in each phase), and,
* a method for initializing simple unit models.

More details on the `UnitModelBlockData` class can be found in the :ref:`technical specifications<technical_specs/core/unit_model_block:Unit Model Class>`.

Unit Model Configuration
------------------------

Configuration arguments in unit models are used to allow users of the model to configure the model to suit the needs of the flowsheet they are simulating. The most common aspects that need to be configured are:

* whether the model should be dynamic or steady-state, and
* what property package to use when calculating thermophysical properties.

The `UnitModelBlockData` class contains a simple configuration block which includes two configuration arguments; “dynamic” and “has_holdup”. These arguments are required for any model which is expected to be used in both steady-state and dynamic flowsheets and are used to determine whether accumulation and holdup terms should be constructed and included in the material balance equations. There are some situations whoever where a model is inherently steady-state (even if it is included in a dynamic flowsheet), notably those where outlet conditions are a function solely of the inlet conditions. Examples of these include:

* unit operations involving equilibrium where the outlet condition can be calculated directly from the inlet condition.
* unit operations with negligible holdup where (total) flow out of the unit is always equal to the (total) flow in.

In general, a unit model is not written with a specific flowsheet or set of thermophysical property calculations in mind, thus it is necessary to provide a configuration argument (or arguments in cases where multiple streams interact) to allow the end-user to specify a property package to use with the model.  The example below shows a typical example of declaring a configuration argument for a single property package, along with a second argument that allows users to pass configuration arguments to the instances of the property packages when they are created.


.. code-block:: python

@declare_process_block_class("NewUnit")
class NewUnitDataData(UnitModelBlockData):

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations."""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s) and used when constructing these."""))

For unit models involving multiple property packages, or those that include reaction packages, additional pairs of configuration arguments are required for each of these. Model developers must provide unique names for each configuration argument, and are encourage to use meaningful names ot assist users in understanding what package should be linked to each argument.

Model developers may also declare additional configuration arguments to give users of their model the ability to change the behavior of different parts of the model. For example, the core iDAES Unit Model Library makes use of these to  provide flexibility in the form of the balance equations. Use of additional configuration arguments is entirely optional.

Unit Model `build` Method
-------------------------

Use of Control Volumes
^^^^^^^^^^^^^^^^^^^^^^

Unit Model Initialization
-------------------------

Whilst the `UnitModelBlockData` class contains a pre-built `initialize` method, this method is relatively simple and is unlikely to work for more complex models. For these situations, model developers will need to write their own `initialize` methods as part of their new unit model.

To create a custom initialization routine, model developers must create an `initialize` method as part of their model, and provide a sequence of steps intended to build up a feasible solutions. Developing initialization routines is one of the hardest aspects of model development, and generally involves starting with a simplified form of the model and progressively adding complexity. Initialization routines generally make use of Pyomo’s tools for activating and deactivating constraints and often involve solving multiple sub-problems whilst building up an initial state.

The example below shows the general form used when declaring a new initialization method:

.. code-block:: python

    def initialize(blk, state_args=None, outlvl=idaeslog.NOTSET,
                   solver='ipopt', optarg={'tol': 1e-6}):

* blk – local name for the block to be initialized.
* state_args – initial guesses for the state variables. The form of this may vary depending on the number and type of inlets to the unit model.
* outlvl – optional argument to allow users to control the amount of diagnostic output from the initialization procedure. His requires the use of the IDAES logger tools to function.
* solver – allows the user to set a solver to use for initialization.
* optarg – dict of options to pass to the solver; used to adjust solver behavior.

Unit Model Report
-----------------

.. code-block:: python

    def _get_performance_contents(self, time_point=0):
        var_dict = {"Volume": self.volume[time_point]}
        if hasattr(self, "heat_duty"):
            var_dict["Heat Duty"] = self.heat_duty[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        return {"vars": var_dict}

Tutorials
---------
Tutorials demonstrating how to create custom unit models are found
:ref:`here<advanced_user_guide/learning_materials/unit_tutorials/index:Unit Model Tutorials>`.    

    
