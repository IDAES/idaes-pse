Custom Unit Models
==================

.. contents:: :local:

UnitModelBlockData
------------------

The starting point for all unit models within the IDAES Process Modeling Framework is the `UnitModelBlockData` base class. This class contains a number of methods to assist users with creating new unit models including:

* checking and validating the “dynamic” and “has_holdup” configuration arguments to ensure consistency,
* adding `Port` objects to the model,
* creating simple material balances between states (equal flow of each component in each phase), and,
* a method for initializing simple unit models.

More details on the `UnitModelBlockData` class can be found in the :ref:`technical specifications<reference_guides/core/unit_model_block:Unit Model Class>`.

Unit Model Configuration
------------------------

Configuration arguments in unit models allow model developers to provide the end-user with the ability to configure the model to suit the needs of the flowsheet they are simulating. The most common aspects that need to be configured are:

* whether the model should be dynamic or steady-state, and
* what property package to use when calculating thermophysical properties.

The `UnitModelBlockData` class contains a simple configuration block which includes two configuration arguments; “dynamic” and “has_holdup”. These arguments are required for any model which is expected to be used in both steady-state and dynamic flowsheets and are used to determine whether accumulation and holdup terms should be constructed and included in the material balance equations. There are some situations, however, where a model is inherently steady-state (even if it is included in a dynamic flowsheet), notably those where outlet conditions are a function solely of the inlet conditions. Examples of these include:

* unit operations involving equilibrium where the outlet condition can be calculated directly from the inlet condition.
* unit operations with negligible holdup where (total) flow out of the unit is always equal to the (total) flow in.

In general, a unit model is not written with a specific flowsheet or set of thermophysical property calculations in mind, thus it is necessary to provide a configuration argument (or arguments in cases where multiple streams interact) to allow the end-user to specify a property package to use with the model. The example below shows how to declare a configuration argument for a single property package, along with a second argument that allows users to pass configuration arguments to the instances of the property packages when they are created.

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

For unit models involving multiple property packages, or those that include reaction packages, additional pairs of configuration arguments are required for each of these. Model developers must provide unique names for each configuration argument, and are encourage to use meaningful names to assist end-users in understanding what package should be linked to each argument.

Model developers may also declare additional configuration arguments to give end-users the ability to change the behavior of different parts of the model. For example, the core IDAES Unit Model Library makes use of these to  provide flexibility in the form of the balance equations. Use of additional configuration arguments is entirely optional.

Unit Model `build` Method
-------------------------

The `build` method for a unit model must include all the instructions necessary for constructing the representation of the unit operation. This generally involves the following steps:

1. Calling `super().build()` to trigger the behind-the-scenes code.
2. Adding any variables and constraints required to describe the system geometry.
3. Adding `State Blocks` to the model to represent each of the material states in the system.
4. Adding the necessary material balances and associated variable to describe the flow of material between each state.
5. Adding the necessary energy balances and associated variable to describe the flow of energy between each state.
6. Adding the necessary momentum balances and associated variable to describe the flow of momentum between each state.
7. Adding any additional performance equations and associated variable that govern the behavior of the unit operation.
8. Adding the required inlet and outlet Ports to allow the unit model to be included in a flowsheet.

For some applications, not all of these steps will be required (e.g. a process in which pressure drop is negligible may be able to skip adding momentum balances).

The above steps represent a significant amount of work, and in many cases require a detailed understanding of how the IDAES framework is structured. To reduce the effort and knowledge required to create new models, the framework provides a number of tools to automate these steps for common cases. Users are encouraged to familiarize themselves with the methods available in :ref:`UnitModelBlockData<reference_guides/core/unit_model_block:Unit Model Class>` and the use of control volumes.

Control Volumes
^^^^^^^^^^^^^^^

The IDAES Process Modeling Framework includes tools to assist users with creating new models in the form of the Control Volume libraries. These libraries contain methods for performing the common task associated with building unit models, such as creating material, energy and momentum balances. Users are free to choose whether or not to use these libraries, but are encouraged to understand what is available in these as they can greatly reduce the amount of effort required by the user.

The IDAES Process Modeling Framework currently includes two types of Control Volumes:

1. :ref:`ControlVolume0D<reference_guides/core/control_volume_0d:0D Control Volume Class>` for inlet-outlet type models where spatial variation are not significant.
2. :ref:`ControlVolume1D<reference_guides/core/control_volume_1d:1D Control Volume Class>` for models where spatial variation in one-dimension are required.

Unit Model Initialization
-------------------------

Whilst the `UnitModelBlockData` class contains a pre-built `initialize` method, this method is relatively simple and is unlikely to work for more complex models. For these situations, model developers will need to write their own `initialize` methods as part of their new unit model.

To create a custom initialization routine, model developers must create an `initialize` method as part of their model, and provide a sequence of steps intended to build up a feasible solution. Developing initialization routines is one of the hardest aspects of model development, and generally involves starting with a simplified form of the model and progressively adding complexity. Initialization routines generally make use of Pyomo’s tools for activating and deactivating constraints and often involve solving multiple sub-problems whilst building up an initial state.

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

Users are likely already aware of the `report` method which is available in all IDAES models and prints a summary of the current state of a given model. This functionality is also part of `UnitModelBlockData` and is thus included in all custom unit models, however model developers need to define what information should be included in the output.

The `report` method will automatically search for and identify all `Ports` in the model to be included in the summary stream table, however model developers must identify any performance variables they wish to include in the summary. This is done by declaring a `_get_performance_contents` method as shown in the example below:

.. code-block:: python

    def _get_performance_contents(self, time_point=0):
        var_dict = {"display name": self.var[time_point]}
        expr_dict = {"display name": self.expr[time_point]}
        param_dict = {"display name": self.param[time_point]}

        return {"vars": var_dict, "exprs": expr_dict, "params": param_dict}

The `_get_performance_contents` method should take two arguments, the first being the model object and the second being a time point at which to report the model state. The method should return a dictionary-of-dictionaries with one to three keys; "vars", "exprs" and "params". The entries from these will be included in the model summary under the headings of Variables, Expressions and Parameters respectively.

Tutorials
---------
Tutorials demonstrating how to create custom unit models are found
:ref:`here<tutorials/tutorials_examples:Tutorials and Examples>`.


