Custom Property Packages
========================

.. contents:: :local:

.. warning:: This section is currently being developed

Physical property packages form the core of all IDAES process models, and the ability for users to develop their own property formulations is a key aspect of the IDAES modeling paradigm. In order to support the flexibility of the IDAES Process Modeling Framework, property packages define a number of key aspects that inform the eventual structure of the final process model. This however places the burden of making these decisions on the developer of each property package, which are implemented as part of the property package classes.

Property Package Classes
------------------------

As users of the IDAES Process Modeling Framework, you are likely already aware that property packages (of both types) consist of two related model components; in this case the Physical Parameter Block and the State Block. However, when creating a new thermophysical property package, developers need to define three (rather than two) new classes.

The first of these classes is a `PhysicalParameterBlock` class, which is responsible for constructing the Physical Parameter Block. However, two classes are required for defining the State Block; a `StateBlockMethods` class and a `StateBlockData` class. The reason for this is because State Blocks are always indexed (by time and occasionally by space) . The `StateBlockData` class represents an individual state at a point in space and time (i.e. one element of the indexed `StateBlock`), and as such contains a set of state variables and the constraints necessary for calculating the desired thermophysical properties at that state. However, we often want to perform actions on the entire set of states (i.e. `StateBlockDatas`) in one go, such as during initialization. Whilst we could initialize each state individually, as the process for each state is generally identical except for values, it is much more efficient to perform the same set of instructions on all states simultaneously. The StateBlockMethods serves this purpose by defining the methods that should be applied to multiple StateBlockDatas simultaneously. When a model requires a State Block, these two classes are combined to produce the final model. The distinction and use of the `StateBlockMethods` and `StateBlockData` classes will be discussed further later in this documentation.

Build-on-Demand Properties
--------------------------

Before moving onto a discussion of the contents of each of the three classes, it is important to introduce the concept of build-on-demand properties. Property packages generally tend to include methods for a number of properties, but not all of these will be required by every unit model. In order to reduce model complexity and avoid calculating properties which are not required for a given unit operation, the IDAES framework supports the concept of build-on-demand properties, where the variables and constraints related to a given property are only constructed if called for in a given state.

It must be noted that this is an advanced feature and is entirely optional. Whilst it can reduce the complexity of individual models, it also increases the complexity of the model instructions and can increase the chance of errors during model constructions. Property package developers should decide up front if they wish to implement build-on-demand properties for their property packages, and which properties this will be implemented for (i.e. it is possible to use the build-on-demand infrastructure) for a subset of the properties within a package.

The Physical Parameter Block
----------------------------

The first part of the physical property package is the `PhysicalParameterBlock`, which defines the global parameters and components of the property package. This includes:

* a reference to the `StateBlock` class associated with the `PhysicalParameterBlock`,
* the chemical species or components in the material,
* the thermodynamic phases of interest,
* the base units of measurement for the property packages,
* the thermophysical properties supported by the property package, and
* the parameters required to calculate the thermophysical properties.

The starting point for creating a new `PhysicalParameterBlock` is shown in the example below. The model developer needs to declare a new class which inherits from the `PhysicalParameterBlock` base class, decorated using the `declare_process_block_class` decorator.

.. code-block:: python

    @declare_process_block_class("NewPhysicalParameterBlock")
    class NewPhysicalParameterData(PhysicalParameterBlock):

        def build(self):
            super().build()

        @classmethod
        def define_metadata(cls, obj):
            obj.add_properties({# properties}})
            obj.add_default_units({# units})

The `NewPhysicalParameterData` class needs to contain a `build` method, and may also include a configuration block and a `define_metadata` classmethod as shown above. These methods and their contents will be explained below.

Physical Parameter Block Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Like all IDAES models, Physical Parameter Blocks can have configuration arguments which can be used to adjust the form of the resulting model. The default configuration block which comes from the `PhysicalParameterBlock` base class contains a single configuration argument:

* “default_arguments” - this configuration argument allows users to specify a set of default configuration arguments that will be passed to all `StateBlocks` created from an instance of a parameter block.

The Physical Parameter Block `build` Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `build` method in the `NewPhysicalParameterBlock` class is responsible for constructing the various modeling components that will be required by the associated `StateBlocks`, such as the sets components and phases that make up the material, and the various parameters required by the property calculations. The `build` method is also responsible for setting up the underlying infrastructure of the property package and making a link to the associated `StateBlock` class so that the modeling framework can automate the construction and linking of these.

The first step in the `build` method is to call `super().build()` to trigger the construction of the underlying infrastructure using the base class' `build` method.

Next, the user must declare an attribute named “_state_block_class” which is a pointer to the associated `StateBlock` class (creation of this will be discussed later). An example of this is shown below, where the associated State Block is named `NewStateBlock`.

.. code-block:: python

    def build(self):

        super().build()
        self._state_block_class = NewStateBlock

The next step in the `build` method is to define the chemical species and phases necessary to describe the material of interest. This is done by adding :ref:`Component<reference_guides/core/comp:Component Class>` and :ref:`Phase<reference_guides/core/phase:Phase Class>` objects, as shown below.

.. code-block:: python

    def build(self):
        self.benzene = Component()
        self.toluene = Component()

        self.liquid = LiquidPhase()
        self.vapor = VaporPhase()

.. note::

    The IDAES Process Modeling Framework supports a number of different types of `Component` and `Phases` objects, as discussed in the associated documentation. Users should use the type most appropriate for their applications. Also note that whilst `Component` and `Phase` objects contain configuration arguments, these are primarily for use by the Generic Property Package framework, and are not required for custom property packages.

Finally, the `build` method needs to declare all the global parameters that will be used by the property calculations. By declaring these in a single central location rather than in each State Block, this reduces the number of parameters present in the model (thus reducing memory requirements) and also facilitates parameter estimation studies using these parameters.

.. note::

    Whilst we generally use the term “parameters” to describe these global coefficients used in property correlations, it is often better to declare these as Pyomo `Var` objects with fixed values (rather than as `Param` objects). The reason for this is because, despite the name, it is not possible to estimate the value of `Params` using parameter estimation tools (as their value is concrete and cannot be changed).

Defining Property Metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^

The last part of creating a new Physical Parameter block is to define the metadata associated with it. The properties metadata serves three purposes:

1. The default units metadata is used by the framework to automatically determine the units of measurement of the resulting property model, and automatically convert between different unit sets where appropriate.
2. The properties metadata is used to set up any build-on-demand properties,
3. The metadata is also used by the Data Management Framework to index the available property packages to create a searchable index for users.

Setting Default Units
"""""""""""""""""""""

The most important part of defining the metadata for a property package is to set the default units of measurement for each of the 7 base quantities (time, length, mass, amount, temperature, current (optional) and luminous intensity (optional)). These units are used by the modeling framework to determine the units of measurement for all other quantities in the process that are related to this property package. Units must be defined using Pyomo `Units` components, as shown in the example below:

.. code-block:: python

    from pyomo.environ import units

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': units.s,
                               'length': units.m,
                               'mass': units.kg,
                               'amount': units.mol,
                               'temperature': units.K})

Setting Properties Metadata
"""""""""""""""""""""""""""

The primary purpose of the properties metadata is to set up the build-on-demand system used to selectively construct only those properties required by a given unit operation. In order to do this, the user needs to add each property they wish to build-on-demand along with the name of a method that will be called whenever the property is required (this method will be created later as part of the `StateBlockData` class). Users are also encouraged to list *all* properties supported by their property packages here, setting `None` as the method associated with the property for those which are always constructed. An example for both uses is shown below:

.. code-block:: python

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'property_1': {'method': method_name},  # a build-on-demand property
                'property_2': {'method': None}})  # a property that will always be constructed

.. note::

  The name of a property in the metadata dictionary must match the name of the property component (normally a variable) that will be called for. These names should be drawn form the :ref:`standard naming conventions<explanations/conventions:Standard Variable Names>`.

The State Block
---------------

The second part of a thermophysical property package is the `StateBlock` class, which as mentioned earlier is defined using two user-written classes; the `StateBlockData` class and the `StateBlockMethods` class. Declaration of the `StateBlock` class is similar to that of other modeling classes, but makes use of a special aspect of the `declare_process_block_class` decorator as shown in the example below.

.. code-block:: python

    @declare_process_block_class("NewStateBlock",
                                 block_class=NewStateBlockMethods)
    class NewStateBlockData(StateBlockData):

        def build(self):
            super().build()

As can be seen, the declaration of the new `StateBlock` class (`NewStateBlock`) looks similar to that of other modeling class declarations, where the `declare_process_block_class` is applied to a user defined `NewStateBlockData` class. However, in this case we also provide an additional argument to the decorator; the "block_class" argument allows us to attach a set of methods declared in a user-defined class (in this case `NewStateBlockMethods`) to the `NewStateBlock` class, which can be applied across all members of an indexed `NewStateBlock` (methods in the `NewStateBlockData` class can only be applied to a single indexed element).

The State Block Data class
^^^^^^^^^^^^^^^^^^^^^^^^^^

As part of the core of the IDAES Process Modeling Framework, the `StateBlockData` class is responsible not only for defining the variables, expressions and constraints which describe the thermophysical properties of the material in question, but also providing information to the rest of the Process Modeling Framework on how the higher levels models should be formulated. As such, `StateBlockData` classes need to define more methods than any other component class. The base class for developing new `StateBlockData` classes is `StateBlockData`, which includes a configuration block with a number of critical configuration arguments as well as the code necessary for supporting “build-on-demand properties”.

State Block Data Configuration Arguments
""""""""""""""""""""""""""""""""""""""""

The `StateBlockData` base class configuration contains three configuration arguments that are expected by the modeling framework and must be included in and user defined `StateBlockData`. These configuration arguments are:

* "parameters" – this argument is used to provide a link back to the associated `PhysicalParameterBlock`, and is generally automatically passed to the `StateBlock` when it is constructed.
* "defined_state" – this argument is used to indicate whether this state represents a point in the process where all state variables are defined. The most common case for this is for inlets to unit models, where all inlets states are known from the outlet of the previous unit model. In these cases, it is not possible to write certain constraints, such as the sum of mole fractions, without over specifying the system of equations; this argument identifies these cases so that generation of these constraints can be automatically skipped.
* "has_phase_equilibrium" – this argument indicates whether phase equilibrium will be considered for this state. Phase equilibrium constraints decrease the degrees of freedom in the system thus it is important to determine when and where these constraints should be written. Note that equilibrium constraints can never be written for cases where the state is fully defined (as above), thus both this and the `defined_state` arguments must be considered when determining whether to include equilibrium constraints.

The State Block `build` Method
""""""""""""""""""""""""""""""

As with all IDAES components, the `build` method forms the core of a `StateBlockData` class, and contains the instructions on how to construct the variables, expressions and constraints required by the thermophysical model. As usual, the first step in the `build` method should be to call `super().build()` to trigger the construction of the underlying components required for State Blocks to function.

State Variables and Properties
""""""""""""""""""""""""""""""

The most important part of the construction of a State Block is defining the necessary set of variables, expression and constraints that make up the property model. There are many different ways in which these can be defined and formulated, and there is no single “best” way to do this; different approaches may work better for different applications. However, there are some general rules that should be followed when defining the variables which make up a State Block.

1. All state variables and properties should use the IDAES naming conventions. Standard names allow linking between different types of models to be automated, as no cross-referencing of names is required.
2. All properties within a property package should use a consistent set of base units. This is most easily accomplished by selecting a set of units for the 7 base SI quantities (time, length, mass, amount, temperature, current and luminous intensity) and deriving units for all quantities from these. Modelers should also select units based solely on convenience or ease of use – scaling of variables and equations is better handled separately using the :ref:`IDAES scaling tools<reference_guides/core/util/scaling:Scaling Methods>`.

Beyond these requirements, modelers are free to choose the form of their model to best suit theirs needs and make the most tractable problem possible. Modelers are also free to combine variable and constraints with expression for some quantities as needed. The IDAES Process Modeling Framework is concerned only that the expected quantities are present (i.e. the expected variable/expression names), not their exact form or how they are calculated.

As described throughout this page, IDAES supports "build-on-demand" for property correlations. Details on how to define methods for building properties on demand is demonstrated in the tutorials (see link at bottom of page).

Required Methods
""""""""""""""""

As the foundation of the entire Process Modeling Framework, the definition of a new `StateBlockData` class needs to include a number of methods that the framework relies on for determining the formulation of the higher level models.

Below is a list of the required methods, along with a short description.

* `get_material_flow_basis(block)` – this method is used to define the basis on which material balance terms will be expressed. This is used by the framework to automatically convert between mass and mole basis if required, and the method needs to return a `MaterialFlowBasis` `Enum`.
* `get_material_flow_terms(block, phase, component)` – this method is used to determine the form of the material flow terms that are constructed as part of the material balance equations in each unit model. This method needs to take three arguments; a reference to the current state block, a phase name and a component name, and must return an expression for the material flow term for the given phase and component.
* `get_material_density_terms(block, phase, component)` – similar to the `get_material_flow_terms` method, this method is used to determine the form of the density term which should be used when constructing material holdup terms in the material balances. This method also needs to take three arguments; a reference to the current state block, a phase name and a component name, and must return an expression for the material density term for the given phase and component.
* `get_material_diffusion_terms(block, phase, component)` – Support for this is not currently implemented.
* `get_enthalpy_flow_terms(block, phase)` – this method is used to determine the form of the enthalpy flow terms that are constructed as part of the energy balance equations in each unit model.  This method needs to take two arguments; a reference to the current state block and a phase name, and must return an expression for the enthalpy flow term for the given phase and component.
* `get_energy_density_terms(block, phase)` – this method is used to determine the form of the energy density terms that are required for the holdup terms in the energy balance equations. This method needs to take two arguments; a reference to the current state block and a phase name, and must return an expression for the energy density term for the given phase and component. Note that the holdup/density term needs to be in terms of internal energy, not enthalpy.
* `get_energy_diffusion_terms(block, phase)` – Support for this is not currently implemented.
* `default_material_balance_type(block)` – this method is used to set a default for the type of material balance to be written by a Control Volume if the user does not specify which type to use. This method needs to return a `MaterialBalanceType` `Enum`.
* `default_energy_balance_type(block)` – this method is used to set a default for the type of energy balance to be written by a Control Volume if the user does not specify which type to use. This method needs to return a `EnergyBalanceType` `Enum`.
* `define_state_vars(block)` – this method is used to define the set of state variables which should be considered the state variables for the property package, and is used in a number of methods associated with model initialization to determine which variables should be fixed. This method must return a Python dict, where the keys are the variable name as a string, and the values are the variables.
* `define_port_members(block)` – similar to the `define_state_vars` method, this method is used to define what variables should be part of the inlet/outlet ports of a unit model. In many cases, these variables are equivalent to the state variables of the property package and if so this method can be skipped (if undefined `define_state_vars` is called instead). This method is similar to the one in the above method, however in this case the key names can be defined by the user for improved readability (instead of having to be the variable name).
* `define_display_vars(block)` – similar again to the `define_state_vars` method, this method is used to define a set of variables which should be used when generating the output of the `report` method for this property package. Again, this is often the same as the state variables, but allows modelers to include additional variables beyond just the state variables (or port members). Similarly to the `define_port_members` method, this method can be skipped (in which case it defaults to `define_state_vars`) and the key names in the dict can be defined by the user.

The State Block Methods class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The purpose of the `StateBlockMethods` class is to define methods which can be applied to to an entire set of indexed `StateBlocks` simultaneously. Whilst the `StateBlockData` class contain the instructions for how to build the variables and constraints that describe the state of a material at a single point in space and time, the `StateBlockMethods` class defined methods for interacting with multiple states across space and time simultaneously. The most common application for this is during initialization of `StateBlocks`, where the same set of instructions needs to to be performed on each indexed state; whilst this could be done by iterating over each state and performing the set of instructions, it is generally more efficient to apply the instructions simultaneously across all states.

Declaration and Base Class
""""""""""""""""""""""""""

Due to the way the `StateBlockMethods` class is provided to the `declare_process_block_class` decorator on the `NewStateBlockData` class, this is one of the few cases where the decorator is not required when declaring a class within IDAES. An example of declaring a new `StateBlockMethods` class is shown below, using the `StateBlock` base class:

.. code-block:: python

    class NewStateBlockMethods(StateBlock):

As the `StateBlockMethods` class is designed to contain methods that can be applied to multiple existing `StateBlockData` object, rather than construct a model itself, the `StateBlockMethods` class does not need a `build` method either, nor is it necessary to call `super().build()` as is normal for other modeling components.

Instead, the `StateBlockMethods` class should contain a set of methods which can be called and applied to an indexed `StateBlock` as required. The two methods that must be declared are:

* `initialize`
* `release_state`

The `initialize` and `release_state` Methods
""""""""""""""""""""""""""""""""""""""""""""

When initializing a unit model, most IDAES models use a hierarchical approach where each state in the model (i.e. each `StateBlockData`) is first initialized at some initial state, after which the unit model attempts to build up and solve the material, energy and momentum balances, etc. The purpose of the `initialize` method is to provide a set of instructions which can take a state from its initial state to a solvable final state at the set of initial conditions (provided as arguments to the `initialize` method). This is generally done by:

1. fixing the state variables at the initial conditions,
2. performing a series of steps to build up the final solution,
3. solving the full state model, and
4. unfixing the state variables (unless they were already fixed when the process began).

However, in order to fully initialize the unit operation (which contains these material state) it is necessary for the unit model to be fully defined (with zero degrees of freedom, i.e. a square model). In order for this to be true however, it is necessary for the inlet states to remain fixed until the unit model has finished initializing. This requires step 4 above to be postponed for inlet states until the unit model has finished initializing, thus the above process is broken into two methods.

1. The `initialize` method covers steps 1-3 above, and is called at the beginning of the unit model initialization process.
2. The `release_state` method covers step 4; for inlet states this is called when the unit model has finished initialization, whilst for all other states it is called immediately by the `initialize` methods when it finishes.

More details on writing initialization methods will be provided elsewhere in the documentation of tutorials.

Tutorials
---------

Tutorials demonstrating how to create custom property packages are being developed. Once they are created, they will be found :ref:`here<tutorials/tutorials_examples:Tutorials and Examples>`.
