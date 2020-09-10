Custom Reaction Packages
========================

.. contents:: :local:

.. warning:: This section is currently being developed

Chemical reactions are a fundamental part of most processes, and models for these come in a wide range of different forms. Much like thermophysical property packages, the ability for users to define custom reaction formulations is a key aspect of the IDAES modeling paradigm.

Reaction packages within IDAES share many similarities with thermophysical property packages, both in form and content. Rather than repeat much of that documentation here, users should start by reading the :ref:`thermophysical property package documentation<advanced_user_guide/custom_models/property_package_development:Custom Property Packages>`, as this document will focus on the content of the reaction package.

What Belongs in a Reaction Package?
-----------------------------------

Chemical reactions are fundamentally governed by the same laws of thermodynamics as thermophysical properties, thus the separation of these into thermophysical and reaction packages is somewhat arbitrary. In IDAES, this separation between thermophysical and reaction packages was based on their expected frequency of use. All unit operations require some set of property calculations (e.g. enthalpy) and these types of calculations were grouped in the thermophysical package, whereas only a small subset of unit operations have chemical reactions and these types of calculations were grouped in the reaction package. This separation benefits the user in that they only need to be concerned about reactions in the unit operations that require them.

.. important::

    For the context of IDAES, chemical reactions are defined as phenomena where one chemical species is converted into another. This includes both rate limited and equilibrium reactions.

    On the other hand, phase equilibrium phenomena (where a chemical species changes phase) are handled via the thermophysical property package. 

However, users should note that reaction properties are fundamentally linked to the thermophysical properties, and that a reaction package should only be used with the thermophysical property package they were developed with (in theory at least). Due to this, when a reactions package is added to a model it must be coupled to a thermophysical property package. The modeling framework performs some limited checks to ensure the two packages are compatible (e.g. same set of base units) and that each reaction packages is only used in conjunction with its coupled thermophysical property package in unit models.

Reaction Package Classes
------------------------

Like thermophysical property packages, reaction property packages consist of two related model components; the Reaction Parameter Block and the Reaction Block, which are analogous to the Physical Parameter Block and State Block components. Similarly, when creating a custom reaction package users need to declare three new classes; the Reaction Parameter Block Data class, the Reaction Block Data class and the Reaction Block Methods class.

Build-on-Demand Properties
--------------------------

IDAES reaction packages also support build-on-demand properties using the same approach as for thermophysical properties.

The Reaction Parameter Block
----------------------------

The first part of the reaction package is the `ReactionParameterBlock`, which defines the global parameters and components of the property package. This includes:

* a reference to the `ReactionBlock` class associated with the `ReactionParameterBlock`,
* a Pyomo Set listing names for all rate-based reactions,
* a dict defining stoichiometric coefficients for all rate-based reactions,
* a Pyomo Set listing names for all equilibrium reactions,
* a dict defining stoichiometric coefficients for all equilibrium reactions,
* the base units of measurement for the property packages,
* the reaction properties supported by the property package, and
* the parameters required to calculate the reaction properties.

The starting point for creating a new `ReactionParameterBlock` is shown in the example below. The model developer needs to declare a new class which inherits from the `ReactionParameterBlock` base class, decorated using the `declare_process_block_class` decorator.

.. code-block:: python

    @declare_process_block_class("NewReactionParameterBlock")
    class NewReactionParameterData(ReactionParameterBlock):

        def build(self):
            super().build()

        @classmethod
        def define_metadata(cls, obj):
            obj.add_properties({# properties}})
            obj.add_default_units({# units})

The `NewReactionParameterData` class needs to contain a `build` method, and may also include a configuration block and a `define_metadata` classmethod as shown above. These methods and their contents will be explained below.

Reaction Parameter Block Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `ReactionParameterBlock` configuration block must contain the following two arguments:

* "property_package" - this configuration argument contains a pointer to the associated thermophysical property package (via an instance of a `PhysicalParameterBlock`), and is used for validating the link between thermophysical and reaction properties (e.g. confirming that both packages use the same set of base units).
* "default_arguments" - this configuration argument allows users to specify a set of default configuration arguments that will be passed to all `ReactionBlocks` created from an instance of a parameter block.

The Reaction Parameter Block `build` Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `build` method in the `NewReactionParameterBlock` class is responsible for constructing the various modeling components that will be required by the associated `ReactionBlocks`. This includes the indexing sets which will be used to identify individual reactions and the stoichiometry of each of these. The `build` method is also responsible for setting up the underlying infrastructure of the property package and making a link to the associated `ReactionBlock` class so that the modeling framework can automate the construction and linking of these.

The first step in the `build` method is to call `super().build()` to trigger the construction of the underlying infrastructure using the base class’ `build` method.

Next, the user must declare an attribute named "_reaction_block_class" which is a pointer to the associated `ReactionBlock` class (creation of this will be discussed later). An example of this is shown below, where the associated Reaction Block class is named `NewReactionBlock`.

.. code-block:: python

    def build(self):

        super().build()
        self._reaction_block_class = NewReactionBlock

Next, the `build` method must create two indexing sets which provide names for the rate- and equilibrium-based reaction respectively. These indexing sets must be named `rate_reaction_idx` and `equilibrium_reaction_idx`. These indexing sets will be used by the unit models and control volumes when creating reaction terms in material balance equations.

.. code-block:: python

    self.rate_reaction_idx = Set(initialize=["rate_rxn_1", "rate_rxn_2"])
    self.equilibrium_reaction_idx = Set(initialize=["equil_rxn_1", "equil_rxn_2"])

.. note::

    Users only need to define indexing sets and stoichiometry dicts for the types of reaction which they wish to model. E.g. users do not need to declare `rate_reaction_idx` and `rate_reaction_stoichiometry` if there are no rate-based reactions in their system.

The `build` method also needs to create stoichiometry `dicts` for the rate- and equilibrium-based reactions present in the system. These `dicts` should be named "rate_reaction_stoichiometry" and "equilibrium_reaction_stoichiometry" and have keys with the form (reaction_index, phase, component) and values equal to the stoichiometric coefficient for the given reaction, phase and component. A positive stoichiometric coefficient indicates a product of the reaction (i.e. generation) whilst a negative coefficient indicates a reactant (i.e. consumption). An example for defining the stoichiometry for rate-based reactions is shown below.

.. code-block:: python

    self.rate_reaction_stoichiometry = {
        ("rate_rxn_1", "phase_1", "component_1"): -1,  # Component 1 in phase 1 is a reactant
        ("rate_rxn_1", "phase_2", "component_1"): 0,  # Reaction 1 does no occur is phase 2
        ("rate_rxn_1", "phase_1", "component_2"): 2,  # Component 2 in phase 1 is a product
        ("rate_rxn_1", "phase_2", "component_2"): 0,
        ("rate_rxn_2", "phase_1", "component_1"): 0,
        ("rate_rxn_2", "phase_2", "component_1"): -1,
        ("rate_rxn_2", "phase_1", "component_2"): 0,
        ("rate_rxn_2", "phase_2", "component_2"): -1}  # etc.

.. important::

    Stoichiometry `dicts` must contain a key for every reaction-phase-component combination, even if the stoichiometric coefficient is zero.

Finally, the `build` method needs to declare all the global parameters that will be used by the reaction calculations. Similar to thermophysical property parameters, users are encouraged to declare these as Pyomo `Vars` rather than `Params` to facilitate parameter estimation studies.

Defining Reaction Metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^

The last part of creating a new Reaction Parameter block is to define the metadata associated with it. The reactions metadata serves four purposes:

1. The default units metadata is used by the framework to automatically determine the units of measurement of the resulting property model, and automatically convert between different unit sets where appropriate.
2. The properties metadata is used to set up any build-on-demand properties,
3. The metadata is also used by the Data Management Framework to index the available property packages to create a searchable index for users.
4. The units metadata is compared to that of the associated thermophysical property package (when an instance of the Reaction Parameter Block is declared), and an exception is raised if they do not match.

Setting Default Units
"""""""""""""""""""""

As with thermophysical property packages, the most important part of defining the metadata for a property package is to set the default units of measurement for each of the 7 base quantities (time, length, mass, amount, temperature, current (optional) and luminous intensity (optional)). These units are used by the modeling framework to determine the units of measurement for all other quantities in the process that are related to this property package. More importantly, the units metadata is used to determine if a reaction package is comparable with a given thermophysical property package when they are declared – if the units metadata does not match, an exception will be raised and the two packages cannot be used together. 

Units must be defined using Pyomo `Units` components, as shown in the example below:

.. code-block:: python

    from pyomo.environ import units

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': units.s,
                               'length': units.m,
                               'mass': units.kg,
                               'amount': units.mol,
                               'temperature': units.K})

Setting Reaction Metadata
"""""""""""""""""""""""""

Similar to thermophysical property packages, reaction packages allow users to specify the set of reaction properties supported by a given reaction package. This is also used to set up the build-on-demand properties system in the same way as thermophysical properties. For more information, see the documentation for :ref:`thermophysical properties metadata<advanced_user_guide/custom_models/property_package_development:Setting Properties Metadata>`. 

The Reaction Block
------------------

The second part of a reaction property package is the `ReactionBlock` class. Similarly to `StateBlock` classes this is defined using two user-written classes; the `ReactionBlockData` class and the `ReactionBlockMethods` class.

.. code-block:: python

    @declare_process_block_class("NewReactionBlock",
                                 block_class=NewReactionBlockMethods)
    class NewReactionBlockData(ReactionBlockData):

        def build(self):
            super().build()

The Reaction Block Data Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One important difference between Reaction Blocks and State Blocks is that while State Blocks are fully self-contained and can be solved in isolation, Reaction Blocks depend upon the State Block for the definition of the state variables. This means that Reaction Blocks do not need to redefine the state variables (which are needed for the reaction properties), but at the cost of not being independent, self-contained models. This is one of the reasons why reaction packages are so closely tied to thermophysical property packages within IDAES.

The purpose of the Reaction Block Data class is to define the reaction properties that will be required by the unit models using this package. The three main properties required for material and energy balances are:

* rate terms for rate-based reactions,
* equilibrium constraints for equilibrium-based reactions, and
* heats of reaction (if required, see note below).

These properties may in turn depend on other reaction properties such as equilibrium and rate constants. All of these properties may be constructed using the build-on-demand framework.

All reaction properties depend upon the state of the material, which is defined in the State Block; thus it is necessary to reference the associated State Block whenever these are needed. In order to facilitate this, the ReactionBlockData base class establishes a reference to the associated State Block which users can use to obtain state variables and properties from the State Block. For example, temperature can be referenced from the state block as shown below:

.. code-block:: python

    temperature = self.state_ref.temperature

.. note::

    There are multiple ways in which heat of reaction may be included in a model, and users should consider which is most suitable for their application. The two most common approaches are to include an explicit heat of reaction term in the energy balance equations, or to incorporate heat of reaction into the specific enthalpy terms (generally via heats of formation). The IDAES Process Modeling Framework supports both of these approaches.

Reaction Block Data Configuration Arguments
"""""""""""""""""""""""""""""""""""""""""""

The ReactionBlockData base class defines three configuration arguments that are required for all Reaction Block Data classes.

* "parameters" – this argument is used to provide a link back to the associated `ReactionParameterBlock`, and is generally automatically passed to the `ReactionBlock` when it is constructed.
* "state_block" – this argument is used to provide a link to the State Block associated with this Reaction Block, as is generally passed to the `ReactionBlock` by the unit model when it is constructed. This argument is used to the `state_ref` attribute shown above for referencing properties from the State Block.
* "has_equilibrium" – this argument indicates whether equilibrium reaction will be considered for this state. In most cases, this argument will always be True, however this allows users the ability to turn off equilibrium reactions if they desire.

The Reaction Block `build` Method
"""""""""""""""""""""""""""""""""

As with all IDAES components, the `build` method forms the core of a `ReactionBlockData` class, and contains the instructions on how to construct the variables, expressions and constraints required by the reaction model. As usual, the first step in the `build` method should be to call `super().build()` to trigger the construction of the underlying components required for Reaction Blocks to function.

Variables and Properties
""""""""""""""""""""""""

The same set of guidelines for defining thermophysical properties apply to reaction properties, :ref:`which can be found here<advanced_user_guide/custom_models/property_package_development:State Variables and Properties>`. 

Required Methods
""""""""""""""""

In addition to the `build` method, Reaction Blocks require one additional method which is used to define the basis for the reaction terms.

* `get_reaction_rate_basis` - must return a `MaterialFlowBasis` `Enum`, and is used to automatically convert reaction terms between mass and mole basis in control volumes.

The Reaction Block Methods Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Reaction Block Methods class is very similar to the :ref:`State Block Methods class<advanced_user_guide/custom_models/property_package_development:The State Block Methods class>`. The Reaction Block Methods class needs to contain an `initialize` method (however a `release_state` method is not required as Reaction Blocks do not contain state variables).

The `initialize` Method
"""""""""""""""""""""""

Initialization of Reaction Blocks is complicated by the fact that they depend upon the State Block for the state variables, and thus cannot be solved as a stand-alone model. Within the wider IDAES modeling framework, this is handled by initializing the Reaction Block after the State Block `initialization` method has been called (and thus all state variables and properties are initialized) but before the `release_State` method is called (thus all state variables are fixed). Thus, the Reaction Block can assume that the state is fully defined and initialized (although it may not be possible to use a solver as part of the Reaction Block’s initialization procedure).

However, Reaction Blocks also tend to be much simpler than State Blocks, involving fewer properties which are generally much less tightly coupled (most reaction properties are functions solely of the state variables), which simplifies the requirements of initializing the sub-model.

Tutorials
---------

Tutorials demonstrating how to create custom reaction packages are being developed. Once they are created, they will be found :ref:`here<tutorials_examples:Tutorials and Examples>`.

