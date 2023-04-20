Physical Property Package Classes
=================================

.. contents:: Contents
    :depth: 2

Physical property packages represent a collection of calculations necessary to determine the state properties of a given material. Property calculations form a critical part of any process model, and thus property packages form the core of the IDAES modeling framework.

Physical property packages consist of two parts:

* PhysicalParameterBlocks, which contain a set of parameters associated with the specific material(s) being modeled, and
* StateBlocks, which contain the actual calculations of the state variables and functions.

Physical Parameter Blocks
-------------------------

Physical Parameter blocks serve as a central location for linking to a property package, and contain all the parameters and indexing sets used by a given property package.

PhysicalParameterBlock Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The role of the PhysicalParameterBlock class is to set up the references required by the rest of the IDAES framework for constructing instances of StateBlocks and attaching these to the PhysicalParameter block for ease of use. This allows other models to be pointed to the PhysicalParameter block in order to collect the necessary information and to construct the necessary StateBlocks without the need for the user to do this manually.

Physical property packages form the core of any process model in the IDAES modeling framework, and are used by all of the other modeling components to inform them of what needs to be constructed. In order to do this, the IDAES modeling framework looks for a number of attributes in the PhysicalParameter block which are used to inform the construction of other components.

* state_block_class - a pointer to the associated class that should be called when constructing StateBlocks. This should only be set by the property package developer.
* phase_list - a Pyomo Set object defining the valid phases of the mixture of interest.
* component_list - a Pyomo Set defining the names of the chemical species present in the mixture.
* element_list - (optional) a Pyomo Set defining the names of the chemical elements that make up the species within the mixture. This is used when doing elemental material balances.
* element_comp - (optional) a dict-like object which defines the elemental composition of each species in component_list. Form: component: {element_1: value, element_2: value, ...}.
* supported properties metadata - a list of supported physical properties that the property package supports, along with instruction to the framework on how to construct the associated variables and constraints, and the units of measurement used for the property. This information is set using the add_properties attribute of the define_metadata class method.

Physical Parameter Configuration Arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Physical Parameter blocks have one standard configuration argument:

* default_arguments - this allows the user to provide a set of default values for construction arguments in associated StateBlocks, which will be passed to all StateBlocks when they are constructed.


.. module:: idaes.core.base.property_base

.. autoclass:: PhysicalParameterBlock
    :members:

State Blocks
------------

State Blocks are used within all IDAES Unit models (generally within ControlVolume Blocks) in order to calculate physical properties given the state of the material. State Blocks are notably different to other types of Blocks within IDAES as they are always indexed by time (and possibly space as well). There are two base Classes associated with State Blocks:

* StateBlockData forms the base class for all StateBlockData objects, which contain the instructions on how to construct each instance of a State Block.
* StateBlock is used for building classes which contain methods to be applied to sets of Indexed State Blocks (or to a subset of these). See the documentation on declare_process_block_class and the IDAES tutorials and examples for more information.

State Block Construction Arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

State Blocks have the following construction arguments:

* parameters - a reference to the associated Physical Parameter block which will be used to make references to all necessary parameters.
* defined_state - this argument indicates whether the State Block should expect the material state to be fully defined by another part of the flowsheet (such as by an upstream unit operation). This argument is used to determine whether constraints such as sums of mole fractions should be enforced.
* has_phase_equilibrium - indicates whether the associated Control Volume or Unit model expects phase equilibrium to be enforced (if applicable).

Constructing State Blocks
^^^^^^^^^^^^^^^^^^^^^^^^^

State Blocks can be constructed directly from the associated Physical Parameter Block by calling the `build_state_block()` method on the Physical Parameter Block. The `parameters` construction argument will be automatically set, and any other arguments (including indexing sets) may be provided to the `build_state_block` method as usual.

StateBlockData Class
^^^^^^^^^^^^^^^^^^^^

StateBlockData contains the code necessary for implementing the as needed construction of variables and constraints.


.. autoclass:: StateBlockData
    :members:

StateBlock Class
^^^^^^^^^^^^^^^^

.. autoclass:: StateBlock
    :members:
