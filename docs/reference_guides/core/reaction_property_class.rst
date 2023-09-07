Reaction Property Package Classes
=================================

.. contents:: Contents
    :depth: 2

Reaction property packages represent a collection of calculations necessary to determine the reaction behavior of a mixture at a given state. Reaction properties depend upon the state and physical properties of the material, and thus must be linked to a StateBlock which provides the necessary state and physical property information.

Reaction property packages consist of two parts:

* ReactionParameterBlocks, which contain a set of parameters associated with the specific reaction(s) being modeled, and
* ReactionBlocks, which contain the actual calculations of the reaction behavior.

Consistency with Thermophysical Properties
------------------------------------------

Within the IDAES modeling framework, all reaction packages are coupled with a thermophysical property package. The thermophysical property package contains the state variables required for calculating reaction properties, and in some cases may also provide thermophysical properties required by reaction calculations. Due to this, reaction packages must be consistent with the thermophysical property package they are linked to and the modeling framework performs some checks to ensure this. Notably, the default units of measurement defined for the reaction package and the thermophysical property package must match.

Reaction Parameter Blocks
-------------------------

Reaction Parameter blocks serve as a central location for linking to a reaction property package, and contain all the parameters and indexing sets used by a given reaction package.

ReactionParameterBlock Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The role of the ReactionParameterBlock class is to set up the references required by the rest of the IDAES framework for constructing instances of ReactionBlocks and attaching these to the ReactionParameter block for ease of use. This allows other models to be pointed to the ReactionParameter block in order to collect the necessary information and to construct the necessary ReactionBlocks without the need for the user to do this manually.

Reaction property packages are used by all of the other modeling components to inform them of what needs to be constructed when dealing with chemical reactions. In order to do this, the IDAES modeling framework looks for a number of attributes in the ReactionParameter block which are used to inform the construction of other components.

* reaction_block_class - a pointer to the associated class that should be called when constructing ReactionBlocks. This should only be set by the property package developer.
* phase_list - a Pyomo Set object defining the valid phases of the mixture of interest.
* component_list - a Pyomo Set defining the names of the chemical species present in the mixture.
* rate_reaction_idx - a Pyomo Set defining a list of names for the kinetically controlled reactions of interest.
* rate_reaction_stoichiometry - a dict-like object defining the stoichiometry of the kinetically controlled reactions. Keys should be tuples of (rate_reaction_idx, phase_list, component_list) and values equal to the stoichiometric coefficient for that index.
* equilibrium_reaction_idx - a Pyomo Set defining a list of names for the equilibrium controlled reactions of interest.
* equilibrium_reaction_stoichiometry - a dict-like object defining the stoichiometry of the equilibrium controlled reactions. Keys should be tuples of (equilibrium_reaction_idx, phase_list, component_list) and values equal to the stoichiometric coefficient for that index.
* supported properties metadata - a list of supported reaction properties that the property package supports, along with instruction to the framework on how to construct the associated variables and constraints, and the units of measurement used for the property. This information is set using the add_properties attribute of the define_metadata class method.
* required properties metadata - a list of physical properties that the reaction property calculations depend upon, and must be supported by the associated StateBlock. This information is set using the add_required_properties attribute of the define_metadata class method.

Reaction Parameter Configuration Arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reaction Parameter blocks have two standard configuration arguments:

* property_package - a pointer to a PhysicalParameterBlock which will be used to construct the StateBlocks to which associated ReactionBlocks will be linked. Reaction property packages must be tied to a single Physical property package, and this is used to validate the connections made later when constructing ReactionBlocks.
* default_arguments - this allows the user to provide a set of default values for construction arguments in associated ReactionBlocks, which will be passed to all ReactionBlocks when they are constructed.

.. module:: idaes.core.base.reaction_base

.. autoclass:: ReactionParameterBlock
    :members:

Reaction Blocks
---------------

Reaction Blocks are used within IDAES Unit models (generally within ControlVolume Blocks) in order to calculate reaction properties given the state of the material (provided by an associated StateBlock). Reaction Blocks are notably different to other types of Blocks within IDAES as they are always indexed by time (and possibly space as well), and are also not fully self contained (in that they depend upon the associated state block for certain variables). There are two bases Classes associated with Reaction Blocks:

* ReactionBlockDataBase forms the base class for all ReactionBlockData objects, which contain the instructions on how to construct each instance of a Reaction Block.
* ReactionBlockBase is used for building classes which contain methods to be applied to sets of Indexed Reaction Blocks (or to a subset of these). See the documentation on declare_process_block_class and the IDAES tutorials and examples for more information.

Reaction Block Construction Arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reaction Blocks have the following construction arguments:

* parameters - a reference to the associated Reaction Parameter block which will be used to make references to all necessary parameters.
* state_block - a reference to the associated StateBlock which will provide the necessary state and physical property information.
* has_equilibrium - indicates whether the associated Control Volume or Unit model expects chemical equilibrium to be enforced (if applicable).

Constructing Reaction Blocks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reaction Blocks can be constructed directly from the associated Reaction Parameter Block by calling the `build_reaction_block()` method on the Reaction Parameter Block. The `parameters` construction argument will be automatically set, and any other arguments (including indexing sets) may be provided to the `build_reaction_block` method as usual.

ReactionBlockDataBase Class
^^^^^^^^^^^^^^^^^^^^^^^^^^^

ReactionBlockDataBase contains the code necessary for implementing the as needed construction of variables and constraints.

.. autoclass:: ReactionBlockDataBase
    :members:

ReactionBlockBase Class
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: ReactionBlockBase
    :members:
