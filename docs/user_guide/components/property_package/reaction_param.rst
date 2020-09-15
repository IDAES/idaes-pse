Reaction Parameter Block
========================

ReactionParameterBlocks serve as a central location for linking to a property package, and 
contain all the parameters and indexing sets used by a given property package.

The role of the :ref:`ReactionParameterBlock Class<technical_specs/core/reaction_property_class:ReactionParameterBlock Class>` 
is to set up the references required by the rest of the IDAES framework for constructing 
instances of :ref:`ReactionBlocks<user_guide/components/property_package/reaction_block:Reaction Block>` 
and attaching these to the ReactionParameterBlock for ease of use. This allows other models to 
be pointed to the ReactionParameterBlock in order to collect the necessary information and to 
construct the necessary ReactionBlocks without the need for the user to do this manually.

Reaction property packages are used by all of the other modeling components to inform them of 
what needs to be constructed when dealing with chemical reactions. In order to do this, the 
IDAES modeling framework looks for a number of attributes in the ReactionParameterBlock which 
are used to inform the construction of other components. These attributes include:

* `reaction_block_class` - a pointer to the associated class that should be called when constructing ReactionBlocks. This should only be set by the property package developer.
* `phase_list` - a Pyomo Set object defining the valid phases of the mixture of interest.
* `component_list` - a Pyomo Set defining the names of the chemical species present in the mixture.
* `rate_reaction_idx` - a Pyomo Set defining a list of names for the kinetically controlled reactions of interest.
* `rate_reaction_stoichiometry` - a dict-like object defining the stoichiometry of the kinetically controlled reactions. Keys should be tuples of (rate_reaction_idx, phase_list, component_list) and values equal to the stoichiometric coefficient for that index.
* `equilibrium_reaction_idx` - a Pyomo Set defining a list of names for the equilibrium controlled reactions of interest.
* `equilibrium_reaction_stoichiometry` - a dict-like object defining the stoichiometry of the equilibrium controlled reactions. Keys should be tuples of (equilibrium_reaction_idx, phase_list, component_list) and values equal to the stoichiometric coefficient for that index.
* supported properties metadata - a list of supported reaction properties that the property package supports, along with instruction to construct the associated variables and constraints, and the units of measurement used for the property. This information is set using the add_properties attribute of the define_metadata class method.
* required properties metadata - a list of physical properties that the reaction property calculations depend upon, and must be supported by the associated StateBlock. This information is set using the add_required_properties attribute of the define_metadata class method.

