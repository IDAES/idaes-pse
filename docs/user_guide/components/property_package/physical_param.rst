Physical Parameter Block
========================

PhysicalParameterBlocks serve as a central location for linking to a property package, and 
contain all the parameters and indexing sets used by a given property package.

The role of the :ref:`PhysicalParameterBlock Class<technical_specs/core/physical_property_class:PhysicalParameterBlock Class>` 
is to set up the references required by the rest of the IDAES Core Modeling Framework for constructing 
instances of :ref:`StateBlocks<user_guide/components/property_package/state_block:State Block>` 
and attaching these to the PhysicalParameterBlock for ease of use. This allows other models to 
be pointed to the PhysicalParameterBlock in order to collect the necessary information and to 
construct the necessary StateBlocks without the need for the user to do this manually.

Several attributes in the PhysicalParameterBlock are used to 
inform the construction of other components. These attributes include:

* `state_block_class` - a pointer to the associated class that should be called when constructing StateBlocks. This should only be set by the property package developer.
* `phase_list` - a Pyomo Set object defining the valid phases of the mixture of interest.
* `component_list` - a Pyomo Set defining the names of the chemical species present in the mixture.
* `element_list` - (optional) a Pyomo Set defining the names of the chemical elements that make up the species within the mixture. This is used when doing elemental material balances.
* `element_comp` - (optional) a dictionary-like object which defines the elemental composition of each species in component_list. Form: component: {element_1: value, element_2: value, ...}.
* supported properties metadata - a dictionary of supported physical properties that the property package supports, along with instruction to construct the associated variables and constraints, and the units of measurement used for the property. This information is set using the `add_properties` attribute of the `define_metadata` class method.


