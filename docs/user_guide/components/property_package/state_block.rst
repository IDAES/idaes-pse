State Block
===========

State Blocks are used within all IDAES Unit models (generally within ControlVolume Blocks) in 
order to calculate physical properties given the state of the material. State Blocks are 
notably different to other types of Blocks within IDAES as they are always indexed by time 
(and possibly space as well). State Blocks consist of two parts:

* StateBlockData forms the base class for all StateBlockData objects, which contain the instructions on how to construct each instance of a State Block.
* StateBlock is used for building classes which contain methods to be applied to sets of Indexed State Blocks (or to a subset of these). See the documentation on declare_process_block_class and the IDAES tutorials and examples for more information.

State Blocks can be constructed directly from the associated Physical Parameter Block by calling 
the `build_state_block()` method on the Physical Parameter Block. The `parameters` construction 
argument will be automatically set, and any other arguments (including indexing sets) may be 
provided to the `build_state_block` method as ususal.

Additional details on :ref:`State Blocks<technical_specs/core/physical_property_class:State Blocks>`
are located in the technical specifications.


