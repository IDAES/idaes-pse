Reaction Block
==============

Reaction Blocks are used within IDAES Unit models (generally within ControlVolume Blocks) in 
order to calculate reaction properties given the state of the material (provided by an 
associated StateBlock). Reaction Blocks are notably different to other types of Blocks within 
IDAES as they are always indexed by time (and possibly space as well), and are also not fully 
self contained (in that they depend upon the associated state block for certain variables). 
Reaction Blocks are composed of two parts:

* ReactionBlockDataBase forms the base class for all ReactionBlockData objects, which contain the instructions on how to construct each instance of a Reaction Block.
* ReactionBlockBase is used for building classes which contain methods to be applied to sets of Indexed Reaction Blocks (or to a subset of these). See the documentation on declare_process_block_class and the IDAES tutorials and examples for more information.

Reaction Blocks can be constructed directly from the associated Reaction Parameter Block by 
calling the `build_reaction_block()` method on the Reaction Parameter Block. The `parameters` 
construction argument will be automatically set, and any other arguments (including indexing 
sets) may be provided to the `build_reaction_block` method as usual.
