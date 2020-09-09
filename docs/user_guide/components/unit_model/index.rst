Unit Model
==========

.. toctree::
    :glob:
    :hidden:
    
    *

Unit models represent pieces of equipment and their processes. They include any unit performance 
constraints and associated variables, such as:

    - constraints relating balance terms to physical phenomena or properties (e.g. relating extent of reaction to reaction rate and volume)
    - constraints describing flow of material into or out of unit (e.g. pressure driven flow constraints)
    - unit level efficiency constraints (e.g. relating mechanical work to fluid work)

Unit models are composed of the following components:

    1. :ref:`ControlVolumeBlocks<user_guide/components/unit_model/control_volume:Control Volume>`, which represent volume of material over which we wish to perform material, energy and/or momentum balances
    2. :ref:`StateBlocks<user_guide/components/property_package/state_block:State Block>` and :ref:`ReactionBlocks<user_guide/components/property_package/reaction_block:Reaction Block>`, which represent the thermophysical, transport and reaction properties of the material at a specific point in space and time
    3. Inlets and Outlets, which allow unit models to connect to other unit models



