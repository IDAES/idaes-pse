Core Library
============

Core Contents
-------------

.. toctree::
    :maxdepth: 1

    configuration
    process_block
    state_block
    reaction_block
    control_volume
    unit_model
    flowsheet_model

Core Overview
-------------

The ProcessBlock class is the base class of IDAES models.

StateBlock supplies information about a materials state including physical
properties and flow rates. ReactionBlock is used in systems where chemical
reactions may take place and supplies reaction rate information based on a
materials state.

ControlVolume is a the basic building block used construct unit models that
contains material and energy holdup and flows in and out. These block contain
energy, mass, and momentum balances. Control volumes contain state and reaction
blocks.

Equipment models are derived from UnitModel. Unit models contain control volumes
and have ports which can be used to connect material and energy flows between
unit models. On top of the balance equations usually contained in control
volumes unit models contain additional performance equations that may calculate
things like heat and mass transfer or efficiency curves.

FlowsheetModel objects contain connected unit models and could even contain
other flowsheet models which are connected by streams.

.. include:: ../global.rst
