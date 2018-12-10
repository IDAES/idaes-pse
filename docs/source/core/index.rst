Core Library
============

Core Contents
-------------

.. toctree::
    :maxdepth: 1

    configuration
    process_block
    flowsheet_model
    state_block
    reaction_block
    unit_model
    control_volume

Core Overview
-------------

All components of the IDAES process modeling framework are built of Pyomo Block components (see Pyomo documentation).

The ProcessBlock class is the base class of IDAES models, and provides the common foundation for all other components.

FlowsheetModel objects represent the top level of the IDAES modeling hierarchy, and contain connected networks of unit models, or even contain other flowsheet models, which are connected by Pyomo Arcs.

Physical property packages supply information about a materials state including physical
properties and flow rates. Reaction property packages are used in systems where chemical
reactions may take place, and supply information on reaction rates and stoichiometry, based on a
materials state.

Equipment models are derived from UnitModel. Unit models contain control volumes
and have ports which can be used to connect material and energy flows between
unit models. On top of the balance equations usually contained in control
volumes unit models contain additional performance equations that may calculate
things like heat and mass transfer or efficiency curves.

ControlVolumes are the basic building block used construct unit models that
contains material and energy holdup and flows in and out. These block contain
energy, mass, and momentum balances, as well as state and reaction
blocks associated with the material within the control volume.

Component References
--------------------

There are many situations in the IDAES modeling framework where a developer may want to make use of a modeling component (e.g. a variable or parameter) from one Block in another Block. The time domain is a good example of this - almost all Blocks within an IDAES model need to make use of the time domain, however the time domain exists only at the top level of the flowsheet structure. In order to make use of the time domain in other parts of the framework, references to the time domain are used instead. By convention, all references within the IDAES modeling framework are indicated by the suffix "_ref" attached to the name of the reference. E.g. all references to the time domain within the framework are called "time_ref".

.. include:: ../global.rst
