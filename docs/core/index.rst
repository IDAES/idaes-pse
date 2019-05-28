Core Library
============

Core Contents
-------------

.. toctree::
    :maxdepth: 1

    configuration
    process_block
    model_concepts
    flowsheet_model
    properties
    unit_model
    control_volume
    util/index

Core Overview
-------------

All components of the IDAES process modeling framework are built of Pyomo Block components (see Pyomo documentation).

The ProcessBlock class is the base class of IDAES models, and provides the common foundation for all other components.

FlowsheetModel objects represent the top level of the IDAES modeling hierarchy, and contain connected networks of unit models, or even contain other flowsheet models, which are connected by Pyomo Arcs.

Physical property packages supply information about a material's state including physical
properties and flow rates. Reaction property packages are used in systems where chemical
reactions may take place, and supply information on reaction rates and stoichiometry, based on a
material's state.

Equipment models are derived from UnitModel. Unit models contain control volumes
and have ports which can be used to connect material and energy flows between
unit models. On top of the balance equations usually contained in control
volumes unit models contain additional performance equations that may calculate
things like heat and mass transfer or efficiency curves.

ControlVolumes are the basic building block used to construct unit models that
contain material and energy holdup and flows in and out. These blocks contain
energy, mass, and momentum balances, as well as state and reaction
blocks associated with the material within the control volume.

More detail on the different types of modeling objects is available in the Modeling Concepts section.
