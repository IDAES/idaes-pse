IDAES Modeling Concepts
=======================

.. contents:: Contents 
    :depth: 2

Introduction
------------

The purpose of this section of the documentation to explain the different parts of the IDAES modeling framework, and what components belong in each part for the hierarchy. Each component is described in greater detail later in the documentation, however this section provides a general introduction to different types of components.

Time Domain
-----------

Before starting on the different types of models present in the IDAES framework, it is important to discuss how time is handled by the framework. When a user first declares a Flowsheet model a time domain is created, the form of which depends on whether the Flowsheet is declared to be dynamic or steady-state (see FlowsheetBlock documentation). When other models are added to the Flowsheet, or to a descendant of the Flowsheet, that model makes a reference to the original time domain. This is handled automatically by the framework so that the user does not need to worry about it. This ensures that all models within the Flowsheet have a consistent time domain.

Different models may handle the time domain differently, but in general all IDAES models contain a component named time, which is a reference to the original time domain. The only exception to this are blocks associated with Property calculations. PropertyBlocks represent the state of the material at a single point in space and time, and thus do not contain the time domain. Instead, PropertyBlocks are index by time (and space where applicable) - i.e. there is a separate PropertyBlock for each point in time. The user should keep this in mind when working with IDAES models, as it is important for understanding where the time index appears within a model.

Another important thing to note is that steady-state models do contain a time domain, however this is generally a single point at time = 0.0. However, models still contain a reference to the time domain, and any components are still indexed by time even in a steady-state model (e.g. PropertyBlocks).

Flowsheets
----------

The top level of the IDAES modeling framework is the Flowsheet model. Flowsheet models represent traditional process flowsheets, containing a number of Unit models representing process unit operations connected together into a flow network. Flowsheets generally contain three types of component:

1. Unit models, representing unit operations,
2. Arcs, representing connections between Unit models, and,
3. Property Parameter blocks, representing the parameters associated with different materials present within the flowsheet.

Flowsheet models may also contain additional constraints relating to how different Unit models behave and interact, such as control and operational constraints. Generally speaking, if a Constraint is purely internal to a single unit, and does not depend on information from other units in the flowsheet, then the Constraint should be placed inside the relevant Unit model. Otherwise, the Constraint should be placed at the Flowsheet level.

Unit Models
-----------

Unit models generally represent individual pieces of equipment present within a process which perform a specific task. Unit models in turn are generally composed of two main types of component:

1. Control Volume Blocks, which represent volume of material over which we wish to perform material, energy and/or momentum balances, and,
2. StateBlocks and ReactionBlocks, which represent the thermophysical, transport and reaction properties of the material at a specific point in space and time.
3. Inlets and Outlets, which allow Unit models to connect to other Unit models.

Unit models will also contain Constraints describing the performance of the unit, which will relate terms in the balance equations to different phenomena.

Control Volumes
^^^^^^^^^^^^^^^

A key feature of the IDAES modeling framework is the use of Control Volume Blocks. As mentioned above, Control Volumes represent a volume of material over which material, energy and/or momentum balances can be performed. Control Volume Blocks contain methods to automate the task of writing common forms of these balance equations. Control Volumes Blocks can also automate the creation of StateBlocks and ReactionBlocks associated with the control volume.

Property Blocks
^^^^^^^^^^^^^^^

Property blocks represent the state of a material at a given point in space and time within the process flowsheet, and contain the state variables, thermophsyical, transport and reaction properties of a material (which are functions solely of the local state of the material). Within the IDAES proces modeling framework, properties are divinded into two types:

* Physcial properties (StateBlocks), including thermophysical and transport properties, and
* Reaction properties (ReactionBlocks), which includes all properties assoicated with chemical reactions.

Additionally, StateBlocks contain information on the extensive flow of material at that point in space and time, which is a departure from how engineers generally think about properties. This is required to facilitate the flexible formulation of the IDAES Framework by allowing the property package to dictate what form the balance equations will take, which requires the StateBlock to know the extensive flow information.

The calculations involved in property blocks of both types generally require a set of parameters which are constant across all instances of that type of property block. Rather than each property block containing its own copy of each of these parameters (thus duplicating parameters between blocks), each type of property block is associated with a Property Parameter Block (PhysicalParameterBlock or ReactionParameterBlock). Property Parameter Blocks serve as a centralized location for the constant parameters involved in property calculations, and all property blocks of the associated type link to the parameters contain in the parameter block.

Component References
--------------------

There are many situations in the IDAES modeling framework where a developer may want to make use of a modeling component (e.g. a variable or parameter) from one Block in another Block. The time domain is a good example of this - almost all Blocks within an IDAES model need to make use of the time domain, however the time domain exists only at the top level of the flowsheet structure. In order to make use of the time domain in other parts of the framework, references to the time domain are used instead. By convention, all references within the IDAES modeling framework are indicated by the suffix "_ref" attached to the name of the reference. E.g. all references to the time domain within the framework are called "time_ref".

What Belongs in Each Type of Block?
-----------------------------------

A common question with the hierarchical structure of the IDAES framework is where does a specific variable or constraint belong (or conversely, where can I find a specific variable or constraint). In general, variables and constraints are divided based on the following guidelines:

1. Property Parameter Blocks - any parameter or quantity that is consistent across all instances of a Property Block belongs in the Property Parameter Block. This includes:

    - component lists,
    - lists of valid phases,
    - universal constants (e.g. R, :math:`\pi`),
    - constants used in calculating properties (e.g. coefficients for calculating :math:`c_p`,
    - reference states (e.g. :math:`P_{ref}` and :math:`T_{ref}`),
    - lists of reaction identifiers,
    - reaction stoichiometry.

2. Property Blocks - all state variables (including extensive flow information) and any quantity that is a function only of state variables plus the constraints required to calculate these. These include:

    - flow rates (can be of different forms, e.g. mass or molar flow, on a total or component basis),
    - temperature,
    - pressure,
    - intensive and extensive state functions (e.g. enthalpy); both variables and constraints.

3. Control Volume Blocks - material, energy and momentum balances and the associated terms. These include:

    - balance equations,
    - holdup volume,
    - material and energy holdups; both variables and constraints,
    - material and energy accumulation terms (Pyomo.dae handles the creation of the associated derivative constraints),
    - material generation terms (kinetic reactions, chemical and phase equilibrium, mass transfer),
    - extent of reaction terms and constraints relating these to the equivalent generation terms,
    - phase fraction within the holdup volume and constrain on the sum of phase fractions,
    - heat and work transfer terms,
    - pressure change term
    - diffusion and conduction terms (where applicable) and associated constraints,
    - Mixer and Splitter blocks for handling multiple inlets/outlets.

4. Unit Model - any unit performance constraints and associated variables, such as:

    - constraints relating balance terms to physical phenomena or properties (e.g. relating extent of reaction to reaction rate and volume),
    - constraints describing flow of material into or out of unit (e.g. pressure driven flow constraints),
    - unit level efficiency constraints (e.g. relating mechanical work to fluid work).

5. Flowsheet Model - any constraints related to interaction of unit models and associated variables. Examples include:

    - control constraints relating behavior between different units (e.g. a constraint on valve opening based on the level in another unit).

