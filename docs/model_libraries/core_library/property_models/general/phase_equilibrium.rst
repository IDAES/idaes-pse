Defining Phase Equilibria
=========================

Phase equilibrium and separation is a key part of almost all chemical processes, and also represent some of the most complex and non-linear constraints in a model, especially when dealing with systems which may cross phase boundaries. Systems may also include multiple interacting phases with equilibrium, which further complicates the problem. As such, good formulations of these constraints is key to a robust and tractable model.

The IDAES Generic Property Package framework supports a range of phase equilibrium behaviors, including multiple phases in equilibrium and different formulations for describing the equilibria. These are all optional, and users do not need to define phase equilibria if it is not required for their system.

Setting up phase equilibrium within the framework is done using three configuration arguments as discussed below. However, users should be aware that some of these options require definition of further properties, such as bubble and dew point calculations.

Define Phases in Equilibrium
----------------------------

The first step in setting up phase equilibrium in the framework is to describe which phases are in equilibrium with each other. In general, for phases to be in equilibrium with each other, the following conditions need to be met:

1. Phases must be in direct contact with each other, and
2. At least one component must appear in both phases (equilibria involving chemical reactions are handled by `ReactionBlocks`).

In order to describe which phases are in equilibrium, the user needs to set the `phase_in_equilibrium` construction argument, which should be a list of 2-tuples where each tuple describes a pair of phases which are in equilibrium. Any component which appears in both phases in a pair is assumed to be in equilibrium.

A simple example for a VLE system is shown below.

.. code-block:: python

    "phases_in_equilibrium" = [("Vap", "Liq")]

.. note::
    Users should take care not to over define their system. For example, in a VLSE system a user could potentially write three sets of equilibrium constraints (VL, LS and VS). However, this would result in an over defined system, as only two of these three are independent. For most situations, a user would consider only the VL and LS equilibria, with the VS being implicitly  defined.

Define Equilibrium State Formulation
------------------------------------

Next, for each pair of phases in equilibrium, the user must define a formulation for the equilibrium state. To handle the complexities of disappearing phases, the IDAES Generic Property Package Framework allows for phase equilibrium to be solved at a separate equilibrium state rather than the actual state of the material. This allows for formulations which avoid disappearing phases by limiting the equilibrium state to exist within the valid two-phase region, whilst returning a negligible amount of any phase which is not valid at the actual material state.

The equilibrium state formulation is set using the `phase_equilibrium_state` configuration argument. This should be a `dict` where the keys are 2-tuples of phases in equilibrium (matching those defined in the `phases_in_equilibrium` argument) and values are a phase equilibrium state formulation method. The IDAES Generic Property Package Framework contains a library of methods for the formulation of the phase equilibrium state, which is shown below.

Phase Equilibrium State Libraries
"""""""""""""""""""""""""""""""""

.. toctree::
    :maxdepth: 1

    pe/smooth_flash

Necessary Properties
--------------------

Next, any component which is involved in a phase equilibrium interaction (i,e, appears in both phases of an interacting pair) must define a form for the required equilibrium constraint. There are a number of ways these constraints can be written depending on the equation of state and scaling of the problem. This is set using the `phase_equilibrium_form` configuration argument in the `Component` objects, and takes the form of a `dict` where the keys are 2-tuples of interacting phases and the value is the formulation to use for the current component across the given phase pair. For example:

.. code-block:: python

    parameters.component_1.config.phase_equilibrium_form = {(phase_1, phase_2): formulation}

A library of common forms for equilibrium constraints is available, and is shown below.

.. toctree::
    :maxdepth: 2

    pe/pe_forms

Bubble and Dew Point Calculations
---------------------------------

Bubble and dew points are often of interest to process engineers for designing process equipment, and appear in some calculations of other thermodynamic properties. They are also useful in getting initial guesses for states in phase equilibrium problems, and some equilibrium state formulations rely on these properties.

Whilst calculation of the saturation pressure for single components is relatively simple, calculating the bubble and dew points of mixtures is more challenging due to the non-linear nature of the equations. Calculation of these properties is generally done through calculations based on the equations of state for the liquid and vapor phases, however these calculations can be greatly simplified if ideal behavior is assumed for both phases (i.e. ideal gas and Raoult's law). To allow for both cases, the IDAES Generic Property Package Framework provides a library of different formulations for the bubble and dew point calculations, which can be set using the following arguments:

* `bubble_dew_method`

A list of available methods is given below:

.. toctree::
    :maxdepth: 3

    bubble_dew
