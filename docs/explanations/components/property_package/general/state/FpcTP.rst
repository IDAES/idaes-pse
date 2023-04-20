``FpcTP``
=========

.. contents:: Contents 
    :depth: 2

State Definition
----------------

This approach describes the material state in terms of phase-component flow (:math:`F_{p, j}`: `flow_mol_phase_comp`), temperature (:math:`T`: `temperature`) and pressure (:math:`P`: `pressure`). As such, there are :math:`2 + phases*N_{components}` state variables.

Application
-----------

This approach required knowledge of the phase-equilibrium of the material in order to define the state variables. Compared to using total flow rate and mole fractions, this approach contains full information on the phase equilibria within the state variables, and thus avoids the needs for flash calculations in many cases. This can greatly reduce the complexity of the problem, and results can significantly affect the tractablity of the problem. However, this approach has a number of limitations which the user should be aware of:

* Users must have knowledge of , or calculate, the phase-component flows of all inlet streams. For single phase flows this is often known, but for streasm with potetnial two-phase behaviour this can reqruire a set of flash calculations for the feed streasm (users can make use of Feed blocks to assist with this).
* State becomes ill-defined when only one component is present and multiphase behavior can occur, as temperature and pressure are insufficient to fully define the thermodynamic state under these conditions.

Bounds
------

The FpcTP module supports bounding of the following variables through the `state_bounds` configuration argument:

* `flow_mol_phase_comp`
* `temperature`
* `pressure`

Supporting Variables and Constraints
------------------------------------

In addition to the state variables, this definition of state creates a number of supporting variables and constraints.

Variables
"""""""""

* `mole_frac_phase_comp` (:math:`x_{p, j}`)

Expressions
"""""""""""

* `flow_mol` (:math:`F = \sum{F_{p,j}}`)
* `flow_mol_phase` (:math:`F_p = \sum{F_{p,j}}_j`)
* `flow_mol_comp` (:math:`F_j = \sum{F_{p,j}}_p`)
* `mole_frac_comp` (:math:`x_j = \frac{\sum{F_{p,j}}}{F}`)
* `phase_frac` (:math:`\psi_p = \frac{F_p}{F}` or :math:`\psi_p = 1` if only single phase)

Constraints
"""""""""""

A set of constraints is created to calculate phase-component mole fractions from the phase-component flow rates.

.. math:: F_j = x_{p, j} \times \sum{F_p} = F_{p, j}

Default Balance Types and Flow Basis
------------------------------------

The following defaults are specified for Unit Models using this state definition:

* Material balances: total component balances
* Material flow basis: molar flow
* Energy balances: total enthalpy
