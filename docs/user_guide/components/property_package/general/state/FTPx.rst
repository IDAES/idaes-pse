``FTPx``
========

.. contents:: Contents 
    :depth: 2

State Definition
----------------

This approach describes the material state in terms of total flow (:math:`F`: `flow_mol`), overall (mixture) mole fractions (:math:`x_j`: `mole_frac_comp`), temperature (:math:`T`: `temperature`) and pressure (:math:`P`: `pressure`). As such, there are :math:`3 + N_{components}` state variables, however only :math:`2 + N_{components}` are independent as the mole fraction must sum to 1.

Application
-----------

This is the simplest approach to fully defining the state of a material, and one of the most easily accessible to the user as it is defined in terms of variables that are easily measured and understood. However, this approach has a number of limitations which the user should be aware of:

* If the property package is set up for multiphase flow, an equilibrium calculation is required at the inlet of each unit, as the state definition does not contain information on multiphase flow. This increases the number of complex equilibrium calculations that must be performed, which could be avoided by using a different state definition. 
* State becomes ill-defined when only one component is present and multiphase behavior can occur, as temperature and pressure are insufficient to fully define the thermodynamic state under these conditions.

Bounds
------

The FTPx module supports bounding of the following variables through the `state_bounds` configuration argument:

* `flow_mol`
* `temperature`
* `pressure`

Note that mole fractions are automatically assigned a lower bound of 0, but the upper bound is left free as this is implicitly defined by the sum of mole fractions constraint. 

Supporting Variables and Constraints
------------------------------------

In addition to the state variables, this definition of state creates a number of supporting variables and constraints.

Variables
"""""""""

* `flow_mol_phase` (:math:`F_{mol, p}`)
* `mole_frac_phase_comp` (:math:`x_{p, j}`)
* `phase_frac` (:math:`\psi_p`)

Constraints
"""""""""""

In all cases, a constraint is written for the sum of the overall mole fractions.

.. math:: \sum_j{x_j} = 1

.. note::
   The sum of mole fractions constraint is not written at inlet states, as all mole fractions should be defined in the inlet stream.

If the property package supports only one phase:

.. math:: F_{mol, p} = F_{mol}
.. math:: x_{p, j} = x_{j} \text{ for all }j
.. math:: \psi_p = 1

If the property package supports only two phases, the Rachford-Rice formulation is used:

.. math:: \sum_p{F_{mol, p}} = F_{mol}
.. math:: F_{mol} \times x_{j} = \sum_p{F_{mol, p} \times x_{p, j}} \text{ for all }j
.. math:: \sum_j{x_{\text{phase 1}, j}} - \sum_j{x_{\text{phase 2}, j}} = 0
.. math:: \psi_p \times F_{mol} = F_{mol, p} \text{ for all }p

If the property package supports more than two phases, the following general formulation is used:

.. math:: F_{mol} \times x_{j} = \sum_p{F_{mol, p} \times x_{p, j}} \text{ for all }j
.. math:: \sum_j{x_{p, j}} = 1 \text{ for all }p
.. math:: \psi_p \times F_{mol} = F_{mol, p} \text{ for all }p

Default Balance Types and Flow Basis
------------------------------------

The following defaults are specified for Unit Models using this state definition:

* Material balances: total component balances
* Material flow basis: molar flow
* Energy balances: total enthalpy
