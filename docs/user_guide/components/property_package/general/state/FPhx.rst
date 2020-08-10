``FPhx``
========

.. contents:: Contents 
    :depth: 2

State Definition
----------------

This approach describes the material state in terms of total flow (:math:`F`: `flow_mol`), overall (mixture) mole fractions (:math:`x_j`: `mole_frac_comp`), total molar enthalpy (:math:`h`: `enth_mol`) and pressure (:math:`P`: `pressure`). As such, there are :math:`3 + N_{components}` state variables, however only :math:`2 + N_{components}` are independent as the mole fraction must sum to 1.

Application
-----------

This approach is commonly used by other process simulation tools as it avoids the issues associated with using temperature and pressure as state variables in single component systems. However, as the user generally does not know the specific enthalpy of their feed streams, this approach requires some method to calculate this for feed streams. This can generally be done by specifying temperature of the feed, and then solving for the specific enthalpy.

This approach suffers from the following limitation which the user should be aware of:

* If the property package is set up for multiphase flow, an equilibrium calculation is required at the inlet of each unit, as the state definition does not contain information on multiphase flow. This increases the number of complex equilibrium calculations that must be performed, which could be avoided by using a different state definition.

Bounds
------

The FPhx module supports bounding of the following variables through the `state_bounds` configuration argument:

* `flow_mol`
* `enth_mol`
* `pressure`
* `temperature`

Supplying bounds for temperature is supported as these are often known to greater accuracy than the enthalpy bounds, and specifying these can help the solver find a feasible solution.

Note that mole fractions are automatically assigned a lower bound of 0, but the upper bound is left free as this is implicitly defined by the sum of mole fractions constraint. 

Supporting Variables and Constraints
------------------------------------

In addition to the state variables, this definition of state creates a number of supporting variables and constraints.

Variables
"""""""""

* `flow_mol_phase` (:math:`F_{mol, p}`)
* `mole_frac_phase_comp` (:math:`x_{p, j}`)
* `temperature` (:math:`T`)
* `phase_frac` (:math:`\psi_p`)

Constraints
"""""""""""

In all cases, a constraint is written for the sum of the overall mole fractions.

.. math:: \sum_j{x_j} = 1

.. note::
   The sum of mole fractions constraint is not written at inlet states, as all mole fractions should be defined in the inlet stream.

Additionally, a constraint relating the total specific enthalpy to the specific enthalpy of each phase is written.

.. math:: h_{mol} = \sum_j{\psi_p \times h_{mol, p}}

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
