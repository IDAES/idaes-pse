Equilibrium Constant Forms
==========================

.. contents:: Contents 
    :depth: 2

van 't Hoff Equation (van_t_hoff)
---------------------------------

The method uses the van 't Hoff equation to calculate the equilibrium constant.

.. math:: k_{eq} = k_{eq, ref} e^{-(\frac{\Delta H_{rxn}}{R})(\frac{1}{T} - \frac{1}{T_{ref, eq}})}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`k_{eq, ref}`", "k_eq_ref", "Varies, depends on equilibrium form", "Equlibrium constant at reference temperature"
   ":math:`T_{ref, eq}`", "temperature_eq_ref", "Base units", "Reference temperature for calculating equilibrium constant"

Gibbs Energy (gibbs_energy)
---------------------------

The method uses the thermodynamic relationship with Gibbs energy to calculate the equilibrium constant.

.. math:: k_{eq} = e^{-\frac{\Delta H_{rxn, ref}}{R T} + \frac{\Delta S_{rxn, ref}}{R}}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`\Delta H_{rxn, ref}`", "dh_rxn_ref", "Base units", "Heat of reaction at reference temperature"
   ":math:`\Delta S_{rxn, ref}`", "ds_rxn_ref", "Base units", "Entropy of reaction at reference temperature"
   ":math:`T_{ref, eq}`", "temperature_eq_ref", "Base units", "Reference temperature, used for calculating scaling factors"

Note that the above form is strictly unitless - units for the equilibrium constant are calculated separately and appended to the expression.
