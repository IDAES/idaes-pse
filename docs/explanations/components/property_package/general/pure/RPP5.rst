Properties of Gases and Liquids 5th edition (``RPP5``)
======================================================

.. contents:: Contents
    :depth: 2

Source
------

Methods for calculating pure component properties from:

The Properties of Gases & Liquids, 5th Edition
Reid, Prausnitz and Polling, 2001, McGraw-Hill

All methods use SI units.

Ideal Gas Molar Heat Capacity (Constant Pressure)
-------------------------------------------------

Properties of Gases and Liquids uses the following correlation for the ideal gas molar heat capacity:

.. math:: \frac{c_{\text{p ig}}}{R} = a0 + a1 \times T + a2 \times T^2 + a3 \times T^3 + a4 \times T^4

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`a0`", "cp_mol_ig_comp_coeff_a0", "None", ""
   ":math:`a1`", "cp_mol_ig_comp_coeff_a1", ":math:`\text{K}^-1`", ""
   ":math:`a2`", "cp_mol_ig_comp_coeff_a2", ":math:`\text{K}^-2`", ""
   ":math:`a3`", "cp_mol_ig_comp_coeff_a3", ":math:`\text{K}^-3`", ""
   ":math:`a4`", "cp_mol_ig_comp_coeff_a4", ":math:`\text{K}^-4`", ""
   ":math:`R`", "gas_constant", "Same as heat capacity", "Universal gas constant"

Ideal Gas Molar Enthalpy
------------------------

The correlation for the ideal gas molar enthalpy is derived from the correlation for the molar heat capacity and is given below:

.. math:: \frac{h_{\text{ig}} - h_{\text{ig ref}}}{R} = a0 \times (T-T_{ref}) + \frac{a1}{2} \times (T^2 - T_{ref}^2) + \frac{a2}{3} \times (T^3 - T_{ref}^3) + \frac{a3}{4} \times (T^4 - T_{ref}^4) + \frac{a4}{5} \times (T^5 - T_{ref}^5) + \Delta h_{\text{form, Vap}}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`a0`", "cp_mol_ig_comp_coeff_a0", "None", ""
   ":math:`a1`", "cp_mol_ig_comp_coeff_a1", ":math:`\text{K}^-1`", ""
   ":math:`a2`", "cp_mol_ig_comp_coeff_a2", ":math:`\text{K}^-2`", ""
   ":math:`a3`", "cp_mol_ig_comp_coeff_a3", ":math:`\text{K}^-3`", ""
   ":math:`a4`", "cp_mol_ig_comp_coeff_a4", ":math:`\text{K}^-4`", ""
   ":math:`\Delta h_{\text{form, Vap}}`", "enth_mol_form_vap_comp_ref", ":math:`\text{J/mol}`", "Molar heat of formation at reference state"

.. note::
    This correlation uses the same parameters as the ideal gas heat capacity correlation.

Ideal Gas Molar Entropy
------------------------

The correlation for the ideal gas molar entropy is derived from the correlation for the molar heat capacity and is given below:

.. math:: \frac{s_{\text{ig}}}{R}= a0 \times ln(T/T_{ref}) + a1 \times (T - T_{ref}) + \frac{a2}{2} \times (T^2 - T_{ref}^2) + \frac{a3}{3} \times (T^3 - T_{ref}^3) + \frac{a4}{4} \times (T^4 - T_{ref}^4) + s_{\text{form, Vap}}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`a0`", "cp_mol_ig_comp_coeff_a0", "None", ""
   ":math:`a1`", "cp_mol_ig_comp_coeff_a1", ":math:`\text{K}^-1`", ""
   ":math:`a2`", "cp_mol_ig_comp_coeff_a2", ":math:`\text{K}^-2`", ""
   ":math:`a3`", "cp_mol_ig_comp_coeff_a3", ":math:`\text{K}^-3`", ""
   ":math:`a4`", "cp_mol_ig_comp_coeff_a4", ":math:`\text{K}^-4`", ""
   ":math:`s_{\text{form, Vap}}`", "entr_mol_form_vap_comp_ref", ":math:`\text{J/mol}\cdotp\text{K}`", "Standard molar entropy of formation at reference state"

.. note::
    This correlation uses the same parameters as the ideal gas heat capacity correlation.

Saturation (Vapor) Pressure
---------------------------

Properties of Gases and Liquids 5th edition uses the following correlation to calculate the vapor pressure of a component:

.. math:: Log{(P_{sat}) = A - \frac{B}{T+C}}

Units are bar and Kelvin.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "pressure_sat_comp_coeff_A", "None", ""
   ":math:`B`", "pressure_sat_comp_coeff_B", ":math:`\text{K}`", ""
   ":math:`C`", "pressure_sat_comp_coeff_C", ":math:`\text{K}`", ""
