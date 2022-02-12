Properties of Gases and Liquids 3rd edition (``RPP3``)
======================================================

.. contents:: Contents 
    :depth: 2

Source
------

Methods for calculating pure component properties from:

The Properties of Gases & Liquids, 3rd Edition
Reid, Prausnitz and Polling, 1977, McGraw-Hill

Ideal Gas Molar Heat Capacity (Constant Pressure)
-------------------------------------------------

Properties of Gases and Liquids uses the following correlation for the ideal gas molar heat capacity:

.. math:: c_{\text{p ig}} = A + B \times T + C \times T^2 + D \times T^3

Units are calories per gram-mole kelvin and Kelvin.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "cp_mol_ig_comp_coeff_A", ":math:`\text{cal/mol}\cdotp\text{K}`", ""
   ":math:`B`", "cp_mol_ig_comp_coeff_B", ":math:`\text{cal/mol}\cdotp\text{K}^2`", ""
   ":math:`C`", "cp_mol_ig_comp_coeff_C", ":math:`\text{cal/mol}\cdotp\text{K}^3`", ""
   ":math:`D`", "cp_mol_ig_comp_coeff_D", ":math:`\text{cal/mol}\cdotp\text{K}^4`", ""

Ideal Gas Molar Enthalpy
------------------------

The correlation for the ideal gas molar enthalpy is derived from the correlation for the molar heat capacity and is given below:

.. math:: h_{\text{ig}} - h_{\text{ig ref}} = A \times (T-T_{ref}) + \frac{B}{2} \times (T^2 - T_{ref}^2) + \frac{C}{3} \times (T^3 - T_{ref}^3) + \frac{D}{4} \times (T^4 - T_{ref}^4) + \Delta h_{\text{form, Vap}}

Units are calories per gram-mole kelvin and Kelvin.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "cp_mol_ig_comp_coeff_A", ":math:`\text{cal/mol}\cdotp\text{K}`", ""
   ":math:`B`", "cp_mol_ig_comp_coeff_B", ":math:`\text{cal/mol}\cdotp\text{K}^2`", ""
   ":math:`C`", "cp_mol_ig_comp_coeff_C", ":math:`\text{cal/mol}\cdotp\text{K}^3`", ""
   ":math:`D`", "cp_mol_ig_comp_coeff_D", ":math:`\text{cal/mol}\cdotp\text{K}^4`", ""
   ":math:`\Delta h_{\text{form, Vap}}`", "enth_mol_form_vap_comp_ref", "Base units", "Molar heat of formation at reference state"

.. note::
    This correlation uses the same parameters as the ideal gas heat capacity correlation.
    Units of molar heat of formation will be derived from the base units defined for the property package.

Ideal Gas Molar Entropy
------------------------

The correlation for the ideal gas molar entropy is derived from the correlation for the molar heat capacity and is given below:

.. math:: s_{\text{ig}} = A \times ln(T/T_{ref}) + B \times (T - T_{ref}) + \frac{C}{2} \times (T^2 - T_{ref}^2) + \frac{D}{3} \times (T^3 - T_{ref}^3) + s_{\text{form, Vap}}

Units are calories per gram-mole kelvin and Kelvin.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "cp_mol_ig_comp_coeff_A", ":math:`\text{cal/mol}\cdotp\text{K}`", ""
   ":math:`B`", "cp_mol_ig_comp_coeff_B", ":math:`\text{cal/mol}\cdotp\text{K}^2`", ""
   ":math:`C`", "cp_mol_ig_comp_coeff_C", ":math:`\text{cal/mol}\cdotp\text{K}^3`", ""
   ":math:`D`", "cp_mol_ig_comp_coeff_D", ":math:`\text{cal/mol}\cdotp\text{K}^4`", ""
   ":math:`s_{\text{form, Vap}}`", "entr_mol_form_vap_comp_ref", "Base units", "Standard molar entropy of formation at reference state"

.. note::
    This correlation uses the same parameters as the ideal gas heat capacity correlation .
    Units of molar entropy of formation will be derived from the base units defined for the property package.

Saturation (Vapor) Pressure
---------------------------

Properties of Gases and Liquids 3rd edition uses the following correlation to calculate the vapor pressure of a component:

.. math:: Ln{(P_{sat}) = A - \frac{B}{T+C}}

Units are mmHg and Kelvin.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "pressure_sat_comp_coeff_A", "None", ""
   ":math:`B`", "pressure_sat_comp_coeff_B", ":math:`\text{K}`", ""
   ":math:`C`", "pressure_sat_comp_coeff_C", ":math:`\text{K}`", ""

