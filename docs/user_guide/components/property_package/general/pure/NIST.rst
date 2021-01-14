NIST Webbook (``NIST``)
=======================

.. contents:: Contents 
    :depth: 2

Source
------

Pure component properties as used by the NIST WebBook, `<https://webbook.nist.gov/chemistry/>`_ Retrieved: September 13th, 2019

Ideal Gas Molar Heat Capacity (Constant Pressure)
-------------------------------------------------

NIST uses the Shomate equation for the ideal gas molar heat capacity, which is shown below:

.. math:: c_{\text{p ig}} = A + B \times t + C \times t^2 + D \times t^3 + \frac{E}{t^2}

where :math:`t = \frac{T}{1000}`. Units are :math:`\text{J/mol}\cdotp\text{K}`.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "cp_mol_ig_comp_coeff_A", ":math:`\text{J/mol}\cdotp\text{K}`", ""
   ":math:`B`", "cp_mol_ig_comp_coeff_B", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}`", ""
   ":math:`C`", "cp_mol_ig_comp_coeff_C", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}^2`", ""
   ":math:`D`", "cp_mol_ig_comp_coeff_D", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}^3`", ""
   ":math:`E`", "cp_mol_ig_comp_coeff_E", ":math:`\text{J}\cdotp\text{kK}^2\text{/mol}\cdotp\text{K}`", ""
   ":math:`F`", "cp_mol_ig_comp_coeff_F", ":math:`\text{kJ/mol}`", ""
   ":math:`G`", "cp_mol_ig_comp_coeff_G", ":math:`\text{J/mol}\cdotp\text{K}`", ""
   ":math:`H`", "cp_mol_ig_comp_coeff_H", ":math:`\text{kJ/mol}`", ""

.. note::
    Due to the division of temperature by 1000 in the expression form, most temperature units are in kilo-Kelvins and reference enthalpies (F and H) are in kJ/mol.
    The parameter `cp_mol_ig_comp_coeff` is also used when calculating specific enthalpy and entropy and parameters 'F', 'G' and 'H' are only required for these properties.

Ideal Gas Molar Enthalpy
------------------------

The correlation for the ideal gas molar enthalpy is derived from the correlation for the molar heat capacity and is given below:

.. math:: \frac{h_{\text{ig}} - h_{\text{ig ref}}}{1000} = A \times (t-t_{ref}) + \frac{B}{2} \times (t^2 - t_{ref}^2) + \frac{C}{3} \times (t^3 - t_{ref}^3) + \frac{D}{4} \times (t^4 - t_{ref}^4) + E \times (\frac{1}{t} - \frac{1}{t_{ref}}) + F - H

Units are :math:`\text{J/mol}`.

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "cp_mol_ig_comp_coeff_A", ":math:`\text{J/mol}\cdotp\text{K}`", ""
   ":math:`B`", "cp_mol_ig_comp_coeff_B", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}`", ""
   ":math:`C`", "cp_mol_ig_comp_coeff_C", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}^2`", ""
   ":math:`D`", "cp_mol_ig_comp_coeff_D", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}^3`", ""
   ":math:`E`", "cp_mol_ig_comp_coeff_E", ":math:`\text{J}\cdotp\text{kK}^2\text{/mol}\cdotp\text{K}`", ""
   ":math:`F`", "cp_mol_ig_comp_coeff_F", ":math:`\text{kJ/mol}`", ""
   ":math:`G`", "cp_mol_ig_comp_coeff_G", ":math:`\text{J/mol}\cdotp\text{K}`", ""
   ":math:`H`", "cp_mol_ig_comp_coeff_H", ":math:`\text{kJ/mol}`", ""

.. note::
    This correlation uses the same parameters as for the ideal gas heat capacity with additional parameters `F` and `H`. These parameters account for the enthalpy at the reference state defined by NIST, where `F` is the constant of integration and `H` is the standard molar heat of formation. Note that the default form of the expression used by NIST subtracts the heat of formation from the specific enthalpy. This behavior can be controlled using the global configuration argument `include_enthalpy_of_formation` - if this is set to `True` (the default setting), then the `H` term is not used when calculating specific enthalpies.
    Due to the division of temperature by 1000 in the expression form, most temperature units are in kilo-Kelvins and reference enthalpies (F and H) are in kJ/mol.

Ideal Gas Molar Entropy
------------------------

The correlation for the ideal gas molar entropy is derived from the correlation for the molar heat capacity and is given below:

.. math:: s_{\text{ig}} = A \times ln(t) + B \times t + \frac{C}{2} \times t^2 + \frac{D}{3} \times t^3 + \frac{E}{2 \times t^2} + G 

Units are :math:`\text{J/mol}\cdotp\text{K}`.

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "cp_mol_ig_comp_coeff_A", ":math:`\text{J/mol}\cdotp\text{K}`", ""
   ":math:`B`", "cp_mol_ig_comp_coeff_B", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}`", ""
   ":math:`C`", "cp_mol_ig_comp_coeff_C", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}^2`", ""
   ":math:`D`", "cp_mol_ig_comp_coeff_D", ":math:`\text{J/mol}\cdotp\text{K}\cdotp\text{kK}^3`", ""
   ":math:`E`", "cp_mol_ig_comp_coeff_E", ":math:`\text{J}\cdotp\text{kK}^2\text{/mol}\cdotp\text{K}`", ""
   ":math:`F`", "cp_mol_ig_comp_coeff_F", ":math:`\text{kJ/mol}`", ""
   ":math:`G`", "cp_mol_ig_comp_coeff_G", ":math:`\text{J/mol}\cdotp\text{K}`", ""
   ":math:`H`", "cp_mol_ig_comp_coeff_H", ":math:`\text{kJ/mol}`", ""

.. note::
    This correlation uses the same parameters as for the ideal gas heat capacity with additional parameter `G`, which accounts for the standard entropy at the reference state defined by NIST. Users wanting to use a different reference state will need to update `G`.
    Due to the division of temperature by 1000 in the expression form, most temperature units are in kilo-Kelvins and reference enthalpies (F and H) are in kJ/mol.

Saturation (Vapor) Pressure
---------------------------

NIST uses the Antoine equation to calculate the vapor pressure of a component, which is given below:

.. math:: log_{10}(P_{sat}) = A - \frac{B}{T+C}

Units are bar and Kelvin.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`A`", "pressure_sat_comp_coeff_A", "None", ""
   ":math:`B`", "pressure_sat_comp_coeff_B", ":math:`\text{K}`", ""
   ":math:`C`", "pressure_sat_comp_coeff_C", ":math:`\text{K}`", ""

