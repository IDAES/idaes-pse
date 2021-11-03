Perry's Chemical Engineers' Handbook (``Perrys``)
=================================================

.. contents:: Contents 
    :depth: 2

Source
------

Methods for calculating pure component properties from:

Perry's Chemical Engineers' Handbook, 7th Edition
Perry, Green, Maloney, 1997, McGraw-Hill

.. note::
    Currently, only the most common correlation forms from Perry's Handbook are implemented. Some components use different forms which are not yet supported.

Ideal Liquid Molar Heat Capacity (Constant Pressure)
----------------------------------------------------

Perry's Handbook uses the following correlation for ideal liquid molar heat capacity:

.. math:: c_{\text{p liq}} = C_1 + C_2 \times T + C_3 \times T^2 + C_4 \times T^3 + C_5 \times T^4

Units are :math:`\text{J/kmol}\cdotp\text{K}`.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`C_1`", "cp_mol_ig_comp_coeff_1", ":math:`\text{J/kmol}\cdotp\text{K}`", ""
   ":math:`C_2`", "cp_mol_ig_comp_coeff_2", ":math:`\text{J/kmol}\cdotp\text{K}^2`", ""
   ":math:`C_3`", "cp_mol_ig_comp_coeff_3", ":math:`\text{J/kmol}\cdotp\text{K}^3`", ""
   ":math:`C_4`", "cp_mol_ig_comp_coeff_4", ":math:`\text{J/kmol}\cdotp\text{K}^4`", ""
   ":math:`C_5`", "cp_mol_ig_comp_coeff_5", ":math:`\text{J/kmol}\cdotp\text{K}^5`", ""

Ideal Liquid Molar Enthalpy
---------------------------

The correlation for the ideal liquid molar enthalpy is derived from the correlation for the molar heat capacity and is given below:

.. math:: h_{\text{liq}} - h_{\text{liq ref}} = C_1 \times (T-T_{ref}) + \frac{C_2}{2} \times (T^2 - T_{ref}^2) + \frac{C_3}{3} \times (T^3 - T_{ref}^3) + \frac{C_4}{4} \times (T^4 - T_{ref}^4) + \frac{C_5}{5} \times (T^5 - T_{ref}^5) + \Delta h_{\text{form, Liq}}

Units are :math:`\text{J/kmol}`.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`C_1`", "cp_mol_ig_comp_coeff_1", ":math:`\text{J/kmol}\cdotp\text{K}`", ""
   ":math:`C_2`", "cp_mol_ig_comp_coeff_2", ":math:`\text{J/kmol}\cdotp\text{K}^2`", ""
   ":math:`C_3`", "cp_mol_ig_comp_coeff_3", ":math:`\text{J/kmol}\cdotp\text{K}^3`", ""
   ":math:`C_4`", "cp_mol_ig_comp_coeff_4", ":math:`\text{J/kmol}\cdotp\text{K}^4`", ""
   ":math:`C_5`", "cp_mol_ig_comp_coeff_5", ":math:`\text{J/kmol}\cdotp\text{K}^5`", ""
   ":math:`\Delta h_{\text{form, Liq}}`", "enth_mol_form_liq_comp_ref", "Base units", "Molar heat of formation at reference state"

.. note::
    This correlation uses the same parameters as the ideal liquid heat capacity.
    Units of molar heat of formation will be derived from the base units defined for the property package.

Ideal Liquid Molar Entropy
---------------------------

The correlation for the ideal liquid molar entropy is derived from the correlation for the molar heat capacity and is given below:

.. math:: s_{\text{liq}} - s_{\text{liq ref}} = C_1 \times ln(T/T_{ref}) + C_2 \times (T-T_{ref}) + \frac{C_3}{2} \times (T^2-T_{ref}^2) + \frac{C_4}{3} \times (T^3-T_{ref}^3) + \frac{C_5}{4} \times (T^4-T_{ref}^4) + s_{\text{form, Liq}}

Units are :math:`\text{J/kmol}\cdotp\text{K}`.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`C_1`", "cp_mol_ig_comp_coeff_1", ":math:`\text{J/kmol}\cdotp\text{K}`", ""
   ":math:`C_2`", "cp_mol_ig_comp_coeff_2", ":math:`\text{J/kmol}\cdotp\text{K}^2`", ""
   ":math:`C_3`", "cp_mol_ig_comp_coeff_3", ":math:`\text{J/kmol}\cdotp\text{K}^3`", ""
   ":math:`C_4`", "cp_mol_ig_comp_coeff_4", ":math:`\text{J/kmol}\cdotp\text{K}^4`", ""
   ":math:`C_5`", "cp_mol_ig_comp_coeff_5", ":math:`\text{J/kmol}\cdotp\text{K}^5`", ""
   ":math:`s_{\text{form, Liq}}`", "entr_mol_form_liq_comp_ref", "Base units", "Standard molar entropy of formation at reference state"

.. note::
    This correlation uses the same parameters as the ideal liquid heat capacity.
    Units of molar entropy of formation will be derived from the base units defined for the property package.

Liquid Molar Density
--------------------

Perry's Handbook uses the following correlations for liquid molar density:

.. math:: \rho_{liq} = \frac{C_1}{C_2^{1 + (1-\frac{T}{C_3})^{C_4}}} \space (1)
.. math:: \rho_{liq} = {C_1} + {C_2} \times {T} + {C_3} \times {T^2} + {C_4} \times {T^3} \space (2)

Units are :math:`\text{kmol/}\text{m}^3`.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units (form 1)", "Units (form 2)", "Description"

   "`eqn_type`", "dens_mol_comp_liq_coeff_eqn_type", "None", "None", "Flag in the set [1, 2]"
   ":math:`C_1`", "dens_mol_comp_liq_coeff_1", ":math:`\text{kmol/}\text{m}^3`", ":math:`\text{kmol/}\text{m}^3`", ""
   ":math:`C_2`", "dens_mol_comp_liq_coeff_2", "None", ":math:`\text{kmol/}\text{m}^3\text{K}`", ""
   ":math:`C_3`", "dens_mol_comp_liq_coeff_3", ":math:`\text{K}`", ":math:`\text{kmol/}\text{m}^3\text{K}^2`", ""
   ":math:`C_4`", "dens_mol_comp_liq_coeff_4", "None", ":math:`\text{kmol/}\text{m}^3\text{K}^3`", ""
    
.. note::
    When Perry's methods are used, an equation form for liquid molar density must be specified as an additional coefficient 'eqn_type.' This parameter
    may be either '1' or '2' to select an equation form. The second correlation form for liquid molar density is most often used for water or o-terphenyl
    (values exist in Perry's Handbook). Fitted coefficients must be entered with the correct Pyomo units to use a specific correlation form.
