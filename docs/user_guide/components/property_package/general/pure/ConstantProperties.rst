Constant Properties
=================================================

.. contents:: Contents 
    :depth: 2

Source
------

Methods for calculating pure component properties independent of temperature from:

Introduction to Chemical Engineering Thermodynamics, 8th Edition
J.M. Smith, Hendrick Van Ness, Michael Abbott, and Mark Swihart, 2018, McGraw-Hill

Ideal Liquid Molar Heat Capacity (Constant Pressure)
----------------------------------------------------

The ideal liquid molar heat capacity is defined as follows:

.. math:: c_{\text{p liq}} = C_1

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`C_1`", "cp_mol_ig_comp_coeff", "Units are defined based on the user's input", ""


Ideal Liquid Molar Enthalpy
---------------------------

The equation for the ideal liquid molar enthalpy is given below:

.. math:: h_{\text{liq}} - h_{\text{liq ref}} = c_{\text{p liq}} \times (T-T_{ref}) + \Delta h_{\text{form, liq}}

Units are defined based on the user's input.

.. note::
    This correlation uses the same parameters as the ideal liquid heat capacity.
    Units of molar heat of formation will be derived from the base units defined for the property package.


Ideal Liquid Molar Entropy
---------------------------

The correlation for the ideal liquid molar entropy is derived from the correlation for the molar heat capacity and is given below:

.. math:: s_{\text{liq}} - s_{\text{liq ref}} = c_{\text{p liq}} \times ln(T/T_{ref}) + s_{\text{form, Liq}}

Units are defined based on the user's input.

.. note::
    This correlation uses the same parameters as the ideal liquid heat capacity.
    Units of molar entropy of formation will be derived from the base units defined for the property package.


Liquid Molar Density
--------------------

The liquid molar density is defined as follows:

.. math:: \rho_{liq} = C_1

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`C_1`", "dens_mol_liq_comp_coeff", "Units are defined based on the user's input", ""


Ideal Gas Molar Heat Capacity (Constant Pressure)
----------------------------------------------------

The ideal gas molar heat capacity is defined as follows:

.. math:: c_{\text{p ig}} = C_1

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Units", "Description"

   ":math:`C_1`", "cp_mol_ig_comp_coeff", "Units are defined based on the user's input", ""


Ideal Gas Molar Enthalpy
---------------------------

The equation for the ideal gas molar enthalpy is given below:

.. math:: h_{\text{ig}} - h_{\text{ig ref}} = c_{\text{p ig}} \times (T-T_{ref}) + \Delta h_{\text{form, ig}}

Units are defined based on the user's input.

.. note::
    This correlation uses the same parameters as the ideal gas heat capacity.
    Units of molar heat of formation will be derived from the base units defined for the property package.


Ideal Gas Molar Entropy
---------------------------

The correlation for the ideal gas molar entropy is derived from the correlation for the molar heat capacity and is given below:

.. math:: s_{\text{ig}} - s_{\text{ig ref}} = c_{\text{p ig}} \times ln(T/T_{ref}) + s_{\text{form, ig}}

Units are defined based on the user's input.

.. note::
    This correlation uses the same parameters as the ideal gas heat capacity.
    Units of molar entropy of formation will be derived from the base units defined for the property package.
