Perry's Chemical Engineer's Handbook
====================================

.. contents:: Contents 
    :depth: 2

Source
------

Methods for calculating pure component properties from:

Perry's Chemical Engineers' Handbook, 7th Edition
Perry, Green, Maloney, 1997, McGraw-Hill

Ideal Liquid Molar Heat Capacity (Constant Pressure)
----------------------------------------------------

Perry's Handbook uses the following correlation for ideal liquid molar heat capacity:

.. math:: \frac{c_{\text{p liq}, j}}{1000} = C_1 + C_2 \times T + C_3 \times T^2 + C_4 \times T^3 + C_5 \times T^4

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`C_1, C_2, C_3, C_4, C_5`", "cp_mol_liq_comp_coeff", "component, `[1, 2, 3, 4, 5]`", ""

Ideal Liquid Molar Enthalpy
---------------------------

The correlation for the ideal liquid molar enthalpy is derived from the correlation for the molar heat capacity and is given below:

.. math:: \frac{h_{\text{liq}, j} - h_{\text{liq ref}, j}}{1000} = C_1 \times (T-T_{ref}) + \frac{C_2}{2} \times (T^2 - T_{ref}^2) + \frac{C_3}{3} \times (T^3 - T_{ref}^3) + \frac{C_4}{4} \times (T^4 - T_{ref}^4) + \frac{C_5}{5} \times (T^5 - T_{ref}^5) + \Delta h_{\text{form, Liq}, j}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`C_1, C_2, C_3, C_4, C_5`", "cp_mol_liq_comp_coeff", "component, `[1, 2, 3, 4, 5]`", ""
   ":math:`\Delta h_{\text{form}}`", "enth_mol_form_phase_comp_ref", "phase, component", "Molar heat of formation at reference state"
   ":math:`T_{ref}`", "temperature_ref", "None", "Temperature at reference state"

.. note::
    This correlation uses the same parameters as for the ideal gas heat capacity.

Ideal Liquid Molar Entropy
---------------------------

The correlation for the ideal liquid molar entropy is derived from the correlation for the molar heat capacity and is given below:

.. math:: s_{\text{liq}, j} = C_1 \times ln(T) + C_2 \times T + \frac{C_3}{2} \times T^2 + \frac{C_4}{3} \times T^3 + \frac{C_5}{4} \times T^4 + s_{\text{ref, Vap}, j}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`C_1, C_2, C_3, C_4, C_5`", "cp_mol_liq_comp_coeff", "component, `[1, 2, 3, 4, 5]`", ""
   ":math:`s_{\text{ref}}`", "entr_mol_phase_comp_ref", "phase, component", "Standard molar entropy at reference state"

.. note::
    This correlation uses the same parameters as for the ideal gas heat capacity.

Liquid Molar Density
--------------------

Perry's Handbook uses the following correlation for liquid molar density:

.. math:: \rho_{liq, j} = \frac{C_1}{C_2^{1 + (1-\frac{T}{C_3})^{C_4}}}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`C_1, C_2, C_3, C_4`", "dens_mol_comp_liq_coeff", "component, `[1, 2, 3, 4]`", ""

.. note::
    Currently, only the most common correlation form from Perry's Handbook is implemented. Some components use different forms which are not yet supported.
