Properties of Gases and Liquids
===============================

.. contents:: Contents 
    :depth: 2

Source
------

Methods for calculating pure component properties from:

The Properties of Gases & Liquids, 4th Edition
Reid, Prausnitz and Polling, 1987, McGraw-Hill

Ideal Gas Molar Heat Capacity (Constant Pressure)
-------------------------------------------------

Properties of Gases and Liquids uses the following correlation for the ideal gas molar heat capacity:

.. math:: c_{\text{p ig}, j} = A + B \times T + C \times T^2 + D \times T^3

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`A, B, C, D`", "cp_ig_coeff", "component, `['A', 'B', 'C', 'D']`", ""

Ideal Gas Molar Enthalpy
------------------------

The correlation for the ideal gas molar enthalpy is derived from the correlation for the molar heat capacity and is given below:

.. math:: h_{\text{ig}, j} - h_{\text{ig ref}, j} = A \times (T-T_{ref}) + \frac{B}{2} \times (T^2 - T_{ref}^2) + \frac{C}{3} \times (T^3 - T_{ref}^3) + \frac{D}{4} \times (T^4 - T_{ref}^4) + \Delta h_{\text{form, Vap}, j}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`A, B, C, D`", "cp_ig_coeff", "component, `['A', 'B', 'C', 'D']`", ""
   ":math:`\Delta h_{\text{form}}`", "enth_mol_form_ref", "phase, component", "Molar heat of formation at reference state"
   ":math:`T_{ref}`", "temperature_ref", "None", "Temperature at reference state"

.. note::
    This correlation uses the same parameters as the ideal gas heat capacity correlation.

Ideal Gas Molar Entropy
------------------------

The correlation for the ideal gas molar entropy is derived from the correlation for the molar heat capacity and is given below:

.. math:: s_{\text{ig}, j} = A \times ln(T) + B \times T + \frac{C}{2} \times T^2 + \frac{D}{3} \times T^3 + s_{\text{ref, Vap}, j}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`A, B, C, D`", "cp_ig_coeff", "component, `['A', 'B', 'C', 'D']`", ""
   ":math:`s_{\text{ref}}`", "entr_mol_ref", "phase, component", "Standard molar entropy at reference state"

.. note::
    This correlation uses the same parameters as the ideal gas heat capacity correlation .

Saturation (Vapor) Pressure
---------------------------

Properties of Gases and Liquids uses the following correlation to calculate the vapor pressure of a component:

.. math:: ln(\frac{P_{sat, j}}{P_{crit}}) \times (1-x) = A \times x + B \times x^1.5 + C \times x^3 + D \times x^6

where :math:`x = 1 - \frac{T}{T_{crit}}`.

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`A, B, C, D`", "pressure_sat_coeff", "component, `['A', 'B', 'C', 'D']`", ""
   ":math:`P_{crit}`", "temperature_crit", "None", "Critical pressure"
   ":math:`T_{crit}`", "temperature_crit", "None", "Critical temperature"

.. note::
    This correlation is only valid at temperatures **below** the critical temperature. Above this point, there is no real solution to the equation.
