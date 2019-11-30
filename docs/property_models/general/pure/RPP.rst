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

Parameters should be provided as a single `Param` named `cp_ig_coeff` indexed by component and parameter names (e.g. `[component, "A"]`).

Ideal Gas Molar Enthalpy
------------------------

The correlation for the ideal gas molar enthalpy is derived from the correaltion for the molar heat capacity and is given below:

.. math:: h_{\text{ig}, j} - h_{\text{ig ref}, j} = A \times (T-T_{ref}) + \frac{B}{2} \times (T^2 - T_{ref}^2) + \frac{C}{3} \times (T^3 - T_{ref}^3) + \frac{D}{4} \times (T^4 - T_{ref}^4)

This correlation uses the same paramters as for the ideal gas heat capacity, which should be provided as a single `Param` named `cp_ig_coeff` indexed by component and parameter names (e.g. `[component, "A"]`). Additionally, the ideal gas molar enthalpy requires a reference temperature (`temperature_ref`) which must also be defined in the parameter block.

Ideal Gas Molar Entrolpy
------------------------

The correlation for the ideal gas molar entropy is derived from the correaltion for the molar heat capacity and is given below:

.. math:: s_{\text{ig}, j} = A \times ln(T) + B \times T + \frac{C}{2} \times T^2 + \frac{D}{3} \times T^3

This correlation uses the same paramters as for the ideal gas heat capacity, which should be provided as a single `Param` named `cp_ig_coeff` indexed by component and parameter names (e.g. `[component, "A"]`).

Saturation (Vapor) Pressure
---------------------------

Properties of Gases and Liquids uses the following correlation to calculate the vapor pressure of a component:

.. math:: ln(\frac{P_{sat, j}}{P_{crit}}) \times (1-x) = A \times x + B \times x^1.5 + C \times x^3 + D \times x^6

where :math:`x = 1 - \frac{T}{T_{crit}}`. This correlation requires the critical pressure and temeprature of each component be defined (`pressure_crit` and `temperature_crit`) along with a single `Param` named `pressure_sat_coeff` indexed by component and parameter names (e.g. `[component, "A"]`).

.. note::
    This correlation is only valid at temperatrues **below** the critical temperture. Above this point, there is no real solution to the equation.

