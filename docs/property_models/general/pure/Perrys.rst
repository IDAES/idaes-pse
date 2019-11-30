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

Perry's Handbook uses the follwoing correlation for ideal liquid molar heat capacity:

.. math:: \frac{c_{\text{p liq}, j}}{1000} = C_1 + C_2 \times T + C_3 \times T^2 + C_4 \times T^3 + C_5 \times T^4

Parameters should be provided as a single `Param` named `cp_liq_coeff` indexed by component and parameter number (e.g. `[component, 1]`).

Ideal Liquid Molar Enthalpy
---------------------------

The correlation for the ideal liquid molar enthalpy is derived from the correaltion for the molar heat capacity and is given below:

.. math:: \frac{h_{\text{liq}, j} - h_{\text{liq ref}, j}}{1000} = C_1 \times (T-T_{ref}) + \frac{C_2}{2} \times (T^2 - T_{ref}^2) + \frac{C_3}{3} \times (T^3 - T_{ref}^3) + \frac{C_4}{4} \times (T^4 - T_{ref}^4) + \frac{C_5}{5} \times (T^5 - T_{ref}^5)

This correlation uses the same paramters as for the ideal liquid heat capacity, which should be provided as a single `Param` named `cp_liq_coeff` indexed by component and parameter number (e.g. `[component, 1]`). Additionally, the ideal liquid molar enthalpy requires a reference temperature (`temperature_ref`) which must also be defined in the parameter block.

Ideal Liquid Molar Entrolpy
---------------------------

The correlation for the ideal liquid molar entropy is derived from the correaltion for the molar heat capacity and is given below:

.. math:: s_{\text{liq}, j} = C_1 \times ln(T) + C_2 \times T + \frac{C_3}{2} \times T^2 + \frac{C_4}{3} \times T^3 + \frac{C_5}{4} \times T^4

This correlation uses the same paramters as for the ideal liquid heat capacity, which should be provided as a single `Param` named `cp_liq_coeff` indexed by component and parameter number (e.g. `[component, 1]`).


