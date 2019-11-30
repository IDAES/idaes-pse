NIST Webbook
============

.. contents:: Contents 
    :depth: 2

Source
------

Pure component properties as used by the NIST WebBook

`<https://webbook.nist.gov/chemistry/>`_ Retrieved: September 13th, 2019

Ideal Gas Molar Heat Capacity (Constant Pressure)
-------------------------------------------------

NIST uses the Shomate equation for the ideal gas molar heat capacity, which is shown below:

.. math:: c_{\text{p ig}, j} = A + B \times t + C \times t^2 + D \times t^3 + \frac{E}{t^2}

where :math:`t = \frac{T}{1000}`. Parameters should be provided as a single `Param` named `cp_ig_coeff` indexed by component and parameter names (e.g. `[component, "A"]`).

Ideal Gas Molar Enthalpy
------------------------

The correlation for the ideal gas molar enthalpy is derived from the correaltion for the molar heat capacity and is given below:

.. math:: \frac{h_{\text{ig}, j} - h_{\text{ig ref}, j}}{1000} = A \times (t-t_{ref}) + \frac{B}{2} \times (t^2 - t_{ref}^2) + \frac{C}{3} \times (t^3 - t_{ref}^3) + \frac{D}{4} \times (t^4 - t_{ref}^4) + E \times (\frac{1}{t} - \frac{1}{t_{ref}})

where :math:`t = \frac{T}{1000}` and :math:`t_{ref} = \frac{T_{ref}}{1000}`. This correlation uses the same paramters as for the ideal gas heat capacity, which should be provided as a single `Param` named `cp_ig_coeff` indexed by component and parameter names (e.g. `[component, "A"]`). Additionally, the ideal gas molar enthalpy requires a reference temperature (`temperature_ref`) which must also be defined in the parameter block.

.. note::
    The NIST Webbook includes additional parameters `F` and `H` in their correlation, however these cancel out when writing :math:`h-h_{ref}`.

Ideal Gas Molar Entrolpy
------------------------

The correlation for the ideal gas molar entropy is derived from the correaltion for the molar heat capacity and is given below:

.. math:: s_{\text{ig}, j} = A \times ln(t) + B \times t + \frac{C}{2} \times t^2 + \frac{D}{3} \times t^3 + \frac{E}{2 \times t^2} + G 

This correlation uses the same paramters as for the ideal gas heat capacity, which should be provided as a single `Param` named `cp_ig_coeff` indexed by component and parameter names (e.g. `[component, "A"]`).

Saturation (Vapor) Pressure
---------------------------

NIST uses the Antoine equation to calculate the vapor pressure of a component, which is given below:

.. math:: log_{10}(P_{sat, j}) = A - \frac{B}{T+C}

This correlation requires a single `Param` named `antoine_coeff` indexed by component and parameter names (e.g. `[component, "A"]`).

.. note::
    The Antoine equation is generally written with saturation pressure expressed in bars. The units of the correlation can be converted to Pascals by adding 5 to :math:`A`.

