Global Options
==============

The IDAES Generic Properties Framework also includes a number of general options for further customizing how certain properties are calculated. These options are listed below. Note that these are generally considered to be global options and will affect all phases and components in the system. Global configuration arguments can be declared at the top level of the property package configuration dict.

.. contents:: Contents 
    :depth: 2

Enthalpy of Formation
---------------------

In most process applications, it is necessary to include calculation of the enthalpy of formation for each component in the system in order to account for latent heat due to chemical reactions or phase changes. However, in certain circumstances, such as systems with no reactions or phase equilibria or cases where the user wishes to specify latent heats as separate parameters, it may be desirable to exclude the heat of formation from the specific enthalpy calculations.

Users can define whether or not the heat of formation should be included in the calculation of the specific enthalpy using the `include_enthalpy_of_formation` configuration argument. If this is set to `True` (the default behavior) then heat of formation will be included in the calculation of the specific enthalpy, if set to `False` then heat of formation will be excluded.



