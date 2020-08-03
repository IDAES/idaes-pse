State Definition
================

Defining State Variables
------------------------

An important part of defining a set of property calculations is choosing the set of variables which will describe the state of the material. The set of state variables needs to include information on the extensive flow, composition and thermodynamic state of the material. However, there are many ways in which this information can be described, and the best choice of state variables depends on many factors.

Within the IDAES Generic Property Package Framework, the definition of state variables is done using sub-modules which create the necessary variables supporting information for the property package. A state definition sub-module may define any set of state variables the user feel appropriate, but must define the following components as either state variables or functions of the state variables:

* `temperature` (must be a Pyomo Var)
* `pressure`
* `mole_frac_phase_comp`
* `phase_frac`

The IDAES Generic Property Package Framework has a library of prebuilt state definition sub-modules for users to use which are listed below.

State Definition Libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 1

    state/FTPx
    state/FcTP
    state/FPhx
    state/FcPh
    state/FpcTP

Setting Bounds on State Variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For optimization applications, it is important to specify a good initial guess and bounds on the state variables in order to improve the robustness of the problem. Further, due to the empirical nature of most thermophysical correlations these correlations are only valid in specific range of states. Users should set the `state_bounds` configuration argument to define the bounds on the state variables of their property package.

The `state_bounds` configuration argument should be a `dict` where the keys are the names of the state variables (using the standard naming convention) and the values should be a tuple with the form `(lower, nominal, upper, units)`. The `lower` and `upper` values are used to set the lower and upper bounds respectively, whilst the `nominal` value is used to set the initial value for the state variable. The `units` value is optional, and is used to specify the units of measurement for the values provided, which will be used to automatically convert these values to the base set of units defined for the property package if required. If the `units` value is omitted, it is assumed that the values provided are in the base unit set for the property package.

.. note::
    Some state definitions allow for setting on additional variables beyond the chosen state variables (temperature is a common example). See the documentation for your state definition for more information on what bounds can be set using the `state_bounds` argument.

Reference State
---------------

Many thermophysical properties are relative quantities, and require the definition of a thermodynamic reference state. Whilst some simpler models and correlations forego this or define the reference state implicitly, the IDAES Generic Property Package requires the user to specify the thermodynamic reference state (even if it is not used explicitly).

As such, users must provide the following two configuration arguments:

* `pressure_ref` - pressure at reference state
* `temperature_ref` - temperature at reference state

