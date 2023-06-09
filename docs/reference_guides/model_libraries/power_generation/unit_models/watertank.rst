Water Tank
==========

The IDAES water tank model represents a unit operation for storing water. The water tank model 
supports several shapes including rectangular, vertical and horizontal cylindrical.

Model Structure
---------------

The water tank unit model consists of a single ``ControlVolume0D`` (named ``control_volume``) with one 
Inlet Port (named ``inlet``) and one Outlet Port (named ``outlet``).


Construction Arguments
----------------------

Similar to other IDAES unit models, the water tank has the following construction arguments:

========================= =================
Argument                  Default Value
========================= =================
dynamic                   False
include_holdup            False
material_balance_type     MaterialBalanceType.componentPhase
energy_balance_type       EnergyBalanceType.enthalpyTotal
momentum_balance_type     MomentumBalanceType.pressureTotal
has_heat_transfer         True
has_pressure_change       True
property_package          Parent value
property_package_args     \--
========================= =================

Additionally, the water tank model has one specific construction argument to declare the tank shape:

* tank_type: configuration argument to define the shape of the tank to be modeled, and accordingly calculate the volume of the filled level. Currently, the supported values are: simple_tank, rectangular_tank, vertical_cylindrical_tank, and horizontal_cylindrical_tank. Being simple_tank the default value.

Variables
---------

The following variables are added to the model independently of the tank type selected:

Model Inputs (variable name) - symbol:

* water inlet (inlet: flow_mol, enth_mol, pressure)
* tank level (tank_level) - :math:`l`
* heat duty (heat_duty) - :math:`Q`

Model Outputs (variable name):

* water outlet (outlet: flow_mol, enth_mol, pressure)
* pressure drop (deltaP) - :math:`\Delta P`

Additionally, some variables are added to the model based on the tank type as indicated below:

=========================== ==================================================
tank_type                   Variable added  
=========================== ==================================================
simple_tank                 tank_cross_sect_area - :math:`A_{c}`
rectangular_tank            tank_width - :math:`W`, tank_length - :math:`L`
vertical_cylindrical_tank   tank_diameter - :math:`d`
horizontal_cylindrical_tank tank_diameter - :math:`d`, tank_length - :math:`L`
=========================== ==================================================


Constraints
-----------

The main assumptions used in the water tank unit model are:

1) Heat loss is a variable given by the user (zero heat loss can be specified if adiabatic)
2) Calculate pressure change due to gravity based on water level
3) Water level is either fixed for steady-state model or calculated for dynamic model
4) Assume enthalpy_in == enthalpy_out + heat loss

In addition to the constraints written by the control volume, the water tank model adds two constraints
for the pressure drop and the volume of the liquid level in the unit. 

Pressure drop constraint:

.. math::
  \Delta P = \Delta P_{gravity} = \rho_{liq} * g * l

Volume of the liquid constraint:

1) for simple_tank, rectangular_tank and vertical_cylindrical_tank:

.. math::
  V_{liq} = l * A_{c}

2) for horizontal_cylindrical_tank:

.. math::
  V_{liq} = L * A_{t}

where:

* :math:`\rho_{liq}`: liquid density
* :math:`l`: level filled by liquid in the unit
* :math:`g`: acceleration gravity
* :math:`L`: tank length
* :math:`A_{c}`: cross sectional area of the tank, which for the simple_tank is an input variable, while for rectangular_tank and vertical_cylindrical_tank is an expression calculated by the model
* :math:`A_{t}`: area of the circular segment covered by the liquid level at one end of the tank. This is an expression calculated by the model and is only valid for the horizontal_cylindrical_tank

The following expressions were used to calculate the tank cross sectional area, and tank area:

* for rectangular_tank: :math:`A_{c} = W * L`

* for vertical_cylindrical_tank: :math:`A_{c} = \pi * r^{2}`, tank_radius (r) is an expression calculated by the model

* for horizontal_cylindrical_tank: :math:`A_{t} = cos^{-1}(1-l/r)*r^{2}-(r-l)*(2rl-l^{2})^{0.5}`

Degrees of Freedom
------------------

The degrees of freedom depend on the tank type as the dimension variables are different for each type, 
but once the dimensions for a specific tank type have been fixed, the model generally has 3-5 degrees of freedom: 
the inlet state (flow_mol, enth_mol, and pressure), the heat duty whether the config argument has_pressure_change is set to True, 
and the tank level whether for steady state simulations


Dynamic Model
-------------

The dynamic model version of the tank model can be constructed by selecting dynamic=True. 
If dynamic = True, material accumulation, energy accumulation, and tank level must be calculated. Therefore, a dynamic initialization method has been developed `set_initial_conditions` to initialize the holdup terms.
