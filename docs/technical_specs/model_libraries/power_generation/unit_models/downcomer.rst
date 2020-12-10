Downcomer Model
===============

.. index::
  pair: idaes.power_generation.unit_models.downcomer;Downcomer

.. currentmodule:: idaes.power_generation.unit_models.downcomer

Introduction
------------
The Downcomer model consists of a simple pipe model (or a set of pipes) where the inlet stream is the Drum outlet and the outlet stream connects with the WaterWall (section 1).
The model simply calculates the pressure change due to friction and gravity, which involves the calculation of fluid velocity, Reynolds number, and friction factor (using Darcy’s correlation). 

Property package: This model requires the Helmholtz EoS (IAPWS95) property package with the mixed phase option, therefore, the phase equilibrium calculations are handled by the property package.

Model inputs (variable name):

* inlet stream (flow_mol, enth_mol, pressure)
* number of downcomer pipes (number_downcomers) same as Drum model
* height of the tubes (height)
* inner diameter of the tubes (diameter)
* heat duty (heat_duty), heat_duty = 0 if adiabatic

Model Outputs:

* outlet stream (flow_mol, enth_mol, pressure)
* pressure change (deltaP) due to gravity and friction


Degrees of Freedom
------------------

By specifying the inlet conditions and downcomer dimensions, the model will be fully specified. Things that are frequently fixed are:

* inlet state vars (generally flow_mol, enth_mol, pressure)
* heat_duty to the downcomer (if applicable)
* number_downcomers, height, diameter

Variables
---------

=========================== ================== =========== =============================================================================
Variable                    Symbol             Index Sets  Doc
=========================== ================== =========== =============================================================================
heat_duty                   :math:`Q`          time        Heat transferred from flue gas to tube side fluid
deltaP                      :math:`deltaP`     time        Pressure change in the unit
=========================== ================== =========== =============================================================================


Constraints
-----------

The main constraints in the model calculate the pressure drop, which is given by deltaP_friction and deltaP_gravity. 

Pressure drop:

.. math::
  \Delta P = \Delta P_{friction} + \Delta P_{gravity}

Friction:

.. math::
  \Delta P_{friction} = f (Friction_{factor}, \rho_{liquid}, V, diameter)
  
Friction factor (Darcy's correlation'):

.. math::
  Friction_{factor} = \frac{0.3164}{Re^{0.25}}

deltaP gravity:

.. math::
  \Delta P_{gravity} = \rho_{liquid} * acceleration_{gravity} * height

Reynolds number:

.. math:: Re = \frac{tube_{di} V \rho}{\mu}


where:

* Re : Reynolds number (liquid)
* V: fluid velocity (m/s, liquid)
* :math:`\rho_{liquid}`: mass density of liquid (kg/m3)
* :math:`\mu` : viscocity (kg/m/s)

Dynamic Model
-------------
The downcomer dynamic model can be constructed by selecting dynamic=True and hold_up=True. 
If dynamic and holdup = True, the energy accumulation and material accumulation variables are constructed, 
which are the derivatives of the corresponding holdup terms with respect to time and are included in the material and energy conservation equations. 
A dynamic model also requires the specification of initial conditions related to the accumulation variables. 
The user needs to provide the initial values for the accumulation variables at all time points and fix the initial conditions to solve the dynamic problem. 
Therefore, a dynamic initialization method has been developed `set_initial_conditions` to initialize the values of the time-indexed accumulation variables to zero 
and fix the variables at the first time point to zero.
