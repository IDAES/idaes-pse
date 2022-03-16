Water Pipe Model
================

.. index::
  pair: idaes.power_generation.unit_models.waterpipe;WaterPipe

.. currentmodule:: idaes.power_generation.unit_models.waterpipe

Introduction
------------

The water pipe model is used to model a water or steam pipe connecting two units in a power plant.  It calculates the pressure change between the pipe inlet and outlet due to friction, gravity, and optional expansion or contraction at the end of the pipe.  The water pipe model does not provide the equations to calculate the heat loss.  However, the user can specify the heat duty if configuration variable “has_heat_transfer” is set to True.
When declaring the water pipe model, the user needs to provide typical configuration variables for a control volume, the base class the model is derived from, including “dynamic”, “has_holdup”, “has_heat_transfer”, “has_pressure_change”, etc.  While most of configuration variables have default values, the configuration variable for “property_package” has to be given as the IDAES property package for water implemented based on IAPWS water property table.  The user needs to set the variable “water_phase" either as “Liq” for liquid water or “Vap” for water vapor.  Currently the model does not support the pipe with two phase flow since the two-phase flow is usually unstable.  The user also needs to specify the configuration variable “contraction_expansion_at_end”. If there is a contraction at the end of the pipe, the value for the variable should be “contraction”.  If there is an expansion at the end, the value should be “expansion”. The value of “None” is used if there is no contraction or expansion at the end of the pipe.


Model inputs (variable name):

* number of pipes (number_of_pipes)
* tube dimensions (length, inner diameter, elevation change) (length, diameter, elevation_change)
* water/steam flow rate and states at inlet (flow_mol, enth_mol, pressure)
* pressure drop correction factor (fcorrection_dp)
* heat duty (usually fixed equal to 0)
* if expansion at the end of pipe is True, user needs to specify area ratio at the end (area_ratio), which is the area after contraction or expansion divided by the cross sectional area of the pipe)

Model Outputs:

* pressure drop (deltaP)
* water/steam flow rate and states at outlet (flow_mol, enth_mol, pressure)


Degrees of Freedom
------------------

As mentioned above, the waterpipe model includes rigorous pressure drop, therefore, detailed pipe dimensions are required. Aside from the inlet conditions and tube dimensions, the waterpipe model usually has one degree
of freedom, the heat duty, which can be fixed for it to be fully specified. 

Variables
---------

=========================== ================== =========== =============================================================================
Variable                    Symbol             Index Sets  Doc
=========================== ================== =========== =============================================================================
heat_duty                   :math:`Q`          time        Heat transferred from flue gas to tube side fluid
deltaP                      :math:`\DeltaP`    time        Pressure drop
=========================== ================== =========== =============================================================================

Constraints
-----------

The main constraints in this model are used to compute the pressure drop. 
Three types of pressure changes are considered in the pipe model, including the pressure change due to friction along the pipe length, the pressure change due to gravity if there is an elevation change from the inlet to the outlet, and the pressure change due to the contraction or expansion at the end.


Pressure drop:

.. math::
  \Delta P = \Delta P_{friction} + \Delta P_{gravity} + \Delta P_{contraction}

.. math:: \Delta P_{tube friction} = f(diameter \rho, V, number of pipes, length, Re, fcorrection_{dp})

.. math::
  deltaP_{gravity} = f(\rho, acceleration gravity, elevation_{change})

.. math::
  deltaP_{contraction} = f(\rho, V, K_{loss})

Reynolds number:

.. math:: Re_{tube} = \frac{tube_{di} V \rho}{\mu}


where:

* fcorrection_dp: correction factor if the pipe is not smooth
* diameter: inner diameter (m)
* :math:`\rho` : density (kg/m^3)
* Re : Reynolds number
* V: fluid velocity (m/s)
* Kloss: loss coefficient due to the contraction (function of the area ratio)

Dynamic Model
-------------

The dynamic model version of the steam heater model can be constructed by selecting dynamic=True. 
If dynamic = True, the material accumulation holdups are constructed. While, the metal energy holdup is not considered.