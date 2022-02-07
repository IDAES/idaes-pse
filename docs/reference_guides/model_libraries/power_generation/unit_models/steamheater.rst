Steam Heater Model
==================

.. index::
  pair: idaes.power_generation.unit_models.steamheater;SteamHeater

.. currentmodule:: idaes.power_generation.unit_models.steamheater

Introduction
------------

The steam heater model consists of a heater model with rigorous heat transfer calculations on the tube side, while the heat duty from fire side is either fixed or provided by the boiler fire side model. 
The model is usually coupled with the IDAES 1-D fire-side model to solve the wall temperatures and heat transfer rate.
If coupled with the fire side model, this model is similar to the water_wall section model. The sum of the net radiation and convective heat fluxes (:math:`q_{rad}^{fire}` and :math:`q_{conv}^{fire}`) at the slag outer layer 
is an output of the fire-side model and an input of the steam heater model (the fluid-side model). While, the temperature of the outer slag layer :math:`T_{w,slag}` is an output of the fluid-side model 
and an input (boundary condition) of the fire-side model. The heat conduction through the slag and tube layers is a part of the fluid-side model. At a steady state, the amount of the heat transferred at the outer slag surface 
(:math:`q_{rad}^{fire}` and :math:`q_{conv}^{fire}`) is equal to the heat conducted through the slag and tube layers, which is equal to the heat convected to the fluid :math:`q_{conv}^{fluid}`.

Model inputs (variable name):

* number of tubes (number_tubes)
* heat duty from fire-side model (sum of net radiation and convection) (heat_fireside)
* tube dimensions (length, inside diameter and thickness) (tube_length, tube_diameter, tube_thickness)
* fin dimension of membrane wall (width and thickness) (fin_length, fin_thickness)
* slag layer thickness (slag_thickness)
* water/steam flow rate and states at inlet (flow_mol, enth_mol, pressure)
* properties of slag and tube metal (thermal conductivity, heat capacity, density) (therm_cond_slag, therm_cond_metal, dens_metal, dens_slag)
* pressure drop correction factor (fcorrection_dp)

Model Outputs:

* temperatures of tube metal at inner wetted surface and at center of the tube thickness (temp_tube_boundary, temp_tube_center)
* temperatures of slag layer at outer surface and at the center of the slag layer (temp_slag_boundary, temp_slag_center)
* pressure drop through each section and heat added to the tube (deltaP and heat_duty, respectively)
* water/steam flow rate and states at outlet (flow_mol, enth_mol, pressure)


Degrees of Freedom
------------------

As mentioned above, the steam heater model includes rigorous heat transfer, therefore, detailed tube and unit dimensions are required. Aside from the inlet conditions and tube dimensions, the steam heater model usually has two degrees
of freedom, heat flux from fire side and slag thickness, which can be fixed for it to be fully specified. 

Variables
---------

=========================== ================== =========== =============================================================================
Variable                    Symbol             Index Sets  Doc
=========================== ================== =========== =============================================================================
heat_duty                   :math:`Q`          time        Heat transferred from flue gas to tube side fluid
hconv                       :math:`h_{conv}`   time        Overal convective heat transfer coefficient
temp_slag_boundary          :math:`T_{w,slag}` time        Temperature of the slag
projected_area              :math:`A`          None        Heat transfer area (total projected area based on tube shape)
=========================== ================== =========== =============================================================================

Constraints
-----------

The main constraints here show the heat flux, convective heat transfer model, and pressure drop. 
This model calculates the slag temperature, slag center temperature, tube boundary temperature, tube center temperatures, and heat flux from fire side to the water/steam side. 


Heat flux equation:

.. math::
  heat_{flux} = Q*pitch/(projected_{area}*perimeter\_slag)


Temperature of slag:

.. math::
  T_{w,slag} - T_{c,slag} = heat_{flux} * slag_{resistance}


Heat flux interface equation:

.. math::
  heat_{flux\_int} * (slag_{resistance} + metal_{resistance}) = (T_{c,slag} - T_{c,tube})


Convective heat flux eqn at tube boundary:

.. math::
  heat_{flux\_conv} * fshape_{conv} * tube_{perimeter} = pitch * h_{conv} * (T_{w,tube} - T_{fluid,in})


Tube boundary wall temperature: 

.. math::
  heat_{flux\_conv} * metal_{resistance} * tube_{perimeter} = interface_{perimeter} * (T_{c,tube} - T_{w,tube})


Heat equation:

.. math:: heat_{duty} = number_{tubes} * heat_{flux\_conv} * tube_{length} * tube_{perimeter}


Pressure drop:

.. math::
  \Delta P = \Delta P_{friction} + \Delta P_{gravity}

Convective heat transfer:

.. math:: h_{conv} = f (tube_{diameter}, N_Re, N_Pr, k)


Prandtl number:

.. math:: Pr_{tube} = \frac{Cp  \mu}{ k  Mw}

Reynolds number:

.. math:: Re_{tube} = \frac{tube_{di} V \rho}{\mu}


where:

* hconv : convective heat transfer coefficient tube side (fluid water/steam) (W / m2 / K)
* projected_area : total projected wall area (m2)
* Pr : Prandtl number
* Re : Reynolds number
* V: fluid velocity (m/s)
* k : thermal conductivity of the fluid (W / m / K)
* MW: molecular weigth of water/steam (kmol/kg)

Note that at the flowsheet level first waterwall section is connected to the economizer, arcs connecting section 2 to n-1 have to be constructed by the user, and the outlet of section n is connected to the drum model or superheater (subcritical and supercritical plant, respectively)

Dynamic Model
-------------

The dynamic model version of the steam heater model can be constructed by selecting dynamic=True. 
If dynamic = True, the energy accumulation of slag and metal, material accumulation holdups are constructed. Therefore, a dynamic initialization method has been developed `set_initial_conditions` to initialize the holdup terms.
