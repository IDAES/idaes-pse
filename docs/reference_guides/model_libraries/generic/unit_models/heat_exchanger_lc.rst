Heat Exchanger (Lumped Capacitance)
===================================

.. index::
   pair: idaes.models.unit_models.heat_exchanger_lc;HeatExchangerLumpedCapacitance

.. currentmodule:: idaes.models.unit_models.heat_exchanger_lc

Author: `Rusty Gentile <https://github.com/rustygentile>`_

The ``HeatExchangerLumpedCapacitance`` model can be imported from :code:`idaes.models.unit_models`.
This model is an extension of ``idaes.models.unit_models.heat_exchanger``,
with wall temperature and heat holdup terms added for use in transient simulations. Using the electric
circuit analogy, heat stored in the wall material is similar to charge in a capacitor.

Degrees of Freedom
------------------

``HeatExchangerLumpedCapacitance`` has three degrees of freedom. For most scenarios, these will be the wall
temperature and exit temperatures on the hot and cold sides. The constraints for this model
should be similar to those of the standard
:ref:`Heat Exchanger <reference_guides/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>`
model. In lieu of "heat_transfer_coefficient", however, the following should be fixed:

* ua_cold_side
* ua_hot_side

The user may also add constraints to make these functions of fluid state.

Model Structure
---------------

Structure is derived from the base 0D :ref:`Heat Exchanger <reference_guides/model_libraries/generic/unit_models/heat_exchanger:Model Structure>`
model. Aside from adding new variables and constraints, the structure of this model is unchanged.

Variables
---------

All variables from the base 0D :ref:`Heat Exchanger <reference_guides/model_libraries/generic/unit_models/heat_exchanger:Variables>`
model are included. This model contains the following variables in addition. Each of these is indexed in the time domain.

=========================== ============================= ==========================================================================================
Variable                    Symbol                        Doc
=========================== ============================= ==========================================================================================
temperature_wall            :math:`T_{wall}`              Average wall temperature
dT_wall_dt                  :math:`\frac{dT_{wall}}{dt}`  Derivative of wall temperature with respect to time
ua_cold_side                :math:`UA_{cold}`             Overall heat transfer coefficient from the cold side
ua_hot_side                 :math:`UA_{hot}`              Overall heat transfer coefficient from the hot side
ua_hot_side_to_wall         :math:`UA_{hot, wall}`        Overall heat transfer coefficient between the hot side and the center of the wall material
hot_side_heat               :math:`Q_{hot}`               Heat entering the hot side control volume (value is typically negative)
cold_side_heat              :math:`Q_{cold}`              Heat entering the cold side control volume (value is typically positive)
=========================== ============================= ==========================================================================================

Note: "ua" variables are of the form :math:`hconv \times area`. So the units of measure of these should be :math:`\frac{energy}{time \times temperature}`

Parameters
----------

The model uses the following parameters. These are not indexed as we assume they are constant.

=========================== ====================== ==========================================================================================
Parameter                    Symbol                 Doc
=========================== ====================== ==========================================================================================
heat_capacity_wall          :math:`C_{wall}`       Total heat capacity of heat exchanger material
thermal_resistance_wall     :math:`R_{wall}`       Total thermal resistance of heat exchanger material
thermal_fouling_cold_side   :math:`R_{foul, cold}` Total thermal resistance due to fouling on the cold side
thermal_fouling_hot_side    :math:`R_{foul, hot}`  Total thermal resistance due to fouling on the hot side
=========================== ====================== ==========================================================================================

Constraints
-----------

Note: all constraints from the base 0D :ref:`Heat Exchanger <reference_guides/model_libraries/generic/unit_models/heat_exchanger:Constraints>`
model are also included.

Total overall heat transfer coefficient:

.. math::
  \frac{1}{UA} = \frac{1}{UA_{hot}} + R_{foul, hot} + R_{wall} + R_{foul, cold} + \frac{1}{UA_{cold}}

Overall heat transfer coefficient from the hot side to the center of the wall:

.. math::
  \frac{1}{UA_{hot, wall}} = \frac{1}{UA_{hot}} + R_{foul, hot} + \frac{1}{2}R_{wall}

Wall temperature equation:

.. math::
  T_{wall} = \frac{1}{2}(T_{hot, in} + T_{hot, out}) + \frac{Q_{hot}}{UA_{hot, wall}}

Dynamic heat balance (if ``dynamic_heat_balance`` is set to ``True``):

.. math::
  Q_{hot} + Q_{cold} + \frac{dT_{wall}}{dt}C_{wall} = 0

Standard heat balance (if ``dynamic_heat_balance`` is set to ``False``):

.. math::
  Q_{hot} + Q_{cold} = 0

Class Documentation
-------------------

.. autoclass:: HeatExchangerLumpedCapacitance
   :members:

.. autoclass:: HeatExchangerLumpedCapacitanceData
   :members:
