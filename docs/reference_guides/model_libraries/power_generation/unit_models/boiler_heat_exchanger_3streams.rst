Heat Exchanger With Three Streams
=================================

.. index::
   pair: idaes.power_generation.unit_models.heat_exchanger_3streams; HeatExchangerWith3Streams

.. currentmodule:: idaes.power_generation.unit_models.heat_exchanger_3streams

The HeatExchangerWith3Streams model consists of a heat exchanger with three inlets, `side_1` represents the hot stream, while `side_2` and `side_3` are cold streams.
This model is a simplified generic heat exchanger model with lumped UA (the product of the overall heat transfer coefficient and the heat transfer area).

In a power plant flowsheet this model is used to represent an air preheater unit. This is because modeling the Ljungström type preheater is quite challenging since it involves not only the hot and cold gas streams but also the energy stored in and relased from the metal parts.


Degrees of Freedom
------------------

Aside from the inlet conditions, a 3 inlet heat exchanger model usually has six degrees
of freedom, which must be fixed for it to be fully specified. Things that are
frequently fixed are two of:

* UA_side_2 - lumped overall heat transfer and heat transfer area of side 2
* UA_side_3 - lumped overall heat transfer and heat transfer area of side 3
* frac_heatloss - fraction of heat loss in the system
* deltaP_side_1 - pressure drop in side 1
* deltaP_side_2 - pressure drop in side 2
* deltaP_side_3 - pressure drop in side 3


Model Structure
---------------

The ``HeatExchangerWith3Streams`` model contains three ``ControlVolume0DBlock`` blocks. The
hot side is named ``side_1`` and two cold sides are named ``side_2`` and ``side_3``. These names are not configurable.
The sign convention is that duty is positive for heat flowing from the hot side to the cold
side.

The control volumes are configured the same as the ``ControlVolume0DBlock`` in the
:ref:`Heater model <reference_guides/model_libraries/generic/unit_models/heater:Heater>`.
The ``HeatExchangerWith3Streams`` model contains additional constraints that calculate the amount
of heat transferred from the hot side to the cold side.

The ``HeatExchangerWith3Streams`` has three inlet ports and three outlet ports. By default these are
``side_1_inlet``, ``side_2_inlet``, ``side_3_inlet``, ``side_1_outlet``, ``side_2_outlet``, ``side_3_outlet``.

Variables
---------

=========================== ================== =========== =============================================================================
Variable                    Symbol             Index Sets  Doc
=========================== ================== =========== =============================================================================
heat_duty                   :math:`Q`          time        Heat transferred (model includes 3 variables, one for each side)
UA                          :math:`UA`         None        lumped Heat transfer area and overall heat transfer coefficient
LMTD                        :math:`LMTD`       time        Log Mean Temperature difference, LMTD
=========================== ================== =========== =============================================================================


Constraints
-----------

The default constraints can be overridden by providing :ref:`alternative rules
<reference_guides/model_libraries/generic/unit_models/heat_exchanger:Callbacks>` for
the heat transfer equation, temperature difference, heat transfer coefficient, shell
and tube pressure drop. This section describes the default constraints.

Heat transfer from hot to cold sides:

.. math::
  Q_{side\_1} * (1-frac_{heat\_loss}) = Q_{side\_2} + Q_{side\_3}

.. math::
  Q_{side\_2} = UA_{side\_2}\Delta T_2

.. math::
  Q_{side\_3} = UA_{side\_3}\Delta T_3

Temperature difference is:

.. math::
  \Delta T = \frac{\Delta T_1 - \Delta T_2}{\log_e\left(\frac{\Delta T_1}{\Delta T_2}\right)}

.. note:: DeltaT2 is a function of hot stream `side 1` and cold stream `side 2`, and DeltaT3 is a function of hot side and cold stream `side 3`.

Class Documentation
-------------------

.. autoclass:: HeatExchangerWith3Streams
   :members:

.. autoclass:: HeatExchangerWith3StreamsData
   :members:
