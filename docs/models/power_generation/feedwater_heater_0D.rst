Feedwater Heater (0D)
=====================

.. index::
  pair: idaes.unit_models.power_generation.feedwater_heater_0D;FWH0D

.. module:: idaes.unit_models.power_generation.feedwater_heater_0D

The FWH0D model is a 0D feedwater heater model suitable for steady state modeling.  It is intended to be used primarily used with the :ref:`IAWPS95 <property_models/water:Water/Steam>` property package. The feedwater heater is split into three sections the condensing section is required while the desuperheating and drain cooling sections are optional. There is also an optional mixer for adding a drain stream from another feedwater heater to the condensing section.  The figure below shows the layout of the feedwater heater.  All but the condensing section are optional.

.. figure:: feedwater_heater_0D.svg
  :width: 800
  :align: center

  Feedwater Heater


Model Structure
---------------

The condensing section uses the :ref:`FWHCondensing0D <models/power_generation/feedwater_heater_condensing_0D:Feedwater Heater (Condensing Section 0D)>` model to calculate a steam flow rate such that all steam is condensed in the condensing section.  This allows turbine steam extraction rates to be calculated. The other sections are regular  :ref:`HeatExchanger <models/heat_exchanger:HeatExchanger (0D)>` models.  The table below shows the unit models which make up the feedwater heater, and the option to include or exclude them.

=========================== ====================== ====================================================================================================================================================================
Unit                        Option                 Doc
=========================== ====================== ====================================================================================================================================================================
``condense``                --                     Condensing section (:ref:`FWHCondensing0D <models/power_generation/feedwater_heater_condensing_0D:Feedwater Heater (Condensing Section 0D)>`)
``desuperheat``             ``has_desuperheat``    Desuperheating section (:ref:`HeatExchanger <models/heat_exchanger:HeatExchanger (0D)>`)
``cooling``                 ``has_drain_cooling``  Drain cooling section (:ref:`HeatExchanger <models/heat_exchanger:HeatExchanger (0D)>`)
``drain_mix``               ``has_drain_mixer``    Mixer for steam and other FWH drain (:ref:`Mixer <models/mixer:Mixer>`)
=========================== ====================== ====================================================================================================================================================================


Degrees of Freedom
------------------

The ``area`` and ``overall_heat_transfer_coefficient`` should be fixed or constraints should be provided to calculate ``overall_heat_transfer_coefficient``.  If the inlets are also fixed except for the inlet steam flow rate (``inlet_1.flow_mol``), the model will have 0 degrees of freedom.

FWH0D Class
-----------

.. autoclass:: FWH0D
  :members:

FWH0DData Class
---------------

.. autoclass:: FWH0DData
  :members:
