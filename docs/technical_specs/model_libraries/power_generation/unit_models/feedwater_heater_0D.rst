Feedwater Heater (0D)
=====================

.. index::
    pair: idaes.power_generation.unit_models.feedwater_heater_0D;FWH0D

.. module:: idaes.power_generation.unit_models.feedwater_heater_0D
  :noindex:

The FWH0D model is a 0D feedwater heater model suitable for steady state modeling.
It is intended to be used primarily with the
:ref:`IAWPS95 <technical_specs/model_libraries/generic/property_models/iapws95:International Association of the Properties of Water and Steam IAPWS-95>` property package.
The feedwater heater is split into three sections the condensing section is required while
the desuperheating and drain cooling sections are optional. There is also an optional mixer
for adding a drain stream from another feedwater heater to the condensing section.  The figure
below shows the layout of the feedwater heater.  All but the condensing section are optional.

.. figure:: feedwater_heater_0D.svg
  :width: 800
  :align: center

  Feedwater Heater


Example
-------

The example below shows how to setup a feedwater heater with all tree sections.  The feedwater flow rate, steam conditions, heat transfer coefficients and areas are not necessarily realistic.

.. testcode::

  import pyomo.environ as pyo
  from idaes.core import FlowsheetBlock
  from idaes.generic_models.unit_models.heat_exchanger import (delta_temperature_underwood_callback,
      delta_temperature_lmtd_callback)
  from idaes.generic_models.properties import iapws95
  from idaes.power_generation.unit_models import FWH0D

  def make_fwh_model():
      model = pyo.ConcreteModel()
      model.fs = FlowsheetBlock(default={
          "dynamic": False,
          "default_property_package": iapws95.Iapws95ParameterBlock()})
      model.fs.properties = model.fs.config.default_property_package
      model.fs.fwh = FWH0D(default={
          "has_desuperheat":True,
          "has_drain_cooling":True,
          "has_drain_mixer":True,
          "property_package":model.fs.properties})

      model.fs.fwh.desuperheat.inlet_1.flow_mol.fix(100)
      model.fs.fwh.desuperheat.inlet_1.flow_mol.unfix()
      model.fs.fwh.desuperheat.inlet_1.pressure.fix(201325)
      model.fs.fwh.desuperheat.inlet_1.enth_mol.fix(60000)
      model.fs.fwh.drain_mix.drain.flow_mol.fix(1)
      model.fs.fwh.drain_mix.drain.pressure.fix(201325)
      model.fs.fwh.drain_mix.drain.enth_mol.fix(20000)
      model.fs.fwh.cooling.inlet_2.flow_mol.fix(400)
      model.fs.fwh.cooling.inlet_2.pressure.fix(101325)
      model.fs.fwh.cooling.inlet_2.enth_mol.fix(3000)
      model.fs.fwh.condense.area.fix(1000)
      model.fs.fwh.condense.overall_heat_transfer_coefficient.fix(100)
      model.fs.fwh.desuperheat.area.fix(1000)
      model.fs.fwh.desuperheat.overall_heat_transfer_coefficient.fix(10)
      model.fs.fwh.cooling.area.fix(1000)
      model.fs.fwh.cooling.overall_heat_transfer_coefficient.fix(10)

      model.fs.fwh.initialize()
      return(model)

  # create a feedwater heater model with all optional units and initialize
  model = make_fwh_model()

Model Structure
---------------

The condensing section uses the
:ref:`FWHCondensing0D <technical_specs/model_libraries/power_generation/unit_models/feedwater_heater_condensing_0D:Feedwater Heater (Condensing Section 0D)>`
model to calculate a steam flow rate such that all steam is condensed in the condensing
section.  This allows turbine steam extraction rates to be calculated. The other sections
are regular
:ref:`HeatExchanger <technical_specs/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>` models.
The table below shows the unit models which make up the feedwater heater, and the option to
include or exclude them.

=========================== ====================== ====================================================================================================================================================================
Unit                        Option                 Doc
=========================== ====================== ====================================================================================================================================================================
``condense``                --                     Condensing section (:ref:`FWHCondensing0D <technical_specs/model_libraries/power_generation/unit_models/feedwater_heater_condensing_0D:Feedwater Heater (Condensing Section 0D)>`)
``desuperheat``             ``has_desuperheat``    Desuperheating section (:ref:`HeatExchanger <technical_specs/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>`)
``cooling``                 ``has_drain_cooling``  Drain cooling section (:ref:`HeatExchanger <technical_specs/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>`)
``drain_mix``               ``has_drain_mixer``    Mixer for steam and other FWH drain (:ref:`Mixer <technical_specs/model_libraries/generic/unit_models/mixer:Mixer>`)
=========================== ====================== ====================================================================================================================================================================


Degrees of Freedom
------------------

The ``area`` and ``overall_heat_transfer_coefficient`` should be fixed or constraints should be provided to calculate ``overall_heat_transfer_coefficient``.  If the inlets are also fixed except for the inlet steam flow rate (``inlet_1.flow_mol``), the model will have 0 degrees of freedom.

See :class:`FWH0D` and :class:`FWH0DData` for full Python class details.
