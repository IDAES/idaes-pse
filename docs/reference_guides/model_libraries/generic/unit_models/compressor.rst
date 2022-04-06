Compressor
==========

The Compressor model is a
:ref:`PressureChanger <reference_guides/model_libraries/generic/unit_models/pressure_changer:Pressure Changer>`,
where the configuration is set so that the "compressor" option can only be True,
and the default "thermodynamic_assumption" is "isentropic."  See the
:ref:`PressureChanger documentation <reference_guides/model_libraries/generic/unit_models/pressure_changer:Pressure Changer>`
for details.

Example
-------

The example below demonstrates the basic Compressor model usage:

.. testcode::

  import pyomo.environ as pyo
  from idaes.core import FlowsheetBlock
  from idaes.models.unit_models import Compressor
  from idaes.models.properties import iapws95

  m = pyo.ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = iapws95.Iapws95ParameterBlock()
  m.fs.unit = Compressor(default={"property_package": m.fs.properties})

  m.fs.unit.inlet.flow_mol[0].fix(100)
  m.fs.unit.inlet.enth_mol[0].fix(4000)
  m.fs.unit.inlet.pressure[0].fix(101325)

  m.fs.unit.deltaP.fix(50000)
  m.fs.unit.efficiency_isentropic.fix(0.9)
