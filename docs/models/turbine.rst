Turbine
=======

The Turbine model is a
:ref:`PressureChanger <models/pressure_changer:Pressure Changer>`,
where the configuration is set so that the "compressor" option can only be False,
and the default "thermodynamic_assumption" is "isentropic."  See the
:ref:`PressureChanger documentation <models/pressure_changer:Pressure Changer>`
for details.

Example
-------

The example below demonstrates the basic Turbine model usage:

.. testcode::

  from idaes.unit_models import Compressor
  from idaes.property_models import iapws95

  m = ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = iapws95.Iapws95ParameterBlock()
  m.fs.unit = Turbine(default={"property_package": m.fs.properties})

  m.fs.unit.inlet.flow_mol[0].fix(1000)
  m.fs.unit.inlet.enth_mol[0].fix(hin = iapws95.htpx(T=800, P=1e7))
  m.fs.unit.inlet.pressure[0].fix(1e7)
  m.fs.unit.deltaP.fix(-2e6)
  m.fs.unit.efficiency_isentropic.fix(0.9)
