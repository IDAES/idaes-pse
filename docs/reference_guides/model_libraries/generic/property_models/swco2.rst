Span-Wager CO2
==============

.. index::
  pair: idaes.models.properties.swco2; SWCO2StateBlock

.. module:: idaes.models.properties.swco2

This implements the Span-Wagner equation of state for CO2 
(:ref:`"Span-Wagner equation of state for CO2" <span-1996>`), and is the same as the generic Helmholtz
equation of state with the component set to ``"co2"``. Please see the 
:ref:`general Helmholtz documentation <reference_guides/model_libraries/generic/property_models/helmholtz:Pure Component Helmholtz EoS>`
for more information.

Example
-------

The Heater unit model
:ref:`example <reference_guides/model_libraries/generic/unit_models/heater:Example>`,
provides a simple example for using water properties.

.. testcode::

  from idaes.models.properties import swco2
  from pyomo.environ import ConcreteModel, units as pyunits, SolverFactory # Pyomo environment
  from idaes.models.unit_models import Compressor
  from idaes.core import FlowsheetBlock

  model = ConcreteModel()
  model.fs = FlowsheetBlock(dynamic=False)
  model.fs.properties = swco2.SWCO2ParameterBlock()
  model.fs.unit = Compressor(property_package=model.fs.properties)
  F = 1000
  Tin = 500
  Pin = 10000
  Pout = 20000
  hin = swco2.htpx(T=Tin*pyunits.K, P=Pin*pyunits.Pa)

  model.fs.unit.inlet.flow_mol[0].fix(F)
  model.fs.unit.inlet.enth_mol[0].fix(hin)
  model.fs.unit.inlet.pressure[0].fix(Pin)
  model.fs.unit.deltaP.fix(Pout - Pin)
  model.fs.unit.efficiency_isentropic.fix(0.9)
  model.fs.unit.initialize(optarg={'tol': 1e-6})

  solver = SolverFactory("ipopt")
  solver.solve(model)


For more information about how StateBlocks and PropertyParameterBlocks work see
the :ref:`StateBlock documentation <reference_guides/core/physical_property_class:Physical Property
Package Classes>`.

SWCO2StateBlock Class
------------------------

.. autoclass:: SWCO2StateBlock
  :members:

SWCO2StateBlockData Class
---------------------------

.. autoclass:: SWCO2StateBlockData
  :members:

SWCO2ParameterBlock Class
---------------------------

.. autoclass:: SWCO2ParameterBlock
  :members:

SWCO2ParameterBlockData Class
-------------------------------

.. autoclass:: SWCO2ParameterBlockData
  :members:

References
----------

.. _span-1996:

Span, R., W. Wagner, W. (1996). "A New Equation of State for Carbon Dioxide
Covering the Fluid Region from the Triple-Point Temperature to 1100 K at
Pressures up to 800 MPa." J. Phys. Chem. Ref. Data, 25(6), 1509-1596.

.. _akasaka-2008:

Akasaka, R. (2008). "A Reliable and Useful Method to Determine the
Saturation State from Helmholtz Energy Equations of State." Journal of
Thermal Science and Technology, 3(3), 442-451.

.. _transport_ref:

Vesovic, V., W.A. Wakeham, G.A. Olchowy, J.V. Sengers, J.T.R. Watson, J.
Millat, (1990). "The transport properties of carbon dioxide." J. Phys.
Chem. Ref. Data, 19, 763-808.

.. _viscosity_ref:

Fenghour, A., W.A. Wakeham, V. Vesovic, (1998). "The Viscosity of Carbon
Dioxide." J. Phys. Chem. Ref. Data, 27, 31-44.