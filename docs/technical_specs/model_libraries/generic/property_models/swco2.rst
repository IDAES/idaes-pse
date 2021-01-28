Span-Wager CO2
==============

.. index::
  pair: idaes.generic_models.properties.swco2; SWCO2StateBlock

.. module:: idaes.generic_models.properties.swco2

This implements the Span-Wagner equation of state for CO2   :ref:`"Span-Wagner equation of state for CO2" <span-1996>`
Please see the :ref:`general Helmholtz documentation <technical_specs/model_libraries/generic/property_models/helmholtz:Pure Component Helmholtz EoS>`
for more information.

Example
-------

The Heater unit model
:ref:`example <technical_specs/model_libraries/generic/unit_models/heater:Example>`,
provides a simple example for using water properties.

.. testcode::

  from idaes.generic_models.properties import swco2
  import pyomo.environ as pe # Pyomo environment
  from idaes.generic_models.unit_models import Compressor
  from idaes.core import FlowsheetBlock

  m = ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = swco2.SWCO2ParameterBlock()
  m.fs.unit = Compressor(default={"property_package": m.fs.properties})
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
  solver.solve(model)


For more information about how StateBlocks and PropertyParameterBlocks work see
the :ref:`StateBlock documentation <technical_specs/core/physical_property_class:Physical Property
Package Classes>`.

Expressions
-----------

The Span-Wager property package contains the standard expressions described in
the :ref:`general Helmholtz documentation <technical_specs/model_libraries/generic/property_models/helmholtz:Pure Component Helmholtz EoS>`,
but it also defines expressions for transport properties.

==================================== =====================================================
Expression                           Description
==================================== =====================================================
``therm_cond_phase[phase]``          Thermal conductivity of phase (W/K/m)
``visc_d_phase[phase]``              Viscosity of phase (Pa/s)
``visc_k_phase[phase]``              Kinimatic viscosity of phase (m\ :superscript:`2`/s)
==================================== =====================================================

Convenience Functions
---------------------

.. autofunction:: htpx

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
