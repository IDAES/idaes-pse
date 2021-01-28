International Association of the Properties of Water and Steam IAPWS-95
=======================================================================

.. index::
  pair: idaes.generic_models.properties.iapws95; Iapws95StateBlock

.. module:: idaes.generic_models.properties.iapws95

Accurate and thermodynamically consistent steam properties are provided for the
IDAES framework by implementing the International Association for the Properties
of Water and Steam's :ref:`"Revised Release on the IAPWS Formulation 1995 for
the Thermodynamic Properties of Ordinary Water Substance for General and
Scientific Use." <iapws-2016>` Non-analytic terms designed to improve accuracy
very near the critical point were omitted, because they cause a singularity at
the critical point, a feature which is undesirable in optimization problems. The
IDAES implementation provides features which make the water and steam property
calculations amenable to rigorous mathematical optimization.

Please see the :ref:`general Helmholtz documentation <technical_specs/model_libraries/generic/property_models/helmholtz:Pure Component Helmholtz EoS>`
for more information.

Example
-------

The Heater unit model
:ref:`example <technical_specs/model_libraries/generic/unit_models/heater:Example>`,
provides a simple example for using water properties.

.. testcode::

  from idaes.generic_models.properties import iapws95
  import pyomo.environ as pe # Pyomo environment
  from idaes.core import FlowsheetBlock, MaterialBalanceType
  from idaes.generic_models.unit_models import Heater

  # Create an empty flowsheet and steam property parameter block.
  model = pe.ConcreteModel()
  model.fs = FlowsheetBlock(default={"dynamic": False})
  model.fs.properties = iapws95.Iapws95ParameterBlock(default={
    "phase_presentation":iapws95.PhaseType.LG,
    "state_vars":iapws95.StateVars.PH})

  # Add a Heater model to the flowsheet.
  model.fs.heater = Heater(default={
    "property_package": model.fs.properties,
    "material_balance_type": MaterialBalanceType.componentTotal})

  # Setup the heater model by fixing the inputs and heat duty
  model.fs.heater.inlet[:].enth_mol.fix(4000)
  model.fs.heater.inlet[:].flow_mol.fix(100)
  model.fs.heater.inlet[:].pressure.fix(101325)
  model.fs.heater.heat_duty[:].fix(100*20000)

  # Initialize the model.
  model.fs.heater.initialize()

Since all properties except the state variables are Pyomo Expressions in the
water properties module, after solving the problem any property can be
calculated in any state block. Continuing from the heater example, to get the
viscosity of both phases, the lines below could be added.

.. testcode::

  mu_l = pe.value(model.fs.heater.control_volume.properties_out[0].visc_d_phase["Liq"])
  mu_v = pe.value(model.fs.heater.control_volume.properties_out[0].visc_d_phase["Vap"])

For more information about how StateBlocks and PropertyParameterBlocks work see
the :ref:`StateBlock documentation <technical_specs/core/physical_property_class:Physical Property
Package Classes>`.

Expressions
-----------

The IAPWS-95 property package contains the standard expressions described in
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

Iapws95StateBlock Class
------------------------

.. autoclass:: Iapws95StateBlock
  :members:

Iapws95StateBlockData Class
---------------------------

.. autoclass:: Iapws95StateBlockData
  :members:

Iapws95ParameterBlock Class
---------------------------

.. autoclass:: Iapws95ParameterBlock
  :members:

Iapws95ParameterBlockData Class
-------------------------------

.. autoclass:: Iapws95ParameterBlockData
  :members:

References
----------

.. _iapws-2016:

International Association for the Properties of Water and Steam (2016).
IAPWS R6-95 (2016), "Revised Release on the IAPWS Formulation 1995 for
the Properties of Ordinary Water Substance for General Scientific Use,"
URL: http://iapws.org/relguide/IAPWS95-2016.pdf

.. _wagner-2002:

Wagner, W.,  A. Pruss (2002). "The IAPWS Formulation 1995 for the
Thermodynamic Properties of Ordinary Water Substance for General and
Scientific Use." J. Phys. Chem. Ref. Data, 31, 387-535.

.. _wagner-2000:

Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997 for the
Thermodynamic Properties of Water and Steam," ASME J. Eng. Gas Turbines
and Power, 122, 150-182.

.. _akasaka-2008:

Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
State from Helmholtz Energy Equations of State." Journal of Thermal
Science and Technology, 3(3), 442-451.

.. _iapws-2011:

International Association for the Properties of Water and Steam (2011).
IAPWS R15-11, "Release on the IAPWS Formulation 2011 for the
Thermal Conductivity of Ordinary Water Substance,"
URL: http://iapws.org/relguide/ThCond.pdf.

.. _iapws-2008:

International Association for the Properties of Water and Steam (2008).
IAPWS R12-08, "Release on the IAPWS Formulation 2008 for the Viscosity of
Ordinary Water Substance,"
URL: http://iapws.org/relguide/visc.pdf.
