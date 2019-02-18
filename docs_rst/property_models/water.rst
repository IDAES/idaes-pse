Water/Steam
===========

Two property modules are available for pure water and steam properties.  The
property modules use the same calculations and yield consistent results, but one
uses pressure and molar enthalpy as state variables and the other uses
temperature, pressure, and vapor fraction as the state variables.

Theses modules can be imported as:

.. testcode::

  from idaes.property_models import iapws95_ph
  # from idaes.property_models import iapws95_tpx

Example
-------

The Heater unit model :ref:`example <models/heater:Example>`, provides a simple
example for using water properties.

Since all properties except the state variables are Pyomo Expressions in the
water properties module, after solving the problem any property can be
calculated in any state block after the problem is solved. To get the viscosity
of the water at the heater outlet, for example the line below could be added.

.. testsetup::

  import pyomo.environ as pe # Pyomo environment
  from idaes.core import FlowsheetBlock, StateBlock
  from idaes.unit_models import Heater
  from idaes.property_models import iapws95_ph

  # Create an empty flowsheet and steam property parameter block.
  model = pe.ConcreteModel()
  model.fs = FlowsheetBlock(default={"dynamic": False})
  model.fs.properties = iapws95_ph.Iapws95ParameterBlock()

  # Add a Heater model to the flowsheet.
  model.fs.heater = Heater(default={"property_package": model.fs.properties})

  # Setup the heater model by fixing the inputs and heat duty
  model.fs.heater.inlet[:].enth_mol.fix(4000)
  model.fs.heater.inlet[:].flow_mol.fix(100)
  model.fs.heater.inlet[:].pressure.fix(101325)
  model.fs.heater.heat_duty[:].fix(100*20000)

  # Initialize the model.
  model.fs.heater.initialize()


.. testcode::

  import pyomo.environ as pe
  mu_l = pe.value(model.fs.heater.control_volume.properties_out[0].visc_d_phase["Liq"])
  mu_v = pe.value(model.fs.heater.control_volume.properties_out[0].visc_d_phase["Vap"])

For more information about how StateBlocks and PropertyParameterBlocks
work see the :ref:`StateBlock documentation <core/state_block:Physical
Property Package Classes>`.

Units
-----

The water property modules are in SI units (m, kg, s, mol). Temperature is in K.

Methods
-------

These methods use the IAPWS-95 formulation for scientific use for thermodynamic
properties (:ref:`Wagner and Pruss, 2002 <wagner-2002>`; :ref:`IAPWS, 2016
<iapws-2016>`). To solve the phase equilibrium, the method of :ref:`Akasaka
(2008) <akasaka-2008>` was used. For solving these equations, some relations from
the IAPWS-97 formulation for industrial use are used as initial values
(:ref:`Wagner et al., 2002 <wagner-2002>`). The industrial formulation is
slightly discontinuous between different regions, so it may not be suitable for
optimization. In addition to thermodynamic quantities, viscosity and thermal
conductivity are calculated (:ref:`IAPWS, 2008 <iapws-2008>`;
:ref:`IAPWS, 2011 <iapws-2011>`).


External Functions
------------------

The IAPWS-95 formulation uses density and temperature as state variables. For
most applications those sate variables are not the most convenient choices. Using
other state variables requires solving an equation or two to get density and
temperature from the chosen state variables. This can have numerous solutions
only one of which is physically meaningful. Rather than solve these equations as
part of the full process simulation, external functions were developed that can
solve the equations required to change state variables and guarantee the correct
roots.

The external property functions are written in C++ and complied such that they
can be called by AMPL solvers.  See the :ref:`installation instructions
<install:Installation Instructions>` for information about compiling these
functions. The external functions provide both first and second derivatives for
all property function calls, however at phase transitions some of these functions
may be non-smooth.

IDAES Framework Wrappers
------------------------

Wrappers for these function are provided for compatibility with the IDAES
framework. Some methods for dealing with non-smoothness may also be included in
the IDAES wrappers. The wrappers provide most properties in the form of Pyomo
Expressions, with only the chosen set of state variables being Pyomo Vars. The
expressions pass the state variables to the external functions and do unit
conversion to put the results in SI units. This means that only the state
variables can be fixed and other quantities cannot be fixed, but constraints
can be added to set them to a specific value.

Since state variables are calculated when solving a model, and the rest of the
properties are Expressions, any property available can be easily calculated
after the model is solved, whether is was needed in the model or not.

Although not generally used the wrappers provide direct access to the
ExternalFunctions also. For more information see section :ref:`ExternalFunctions
<property_models/water:ExternalFunctions>`

Pressure-Enthalpy
~~~~~~~~~~~~~~~~~

Although Expressions for properties of different phases are available, the
pressure-enthalpy formulation treats the fluid as a single mixed phase with a
vapor fraction.  This bypasses some of the IDAES framework phase equilibrium
mechanisms and phase equilibrium is always calculated.

The advantage of this choice of state variables is that it is very robust when
phase changes occur, and is especially useful when it is not known if a phase
change will occur.  The disadvantage of this choice of state variables is that
for equations like heat transfer equations that are highly dependent on
temperature, a model could be harder to solve near regions with phase change.
Temperature is a non-smooth function with a zero derivative with respect to
enthalpy in the two-phase region.

The variables for this form are ``flow_mol`` (mol/s), ``pressure`` (Pa), and
``enth_mol`` (J/mol).

Since temperature and vapor fraction are not state variables in this formulation,
they are provided by expressions, and cannot be fixed.  For example, to set a
temperature to a specific value, a constraint could be added which says the
temperature expression equals a fixed value.

These expressions are specific to the P-H formulation:

``temperature``
  Expression that calculates temperature by calling an ExternalFunction of
  enthalpy and pressure. This expression is non-smooth in the transition from
  single-phase to two-phase and has a zero derivative with respect to enthalpy
  in the two-phase region.
``vapor_frac``
  Expression that calculates vapor fraction by calling an ExternalFunction of
  enthalpy and pressure. This expression is non-smooth in the transition from
  single-phase to two-phase and has a zero derivative with respect to enthalpy
  in the single-phase region, where the value is 0 (liquid) or 1 (vapor).

Temperature-Pressure-Vapor Fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coming soon.


Expressions
~~~~~~~~~~~

Unless otherwise noted above, the property expressions are common to both the
T-P-x and P-H formulations. For phase specific properties, valid phase indexes
are ``"Liq"`` and ``"Vap"``

================================ ===============================================================================================
Expression                       Description
================================ ===============================================================================================
``mw``                           Molecular weight (kg/mol)
``tau``                          Critical temperature divided by temperature (unitless)
``temperature_red``              Reduced temperature, temperature divided by critical temperature (unitless)
``temperature_sat``              Saturation temperature (K)
``tau_sat``                      Critical temperature divided by saturation temperature (unitless)
``pressure_sat``                 Saturation pressure (Pa)
``dens_mass_phase[phase]``       Density phase (kg/m\ :superscript:`3`)
``dens_phase_red[phase]``        Phase reduced density (:math:`\delta`), mass density divided by critical density (unitless)
``dens_mass``                    Total mixed phase mass density (kg/m\ :superscript:`3`)
``dens_mol``                     Total mixed phase mole density (kg/m\ :superscript:`3`)
``flow_vol``                     Total volumetric flow rate (m\ :superscript:`3`/s)
``enth_mass``                    Mass enthalpy (kJ/kg)
``enth_mol_sat_phase[phase]``    Saturation enthalpy of phase, enthalpy at P and T\ :subscript:`sat` (kJ/mol)
``enth_mol_phase[phase]``        Molar enthalpy of phase (kJ/mol)
``entr_mol_phase``               Molar entropy of phase (kJ/mol/K)
``entr_mol``                     Total mixed phase entropy (kJ/mol/K)
``cp_mol_phase[phase]``          Constant pressure molar heat capacity of phase (kJ/mol/K)
``cv_mol_phase[phase]``          Constant pressure volume heat capacity of phase (kJ/mol/K)
``cp_mol``                       Total mixed phase constant pressure heat capacity (kJ/mol/K)
``cv_mol``                       Total mixed phase constant volume heat capacity (kJ/mol/K)
``heat_capacity_ratio``          :code:`cp_mol/cv_mol`
``speed_sound_phase[phase]``     Speed of sound in phase (m/s)
``dens_mol_phase[phase]``        Mole density of phase (mol/m\ :superscript:`3`)
``therm_cond_phase[phase]``      Thermal conductivity of phase (W/K/m)
``visc_d_phase[phase]``          Viscosity of phase (Pa/s)
``visc_k_phase[phase]``          Kinimatic viscosity of phase (m\ :superscript:`2`/s)
``phase_frac[phase]``            Phase fraction
``flow_mol_comp["H2O"]``         Same as total flow since only water (mol/s)
================================ ===============================================================================================

ExternalFunctions
~~~~~~~~~~~~~~~~~

This provides a list of ExternalFuctions available in the wrappers.  These
functions do not use SI units and are not usually called directly.  If these
functions are needed, they should be used with caution. Some of these are used
in the property expressions, some are just provided to allow easier testing with
a Python framework.

All of these functions provide first and second derivative and are generally
suited to optimization. Some functions may have non-smoothness at phase
transitions.  The delta_vap and delta_liq functions return the same values in
the critical region.  They will also return real values when a phase doesn't
exist, but those values do not necessarily have physical meaning.

There are a few variables that are common to a lot of these functions, so they
are summarized here :math:`tau` is the critical temperature divided by the
temperature :math:`\frac{T_c}{T}`, :math:`\delta` is density divided by the
critical density :math:`\frac{\rho}{\rho_c}`, and :math:`\phi` is Helmholtz free
energy divided by the ideal gas constant and temperature :math:`\frac{f}{RT}`.

==================== ============== ================================================================ ===========================
Pyomo Function       C Function     Returns                                                          Arguments
==================== ============== ================================================================ ===========================
func_p               p              pressure (kPa)                                                   :math:`\delta, \tau`
func_u               u              internal energy (kJ/kg)                                          :math:`\delta, \tau`
func_s               s              entropy (kJ/K/kg)                                                :math:`\delta, \tau`
func_h               h              enthalpy (kJ/kg)                                                 :math:`\delta, \tau`
func_hvpt            hvpt           vapor enthalpy (kJ/kg)                                           P (kPa), :math:`\tau`
func_hlpt            hlpt           liquid enthalpy (kJ/kg)                                          P (kPa), :math:`\tau`
func_tau             tau            :math:`\tau` (unitless)                                          h (kJ/kg), P (kPa)
func_vf              vf             vapor fraction (unitless)                                        h (kJ/kg), P (kPa)
func_g               g              Gibbs free energy (kJ/kg)                                        :math:`\delta, \tau`
func_f               f              Helmholtz free energy (kJ/kg)                                    :math:`\delta, \tau`
func_cv              cv             const. volume heat capacity (kJ/K/kg)                            :math:`\delta, \tau`
func_cp              cp             const. pressure heat capacity (kJ/K/kg)                          :math:`\delta, \tau`
func_w               w              speed of sound (m/s)                                             :math:`\delta, \tau`
func_delta_liq       delta_liq      liquid :math:`\delta` (unitless)                                 P (kPa), :math:`\tau`
func_delta_vap       delta_vap      vapor :math:`\delta` (unitless)                                  P (kPa), :math:`\tau`
func_delta_sat_l     delta_sat_l    sat. liquid :math:`\delta` (unitless)                            :math:`\tau`
func_delta_sat_v     delta_sat_v    sat. vapor :math:`\delta` (unitless)                             :math:`\tau`
func_p_sat           p_sat          sat. pressure (kPa)                                              :math:`\tau`
func_tau_sat         tau_sat        sat. :math:`\tau` (unitless)                                     P (kPa)
func_phi0            phi0           :math:`\phi` idaes gas part (unitless)                           :math:`\delta, \tau`
func_phi0_delta      phi0_delta     :math:`\frac{\partial \phi_0}{\partial \delta}`                  :math:`\delta`
func_phi0_delta2     phi0_delta2    :math:`\frac{\partial^2 \phi_0}{\partial \delta^2}`              :math:`\delta`
func_phi0_tau        phi0_tau       :math:`\frac{\partial \phi_0}{\partial \tau}`                    :math:`\tau`
func_phi0_tau2       phi0_tau2      :math:`\frac{\partial^2 \phi_0}{\partial \tau^2}`                :math:`\tau`
func_phir            phir           :math:`\phi` real gas part (unitless)                            :math:`\delta, \tau`
func_phir_delta      phir_delta     :math:`\frac{\partial \phi_r}{\partial \delta}`                  :math:`\delta, \tau`
func_phir_delta2     phir_delta2    :math:`\frac{\partial^2 \phi_r}{\partial \delta^2}`              :math:`\delta, \tau`
func_phir_tau        phir_tau       :math:`\frac{\partial \phi_r}{\partial \tau}`                    :math:`\delta, \tau`
func_phir_tau2       phir_tau2      :math:`\frac{\partial^2 \phi_r}{\partial \tau^2}`                :math:`\delta, \tau`
func_phir_delta_tau  phir_delta_tau :math:`\frac{\partial^2 \phi_r}{\partial \delta \partial \tau}`  :math:`\delta, \tau`
==================== ============== ================================================================ ===========================


Initialization
--------------

The IAPWS-95 property functions do provide initilaization functions for general
compatibility with the IDAES framework, but as long as the state variables are
specified to some resonalbe value, initialization is not required. All required
solves are handled by external functions.

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
