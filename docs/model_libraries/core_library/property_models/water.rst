Water/Steam - IAPWS95
======================

.. index::
  pair: idaes.core_lib.properties.iapws95; Iapws95StateBlock

.. module:: idaes.core_lib.properties.iapws95

Accurate and thermodynamically consistent steam properties are provided for the
IDAES framework by implementing the International Association for the Properties
of Water and Steam's :ref:`"Revised Release on the IAPWS Formulation 1995 for
the Thermodynamic Properties of Ordinary Water Substance for General and
Scientific Use." <iapws-2016>` Non-analytic terms designed to improve accuracy
very near the critical point were omitted, because they cause a singularity at
the critical point, a feature which is undesirable in optimization problems. The
IDAES implementation provides features which make the water and steam property
calculations amenable to rigorous mathematical optimization.

Example
-------

Theses modules can be imported as:

.. testcode::

  from idaes.core_lib.properties import iapws95

The Heater unit model :ref:`example <models/heater:Example>`, provides a simple
example for using water properties.

.. testcode::

  import pyomo.environ as pe # Pyomo environment
  from idaes.core import FlowsheetBlock, MaterialBalanceType
  from idaes.unit_models import Heater
  from idaes.core_lib.properties import iapws95

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
the :ref:`StateBlock documentation <core/state_block:Physical Property Package
Classes>`.

Units
-----

The iapws95 property module uses SI units (m, kg, s, J, mol) for all public
variables and expressions. Temperature is in K. Note that this means molecular
weight is in the unusual unit of kg/mol.

A few expressions intended to be used internally and all external function calls
use units of kg, kJ, kPa, and K.  These generally are not needed by the end user.

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
most applications those state variables are not the most convenient choices.
Using other state variables requires solving equations to get density and
temperature from the chosen state variables. These equations can have numerous
solutions only one of which is physically meaningful. Rather than solve these
equations as part of the full process simulation, external functions were
developed that can solve the equations required to change state variables and
guarantee the correct roots.

The external property functions are written in C++ and complied such that they
can be called by AMPL solvers.  See the :ref:`idaes_installation` page for
information about compiling these functions. The external functions provide
both first and second derivatives for all property function calls, however at
phase transitions some of these functions may be non-smooth.

IDAES Framework Wrapper
-----------------------

A wrapper for the external functions is provided for compatibility with the IDAES
framework. Most properties are available as Pyomo Expressions from the wrapper.
Only the state variables are model variables. Benefits of using mostly
expressions in the property package are: no initialization is required
specifically for the property package, the model has fewer equations, and
all properties can be easily calculated after the model is solved from the
state variable values even if they were not used in the model. Calls to the
external functions are used within expressions so users do not need to directly
call any functions. The potential downside of the extensive use of expressions
here is that combining the expressions to form constraints could yield equations
that are more difficult to solve than, they would have been if an equivalent
system of equations was written with more variables and simpler equations.
Quantifying the effect of writing larger equations with fewer variables is
difficult. Experience suggests in this particular case more expressions and fewer
variables is better.

Although not generally used, the wrapper provides direct access to the
ExternalFunctions, including intermediate functions. For more information see
section :ref:`ExternalFunctions <property_models/water:ExternalFunctions>`.
These are mostly available for testing purposes.

Phase Presentation
~~~~~~~~~~~~~~~~~~

The property package wrapper can present fluid phase information to the
IDAES framework in different ways.  See the
:ref:`class reference <property_models/water:Iapws95ParameterBlock Class>` for details
on how to set these options.  The ``phase_presentation=PhaseType.MIX`` option
looks like one phase called "Mix" to the IDAES framework. The property
package will calculate a phase fraction. This will bypass any two phase
handling equations written for unit models, and should work with any unit model
options as long as you do not want to separate the phases. The benefit of this
option is that it can potentially lead to a simpler set of equations.

The ``phase_presentation=PhaseType.LG`` option appears to the IDAES framework to
be two phases "Vap" and "Liq".  This option requires one of two unit model
options to be set.  You can use the total material balance option for unit
models, to specify that only one material balance equation should be written
not one per phase. The other possible option is to specify
``has_phase_equlibrium=True``. This will still write a material balance
per phase, but will add a phase generation term to the model. For the IAPWS-95
package, it is generally recommended that specifying total material balances is
best because it results in a problem with fewer variables.

There are also two single phase options ``phase_presentation=PhaseType.L`` and
``phase_presentation=PhaseType.G``, these present a single phase "Liq" or "Vap"
to the framework.  The vapor fraction will also always return 0 or 1 as
appropriate. These options can be used when the phase of a fluid is know for
certain to only be liquid or only be vapor. For the temperature-pressure-vapor
fraction formulation, this eliminates the complementarity constraint, but for the
enthalpy-pressure formulation, where the vapor fraction is always calculated,
the single phase options probably do not provide any real benefit.

Pressure-Enthalpy Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The advantage of this choice of state variables is that it is very robust when
phase changes occur, and is especially useful when it is not known if a phase
change will occur.  The disadvantage of this choice of state variables is that
for equations like heat transfer equations that are highly dependent on
temperature, a model could be harder to solve near regions with phase change.
Temperature is a non-smooth function with non-smoothness when transitioning
from the single-phase to the two-phase region. Temperature also has a zero
derivative with respect to enthalpy in the two-phase region, so near the
two-phase region solving a constraint that specifies a specific temperature
may not be possible.

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

This formulation uses temperature (K), pressure (Pa), and vapor fraction as
state variables.  When a single phase option is given, the vapor fraction is
fixed to the appropriate value and not included in the state variable set. For
single phase, the complementarity constraint is also deactivated.

A complementarity constraint is required for the T-P-x formulation.  First, two
expressions are defined below where :math:`P^-` is pressure under saturation
pressure and :math:`P^+` is pressure over saturation pressure. The max function
is provided by an IDAES utility function which provides a smooth max expression.

.. math::

  P^- = \max(0, P_{\text{sat}} - P)

.. math::

  P^+ = \max(0, P - P_{\text{sat}})

With the pressure over and pressure under saturated pressure expressions a
complementarity constraint can be written.  If the pressure under saturation is
more than zero, only vapor exists.  If the pressure over saturation is greater
than zero only a liquid exists.  If both are about zero two phases can exist.
The saturation pressure function maxes out at the critical pressure and any
temperature above the critical temperature will yield a saturation pressure that
is the critical pressure, so supercritical fluids will be classified as liquids
as the convention for this property package.

.. math::

  0 = xP^+  - (1 - x)P^-

Assuming the vapor fraction (:math:`x`) is positive and noting that only one of
:math:`P^+` and :math:`P^-` can be nonzero (approximately), the complementarity
equation above requires :math:`x` to be 0 when :math:`P^+` is not zero (liquid)
or :math:`x` to be 1 when :math:`P^-`` is not zero (vapor).  When both
:math:`P^+` and :math:`P^-`` are about 0, the complementarity constraint says
nothing about x, but it does provide another constraint, that
:math:`P=P_{\text{sat}}`. When two phases are present :math:`x` can be found
by the unit model energy balance and the temperature will be
:math:`T_{\text{sat}}`.

An alternative approach is sometimes useful. If you know for certain that you
have two phases, the complementarity constraint can be deactivated and a
:math:`P=P_{\text{sat}}` or :math:`T=T_{\text{sat}}` constraint can be added.

Using the T-P-x formulation requires better initial guesses than the P-H form.
It is not strictly necessary but it is best to try to get an initial guess that
is in the correct phase region for the expected result model.


Expressions
~~~~~~~~~~~

Unless otherwise noted, the property expressions are common to both the
T-P-x and P-H formulations. For phase specific properties, valid phase indexes
are ``"Liq"`` and ``"Vap"``

==================================== ===============================================================================================
Expression                           Description
==================================== ===============================================================================================
``mw``                               Molecular weight (kg/mol)
``tau``                              Critical temperature divided by temperature (unitless)
``temperature``                      Temperature (K) if PH form
``temperature_red``                  Reduced temperature, temperature divided by critical temperature (unitless)
``temperature_sat``                  Saturation temperature (K)
``tau_sat``                          Critical temperature divided by saturation temperature (unitless)
``pressure_sat``                     Saturation pressure (Pa)
``dens_mass_phase[phase]``           Density phase (kg/m\ :superscript:`3`)
``dens_phase_red[phase]``            Phase reduced density (:math:`\delta`), mass density divided by critical density (unitless)
``dens_mass``                        Total mixed phase mass density (kg/m\ :superscript:`3`)
``dens_mol``                         Total mixed phase mole density (kg/m\ :superscript:`3`)
``flow_vol``                         Total volumetric flow rate (m\ :superscript:`3`/s)
``enth_mass``                        Mass enthalpy (J/kg)
``enth_mol_sat_phase[phase]``        Saturation enthalpy of phase, enthalpy at P and T\ :subscript:`sat` (J/mol)
``enth_mol``                         Molar enthalpy (J/mol) if TPx form
``enth_mol_phase[phase]``            Molar enthalpy of phase (J/mol)
``energy_internal_mol``              molar internal energy (J/mol)
``energy_internal_mol_phase[phase]`` Molar internal energy of phase (J/mol)
``entr_mol_phase``                   Molar entropy of phase (J/mol/K)
``entr_mol``                         Total mixed phase entropy (J/mol/K)
``cp_mol_phase[phase]``              Constant pressure molar heat capacity of phase (J/mol/K)
``cv_mol_phase[phase]``              Constant pressure volume heat capacity of phase (J/mol/K)
``cp_mol``                           Total mixed phase constant pressure heat capacity (J/mol/K)
``cv_mol``                           Total mixed phase constant volume heat capacity (J/mol/K)
``heat_capacity_ratio``              :code:`cp_mol/cv_mol`
``speed_sound_phase[phase]``         Speed of sound in phase (m/s)
``dens_mol_phase[phase]``            Mole density of phase (mol/m\ :superscript:`3`)
``therm_cond_phase[phase]``          Thermal conductivity of phase (W/K/m)
``vapor_frac``                       Vapor fraction, if PH form
``visc_d_phase[phase]``              Viscosity of phase (Pa/s)
``visc_k_phase[phase]``              Kinimatic viscosity of phase (m\ :superscript:`2`/s)
``phase_frac[phase]``                Phase fraction
``flow_mol_comp["H2O"]``             Same as total flow since only water (mol/s)
``P_under_sat``                      Pressure under saturation pressure (kPA)
``P_over_sat``                       Pressure over saturation pressure (kPA)
==================================== ===============================================================================================

ExternalFunctions
~~~~~~~~~~~~~~~~~

This provides a list of ExternalFuctions available in the wrapper.  These
functions do not use SI units and are not usually called directly.  If these
functions are needed, they should be used with caution. Some of these are used
in the property expressions, some are just provided to allow easier testing with
a Python framework.

All of these functions provide first and second derivative and are generally
suited to optimization (including the ones that return derivatives of Helmholtz
free energy). Some functions may have non-smoothness at phase transitions.  The
``delta_vap`` and ``delta_liq`` functions return the same values in the critical
region.  They will also return real values when a phase doesn't exist, but those
values do not necessarily have physical meaning.

There are a few variables that are common to a lot of these functions, so they
are summarized here :math:`\tau` is the critical temperature divided by the
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

The IAPWS-95 property functions do provide initialization functions for general
compatibility with the IDAES framework, but as long as the state variables are
specified to some reasonable value, initialization is not required. All required
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
