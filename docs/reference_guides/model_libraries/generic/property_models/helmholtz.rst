Pure Component Helmholtz EoS
============================

.. index::
  pair: idaes.generic_models.properties.helmholtz.helmholtz; HelmholtzStateBlock

.. module:: idaes.generic_models.properties.helmholtz.helmholtz

The Helmholtz Equation of State (EoS) classes serve as a common core for pure
component property packages where very accurate and thermodynamically consistent
pure component properties are required. This contains general information.
Thermodynamic properties for all Helmholtz EoS packages are calculated by the core
class only the parameters differ between specific component implementation.
Specific implementations may also contain additional properties such as viscosity
and thermal conductivity. For specific property packages details see the pages below.

.. toctree::
    :maxdepth: 1

    iapws95
    swco2

The basic Helmholtz EoS is described :ref:`"Revised Release on the IAPWS Formulation
1995 for the Thermodynamic Properties of Ordinary Water Substance for General and
Scientific Use." <iapws-2016>`.  The Helmholtz EoS as used in the IAPWS-95 contains
non-analytic terms to improve accuracy near the critical point.  These terms, however
cause a singularity at the critical point and can causes computational difficulty,
so the non-analytic where omitted in the IDAES implementation.

The IDAES implementation is of the Helmholtz EoS makes use of external function
for many of the properties.  Solving the VLE and changing state variables require
solution of non-linear equations with multiple solutions, so solving then externally
provides a method of decomposition where it can be guaranteed that the nonlinear
equations associated with the Helmholtz EoS are solved correctly. The external
functions provide first and second derivatives, and are compatible with
advanced optimization solvers.  Phase change does cause problems to be non-smooth,
but, as a practical matter, problem using the IDAES implementation of the Helmholtz
EoS, still seem to solve well even with phase change.

Units
-----

The iapws95 property module uses SI units (m, kg, s, J, mol) for all public
variables and expressions. Temperature is in K. Note that this means molecular
weight is in the unusual unit of kg/mol.

A few expressions intended used internally and all external function calls
use units of kg, kJ, kPa, and K. These generally are not needed by the end user.


Phase Presentation
------------------

The property package wrapper can present fluid phase information to the
IDAES framework in different ways.  For specifics on how to set the options
see a specific implementation page.

The ``PhaseType.MIX`` option causes to modeling framework to view water and
steam as a single mixed liquid and vapor phase. This generally reduces model
complexity. Phase equilibrium is still calculated and ``vapor_frac`` and
individual phase properties are available, just are they would be with the
two-phase presentation.  The mixed-phase presentation can be used with most
standard unit models that do not provide phase separation.  If phase separation
is required, either use the two-phase presentation or create a custom model.

The ``PhaseType.LG`` option appears to the IDAES framework to be two phases "Vap"
and "Liq".  This option requires one of two unit model options to be set.  You
can use the total material balance option for unit models, to specify that only
one material balance equation should be written not one per phase. The other
possible option is to specify ``has_phase_equlibrium=True``. This will
write a material balance per phase, but will add a phase generation term to the
model. For Helmholtz EoS packages, it is generally recommended that specifying
total material balances is best because it results in a problem with fewer
variables, and phase equilibrium is always calculated by the property package.

There are two single phase options ``PhaseType.L`` and ``PhaseType.G``; these
present a single phase "Liq" or "Vap" to the framework. The vapor fraction will
also always return 0 or 1 as appropriate. These options can be used when the phase
of a fluid is know for certain to only be liquid or only be vapor. For the
temperature-pressure-vapor fraction formulation, this eliminates the
complementarity constraint, but for the enthalpy-pressure formulation, where the
vapor fraction is always calculated, the single phase options probably do not
provide any real benefit.

State Variables
---------------

There is a choice of state variables, pressure-enthalpy and
temperature-pressure-vapor fraction.  In general the enthalpy-pressure form is
preferable. Both the pressure and enthalpy variables are smooth and sufficient
to define the fluid state. For systems where two-phases may be present, it is
expected that pressure-enthalpy is the best choice of state variables.

The temperature-pressure-vapor fraction form is more convenient, since temperature
is directly measurable and more familiar than enthalpy. Complementarity
constraints are used to deal with the vapor fraction variable, but the additional
complimentary constraints may make the problem less robust.  Temperature-pressure
is often a good choice of state variables where there is only a single known phase.

Pressure-Enthalpy Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The advantage of this choice of state variables is that it is more robust when
phase changes occur, and is especially useful when it is not known if a phase
change will occur. The disadvantage of this choice of state variables is that
for equations like heat transfer equations that are highly dependent on
temperature, a model could be harder to solve near regions with phase change.
Temperature is a non-smooth function with non-smoothness when transitioning
from the single-phase to the two-phase region. Temperature also has a zero
derivative with respect to enthalpy in the two-phase region, so near the
two-phase region solving a constraint that specifies a specific temperature
may be difficult.

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
fixed to the appropriate value and the complementarity constraint is deactivated.

A complementarity constraint is required for the T-P-x formulation when two-phases
may be present.  First, two expressions are defined below where :math:`P^-` is
pressure under saturation pressure and :math:`P^+` is pressure over saturation
pressure. The :math:`\max()` function is provided as an IDAES utility which
provides a smooth max expression.

.. math::

  P^- = \max(0, P_{\text{sat}} - P)

.. math::

  P^+ = \max(0, P - P_{\text{sat}})

With the "pressure over" and "pressure under" expressions a complementarity
constraint can be written.  If the pressure under saturation is more than zero,
only vapor exists.  If the pressure over saturation is greater than zero only a
liquid exists.  If both are about zero two phases can exist. The saturation
pressure function maxes out at the critical pressure and any temperature above
the critical temperature will yield a saturation pressure that is the critical
pressure, so supercritical fluids will be classified as liquids as is the
convention for this property package.

.. math::

  0 = xP^+  - (1 - x)P^-

Assuming the vapor fraction (:math:`x`) is positive and noting that only one of
:math:`P^+` and :math:`P^-` can be nonzero (approximately), the complementarity
equation above requires :math:`x` to be 0 when :math:`P^+` is not zero (liquid)
or :math:`x` to be 1 when :math:`P^-`` is not zero (vapor).  When both
:math:`P^+` and :math:`P^-`` are about 0, the complementarity constraint says
nothing about x, but it basically reduces another constraint, that
:math:`P=P_{\text{sat}}`. When two phases are present :math:`x` is found
by the unit model energy balance, where the temperature will be
:math:`T_{\text{sat}}` (because :math:`P=P_{\text{sat}}`).

An alternative approach is sometimes useful to simplify the problem when it is
certain that there are two phases. The complementarity constraint can be
deactivated and a :math:`P=P_{\text{sat}}` or :math:`T=T_{\text{sat}}` constraint
can be added.

Using the T-P-x formulation requires better initial guesses than the P-H form.
It is not strictly necessary but it is best to try to get an initial guess that
is in the correct phase region for the expected result model.

Expressions
~~~~~~~~~~~

Unless otherwise noted, the property expressions are common to both the
T-P-x and P-H formulations. For phase specific properties, valid phase indexes
are ``"Liq"`` and ``"Vap"``.  Even when using the mixed phase version of the
property package, both liquid and vapor properties are available.

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
``vapor_frac``                       Vapor fraction, if PH form
``phase_frac[phase]``                Phase fraction
``flow_mol_comp["H2O"]``             Same as total flow since only water (mol/s)
``P_under_sat``                      Pressure under saturation pressure (kPa)
``P_over_sat``                       Pressure over saturation pressure (kPa)
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
Object               C Function     Returns                                                          Arguments
==================== ============== ================================================================ ===========================
func_p               p              pressure (kPa)                                                   :math:`\delta, \tau`
func_p_stau          p_stau         pressure (kPa)                                                   s (kJ/kg/K), :math:`\tau`
func_u               u              internal energy (kJ/kg)                                          :math:`\delta, \tau`
func_s               s              entropy (kJ/K/kg)                                                :math:`\delta, \tau`
func_h               h              enthalpy (kJ/kg)                                                 :math:`\delta, \tau`
func_hvpt            hvpt           vapor enthalpy (kJ/kg)                                           P (kPa), :math:`\tau`
func_hlpt            hlpt           liquid enthalpy (kJ/kg)                                          P (kPa), :math:`\tau`
func_svpt            svpt           vapor entropy (kJ/kg/K)                                          P (kPa), :math:`\tau`
func_slpt            slpt           liquid entropy (kJ/kg/K)                                         P (kPa), :math:`\tau`
func_uvpt            uvpt           vapor internal energy (kJ/kg)                                    P (kPa), :math:`\tau`
func_ulpt            ulpt           liquid internal energy (kJ/kg)                                   P (kPa), :math:`\tau`
func_tau             tau            :math:`\tau` (unitless)                                          h (kJ/kg), P (kPa)
func_tau_sp          tau_sp         :math:`\tau` (unitless)                                          s (kJ/kg/K), P (kPa)
func_tau_up          tau_up         :math:`\tau` (unitless)                                          u (kJ/kg), P (kPa)
func_vf              vf             vapor fraction (unitless)                                        h (kJ/kg), P (kPa)
func_vfs             vfs            vapor fraction (unitless)                                        s (kJ/kg/K), P (kPa)
func_vfu             vfu            vapor fraction (unitless)                                        u (kJ/kg), P (kPa)
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

Although a general Helmholtz EoS was developed, the equations where taken from
the IAPWS-95 standard. For specific parameter sources see specific implementation
documentation.

.. _iapws-2016:

International Association for the Properties of Water and Steam (2016).
IAPWS R6-95 (2016), "Revised Release on the IAPWS Formulation 1995 for
the Properties of Ordinary Water Substance for General Scientific Use,"
URL: http://iapws.org/relguide/IAPWS95-2016.pdf

.. _wagner-2002:

Wagner, W.,  A. Pruss (2002). "The IAPWS Formulation 1995 for the
Thermodynamic Properties of Ordinary Water Substance for General and
Scientific Use." J. Phys. Chem. Ref. Data, 31, 387-535.

.. _akasaka-2008:

Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
State from Helmholtz Energy Equations of State." Journal of Thermal
Science and Technology, 3(3), 442-451.
