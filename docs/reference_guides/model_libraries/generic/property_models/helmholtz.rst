Pure Component Helmholtz EoS
============================

.. index::
  pair: idaes.models.properties.general_helmholtz.helmholtz_state; HelmholtzStateBlock

.. module:: idaes.models.properties.general_helmholtz

The Helmholtz Equation of State (EoS) classes serve as a common core for pure
component property packages where accurate and thermodynamically consistent
pure component properties are required. New substances can be added by providing 
parameter and expression files. The Helmholtz EoS functions use ExternalFunction, 
so the IDAES binary extensions are required.

This page describes the standard HelmholtzParameterBlock and HelmholtzStateBlock. For 
more information on accessing property expressions using a function-like interface,
or adding and modifying substance parameters see the following pages. There are also
two modules available for backward compatibility ``iapws95`` and ``swco2``, which are
the same as the general module, just with the pure component automatically set 
to ``"h2o"`` or ``"co2"`` respectively.

.. toctree::
  :maxdepth: 1

  helmholtz_expressions
  helmholtz_parameters
  iapws95
  swco2

Defined Components
------------------

If parameter and expression files are available for a component in the parameter file
directory, they will be registered automatically. IDAES comes with some defined components
and more can be added by users. Each component must have an equation of state model, while
transport models for viscosity, thermal conductivity and surface tension are optional.

Functions listed in this section allow you to discover what components and models are 
available and get references for the models.

.. autofunction:: registered_components

.. autofunction:: viscosity_available
    
.. autofunction:: thermal_conductivity_available
    
.. autofunction:: surface_tension_available

.. autofunction:: component_registered

.. autofunction:: clear_component_registry

.. autofunction:: eos_reference

.. autofunction:: viscosity_reference

.. autofunction:: thermal_conductivity_reference

.. autofunction:: surface_tension_reference

Parameter File Location
-----------------------

You can get or set the parameter path with the functions described below.  When
the parameter path is set components will automatically be registered based on the
available parameter files.

.. autofunction:: get_parameter_path

.. autofunction:: set_parameter_path

Classes
-------

The main class used to define a property method is the HelmholtzParameterBlock. Unit models
usually create their own state blocks, but it may be useful to create state block objects
apart from unit models.

.. autoclass:: HelmholtzParameterBlock
  :members:

.. autoclass:: HelmholtzParameterBlockData
  :members:

.. autoclass:: HelmholtzStateBlock
  :members:

.. autoclass:: HelmholtzStateBlockData
  :members:

.. autoclass:: StateVars

.. autoclass:: PhaseType

.. autoclass:: AmountBasis

Helper Function
---------------

Since conditions may be difficult to specify for some choices of state variables, for example, 
fixing temperature and pressure of an inlet when the state variables are enthalpy and pressure,
helper methods are provided by the Parameter block class.  See the ``htpx()``, ``stpx()``, or
``uptx()`` documentation in the parameter block class above.  

Example
-------

The Heater unit model
:ref:`example <reference_guides/model_libraries/generic/unit_models/heater:Example>`,
provides a simple example for using water properties.

.. testcode::

  import pyomo.environ as pe # Pyomo environment
  from idaes.core import FlowsheetBlock, MaterialBalanceType
  from idaes.models.unit_models import Heater
  from idaes.models.properties.general_helmholtz import (
      HelmholtzParameterBlock,
      PhaseType,
      StateVars,
  )

  # Create an empty flowsheet and steam property parameter block.
  model = pe.ConcreteModel()
  model.fs = FlowsheetBlock(dynamic=False)
  model.fs.properties = HelmholtzParameterBlock(
    pure_component="h2o",
    phase_presentation=PhaseType.LG,
    state_vars=StateVars.PH
  )

  # Add a Heater model to the flowsheet.
  model.fs.heater = Heater(
    property_package=model.fs.properties,
    material_balance_type=MaterialBalanceType.componentTotal
  )

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
the :ref:`StateBlock documentation <reference_guides/core/physical_property_class:Physical Property
Package Classes>`.


Units
-----

SI units are used for property variables and expressions (J, Pa, kg, mol, m, s, W).


Phase Presentation
------------------

The property package wrapper can present fluid phase information to the
IDAES framework in different ways. The ``PhaseType.MIX`` option causes 
the modeling framework to view liquid and vapor as a single mixed liquid
and vapor phase. This generally reduces model complexity. Phase equilibrium
is still calculated and ``vapor_frac`` and individual phase properties are 
available, just as they would be with the two-phase presentation. The 
mixed-phase presentation can be used with most standard unit models that do
not provide phase separation.  If phase separation is required, either use 
the two-phase presentation or create a custom model.

.. warning::
    The "has_phase_equilibrium" argument is ignored when constructing Helmholtz
    property packages using mixed phase presentation. However, setting this to
    `True` may cause errors in unit models as it is not possible to construct
    phase equilibrium transfer terms with only one phase present.

The ``PhaseType.LG`` option appears to the IDAES framework to be two phases "Vap"
and "Liq".  This option requires one of two unit model options to be set.  You
can use the total material balance option for unit models, to specify that only
one material balance equation should be written not one per phase. The other
possible option is to specify ``has_phase_equilibrium=True``. This will
write a material balance per phase, but will add a phase generation term to the
model. For Helmholtz EoS packages, it is generally recommended that specifying
total material balances is best because it results in a problem with fewer
variables, and phase equilibrium is always calculated by the property package.

There are two single phase options ``PhaseType.L`` and ``PhaseType.G``; these
present a single phase "Liq" or "Vap" to the framework. The vapor fraction will
also always return 0 or 1 as appropriate. These options can be used when the phase
of a fluid is known for certain to only be liquid or only be vapor. For the
temperature-pressure-vapor fraction formulation, this eliminates the
complementarity constraint, but for the enthalpy-pressure formulation, where the
vapor fraction is always calculated, the single phase options probably do not
provide any real benefit over mixed phase.

State Variables
---------------

There is a choice of state variables, pressure-enthalpy, pressure-entropy, 
pressure-internal energy and temperature-pressure-vapor fraction.  In general
the enthalpy-pressure form is preferable. Both the pressure and enthalpy 
variables are smooth and sufficient to define the fluid state. For systems 
where two-phases may be present, it is expected that pressure-enthalpy is 
the best choice of state variables.

The temperature-pressure-vapor fraction form is more convenient, since temperature
is directly measurable and more familiar than enthalpy. Complementarity
constraints are used to deal with the vapor fraction variable, but the additional
complimentary constraints may make the problem less robust.  Temperature-pressure
is often a good choice of state variables where there is only a single known phase.

Pressure-Enthalpy, Entropy, or Internal Energy Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The advantage of this choice of state variables is that it is more robust when
phase changes occur, and is especially useful when it is not known if a phase
change will occur. The disadvantage of this choice of state variables is that
for equations like heat transfer that are highly dependent on temperature, a 
model could be harder to solve near regions with phase change. Temperature 
is a non-smooth function with non-smoothness when transitioning from the 
single-phase to the two-phase region. Temperature also has a zero derivative
with respect to enthalpy in the two-phase region, so near the two-phase region
solving a constraint that specifies a specific temperature may be difficult.

When a mass basis is used the variables in these forms are flow_mass (kg/s),
pressure (Pa), and one of ``enth_mass`` (J/kg), ``entr_mass`` (J/kg/K), or 
``energy_ineternal_mass`` (J/kg).

When a mass basis is used the variables in these forms are flow_mol (mol/s),
pressure (Pa), and one of ``enth_mol`` (J/mol), ``entr_mol`` (J/mol/K), or 
``energy_ineternal_mole`` (J/mol).

Since temperature and vapor fraction are not state variables in this formulation,
they are provided by expressions, and cannot be fixed.  For example, to set a
temperature to a specific value, a constraint could be added which says the
temperature expression equals a fixed value. 

Temperature-Pressure-Vapor Fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This formulation uses temperature (K), pressure (Pa), and vapor fraction as
state variables.  When a single phase option is given, the vapor fraction does
not need to be specified and is instead an expression with the appropriate value.

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
is in the correct phase region for the expected result.

Non-Existent Phases
~~~~~~~~~~~~~~~~~~~

This section describes the behavior of specific phase property calculations where that phase
does not exist.

Liquid phase properties calculated where a liquid does not exist will solve for the density root
smoothly extending from the saturation curve into the vapor region. Once the temperature 
corresponding to the density root becomes higher than the critical temperature the vapor density
root is used.

Vapor phase properties calculated where a vapor does not exist will solve for the density root
smoothly extending from the saturation curve into the liquid region. Once the pressure
corresponding to the density root becomes higher than the critical pressure the liquid density
root is used.

This method for non-existing phase properties should provide a smooth buffer when solving
equations.  It also returns properties for superheated liquids or sub-cooled vapor. It also
ensures that both the liquid and vapor phase properties are the same for the supercritical
region. 

Variables
~~~~~~~~~

Variables are listed in the table below.  What is a variable and what is an expression 
depends on the selected state variables and amount basis.  Pressure is always a variable.  

========================================= ===============================================================================================
Variable                                  Description
========================================= ===============================================================================================
``pressure``                              Pressure (Pa)
``flow_mass``                             Mass flow (kg/s), if ``amount_basis=AmountBasis.MASS``; otherwise, expression
``flow_mol``                              Mole flow (mol/s) if ``amount_basis=AmountBasis.MOLE``; otherwise, expression
``temperature``                           Temperature (K), is a variable in T-P-x formulation
``vapor_frac``                            Vapor fraction (dimensionless), is a variable in T-P-x two phase formulation
``enth_mass``                             Specific enthalpy (J/kg), is a variable in P-H formulation with mass basis
``enth_mol``                              Molar enthalpy (J/mol), is a variable in P-H formulation with mole basis
``entr_mass``                             Specific enthalpy (J/kg/K), is a variable in P-S formulation with mole basis
``entr_mol``                              Molar enthalpy (J/mol/K), is a variable in P-S formulation with mole basis
``energy_internal_mass``                  Specific internal energy (J/kg), is a variable in P-U formulation with mass basis
``energy_internal_mol``                   Molar internal energy (J/mol), is a variable in P-U formulation with mole basis
========================================= ===============================================================================================


Expressions
~~~~~~~~~~~

Unless otherwise noted, the property expressions are common to both the
T-P-x and P-H formulations. For phase specific properties, valid phase indexes
are ``"Liq"`` and ``"Vap"``.  Even when using the mixed phase version of the
property package, both liquid and vapor properties are available.

========================================= ===============================================================================================
Expression                                Description
========================================= ===============================================================================================
``flow_mass``                             Mass flow (kg/s), variable if ``amount_basis=AmountBasis.MASS``
``flow_mol``                              Mole flow (mol/s) variable if ``amount_basis=AmountBasis.MOLE``
``temperature``                           Temperature (K), is a variable in T-P-x formulation
``vapor_frac``                            Vapor fraction (dimensionless), is a variable in T-P-x two phase formulation
``enth_mass``                             Specific enthalpy (J/kg), is a variable in P-H formulation with mass basis
``enth_mol``                              Molar enthalpy (J/mol), is a variable in P-H formulation with mole basis
``entr_mass``                             Specific enthalpy (J/kg/K), is a variable in P-S formulation with mole basis
``entr_mol``                              Molar enthalpy (J/mol/K), is a variable in P-S formulation with mole basis
``energy_internal_mass``                  Specific internal energy (J/kg), is a variable in P-U formulation with mass basis
``energy_internal_mol``                   Molar internal energy (J/mol), is a variable in P-U formulation with mole basis
``mw``                                    Molecular weight (kg/mol)
``mole_frac_comp[comp]``                  Mole fraction of component (dimensionless), since pure component, returns 1
``mole_frac_phase_comp[phase][comp]``     Mole fraction of component in phase (dimensionless), since pure component, returns 1
``temperature_crit``                      Critical temperature (K)
``temperature_star``                      Reducing temperature :math:`\tau=\frac{T^*}{T}` (K)
``pressure_crit``                         Critical pressure (Pa)
``dens_mass_crit``                        Critical mass density (kg/m\ :superscript:`3`)
``dens_mass_star``                        Reducing mass density :math:`\delta=\frac{\rho}{\rho^*}` (kg/m\ :superscript:`3`)
``dens_mol_crit``                         Critical mole density (mol/m\ :superscript:`3`)
``dens_mol_star``                         Reducing mole density :math:`\delta=\frac{\rho}{\rho^*}` (kg/m\ :superscript:`3`)
``temperature_sat``                       Saturation temperature (K), if supercritical, Tsat=Tcrit
``pressure_sat``                          Saturation pressure (Pa), if supercritical, Psat=Pcrit
``enth_mass_sat_phase[phase]``            Saturation specific enthalpy of phase (J/kg)
``enth_mol_sat_phase[phase]``             Saturation molar enthalpy of phase (J/mol)
``entr_mass_sat_phase[phase]``            Saturation specific entropy of phase (J/kg/K)
``entr_mol_sat_phase[phase]``             Saturation molar entropy of phase (J/mol/K)
``energy_internal_mass_sat_phase[phase]`` Saturation specific internal energy of phase (J/kg)
``energy_internal_mol_sat_phase[phase]``  Saturation molar internal energy of phase (J/mol)
``volume_mass_sat_phase[phase]``          Saturation specific volume of phase (m\ :superscript:`3`/kg)
``volume_mol_sat_phase[phase]``           Saturation molar volume of phase (m\ :superscript:`3`/mol)
``dh_vap_mass``                           Specific enthalpy of vaporization at P (J/kg)
``dh_vap_mol``                            Molar enthalpy of vaporization at P (J/mol)
``ds_vap_mass``                           Specific entropy of vaporization at P (J/kg/K)
``ds_vap_mol``                            Molar entropy of vaporization at P (J/mol/K)
``du_vap_mass``                           Specific internal energy of vaporization at P (J/kg)
``du_vap_mol``                            Molar internal energy of vaporization at P (J/mol)
``phase_frac[phase]``                     Phase fraction (dimensionless) mole and mass fraction same for pure
``enthalpy_mass_phase[phase]``            Specific enthalpy of phase (J/kg) 
``enthalpy_mol_phase[phase]``             Molar enthalpy of phase (J/mol) 
``entropy_mass_phase[phase]``             Specific entropy of phase (J/kg/K) 
``entropy_mol_phase[phase]``              Molar entropy of phase (J/mol/K) 
``energy_internal_mass_phase[phase]``     Specific internal energy of phase (J/kg) 
``energy_internal_mol_phase[phase]``      Molar internal energy of phase (J/mol) 
``cp_mass_phase[phase]``                  Specific isobaric heat capacity for phase (J/kg)       
``cp_mol_phase[phase]``                   Molar isobaric heat capacity for phase (J/mol)      
``cv_mass_phase[phase]``                  Specific isochoric heat capacity for phase (J/kg)
``cv_mol_phase[phase]``                   Molar isochoric heat capacity for phase (J/mol)     
``speed_sound_phase[phase]``              Speed of sound in phase (m/s)
``vol_mass_phase[phase]``                 Specific volume of phase (m\ :superscript:`3`/kg)
``vol_mol_phase[phase]``                  Molar volume of phase (m\ :superscript:`3`/mol)
``dens_mass_phase[phase]``                Mass density of phase (kg/m\ :superscript:`3`)
``dens_mol_phase[phase]``                 Mole density of phase (mol/m\ :superscript:`3`)
``flow_mass_comp[comp]``                  Same as total mass flow since pure (kg/s)
``flow_mol_comp[comp]``                   Same as total mole flow since pure (mol/s)
``cp_mass``                               Specific isobaric heat capacity for mixed phase (J/kg)       
``cp_mol``                                Molar isobaric heat capacity for mixed phase (J/mol)      
``cv_mass``                               Specific isochoric heat capacity for mixed phase (J/kg)
``cv_mol``                                Molar isochoric heat capacity for mixed phase (J/mol)
``dens_mass``                             Mass density of mixed phase (kg/m\ :superscript:`3`)
``dens_mol``                              Mole density of mixed phase (mol/m\ :superscript:`3`)
``flow_vol``                              Total mixed phase volumetric flow (m\ :superscript:`3`/s)
``heat_capacity_ratio``                   Mixed phase cp/cv (dimensionless)
``visc_d_phase[phase]``                   Dynamic viscosity (Pa*s), depending on substance, may not be available
``visc_k_phase[phase]``                   Kinematic viscosity (m\ :superscript:`2`/s), depending on substance, may not be available
``therm_cond_phase[phase]``               Thermal conductivity of phase (W/m/s), depending on substance, may not be available
``surface_tension``                       Surface tension (N/m), depending on substance, may not be available
``P_under_sat``                           Pressure under saturation pressure (Pa)
``P_over_sat``                            Pressure over saturation pressure (Pa)
========================================= ===============================================================================================


Initialization
--------------

The Helmholtz EoS state blocks provide initialization functions for general
compatibility with the IDAES framework, but as long as the state variables are
specified to some reasonable value, initialization is not required. All required
solves are handled by external functions.

However, in order to support a general hierarchical initialization for unit models
which use Helmholtz equation of state properties, a custom ``Initializer`` for these
property packages is available.

.. module:: idaes.models.properties.general_helmholtz.helmholtz_state

.. autoclass:: HelmholtzEoSInitializer
   :members: initialize

