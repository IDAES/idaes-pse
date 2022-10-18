Fixed Bed 1D Reactor
====================

The IDAES Fixed Bed 1D Reactor (FixedBed1D) model represents a unit operation where a gas stream  
passes through a solid phase bed in a linear reactor vessel.
The FixedBed1D mathematical model is a 1-D time variant model with two phases (gas and solid). The model captures the
gas-solid interaction between both phases through reaction/adsorption, mass and heat transfer. 

**Assumptions:**

* The radial concentration and temperature gradients are assumed to be negligible. 
* The reactor is assumed to be adiabatic.

**Requirements:**

* Property package contains temperature and pressure variables.

Degrees of Freedom
------------------

FixedBed1D Reactors generally have at least 2 (or more) degrees of freedom, consisting of design and operating variables. The design variables of reactor length and diameter are typically the minimum variables to be fixed.

Model Structure
---------------

The core FixedBed1D unit model consists of one ControlVolume1DBlock Block named gas_phase, which has one 
Inlet Port (named gas_inlet) and one Outlet Port (named gas_outlet). As there are no flow variables in the solid_phase 
a ControlVolume1DBlock Block is not used, instead indexed state blocks are used to access the solid thermophysical properties whilst the mass and energy balance equations
of the solid are written in the core FixedBed1D unit model.

Construction Arguments
----------------------

The IDAES FixedBed1D model has construction arguments specific to the whole unit and to the individual regions.

**Arguments that are applicable to the FixedBed1D unit as a whole are:**

* finite_elements - sets the number of finite elements to use when discretizing the spatial domains (default = 10).
* length_domain_set - sets the list of point to use to initialize a new ContinuousSet (default = [0.0, 1.0]).
* transformation_method - sets the discretization method to use by the Pyomo TransformationFactory 
  to transform the spatial domain (default = dae.finite_difference):
  
      - dae.finite_difference - finite difference method.
      - dae.collocation - orthogonal collocation method.
  
* transformation_scheme - sets the scheme to use when transforming a domain. 
  Selected schemes should be compatible with the transformation_method chosen (default = None):
  
      - None - defaults to "BACKWARD" for finite difference transformation method and to "LAGRANGE-RADAU" for collocation transformation method
      - BACKWARD - use a finite difference transformation method.
      - FORWARD - use a finite difference transformation method.
      - LAGRANGE-RADAU - use a collocation transformation method.   
  
* collocation_points - sets the number of collocation points to use when discretizing the spatial domains (default = 3, collocation methods only).

* flow_type - indicates the flow arrangement within the unit to be modeled. Options are:

      - 'forward_flow' - (default) gas flows in the forward direction (from x=0 to x=1)
      - 'reverse_flow' - gas flows in the reverse direction (from x=1 to x=0).

* material_balance_type - indicates what type of energy balance should be constructed (default = MaterialBalanceType.componentTotal).

      - MaterialBalanceType.componentTotal - use total component balances.
      - MaterialBalanceType.total - use total material balance.
    
* energy_balance_type - indicates what type of energy balance should be constructed (default = EnergyBalanceType.enthalpyTotal).

      - EnergyBalanceType.none - excludes energy balances.
      - EnergyBalanceType.enthalpyTotal - single enthalpy balance for material.

* momentum_balance_type - indicates what type of momentum balance should be constructed (default = MomentumBalanceType.pressureTotal).

      - MomentumBalanceType.none - exclude momentum balances.
      - MomentumBalanceType.pressureTotal - single pressure balance for material.

* has_pressure_change - indicates whether terms for pressure change should be constructed (default = True).

      - True - include pressure change terms.
      - False - exclude pressure change terms.

* pressure_drop_type - indicates what type of pressure drop correlation should be used (default = "ergun_correlation").

      - "ergun_correlation" - use the Ergun equation.
      - "simple_correlation" - use a simplified pressure drop correlation.

**Arguments that are applicable to the gas phase:**

* property_package - property package to use when constructing gas phase Property Blocks (default = 'use_parent_value'). 
  This is provided as a Physical Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the gas phase Property Blocks when they are created (default = 'use_parent_value').
* reaction_package - reaction package to use when constructing gas phase Reaction Blocks (default = None). 
  This is provided as a Reaction Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* reaction_package_args - set of arguments to be passed to the gas phase Reaction Blocks when they are created (default = None).
* has_equilibrium_reactions - sets flag to indicate if terms of equilibrium controlled reactions should be constructed (default = False).

**Arguments that are applicable to the solid phase:**

* property_package - property package to use when constructing solid phase Property Blocks (default = 'use_parent_value'). 
  This is provided as a Physical Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the Indexed State Blocks will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the solid phase Property Blocks when they are created (default = 'use_parent_value').
* reaction_package - reaction package to use when constructing solid phase Reaction Blocks (default = None). 
  This is provided as a Reaction Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the Indexed State Blocks will try to use the default property package if one is defined.
* reaction_package_args - set of arguments to be passed to the solid phase Reaction Blocks when they are created (default = None).
* has_equilibrium_reactions - sets flag to indicate if terms of equilibrium controlled reactions should be constructed (default = False).

Additionally, FixedBed1D units have the following construction arguments which are passed to all the ControlVolume1DBlock and state Blocks
and are always specified to True.

========================= ==========================
Argument                  Value
========================= ==========================
dynamic                   True
has_holdup                True
========================= ==========================

Constraints
-----------

In the following, the subscripts :math:`g` and :math:`s` refer to the gas and solid phases, respectively. 
In addition to the constraints written by the control_volume Block, FixedBed1D units write the following Constraints:

Geometry Constraints
^^^^^^^^^^^^^^^^^^^^

**Area of the reactor bed:**

.. math:: A_{bed} = \pi \left( \frac{ D_{bed} }{ 2 } \right)^2

**Area of the gas domain:**

.. math:: A_{g,t,x} = \varepsilon A_{bed} 

**Area of the solid domain:**

.. math:: A_{s,t,x} = (1 - \varepsilon) A_{bed}

Hydrodynamic Constraints
^^^^^^^^^^^^^^^^^^^^^^^^

**Superficial velocity of the gas:**

.. math:: u_{g,t,x} = \frac{ F_{mol,g,t,x} }{ A_{bed} \rho_{mol,g,t,x} }

**Pressure drop:**

The constraints written by the FixedBed1D model to compute the pressure drop (if `has_pressure_change` is 'True') in the reactor depend upon the 
construction arguments chosen:

If `pressure_drop_type` is `simple_correlation`:

.. math:: - \frac{ dP_{g,t,x} }{ dx } = 0.2 \left( \rho_{mass,s,t,x} - \rho_{mass,g,t,x} \right) u_{g,t,x}

If `pressure_drop_type` is `ergun_correlation`:

.. math:: - \frac{ dP_{g,t,x} }{ dx } = \frac{ 150 \mu_{g,t,x} {\left( 1 - \varepsilon \right)}^{2} u_{g,t,x} }{ \varepsilon^{3} d_{p}^2 } + \frac{ 1.75 \left( 1 - \varepsilon \right) \rho_{mass,g,t,x} u_{g,t,x}^{2} }{ \varepsilon^{3} d_{p} }

Reaction Constraints
^^^^^^^^^^^^^^^^^^^^

If `gas_phase_config.reaction_package` is not 'None':

    **Gas phase reaction extent:**

    .. math:: \xi_{g,t,x,r} = r_{g,t,x,r} A_{g,t,x}

If `solid_phase_config.reaction_package` is not 'None':

    **Gas phase heterogeneous rate generation/consumption:** 

    .. math:: M_{g,t,x,p,j} = A_{s,t,x} \sum_{r}^{reactions} {\nu_{s,j,r} r_{s,t,x,r}}

Dimensionless numbers, mass and heat transfer coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Sum of component mass fractions in solid phase:**

.. math:: {\Sigma_{j}} {x_{mass,s,t,x,j}} = 1 \space \space \space \forall \space t \forall \space x

**Solid material holdup:**

.. math:: J_{mass,s,t,x,j} = A_{s,t,x} \rho_{mass,s,t,x} x_{mass,s,t,x,j}

If `solid_phase_config.reaction_package` is not 'None':

    **Solid material accumulation:**

    .. math:: {\dot{J}}_{mass,s,t,x,j} = A_{s,t,x} {MW}_{s,j} {\Sigma}_{r} {r_{s,t,x,r} \nu_{s,j,r}}

If `solid_phase_config.reaction_package` is 'None':

    **Solid material accumulation:**

    .. math:: {\dot{J}}_{mass,s,t,x,j} = 0

If `energy_balance_type` not `EnergyBalanceType.none`:

    **Solid energy holdup:**

    .. math:: q_{energy,s,t,x} = A_{s,t,x} \rho_{mass,s,t,x} H_{mass,s,t,x}

    **Solid energy accumulation:**

    .. math:: {\dot{q}}_{energy,s,t,x} = H_{s,t,x} - A_{s,t,x} {\Sigma}_{r} {r_{s,t,x,r} H_{rxn,s,t,x,r}}

If `energy_balance_type` is `EnergyBalanceType.none`:

    **Isothermal solid phase:**

    .. math:: T_{s,t,x} = T_{s,t=0,x}

    **Isothermal gas phase:**

    .. math:: T_{g,t,x} = T_{s,t=0,x}

If `energy_balance_type` not `EnergyBalanceType.none`:

    **Particle Reynolds number:**

    .. math:: Re_{p,t,x} = \frac{ u_{g,t,x} \rho_{mass,g,t,x} }{ \mu_{g,t,x} d_{p}}

    **Prandtl number:**

    .. math:: Pr_{t,x} = \frac{ c_{p,t,x} \mu_{g,t,x} }{ k_{g,t,x} }

    **Particle Nusselt number:**

    .. math:: Nu_{p,t,x} = 2 + 1.1 Pr_{t,x}^{1/3} \left| Re_{p,t,x} \right|^{0.6}

    **Particle to fluid heat transfer coefficient**

    .. math:: h_{gs,t,x} d_{p} = Nu_{p,t,x} k_{g,t,x}

    **Gas phase - gas to solid heat transfer:**

    .. math:: H_{g,t,x} = - \frac{ 6 } { d_{p} } h_{gs,t,x} \left( T_{g,t,x} - T_{s,t,x} \right) A_{s,t,x}

    **Solid phase - gas to solid heat transfer:**

    .. math:: H_{s,t,x} = \frac{ 6 } { d_{p} } h_{gs,t,x} \left( T_{g,t,x} - T_{s,t,x} \right) A_{s,t,x}

List of Variables
-----------------

.. csv-table::
   :header: "Variable", "Description", "Reference to"

   ":math:`A_{bed}`", "Reactor bed cross-sectional area", "``bed_area``"
   ":math:`A_{g,t,x}`", "Gas phase area (interstitial cross-sectional area)", "``gas_phase.area``"
   ":math:`A_{s,t,x}`", "Solid phase area", "``solid_phase.area``"
   ":math:`c_{p,t,x}`", "Gas phase heat capacity (constant :math:`P`)", "``gas_phase.properties.cp_mass``"
   ":math:`D_{bed}`", "Reactor bed diameter", "``bed_diameter``"
   ":math:`F_{mol,g,t,x}`", "Total molar flow rate of gas", "``gas_phase.properties.flow_mol``"
   ":math:`H_{g,t,x}`", "Gas to solid heat transfer term, gas phase", "``gas_phase.heat``"
   ":math:`H_{s,t,x}`", "Gas to solid heat transfer term, solid phase", "``solid_phase.heat``"
   ":math:`H_{rxn,s,t,x,r}`", "Solid phase heat of reaction", "``solids.reactions.dh_rxn``"
   ":math:`h_{gs,t,x}`", "Gas-solid heat transfer coefficient", "``gas_solid_htc``"
   ":math:`J_{mass,s,t,x,j}`", "Material holdup, solid phase", "``solids.solids_material_holdup``"
   ":math:`{\dot{J}}_{mass,s,t,x,j}`", "Material accumulation, solid phase", "``solids.solids_material_accumulation``"
   ":math:`k_{g,t,x}`", "Gas thermal conductivity", "``gas_phase.properties.therm_cond``"
   ":math:`L_{bed}`", "Reactor bed height", "``bed_height``"
   ":math:`M_{g,t,x,p,j}`", "Rate generation/consumption term, gas phase", "``gas_phase.mass_transfer_term``"
   ":math:`Nu_{p,t,x}`", "Particle Nusselt number", "``Nu_particle``"
   ":math:`dP_{g,t,x}`", "Total pressure derivative w.r.t. :math:`x` (axial position)", "``gas_phase.deltaP``"
   ":math:`Pr_{t,x}`", "Prandtl number", "``Pr``"
   ":math:`q_{energy,s,t,x}`", "Energy holdup, solid phase", "``solids.energy_material_holdup``"
   ":math:`{\dot{q}}_{energy,s,t,x}`", "Energy accumulation, solid phase", "``solids.solids_energy_accumulation``"
   ":math:`r_{g,t,x,r}`", "Gas phase reaction rate", "``gas_phase.reactions.reaction_rate``"
   ":math:`r_{s,t,x,r}`", "Solid phase reaction rate", "``solid_phase.reactions.reaction_rate``"
   ":math:`Re_{p,t,x}`", "Particle Reynolds number", "``Re_particle``"
   ":math:`T_{g,t,x}`", "Gas phase temperature", "``gas_phase.properties.temperature``"
   ":math:`T_{s,t,x}`", "Solid phase temperature", "``solid_phase.properties.temperature``"
   ":math:`u_{g,t,x}`", "Superficial velocity of the gas", "``velocity_superficial_gas``"
   "*Greek letters*", " ", " "
   ":math:`\varepsilon`", "Reactor bed voidage", "``bed_voidage``"
   ":math:`\mu_{g,t,x}`", "Dynamic viscosity of gas mixture", "``gas_phase.properties.visc_d``"
   ":math:`\xi_{g,t,x,r}`", "Gas phase reaction extent", "``gas_phase.rate_reaction_extent``"
   ":math:`\rho_{mass,g,t,inlet}`", "Density of gas mixture", "``gas_phase.properties.dens_mass``"
   ":math:`\rho_{mass,s,t,inlet}`", "Density of solid particles", "``solid_phase.properties.dens_mass_particle``"
   ":math:`\rho_{mol,g,t,x}`", "Molar density of the gas", "``gas_phase.properties.dens_mole``"

List of Parameters
------------------

.. csv-table::
   :header: "Parameter", "Description", "Reference to"

   ":math:`d_{p}`", "Solid particle diameter", "``solid_phase.properties._params.particle_dia``"
   ":math:`\nu_{s,j,r}`", "Stoichiometric coefficients", "``solid_phase.reactions.rate_reaction_stoichiometry``"


Initialization
--------------

The initialization method for this model will save the current state of the model before
commencing initialization and reloads it afterwards. The state of the model will be the same after
initialization, only the initial guesses for unfixed variables will be changed.

The model allows for the passing of a dictionary of values of the state variables of the gas and
solid phases that can be used as initial guesses for the state variables throughout the time and
spatial domains of the model. This is optional but recommended. A typical guess could be values
of the gas and solid inlet port variables at time :math:`t=0`.

The model initialization proceeds through a sequential hierarchical method where the model
equations are deactivated at the start of the initialization routine, and the complexity of the model
is built up through activation and solution of various sub-model blocks and equations at each
initialization step. At each step the model variables are updated to better guesses obtained from
the model solution at that step.

The initialization routine proceeds as follows:

*  Step 1:  Initialize the thermo-physical and transport properties model blocks.
*  Step 2:  Initialize the hydrodynamic properties.
*  Step 3a: Initialize mass balances without reactions and pressure drop.
*  Step 3b: Initialize mass balances with reactions and without pressure drop.
*  Step 3c: Initialize mass balances with reactions and pressure drop.
*  Step 4:  Initialize energy balances.


Block Triangularization Initialization
--------------------------------------

The Block Triangularization initialization routine initializes the 1DFixedBed model in the following steps:

*  Step 1:  Convert the system of equations from a partial differential set of equations (PDAE) to an algebraic set of equations (AE) by deactivating
    the discretization equations and sum_component_equations (if present), and fixing the state variables to a guess (typically inlet conditions).
*  Step 2:  Decompose the AE into strongly connected components via block triangularization and solve individual blocks.
*  Step 3:  Revert the fixed variables and deactivated constraints to their original states.

The routine allows for the passing of a dictionary of values of the state variables of the gas and
solid phases that can be used as initial guesses for the state variables throughout the time and
spatial domains of the model. This is optional but recommended. A typical guess could be values
of the gas and solid state variables at time :math:`t=0`.


FixedBed1D Class
----------------

.. module:: idaes.models_extra.gas_solid_contactors.unit_models.fixed_bed_1D

.. autoclass:: FixedBed1D
   :members:

FixedBed1DData Class
--------------------

.. autoclass:: FixedBed1DData
   :members:

