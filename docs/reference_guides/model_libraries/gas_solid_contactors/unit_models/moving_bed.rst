Moving Bed Reactor
==================

The IDAES Moving Bed Reactor (MBR) model represents a unit operation where two material streams – a 
solid phase and a gas phase – pass through a linear reactor vessel while undergoing chemical reaction(s). 
The two streams have opposite flow directions (counter-flow).
The MBR mathematical model is a 1-D rigorous first-principles model consisting of a set of differential 
equations obtained by applying the mass, energy (for each phase) and momentum balance equations. 

**Assumptions:**

* The radial concentration and temperature gradients are assumed to be negligible. 
* The reactor is assumed to be adiabatic. 
* The solid phase is assumed to be moving at a constant velocity determined by the solids feed rate to the reactor.

**Requirements:**

* Property package contains temperature and pressure variables.
* Property package contains minimum fluidization velocity.

The MBR model is based on:

A. Ostace, A. Lee, C.O. Okoli, A.P. Burgard, D.C. Miller, D. Bhattacharyya, Mathematical modeling of a moving-bed reactor for chemical looping combustion of methane, in: M.R. Eden, M. Ierapetritou, G.P. Towler (Eds.),13th Int. Symp. Process Syst. Eng. (PSE 2018), Computer-Aided Chemical Engineering 2018, pp. 325–330 , San Diego, CA.

Degrees of Freedom
------------------

MBRs generally have at least 2 (or more) degrees of freedom, consisting of design and operating variables. The design variables of reactor length and diameter are typically the minimum variables to be fixed.

Model Structure
---------------

The core MBR unit model consists of two ControlVolume1DBlock Blocks (named gas_phase and solid_phase), each with one 
Inlet Port (named gas_inlet and solid_inlet) and one Outlet Port (named gas_outlet and solid_outlet).

Discretization
--------------

The default spatial discretization is a backward finite difference with respect
to increasing :math:`x` coordinate. For dynamic operation, we recommend a
discretization that is backward with respect to direction of flow, i.e.
backward for the gas phase and forward for the solid phase. Discretizations may
be set for gas and solid phases individually using the
``gas_transformation_scheme`` and ``solid_transformation_scheme`` config options.

Constraints
-----------

In the following, the subscripts :math:`g` and :math:`s` refer to the gas and solid phases, respectively. 
In addition to the constraints written by the control_volume Block, MBR units write the following Constraints:

Geometry Constraints
^^^^^^^^^^^^^^^^^^^^

**Area of the reactor bed:**

.. math:: A_{bed} = \pi \left( \frac{ D_{bed} }{ 2 } \right)^2

**Area of the gas domain:**

.. math:: A_{g,t,x} = \varepsilon A_{bed} 

**Area of the solid domain:**

.. math:: A_{s,t,x} = (1 - \varepsilon) A_{bed}

**Length of the gas domain:**

.. math:: L_{g} = L_{bed}

**Length of the solid domain:**

.. math:: L_{s} = L_{bed}

Hydrodynamic Constraints
^^^^^^^^^^^^^^^^^^^^^^^^

**Superficial velocity of the gas:**

.. math:: u_{g,t,x} = \frac{ F_{mol,g,t,x} }{ A_{bed} \rho_{mol,g,t,x} }

**Superficial velocity of the solids:**

.. math:: u_{s,t} = \frac{ F_{mass,s,t,inlet} }{ A_{bed} \rho_{mass,s,t,inlet} }

**Pressure drop:**

The constraints written by the MBR model to compute the pressure drop (if `has_pressure_change` is 'True') in the reactor depend upon the 
construction arguments chosen:

If `pressure_drop_type` is `simple_correlation`:

.. math:: - \frac{ dP_{g,t,x} }{ dx } = 0.2 \left( \rho_{mass,s,t,x} - \rho_{mass,g,t,x} \right) u_{g,t,x}

If `pressure_drop_type` is `ergun_correlation`:

.. math:: - \frac{ dP_{g,t,x} }{ dx } = \frac{ 150 \mu_{g,t,x} {\left( 1 - \varepsilon \right)}^{2} \left( u_{g,t,x} + u_{s,t} \right) }{ \varepsilon^{3} d_{p}^2 } + \frac{ 1.75 \left( 1 - \varepsilon \right) \rho_{mass,g,t,x} \left( u_{g,t,x} + u_{s,t} \right)^{2} }{ \varepsilon^{3} d_{p} }

Reaction Constraints
^^^^^^^^^^^^^^^^^^^^

**Gas phase reaction extent:**

If `gas_phase_config.reaction_package` is not 'None':

.. math:: \xi_{g,t,x,r} = r_{g,t,x,r} A_{g,t,x}

**Solid phase reaction extent:**

If `solid_phase_config.reaction_package` is not 'None':

.. math:: \xi_{s,t,x,r} = r_{s,t,x,r} A_{s,t,x}

**Gas phase heterogeneous rate generation/consumption:** 

.. math:: M_{g,t,x,p,j} = A_{s,t,x} \sum_{r}^{reactions} {\nu_{s,j,r} r_{s,t,x,r}}

Dimensionless numbers, mass and heat transfer coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Particle Reynolds number:**

.. math:: Re_{p,t,x} = \frac{ u_{g,t,x} \rho_{mass,g,t,x} }{ \mu_{g,t,x} d_{p}}

**Prandtl number:**

.. math:: Pr_{t,x} = \frac{ c_{p,t,x} \mu_{g,t,x} }{ k_{g,t,x} }

**Particle Nusselt number:**

.. math:: Nu_{p,t,x} = 2 + 1.1 Pr_{t,x}^{1/3} \left| Re_{p,t,x} \right|^{0.6}

**Particle to fluid heat transfer coefficient**

.. math:: h_{gs,t,x} d_{p} = Nu_{p,t,x} k_{g,t,x}


If `energy_balance_type` not `EnergyBalanceType.none`:

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
   ":math:`F_{mass,s,t,inlet}`", "Total mass flow rate of solids, at inlet (:math:`x=1`)", "``solid_phase.properties.flow_mass``"
   ":math:`F_{mol,g,t,x}`", "Total molar flow rate of gas", "``gas_phase.properties.flow_mol``"
   ":math:`H_{g,t,x}`", "Gas to solid heat transfer term, gas phase", "``gas_phase.heat``"
   ":math:`H_{s,t,x}`", "Gas to solid heat transfer term, solid phase", "``solid_phase.heat``"
   ":math:`h_{gs,t,x}`", "Gas-solid heat transfer coefficient", "``gas_solid_htc``"
   ":math:`k_{g,t,x}`", "Gas thermal conductivity", "``gas_phase.properties.therm_cond``"
   ":math:`L_{bed}`", "Reactor bed height", "``bed_height``"
   ":math:`L_{g}`", "Gas domain length", "``gas_phase.length``"
   ":math:`L_{s}`", "Solid domain length", "``solid_phase.length``"
   ":math:`M_{g,t,x,p,j}`", "Rate generation/consumption term, gas phase", "``gas_phase.mass_transfer_term``"
   ":math:`Nu_{p,t,x}`", "Particle Nusselt number", "``Nu_particle``"
   ":math:`dP_{g,t,x}`", "Total pressure derivative w.r.t. :math:`x` (axial position)", "``gas_phase.deltaP``"
   ":math:`Pr_{t,x}`", "Prandtl number", "``Pr``"
   ":math:`r_{g,t,x,r}`", "Gas phase reaction rate", "``gas_phase.reactions.reaction_rate``"
   ":math:`r_{s,t,x,r}`", "Solid phase reaction rate", "``solid_phase.reactions.reaction_rate``"
   ":math:`Re_{p,t,x}`", "Particle Reynolds number", "``Re_particle``"
   ":math:`T_{g,t,x}`", "Gas phase temperature", "``gas_phase.properties.temperature``"
   ":math:`T_{s,t,x}`", "Solid phase temperature", "``solid_phase.properties.temperature``"
   ":math:`u_{g,t,x}`", "Superficial velocity of the gas", "``velocity_superficial_gas``"
   ":math:`u_{s,t}`", "Superficial velocity of the solids", "``velocity_superficial_solid``"
   "*Greek letters*", " ", " "
   ":math:`\varepsilon`", "Reactor bed voidage", "``bed_voidage``"
   ":math:`\mu_{g,t,x}`", "Dynamic viscosity of gas mixture", "``gas_phase.properties.visc_d``"
   ":math:`\xi_{g,t,x,r}`", "Gas phase reaction extent", "``gas_phase.rate_reaction_extent``"
   ":math:`\xi_{s,t,x,r}`", "Solid phase reaction extent", "``solid_phase.rate_reaction_extent``"
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

The initialization routine proceeds in as follows:

*  Step 1:  Initialize the thermo-physical and transport properties model blocks.
*  Step 2:  Initialize the hydrodynamic properties.
*  Step 3a: Initialize mass balances without reactions and pressure drop.
*  Step 3b: Initialize mass balances with reactions and without pressure drop.
*  Step 3c: Initialize mass balances with reactions and pressure drop.
*  Step 4:  Initialize energy balances.


MBR Class
---------

.. module:: idaes.models_extra.gas_solid_contactors.unit_models.moving_bed

.. autoclass:: MBR
   :members:

MBRData Class
-------------

.. autoclass:: MBRData
   :members:

