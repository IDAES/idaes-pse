Fixed Bed 0D Reactor
====================

The IDAES Fixed Bed 0D Reactor (FixedBed0D) model represents a unit operation where a gas stream passes through a 
solid phase bed.
The FixedBed0D mathematical model is a 0-D time variant model with two phases (gas and solid). The model captures the
gas-solid interaction between both phases through reaction, mass and heat transfer.

**Assumptions:**

* There is assumed to be no axial or radial variation in solid composition. 
* Excess gas flowrate and ideal gas behavior are assumed.
* The reactor is assumed to be isothermal, and the reaction is one-way (no equilibrium).

**Requirements:**

* Property packages contain temperature and pressure variables.

Degrees of Freedom
------------------

FixedBed0D Reactors generally have at least 2 (or more) degrees of freedom, including the reactor bed diameter and height.

Model Structure
---------------

The FixedBed0D model assumes a single control volume, and writes a set of material and energy balances for the solid phase. The gas phase conditions are considered to be fixed at the inlet values.

Construction Arguments
----------------------

The IDAES FixedBed0D model has construction arguments specific to the whole unit and to the individual regions.

**Arguments that are applicable to the FixedBed0D unit as a whole are:**

* dynamic - indicates whether the model will be dynamic or not, must be True (default = True).
* has_holdup - indicates whether holdup terms are constructed or not, must be True (default = True).
* energy_balance_type - indicates what type of energy balance should be constructed (default = enthalpyTotal):
  
  - EnergyBalanceType.none - exclude energy balances
  - EnergyBalanceType.enthalpyTotal - single enthalpy balance for material
  - EnergyBalanceType.enthalpyPhase - enthalpy balances for each phase
  - EnergyBalanceType.energyTotal - single energy balance for material
  - EnergyBalanceType.energyPhase - energy balances for each phase

**Arguments that are applicable to a specific region:**

* gas_property_package - property package to use when calculating properties for the gas phase (default = 'None'). 
  This is provided as a Physical Parameter Block by the Flowsheet when creating the model. 
* gas_property_package_args - set of arguments to be passed to the gas phase Property Blocks when they are created (default = 'None').
* solid_property_package - property package to use when calculating properties for the solid phase (default = 'None'). 
  This is provided as a Physical Parameter Block by the Flowsheet when creating the model. 
* solid_property_package_args - set of arguments to be passed to the solid phase Property Blocks when they are created (default = 'None').
* reaction_package - reaction package to use when calculating properties for solid phase reaction (default = None). 
  This is provided as a Reaction Parameter Block by the Flowsheet when creating the model. 
* reaction_package_args - set of arguments to be passed to the Reaction Blocks when they are created (default = None).

Constraints
-----------

In the following, the subscripts :math:`g` and :math:`s` refer to the gas and solid phases, respectively. 
FixedBed0D units write the following Constraints:

Reaction and Mass/Heat Transfer Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Volume of the bed:

.. math:: V_{bed} = \pi L_{bed} (0.5 D_{bed})^2

Volume of the solid:

.. math:: V_{s,t} = V_{bed} (1 - \varepsilon)

Solid material holdup:

.. math:: J_{mass,s,t,j} = V_{s} \rho_{mass,s,t} x_{mass,s,t,j}

Solid material accumulation:

.. math:: {\dot{J}}_{mass,s,t,j} = V_{s} {MW}_{s,j} {\Sigma}_{r} {r_{s,t,r} \nu_{s,j,r}}

Total mass of solids:

.. math:: M_{s,t} = V_{s} {\Sigma_{j} {\rho_{mass,s,t,j}}}

Sum of component mass fractions:

.. math:: {\Sigma_{j}} {x_{mass,s,t,j}} = 1 \space \space \space \forall \space t

if self.config.energy_balance_type != EnergyBalanceType.none:

    Solid energy holdup:

    .. math:: q_{energy,s,t} = V_{s} \rho_{mass,s,t} H_{mass,s,t}

    Solid energy accumulation:

    .. math:: {\dot{q}}_{energy,s,t} = - V_{s} {\Sigma}_{r} {r_{s,t,r} H_{rxn,s,t,r}}

if self.config.energy_balance_type == EnergyBalanceType.none:

    Isothermal solid phase:

    .. math:: T_{s,t} = T_{s,t=0}

List of Variables
-----------------

.. csv-table::
   :header: "Variable", "Description", "Reference to"

   ":math:`D_{bed}`", "Reactor bed diameter", "``bed_diameter``"
   ":math:`H_{mass,s,t}`", "Solid phase mass enthalpy", "``solids.enth_mass``"
   ":math:`H_{rxn,s,t,r}`", "Solid phase reaction enthalpy", "``solids.reactions.dh_rxn``"
   ":math:`J_{mass,s,t,j}`", "Material holdup, solid phase", "``solids.solids_material_holdup``"
   ":math:`{\dot{J}}_{mass,s,t,j}`", "Material accumulation, solid phase", "``solids.solids_material_accumulation``"
   ":math:`L_{bed}`", "Reactor bed height", "``bed_height``"
   ":math:`M_{s,t}`", "Total mass of solids", "``solids.mass_solids``"
   ":math:`q_{energy,s,t}`", "Energy holdup, solid phase", "``solids.energy_material_holdup``"
   ":math:`{\dot{q}}_{energy,s,t}`", "Energy accumulation, solid phase", "``solids.solids_energy_accumulation``"
   ":math:`r_{s,t,r}`", "Solid phase reaction rate of reaction r", "``solids.reactions.reaction_rate``"
   ":math:`T_{s,t}`", "Solid phase temperature", "``solids.temperature``"
   ":math:`x_{mass,s,t,j}`", "Mass fraction of component j, solid phase", "``solids.mass_frac_comp``"
   ":math:`V_{bed}`", "Total volume of the bed", "``volume_bed``"
   ":math:`V_{s,t}`", "Total volume of solids", "``solids.volume_solid``"
   "*Greek letters*", " ", " "
   ":math:`\rho_{mass,s,t}`", "Density of solid particles", "``solids.dens_mass_particle``"

List of Parameters
------------------

.. csv-table::
   :header: "Parameter", "Description", "Reference to"

   ":math:`\varepsilon`", "Reactor bed voidage", "``bed_voidage``"
   ":math:`{MW}_{s,j}`", "Molecular weight of solid component j", "``solids._params.mw_comp``"
   ":math:`\nu_{s,j,r}`", "Stoichiometric coefficients", "``solids.reaction_package.rate_reaction_stoichiometry``"

Initialization
--------------

The initialization method for this model will save the current state of the model before
commencing initialization and reloads it afterwards. The state of the model will be the same after
initialization, only the initial guesses for unfixed variables will be changed. For the FixedBed0D
model, the initialization method is a general purpose routine for simple unit models where a single
control volume is initialized followed by attempting to solve the entire unit.

The model allows for the passing of a dictionary of values of the state variables of the gas and
solid phases that can be used as initial guesses for the state variables throughout the time domain
of the model. This is optional but recommended. A typical guess could be values of the gas and solid
inlet port variables at time :math:`t=0`.

The model initialization proceeds through a sequential hierarchical method where the model
equations are deactivated at the start of the initialization routine, and the complexity of the model
is built up through activation and solution of various sub-model blocks and equations at each
initialization step. At each step the model variables are updated to better guesses obtained from
the model solution at that step.

The initialization routine proceeds as follows:

*  Step 1:  Initialize the thermo-physical model blocks.
*  Step 2:  Initialize the reaction properties.


FixedBed0D Class
----------------

.. module:: idaes.models_extra.gas_solid_contactors.unit_models.fixed_bed_0D

.. autoclass:: FixedBed0D
   :members:

FixedBed0DData Class
--------------------

.. autoclass:: FixedBed0DData
   :members:

