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

The FixedBed0D model equations are derived from a similar moving-bed model:

A. Ostace, A. Lee, C.O. Okoli, A.P. Burgard, D.C. Miller, D. Bhattacharyya, Mathematical modeling of a moving-bed reactor for chemical looping combustion of methane, in: M.R. Eden, M. Ierapetritou, G.P. Towler (Eds.),13th Int. Symp. Process Syst. Eng. (PSE 2018), Computer-Aided Chemical Engineering 2018, pp. 325–330 , San Diego, CA.

Degrees of Freedom
------------------

FixedBed0D Reactors generally have at least 0 (or more) degrees of freedom, consisting of design and operating variables. At a minimum, the inlet stream properties must be specified (no need to specify geometric or design parameters).

Model Structure
---------------

The FixedBed0D model does not use ControlVolume blocks, and instead writes a set of material and energy balances for the solid phase. The gas phase conditions are considered to be fixed at the inlet values.

Construction Arguments (need to update)
---------------------------------------

placeholder text

Constraints (need to update)
----------------------------

In the following, the subscripts :math:`g` and :math:`s` refer to the gas and solid phases, respectively. 
FixedBed0D units write the following Constraints:

Reaction and Mass/Heat Transfer Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Solid material holdup:

.. math:: h_{mass,s,t,j} = V_{s} \rho_{mass,s,t,j} x_{mass,s,t,j}

Solid material accumulation:

.. math:: {\dot{h}}_{mass,s,t,j} = V_{s} {MW}_{s,j} {\Sigma}_{r} {r_{s,t,r} \nu_{s,j,r}}

Total mass of solids

.. math:: S_{s,t} = V_{s} {\Sigma_{j} {\rho_{mass,s,t,j}}}

List of Variables (need to update)
----------------------------------

.. csv-table::
   :header: "Variable", "Description", "Reference to"

   ":math:`h_{mass,s,t,j}`", "Material holdup, solid phase", "``solids.solids_material_holdup``"
   ":math:`{\dot{h}}_{mass,s,t,j}`", "Material accumulation, solid phase", "``solids.solids_material_accumulation``"
   ":math:`r_{s,t,r}`", "Solid phase reaction rate of reaction r", "``solids.reactions.reaction_rate``"
   ":math:`S_{s,t}`", "Total mass of solids", "``solids.mass_solids``"
   ":math:`x_{mass,s,t,j}`", "Mass fraction of component j, solid phase", "``solids.mass_frac_comp``"
   ":math:`V_{s}`", "Total volume of solids", "``solids.volume_solid``"
   "*Greek letters*", " ", " "
   ":math:`\rho_{mass,s,t,j}`", "Density of solid particles of component j", "``solids.dens_mass_particle``"
   






   
   
   ":math:`h_{energy,s,t}`", "Energy holdup, solid phase", "``solid_phase.heat``"
   ":math:`a_{mass,s,t}`", "Material accumulation, solid phase", "``gas_phase.heat``"
   ":math:`a_{energy,s,t}`", "Energy accumulation, solid phase", "``solid_phase.heat``"
   
   ":math:`H_{rxn,s,t,r}`", "Solid phase reaction enthalpy", "``solid_phase.reactions.reaction_rate``"
   ":math:`T_{g,t}`", "Gas phase temperature", "``gas_phase.properties.temperature``"
   ":math:`T_{s,t}`", "Solid phase temperature", "``solid_phase.properties.temperature``"
   
   
   ":math:`\phi_{s,t}`", "Bed porosity", "``solid_phase.properties.porosity``"

List of Parameters (need to update)
-----------------------------------

.. csv-table::
   :header: "Parameter", "Description", "Reference to"

   ":math:`{MW}_{s,j}`", "Molecular weight of solid component j", "``solids._params.mw_comp``"
   ":math:`\nu_{s,j,r}`", "Stoichiometric coefficients", "``solids.reaction_package.rate_reaction_stoichiometry``"

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

*  Step 1:  Initialize the thermo-physical model blocks.
*  Step 2:  Initialize the reaction properties.


FixedBed0D Class
----------------

.. module:: idaes.gas_solid_contactors.unit_models.fixed_bed_0D

.. autoclass:: FixedBed0D
   :members:

FixedBed0DData Class
--------------------

.. autoclass:: FixedBed0DData
   :members:

