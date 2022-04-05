Flash Unit
==========

The IDAES Flash model represents a unit operation where a single stream undergoes a flash separation into two phases. The Flash model supports multiple types of flash operations, including pressure changes and addition or removal of heat.

Degrees of Freedom
------------------

Flash units generally have 2 degrees of freedom.

Typical fixed variables are:

* heat duty or outlet temperature (see note),
* pressure change or outlet pressure.

Note: When setting the outlet temperature of a Flash unit, it is best to set control_volume.properties_out[t].temperature. Setting the temperature in one of the outlet streams directly results in a much harder problem to solve, and may be degenerate or unbounded in some cases.

Model Structure
---------------

The core Flash unit model consists of a single ControlVolume0DBlock (named control_volume) with one Inlet Port (named inlet) connected to a Separator unit model with two outlet Ports named 'vap_outlet' and 'liq_outlet'. The Flash model utilizes the separator unit model in IDAES to split the outlets by phase flows to the liquid and vapor outlets respectively.

The Separator unit model supports both direct splitting of state variables and writing of full splitting constraints via the `ideal_separation` construction argument. Full details on the Separator unit model can be found in the documentation for that unit. To support direct splitting, the property package must use one of a specified set of state variables and support a certain set of property calacuations, as outlined in the table below.

==================================== ===================================
State Variables                      Required Properties
==================================== ===================================
**Material flow and composition**
------------------------------------------------------------------------
flow_mol  &  mole_frac               flow_mol_phase  &  mole_frac_phase
flow_mol_phase  &  mole_frac_phase   flow_mol_phase  &  mole_frac_phase
flow_mol_comp                        flow_mol_phase_comp
flow_mol_phase_comp                  flow_mol_phase_comp
flow_mass  &  mass_frac              flow_mass_phase  &  mass_frac_phase
flow_mass_phase  &  mass_frac_phase  flow_mass_phase  &  mass_frac_phase
flow_mass_comp                       flow_mass_phase_comp
flow_mass_phase_comp                 flow_mass_phase_comp
**Energy state**
------------------------------------------------------------------------
temperature                          temperature
enth_mol                             enth_mol_phase
enth_mol_phase                       enth_mol_phase
enth_mass                            enth_mass_phase
enth_mass_phase                      enth_mass_phase
**Pressure state**
------------------------------------------------------------------------
pressure                             pressure
==================================== ===================================

Construction Arguments
----------------------

Flash units have the following construction arguments:

* property_package - property package to use when constructing Property Blocks (default = 'use_parent_value'). This is provided as a Property Parameter Block by the Flowsheet when creating the model. If a value is not provided, the Holdup Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the Property Blocks when they are created.

Additionally, Flash units have the following construction arguments which are passed to the Holdup Block for determining which terms to construct in the balance equations.

========================= =================
Argument                  Default Value
========================= =================
dynamic                   False
include_holdup            False
material_balance_type     MaterialBalanceType.componentPhase
energy_balance_type       EnergyBalanceType.enthalpyTotal
momentum_balance_type     MomentumBalanceType.pressureTotal
has_phase_equilibrium     True
has_heat_transfer         True
has_pressure_change       True
========================= =================

Finally, Flash units also have the following arguments which are passed to the Separator block for determining how to split to two-phase mixture.

=================== =================
Argument            Default Value
=================== =================
ideal_separation    True
energy_split_basis  EnergySplittingType.equal_temperature
=================== =================

Additional Constraints
----------------------

Flash units write no additional Constraints beyond those written by the ControlVolume0DBlock and the Separator block.

Variables
---------

============== ==================================================
Name           Notes
============== ==================================================
heat_duty      Reference to control_volume.heat
deltaP         Reference to control_volume.deltaP
============== ==================================================

Flash Class
-----------

.. module:: idaes.models.unit_models.flash

.. autoclass:: Flash
    :members:

FlashData Class
---------------

.. autoclass:: FlashData
    :members:
