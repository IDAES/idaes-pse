Control Volume Classes
======================

.. toctree::
    :maxdepth: 1

    control_volume_0d
    control_volume_1d

Control Volumes are the center of the IDAES process modeling framework, and serve as the fundamental building block of all unit operations. Control Volumes represent a single, well-defined volume of material over which material, energy and/or momentum balances will be performed.

The IDAES Control Volume classes are designed to facilitate the construction of these balance equations by providing the model developer with a set of pre-built methods to perform the most common tasks in developing models of unit operations. The Control Volume classes contain methods for creating and linking the necessary property calculations and writing common forms of the balance equations so that the model developer can focus their time on the aspects that make each unit model unique.

The IDAES process modeling framework currently supports two types of Control Volume:

* ControlVolume0D represents a single well-mixed volume of material with a single inlet and a single outlet. This type of control volume is sufficient to model most inlet-outlet type unit operations which do not require spatial discretization.
* ControlVolume1D represents a volume with spatial variation in one dimension parallel to the material flow. This type of control volume is useful for representing flow in pipes and simple 1D flow reactors.

Common Control Volume Tasks
---------------------------

All of the IDAES Control Volume classes are built on a common core class which defines a set of common tasks required for all Control Volumes. The more specific Control Volume classes then build upon these common tasks to provide tools appropriate for their specific application.

All Control Volume classes begin with the following tasks:

* Determining if the Control Volume should be steady-state or dynamic, getting the time domain, and determining whether material and energy holdups should be calculated.
* Collecting the information necessary for creating StateBlocks and ReactionBlocks.
* Creating references to the phase_list and component_list Sets from the PhysicalParameterBlock.

More details on these steps is provided later.

Control Volume Configuration Arguments
--------------------------------------

All Control Volume Blocks have the following configuration arguments:

* dynamic - indicates whether this will be a dynamic (True) or steady-state (False) control volume.
* has_holdup - indicates whether material and energy holdup terms should be calculated. If dynamic = True, this must also be True, and a ConfigurationError will be returned if it is False.
* property_package - a pointer to the PhysicalParameterBlock instance to be used when constructing StateBlocks within this Control Volume.
* propery_package_args - a ConfigBlock of arguments to be passed to all StateBlocks as they are instantiated.
* reaction_package - a pointer to the ReactionParameterBlock instance to be used when constructing ReactionBlocks within this Control Volume (if necessary). This argument is only necessary when reaction terms are to be added to the balance equations.
* reaction_package_args - a ConfigBlock of arguments to be passed to all ReactionBlocks as they are instantiated.
* auto_construct - flag indicating whether the Control Volume should attempt to automatically construct a set of material, energy and momentum balances based on unit model configuration arguments (True), or if the unit model will explicitly call the methods to construct the Control Volume (False). Default value is False.


Setting up the time domain
--------------------------

The first common task the Control Volume block performs is to determine if it should be dynamic or steady-state and to collect the time domain from the UnitModel. Control Volume blocks have an argument `dynamic` which can be provided during construction which specifies if the Control Volume should be dynamic (dynamic = True) or steady-state (dynamic = False). If the argument is not provided, the Control Volume block will inherit this argument from its parent UnitModel.

After setting the dynamic argument, the Control Volume block then creates a reference to the time domain (time_ref) from the UnitModel. If the block containing the Control Volume block does not have an attribute named time a DynamicsError will be returned.

Finally, the Control Volume checks that the `has_holdup` argument is consistent with the `dynamic` argument, and raises a ConfigurationError if it is not.

Getting Property Package Information
-------------------------------------

If a reference to a property package was not provided by the UnitModel as an argument, the Control Volume first checks to see if the UnitModel has a property_package argument set, and uses this if present. Otherwise, the Control Volume block begins searching up the model tree looking for an argument named default_property_package and uses the first of these that it finds. If not default_property_package is found, a ConfigurationError is returned.

Collecting Indexing Sets for Property Package
---------------------------------------------

The final common step for all Control Volumes is to collect any required indexing sets from the physical property package (for example component and phase lists). These are used by the Control Volume for determining what balance equations need to be written, and what terms to create.

The indexing sets the Control Volume looks for are:

* component_list - used to determine what components are present, and thus what material balances are required
* phase_list - used to determine what phases are present, and thus what balance equations are required

ControlVolumeBase Class
-----------------------

A key purpose of Control Volumes is to automate as much of the task of writing a unit model as possible. For this purpose, Control Volumes support a number of methods for common tasks model developers may want to perform. The specifics of these methods will different between different types of Control Volumes, and certain methods may not be applicable to some types of Control Volumes (in which case a NotImplementedError will be returned). A full list of potential methods is provided here, however users should check the documentation for the specific Control Volume they are using for more details on what methods are supported in that specific Control Volume.

.. module:: idaes.core.control_volume_base

.. autoclass:: ControlVolumeBase
    :members:

Auto-Construct Method
---------------------

To reduce the demands on unit model developers even further, Control Volumes have an optional auto-construct feature that will attempt to populate the Control Volume based on a set of instructions provided at the Unit Model level. If the `auto_construct` configuration argument is set to True, the following methods are called automatically in the following order when instantiating the Control Volume.

* add_geometry
* add_state_blocks
* add_reaction_blocks
* add_meterial_balances
* add_energy_balances
* add_momentum_balances
* apply_transformation

To determine what terms are required for the balance equations, the Control Volume expects the Unit Model to have the following configuration argumnets, which are used as arguments to the methods above.

* dynamic
* has_holdup
* material_balance_type
* energy_balance_type
* momentum_balance_type
* has_rate_reactions
* has_equilibrium_reactions
* has_phase_equilibrium
* has_mass_transfer
* has_heat_of_reaction
* has_heat_transfer
* has_work_transfer
* has_pressure_change
* property_package
* property_package_args
* reaction_package
* reaction_package_args

For convenience, a template ConfigBlock (named CONIG_Template) is available in the control_volume_base.py module which contains all the necessary arguments which can be inherited by unit models wishingto use the auto-construct feature.
