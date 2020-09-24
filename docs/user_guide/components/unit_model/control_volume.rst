Control Volume
==============

.. contents:: :local:

Overview
--------

Control volumes serve as the 
fundamental building block of all unit operations. Control Volumes represent a single, 
well-defined volume of material over which material, energy and/or momentum balances will 
be performed.

The IDAES ControlVolume classes are designed to facilitate the construction of these balance 
equations by providing the model developer with a set of pre-built methods to perform the most 
common tasks in developing models of unit operations. The ControlVolume classes contain methods 
for creating and linking the necessary property calculations and writing common forms of the 
balance equations so that the model developer can focus their time on the aspects that make each 
unit model unique.

The IDAES process modeling framework currently supports two types of ControlVolumes:

* :ref:`ControlVolume0DBlock<technical_specs/core/control_volume_0d:0D Control Volume Class>` represents a single well-mixed volume of material with a single inlet and a single outlet. This type of control volume is sufficient to model most inlet-outlet type unit operations which do not require spatial discretization.
* :ref:`ControlVolume1DBlock<technical_specs/core/control_volume_1d:1D Control Volume Class>` represents a volume with spatial variation in one dimension parallel to the material flow. This type of control volume is useful for representing flow in pipes and simple 1D flow reactors.

Common Control Volume Tasks
---------------------------

All of the IDAES ControlVolume classes are built on a common core ControlVolumeBlockData which 
defines a set of common tasks required for all Control Volumes. The more specific ControlVolume classes 
then build upon these common tasks to provide tools appropriate for their specific application.

All ControlVolume classes begin with the following tasks:

* Determine if the ControlVolume should be steady-state or dynamic.
* Get the time domain.
* Determine whether material and energy holdups should be calculated.
* Collect information necessary for creating StateBlocks and ReactionBlocks.
* Create references to phase_list and component_list Sets in the PhysicalParameterBlock

Setting up the time domain
--------------------------

The first common task the ControlVolumeBlock performs is to determine if it should be dynamic 
or steady-state and to collect the time domain from the UnitModel. ControlVolumeBlocks have 
an argument `dynamic` which can be provided during construction which specifies if the 
Control Volume should be dynamic (`dynamic=True`) or steady-state (`dynamic=False`). If the 
argument is not provided, the ControlVolumeBlock will inherit this argument from its parent 
Unit model.

Finally, the ControlVolume checks that the `has_holdup` argument is consistent with the 
`dynamic` argument, and raises a `ConfigurationError` if it is not.

Getting Property Package Information
------------------------------------

If a reference to a property package was not provided by the UnitModel as an argument, 
the Control Volume first checks to see if the unit model has a `property_package` argument 
set, and uses this if present. Otherwise, the ControlVolumeBlock begins searching up the model 
tree looking for an argument named `default_property_package` and uses the first of these 
that it finds. If no `default_property_package` is found, a `ConfigurationError` is returned.

Collecting Indexing Sets for Property Package
---------------------------------------------

The final common step for all ControlVolumes is to collect any required indexing sets from the physical property package (for example component and phase lists). These are used by the Control Volume for determining what balance equations need to be written, and what terms to create.

The indexing sets the ControlVolume looks for are:

* `component_list` - used to determine what components are present, and thus what material balances are required
* `phase_list` - used to determine what phases are present, and thus what balance equations are required

ControlVolume and ControlVolumeBlockData Classes
------------------------------------------------

A key purpose of ControlVolumes is to automate as much of the task of writing a unit model as 
possible. For this purpose, ControlVolumes support a number of methods for common tasks model 
developers may want to perform. The specifics of these methods will be different between 
different types of ControlVolumes, and certain methods may not be applicable to some types of 
Control Volumes (in which case a `NotImplementedError` will be returned). A full list of 
potential methods is provided here, however users should check the documentation for the 
specific Control Volume they are using for more details on what methods are supported in that 
specific Control Volume.

A key feature of the IDAES Core Modeling Framework is the use of ControlVolumeBlocks. ControlVolumes 
represent a volume of material over which material, energy and/or momentum balances 
can be performed. ControlVolumeBlocks contain methods to automate the task of writing common 
forms of these balance equations. ControlVolumeBlocks can also automate the creation of 
StateBlocks and ReactionBlocks associated with the control volume.


