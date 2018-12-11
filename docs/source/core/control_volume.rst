Control Volumes
===============

.. toctree::
    :maxdepth: 1

    control_volume_0d
    control_volume_1d

Control Volumes are the center of the IDAES process modeling framework, and serve as the fundamental building block of all unit operations. Control Volumes represent a single, well-defined volume of material over which material, energy and/or momentum balances will be performed.

The IDAES Control Volume classes are designed to facilitate the construction of these balance equations by providing the model developer with a set of pre-built methods to perform the most common tasks in developing models of unit operations. The Control Volume classes contain methods for creating and linking the necessary property calculations and writing common forms of the balance equations so that the model developer can focus their time on the aspects that make each unit model unique.

The IDAES process modeling framework currently supports two types of Control Volume:

* ControlVolume0D represents a single well-mixed volume of material with a single inlet and a single outlet. This type of control volume is sufficient to model most inlet-outlet type unit operations which do not require spatial discretization.
* ControlVolume1D represents a volume with spatial variation in one dimension parallel to the material flow. This type of control volume is useful for representing flow in pipes and simple 1D flow reactors.

ControlVolumeBase Class
-----------------------

.. module:: idaes.core.control_volume_base

.. autoclass:: ControlVolumeBase
    :members:

