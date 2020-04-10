Phase Objects
=============

.. contents:: Contents 
    :depth: 2

Introduction
------------

Phase objects are used in the IDAES Process Modeling Framework to identify the thermodynamic phases of interest in a property package and to contain information describing the behavior of that phase (for example the equation of state which describes that phase).

The IDAES Process Modeling Framework currently supports the following types of phases, along with a generic Phase object.

* LiquidPhase
* SolidPhase
* VaporPhase

Phase Objects
-------------

In a number of unit operations, different phases behave in different ways. For example, in a Flash operation, the vapor phase exits through the top outlet whilst liquid phase(s) (and any solids) exit through the bottom outlet. In order to determine how a given phase should behave in these situations, each Phase object implements the following three methods:

* `is_liquid_phase()`
* `is_solid_phase()`
* `is_vapor_phase()`

These methods return a boolean (`True` or `False`) indicating whether the unit operation should treat the phase as being of the specified type in order to decide on how it should behave. Each type of phase returns `True` for its type and `False` for all other types (e.g. LiquidPhase returns `True` for `is_liquid_phase()` and `False` for `is_solid_phase()` and `is_vapor_phase()`.

The generic Phase object determines what to return for each method based on the user-provided name for the instance of the Phase object as shown below:

* `is_liquid_phase()` returns `True` if the Phase name contains the string `Liq`, otherwise it returns `False`.
* `is_solid_phase()` returns `True` if the Phase name contains the string `Sol`, otherwise it returns `False`.
* `is_vapor_phase()` returns `True` if the Phase name contains the string `Vap`, otherwise it returns `False`.

Users should avoid using the generic `Phase` object, as this is primarily intended as a base class for the specific phase classes and for backwards compatibility.

Phase Class
^^^^^^^^^^^

.. module:: idaes.core.phases

.. autoclass:: Phase
    :members:

Phase Type Enum
---------------

In some cases, it is useful to be able to indicate a given type of phase, rather than an instance specific `Phase` class; an example would be indicating the set of valid phases for a given chemical species. In these cases, the `PhaseType` `Enum` can be used, which enumerates the different types of phases recognized by the IDAES framework.

The `PhaseType` `Enum` has the following possible values:

* `liquidPhase` (1)
* `vaporPhase` (2)
* `solidPhase` (3)
