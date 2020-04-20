Phase Objects
=============

.. contents:: Contents 
    :depth: 2

Introduction
------------

Phase objects are used in the IDAES Process Modeling Framework to identify the thermodynamic phases of interst in a property package and to contain information describing the behavior of that phase (for example the equation of state which describes that phase).

The IDAES Process Mdoeling Framework currently supports the follwoing types of phases, along with a generic Phase object.

* LiquidPhase
* SolidPhase
* VaporPhase

Phase Type
----------

In a number of unit operations, differnt phases behave in different ways. For example, in a Flash operation, the vapor phase exits throguh the top outlet whilst liquid phase(s) (and any solids) exit throguh the bottom outlet. In order to determine how a given phase should behave in these situations, each Phase object implements the following three methods:

* `is_liquid_phase()`
* `is_solid_phase()`
* `is_vapor_phase()`

These methods return a boolean (`True` or `False`) indicating whether the unit operation should treat the phase as being of the specified type in order to decide on how it should behave. Each type of phase returns `True` for its type and `False` for all other types (e.g. LiquidPhase returns `True` for `is_liquid_phase()` and `False` for `is_solid_phase()` and `is_vapor_phase()`.

The generic Phase object determines what to return for each method based on the user-provided name for the instance of the Phase object as shown below:

* `is_liquid_phase()` returns `True` if the Phase name contains the string `Liq`, otherwise it returns `False`.
* `is_solid_phase()` returns `True` if the Phase name contains the string `Sol`, otherwise it returns `False`.
* `is_vapor_phase()` returns `True` if the Phase name contains the string `Vap`, otherwise it returns `False`.

Users should avoid using the generic Phase object, as this is primarily intended as a base class for the specific phase classes and for backwards compatability.

Phase Class
^^^^^^^^^^^

.. module:: idaes.core.phases

.. autoclass:: Phase
    :members:


