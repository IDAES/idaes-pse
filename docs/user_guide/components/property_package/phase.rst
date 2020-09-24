Phase Object
============

Phase objects are used to identify the thermodynamic 
phases of interest in a property package and to contain information describing the behavior of 
that phase (for example the equation of state which describes that phase). Additional 
information on the :ref:`Phase Class<technical_specs/core/phase:Phase Class>` is 
provided in the technical specifications.

TThe following types of phases, along with a generic Phase object, are supported:

* `LiquidPhase`
* `SolidPhase`
* `VaporPhase`

In a number of unit operations, different phases behave in different ways. For example, in a 
Flash operation, the vapor phase exits through the top outlet whilst liquid phase(s) 
(and any solids) exit through the bottom outlet. In order to determine how a given phase should 
behave in these situations, each Phase object implements the following three methods:

* `is_liquid_phase()`
* `is_solid_phase()`
* `is_vapor_phase()`

These methods return a boolean (`True` or `False`) indicating whether the unit operation should 
treat the phase as being of the specified type in order to decide on how it should behave. Each 
type of phase returns `True` for its type and `False` for all other types (e.g. LiquidPhase 
returns `True` for `is_liquid_phase()` and `False` for `is_solid_phase()` and `is_vapor_phase()`.

The generic Phase object determines what to return for each method based on the user-provided 
name for the instance of the Phase object as shown below:

* `is_liquid_phase()` returns `True` if the Phase name contains the string `Liq`, otherwise it returns `False`.
* `is_solid_phase()` returns `True` if the Phase name contains the string `Sol`, otherwise it returns `False`.
* `is_vapor_phase()` returns `True` if the Phase name contains the string `Vap`, otherwise it returns `False`.

Users should avoid using the generic `Phase` object, as this is primarily intended as a base 
class for the specific phase classes and for backwards compatibility.
