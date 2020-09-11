Component Object
================

Component objects are used to identify the chemical 
species of interest in a property package and to contain information describing the behavior 
of that component (such as properties of that component). Additional information on the 
:ref:`Component Class<technical_specs/core/comp:Component Class>` is provided in the technical
specifications.

The following types of components are currently supported.

* `Component` - general purpose object for representing chemical species.
* `Solute` - component object for representing species which should be treated as a solute in a `LiquidPhase`.
* `Solvent` - component object for representing species which should be treated as a solvent in a `LiquidPhase`.
* `Ion` - general purpose component object for representing ion species (`LiquidPhase` only). Users should generally use the `Anion` or `Cation` components instead.
* `Anion` - component object for representing ion species with a negative charge (`LiquidPhase` only).
* `Cation` - component object for representing ion species with a positive charger(`LiquidPhase` only).

Component objects are intended to store all the necessary information regarding a given 
chemical species for use within a process model. Examples of such information include the 
methods and parameters required for calculating thermophysical properties. Additionally, 
certain unit operations handle components in different ways depending on certain criteria. 
An example of this is Reverse Osmosis, where the driving force across the membrane is calculated 
differently for solvent species and solute species.

Component objects implement the following methods for determining species behavior:

* `is_solute()` - returns `True` if species is a solute (`Solute`, `Ion`, `Anion` or `Cation` component objects), otherwise `False`.
* `is_solvent()` - returns `True` if species is a solvent (`Solvent` component object), otherwise `False`.

.. note:: The general purpose `Component` object does not distinguish solutes and solvents, and these methods will will raise a `TypeError` instead.
