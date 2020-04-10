Component Objects
=================

.. contents:: Contents 
    :depth: 2

Introduction
------------

Component objects are used in the IDAES Process Modeling Framework to identify the chemical species of interest in a property package and to contain information describing the behavior of that component (such as properties of that component).

The IDAES Process Modeling Framework currently supports the following types of components.

* `Component` - general purpose object for representing chemical species.
* `Solute` - component object for representing species which should be treated as a solute in a `LiquidPhase`.
* `Solvent` - component object for representing species which should be treated as a solvent in a `LiquidPhase`.
* `Ion` - general purpose component object for representing ion species (`LiquidPhase` only). Users should generally use the `Anion` or `Cation` components instead.
* `Anion` - component object for representing ion species with a negative charge (`LiquidPhase` only).
* `Cation` - component object for representing ion species with a positive charger(`LiquidPhase` only).

Component Object Methods
------------------------

Component objects are intended to store all the necessary information regarding a given chemical species for use within a process model. Examples of such information include the methods and parameters required for calculating thermophysical properties. Additionally, certain unit operations handle components in different ways depending on certain criteria. An example of this is Reverse Osmosis membranes, where the driving force across the member is calculated differently for solvent species and solute species.

Component objects implement the following methods for determining species behavior:

* `is_solute()` - returns `True` if species is a solute (`Solute`, `Ion`, `Anion` or `Cation` component objects), otherwise `False`.
* `is_solvent()` - returns `True` if species is a solvent (`Solvent` component object), otherwise `False`.

.. note:: The general purpose `Component` object does not distinguish solutes and solvents, and these methods will will raise a `TypeError` instead.

Types of Components
-------------------

Component Class
^^^^^^^^^^^^^^^

This is a general purpose Component object, and is suitable for general cases where the user is not concerned about distinguishing solutes from solvents (`is_solute()` and `is_solvent()` will both raise `TypeErrors`). The also forms the base class for all other Component types.

.. module:: idaes.core.components

.. autoclass:: Component
    :members:

Solute Class
^^^^^^^^^^^^

The component object is suitable for species which should be treated as solutes in a `LiquidPhase`. The only difference between this and a general `Component` is that `is_solute()` returns `True` and `is_solvent()` returns `False`.

Solvent Class
^^^^^^^^^^^^^

The component object is suitable for species which should be treated as solvents in a `LiquidPhase`. The only difference between this and a general `Component` is that `is_solute()` returns `False` and `is_solvent()` returns `True`.

Ion Class
^^^^^^^^^

The `Ion` class is suitable for ionic species which appear in `LiquidPhases`. This is similar to the `Solute` class, in that `is_solute()` returns `True` and `is_solvent()` returns `False`. Additionally, `Ion` objects have a `charge` configuration argument for recording the charge on the ion (must be an integer) and do not have a `valid_phase_types` argument (as it is assumed they can only exist in `LiquidPhases`).

.. note:: Users are encouraged to use the `Anion` and `Cation` classes instead of the generic `Ion` class, as these validate that sign of the `charge` configuration argument.

Anion Class
^^^^^^^^^^^

The `Anion` class is suitable for anionic species (i.e. negatively charged) which appear in `LiquidPhases`. This is a subclass of `Ion`, which enforces that the sign on the `charge` configuration argument be negative.

Cation Class
^^^^^^^^^^^^

The `Cation` class is suitable for cationic species (i.e. positively charged) which appear in `LiquidPhases`. This is a subclass of `Ion`, which enforces that the sign on the `charge` configuration argument be positive.
