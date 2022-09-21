Component Class
===============

This is a general purpose Component object, and is suitable for general cases where the user is 
not concerned about distinguishing solutes from solvents (`is_solute()` and `is_solvent()` will 
both raise `TypeErrors`). This also forms the base class for all other Component types.

.. module:: idaes.core.base.components

.. autoclass:: Component
    :members:

Solute Class
^^^^^^^^^^^^

The component object is suitable for species which should be treated as solutes in a 
`LiquidPhase`. The only difference between this and a general `Component` is that `is_solute()` 
returns `True` and `is_solvent()` returns `False`.

Solvent Class
^^^^^^^^^^^^^

The component object is suitable for species which should be treated as solvents in a 
`LiquidPhase`. The only difference between this and a general `Component` is that `is_solute()` 
returns `False` and `is_solvent()` returns `True`.

Ion Class
^^^^^^^^^

The `Ion` class is suitable for ionic species which appear in `LiquidPhases`. This is similar to 
the `Solute` class, in that `is_solute()` returns `True` and `is_solvent()` returns `False`. 
Additionally, `Ion` objects have a `charge` configuration argument for recording the charge on 
the ion (must be an integer) and do not have a `valid_phase_types` argument (as it is assumed 
they can only exist in `LiquidPhases`).

.. note:: Users are encouraged to use the `Anion` and `Cation` classes instead of the generic `Ion` class, as these validate that sign of the `charge` configuration argument.

Anion Class
^^^^^^^^^^^

The `Anion` class is suitable for anionic species (i.e. negatively charged) which appear in 
`LiquidPhases`. This is a subclass of `Ion`, which enforces that the sign on the `charge` 
configuration argument be negative.

Cation Class
^^^^^^^^^^^^

The `Cation` class is suitable for cationic species (i.e. positively charged) which appear in 
`LiquidPhases`. This is a subclass of `Ion`, which enforces that the sign on the `charge` 
configuration argument be positive.
