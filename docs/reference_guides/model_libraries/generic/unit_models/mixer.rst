Mixer
=====

The IDAES Mixer unit model represents operations where multiple streams of material are combined into a single flow. The Mixer class can be used to create either a stand-alone mixer unit, or as part of a unit model where multiple streams need to be mixed.

Degrees of Freedom
------------------

Mixer units have zero degrees of freedom.

Model Structure
---------------

The IDAES Mixer unit model does not use ControlVolumes, and instead writes a set of material, energy and momentum balances to combine the inlet streams into a single mixed stream. Mixer models have a user-defined number of inlet Ports (by default named inlet_1, inlet_2, etc.) and one outlet Port (named outlet).

**Mixed State Block**

If a mixed state block is provided in the construction arguments, the Mixer model will use this as the StateBlock for the mixed stream in the resulting balance equations. This allows a Mixer unit to be used as part of a larger unit operation by linking multiple inlet streams to a single existing StateBlock.

Variables
---------

Mixer units have the following variables (:math:`i` indicates index by inlet):

============================ ===================== ===========================================================================
Variable Name                Symbol                Notes
============================ ===================== ===========================================================================
phase_equilibrium_generation :math:`X_{eq, t, r}`  Only if has_phase_equilibrium = True, Generation term for phase equilibrium
minimum_pressure             :math:`P_{min, t, i}` Only if momentum_mixing_type = MomemntumMixingType.minimize
============================ ===================== ===========================================================================

Parameters
----------

Mixer units have the following parameters:

============= ================ =====================================================================================
Variable Name Symbol           Notes
============= ================ =====================================================================================
eps_pressure  :math:`\epsilon` Only if momentum_mixing_type = MomemntumMixingType.minimize, smooth minimum parameter
============= ================ =====================================================================================

Constraints
-----------

The constraints written by the Mixer model depend upon the construction arguments chosen.

If `material_mixing_type` is `extensive`:

* If `material_balance_type` is `componentPhase`:

`material_mixing_equations(t, p, j)`:

.. math:: 0 = \sum_i{F_{in, i, p, j}} - F_{out, p, j} + \sum_r {n_{r, p, j} \times X_{eq, t, r}}

* If `material_balance_type` is `componentTotal`:

`material_mixing_equations(t, j)`:

.. math:: 0 = \sum_p{(\sum_i{F_{in, i, p, j}} - F_{out, p, j} + \sum_r {n_{r, p, j} \times X_{eq, t, r}})}

* If `material_balance_type` is `total`:


`material_mixing_equations(t)`:

.. math:: 0 = \sum_p{\sum_j{(\sum_i{F_{in, i, p, j}} - F_{out, p, j} + \sum_r {n_{r, p, j} \times X_{eq, t, r}})}}

where :math:`n_{r, p, j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` in reaction :math:`r`.

If 'energy_mixing_type` is `extensive`:

`enthalpy_mixing_equations(t)`:

.. math:: 0 = \sum_i{\sum_p{H_{in, i, p}}} - \sum_p{H_{out, p}}

If 'momentum_mixing_type` is `minimize`, a series of smooth minimum operations are performed:

`minimum_pressure_constraint(t, i)`:

   For the first inlet:

      .. math:: P_{min, t, i} = P_{t, i}

   Otherwise:

      .. math:: P_{min, t, i} = smin(P_{min, t, i-1}, P_{t, i}, eps)

Here, :math:`P_{t, i}` is the pressure in inlet :math:`i` at time :math:`t`, :math:`P_{min, t, i}` is the minimum pressure in all inlets up to inlet :math:`i`, and :math:`smin` is the smooth minimum operator (see IDAES Utility Function documentation).

The minimum pressure in all inlets is then:

`mixture_pressure(t)`:

.. math:: P_{mix, t} = P_{min, t, i=last}

If `momentum_mixing_type` is `equality`, the pressure in all inlets and the outlet are equated.

.. note:: This may result in an over-specified problem if the user is not careful.

`pressure_equality_constraints(t, i)`:

.. math:: P_{mix, t} = P_{t, i}

Often the minimum inlet pressure constraint is useful for sequential modular type initialization, but the equal pressure constants are required for pressure-driven flow models.  In these cases it may be convenient to use the minimum pressure constraint for some initialization steps, them deactivate it and use the equal pressure constraints.  The `momentum_mixing_type` is `minimum_and_equality` this will create the constraints for both with the minimum pressure constraint being active.

The `mixture_pressure(t)` and `pressure_equality_constraints(t, i)` can be directly activated and deactivated, but only one set of constraints should be active at a time. The ``use_minimum_inlet_pressure_constraint()`` and ``use_equal_pressure_constraint()`` methods are also provided to switch between constant sets. 

Initialization
--------------

.. module:: idaes.models.unit_models.mixer

.. autoclass:: MixerInitializer
   :members: initialization_routine

Mixer Class
-----------

.. autoclass:: Mixer
  :members:

MixerData Class
---------------

.. autoclass:: MixerData
  :members:
