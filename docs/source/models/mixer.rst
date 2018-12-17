Mixer
=====

The IDAES Mixer unit model represents operations where multiple streams of material are combined into a single flow. The Mixer class can be used to create either a stand-alone mixer unti, or as part of a unit model where mutliple streams need to be mixed.

**Degrees of Freedom**

Mixer units have zero degrees of freedom.

**Model Structure**

The IDAES Mixer unit model does not use ControlVolumes, and instead writes a set of material, energy and momentum balances to combine the inlet streams into a singl mixed stream. Mixer modles have a user-defined number of inlet Ports (by default names inlet_1, inlet_2, etc.) and one outlet Port (by default names outlet).

**Construction Arguments**

Mixer units have the following construction arguments:

* property_package - property package to use when constructing StateBlocks (default = 'use_parent_value'). This is provided as a Physical Parameter Block by the Flowsheet when creating the model. If a value is not provided, the uunit model will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the StateBlocks when they are created.

Mixers also have the following configuration arguments for determining the form of the balance equations. Onlt one of `inlet_list` and `num_inlets` needs to be provided, and if both are provided and do not agree (i.e. the lenght of `inlet_list` is not equal to `num_inlets`) a ConfigurationError will be returned.

=========================== =========================== ======================================================
Argument                    Default Value               Description
=========================== =========================== ======================================================
inlet_list                  None                        A list of names to use for inlets
num_inlets                  2                           Number of inlets to Mixer
calculate_phase_equilibrium False                       Whether mixed stream should be in phase equilibrium
material_mixing_type        MixingType.extensive
energy_mixing_type          MixingType.extensive
momentum_mixing_type        MomentumMixingType.minimize
mixed_state_block           None                        Optional - existing StateBlock to use for mixed stream
construct_ports             True                        If False, Port objects will not be constructed
=========================== =========================== ======================================================

**Mixed State Block**

If a mixed state block is provided in the construction arguments, the Mixer model will use this as the StateBlock for the mixed stream in the resulting balance equations. This allows a Mixer unit to be used as part of a larger unit operation by linking multiple inlet streams to a single exisiting StateBlock.

**Variables**

Mixer units have the following variables (:math:`i` indicates index by inlet):

============================ ===================== ===========================================================================
Variable Name                Symbol                Notes
============================ ===================== ===========================================================================
phase_equilibrium_generation :math:`X_{eq, t, r}`  Only if has_phase_equilibrium = True, Generation term for phase equilibrium
minimum_pressure             :math:`P_{min, t, i}` Only if momentum_mixing_type = MomemntumMixingType.minimize
============================ ===================== ===========================================================================

**Parameters**

Mixer units have the following parameters:

============= ================ =====================================================================================
Variable Name Symbol           Notes
============= ================ =====================================================================================
eps_pressure  :math:`\epsilon` Only if momentum_mixing_type = MomemntumMixingType.minimize, smooth minimum parameter
============= ================ =====================================================================================

**Constraints**

The constraints written by the Mixer model depend upon the construction arguments chosen.

If `material_mixing_type` is `extensive`:

`material_mixing_equations(t, p, j)`:

.. math:: 0 = \sum_i{F_{in, i, p, j}} - F_{out, p, j} + \sum_r {n_{r, p, j} \times X_{eq, t, r}}

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

If `momentum_mixing_type` is `equality`, the pressure in all inlets and the outlet are equated. **N.B.** This may result in an overspecified problem if the user is not careful.

`pressure_equality_constraints(t, i)`:

.. math:: P_{mix, t} = P_{t, i}

Mixer Class
-----------

.. module:: idaes.unit_models.mixer

.. autoclass:: Mixer
  :members:

.. autoclass:: MixerData
  :members:
