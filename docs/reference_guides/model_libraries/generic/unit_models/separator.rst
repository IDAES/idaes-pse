Separator
=========

The IDAES Separator unit model represents operations where a single stream is split into multiple flows. The Separator model supports separation using split fractions, or by ideal separation of flows. The Separator class can be used to create either a stand-alone separator unit, or as part of a unit model where a flow needs to be separated.

Degrees of Freedom
------------------

Separator units have a number of degrees of freedom based on the separation type chosen.

* If `split_basis` = 'phaseFlow', degrees of freedom are generally :math:`(no. outlets-1) \times no. phases`
* If `split_basis` = 'componentFlow', degrees of freedom are generally :math:`(no. outlets-1) \times no. components`
* If `split_basis` = 'phaseComponentFlow', degrees of freedom are generally :math:`(no. outlets-1) \times no. phases \times no. components`
* If `split_basis` = 'totalFlow', degrees of freedom are generally :math:`(no. outlets-1) \times no. phases \times no. components`

Typical fixed variables are:

* split fractions.

Model Structure
---------------

The IDAES Separator unit model does not use ControlVolumes, and instead writes a set of material, energy and momentum balances to split the inlet stream into a number of outlet streams. Separator models have a single inlet Port (named inlet) and a user-defined number of outlet Ports (by default named outlet_1, outlet_2, etc.).

**Mixed State Block**

If a mixed state block is provided in the construction arguments, the Mixer model will use this as the StateBlock for the mixed stream in the resulting balance equations. This allows a Mixer unit to be used as part of a larger unit operation by linking to an existing StateBlock.

Ideal Separation
----------------

The IDAES Separator model supports ideal separations, where all of a given subset of the mixed stream is sent to a single outlet (i.e. split fractions are equal to zero or one). In these cases, no Constraints are necessary for performing the separation, as the mixed stream states can be directly partitioned to the outlets.

Ideal separations will not work for all choices of state variables, and thus will not work for all property packages. To use ideal separations, the user must provide a map of what part of the mixed flow should be partitioned to each outlet. The `ideal_split_map` should be a dict-like object with keys as tuples matching the `split_basis` argument and values indicating which outlet this subset should be partitioned to.

Variables
---------

Separator units have the following variables (:math:`o` indicates index by outlet):

=============== ====================== ===========================================
Variable Name   Symbol                 Notes
=============== ====================== ===========================================
split_fraction  :math:`\phi_{t, o, *}` Indexing sets depend upon `split_basis`
=============== ====================== ===========================================

Constraints
-----------

Separator units have the following Constraints, unless `ideal_separation` is True.

* If `material_balance_type` is `componentPhase`:

`material_splitting_eqn(t, o, p, j)`:

.. math:: F_{in, t, p, j} = \phi_{t, p, *} \times F_{t, o, p, j}

* If `material_balance_type` is `componentTotal`:

`material_splitting_eqn(t, o, j)`:

.. math:: \sum_p{F_{in, t, p, j}} = \sum_p{\phi_{t, p, *} \times F_{t, o, p, j}}

* If `material_balance_type` is `total`:

`material_splitting_eqn(t, o)`:

.. math:: \sum_p{\sum_j{F_{in, t, p, j}}} = \sum_p{\sum_j{\phi_{t, p, *} \times F_{t, o, p, j}}}

If `energy_split_basis` is `equal_temperature`:

`temperature_equality_eqn(t, o)`:

.. math:: T_{in, t} = T_{t, o}

If `energy_split_basis` is `equal_molar_enthalpy`:

`molar_enthalpy_equality_eqn(t, o)`:

.. math:: h_{in, t} = h_{t, o}

If `energy_split_basis` is `enthalpy_split`:

`molar_enthalpy_splitting_eqn(t, o)`:

.. math:: sum_p{h_{in, t, p}*sf_{t, o, p}} = sum_p{h_{t, o, p}}

If `momentum_balance_type` is `pressureTotal`:

`pressure_equality_eqn(t, o)`:

.. math:: P_{in, t} = P_{t, o}

Separators do not support momentum balances using the `pressurePhase`, `momentumTotal` or `momentumPhase` options.

Separator Class
---------------

.. module:: idaes.models.unit_models.separator

.. autoclass:: Separator
  :members:

SeparatorData Class
-------------------

.. autoclass:: SeparatorData
  :members:
