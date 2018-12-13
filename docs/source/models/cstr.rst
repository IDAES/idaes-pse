Continuous Stirred Tank Reactor
===============================

The IDAES CSTR model represents a unit operation where a material stream undergoes some chemical reaction(s) in a well-mixed vessel.

Degrees of Freedom
------------------

CSTRs generally have one degree of freedom. Typically, the fixed variable is reactor volume.

Model Structure
---------------

The core CSTR unit model consists of a single ``ControlVolumn0D`` (named ``control_volume``) with one Inlet Port (named ``inlet``) and one Outlet Port (named ``outlet``).

Construction Arguments
----------------------

CSTR units have the following construction arguments:

* ``property_package`` - property package to use when constructing State Blocks (default = ``useDefault``). This is provided as a ``PhysicalParameterBlock`` by the Flowsheet or the parent model when creating the model. If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* ``property_package_args`` - set of arguments to be passed to the State Blocks when they are created.
* ``reaction_package`` - reaction package to use when constructing Reaction Blocks (default = ``None``).
* ``reaction_package_args`` - set of arguments to be passed to the Reaction Blocks when they are created.

Additionally, CSTR units have the following construction arguments which are passed to the ControlVolume Block for determining which terms to construct in the balance equations.

============================= ======================================
Argument                      Default Value
============================= ======================================
``material_balance_type``     ``MaterialBalanceType.componentPhase``
``energy_balance_type``       ``EnergyBalanceType.enthalpyTotal``
``momentum_balance_type``     ``MomentumBalanceType.pressureTotal``
``dynamic``                   ``False``
``has_holdup``                ``False``
``has_rate_reactions``        ``True``
``has_equilibrium_reactions`` ``True``
``has_heat_of_reaction``      ``False``
``has_mass_transfer``         ``False``
``has_heat_transfer``         ``False``
``has_pressure_change``       ``False``
============================= ======================================

Additional Constraints
----------------------

CSTR units write the following additional Constraints beyond those written by the ControlVolume Block.

.. math:: X_{t,r} = V_t \times r_{t,r}

where :math:`X_{t,r}` is the extent of reaction of reaction :math:`r` at time :math:`t`, :math:`V_t` is the volume of the reacting material at time :math:`t` (allows for varying reactor volume with time) and :math:`r_{t,r}` is the volumetric rate of reaction of reaction :math:`r` at time :math:`t` (from the outlet property package).

Variables
---------

CSTR units add the following additional Variables beyond those created by the ControlVolume Block.

================ ====== =========================================================================================================================
Variable         Name   Notes
================ ====== =========================================================================================================================
:math:`V_t`      volume If ``has_holdup`` = ``True`` this is a reference to ``control_volume.volume``, otherwise a Var attached to the Unit Model
:math:`Q_t`      heat   Only ``if has_heat_transfer`` = ``True``, reference to ``control_volume.heat``
================ ====== =========================================================================================================================

CSTRData Class
--------------

.. module:: idaes.unit_models.cstr

.. autoclass:: CSTR
  :members:

.. autoclass:: CSTRData
  :members:
