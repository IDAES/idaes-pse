Stoichiometric (Yield) Reactor
==============================

The IDAES Stoichiometric reactor model represents a unit operation where a single material stream undergoes some chemical reaction(s) subject to a set of extent or yield specifications.

Degrees of Freedom
------------------

Stoichiometric reactors generally have degrees of freedom equal to the number of reactions + 1.

Typical fixed variables are:

* reaction extents or yields (1 per reaction),
* reactor heat duty (has_heat_transfer = True only).

Model Structure
---------------

The core Stoichiometric reactor unit model consists of a single ControlVolume0D (named control_volume) with one Inlet Port (named inlet) and one Outlet Port (named outlet).

Construction Arguments
----------------------

Stoichiometric reactor units have the following construction arguments:

* dynamic - indicates whether this model will be dynamic or not. The default = False which implies the model is steady-state.
* material_balance_type - indicates what type of mass balance should be constructed. The default is componentPhase, which uses phase component balances.
* energy_balance_type - indicates what type of energy balance should be constructed. The default is enthalpyTotal, which uses a single energy balance for the inlet stream.
* momentum_balance_type - indicates what type of momentum balance should be constructed. The default is None, which excludes momentum balance.
* has_heat_of_reaction - indicates whether there is heat of reaction. The default is False.
* has_heat_transfer - indicates whether there is heat transfer. The default is False.
* has_pressure_change - indicates whether there is pressure change. The default is False.
* property_package - property package to use when constructing Property Blocks (default = None). This is provided as a Property Parameter Block by the Flowsheet when creating the model. If a value is not provided, the control_volume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the Property Blocks (default = None) when they are created.
* reaction_package - reaction package to use for control volume (default = None). This is provided as a Reaction Parameter Block.
* reaction_package_args - set of arguments to be passed to the Reactor Blocks (default = None) when they are created.

Additionally, Stoichiometric reactor units have the following construction arguments which are passed to the control_volume Block for determining which terms to construct in the balance equations.

========================= =================
Argument                  Default Value
========================= =================
has_rate_reactions        True
========================= =================

Additional Constraints
----------------------

Stoichiometric reactor units write no additional Constraints beyond those written by the control_volume Block.

Variables
---------

Stoichiometric reactors units add the following additional variables beyond those created by the control_volume Block.

================ ==================== ===========================================================================
Variable         Name                 Notes
================ ==================== ===========================================================================
:math:`Q_t`      heat                 Only if has_heat_transfer = True, reference to control_volume.heat
:math:`deltaP`    pressure change     Only if has_pressure_change = True, reference to control_volume.deltaP
================ ==================== ===========================================================================

StoichiometricReactorData Class
-------------------------------

..module:: idaes.unit_models.stoichiometric_reactor

..autoclass:: StoichiometricReactorData
    :members:
