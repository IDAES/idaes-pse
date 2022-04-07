Stoichiometric (Yield) Reactor
==============================

The IDAES Stoichiometric reactor model represents a unit operation where a single material stream undergoes some chemical reaction(s) subject to a set of extent or yield specifications.

Degrees of Freedom
------------------

Stoichiometric reactors generally have degrees of freedom equal to the number of reactions + 1.

Typical fixed variables are:

* reaction extents or yields (1 per reaction),
* reactor heat duty (has_heat_transfer = True only).

Alternative Specifications
==========================

Alternative approaches for specifying the performance of stoichiometric reactors are parameters such as conversion or selectivity. However, these parameters are defined with respect to a specific key component within the system, which is difficult to implement in a general fashion. Users are encouraged to implement additional variables and constraints to represent these parameters if desired.

Model Structure
---------------

The core Stoichiometric reactor unit model consists of a single ControlVolume0DBlock (named control_volume) with one Inlet Port (named inlet) and one Outlet Port (named outlet).

Variables
---------

Stoichiometric reactors units add the following variables:

================ ===================== ===========================================================================
Variable         Name                  Notes
================ ===================== ===========================================================================
:math:`X_t`      rate_reaction_extent  Extent of reaction, indexed by time and reactions
:math:`Q_t`      heat                  Only if has_heat_transfer = True, reference to control_volume.heat
:math:`deltaP_t` pressure change       Only if has_pressure_change = True, reference to control_volume.deltaP
================ ===================== ===========================================================================

Constraints
-----------

Stoichiometric reactor units write no additional Constraints beyond those written by the control_volume Block.

StoichiometricReactor Class
---------------------------

.. module:: idaes.models.unit_models.stoichiometric_reactor

.. autoclass:: StoichiometricReactor
  :members:

StoichiometricReactorData Class
-------------------------------

.. autoclass:: StoichiometricReactorData
    :members:
