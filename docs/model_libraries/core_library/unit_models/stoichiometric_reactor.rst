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

The core Stoichiometric reactor unit model consists of a single ControlVolume0DBlock (named control_volume) with one Inlet Port (named inlet) and one Outlet Port (named outlet).

Variables
---------

Stoichiometric reactors units add the following variables:

================ ==================== ===========================================================================
Variable         Name                 Notes
================ ==================== ===========================================================================
:math:`Q_t`      heat                 Only if has_heat_transfer = True, reference to control_volume.heat
:math:`deltaP_t` pressure change      Only if has_pressure_change = True, reference to control_volume.deltaP
================ ==================== ===========================================================================

Constraints
-----------

Stoichiometric reactor units write no additional Constraints beyond those written by the control_volume Block.

StoichiometricReactor Class
---------------------------

.. module:: idaes.generic_models.unit_models.stoichiometric_reactor

.. autoclass:: StoichiometricReactor
  :members:

StoichiometricReactorData Class
-------------------------------

.. autoclass:: StoichiometricReactorData
    :members:
