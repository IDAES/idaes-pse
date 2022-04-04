Equilibrium Reactor
===================

The IDAES Equilibrium reactor model represents a unit operation where a material stream undergoes some chemical reaction(s) to reach an equilibrium state. This model is for systems with reaction with equilibrium coefficients - for Gibbs energy minimization see Gibbs reactor documentation.

Degrees of Freedom
------------------

Equilibrium reactors generally have 1 degree of freedom.

Typical fixed variables are:

* reactor heat duty (has_heat_transfer = True only).

Model Structure
---------------

The core Equilibrium reactor unit model consists of a single ``ControlVolume0D`` (named ``control_volume``) with one Inlet Port (named ``inlet``) and one Outlet Port (named ``outlet``).

Additional Constraints
----------------------

Equilibrium reactors units write the following additional Constraints beyond those written by the Control Volume if rate controlled reactions are present.

.. math:: r_{t,r} = 0

where :math:`r_{t,r}` is the rate of reaction for reaction :math:`r` at time :math:`t`. This enforces equilibrium in any reversible rate controlled reactions which are present. Any non-reversible reaction that may be present will proceed to completion.

Variables
---------

Equilibrium reactor units add the following additional Variables beyond those created by the Control Volume.

================ ====== =========================================================================================================================
Variable         Name   Notes
================ ====== =========================================================================================================================
:math:`V_t`      volume If ``has_holdup`` = ``True`` this is a reference to ``control_volume.volume``, otherwise a Var attached to the Unit Model
:math:`Q_t`      heat   Only if ``has_heat_transfer`` = ``True``, reference to ``control_volume.heat``
================ ====== =========================================================================================================================

.. module:: idaes.models.unit_models.equilibrium_reactor

EquilibriumReactor Class
------------------------

.. autoclass:: EquilibriumReactor
  :members:

EquilibriumReactorData Class
----------------------------

.. autoclass:: EquilibriumReactorData
  :members:
