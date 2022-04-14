Continuous Stirred Tank Reactor
===============================

The IDAES CSTR model represents a unit operation where a material stream undergoes some chemical reaction(s) in a well-mixed vessel.

Degrees of Freedom
------------------

CSTRs generally have one degree of freedom. Typically, the fixed variable is reactor volume.

Model Structure
---------------

The core CSTR unit model consists of a single ``ControlVolume0D`` (named ``control_volume``) with one Inlet Port (named ``inlet``) and one Outlet Port (named ``outlet``).

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
:math:`Q_t`      heat   Only if ``has_heat_transfer`` = ``True``, reference to ``control_volume.heat``
================ ====== =========================================================================================================================

.. module:: idaes.models.unit_models.cstr

CSTR Class
----------

.. autoclass:: CSTR
  :members:

CSTRData Class
--------------

.. autoclass:: CSTRData
  :members:
