Feed Block with Flash
=====================

Feed Blocks are used to represent sources of material in Flowsheets. In some cases, users may have a situation where a feed stream may be in a multi-phase state, but may not know the full details of the equilibrium state. The IDAES Feed Block with Flash (FeedFlash) allows users to define a feed block where the outlet is in phase equilibrium based on calculations from the chosen property package and a sufficient set of state variables prior to being passed to the first unit operation. The phase equilibrium is performed assuming an isobaric and isothermal flash operation.

A Feed Block with Flash is only required in cases where the feed may be in phase equilibrium AND the chosen property package uses a state definition that includes phase separations. Some property packages support phase equilibrium, but use a state definition that involves only total flows - in these cases a flash calculation is performed at the inlet of every unit and thus it is not necessary to perform a flash calculation at the feed block.

Degrees of Freedom
------------------

The degrees of freedom of FeedFlash blocks depends on the property package being used and the number of state variables necessary to fully define the system. Users should refer to documentation on the property package they are using.

Model Structure
---------------

FeedFlash Blocks contain a single ``ControlVolume0D`` (named ``control_volume``) with one Outlet Port (named ``outlet``). FeedFlash Blocks also contain References to the state variables defined within the inlet StateBlock of the ControlVolume (representing the unflashed state of the feed).

FeedFlash Blocks do not write a set of energy balances within the Control Volume - instead a constraint is written which enforces an isothermal flash.

Additional Constraints
----------------------

The FeedFlash Block writes one additional constraint to enforce isothermal behavior.

.. math:: T_{in, t} = T_{out, t}

where :math:`T_{in, t}` and :math:`T_{out, t}` are the temperatures of the material before and after the flash operation.

Variables
---------

FeedFlash blocks add no additional Variables.

.. module:: idaes.models.unit_models.feed_flash

FeedFlash Class
---------------

.. autoclass:: FeedFlash
  :members:

FeedFlashData Class
-------------------

.. autoclass:: FeedFlashData
  :members:

