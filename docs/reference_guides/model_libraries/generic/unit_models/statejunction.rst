StateJunction Block
===================

The IDAES StateJunction block represents a pass-through unit or simple pipe with no holdup. The primary use for this unit is in conceptual design applications for linking Arcs to/from different process alternatives.

Degrees of Freedom
------------------

StateJunctions have no degrees of freedom.

Model Structure
---------------

A StateJunction consists of a single StateBlock with two Ports (inlet and outlet), where the state variables in the state block are simultaneously connected to both Ports.

Additional Constraints
----------------------

StateJunctions write no additional constraints beyond those in the StateBlock.

Variables
---------

StateJunctions have no additional variables.

.. module:: idaes.models.unit_models.statejunction

StateJunction Class
-------------------

.. autoclass:: StateJunction
  :members:

StateJunctionData Class
-----------------------

.. autoclass:: StateJunctionData
  :members:
