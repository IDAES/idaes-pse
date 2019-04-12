Node Block
==========

The IDAES Node block represents a pass-through unit or simple pipe with no holdup. The primary use for this unit is in conceptual design applications for linking Arcs to/from different process alternatives.

Degrees of Freedom
------------------

Nodes have no degrees of freedom.

Model Structure
---------------

A Node consists of a single StateBlock with two Ports (inlet and outlet), where the state variables in the state block are simultaneously connected to both Ports.

Additional Constraints
----------------------

Nodes write no additional constraints beyond those in the StateBlock.

Variables
---------

Nodes have no additional variables.

.. module:: idaes.unit_models.node

Node Class
----------

.. autoclass:: Node
  :members:

NodeData Class
--------------

.. autoclass:: NodeData
  :members:
