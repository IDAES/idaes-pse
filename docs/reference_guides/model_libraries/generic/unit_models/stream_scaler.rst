Stream Scaler Block
===================

Stream Scaler Blocks are used to adjust size of streams to represent, for example, a stream being split across several identical units, which are then all modeled as a single IDAES unit

Degrees of Freedom
------------------

Stream Scaler blocks have one degree of freedom (and an additional variable).

Model Structure
---------------

Stream Scaler Blocks consists of a single StateBlock (named properties), each with an inlet and outlet port.

Additional Constraints
----------------------

Stream Scaler Blocks write no additional constraints* (besides those naturally occurring in StateBlocks).

Variables
---------

Stream Scaler blocks add no additional Variables.

.. module:: idaes.models.unit_models.stream_scaler


Initialization
--------------

.. autoclass:: StreamScalerInitializer
   :members: initialization_routine

StreamScaler Class
------------------

.. autoclass:: StreamScaler
  :members:

StreamScalerData Class
----------------------

.. autoclass:: StreamScalerData
  :members:

