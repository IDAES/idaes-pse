Feed Block
==========

Feed Blocks are used to represent sources of material in Flowsheets. These can be used to determine the full state of a material (including equilibrium) based on a sufficient set of state variables prior to being passed to the first unit operation.

Degrees of Freedom
------------------

The degrees of freedom of Feed blocks depends on the property package being used and the number of state variables necessary to fully define the system. Users should refer to documentation on the property package they are using.

Model Structure
---------------

Feed Blocks consists of a single StateBlock (named properties), each with one Outlet Port (named outlet). Feed Blocks also contain References to the state variables defined within the StateBlock

Additional Constraints
----------------------

Feed Blocks write no additional constraints to the model.

Variables
---------

Feed blocks add no additional Variables.

.. module:: idaes.unit_models.feed

Feed Class
----------

.. autoclass:: Feed
  :members:

FeedData Class
--------------

.. autoclass:: FeedData
  :members:

