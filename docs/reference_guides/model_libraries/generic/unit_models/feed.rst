Feed Block
==========

Feed Blocks are used to represent sources of material in Flowsheets. Feed blocks do not require calculation of the phase equilibrium of the feed stream (some property packages, like the :ref:`Generic Property Package <explanations/components/property_package/general/index:Generic Property Package Framework>` with a :ref:`FTPx <explanations/components/property_package/general/state/FTPx:State Definition>` state definition, always flash regardless), and the composition of the material in the outlet stream will be exactly as specified in the input. For applications where the users wishes the outlet stream to be in phase equilibrium, see the Feed_Flash unit model.

Degrees of Freedom
------------------

The degrees of freedom of Feed blocks depends on the property package being used and the number of state variables necessary to fully define the system. Users should refer to documentation on the property package they are using.

Model Structure
---------------

Feed Blocks consists of a single StateBlock (named properties), each with one Outlet Port (named outlet). Feed Blocks also contain References to the state variables defined within the StateBlock.

Additional Constraints
----------------------

Feed Blocks do not write any additional constraints to the model themselves. However, the StateBlock they contain may write constraints.

Variables
---------

Feed blocks add no additional Variables.

.. module:: idaes.models.unit_models.feed

Feed Class
----------

.. autoclass:: Feed
  :members:

FeedData Class
--------------

.. autoclass:: FeedData
  :members:

