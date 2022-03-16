Product Block
=============

Product Blocks are used to represent sinks of material in Flowsheets. These can be used as a conventient way to mark the final destination of a material stream and to view the state of that material.

Degrees of Freedom
------------------

Product blocks generally have zero degrees of freedom.

Model Structure
---------------

Product Blocks consists of a single StateBlock (named properties), each with one Inlet Port (named inlet). Product Blocks also contain References to the state variables defined within the StateBlock

Additional Constraints
----------------------

Product Blocks write no additional constraints to the model.

Variables
---------

Product blocks add no additional Variables.

.. module:: idaes.generic_models.unit_models.product

Product Class
-------------

.. autoclass:: Product
  :members:

ProductData Class
-----------------

.. autoclass:: ProductData
  :members:

