Translator Block
================

Translator blocks are used in complex flowsheets where the user desires to use different property packages for different parts of the flowsheet. In order to link two streams using different property packages, a translator block is required.

The core translator block provides a general framework for constructing Translator Blocks, however users need to add constraints to map the incoming states to the outgoing states as required by their specific application.

Degrees of Freedom
------------------

The degrees of freedom of Translator blocks depends on the property packages being used, and the user should write a sufficient number of constraints mapping inlet states to outlet states to satisfy these degrees of freedom.

Model Structure
---------------

The core Translator Block consists of two State Blocks, names ``properties_in`` and ``properties_out``, which are linked to two Ports names ``inlet`` and ``outlet`` respectively.

Additional Constraints
----------------------

The core Translator Block writes no additional constraints. Users should add constraints to their instances as required.

Variables
---------

Translator blocks add no additional Variables.

.. module:: idaes.models.unit_models.translator

Translator Class
----------------

.. autoclass:: Translator
  :members:

TranslatorData Class
--------------------

.. autoclass:: TranslatorData
  :members:


