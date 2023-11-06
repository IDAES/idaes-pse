Generic Solid-Liquid Separator
==============================

The ``SLSeparator`` unit model is a general purpose model for representing the separation of mixed streams containing solids and liquids, such as filtration and thickening. The model assumes that the mixed stream is separated into a concentrated solid-liquid stream and a stream of clarified liquid.

Degrees of Freedom
------------------

The ``SLSeparator`` model has 1 degree of freedom, which is generally the liquid recovery fraction. Note that due to the assumption that all solids report to the concentrated outlet, the solid inlet and outlet Ports link to the same variables (and thus cannot be fixed independently).

Model Structure
---------------

The model consists of two parts. The solid phase is represented by a single StateBlock as all solids in the inlet are assumed to report to the concentrated outlet. The liquid phase is represented by a standard :ref:`Separator <reference_guides/model_libraries/generic/unit_models/separator:Separator>` block which splits the incoming liquid between the concentrated and clarified outlets.

The ``SLSeparator`` unit model has a total of 5 ports:

* ``solid_inlet`` representing the solids in the incoming stream,
* ``liquid_inlet`` representing the liquid in the incoming stream,
* ``solid_outlet`` representing the solids in the concentrated stream,
* ``retained_liquid_outlet`` representing the liquid in the concentrated stream, and,
* ``recovered_liquid_outlet`` representing the liquid in the clarified stream.

Additional Constraints
----------------------

The ``SLSeparator`` model adds no additional constraints beyond those written by the :ref:`Separator <reference_guides/model_libraries/generic/unit_models/separator:Separator>` and StateBlocks.

Variables
---------

The ``SLSeparator`` adds one Reference to a variable in the :ref:`Separator <reference_guides/model_libraries/generic/unit_models/separator:Separator>`.

============ ================ ================================================================================================
Variable     Name             Notes
============ ================ ================================================================================================
:math:`R_t`  liquid_recovery  Fraction of liquid which reports to the clarified stream, reference to ``split.split_fraction``
============ ================ ================================================================================================

.. module:: idaes.models.unit_models.solid_liquid.sl_separator

SLSeparator Class
-----------------

.. autoclass:: SLSeparator
  :members:

SLSeparatorData Class
---------------------

.. autoclass:: SLSeparatorData
  :members:
