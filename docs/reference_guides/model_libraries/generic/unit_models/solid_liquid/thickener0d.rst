Thickener (0D)
==============

The ``Thickener0D`` unit model is an extension of the :ref:`SLSeparator <reference_guides/model_libraries/generic/unit_models/solid_liquid/sl_separator:Generic Solid-Liquid Separator>` model which adds constraints to estimate the area and height of a vessel required ot achive the desired separation of solid and liquid based on experimental measurements of the settling velocity. This model is based on corrleations described in:

[1] Couldson & Richardson's Chemical Engineering, Volume 2 Particle Technology & Separation Processes (4th Ed.), Butterworth-Heinemann (2001)

Degrees of Freedom
------------------

The ``Thickener0D`` model has 5 degrees of freedom, which is are generally chosen from:

* the liquid recovery fraction or underflow liquid-to-solid ratio,
* the liquid-to-solid ratio and settling velocity at the pinch (critical) point,
* the cross-sectional area of the thickener,
* the depth of the clarification zone,
* the depth of the sedimetnation zone or the required settling time (experimental).

Model Structure
---------------

The ``Thickener0D`` model has the same strucutre as the :ref:`SLSeparator <reference_guides/model_libraries/generic/unit_models/solid_liquid/sl_separator:Model Structure>`.

Additional Constraints
----------------------

The ``Thickener0D`` model adds the following additional constraints beyond those written by the :ref:`SLSeparator <reference_guides/model_libraries/generic/unit_models/solid_liquid/sl_separator:Additional Constraints>`.

The cross-sectional area of the thickener is calculated using:

.. math:: A = \frac{S_t}{\rho_{liq, t}} \times \frac{(Y_{pinch, t}-Y_{under, t})}{u_{pinch, t}}

where :math:`A` is the cross-section area ofthe thickener, :math:`S_t` is the mass flowrate of solids entering the thickener, :math:`\rho_{liq}` is hte mass density of the liquid phase, :math:`Y_{pinch}` and :math:`Y_{under}` are the mass-based liquid-to-solid ratios atthe pinch point and underflow respectively and :math:`u_{pinch}` is the sedimentaion velocity of the suspension at the pinch point (Eqn 5.54, pg. 198 in [1]).

The liquid-solid ratio at the underflow can be calculated using:

.. math:: Y_{under, t} = \frac{L_{under, t}}{S_t}

where :math:`L_{under}` is hte liquid mass flowrate at the underflow.

The total height of the thickener is calculated using:

.. math:: H = \frac{S_t \tau_{t}}{A\rho_{sol, t}} \times (1+\frac{\rho_{sol}}{\rho_{liq}} \times Y_{avg}) + H_{clarified}

where :math:`H` is the total height of the thickener, :math:`H_{clarified}` is the height of the clarification zone, :math:`\tau` is the empirically measured settling time required to achieve the desired underflow conditions and :math:`Y_avg` is the average liquid-solid ratio in the thickener (assumed to be linear) (Eqn 5.55, pg. 198 in [1]).

Variables
---------

The ``Thickener0D`` adds the following variables in addition to those in the :ref:`SLSeparator <reference_guides/model_libraries/generic/unit_models/solid_liquid/sl_separator:Variables>`.

============ ================ ================================================================================================
Variable     Name             Notes
============ ================ ================================================================================================
:math:`R_t`  liquid_recovery  Fraction of liquid which reports to the clarified stream, reference to ``split.split_fraction``
============ ================ ================================================================================================

.. module:: idaes.models.unit_models.solid_liquid.thickener

SLSeparator Class
-----------------

.. autoclass:: Thickener0D
  :members:

SLSeparatorData Class
---------------------

.. autoclass:: Thickener0DData
  :members:
