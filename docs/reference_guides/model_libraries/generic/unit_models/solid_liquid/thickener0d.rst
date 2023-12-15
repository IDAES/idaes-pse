Thickener (0D)
==============

.. warning::
    The Thickener model is currently in beta status and will likely change in the next release as a more predictive version is developed.

The ``Thickener0D`` unit model is an extension of the :ref:`SLSeparator <reference_guides/model_libraries/generic/unit_models/solid_liquid/sl_separator:Generic Solid-Liquid Separator>` model which adds constraints to estimate the area and height of a vessel required to achieve the desired separation of solid and liquid based on experimental measurements of the settling velocity. This model is based on correlations described in:

[1] Coulson & Richardson's Chemical Engineering, Volume 2 Particle Technology & Separation Processes (4th Ed.), Butterworth-Heinemann (2001)

Sizing Thickeners and Pinch Point
---------------------------------

The approach for sizing the thickener vessel used in this model relies on identifying a pinch point in the thickener which is the limiting condition for the settling velocity of the suspension as described in [1]. The pinch point is described by the condition:

.. math:: max \left( \frac{(Y - Y_{under})}{u(Y)} \right)

where :math:`Y` is the mass based liquid-to-solid ratio, :math:`u(Y)` is the settling velocity of the suspension as a function of :math:`Y` and :math:`Y_{under}` is the liquid-solid ratio at the thickener underflow. Correlations exist which can predict :math:`u(Y)` however in many cases it is necessary to measure this experimentally. Additionally, whilst there are techniques for embedding the maximization operation within an equation-oriented model, these approaches tend to be highly non-linear and require the user to provide good values for scaling parameters and thus have not been implemented yet.

In the current implementation of the model, :math:`Y_{pinch}` and :math:`u_{pinch}` are created as user-defined input variables which should generally be fixed at appropriate values. Users may choose to add a correlation for :math:`u(Y)` if they choose.

Degrees of Freedom
------------------

The ``Thickener0D`` model has 5 degrees of freedom, which are generally chosen from:

* the liquid recovery fraction or underflow liquid-to-solid ratio,
* the liquid-to-solid ratio and settling velocity at the pinch (critical) point,
* the cross-sectional area of the thickener,
* the depth of the clarification zone,
* the depth of the sedimentation zone or the required settling time (experimental).

Model Structure
---------------

The ``Thickener0D`` model has the same structure as the :ref:`SLSeparator <reference_guides/model_libraries/generic/unit_models/solid_liquid/sl_separator:Model Structure>`.

Additional Constraints
----------------------

The ``Thickener0D`` model adds the following additional constraints beyond those written by the :ref:`SLSeparator <reference_guides/model_libraries/generic/unit_models/solid_liquid/sl_separator:Additional Constraints>`.

The cross-sectional area of the thickener is calculated using:

.. math:: A = \frac{S_t}{\rho_{liq, t}} \times \frac{(Y_{pinch, t}-Y_{under, t})}{u_{pinch, t}}

where :math:`A` is the cross-section area of the thickener, :math:`S_t` is the mass flowrate of solids entering the thickener, :math:`\rho_{liq}` is the mass density of the liquid phase, :math:`Y_{pinch}` and :math:`Y_{under}` are the mass-based liquid-to-solid ratios at the pinch point and underflow respectively and :math:`u_{pinch}` is the sedimentation velocity of the suspension at the pinch point (Eqn 5.54, pg. 198 in [1]).

The liquid-solid ratio at the underflow can be calculated using:

.. math:: Y_{under, t} = \frac{L_{under, t}}{S_t}

where :math:`L_{under}` is the liquid mass flowrate at the underflow.

The total height of the thickener is calculated using:

.. math:: H = \frac{S_t \tau_{t}}{A\rho_{sol, t}} \times \left( 1+\frac{\rho_{sol}}{\rho_{liq}} \times Y_{avg} \right) + H_{clarified}

where :math:`H` is the total height of the thickener, :math:`H_{clarified}` is the height of the clarification zone, :math:`\tau` is the empirically measured settling time required to achieve the desired underflow conditions and :math:`Y_avg` is the average liquid-solid ratio in the thickener (assumed to be linear) (Eqn 5.55, pg. 198 in [1]).

Variables
---------

The ``Thickener0D`` adds the following variables in addition to those in the :ref:`SLSeparator <reference_guides/model_libraries/generic/unit_models/solid_liquid/sl_separator:Variables>`.

===================== ======================== ====== =====================================================================
Variable              Name                     Index  Notes
===================== ======================== ====== =====================================================================
:math:`A`             area                     None   Cross-sectional area of thickener
:math:`H`             height                   None   Total height of thickener
:math:`H_{clarified}` height_clarified         None   Height of clarification zone in thickener (height above feed point)
:math:`u_{pinch}`     settling_velocity_pinch  time   Settling velocity of suspension at pinch point
:math:`Y_{pinch}`     liquid_solid_pinch       time   Liquid-solid ratio at pinch point
:math:`Y_{under}`     liquid_solid_underflow   time   Liquid-solid ratio at underflow
:math:`\tau`          settling_time            time   Settling time in thickener
===================== ======================== ====== =====================================================================

.. module:: idaes.models.unit_models.solid_liquid.thickener

SLSeparator Class
-----------------

.. autoclass:: Thickener0D
  :members:

SLSeparatorData Class
---------------------

.. autoclass:: Thickener0DData
  :members:
