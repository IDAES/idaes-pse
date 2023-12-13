Thickener (0D)
==============

The ``Thickener0D`` unit model is a predictive model for clarifiers and thickeners based on the following references:

[1] R. Burger, F. Concha, K.H. Karlsen, A. Narvaez, Numerical simulation of clarifier-thickener units treating ideal suspensions with a flux density function having two inflection points, Mathematical and Computer Modelling 44 (2006) 255–275, doi:10.1016/j.mcm.2005.11.008

[2] N.G. Barton, C.-H. Li, S.J. Spencer, Control of a surface of discontinuity in continuous thickeners, J. Aust. Math. Soc. Ser. B 33 (1992) 269–289


Degrees of Freedom
------------------

The ``Thickener0D`` model has 6 degrees of freedom, which are generally chosen from:

* the cross-sectional area of the settler ``area``,
* the solid particle size ``particle_size``,
* the maximum achievable solid fraction ``solid_fraction_max``,
* two empirical coefficients for calcuating the settling velocity ``v1``, ``C``,
* the volumetric flowrate or volume fraction of solids at the underflow.

Model Structure
---------------

The ``Thickener0D`` model contains two separators, ``solid_split`` and ``liquid_spit`` to separate the solid and liquid flows. The separations are assumed to be based on total flow with no change in composition, temperature or pressure. The ``Thickener0D`` model has 6 Ports (2 inlet and 4 outlet):

* solid_inlet,
* liquid_inlet,
* solid_overflow,
* liquid_overflow,
* solid_underflow,
* liquid_underflow.

Additional Constraints
----------------------

The ``Thickener0D`` model adds the following constraints to calcuate the split fractions of the solid and liquid stream.

The flux density at the overflow and underflow are calcuated using the following constraint:

If :math:`0 \leqslant \epsilon_{x,t} \leqslant \epsilon_{max}`:

.. math:: F_{x,t} = v_0 \times \epsilon_{x, t} \times(1- \frac{\epsilon_{x,t}}{\epsilon_{max}})^C + v_1 \times \epsilon_{x,t}^2 \times(\epsilon_{max}-\epsilon_{x,t})

otherwise:

.. math:: F_{x,t} = 0

where :math:`x` represents the overflow or underflow, :math:`F` is the flux density, :math:`v_0` is the Stokes velocity of a single particle, :math:`\epsilon` is the solid volume fraction, :math:`\epsilon_{max}` os the maximum attainable solid volume fraction, and :math:`v_1` and :math:`C` are empirical constants [2].




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
