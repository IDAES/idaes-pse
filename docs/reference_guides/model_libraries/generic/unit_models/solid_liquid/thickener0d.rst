Thickener (0D)
==============

The ``Thickener0D`` unit model is a predictive model for clarifiers and thickeners based on the following references:

[1] R. Burger, F. Concha, K.H. Karlsen, A. Narvaez, Numerical simulation of clarifier-thickener units treating ideal suspensions with a flux density function having two inflection points, Mathematical and Computer Modelling 44 (2006) 255–275, doi:10.1016/j.mcm.2005.11.008

[2] N.G. Barton, C.-H. Li, S.J. Spencer, Control of a surface of discontinuity in continuous thickeners, J. Aust. Math. Soc. Ser. B 33 (1992) 269–289


Degrees of Freedom
------------------

The ``Thickener0D`` model has 6 degrees of freedom. Four of these are used to define the solids settling velocity and must be provided by the user.

* the solid particle size ``particle_size``,
* the maximum achievable solid fraction ``solid_fraction_max``,
* two empirical coefficients for calculating the settling velocity ``v1``, ``C``.

Additionally, the user must provide 2 additional degrees of freedom regarding the design and/or operation of the thickener such as:

* the cross-sectional area of the settler ``area``,
* the total volumetric flowrate or volume fraction of solids at the underflow, or
* the total volumetric flowrate or volume fraction of solids at the overflow.

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

The ``Thickener0D`` model adds the following constraints to calculate the split fractions of the solid and liquid stream.

Total volumetric flowrates at the feed  overflow and underflow are calculated as the sum of the volumetric flowrates of the solid and liquid streams at each location.

.. math:: Q_{x,t} = Q_{solid,x,t} + Q_{liquid,x,t}

where :math:`Q` represents volumetric flowrate and :math:`x` represents the feed, overflow or underflow point.

The flux density at the overflow and underflow are calculated using the following constraint:

If :math:`0 \leqslant \epsilon_{x,t} \leqslant \epsilon_{max}`:

.. math:: F_{x,t} = v_0 \times \epsilon_{x, t} \times(1- \frac{\epsilon_{x,t}}{\epsilon_{max}})^C + v_1 \times \epsilon_{x,t}^2 \times(\epsilon_{max}-\epsilon_{x,t})

otherwise:

.. math:: F_{x,t} = 0

where :math:`F` is the flux density, :math:`v_0` is the Stokes velocity of a single particle, :math:`\epsilon` is the solid volume fraction, :math:`\epsilon_{max}` is the maximum attainable solid volume fraction, and :math:`v_1` and :math:`C` are empirical constants [2].

The solids fraction at the overflow and underflow is described using the following equation (derived from [1]).

.. math:: Q_{feed, t} \times \epsilon_{feed,t} = A \times (F_{overflow,t} + F_{underflow,t}) - Q_{overflow,t} \times (\epsilon_{overflow,t} - \epsilon_{feed,t}) + Q_{underflow,t} \times (\epsilon_{underflow,t} - \epsilon_{feed,t})

where :math:`A` is the cross-sectional area of the thickener (assumed constant in space and time).

Conservation of solids is enforced using the following constraint:

.. math:: Q_{feed,t} \times \epsilon_{feed,t} = Q_{overflow,t} \times \epsilon_{overflow,t} + Q_{underflow,t} \times \epsilon_{underflow,t}

Two constraints are written to define the solid volume fraction at the feed and underflow point (no constraint is written for the overflow point).

.. math:: Q_{solid,x,t} = \epsilon_{x,t} \times (Q_{solid,x,t} + Q_{liquid,x,t})

The solids volume fraction at the overflow and underflow are bounded by the following inequality constraints:

.. math:: \epsilon_{t,x} \leqslant \epsilon_{max}

Finally, the Stokes velocity is calculated using Stokes Law:

.. math:: 18 \times v0_t \times \mu_{liquid,t} = (\rho_{solid,t} - \rho_{liquid,t}) \times g \times d_p^2

where :math:`\mu_{liquid}` is the viscosity of the liquid phase, :math:`\rho_{liquid}` and :math:`\rho_{solid}` are the densities of the liquid and solid phases, :math:`g` is the acceleration due to gravity and :math:`d_p` is the solid particle diameter.

Variables
---------

The ``Thickener0D`` adds the following variables to those contained in the Separators.

============================= ========================= ====== ========================================================
Variable                      Name                      Index  Notes
============================= ========================= ====== ========================================================
:math:`A`                     area                      None   Cross-sectional area of thickener
:math:`Q_{feed}`              flow_vol_feed             time   Total volumetric flowrate at feed point
:math:`Q_{overflow}`          flow_vol_overflow         time   Total volumetric flowrate at overflow point
:math:`Q_{underflow}`         flow_vol_underflow        time   Total volumetric flowrate at underflow point
:math:`\epsilon_{feed}`       solid_fraction_feed       time   Volumetric solids fraction at feed point
:math:`\epsilon_{overflow}`   solid_fraction_overflow   time   Volumetric solids fraction at overflow point
:math:`\epsilon_{underflow}`  solid_fraction_underflow  time   Volumetric solids fraction at underflow point
:math:`F_{overflow}`          flux_density_overflow     time   Solids flux density at overflow point
:math:`F_{underflow}`         flux_density_underflow    time   Solids flux density at underflow point
:math:`d_p`                   particle size             time   Solid particle size
:math:`v0`                    v0                        time   Stokes velocity of isolated particle
:math:`v1`                    v1                        None   Empirical parameter in settling velocity correlation
:math:`C`                     C                         None   Empirical parameter in settling velocity correlation
math:`\epsilon_{max}`         solid_fraction_max        None   Maximum achievable volumetric solids fraction
============================= ========================= ====== ========================================================

.. module:: idaes.models.unit_models.solid_liquid.thickener

Thickener0D Class
-----------------

.. autoclass:: Thickener0D
  :members:

Thickener0DData Class
---------------------

.. autoclass:: Thickener0DData
  :members:
