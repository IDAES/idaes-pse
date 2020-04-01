Moving Bed Reactor
==================

The IDAES Moving Bed Reactor (MBR) model represents a unit operation where two material streams – a solid phase and a gas phase – pass through a linear reactor vessel while undergoing chemical reaction(s). The two streams have opposite flow directions (counter-flow).
The MBR mathematical model is a 1-D rigorous first-principles model consisting of a set of differential equations obtained by applying the mass, energy and momentum balance equations. The radial concentration and temperature gradients are assumed to be negligible. A mass and an energy balance equation is written for each phase. The reactor is assumed to be adiabatic. The solid phase is assumed to be moving at a constant velocity determined by the solids feed rate to the reactor, while the gas phase velocity through the reactor bed is computed using a simplified pressure balance equation, or the Ergun equation.

The MBR model is based on:

A. Ostace, A. Lee, C.O. Okoli, A.P. Burgard, D.C. Miller, D. Bhattacharyya, Mathematical modeling of a moving-bed reactor for chemical looping combustion of methane, in: M.R. Eden, M. Ierapetritou, G.P. Towler (Eds.),13th Int. Symp. Process Syst. Eng. (PSE 2018), Computer-Aided Chemical Engineering 2018, pp. 325–330 , San Diego, CA.

Degrees of Freedom
------------------

MBRs generally have at least 2 (or more) degrees of freedom. Typically fixed variables are reactor length and diameter.

Model Structure
---------------

The core MBR unit model consists of two ControlVolume1DBlock Blocks (named gas_phase and solid_phase), each with one Inlet Port (named gas_inlet and solid_inlet) and one Outlet Port (named gas_outlet and solid_outlet).

Variables
---------

PFR units add the following additional Variables:

====================== ======= ===============================================================
Variable               Name    Notes
====================== ======= ===============================================================
:math:`L`              length  Reference to control_volume.length
:math:`A`              area    Reference to control_volume.area
:math:`V`              volume  Reference to control_volume.volume
:math:`Q_{t,x}`        heat    Only if has_heat_transfer = True, reference to holdup.heat
:math:`\Delta P_{t,x}` deltaP  Only if has_pressure_change = True, reference to holdup.deltaP
====================== ======= ===============================================================

Constraints
-----------

In addition to the constraints written by the control_volume Block, MBR units write the following Constraints at all points in the spatial domain:

.. math:: X_{t,x,r} = A \times r_{t,x,r}

where :math:`X_{t,x,r}` is the extent of reaction of reaction :math:`r` at point :math:`x` and time :math:`t`, :math:`A` is the cross-sectional area of the reactor and :math:`r_{t,r}` is the volumetric rate of reaction of reaction :math:`r` at point :math:`x` and time :math:`t` (from the outlet StateBlock).

The constraints written by the MBR model to compute the pressure drop in the reactor depend upon the construction arguments chosen:

If `pressure_drop_type` is `simple_correlation`:

.. math:: - \frac{ dP }{ dz } = \left( \rho_{p} - \rho_{g} \right) a_{E} u_{g}

where :math:`P` is the system pressure, :math:`z` is the spatial domain, :math:`\rho_{p}` is the density of the particles, :math:`\rho_{g}` is the density of the gas, :math:`u_{g}` is the superficial velocity of the gas, and :math:`a_{E}` is a parameter having the value of 0.2.

If `pressure_drop_type` is `ergun_correlation`:

.. math:: - \frac{ dP }{ dz } = \frac{ 150 \mu_{g} {\left( 1 - \epsilon \right)}^{2} \left( u_{g} + u_{s} \right) }{ \epsilon^{3} d_{p}^2 } + \frac{ 1.75 \left( 1 - \epsilon \right) \rho_{g} \left( u_{g} + u_{s} \right)^{2} }{ \epsilon^{3} d_{p} }

where :math:`\mu_{g}` is the gas viscosity, :math:`\epsilon` is the reactor bed voidage, :math:`u_{s}` is the velocity of the solids, and :math:`d_{p}` is the diameter of the solid particles.

Initialization
--------------

The initialization method for this model uses a multi-step, sequential, hierarchical initialization approach.

MBR Class
---------

.. module:: idaes.gas_solid_contactors.unit_models.moving_bed

.. autoclass:: MBR
    :members:

MBRData Class
-------------

.. autoclass:: MBRData
    :members:
