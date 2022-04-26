Plug Flow Reactor
=================

The IDAES Plug Flow Reactor (PFR) model represents a unit operation where a material stream passes through a linear reactor vessel whilst undergoing some chemical reaction(s). This model requires modeling the system in one spatial dimension.

Degrees of Freedom
------------------

PFRs generally have at least 2 degrees of freedom.

Typical fixed variables are:

* 2 of reactor length, area and volume.

Model Structure
---------------

The core PFR unit model consists of a single ControlVolume1DBlock (named control_volume) with one Inlet Port (named inlet) and one Outlet Port (named outlet).

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

PFR units write the following additional Constraints at all points in the spatial domain:

.. math:: X_{t,x,r} = A \times r_{t,x,r}

where :math:`X_{t,x,r}` is the extent of reaction of reaction :math:`r` at point :math:`x` and time :math:`t`, :math:`A` is the cross-sectional area of the reactor and :math:`r_{t,r}` is the volumetric rate of reaction of reaction :math:`r` at point :math:`x` and time :math:`t` (from the outlet StateBlock).

PFR Class
---------

.. module:: idaes.models.unit_models.plug_flow_reactor

.. autoclass:: PFR
    :members:

PFRData Class
-------------

.. autoclass:: PFRData
    :members:
