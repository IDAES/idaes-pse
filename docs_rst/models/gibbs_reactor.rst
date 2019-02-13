Gibbs Reactor
=============

The IDAES Gibbs reactor model represents a unit operation where a material stream undergoes some set of reactions such that the Gibbs energy of the resulting mixture is minimized. Gibbs reactors rely on conservation of individual elements within the system, and thus require element balances, and make use of Lagrange multipliers to find the minimum Gibbs energy state of the system.

Degrees of Freedom
------------------

Gibbs reactors generally have between 0 and 2 degrees of freedom, depending on construction arguments.

Typical fixed variables are:

* reactor heat duty (has_heat_transfer = True only).
* reactor pressure change (has_pressure_change = True only).

Model Structure
---------------

The core Gibbs reactor unit model consists of a single ControlVolume0DBlock (named control_volume) with one Inlet Port (named inlet) and one Outlet Port (named outlet).

Variables
---------

Gibbs reactor units add the following additional Variables beyond those created by the Control Volume Block.

=============== ================== =============================================
Variable Name   Symbol             Notes
=============== ================== =============================================
lagrange_mult   :math:`L_{t,e}`    Lagrange multipliers
heat_duty       :math:`Q_t`        Only if has_heat_transfer = True, reference
deltaP          :math:`\Delta P_t` Only if has_pressure_change = True, reference
=============== ================== =============================================

Constraints
-----------

Gibbs reactor models write the following additional constraints to calculate the state that corresponds to the minimum Gibbs energy of the system.

`gibbs_minimization(time, phase, component)`:

.. math :: 0 = g_{partial,t,j} + \sum_e{L_{t,e} \times \alpha_{j,e}})

where :math:`g_{partial,t,j}` is the partial molar Gibbs energy of component :math:`j` at time :math:`t`, :math:`L_{t,e}` is the Lagrange multiplier for element :math:`e` at time :math:`t` and :math:`\alpha_{j,e}` is the number of moles of element :math:`e` in one mole of component :math:`j`. :math:`g_{partial,t,j}` and :math:`\alpha_{j,e}` come from the outlet StateBlock.

GibbsReactor Class
------------------

.. module:: idaes.unit_models.gibbs_reactor

.. autoclass:: GibbsReactor
  :members:

GibbsReactorData Class
----------------------

.. autoclass:: GibbsReactorData
  :members:
