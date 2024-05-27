Cubic Smooth Vapor-Liquid Equilibrium Formulation (SmoothVLE2)
==============================================================

.. contents:: Contents 
    :depth: 2

.. note::
  This formulation for vapor-liquid equilibrium is only valid if using a cubic Equation of State. For other equations of state, use :ref:`SmoothVLE <explanations/components/property_package/general/pe/smooth_flash:Smooth Vapor-Liquid Equilibrium Formulation (SmoothVLE)>`.

Source
------

Dabadghao, V., Ghouse, J., Eslick, J., Lee, A., Burgard, A., Miller, D., and Biegler, L., A Complementarity-based Vapor-Liquid Equilibrium Formulation for Equation-Oriented Simulation and Optimization, AIChE Journal, 2023, Volume 69(4), e18029. https://doi.org/10.1002/aic.18029

Introduction
------------

Typically, equilibrium calculations are only used when the user knows the current state is within the two-phase envelope. For simulation only studies, the user may know a priori the condition of the stream but when the same set of equations are used for optimization, there is a high probability that the specifications can transcend the phase envelope. In these situations, the equilibrium calculations become trivial, thus it is necessary to find a formulation that has non-trivial solutions at all states.

To address this, the cubic smooth vapor-liquid equilibrium (VLE) formulation always solves the equilibrium calculations at a condition where a valid two-phase solution exists. In situations where only a single phase is present, the phase equilibrium is solved at the either the bubble or dew point, where the non-existent phase exists but in negligible amounts. In this way, a non-trivial solution is guaranteed but still gives near-zero material in the non-existent phase in the single phase regions. Rather than explicitly calculate the bubble and dew points (as is done in the non-cubic formulation), this formulation leverages properties of the  cubic equation of state to identify the "equilibrium temperature".

Formulation
-----------

.. note::
  For the full derivation of the cubic smooth VLE formulation, see the reference above.

.. note::
  For consistency of naming between the cubic and non-cubic formulations, :math:`\bar{T}` is referred to as :math:`T_{eq}` in this document and the resulting model.

The approach used by the smooth VLE formulation is to define an "equilibrium temperature" (:math:`T_{eq}`) at which the equilibrium calculations will be performed. The equilibrium temperature is defined such that:

.. math:: T = T_{eq} - s_{vap} + s_{liq}

where :math:`T` is the actual state temperature, and :math:`s_{liq}` and :math:`s_{vap}` are non-negative slack variables. For systems existing the the liquid-only region, :math:`s_{liq}` will be non-zero whilst :math:`s_{vap}=0` (indicating that the system is below the bubble point and thus :math:`T_{eq}>T`). Similarly, for systems in the vapor-only region, :math:`s_{vap}` will be non-zero whilst :math:`s_{liq}=0`. Finally, in the two-phase region, :math:`s_{liq}=s_{vap}=0`, indicating that :math:`T_{eq}=T`.

In order to determine the values of :math:`s_{liq}` and :math:`s_{vap}`, the following complementarity constraints are written:

.. math:: 0 = \min(s_{liq}, F_{liq})
.. math:: 0 = \min(s_{vap}, F_{vap})

where :math:`F_{p}` is the flow rate of each phase :math:`p`. That is, for each phase (liquid and vapor), if thee is any flowrate associated with that phase (i.e., the phase exists), its slack variable must be equal to zero.

Additionally, the follow complementarities are written to constraint the roots of the cubic equation of state.

.. math:: 0 = \min(g^{+}_{liq}, F_{liq})
.. math:: 0 = \min(g^{-}_{vap}, F_{vap})

where :math:`g^{+}_p` and :math:`g^{-}_p` are another pair of non-negative slack variables associated with each phase :math:`p`. These slack variables are defined such that:

.. math:: f''(Z_p) = g^{+}_{p} - g^{-}_{p}

where :math:`f''(Z_p)` is the second derivative of the cubic equation of state written in terms of the compressibility factor :math:`Z_p` for each phase :math:`p`.

Smooth Approximation
''''''''''''''''''''

In order to express the minimum operators in a tractable form, these equations are reformulated using the IDAES `smooth_min` function:

.. math:: \min(a, b) =  0.5{\left[a + b - \sqrt{(a-b)^2 + \epsilon^2}\right]}

Each complementarity requires a smoothing parameter, named :math:`\epsilon_T` and :math:`\epsilon_Z` for the temperature and cubic root constraints respectively. Within the IDAES model, these are rendered as ``eps_t_phase1_phase2`` and ``eps_z_phase1_phase2``, where ``phase1`` and ``phase2`` are the names assigned to the liquid and vapor phases in the property package (order will depend on the order these are declared).

The tractability of the VLE problem depends heavily upon the values chosen for :math:`\epsilon_T` and :math:`\epsilon_Z`, with larger values resulting in smoother transitions at the phase boundaries (and thus increased tractability) at the expense of decreased accuracy near these points. It is recommended that users employ a 2-stage approach to solving these problems, starting with a larger value of :math:`\epsilon_T` and :math:`\epsilon_Z` initially to determine which region the solution lies in, followed by a second solve using smaller values to refine the solution.

As a rule of thumb, the values of :math:`\epsilon_T` and :math:`\epsilon_Z` should be between 2 and 4 orders of magnitude smaller than the largest quantify involved in the smooth maximum operation. This means the value of :math:`\epsilon_T` should be based on the larger of :math:`T` and :math:`F_p`, whilst :math:`\epsilon_Z` should be based on the larger of :math:`f''(Z_p)` and :math:`F_p`. The value of :math:`f''(Z_p)` may be difficult to determine *a priori*, however :math:`F_p` is likely to dominate in most cases unless :math:`F_p` is small or :math:`P` is large.

