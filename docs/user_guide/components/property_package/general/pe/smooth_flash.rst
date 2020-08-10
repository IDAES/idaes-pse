Smooth Vapor-Liquid Equilibrium Formulation (``smooth_VLE``)
============================================================

.. contents:: Contents 
    :depth: 2

Source
------

Burgard, A.P., Eason, J.P., Eslick, J.C., Ghouse, J.H., Lee, A., Biegler, L.T., Miller, D.C., 2018, A Smooth, Square Flash Formulation for Equation-Oriented Flowsheet Optimization. Proceedings of the 13th International Symposium on Process Systems Engineering – PSE 2018, July 1-5, 2018, San Diego.

Introduction
------------

Typically, equilibrium calculations are only used when the user knows the current state is within the two-phase envelope. For simulation only studies, the user may know a priori the condition of the stream but when the same set of equations are used for optimization, there is a high probability that the specifications can transcend the phase envelope. In these situations, the equilibrium calculations become trivial, thus it is necessary to find a formulation that has non-trivial solutions at all states.

To address this, the smooth vapor-liquid equilibrium (VLE) formulation always solves the equilibrium calculations at a condition where a valid two-phase solution exists. In situations where only a single phase is present, the phase equilibrium is solved at the either the bubble or dew point, where the non-existent phase exists but in negligible amounts. In this way, a non-trivial solution is guaranteed but still gives near-zero material in the non-existent phase in the single phase regions.

Formulation
-----------

The approach used by the smooth VLE formulation is to define an "equilibrium temperature" (:math:`T_{eq}`) at which the equilibrium calculations will be performed. The equilibrium temperature is computed as follows:

.. math:: T_{1} = max(T_{bubble}, T) 
.. math:: T_{eq} = min(T_{1}, T_{dew})

where :math:`T` is the actual stream temperature, :math:`T_{1}` is an intermediate temperature variable and :math:`T_{bubble}` and :math:`T_{dew}` are the bubble and dew point temperature of mixture. In order to express the maximum and minimum operators in a tractable form, these equations are reformulated using the IDAES `smooth_max` and `smooth_min` operators which results in the following equations:

.. math:: T_{1} = 0.5{\left[T + T_{bubble} + \sqrt{(T-T_{bubble})^2 + \epsilon_{1}^2}\right]}
.. math:: T_{eq} = 0.5{\left[T_{1} + T_{dew} - \sqrt{(T-T_{dew})^2 + \epsilon_{2}^2}\right]}

where :math:`\epsilon_1` and :math:`\epsilon_2` are smoothing parameters(mutable `Params` named `eps_1` and `eps_2`). The default values are 0.01 and 0.0005 respectively, and it is recommended that :math:`\epsilon_1` > :math:`\epsilon_2`. It can be seen that if the stream temperature is less than that of the bubble point temperature, the VLE calculations will be computed at the bubble point. Similarly, if the stream temperature is greater than the dew point temperature, then the VLE calculations are computed at the dew point temperature. For all other conditions, the equilibrium calculations will be computed at the actual temperature.

Finally, the phase equilibrium is expressed using the following equation:

.. math:: \Phi_{\text{Vap}, j}(T_{eq}) = \Phi_{\text{Liq}, j}(T_{eq})

where :math:`\Phi_{p, j}(T_{eq})` is the fugacity of component :math:`j` in the phase :math:`p` calculated at :math:`T_{eq}`. The fugacities are calculated using methods defined by the equation of state chosen by the user for each phase.
