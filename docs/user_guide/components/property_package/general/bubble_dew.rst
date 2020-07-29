Bubble and Dew Point Methods
============================

.. contents:: Contents 
    :depth: 3

Ideal Assumptions (``IdealBubbleDew``)
--------------------------------------

In the case where ideal behavior can be assumed, i.e. ideal gas assumption and Raoult's Law holds, the bubble and dew points can be calculated directly from the saturation pressure using the following equations.

Ideal Bubble Pressure
^^^^^^^^^^^^^^^^^^^^^

.. math:: P_{bub} = \sum_j{x_j \times P_{sat, j}(T)}
.. math:: x_j(P_{bub}) \times P_{bub} = x_j \times P_{sat, j}(T)

where :math:`P_{bub}` is the bubble pressure of the mixture, :math:`P_{sat, j}(T)` is the saturation pressure of component :math:`j` at the system temperature, :math:`T`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(P_{bub})` is the mole fraction of the vapor phase at the bubble pressure.

Ideal Bubble Temperature
^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: \sum_j{\left(x_j \times P_{sat, j}(T_{bub})\right)} - P = 0
.. math:: x_j(T_{bub}) \times P = x_j \times P_{sat, j}(T_{bub})

where :math:`P` is the system pressure, :math:`P_{sat, j}(T_{bub})` is the saturation pressure of component :math:`j` at the bubble temperature, :math:`T_{bub}`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(T_{bub})` is the mole fraction of the vapor phase at the bubble temperature.

Ideal Dew Pressure
^^^^^^^^^^^^^^^^^^

.. math:: 0 = 1 - P_{dew} \times \sum_j{x_j \times P_{sat, j}(T)}
.. math:: x_j(P_{dew}) \times P_{sat, j}(T) = x_j \times P_{dew}

where :math:`P_{dew}` is the dew pressure of the mixture, :math:`P_{sat, j}(T)` is the saturation pressure of component :math:`j` at the system temperature, :math:`T`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(P_{dew})` is the mole fraction of the liquid phase at the dew pressure.

Ideal Dew Temperature
^^^^^^^^^^^^^^^^^^^^^

.. math:: P \times \sum_j{\left(x_j \times P_{sat, j}(T_{dew})\right)} - 1 = 0
.. math:: x_j(T_{dew}) \times P_{sat, j}(T_{dew}) = x_j \times P

where :math:`P` is the system pressure, :math:`P_{sat, j}(T_{dew})` is the saturation pressure of component :math:`j` at the dew temperature, :math:`T_{bub}`, :math:`x_j` is the overall mixture mole fraction and :math:`y_j(T_{dew})` is the mole fraction of the liquid phase at the dew temperature.

Equal Fugacity (log form) (``LogBubbleDew``)
--------------------------------------------

For cases where ideal behavior is insufficient, it is necessary to calculate the fugacity of each component at the relevant transition point and enforce equality of the fugacity in each phase. As such, this methods depends upon the definition of fugacity for each phase and component. In this formulation, the logarithm of the phase equilibrium constraint is used.

Bubble Pressure (log form)
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: ln(x_j) + ln(f_{liquid, j}(P_{bub})) = ln(x_j(P_{bub})) + ln(f_{vapor, j}(P_{bub}))
.. math:: 1 = \sum_j{x_j(P_{bub})}

where :math:`P_{bub}` is the bubble pressure of the mixture, :math:`f_{p, j}(P_{bub})` is the fugacity of component :math:`j` in phase :math:`p` at :math:`P_{bub}`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(P_{bub})` is the mole fraction of the vapor phase at the bubble pressure. 

Bubble Temperature (log form)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: ln(x_j) + ln(f_{liquid, j}(T_{bub})) = ln(x_j(T_{bub})) + ln(f_{vapor, j}(T_{bub}))
.. math:: 1 = \sum_j{x_j(T_{bub})}

where :math:`T_{bub}` is the bubble temperature of the mixture, :math:`f_{p, j}(T_{bub})` is the fugacity of component :math:`j` in phase :math:`p` at :math:`T_{bub}`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(T_{bub})` is the mole fraction of the vapor phase at the bubble temperature. 

Dew Pressure (log form)
^^^^^^^^^^^^^^^^^^^^^^^

.. math:: ln(x_j(P_{dew})) + ln(f_{liquid, j}(P_{dew})) = ln(x_j) + ln(f_{vapor, j}(P_{dew}))
.. math:: 1 = \sum_j{x_j(P_{dew})}

where :math:`P_{dew}` is the dew pressure of the mixture, :math:`f_{p, j}(P_{dew})` is the fugacity of component :math:`j` in phase :math:`p` at :math:`P_{dew}`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(P_{dew})` is the mole fraction of the vapor phase at the dew pressure. 

Dew Temperature (log form)
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: ln(x_j(T_{dew})) + ln(f_{liquid, j}(T_{dew})) = ln(x_j) + ln(f_{vapor, j}(T_{dew}))
.. math:: 1 = \sum_j{x_j(T_{dew})}

where :math:`T_{dew}` is the dew temperature of the mixture, :math:`f_{p, j}(T_{dew})` is the fugacity of component :math:`j` in phase :math:`p` at :math:`T_{dew}`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(T_{dew})` is the mole fraction of the vapor phase at the dew temperature. 
