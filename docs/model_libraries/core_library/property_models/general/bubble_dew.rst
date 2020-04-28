Bubble and Dew Point Methods
============================

.. contents:: Contents 
    :depth: 3

Ideal Assumptions
-----------------

In the case where ideal behavior can be assumed, i.e. Raoult's Law holds, the bubble and dew points can be calculated directly from the saturation pressure using the following equations.

Ideal Bubble Pressure
^^^^^^^^^^^^^^^^^^^^^

This method is implemented as `bubble_press_ideal`.

.. math:: P_{bub} = \sum_j{x_j \times P_{sat, j}(T)}
.. math:: x_j(P_{bub}) \times P_{bub} = x_j \times P_{sat, j}(T)

where :math:`P_{bub}` is the bubble pressure of the mixture, :math:`P_{sat, j}(T)` is the saturation pressure of component :math:`j` at the system temperature, :math:`T`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(P_{bub})` is the mole fraction of the vapor phase at the bubble pressure.

Ideal Bubble Temperature
^^^^^^^^^^^^^^^^^^^^^^^^

This method is implemented as `bubble_temp_ideal`.

.. math:: \sum_j{\left(x_j \times P_{sat, j}(T_{bub})\right)} - P = 0
.. math:: x_j(T_{bub}) \times P = x_j \times P_{sat, j}(T_{bub})

where :math:`P` is the system pressure, :math:`P_{sat, j}(T_{bub})` is the saturation pressure of component :math:`j` at the bubble temperature, :math:`T_{bub}`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(T_{bub})` is the mole fraction of the vapor phase at the bubble temperature.

Ideal Dew Pressure
^^^^^^^^^^^^^^^^^^

This method is implemented as `dew_press_ideal`.

.. math:: 0 = 1 - P_{dew} \times \sum_j{x_j \times P_{sat, j}(T)}
.. math:: x_j(P_{dew}) \times P_{sat, j}(T) = x_j \times P_{dew}

where :math:`P_{dew}` is the dew pressure of the mixture, :math:`P_{sat, j}(T)` is the saturation pressure of component :math:`j` at the system temperature, :math:`T`, :math:`x_j` is the overall mixture mole fraction and :math:`x_j(P_{dew})` is the mole fraction of the liquid phase at the dew pressure.

Ideal Dew Temperature
^^^^^^^^^^^^^^^^^^^^^

This method is implemented as `dew_temp_ideal`.

.. math:: P \times \sum_j{\left(x_j \times P_{sat, j}(T_{dew})\right)} - 1 = 0
.. math:: x_j(T_{dew}) \times P_{sat, j}(T_{dew}) = x_j \times P

where :math:`P` is the system pressure, :math:`P_{sat, j}(T_{dew})` is the saturation pressure of component :math:`j` at the dew temperature, :math:`T_{bub}`, :math:`x_j` is the overall mixture mole fraction and :math:`y_j(T_{dew})` is the mole fraction of the liquid phase at the dew temperature.
