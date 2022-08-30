Coordinator
==================
To ensure the energy supply equals demand at every moment across the transmission
network, Independent System Operators (ISOs) clear the markets seeks to minimize
the overall operation and production costs by coordinating the operation schedules
of generation units in the network while considering operational and physical
constraints of power systems. As mentioned above, the market is cleared by a
two-settlement system: the day-ahead market and the real-time market. In the framework,
the ISOs are modeled by the PCM Prescient. The ``DoubleLoopCoordinator`` will
coordinate the operations among ``Bidder``, ``Tracker``, and Prescient for the
simulations.

.. module:: idaes.apps.grid_integration.coordinator

.. autoclass:: DoubleLoopCoordinator
  :members:
