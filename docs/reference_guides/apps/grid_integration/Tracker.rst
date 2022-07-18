Tracker
============
Production Cost Models (PCM) computes the time-varying dispatch schedules for each resource using
simplified models. Emerging resources including IESs and industrial demand response
need to determine the optimal operations strategy to track their market dispatch
signal. The ``Tracker`` formulates these decisions as a model predictive control
(MPC) problem. The figure below shows an example of the optimal tracking from an
integrated energy system which consists of a thermal generator and an energy storage.
The figure shows that to track the dispatch (load) the energy system can optimally
use power output from charging and discharging cycle.

.. |tracking_example| image:: images/tracking_example.png
  :width: 1200
  :alt: Alternative text
  :align: middle

|tracking_example|

.. module:: idaes.apps.grid_integration.tracker

.. autoclass:: Tracker
  :members:
