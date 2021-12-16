.. grid_integration documentation master file, created by
   sphinx-quickstart on Wed Dec 15 17:14:39 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to grid_nmpc's documentation!
============================================

.. |doubleloop| image:: images/framework.png
  :width: 400
  :alt: Alternative text
  :align: middle

|doubleloop|

The subpackage implements a new generalized multiscale simulation framework, shown in
figure, that integrates the process and grid modeling paradigms to quantify
the operational (hours to year timescale) interactions between energy systems and
wholesale electricity markets. We integrate optimal design, operations, and
control of an  energy system with market clearing via a high-fidelity Production
Cost Model (PCM) Prescient that optimizes resource dispatch decisions across an
entire transmission network. While optimal bidding and scheduling of energy system,
market clearing and dispatch, and advanced control of the energy systems have been
studied in isolation previously, this framework bridges these process- and
grid-centric modeling paradigms. Specifically, in the day-ahead market loop,
before the market is cleared, the individual energy systems forecast the market
uncertainties, e.g., energy locational marginal prices (LMPs), weather, fuel
prices, renewable generation. Any accurate forecasting method is compatible with
the framework. These forecasts then facilitate the energy systems to compute optimal
market bids under uncertainty. Energy bids are usually marginal production costs
of energy systems, and they are submitted to the market operator and enter the objective
functions of market-clearing problems. The day-ahead market (figure right) is cleared by
solving unit commitment problems. Similarly, in the real-time market cycle
(figure left), economic dispatch problems are solved every
hour to clear the market and set generation targets for each resource. Each
energy system then solves optimal control problems to track these dispatch signals from
real-time markets constrained by (non)linear process models. Finally, the market
calculates the settlements and pays the energy systems.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Bidder
   Tracker
   Coordinator



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
