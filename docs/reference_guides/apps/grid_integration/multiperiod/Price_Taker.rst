Price Taker
===========

Price takers are companies or individuals who must accept market prices since they lack the market share
to directly influence the market price. Likewise, it is assumed that a price taker's resource or energy
system is small enough such that it does not significantly impact the market. When coupled with multi-period modeling,
the price-taker model is able to synthesize grid-centric modeling with steady-state process-centric modeling, as
depicted in figure below.

.. |pricetaker| image:: images/pricetaker.png
  :width: 1200
  :alt: Alternative text
  :align: middle

|pricetaker|

The following equations represent the multi-period price taker model, where :math:`d` are design decisions,
:math:`u` are operating decisions, :math:`δ`: are power decisions ,:math:`s` are scenarios (timeseries/representative days),
:math:`w` is weight/frequency, :math:`R` is revenue, :math:`π` is price data,
:math:`C` is capital and operating costs, :math:`g` is the process model, and :math:`h` is the temporal constraint.

    .. math::

       max_{d,u, δ} = \sum_{s ∈ S} \sum{t ∈ T} w_{s}[R(d,u_{s,t},δ_{s,t},π_{s,t}) - C(d,u_{s,t},δ_{s,t})]

       g(d,u_{s,t},δ_{s,t}) = 0
       ∀_{s} ∈ S, t ∈ T

       h(d,u_{s,t},δ_{s,t},u_{s,t+1},δ_{s,t+1}) = 0
       ∀_{s} ∈ S, t ∈ T


The price taker multi-period modeling workflow involves the integration of multiple software platforms into the IDAES optimization model
and can be broken down into two distinct functions, as shown in the figure below. In part 1, simulated or historical
ISO (International Organization for Standardization) data is used to generate locational marginal price (LMP)
signals, and production cost models (PCMs) are used to compute and optimize the time-varying dispatch schedules for each
resource based on their respective bid curves. Advanced data analytics (RAVEN) reinterpret the LMP signals and PCM
as stochastic realizations of the LMPs in the form of representative days (or simply the full-year price signals).
In part 2, PRESCIENT uses a variety of input parameters (design capacity, minimum power output, ramp rate, minimum up/down time, marginal cost, no load cost, and startup profile)
to generate data for the market surrogates. Meanwhile, IDAES uses the double loop simulation to integrate detailed
process models (b, ii) into the daily (a, c) and hourly (i, iii) grid operations workflow.

.. |hybrid_energy_system| image:: images/hybrid_energy_system.png
  :width: 1200
  :alt: Alternative text
  :align: middle

|hybrid_energy_system|