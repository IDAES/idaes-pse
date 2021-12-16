Bidder
============================================
Market participating resources (e.g., generators, IESs) submit energy bids
(a.k.a., bid curves) to the day-ahead and real-time markets for each trading time
period to communicate their flexibility and marginal costs. An energy bid is a
piecewise constant function described by several energy offer price (\$/MWh) and
operating level (MW) pairs. Bid curves from each resource are inputs (i.e.,
parameters) in the market-clearing optimization problems solved by PCM. Currently,
the ``Bidder`` formulates a two-stage stochastic program to calculate the optimized
time-varying bid curves for thermal generators.

.. autoclass:: bidder.Bidder
  :members:
