One-dimensional membrane class for CO2 gas separation
================================================================

This is a one-dimensional model for gas separation in CO₂ capture applications.
The model will be discretized in the flow direction, and it supports two flow patterns:
counter-current flow and co-current flow. The model was customized for gas-phase separation
in CO₂ capture with a single-layer design. If a multi-layer design is needed, multiple units
can be connected for this application. The two sides of the membrane are called the feed side
and sweep side. The sweep stream inlet is optional. The driving force across the membrane is the
partial pressure difference in this gas separation application. Additionally, the energy balance
assumes that temperature remains constant on each side of the membrane.

Variables
---------

Model Inputs - symbol:

* Membrane length - :math:`L`
* Membrane Area - :math:`A`
* Permeance - :math:`per`
* Feed flowrate - :math:`F_fr`
* Feed compositions - :math:`x`
* Feed pressure - :math:`P`
* Feed temperature - :math:`T`


Model Outputs :

* Permeate compositions
* Permeate flowrate

Degrees of Freedom
------------------

The DOF should be 0 for square problem simulations.




