Feedwater Heater (Condensing Section 0D)
========================================

.. index::
  pair: idaes.power_generation.unit_models.feedwater_heater_0D; FWHCondensing0D

.. module:: idaes.power_generation.unit_models.feedwater_heater_0D

The condensing feedwater heater is the same as the
:ref:`HeatExchanger <technical_specs/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>`
model with one additional constraint to calculate the inlet flow rate such that all the
entering steam is condensed.  This model is suitable for steady state modeling, and is
intended to be used with the :ref:`IAWPS95 <technical_specs/model_libraries/generic/property_models/iapws95:International Association of the Properties of Water and Steam IAPWS-95>` 
property package.  For dynamic modeling, the 1D feedwater heater models should be used
(not yet publicly available).

Degrees of Freedom
------------------

Usually ``area`` and ``overall_heat_transfer_coefficient`` are fixed or constraints are provided to calculate ``overall_heat_transfer_coefficient``.  If the inlets are also fixed except for the inlet steam flow rate (``inlet_1.flow_mol``), the model will have 0 degrees of freedom.

Variables
---------

The variables are the same as :ref:`HeatExchanger <technical_specs/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>`.

Constraints
-----------

In addition to the :ref:`HeatExchanger <technical_specs/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>` constraints, there is one additional constraint to calculate the inlet steam flow such that all steam condenses. The constraint is called ``extraction_rate_constraint``, and is defined below.

.. math::

  h_{steam, out} = h_{sat, liquid}(P)

Where :math:`h` is molar enthalpy, and the saturated liquid enthalpy is a function of pressure.


FWHCondensing0D Class
---------------------

.. autoclass:: FWHCondensing0D
  :members:

FWHCondensing0DData Class
--------------------------

.. autoclass:: FWHCondensing0DData
  :members:
