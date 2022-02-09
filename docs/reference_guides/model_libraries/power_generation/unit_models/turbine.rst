Turbine (Isentropic)
====================

.. index::
  pair: idaes.power_generation.unit_models.helm.turbine_stage;HelmIsentropicTurbine

.. module:: idaes.power_generation.unit_models.helm.turbine

This is a steam power generation turbine model for the basic isentropic turbine
calculations. It is the basis of the :ref:`TurbineInletStage
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Inlet Stage)>`,
`TurbineOutletStage
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Outlet Stage)>`,
and,
`TurbineOutletStage
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Stage)>`
models.

Variables
---------

=========================== ======================== =========== ==========================================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ==========================================================================================
``efficiency_isentropic``   :math:`\eta_{isen}`      time        Isentropic efficiency
``deltaP``                  :math:`\Delta P`         time        Pressure change (:math:`P_{out} - P_{in}`) [Pa]
``ratioP``                  :math:`P_{ratio}`        time        Ratio of discharge pressure to inlet pressure :math:`\left(\frac{P_{out}}{P_{in}}\right)`
=========================== ======================== =========== ==========================================================================================

Expressions
-----------

This model provides two expressions that are not available in the
pressure changer model.

=========================== ========================= =========== ======================================================================
Expression                  Symbol                    Index Sets  Doc
=========================== ========================= =========== ======================================================================
``h_is``                    :math:`h_{is}`            time        Isentropic outlet molar enthalpy [J/mol]
``delta_enth_isentropic``   :math:`\Delta h_{is}`     time        Isentropic enthalpy change (:math:`h_{is} - h_{in}`) [J/mol]
``work_isentropic``         :math:`w_{is}`            time        Isentropic work (W)
=========================== ========================= =========== ======================================================================

Constraints
-----------

In addition to the mass and energy balances provided by the control volume the
following equation is used to calculate the outlet enthalpy, so work comes from
the control volume energy balance.

.. math::

  h_{out} = h_{in} - \eta_{is}\left(h_{in} - h_{is}\right)

Initialization
--------------

To initialize the turbine model, a reasonable guess for the inlet condition and
deltaP and efficiency should be set by setting the appropriate variables.

TurbineStage Class
------------------

.. autoclass:: HelmIsentropicTurbine
  :members:

TurbineStageData Class
----------------------

.. autoclass:: HelmIsentropicTurbineData
  :members:
