Turbine (Stage)
===============

.. index::
  pair: idaes.power_generation.unit_models.helm.turbine_stage;HelmTurbineStage

.. module:: idaes.power_generation.unit_models.helm.turbine_stage

This is a steam power generation turbine model for intermediate stages between the
inlet and outlet.  It inherits `HelmIsentropicTurbine
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Isentropic)>`.

Example
-------

.. testcode::

    from pyomo.environ import ConcreteModel, SolverFactory

    from idaes.core import FlowsheetBlock
    from idaes.power_generation.unit_models.helm import HelmTurbineStage
    from idaes.generic_models.properties import iapws95

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = HelmTurbineStage(default={"property_package": m.fs.properties})
    # set inlet
    m.fs.turb.inlet[:].enth_mol.fix(70000)
    m.fs.turb.inlet[:].flow_mol.fix(15000)
    m.fs.turb.inlet[:].pressure.fix(8e6)
    m.fs.turb.efficiency_isentropic[:].fix(0.8)
    m.fs.turb.ratioP[:].fix(0.7)
    m.fs.turb.initialize()


Variables
---------

=========================== ======================== =========== ==========================================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ==========================================================================================
``efficiency_mech``         :math:`\eta_{mech}`      None        Mechanical efficiency (accounts for losses in bearings...)
``efficiency_isentropic``   :math:`\eta_{isen}`      time        Isentropic efficiency
``deltaP``                  :math:`\Delta P`         time        Pressure change (:math:`P_{out} - P_{in}`) [Pa]
``ratioP``                  :math:`P_{ratio}`        time        Ratio of discharge pressure to inlet pressure :math:`\left(\frac{P_{out}}{P_{in}}\right)`
``shaft_speed``             :math:`s`                time        Shaft speed [hz]
=========================== ======================== =========== ==========================================================================================

The shaft speed is used to calculate specific speed for more advanced turbine models,
the specific speed expression is available, but otherwise has no effect on the model
results.

Expressions
-----------

This model provides two expressions that are not available in the
pressure changer model.

=========================== ========================= =========== ======================================================================
Variable                    Symbol                    Index Sets  Doc
=========================== ========================= =========== ======================================================================
``power_thermo``            :math:`\dot{w}_{thermo}`  time        Turbine stage power output not including mechanical loss [W]
``power_shaft``             :math:`\dot{w}_{shaft}`   time        Turbine stage power output including mechanical loss (bearings...) [W]
``specific_speed``          :math:`n_s`               time        Turbine stage specific speed [dimensionless]
=========================== ========================= =========== ======================================================================

.. math::

  n_s = s\dot{v}^0.5 (w_{isen} / \dot{m}) ** (-0.75)

Where :math:`\dot{m}` is the mass flow rate and :math:`\dot{v}` is the outlet
volumetric flow.

Constraints
-----------

There are no additional constraints.

Initialization
--------------

To initialize the turbine model, a reasonable guess for the inlet condition and
deltaP and efficiency should be set by setting the appropriate variables.

TurbineStage Class
------------------

.. autoclass:: HelmTurbineStage
  :members:

TurbineStageData Class
----------------------

.. autoclass:: HelmTurbineStageData
  :members:
