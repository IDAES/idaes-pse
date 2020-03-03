Turbine (Stage)
===============

.. index::
  pair: idaes.power_generation.unit_models.turbine_stage;TurbineStage

.. module:: idaes.power_generation.unit_models.turbine_stage

This is a steam power generation turbine model for the stages between the inlet
and outlet.
This model inherits the :ref:`PressureChanger model
<model_libraries/core_library/unit_models/pressure_changer:Pressure Changer>` with the isentropic options. The
initialization scheme is the same as the :ref:`TurbineInletStage model
<model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Inlet Stage)>`.

Example
-------

.. testcode::

    from pyomo.environ import ConcreteModel, SolverFactory

    from idaes.core import FlowsheetBlock
    from idaes.power_generation.unit_models import TurbineStage
    from idaes.generic_models.properties import iapws95

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = TurbineStage(default={"property_package": m.fs.properties})
    # set inlet
    m.fs.turb.inlet[:].enth_mol.fix(70000)
    m.fs.turb.inlet[:].flow_mol.fix(15000)
    m.fs.turb.inlet[:].pressure.fix(8e6)
    m.fs.turb.efficiency_isentropic[:].fix(0.8)
    m.fs.turb.ratioP[:].fix(0.7)
    m.fs.turb.initialize()


Variables
---------

This model adds a variable to the base ``PressureChanger model`` to account
for mechanical efficiency .

=========================== ======================== =========== ======================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ======================================================================
``efficiency_mech``         :math:`\eta_{mech}`      None        Mechanical Efficiency (accounts for losses in bearings...)
=========================== ======================== =========== ======================================================================

The table below shows important variables inherited from the pressure changer model.

=========================== ======================== =========== ==========================================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ==========================================================================================
``efficiency_isentropic``   :math:`\eta_{isen}`      time        Isentropic efficiency
``deltaP``                  :math:`\Delta P`         time        Pressure change (:math:`P_{out} - P_{in}`) [Pa]
``ratioP``                  :math:`P_{ratio}`        time        Ratio of discharge pressure to inlet pressure :math:`\left(\frac{P_{out}}{P_{in}}\right)`
=========================== ======================== =========== ==========================================================================================

:math:`\eta_{isentropic,t}` efficiency_isentropic Isentropic assumption only

Expressions
-----------

This model provides two expressions that are not available in the
pressure changer model.

=========================== ========================= =========== ======================================================================
Variable                    Symbol                    Index Sets  Doc
=========================== ========================= =========== ======================================================================
``power_thermo``            :math:`\dot{w}_{thermo}`  time        Turbine stage power output not including mechanical loss [W]
``power_shaft``             :math:`\dot{w}_{shaft}`   time        Turbine stage power output including mechanical loss (bearings...) [W]
=========================== ========================= =========== ======================================================================

Constraints
-----------

There are no additional constraints.

Initialization
--------------

This just calls the initialization routine from ``PressureChanger``, but it is wrapped in
a function to ensure the state after initialization is the same as before initialization.
The arguments to the initialization method are the same as PressureChanger.

TurbineStage Class
------------------

.. autoclass:: TurbineStage
  :members:

TurbineStageData Class
----------------------

.. autoclass:: TurbineStageData
  :members:
