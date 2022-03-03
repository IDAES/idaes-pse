Turbine (Inlet Stage)
=====================

.. index::
  pair: idaes.power_generation.unit_models.helm.turbine_inlet;HelmTurbineInletStage

.. module:: idaes.power_generation.unit_models.helm.turbine_inlet

This is a steam power generation turbine model for the inlet stage. It inherits
`HelmIsentropicTurbine
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Isentropic)>`.

The turbine inlet model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission for Use in a Process Simulation Software Environment." Journal of Engineering for Gas Turbines and Power. v136.


Example
-------

.. testcode::

    from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory, units
    from idaes.core import FlowsheetBlock
    from idaes.power_generation.unit_models.helm import HelmTurbineInletStage
    from idaes.generic_models.properties import iapws95

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = HelmTurbineInletStage(default={"property_package": m.fs.properties})
    hin = iapws95.htpx(T=880*units.K, P=2.4233e7*units.Pa)
    # set inlet
    m.fs.turb.inlet[:].enth_mol.fix(hin)
    m.fs.turb.inlet[:].flow_mol.fix(26000/4.0)
    m.fs.turb.inlet[:].pressure.fix(2.4233e7)
    m.fs.turb.eff_nozzle.fix(0.95)
    m.fs.turb.blade_reaction.fix(0.9)
    m.fs.turb.flow_coeff.fix(1.053/3600.0)
    m.fs.turb.blade_velocity.fix(110.0)
    m.fs.turb.efficiency_mech.fix(0.98)

    m.fs.turb.initialize()

Degrees of Freedom
------------------

Usually the inlet stream, or the inlet stream minus flow rate plus discharge
pressure are fixed. There are also a few variables which are turbine parameters
and are usually fixed, like flow coefficients.  See the variables section for
more information.

Model Structure
---------------

The turbine inlet stage model contains one :ref:`ControlVolume0DBlock block
<reference_guides/core/control_volume_0d:0D Control Volume Class>` called control\_volume and
inherits `HelmIsentropicTurbine
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Isentropic)>`.

Variables
---------
The variables below are defined in the HelmIsentropicTurbine model.

=========================== ======================== =========== ======================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ======================================================================
``blade_reaction``          :math:`R`                None        Blade reaction
``eff_nozzle``              :math:`\eta_{nozzle}`    None        Nozzle efficiency
``efficiency_mech``         :math:`\eta_{mech}`      None        Mechanical Efficiency (accounts for losses in bearings...)
``flow_coeff``              :math:`C_{flow}`         None        Turbine stage flow coefficient [kg*C^0.5/Pa/s]
``blade_velocity``          :math:`V_{rbl}`          None        Turbine blade velocity (should be constant while running) [m/s]
``delta_enth_isentropic``   :math:`\Delta h_{isen}`  time        Isentropic enthalpy change through stage [J/mol]
=========================== ======================== =========== ======================================================================

The table below shows important variables inherited from the pressure changer model.

=========================== ======================== =========== ==========================================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ==========================================================================================
``efficiency_isentropic``   :math:`\eta_{isen}`      time        Isentropic efficiency
``deltaP``                  :math:`\Delta P`         time        Pressure change (:math:`P_{out} - P_{in}`) [Pa]
``ratioP``                  :math:`P_{ratio}`        time        Ratio of discharge pressure to inlet pressure :math:`\left(\frac{P_{out}}{P_{in}}\right)`
=========================== ======================== =========== ==========================================================================================

Expressions
-----------

=========================== ========================= =========== ======================================================================
Variable                    Symbol                    Index Sets  Doc
=========================== ========================= =========== ======================================================================
``power_thermo``            :math:`\dot{w}_{thermo}`  time        Turbine stage power output not including mechanical loss [W]
``power_shaft``             :math:`\dot{w}_{shaft}`   time        Turbine stage power output including mechanical loss (bearings...) [W]
``steam_entering_velocity`` :math:`V_0`               time        Steam velocity entering stage [m/s]
=========================== ========================= =========== ======================================================================

The expression defined below provides a calculation for steam velocity entering
the stage, which is used in the efficiency calculation.

.. math::

  V_0 = 1.414\sqrt{\frac{-(1 - R)\Delta h_{isen}}{WT_{in}\eta_{nozzel}}}

Constraints
-----------

In addition to the constraints inherited from the `HelmTurbineStage
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Stage)>`,
this model contains two more constraints, one to estimate efficiency and
one pressure-flow relation. From the isentropic pressure changer model, these
constraints eliminate the need to specify efficiency and either inlet flow or
outlet pressure.

The isentropic efficiency is given by:

.. math::

  \eta_{isen} = 2 \frac{V_{rbl}}{V_0}\left[\left(\sqrt{1 - R} - \frac{V_{rbl}}{V_0}\right) +
   \sqrt{\left(\sqrt{1 - R} - \frac{V_{rbl}}{V_0}\right)^2 + R}\right]

The pressure-flow relation is given by:

.. math::

  \dot{m} = C_{flow}\frac{P_{in}}{\sqrt{T_{in}-273.15}}\sqrt{\frac{\gamma}{\gamma-1} \left[
    \left(\frac{P_{out}}{P_{in}}\right)^{\frac{2}{\gamma}} -
    \left(\frac{P_{out}}{P_{in}}\right)^{\frac{\gamma+1}{\gamma}} \right]}

Initialization
--------------

The initialization method for this model will save the current state of the model
before commencing initialization and reloads it afterwards.  The state of the model
will be the same after initialization, only the initial guesses for
unfixed variables will be changed and optionally a flow coefficient value can be
calculated.  To initialize this model, provide a starting value for the inlet port
variables. Then provide a guess for one of: discharge pressure, ``deltaP``, or
``ratioP``.  Since it can be hard to determine a proper flow coefficient, the
``calculate_cf`` argument of the ``initialize()`` method can be set to True, and
the deltaP guess will be used to calculate and set a corresponding flow coefficient.

The model should initialize readily, but it is possible to provide a flow
coefficient that is incompatible with the given flow rate resulting in an
infeasible problem.

TurbineInletStage Class
-----------------------

.. autoclass:: HelmTurbineInletStage
  :members:

TurbineInletStageData Class
---------------------------

.. autoclass:: HelmTurbineInletStageData
  :members:
