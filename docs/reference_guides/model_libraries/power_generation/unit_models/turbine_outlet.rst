Turbine (Outlet Stage)
======================

.. index::
  pair: idaes.power_generation.unit_models.helm.turbine_outlet;HelmTurbineOutletStage

.. module:: idaes.power_generation.unit_models.helm.turbine_outlet

This is a steam power generation turbine model for the outlet stage. The turbine outlet model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission for Use in a Process Simulation Software Environment." Journal of Engineering for Gas Turbines and Power. v136.


Example
-------

.. testcode::

    from pyomo.environ import ConcreteModel, SolverFactory
    from idaes.core import FlowsheetBlock
    from idaes.power_generation.unit_models.helm import HelmTurbineOutletStage
    from idaes.generic_models.properties import iapws95

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = HelmTurbineOutletStage(default={"property_package": m.fs.properties})
    # set inlet
    m.fs.turb.inlet[:].enth_mol.fix(47115)
    m.fs.turb.inlet[:].flow_mol.fix(15000)
    m.fs.turb.inlet[:].pressure.fix(8e4)

    m.fs.turb.initialize()

Degrees of Freedom
------------------

Usually the inlet stream, or the inlet stream minus flow rate plus discharge
pressure are fixed. There are also a few variables which are turbine parameters
and are usually fixed.  See the variables section for more information.

Model Structure
---------------

The turbine outlet stage model contains one :ref:`ControlVolume0DBlock block
<reference_guides/core/control_volume_0d:0D Control Volume Class>` called control\_volume and
inherits the `HelmIsentropicTurbine
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Isentropic)>`.

Variables
---------
The variables below are defined int the TurbineInletStage model. Additional variables
are in inherited from the :`HelmIsentropicTurbine
<reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Isentropic)>` model.

=========================== ======================== =========== ======================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ======================================================================
``eff_dry``                 :math:`\eta_{dry}`       None        Turbine efficiency when no liquid is present.
``efficiency_mech``         :math:`\eta_{mech}`      None        Mechanical Efficiency (accounts for losses in bearings...)
``flow_coeff``              :math:`C_{flow}`         None        Turbine stage flow coefficient [kg*C^0.5/Pa/s]
``design_exhaust_flow_vol`` :math:`V_{des,exhaust}`  None        Design volumetric flow out of stage [m^3/s]
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
``tel``                     :math:`\text{TEL}`        time        Total exhaust loss [J/mol]
=========================== ========================= =========== ======================================================================

The expression defined below provides a total exhaust loss.

.. math::

  \text{TEL} = 1\times 10^6*\left(-0.0035f^5 + 0.022f^4 - 0.0542f^3 + 0.0638f^2 - 0.0328f + 0.0064\right)

Where :math:`f` is the total volumetric flow of the exhaust divided by the design flow.

Constraints
-----------

In addition to the constraints inherited from the :ref:`PressureChanger model
<reference_guides/model_libraries/generic/unit_models/pressure_changer:Pressure Changer>` with the isentropic options, this
model contains two more constraints, one to estimate efficiency and one pressure-flow
relation.  From the isentropic pressure changer model, these constraints eliminate the
need to specify efficiency and either inlet flow or outlet pressure.

The isentropic efficiency is given by:

.. math::

  \eta_{isen} = \eta_{dry}x\left(1 - 0.65(1 - x)\right)*\left(1 + \frac{\text{TEL}}{\Delta h_{isen}}\right)

Where :math:`x` is the steam quality (vapor fraction).


The pressure-flow relation is given by the Stodola Equation:

.. math::

  \dot{m}\sqrt{Tin - 273.15} = C_{flow}P_{in}\sqrt{1 - Pr^2}

Initialization
--------------

The initialization method for this model will save the current state of the model
before commencing initialization and reloads it afterwards.  The state of the model
will be the same after initialization, only the initial guesses for
unfixed variables will be changed except for optional calculation of the flow coefficient.
To initialize this model, provide a starting value for the inlet port variables.  Then
provide a guess for one of: discharge pressure, ``deltaP``, or ``ratioP``.  Since a
good flow coefficient can be difficult to determine, the ``calculate_cf`` option will
calculate and set a flow coefficient based on the specified inlet flow and ``deltaP``.

The model should initialize readily, but it is possible to provide a flow
coefficient that is incompatible with the given flow rate resulting in an
infeasible problem.

TurbineOutletStage Class
------------------------

.. autoclass:: HelmTurbineOutletStage
  :members:

TurbineOutletStageData Class
----------------------------

.. autoclass:: HelmTurbineOutletStageData
  :members:
