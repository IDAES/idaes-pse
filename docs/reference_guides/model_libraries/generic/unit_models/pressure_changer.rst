Pressure Changer
================

The IDAES Pressure Changer model represents a unit operation with a single stream of material which undergoes a change in pressure due to the application of a work. The Pressure Changer model contains support for a number of different thermodynamic assumptions regarding the working fluid.

Degrees of Freedom
------------------

Pressure Changer units generally have one or more degrees of freedom, depending on the thermodynamic assumption used.

Typical fixed variables are:

* outlet pressure, :math:`P_{ratio}` or :math:`\Delta P`,
* unit efficiency (isentropic or pump assumption).

Model Structure
---------------

The core Pressure Changer unit model consists of a single ``ControlVolume0D`` (named ``control_volume``) with one Inlet Port (named ``inlet``) and one Outlet Port (named ``outlet``). Additionally, if an isentropic pressure changer is used, the unit model contains an additional ``StateBlock`` named ``properties_isentropic`` at the unit model level.

Variables
---------

Pressure Changers contain the following Variables (not including those contained within the control volume Block):

=========================== ===================== ===================================================================================
Variable                    Name                  Notes
=========================== ===================== ===================================================================================
:math:`P_{ratio}`           ratioP
:math:`V_t`                 volume                Only if has_rate_reactions = True, reference to control_volume.rate_reaction_extent
:math:`W_{mechanical,t}`    work_mechanical       Reference to control_volume.work
:math:`W_{fluid,t}`         work_fluid            Pump assumption only
:math:`\eta_{pump,t}`       efficiency_pump       Pump assumption only
:math:`W_{isentropic,t}`    work_isentropic       Isentropic assumption only
:math:`\eta_{isentropic,t}` efficiency_isentropic Isentropic assumption only
=========================== ===================== ===================================================================================

Isentropic Pressure Changers also have an additional Property Block named `properties_isentropic` (attached to the Unit Model).

Constraints
-----------

In addition to the Constraints written by the Control Volume block, Pressure Changer writes additional Constraints which depend on the thermodynamic assumption chosen. All Pressure Changers add the following Constraint to calculate the pressure ratio:

.. math:: P_{ratio,t} \times P_{in,t} = P_{out,t}

Isothermal Assumption
~~~~~~~~~~~~~~~~~~~~~

The isothermal assumption writes one additional Constraint:

.. math:: T_{out} = T_{in}

Adiabatic Assumption
~~~~~~~~~~~~~~~~~~~~

The isothermal assumption writes one additional Constraint:

.. math:: H_{out} = H_{in}

Isentropic Assumption
~~~~~~~~~~~~~~~~~~~~~

The isentropic assumption creates an additional set of Property Blocks (indexed by time) for the isentropic fluid calculations (named properties_isentropic). This requires a set of balance equations relating the inlet state to the isentropic conditions, which are shown below:

.. math:: F_{in,t,p,j} = F_{isentropic,t,p,j}
.. math:: s_{in,t} = s_{isentropic,t}
.. math:: P_{in,t} \times P_{ratio,t} = P_{isentropic,t}

where :math:`F_{t,p,j}` is the flow of component :math:`j` in phase :math:`p` at time :math:`t` and :math:`s` is the specific entropy of the fluid at time :math:`t`.

Next, the isentropic work is calculated as follows:

.. math:: W_{isentropic,t} = \sum_p{H_{isentropic,t,p}} - \sum_p{H_{in,t,p}}

where :math:`H_{t,p}` is the total energy flow of phase :math:`p` at time :math:`t`. Finally, a constraint which relates the fluid work to the actual mechanical work via an efficiency term :math:`\eta`.

If compressor is True, :math:`W_{isentropic,t} = W_{mechanical,t} \times \eta_t`

If compressor is False, :math:`W_{isentropic,t} \times \eta_t = W_{mechanical,t}`

Pump (Incompressible Fluid) Assumption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The incompressible fluid assumption writes two additional constraints. Firstly, a Constraint is written which relates fluid work to the pressure change of the fluid.

.. math:: W_{fluid,t} = (P_{out,t}-P_{in,t})\times F_{vol,t}

where :math:`F_{vol,t}` is the total volumetric flowrate of material at time :math:`t` (from the outlet Property Block). Secondly, a constraint which relates the fluid work to the actual mechanical work via an efficiency term :math:`\eta`.

If compressor is True, :math:`W_{fluid,t} = W_{mechanical,t} \times \eta_t`

If compressor is False, :math:`W_{fluid,t} \times \eta_t = W_{mechanical,t}`

Performance Curves
------------------
Isentropic pressure changers support optional performance curve constraints.
The exact form of these constraints is left to the user, but generally the
constraints take the form of one or two equations which provide a correlation
between head, efficiency, or pressure ratio and mass or volumetric flow.
Additional variables such as compressor or turbine speed can be added if needed.

Performance curves should be added to the ``performance_curve`` sub-block rather
than adding them elsewhere because it allows them to be integrated into the
unit model initialization. It also provides standardization for users and
provides a convenient way to turn the performance equations on and off by
activating and deactivating the block.

Usually there are one or two performance curve constraints. Either directly or
indirectly, these curves specify an efficiency and pressure drop, so in adding
performance curves the efficiency and/or pressure drop should be freed as
appropriate.

Performance equations generally are in a simple form (e.g. efficiency = f(flow)),
where no special initialization is needed. Performance curves also are specific
to a particular property package selection and pressure changer, which allows
the performance curve equations to be written in a well-scaled way since units
of measure and magnitudes are known.

To add performance curves to an isentropic pressure changer, simply supply the
``"support_isentropic_performance_curves": True`` options in the pressure
changer config dict.  This will create a ``performance_curve`` sub-block of
the pressure changer model. By default this block will have the expressions
``head`` and ``heat_isentropic`` for convenience, as these quantities often
appear in performance curves.

Two examples are provided below that demonstrate two ways to add performance
curves. The first does not use a callback the second does.

.. testcode::

  from pyomo.environ import ConcreteModel, SolverFactory, units, value
  from idaes.core import FlowsheetBlock
  from idaes.models.unit_models.pressure_changer import Turbine
  from idaes.models.properties import iapws95
  import pytest

  solver = SolverFactory('ipopt')
  m = ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = iapws95.Iapws95ParameterBlock()
  m.fs.unit = Turbine(default={
      "property_package": m.fs.properties,
      "support_isentropic_performance_curves":True})

  # Add performance curves
  @m.fs.unit.performance_curve.Constraint(m.fs.config.time)
  def pc_isen_eff_eqn(b, t):
      # main pressure changer block parent of performance_curve
      prnt = b.parent_block()
      return prnt.efficiency_isentropic[t] == 0.9
  @m.fs.unit.performance_curve.Constraint(m.fs.config.time)
  def pc_isen_head_eqn(b, t):
      # divide both sides by 1000 for scaling
      return b.head_isentropic[t]/1000 == -75530.8/1000*units.J/units.kg

  # set inputs
  m.fs.unit.inlet.flow_mol[0].fix(1000)  # mol/s
  Tin = 500  # K
  Pin = 1000000  # Pa
  Pout = 700000  # Pa
  hin = iapws95.htpx(Tin*units.K, Pin*units.Pa)
  m.fs.unit.inlet.enth_mol[0].fix(hin)
  m.fs.unit.inlet.pressure[0].fix(Pin)

  m.fs.unit.initialize()
  solver.solve(m, tee=False)

  assert value(m.fs.unit.efficiency_isentropic[0]) == pytest.approx(0.9, rel=1e-3)
  assert value(m.fs.unit.deltaP[0]) == pytest.approx(-3e5, rel=1e-3)

The next example shows how to use a callback to add performance curves.

.. testcode::

  from pyomo.environ import ConcreteModel, SolverFactory, units, value
  from idaes.core import FlowsheetBlock
  from idaes.models.unit_models.pressure_changer import Turbine
  from idaes.models.properties import iapws95
  import pytest

  solver = SolverFactory('ipopt')
  m = ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = iapws95.Iapws95ParameterBlock()

  def perf_callback(blk):
      # This callback adds constraints to the performance_cruve block. blk is the
      # performance_curve block, but we also want to use quantities from the main
      # pressure changer model, which is the parent block.
      prnt = blk.parent_block()
      # this is the pressure changer model block
      @blk.Constraint(m.fs.config.time)
      def pc_isen_eff_eqn(b, t):
          return prnt.efficiency_isentropic[t] == 0.9
      @blk.Constraint(m.fs.config.time)
      def pc_isen_head_eqn(b, t):
          return b.head_isentropic[t]/1000 == -75530.8/1000*units.J/units.kg

  m.fs.unit = Turbine(default={
      "property_package": m.fs.properties,
      "support_isentropic_performance_curves":True,
      "isentropic_performance_curves": {"build_callback": perf_callback}})

  # set inputs
  m.fs.unit.inlet.flow_mol[0].fix(1000)  # mol/s
  Tin = 500  # K
  Pin = 1000000  # Pa
  Pout = 700000  # Pa
  hin = iapws95.htpx(Tin*units.K, Pin*units.Pa)
  m.fs.unit.inlet.enth_mol[0].fix(hin)
  m.fs.unit.inlet.pressure[0].fix(Pin)

  m.fs.unit.initialize()
  solver.solve(m, tee=False)

  assert value(m.fs.unit.efficiency_isentropic[0]) == pytest.approx(0.9, rel=1e-3)
  assert value(m.fs.unit.deltaP[0]) == pytest.approx(-3e5, rel=1e-3)


PressureChanger Class
----------------------

.. module:: idaes.models.unit_models.pressure_changer

.. autoclass:: PressureChanger
  :members:

PressureChangerData Class
--------------------------

.. autoclass:: PressureChangerData
  :members:
