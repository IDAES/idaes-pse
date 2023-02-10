Proportional-Integral-Derivative (PID) Controller
=================================================

.. index::
   pair: idaes.models.control.controller; Proportional-Integral-Derivative (PID) Controller

.. currentmodule:: idaes.models.control.controller

The IDAES framework contains a basic PID control implementation, which is described
in this section.

Example
-------

The following code provides a full example for using the PIDController model.
The flowsheet model is a tank with a valve on the inlet and outlet. The dynamics
of the valves are assumed to be very fast and is not considered. The system has
steam flowing through it and the inlet valve is used to maintain the tank
pressure at a given setpoint. The final model solved demonstrates the
controller with a step change in the inlet pressure. Results are saved to a
file "pid_steam_tank_pressure.pdf."

.. testcode::

  import pyomo.environ as pyo
  from pyomo.network import Arc
  from idaes.core import FlowsheetBlock, MaterialBalanceType
  from idaes.models.unit_models import Heater, Valve
  from idaes.models.properties import iapws95
  from idaes.core.util.initialization import propagate_state
  from idaes.models.control.controller import (
      PIDController,
      ControllerType,
      ControllerMVBoundType
  )
  import idaes.core.util.scaling as iscale
  from idaes.core.solvers import get_solver
  from idaes.core.util.plot import plot_grid_dynamic

  def _valve_pressure_flow_cb(b):
      """Callback for the pressure flow relations for the valves"""
      umeta = b.config.property_package.get_metadata().get_derived_units
      b.Cv = pyo.Var(
          initialize=0.1,
          doc="Valve flow coefficent",
          units=umeta("amount") / umeta("time") / umeta("pressure"),
      )
      b.Cv.fix()
      b.flow_var = pyo.Reference(b.control_volume.properties_in[:].flow_mol)
      b.pressure_flow_equation_scale = lambda x: x ** 2
      @b.Constraint(b.flowsheet().time)
      def pressure_flow_equation(b2, t):
          Po = b2.control_volume.properties_out[t].pressure
          Pi = b2.control_volume.properties_in[t].pressure
          F = b2.control_volume.properties_in[t].flow_mol
          Cv = b2.Cv
          fun = b2.valve_function[t]
          return F ** 2 == Cv ** 2 * (Pi ** 2 - Po ** 2) * fun ** 2

  def _add_inlet_pressure_step(m, time=1, value=6.0e5):
      """Add an inlet pressure step change"""
      for t in m.fs.time:
          if t >= time:
              m.fs.valve_1.inlet.pressure[t].fix(value)

  def create_model(time_set=None, time_units=pyo.units.s, nfe=5, tee=False):
      """Build and initialize the flowsheet model

             valve_1   +----+
    steam ----|><|-->--| ta |    valve_2
                       | nk |-----|><|--->--- steam
                       +----+
      """
      m = pyo.ConcreteModel(name="Dynamic Steam Tank with PID Control")
      # Add flowsheet
      m.fs = FlowsheetBlock(dynamic=True, time_set=time_set, time_units=time_units)
      # Add water property parameter block
      m.fs.prop_water = iapws95.Iapws95ParameterBlock(
          phase_presentation=iapws95.PhaseType.LG
      )
      # Add valve 1
      m.fs.valve_1 = Valve(
          dynamic=False,
          has_holdup=False,
          pressure_flow_callback=_valve_pressure_flow_cb,
          material_balance_type=MaterialBalanceType.componentTotal,
          property_package=m.fs.prop_water,
      )
      # Add heater model to represent a tank (close to bare control volume model).
      m.fs.tank = Heater(
          has_holdup=True,
          material_balance_type=MaterialBalanceType.componentTotal,
          property_package=m.fs.prop_water,
      )
      # Add valve 2
      m.fs.valve_2 = Valve(
          dynamic=False,
          has_holdup=False,
          pressure_flow_callback=_valve_pressure_flow_cb,
          material_balance_type=MaterialBalanceType.componentTotal,
          property_package=m.fs.prop_water,
      )
      # Add a controller
      m.fs.ctrl = PIDController(
          process_var=m.fs.tank.control_volume.properties_out[:].pressure,
          manipulated_var=m.fs.valve_1.valve_opening,
          calculate_initial_integral=True,
          mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
          controller_type=ControllerType.PI,
      )
      # The control volume block doesn't assume equilibrium, so I'll make that
      # assumption here. I don't actually expect liquid to form but who knows?
      # The phase_fraction in the control volume is volumetric phase fraction.
      @m.fs.tank.Constraint(m.fs.time)
      def vol_frac_vap(b, t):
          return (
              b.control_volume.properties_out[t].phase_frac["Vap"]
              * b.control_volume.properties_out[t].dens_mol
              / b.control_volume.properties_out[t].dens_mol_phase["Vap"]
          ) == (b.control_volume.phase_fraction[t, "Vap"])

      # Connect the models
      m.fs.v1_to_tank = Arc(source=m.fs.valve_1.outlet, destination=m.fs.tank.inlet)
      m.fs.tank_to_v2 = Arc(source=m.fs.tank.outlet, destination=m.fs.valve_2.inlet)

      # Add the stream constraints
      pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)

      # Do DAE discretization
      pyo.TransformationFactory("dae.finite_difference").apply_to(
          m.fs,
          nfe=nfe,
          wrt=m.fs.time,
          scheme="BACKWARD"
      )

      # Fix the derivative variables to zero at time 0 (steady state assumption)
      m.fs.fix_initial_conditions()

      # Fix the input variables
      m.fs.valve_1.Cv.fix(0.001) # valve 1 flow coefficent
      m.fs.valve_2.Cv.fix(0.001) # valve 2 flow coefficent
      m.fs.ctrl.gain_p.fix(1e-6) # porportional gain
      m.fs.ctrl.gain_i.fix(1e-5) # integral error gain
      m.fs.ctrl.setpoint.fix(3e5) # setpoint is tank pressure of 300 kPa
      m.fs.ctrl.mv_ref.fix(0) # controller bias
      m.fs.ctrl.mv_lb = 0.0 # valve opening lower bound is 0
      m.fs.ctrl.mv_ub = 1.0 # value opening upper bound is 1
      m.fs.valve_1.inlet.enth_mol.fix(50000) # inlet fluid molar enthalpy
      m.fs.valve_1.inlet.pressure.fix(5e5) # inlet pressure of valve 1
      m.fs.valve_2.outlet.pressure.fix(101325) # valve 2 outlet pressure
      m.fs.valve_1.valve_opening.fix(1) # fix valve 1 opening to full open
      m.fs.valve_2.valve_opening.fix(1) # fix valve 2 opening to full open
      m.fs.tank.heat_duty.fix(0) # no heat transfer to tank
      m.fs.tank.control_volume.volume.fix(2.0) # 2 m^3 tank volume

      # Initialize the model
      solver = get_solver(options={"max_iter": 50})

      for t in m.fs.time:
          m.fs.valve_1.inlet.flow_mol[t] = 100  # initial guess on flow
      # simple initialize
      m.fs.valve_1.initialize()
      propagate_state(m.fs.v1_to_tank)
      m.fs.tank.initialize()
      propagate_state(m.fs.tank_to_v2)
      # Can't specify both flow and outlet pressure so free the outlet pressure
      # for initialization and refix it after.  Inlet flow gets fixed in init
      op = {}
      for t in m.fs.time:
          op[t] = pyo.value(m.fs.valve_2.outlet.pressure[t])
          m.fs.valve_2.outlet.pressure[t].unfix()
      m.fs.valve_2.initialize()
      for t in m.fs.time:
          m.fs.valve_2.outlet.pressure[t].fix(op[t])

      # Solve dynamic problem with deactivated controller
      m.fs.ctrl.deactivate()
      m.fs.valve_1.valve_opening.fix()
      solver.solve(m, tee=tee)

      # Solve dynamic problem with controller on
      m.fs.ctrl.activate()
      m.fs.valve_1.valve_opening.unfix()
      m.fs.valve_1.valve_opening[m.fs.time.first()].fix()
      solver.solve(m, tee=tee)

      # Return the model and solver
      return m, solver

  # Create a model for the 0 to 12 sec time period
  m, solver = create_model(time_set=[0, 12], nfe=30, tee=False)

  # Add a step change in inlet pressure at t=6 from 500 kPa to 550 kPa
  _add_inlet_pressure_step(m, time=6, value=5.5e5)

  # Solve with step change
  solver.solve(m, tee=False)

  # Plot some results and save the figure to a file.
  plot_grid_dynamic(
      x=m.fs.time,
      xlabel="time (s)",
      y=[
          m.fs.valve_1.valve_opening,
          m.fs.tank.control_volume.properties_out[:].pressure,
          m.fs.valve_1.control_volume.properties_in[:].pressure,
      ],
      ylabel=[
          "opening (fraction open)",
          "tank pressure (kPa)",
          "inlet pressure (kPa)",
      ],
      yunits=[
          None,
          pyo.units.kPa,
          pyo.units.kPa,
      ],
      cols=3,
      rows=1,
      to_file="pid_steam_tank_pressure.pdf"
  )

  assert abs(
    pyo.value(m.fs.valve_1.valve_opening[m.fs.time.last()]) - 0.61254) < 0.001



Class Documentation
-------------------

.. autoclass:: PIDController
   :members:

.. autoclass:: PIDControllerData
   :members:

.. _PIDVarsSection:

Variables, Parameters, and Expressions
--------------------------------------

.. table::
    :widths: 20 25 40

    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | Symbol                                              | Name in Model                   | Description                                 |
    +=====================================================+=================================+=============================================+
    | :math:`y_{sp}(t)`                                   | ``setpoint[t]``                 | Setpoint variable                           |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`y(t)`                                        | ``process_var[t]``              | Measured process variable (Reference)       |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`u(t)`                                        | ``manipulated_var[t]``          | Manipulated variable (Reference)            |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`K_p(t)`                                      | ``gain_p[t]``                   | Controller gain  (usually fixed)            |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`K_i(t)`                                      | ``gain_i[t]``                   | Integral gain (usually fixed)               |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`K_d(t)`                                      | ``gain_d[t]``                   | Derivative gain (usually fixed)             |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`K_b(t)`                                      | ``gain_b[t]``                   | Back-calculation gain (usually fixed)       |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`u_{ref}`                                     | ``mv_ref[t]``                   | Reference value of manipulated variable     |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`e(t)`                                        | ``error[t]``                    | Error expression or variable (setpoint - pv)|
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`e_d(t) \quad\text{or}\; -\frac{dy(t)}{dt}`   | ``derivative_term[t]``          | Derivative of either error or negative pv   |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`u_i(t)`                                      | ``mv_integral_component[t]``    | Portion of control from integral action     |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`\frac{de_i(t)}{dt}`                          | ``mv_integral_component_dot[t]``| Derivative of integral term w.r.t. time     |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`u_{ub}`                                      | ``mv_ub``                       | Upper limit of output parameter             |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`u_{lb}`                                      | ``mv_lb``                       | Lower limit of output parameter             |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`\epsilon`                                    | ``smooth_eps``                  | Smooth min/max smoothing parameter          |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`k_\text{log}`                                | ``logistic_bound_k``            | Logistic bound steepness parameter          |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+
    | :math:`k_\text{cond}`                               | ``conditional_integration_k``   | Conditional integration steepness parameter |
    +-----------------------------------------------------+---------------------------------+---------------------------------------------+

Formulation
-----------
Textbook PID
~~~~~~~~~~~~

The PID controller equations are given by the following equations (see the
:ref:`Variables, Parameters and Expressions <PIDVarsSection>` section for a
description of the variables).

The PID unbounded controller model is given by the following equation where the
integral and derivative terms are optional.

.. math::

  u(t) = K_p e(t) + K_i e_i(t) + K_d e_d(t) + u_{ref}(t)


The error term is given by the following equation.

.. math::

  e(t) = y_{sp}(t) - y(t)

The integral error is defined as follows, where the initial condition of the
integral error can be calculated from the manipulated variable initial condition
to ensure that the initial conditions are consistent and avoid an initial jump in
the manipulated variable value.

.. math::

  e_i(t) = e_i(t_0) + \int_{t_0}^t e(\tau) \text{d}\tau

The derivative of error is given below.

.. math::

  e_d(t) = \frac{de(t)}{dt}

Actual Implementation
~~~~~~~~~~~~~~~~~~~~~

There is often reason to change the integral gain :math:`K_i` (for example, using different values at different
operating points). However, if :math:`K_i` changes while :math:`e_i` remains unchanged, a large "bump" in :math:`u(t)`
results. Therefore, the product :math:`K_i e_i` is integrated instead.

.. math::

  u_i(t) = u_i(t_0) + \int_{t_0}^t K_i(\tau) e(\tau) d\tau

The integral error is formulated as a differential equation in the model.

.. math::

  \frac{du_i(t)}{dt} = K_i(t) e(t)

A naive implementation of derivative action also runs into problems. If the setpoint experiences a step change, the
derivative of error is undefined, which can result in numerical problems in simulations and a phenomenon called
"derivative kick" in practical implementations. A common way to address this problem is by taking the derivative
of :math:`-y(t)` rather than :math:`y_{sp}(t) - y(t)`. The derivative action is the same when the setpoint is fixed
but no derivative kick happens when the setpoint is changed. The naive derivative-in-error implementation can be
enabled by setting `derivative_in_error=True` in the config, but it is disabled by default.

Therefore, the final control action for a PID controller with `derivative_in_error=False` is given by: 

.. math::

  u(t) = K_p(t) e(t) + u_i(t) - K_d(t) \frac{\text{d}y(t)}{\text{d}t} + u_{ref}(t)



Output Limits
~~~~~~~~~~~~~

There are two available ways to bound the controller output.  The smooth max
method provides the closest approximation, but the logistic method provides a
method that is smoother and potentially more tractable.

Smooth Max
""""""""""

Smooth min and smooth max expressions are used to keep the controller output
between a minimum and maximum value. The exact form of is given by the following.

.. math::

  u_{bound}(t) = \min\left(\max\left(u_{unbound}(t), u_{lb} \right), u_{ub}\right)

To avoid a non-smooth formulation, the min and max functions are replaced by
smooth versions as follows.

.. math::

  \max(a, b) \approx \frac{1}{2}\left[a + b + \left((a-b)^2 + \epsilon^2\right)^{\frac{1}{2}}\right]

.. math::

  \min(a, b) \approx \frac{1}{2}\left[a + b - \left((a-b)^2 + \epsilon^2\right)^{\frac{1}{2}}\right]

Logistic Function
"""""""""""""""""

The logistic bounding function is more smooth and possibly more tractable, but is
less faithful to the unbound manipulated variable values and slower to arrive
at the bounds for reasonable values of :math:`k_\text{log}`.

.. math::

  u_{bound}(t) = u_{lb} + \frac{(u_{ub} - u_{lb})}{1 + \exp \left( -k_\text{log}\frac{u_{unbound}(t) - (u_{ub} + u_{lb})/2}{(u_{ub} - u_{lb})}\right)}


Anti-integral-windup
~~~~~~~~~~~~~~~~~~~~

When the controller output is clipped to a bound the integral term will continue
to grow.  This can cause a delay in the controller output moving away from its bound.
In order to mitigate this phenomenon, two methods of anti-windup have been implemented:
conditional integration and back calculation.

Conditional Integration
"""""""""""""""""""""""

Conditional integration, enabled by setting

.. code-block::

    antiwindup=ControllerAntiwindupType.CONDITIONAL_INTEGRATION

in the controller config, operates on simple logic:

.. math::

  \frac{du_i(t)}{dt} =
  \begin{cases}
        K_i(t) e(t) & \text{if} \quad u_{lb}\le u_\text{unbound}\le u_{ub} \\
        0 & \text{otherwise}
  \end{cases}

This logic is implemented by first writing this differential equation in terms of
Heaviside step functions :math:`H(x)`.

.. math::

    \frac{du_i(t)}{dt} = K_i(t) e(t) \left(H\left(\frac{u_\text{unbound} -u_{lb}}{u_{ub} - u_{lb}}\right)
    - H\left(\frac{u_\text{unbound} -u_{ub}}{u_{ub} - u_{lb}}\right) \right)

and then approximating the (discontinuous) step function using a steep logistic function

.. math::

    H(x) \approx \frac{1}{1+\exp(-2k_\text{cond}x)}

in which :math:`k_\text{cond}` is a parameter governing the steepness of the approximation---the larger it is, the
steeper the transitions between 0 and 1 are.

Back Calculation
""""""""""""""""

Conditional integration is simple to use, but the additional nonlinearity can cause problems for integrators. Back
calculation, enabled by setting

.. code-block::

    antiwindup=ControllerAntiwindupType.BACK_CALCULATION

in the `PIDController` configuration, does not add any additional nonlinearity, but utilizes the nonlinearity
inherent in enforcing controller bounds. The differential equation describing error accumulation is modified:

.. math::

      \frac{du_i(t)}{dt} = K_i(t) e(t) + K_b(t)(u_\text{bound}(t) - u_\text{unbound}(t))

in which :math:`K_b(t) \ge 0` is  a new gain term that must be chosen by the user. The more the unbounded control
input exceeds the bounded one, the more the integral term unwinds as a result.



To do
-----

Example with step-by-step initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dynamic models with a PID controller may be more difficult to initialize.  The
error and integral error variables may be especially difficult.  It is recommended
that these models be initialized with a time-stepping method. Examples of this
type of initialization should be added.
