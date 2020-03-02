Proportional-Integral-Derivative (PID) Controller
=================================================

.. index::
   pair: idaes.unit_models.heat_exchanger; Proportional-Integral-Derivative (PID) Controller

.. currentmodule:: idaes.dynamic.pid_controller

The IDAES framework contains a basic PID control implementation, which is described
in this section.

Example
-------

The following code demonstrated the creation of a PIDBlock, but for simplicity, it does
not create a dynamic process model.  A full example of a dynamic process with PID control
is being prepared for the IDAES examples repository and will be referenced here once completed.

The valve opening is the controlled output variable and pressure "1" is the measured variable.
The controller output for the valve opening is restricted to be between 0 and 1. The measured
and output variables should be indexed only by time.  Fortunately there is no need to create
new variables if the variables are in a property block or not indexed only by time. Pyomo's
Reference objects can be use to create references to existing variables with the proper
indexing as shown in the example.

The ``calculate_initial_integral`` option calculates the integral error in the first time
step to match the initial controller output.  This keeps the controller output from
immediately jumping to a new value.  Unless the initial integral error is known, this option
should usually be True.

The controller should be added after the DAE expansion is done. There are several variables
in the controller that are usually meant to be fixed; as shown in the example, they are
``gain``, ``time_i``, ``time_d``, and ``setpoint``.  For more information about the
variables, expressions, and parameters in the PIDBlock, model see :ref:`PIDVarsSection`.

.. testcode::

  from idaes.dynamic import PIDBlock, PIDForm
  from idaes.core import FlowsheetBlock
  import pyomo.environ as pyo

  m = pyo.ConcreteModel(name="PID Example")
  m.fs = FlowsheetBlock(default={"dynamic":True, "time_set":[0,10]})

  m.fs.valve_opening = pyo.Var(m.fs.time, doc="Valve opening")
  m.fs.pressure = pyo.Var(m.fs.time, [1,2], doc="Pressure in unit 1 and 2")

  pyo.TransformationFactory('dae.finite_difference').apply_to(
      m.fs,
      nfe=10,
      wrt=m.fs.time,
      scheme='BACKWARD',
  )

  m.fs.measured_variable = pyo.Reference(m.fs.pressure[:,1])

  m.fs.ctrl = PIDBlock(
      default={
          "pv":m.fs.measured_variable,
          "output":m.fs.valve_opening,
          "upper":1.0,
          "lower":0.0,
          "calculate_initial_integral":True,
          "pid_form":PIDForm.velocity,
      }
  )

  m.fs.ctrl.gain.fix(1e-6)
  m.fs.ctrl.time_i.fix(0.1)
  m.fs.ctrl.time_d.fix(0.1)
  m.fs.ctrl.setpoint.fix(3e5)

Controller Windup
-----------------

The current PID controller model has no integral windup prevention. This will be added
to the model in the near future.

Class Documentation
-------------------

.. autoclass:: PIDBlock
   :members:

.. autoclass:: PIDBlockData
   :members:

.. _PIDVarsSection:

Variables and Expressions
-------------------------

+-------------------------+--------------------+-------------------------------------------+
| Symbol                  | Name in Model      | Description                               |
+=========================+====================+===========================================+
| :math:`v_{sp}(t)`       | ``setpoint[t]``    | Setpoint variable (usually fixed)         |
+-------------------------+--------------------+-------------------------------------------+
| :math:`v_{mv}(t)`       | ``pv[t]``          | Measured process variable reference       |
+-------------------------+--------------------+-------------------------------------------+
| :math:`u(t)`            | ``output[t]``      | Controller output variable reference      |
+-------------------------+--------------------+-------------------------------------------+
| :math:`K_p(t)`          | ``gain[t]``        | Controller gain  (usually fixed)          |
+-------------------------+--------------------+-------------------------------------------+
| :math:`T_i(t)`          | ``time_i[t]``      | Integral time (usually fixed)             |
+-------------------------+--------------------+-------------------------------------------+
| :math:`T_d(t)`          | ``time_d[t]``      | Derivative time (usually fixed)           |
+-------------------------+--------------------+-------------------------------------------+
| :math:`e(t)`            | ``err[t]``         | Error expression (setpoint - pv)          |
+-------------------------+--------------------+-------------------------------------------+
| --                      | ``err_d[t]``       | Derivative error expression               |
+-------------------------+--------------------+-------------------------------------------+
| --                      | ``err_i[t]``       | Integral error expression (standard form) |
+-------------------------+--------------------+-------------------------------------------+
| --                      | ``err_d0``         | Initial derivative error value (fixed)    |
+-------------------------+--------------------+-------------------------------------------+
| :math:`e_{integral}(0)` | ``err_i0``         | Initial integral error value (fixed)      |
+-------------------------+--------------------+-------------------------------------------+
| --                      | ``err_i_end``      | Last initial integral error expression    |
+-------------------------+--------------------+-------------------------------------------+
| --                      | ``limits["h"]``    | Upper limit of output parameter           |
+-------------------------+--------------------+-------------------------------------------+
| --                      | ``limits["l"]``    | Lower limit of output parameter           |
+-------------------------+--------------------+-------------------------------------------+
| --                      | ``smooth_eps``     | Smooth min/max parameter                  |
+-------------------------+--------------------+-------------------------------------------+

Formulation
-----------

There are two forms of the PID controller equation.  The standard formulation
can result in equations with very large summations.  In the velocity form of the
equation the controller output can be calculated based only on the previous state.

The two forms of the equations are equivalent, but the choice of form will affect
robustness and solution time.  It is not necessarily clear that the velocity form
of the equation is always more numerically favorable, however it would usually be
the default choice. Both forms are provided in case the standard form works better
in some situations.

Standard Formulation
~~~~~~~~~~~~~~~~~~~~

The PID controller equations are given by the following equations

.. math::

  e(t) = v_{sp}(t) - v_{mv}(t)

.. math::

  u(t) = K_p \left[ e(t) + \frac{1}{T_i} \int_0^t e(s) \text{d}s + T_d \frac{\text{d}e(t)}{\text{d}t} \right]

The PID equation now must be discretized.

.. math::

  u(t_i) = K_p \left[
    e(t_i) +
    \frac{e_{integral}(0)}{T_i} + \frac{1}{T_i} \sum_{j=0}^{i-1} \Delta t_j \frac{e(t_j) + e(t_{j+1})}{2} +
    T_d \frac{e(t_i) - e(t_{i-1})}{\Delta t_{i-1}} \right]

Velocity Formulation
~~~~~~~~~~~~~~~~~~~~

The velocity formulation of the PID equation may also be useful.  The way the
equations are written in the PID model, the integral term is a summation expression
and as time increases the integral term will build up an increasing number of terms
potentially becoming very large.  This potentially has two affects, increasing
round off error and computation time.  The velocity formulation allows the controller
output to be calculated based on the previous output.

Frist the usual PID controller equation can be rearranged to solve for the integral
error.

.. math::

  \frac{1}{T_i} \int_0^t e(s) \text{d}s = \frac{u(t)}{K_p} - e(t) - T_d \frac{\text{d}e(t)}{\text{d}t}

The PID equation for some time (:math:`t + \Delta t`) is

.. math::

  u(t + \Delta t) = K_p \left[
    e(t + \Delta t) +
    \frac{1}{T_i} \int_0^{t+\Delta t} e(s) \text{d}s +
    T_d \frac{\text{d}e(t+\Delta t)}{\text{d}t}
  \right]

.. math::

  u(t + \Delta t) = K_p \left[
    e(t + \Delta t) +
    \frac{1}{T_i} \int_t^{t+\Delta t} e(s) \text{d}s +
    \frac{1}{T_i} \int_0^{t} e(s) \text{d}s +
    T_d \frac{\text{d}e(t+\Delta t)}{\text{d}t}
  \right]

.. math::

  u(t + \Delta t) = u(t) + K_p \left[
    e(t + \Delta t) - e(t) +
    \frac{1}{T_i} \int_t^{t+\Delta t} e(s) \text{d}s +
    T_d \left( \frac{\text{d}e(t+\Delta t)}{\text{d}t} - \frac{\text{d}e(t)}{\text{d}t}\right)
  \right]

Now we can discretize the equation using the trapezoid rule for the integral.

.. math::

  u(t + \Delta t) = u(t) + K_p \left[
    e(t + \Delta t) - e(t) +
    \frac{\Delta t}{T_i} \left(\frac{e(t+\Delta t) + e(t)}{2} \right) +
    T_d \left( \frac{\text{d}e(t+\Delta t)}{\text{d}t} - \frac{\text{d}e(t)}{\text{d}t}\right)
  \right]

Since the derivative error term will require the error at the previous time step
to calculate, this form will still result in a large summation being formed since
in the model there is no derivative error variable.  To avoid this problem, the
derivative error term can equivalently be replaced with the derivative of the
negative measured process variable.

.. math::

  u(t + \Delta t) = u(t) + K_p \left[
    e(t + \Delta t) - e(t) +
    \frac{\Delta t}{T_i} \left(\frac{e(t+\Delta t) + e(t)}{2} \right)+
    T_d \left( \frac{\text{d}v_{mv}(t+\Delta t)}{\text{d}t} - \frac{\text{d}v_{mv}(t)}{\text{d}t}\right)
  \right]

Now the velocity form of the PID controller equation can be calculated at each time
point from just the state at the previous time point.

Substitution
~~~~~~~~~~~~

In both the proportional and integral terms, error can be replaced with the negative
measured process variable yielding equivalent results. This substitution is provided
by the PID class and is done by default.

Output Limits
~~~~~~~~~~~~~

Smooth min and smooth max expressions are used to keep the controller output between
a minimum and maximum value.
