Proportional-Integral-Derivative (PID) Controller
=================================================

.. index::
   pair: idaes.unit_models.heat_exchanger; Proportional-Integral-Derivative (PID) Controller

.. currentmodule:: idaes.dynamic.pid_controller

Class Documentation
===================

.. autoclass:: PIDBlock
   :members:

.. autoclass:: PIDBlockData
   :members:

Variables and Expressions
=========================

To setup a PID block, once it's created you generally need to set the setpoint, gain,
integral time, and derivative time.  The controller output limits should be set in the
PID configuration when the object is created.  The ``pv`` and ``output`` variables should be
model variables, expressions, or references pass the PIDBlock in the configuration
when the object is created.

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
===========

There are two forms of the PID controller equation.  The standard formulation
can result in equations with very large summations.  In the velocity form of the
equation the controller output can be calculated based only on the previous state.

The two forms of the equations are equivalent, but the choice of form will affect
robustness and solution time.  It is not necessarily clear that the velocity form
of the equation is always more numerically favorable, however it would usually be
the default choice. Both forms are provided in case the standard form works better
in some situations.

Standard Formulation
--------------------

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
--------------------

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
------------

In both the proportional and integral terms error can be replaced with the negative
measured process variable yielding equivalent results. This substitution is provided
by the PID class and is done by default.

Output Limits
-------------

Smooth min and smooth max expressions are used to keep the controller output between
a minimum and maximum value.
