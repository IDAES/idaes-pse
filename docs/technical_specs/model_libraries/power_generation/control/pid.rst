Proportional-Integral-Derivative (PID) Controller
=================================================

.. index::
   pair: idaes.power_generation.control.pid_controller; Proportional-Integral-Derivative (PID) Controller

.. currentmodule:: idaes.power_generation.control.pid_controller

The IDAES framework contains a basic PID control implementation, which is described
in this section.

Introduction
------------
The PID controller model represents a PID controller in a feedback control loop of a process plant.  Depending on the user specified configuration, this model can be configurated as a proportional only (P), proportional and integral (PI), proportional and derivative (PD), or proportional, integral and derivative (PID) controller.  When declaring a controller model, the user needs to specify the model type through the configuration variable “type”, the process variable to be controlled y(t) through the configuration variable “pv”, and the manipulated variable u(t) through the configuration variable “mv”.  The “pv” and “mv” variables can be any of the time-indexed variable on a dynamic flowsheet.
The setpoint of the process variable r(t) is a variable declared inside the model named as “setpoint”, which is usually fixed or specified as a function of other process variables.
The variable or expression for the error e(t) is defined as the setpoint minus the process variable to be controlled.

.. math::

  e(t) = r(t) - y(t)

Since the proportional part exists in any type of controllers, the gain for the proportional part K_p (t) named as “gain_p” is always included as a variable in the model.  Note that the gain is declared as a time-indexed variable since it could be changed from time to time if a gain scheduling approach is applied.
If the model contains the integral part (for a PI or PID type controller), the user needs to specify the gain for the integral part K_i (t) named as “gain_i”. 
A variable named “integral_of_error” is also declared inside the model for the integral error e_i (t) defined as

.. math::

  e_i(t) = \int_0^t e(t')\text{d}t'


Since Pyomo.DAE does not provide a direct method to calculate the integral of a variable for the discretized equations, the PID controller model declares the integral term as a regular variable and a constraint that sets the derivative of the integral term with respect to time as the error term e(t)

.. math::

  \frac{\text{d}e_i(t)}{\text{d}t}= e(t)

If the model contains the derivative part (for a PD or PID type controller), the user needs to specify the gain for the derivative part K_i (t) named as “gain_d”.  A derivative variable named “derivative_of_error” is also declared inside the model for the derivative error e_d (t) defined as

.. math::

  e_i(t) = \frac{\text{d}e(t)}{\text{d}t}
  
Pyomo.DAE provides a method to calculate the derivative error using its “DerivativeVar” declaration for the discretized equations.

For a general PID controller, the deviation of the manipulated variable from its steady-state value u'(t) can be expressed as

.. math::

  u'(t) = K_p e(t) + K_i e_i (t) + K_d e_d (t)
  
Note that the terms for the integral and derivative part on the right side of the equation could be zero depending on the controller type.
The actual output for the manipulated variable u(t) is calculated as the sum of the deviation term u'(t) and the reference value also known as steady-state bias u_ref

.. math::

  u(t) = u'(t) + u_ref
  
Note that the steady-state bias is reference value that is not time-indexed in the current model.
To account for an actuator saturation condition, the calculated manipulated variable u(t) can optionally be clamped within a range between its lower and upper bounds.  For example, a control valve cannot close to less than 0% or open to more than 100%.  To declare this option, the configuration variable “bounded_output” is set to True.  With this option turned on, the user needs to set the value for two mutable parameters.  Parameter “mv_lb” is for the lower bound and parameter “mv_ub” is for the upper bound.  If they are not set by the user, the default values of 0.05 and 1 for the lower and upper bounds will be used as defaults, respectively.
Different clamping function has been tried to make the output of the manipulated variable u(t) within the lower and upper bounds.  Since Pyomo model requires all functions to be smooth and differentiable, current PID controller model uses a sigmoid function as the clamping function.  If the lower and upper bounds of u are u_l and u_u, respectively, the sigmoid function is defined as

.. math::

  f(u) = u_l + \frac{u_u-u_l}{1 + exp⁡[-\frac{4 (u - \frac{u_l + u_u}{2})} {u_u-u_l} ] }

where u is the calculated manipulated variable before the clamping function is applied and f(u) is the actual output value for the manipulated variable.  This function has the following properties:

.. math::

  f(-\infty) = u_l

.. math::

  f(\infty) = u_u

.. math::
  
  f(\frac{u_l + u_u}{2}) = \frac{u_l + u_u}{2}

.. math::

  \frac{\text{d}f}{\text{d}u} (\frac{u_l + u_u}{2}) = 1

Certainly, it is also smooth and differentiable.  One disadvantage of the function is that when all error terms are zero (e(t)=e_i (t)=e_d (t)=0), the final output for the manipulated variable is not steady-state bias unless the bias is the same as the average of the lower and upper bound values.

.. math::

 f(u) \neq U_ref
 
For a PI controller, if e(t)=0, to make f(u)=u_ref, the integral error e_i (t) should be set to a non-zero value, which can be calculated as

.. math::

  e_(i,o) (t) = \frac{1}{K_i} [\frac{u_l+u_u}{2} - U_ref - \frac{u_u-u_l}{4}  \ln⁡(\frac{u_u-u_l}{u_ref-u_l} - 1) ]

The model provides an expression named “integral_of_error_ref” for the above equation.  There is a similar expression named “integral_of_error_mv” to calculate the required integral error for a given f(u) value when e(t)=0.
While current clamping function is quite robust in solving the dynamic system with PID controllers, a better clamping function is still under development.  If a manipulated variable is unlikely to reach saturation, it is recommended to disable the clamping option.
It needs to be mentioned that the integral error in an PID controller could causes the wind-up issue if it gets larger and larger, especially for a slow process.  To reset the wind-up at certain time, the dynamic simulation can be performed in multiple discretized time periods and the integral error as a variable can be reset to zero or other small values after solving one time period and before solving a subsequent time period.
The PID controller model should be declared only for a dynamic flowsheet.  It needs to be mentioned that when a dynamic flowsheet is discretized in a time domain, the calculation for the manipulated variable at the initial time is skipped.  The user should provide the initial condition for the manipulated variable.  Also note that the current implementation of the PID controller model ignore the measurement delay for the process variable.  The current model simply provides continuous equations for the controller model, which is solved by the discretization through Pyomo.DAE.  This is different from the actual PID installed at a plant where the controller calculates the current maneuver based on the measured process variable from the previous time step.

