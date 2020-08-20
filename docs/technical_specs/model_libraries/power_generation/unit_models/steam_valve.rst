Valve
=====

.. index::
  pair: idaes.power_generation.unit_models.helm.valve_steam;HelmValve

.. module:: idaes.power_generation.unit_models.helm.valve_steam

This documentation covers the valve
(``idaes.power_generation.unit_models.helm.HelmValve``)
unit model intended to be used with Helmholtz Equation of State property
packages.

Examples
--------

The following example shows the basic usage of a valve model.

.. code-block:: python

  from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory

  from idaes.core import FlowsheetBlock
  from idaes.power_generation.unit_models.helm import HelmValve
  from idaes.generic_models.properties import iapws95
  from idaes.ui.report import degrees_of_freedom, active_equalities

  solver = SolverFactory('ipopt')
  solver.options = {'tol': 1e-6}

  m = ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = iapws95.Iapws95ParameterBlock()
  m.fs.valve = HelmValve(default={"property_package": m.fs.properties})

  hin = iapws95.htpx(T=880, P=2.4233e7)
  # set inlet
  m.fs.valve.inlet.enth_mol[0].fix(hin)
  m.fs.valve.inlet.flow_mol[0].fix(26000/4.0)
  m.fs.valve.inlet.pressure[0].fix(2.5e7)
  m.fs.valve.Cv.fix(0.01)
  m.fs.valve.valve_opening.fix(0.5)
  m.fs.valve.initialize(outlvl=1)


The next example is the same as the first, but the custom valve function is used
to add the valve function for an equal-percentage opening valve. Equal-percentage
is an available option, but this example is provided just to demonstrate how to
use a custom valve function.





Variables
---------

The variables in the table below are defined in the valve model, in addition to
the standard control volume variables.

=========================== ======================== =========== ======================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ======================================================================
``Cv``                      :math:`C_v`              None        Valve coefficient for liquid [mol/s/Pa^0.5] for vapor [mol/s/Pa]
``valve_opening``           :math:`x`                time        The fraction that the valve is open from 0 to 1
``alpha``                   :math:`\alpha`           None        Parameter only for equal percentage valve function
=========================== ======================== =========== ======================================================================

Expressions
-----------

The valve model provides a selection of expressions for the valve function.
Users may also provide their own custom valve function by providing a callback
function that takes the model as and argument and adds a ``valve_function``
Expression.  The custom callback can also add and additional model components
(e.g. variables for valve function parameters).

=========================== ========================= =========== ==================================================
Expression                  Symbol                    Index Sets  Doc
=========================== ========================= =========== ==================================================
``valve_function``          :math:`f(x)`              time        Transformation function for the valve opening.
=========================== ========================= =========== ==================================================

Predefined valve function expressions are given below, where :math:`x` is the
valve opening.  See the configurations options in section
:ref:`HelmValve<technical_specs/model_libraries/power_generation/unit_models/steam_valve:HelmValve Class>`
for information on selecting the valve function and providing custom functions.

**Predefined Valve Function Expressions**

============================= =======================================
Linear                        :math:`f(x) = x`
----------------------------- ---------------------------------------
Quick                         :math:`f(x) = \sqrt{x}`
----------------------------- ---------------------------------------
Equal Percentage              :math:`f(x) = \alpha^{(x - 1)}`
============================= =======================================


Constraints
-----------

The pressure flow relation is the only constraint in the valve model beyond the
standard mass and energy balance provided by the control volume block.

If the ``phase`` option is set to ``"Liq"`` the following equation describes the
pressure-flow relation.

.. math::

  F^2 = C_v^2\left(P_{in} - P_{out}\right)f(x)^2

If the ``phase`` option is set to ``"Vap"`` the following equation describes the
pressure-flow relation.

.. math::

  F^2 = C_v^2\left(P_{in}^2 - P_{out}^2\right)f(x)^2


Initialization
--------------

The valve model can be initialized in two ways.  The inlet conditions are fixed.
Then if ```control_volume.deltaP``` or the outlet pressure is fixed, inlet flow
is unfixed and the flow through the valve is calculated.  If the pressure drop
is not fixed the flow is used to calculate the pressure drop.  Mass and energy
balances are used to fill in values for the unfixed outlet variables.

The initialization just fills in values for unfixed variables, and does not
change the model specifications.


HelmValve Class
---------------

.. autoclass:: HelmValve
  :members:

HelmValveData Class
-------------------

.. autoclass:: HelmValveData
  :members:
