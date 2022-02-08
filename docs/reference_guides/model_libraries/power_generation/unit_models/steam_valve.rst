HelmValve
=========

.. index::
  pair: idaes.power_generation.unit_models.helm.valve_steam;HelmValve

.. module:: idaes.power_generation.unit_models.helm.valve_steam

This is a steam power generation turbine model for the stages between the inlet
and outlet.

Example
-------

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

Variables
---------

This model adds a variable to account for mechanical efficiency to the base PressureChanger
model.

=========================== ======================== =========== ======================================================================
Variable                    Symbol                   Index Sets  Doc
=========================== ======================== =========== ======================================================================
``Cv``                      :math:`C_v`              None        Valve coefficient for liquid [mol/s/Pa^0.5] for vapor [mol/s/Pa]
``valve_opening``           :math:`x`                time        The fraction that the valve is open from 0 to 1
=========================== ======================== =========== ======================================================================

Expressions
-----------

Currently this model provides two additional expressions, with are not available
in the pressure changer model.

=========================== ========================= =========== ===========================================================================
Expression                  Symbol                    Index Sets  Doc
=========================== ========================= =========== ===========================================================================
``valve_function``          :math:`f(x)`              time        This is a valve function that describes how the fraction open affects flow.
=========================== ========================= =========== ===========================================================================

Constraints
-----------

The pressure flow relation is added to the inherited constraints from the :ref:`PressureChanger model
<reference_guides/model_libraries/generic/unit_models/pressure_changer:Pressure Changer>`.

If the ``phase`` option is set to ``"Liq"`` the following equation describes the pressure-flow relation.

.. math::

  \frac{1}{s_f^2}F^2 = \frac{1}{s_f^2}C_v^2\left(P_{in} - P_{out}\right)f(x)^2

If the ``phase`` option is set to ``"Vap"`` the following equation describes the pressure-flow relation.

.. math::

  \frac{1}{s_f^2}F^2 = \frac{1}{s_f^2}C_v^2\left(P_{in}^2 - P_{out}^2\right)f(x)^2


Initialization
--------------

This just calls the initialization routine from PressureChanger, but it is wrapped in
a function to ensure the state after initialization is the same as before initialization.
The arguments to the initialization method are the same as PressureChanger.

HelmValve Class
----------------

.. autoclass:: HelmValve
  :members:

HelmValveData Class
---------------------

.. autoclass:: HelmValveData
  :members:
