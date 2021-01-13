Valve
=====

.. index::
  pair: idaes.generic_models.unit_models.valve;Valve

.. module:: idaes.generic_models.unit_models.valve

This section describes the generic adiabatic valve model. By default the model
is based on molar flow, but the pressure-flow equation and the flow basis is
configurable. This model inherits the :ref:`PressureChanger model
<technical_specs/model_libraries/generic/unit_models/pressure_changer:Pressure Changer>`
with the adiabatic options. Beyond the base pressure changer model this provides a pressure
flow relation as a function of the valve opening fraction.

Example
-------

.. code-block:: python

  from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory

  from idaes.core import FlowsheetBlock
  from idaes.generic_models.unit_models import Valve
  from idaes.generic_models.properties import iapws95

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

======================= ======================== =========== =======================================================
Variable                Symbol                   Index Sets  Doc
======================= ======================== =========== =======================================================
``Cv``                  :math:`C_v`              None        Valve coefficient
``valve_opening``       :math:`x`                time        The fraction that the valve is open from 0 to 1
======================= ======================== =========== =======================================================

The ``Cv`` variable is highly recommended but can be omitted in custom pressure-flow relations.

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
<technical_specs/model_libraries/generic/unit_models/pressure_changer:Pressure Changer>`.

The default pressure-flow relation is given below where f is the molar flow.

.. math::

  F^2 = C_v^2\left(P_{in} - P_{out}\right)f(x)^2

Other functions can be specified via callback supplied to the unit configuration.
The callback allows both the form and flow basis of the pressure-flow equation to
be changed.  The callback should add a constraint that calculates flow as a function
of pressure.  It should also add an attribute ``pressure_flow_equation_scale`` which
is a function that provides the form that the flow apears in the pressure-flow
equation.  For example, by default this is ``lambda x : x ** 2``.  This is used to
calculate scale factors for the equation based on the flow variable scale. To
calculate scaling the flow basis is also needed, which is provided by assigning
a reference to the flow basis variable as ``flow_basis``. The reference should
be time-indexed.

Initialization
--------------

This just calls the initialization routine from PressureChanger.  Either and
outlet pressure value or deltaP can be specified to aid the initialization.

Valve Class
-----------

.. autoclass:: Valve
  :members:

ValveData Class
---------------------

.. autoclass:: ValveData
  :members:
