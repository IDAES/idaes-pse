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
  import idaes.core.util.scaling as iscale

  m = ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = iapws95.Iapws95ParameterBlock()
  m.fs.valve = Valve(default={"property_package": m.fs.properties})
  fin = 900 # mol/s
  pin = 200000 # Pa
  pout = 100000 # Pa
  tin = 300 # K
  hin = iapws95.htpx(T=tin*units.K, P=pin*units.Pa) # J/mol
  # Calculate the flow coefficient to give 1000 mol/s flow with given P
  cv = 1000/math.sqrt(pin - pout)/0.5
  # set inlet
  m.fs.valve.inlet.enth_mol[0].fix(hin)
  m.fs.valve.inlet.flow_mol[0].fix(fin)
  m.fs.valve.inlet.flow_mol[0].unfix()
  m.fs.valve.inlet.pressure[0].fix(pin)
  m.fs.valve.outlet.pressure[0].fix(pout)
  m.fs.valve.Cv.fix(cv)
  m.fs.valve.valve_opening.fix(0.5)
  iscale.calculate_scaling_factors(m)
  m.fs.valve.initialize(outlvl=1)

  solver = pyo.SolverFactory("ipopt")
  solver.options = {"nlp_scaling_method": "user-scaling"}
  solver(m, tee=True)


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

=========================== ========================= =========== ===========================================================================
Expression                  Symbol                    Index Sets  Doc
=========================== ========================= =========== ===========================================================================
``valve_function``          :math:`f(x)`              time        This is a valve function that describes how the fraction open affects flow.
=========================== ========================= =========== ===========================================================================

Built-in Valve Functions
~~~~~~~~~~~~~~~~~~~~~~~~

Standard valve functions can be specified by providing a ``ValveFunctionType``
enumerated type to the ``valve_function_callback`` argument.  Standard functions
are given below.

``ValveFunctionType.linear``

.. math::

  f(x) = x

``ValveFunctionType.quick_opening``

.. math::

  f(x) = \sqrt{x}

``ValveFunctionType.equal_percentage``

.. math::

  f(x) = \alpha^{x - 1}

For the equal-percentage valve function an additional variable ``alpha`` is defined
which by default is fixed and set to 100.

Custom Valve Functions
~~~~~~~~~~~~~~~~~~~~~~

In general, the valve opening should be restricted to range from 0 to 1.  The valve
function should be a named expression attached to the valve model called
``valve_function`` which takes the valve opening and computes a value that goes
from approximately zero when valve opening is 0 to 1 when the valve opening is one.
The valve function can have parameters as needed, so custom valve functions are
defined using a callback function.

The callback function should take an object of the ``Valve`` class as an argument
and add the ``valve_function`` named expression.  Any additional parameters can also
be added.  The standard equal-percentage valve function is provided below as an
example.  The callback can be provided for the ``valve_function_callback``
configuration option.

.. code-block:: python

  def equal_percentage_cb(b):
      """
      Equal percentage valve function callback.
      """
      # Parameters can be defined as Var or Param.  If Var is used the parameter
      # can be included in a parameter estimation problem.
      b.alpha = pyo.Var(initialize=100, doc="Valve function parameter")
      b.alpha.fix()
      @b.Expression(b.flowsheet().config.time)
      def valve_function(b2, t):
          return b2.alpha ** (b2.valve_opening[t] - 1)

Constraints
-----------

The pressure flow relation is added to the inherited constraints from the :ref:`PressureChanger model
<technical_specs/model_libraries/generic/unit_models/pressure_changer:Pressure Changer>`.

The default pressure-flow relation is given below where :math:`F` is the molar flow. The default
valve function assumes an incompressible fluid of constant density.  In this case the
fluid specific gravity is included in the flow coefficient.  For rigorous modeling of valves
with gases, it is recommended that a custom pressure-flow equation be specified.

.. math::

  F^2 = C_v^2\left(P_{in} - P_{out}\right)f(x)^2

Custom Pressure Flow Relations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Other pressure-flow equations can be specified via callback supplied to the unit
configuration option ``pressure_flow_callback``. The callback allows both the form
and flow basis of the pressure-flow equation to be specified.

The callback can add parameters and variables as needed. It is recommended that
only the ``pressure_flow_equation`` be specified as additional constraints would
not be scaled by the valve model's scaling routines.  The pressure flow relation
generally should be written in the form below to facilitate scaling where :math:`F`
is flow variable.

.. math::

  f_1(F) = f_2(P_{in}, P_{out})

The callback takes a Valve model object as an argument. There are three attributes
that the ``pressure_flow_callback`` should define:
1. ``flow_var`` a time indexed reference to the flow variable basis,
2. ``pressure_flow_equation_scale`` a function that takes ``flow_var`` and defines the form of the flow term
3. ``pressure_flow_equation`` the pressure flow relation constraint.

The first two items, ``flow_var`` and ``pressure_flow_equation_scale``, are
not directly used in the model, but are used by the model scaling routine.

The example callback below is the model default pressure-flow equation.

.. code-block:: python

  def pressure_flow_default_callback(b):
      """
      Add the default pressure flow relation constraint.  This will be used in the
      valve model, a custom callback is provided.
      """
      umeta = b.config.property_package.get_metadata().get_derived_units

      b.Cv = pyo.Var(
          initialize=0.1,
          doc="Valve flow coefficent",
          units=umeta("amount")/umeta("time")/umeta("pressure")**0.5
      )
      b.Cv.fix()

      b.flow_var = pyo.Reference(b.control_volume.properties_in[:].flow_mol)
      b.pressure_flow_equation_scale = lambda x : x**2

      @b.Constraint(b.flowsheet().config.time)
      def pressure_flow_equation(b2, t):
          Po = b2.control_volume.properties_out[t].pressure
          Pi = b2.control_volume.properties_in[t].pressure
          F = b2.control_volume.properties_in[t].flow_mol
          Cv = b2.Cv
          fun = b2.valve_function[t]
          return F ** 2 == Cv ** 2 * (Pi - Po) * fun ** 2


Initialization
--------------

This just calls the initialization routine from PressureChanger. Either an
outlet pressure value or deltaP can be specified to aid the initialization.

Valve Class
-----------

.. autoclass:: Valve
  :members:

ValveData Class
---------------------

.. autoclass:: ValveData
  :members:
