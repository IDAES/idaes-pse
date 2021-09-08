Tagging Classes
===============

IDAES contains classes for tagging model quantities and grouping them.  The tags
allow a convenient short cut to important model inputs and outputs.  There are
also facilities for numeric formatting and displaying output in desired units.

Examples:
---------

The code below creates a simple model that can be used to demonstrate the use
of tags.

.. testcode::

  import pyomo.environ as pyo
  from idaes.core.util import ModelTag, ModelTagGroup

  def model():
      m = pyo.ConcreteModel()
      m.w = pyo.Var([1, 2, 3], ["a", "b"], initialize=4, units=pyo.units.kg)
      m.x = pyo.Var([1, 2, 3], initialize=5, units=pyo.units.kg)
      m.y = pyo.Var(initialize=6, units=pyo.units.s)
      m.z = pyo.Var(initialize=7)
      m.e = pyo.Expression(expr=m.w[1, "a"] / m.x[1])
      m.f = pyo.Expression(expr=m.x[1] / m.y)

      @m.Expression([1, 2, 3], ["a", "b"])
      def g(b, i, j):
          return m.w[i, j] / m.x[i] * 100

      return m

The next code snippet creates a single tag object for the model variable ``w``.
While the structure of the example model is simple, IDAES models often have a
complex structure. Tags generally provide a much shorter means to reference a
quantity.

.. testcode::

  m = model()
  tag = ModelTag(expr=m.w, format_string="{:.3f}", display_units=pyo.units.g)

Now we an use the tag to set model input and display model output.

.. testcode::

  # set all the elements of w with the second index of "a".
  tag[:,"a"].set(2*pyo.units.kg)
  tag[:,"b"].fix(3*pyo.units.kg)
  assert str(tag[1, "a"]) == "2000.000 g"
  assert str(tag[2, "a"]) == "2000.000 g"
  assert str(tag[3, "a"]) == "2000.000 g"
  assert str(tag[1, "b"]) == "3000.000 g"
  assert str(tag[2, "b"]) == "3000.000 g"
  assert str(tag[3, "b"]) == "3000.000 g"
  assert tag[3, "b"].expression.fixed

  # if no units are provided setting set_in_display_units to True will assume
  # the display units. If it is False, the native units of the quantity will be
  #used.
  tag.set_in_display_units = True
  tag[:,"a"].set(2)
  assert str(tag[1, "a"]) == "2.000 g"
  assert str(tag[2, "a"]) == "2.000 g"
  assert str(tag[3, "a"]) == "2.000 g"

  tags.str_include_units = False
  assert str(tag[1, "a"]) == "2.000"

In addition to creating single tag objects, a tag group can be created.  The
ModelTagGroup class is a dictionary with added methods for dealing with groups
of tags.  When tags are in groups ``set_in_display_units`` and
``str_include_units`` are set for the group as a whole and cannot be set
independently.

.. testcode::

  m = model()
  g = ModelTagGroup()

  g["w"] = ModelTag(expr=m.w, format_string="{:.3f}")
  g["x"] = ModelTag(expr=m.x, format_string="{:.3f}", display_units=pyo.units.g)
  g["y"] = ModelTag(expr=m.y, format_string="{:.3f}")
  g["z"] = ModelTag(expr=m.z, format_string="{:.3f}")
  g["e"] = ModelTag(expr=m.e, format_string="{:.3f}")
  g["f"] = ModelTag(expr=m.f, format_string="{:.3f}")
  g["g"] = ModelTag(expr=m.g, format_string="{:.3f}")

  g.set_in_display_units = True
  g.str_include_units = False

  g["x"].set(2)
  g["x"].setlb(1)
  g["x"].setub(3)

  assert str(g["x"][1]) == "2.000"
  assert abs(g["x"][1].expression.lb - 0.001) < 1e-5 # x is in kg
  assert abs(g["x"][1].expression.ub - 0.003) < 1e-5 # x is in kg


Available Methods
-----------------

.. autoclass:: idaes.core.util.tags.ModelTag
  :members:

.. autoclass:: idaes.core.util.tags.ModelTagGroup
  :members:
