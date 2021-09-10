Tagging Classes
===============

IDAES contains classes for tagging model quantities and grouping them.  The tags
provide a convenient short cut to important model inputs and outputs and
facilities for numeric formatting and displaying output in desired units.

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
complex structure, so tags provide a much shorter means to reference a quantity.

.. testcode::

  m = model()
  tag = ModelTag(expr=m.w, format_string="{:.3f}", display_units=pyo.units.g)

Now we can use the tag to set model input and display model output.

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

  tag.str_include_units = False
  assert str(tag[1, "a"]) == "2.000"

In addition to creating single tag objects, a tag group can be created.  The
ModelTagGroup class is a dictionary with added methods for dealing with groups
of tags.  When tags are in groups ``set_in_display_units`` and
``str_include_units`` are set for the group as a whole and cannot be set
independently.

.. testcode::

  m = model()
  group = ModelTagGroup()

  group["w"] = ModelTag(expr=m.w, format_string="{:.3f}")
  group["x"] = ModelTag(expr=m.x, format_string="{:.3f}", display_units=pyo.units.g)
  group["y"] = ModelTag(expr=m.y, format_string="{:.3f}")
  group["z"] = ModelTag(expr=m.z, format_string="{:.3f}")
  group["e"] = ModelTag(expr=m.e, format_string="{:.3f}")
  group["f"] = ModelTag(expr=m.f, format_string="{:.3f}")
  group["g"] = ModelTag(expr=m.g, format_string="{:.3f}")

  group.set_in_display_units = True
  group.str_include_units = False

  group["x"].set(2)
  group["x"].setlb(1)
  group["x"].setub(3)

  assert str(group["x"][1]) == "2.000"
  assert abs(group["x"][1].expression.lb - 0.001) < 1e-5 # x is in kg
  assert abs(group["x"][1].expression.ub - 0.003) < 1e-5 # x is in kg

When a tagged a quantity can vary over several orders of magnitude, it can be
helpful to provide conditional formatting. To do this a callable can be provided
as the ``format_string`` which takes the quantity value and returns a format
string. A simple example is given below.

.. testcode::

  m = model()

  tagw = ModelTag(
    expr=m.w,
    format_string=lambda x: "{:,.0f}" if x >= 100 else "{:.2f}",
    display_units=pyo.units.g,
  )

  tagw.set(1*pyo.units.g)
  assert str(tagw[1, "a"]) == "1.00 g"
  tagw.set(1*pyo.units.kg)
  assert str(tagw[1,"a"]) == "1,000 g"

Available Classes
-----------------

.. autoclass:: idaes.core.util.tags.ModelTag
  :members:

.. autoclass:: idaes.core.util.tags.ModelTagGroup
  :members:

.. autofunction:: idaes.core.util.tags.svg_tag
