Thermogravimetric Analysis in a Fixed Bed (0D dynamic)
=============================================================

.. currentmodule:: idaes.models_extra.gas_solid_contactors.flowsheets.dyn_TGA_example

Dynamic flowsheet example of the fixed bed 0D model for thermogravimetric analysis.

This model is for demonstration and tutorial purposes only.

The model code is located in the `main flowsheet module <https://github.com/IDAES/idaes-pse/blob/main/idaes/models_extra/gas_solid_contactors/flowsheets/dyn_TGA_example.py>`_ and may be imported as:

.. code:: python

  >>> from pyomo.environ import ConcreteModel
  >>> from idaes.models_extra.gas_solid_contactors.flowsheets.dyn_TGA_example import main
  >>> m = ConcreteModel()
  >>> m = main(m)

Inputs:

* Gas feed  - flowrate, pressure, temperature, composition
* Solid feed - flowrate, particle porosity, temperature, composition