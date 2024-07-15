Oxygen Carrier Oxidation in a Bubbling Fluidized Bed (steady-state)
===================================================================

.. currentmodule:: idaes.models_extra.gas_solid_contactors.flowsheets.ss_BFB_OC_oxidation

Steady-state flowsheet example of the bubbling fluidized bed model for oxidation of an iron-oxide based oxygen carrier with oxygen (air).

This model is for demonstration and tutorial purposes only.

The model code is located in the `main flowsheet module <https://github.com/IDAES/idaes-pse/blob/main/idaes/models_extra/gas_solid_contactors/flowsheets/ss_BFB_OC_oxidation.py>`_ and may be imported as:

.. code:: python

  >>> from idaes.models_extra.gas_solid_contactors.flowsheets.ss_BFB_OC_oxidation import main
  >>> m = main()

Inputs:

* Bed diameter
* Bed height
* Number of orifices
* Gas feed  - flowrate, pressure, temperature, composition
* Solid feed - flowrate, particle porosity, temperature, composition