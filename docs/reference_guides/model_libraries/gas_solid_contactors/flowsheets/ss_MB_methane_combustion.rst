Methane Combustion in a Moving Bed (steady-state)
=================================================

.. currentmodule:: idaes.models_extra.gas_solid_contactors.flowsheets.ss_MB_methane_combustion

Steady-state flowsheet example of the moving bed reactor model for methane combustion with an iron-oxide based oxygen carrier.

This model is for demonstration and tutorial purposes only.

The model code is located in the `main flowsheet module <https://github.com/IDAES/idaes-pse/blob/main/idaes/models_extra/gas_solid_contactors/flowsheets/ss_MB_methane_combustion.py>`_ and may be imported as:

.. code:: python

  >>> from idaes.models_extra.gas_solid_contactors.flowsheets.ss_MB_methane_combustion import main
  >>> m = main()

Inputs:

* Bed diameter
* Bed height
* Gas feed  - flowrate, pressure, temperature, composition
* Solid feed - flowrate, particle porosity, temperature, composition