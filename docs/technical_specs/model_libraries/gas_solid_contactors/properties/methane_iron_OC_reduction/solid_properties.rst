Solid properties
================

.. index::
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_properties;SolidPhaseThermoParameterBlock
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_properties;SolidPhaseThermoParameterData
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_properties;SolidPhaseThermoStateBlock
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_properties;SolidPhaseThermoStateBlockData

.. currentmodule:: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_properties

This property package provides the solid phase properties for the chemical looping combustion of an iron-based oxygen carrier.
The components modeled are Fe2O3, Fe3O4, and Al2O3.

Main parameters:

* molecular weight of components in kg/mol indexed by component list,
* skeletal densities of components in kg/m3 indexed by component list,
* standard heat of formation in kJ/mol indexed by component list,
* constants for specific heat capacity in kJ/(mol K) indexed by component list and parameters 1 to 8,
* particle diameter in m,
* minimum fluidization velocity in m/s,
* minimum fluidization voidage in (dimensionless)
* thermal conductivity of oxygen carrier in kJ/(m K s)

The main methods supported are:

* skeletal density of oxygen carrier in kg/m3,
* particle density of oxygen carrier in kg/m3,
* molar heat capacity of component in kJ/(mol K)
* mass heat capacity of oxygen carrier in kJ/(kg K)
* molar enthalpy of component in kJ/mol
* mass enthalpy of oxygen carrier in kJ/kg

References:

* National Institute of Standards and Technology, NIST Chemistry WebBook, (n.d.). (accessed March 10, 2018).
* R. Stevens, R. Newby, V. Shah, N. Kuehn, D. Keairns, Guidance for NETL’s oxycombustion R&D program: Chemical looping combustion reference plant designs and sensitivity studies, 2014.