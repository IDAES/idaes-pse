Gas properties
==============

.. index::
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_properties;GasPhaseParameterBlock
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_properties;GasPhaseParameterData
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_properties;GasPhaseStateBlock
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_properties;GasPhaseStateBlockData

.. currentmodule:: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_properties

This property package provides the gas phase properties for the chemical looping combustion of methane.
The components modeled are methane, carbon dioxide, and water.

Main parameters:

* gas constant in kJ/(mol K),
* molecular weights in kg/mol indexed by component list,
* standard heat of formation in kJ/mol indexed by component list,
* constants for specific heat capacity in kJ/(mol K) indexed by component list and parameters 1 to 8,
* constants for dynamic viscosity in kg/(m s) indexed by component list and parameters 1 to 4,
* constants for thermal conductivity in kJ/(m K s) indexed by component list and parameters 1 to 4,
* constants for diffusivity of components in cm2/s indexed by component list

The main methods supported are:

* molecular weight of mixture in kg/mol
* molar heat capacity of component in kJ/(mol K)
* molar heat capacity of mixture in kJ/(mol K)
* mass heat capacity of mixture in kJ/(kg K)
* molar density of component in mol/m3
* molar density of mixture in mol/m3
* mass density of mixture in kg/m3
* molar enthalpy of component in kJ/mol
* molar enthalpy of mixture in kJ/mol
* dynamic viscosity of mixture in kg/(m s)
* thermal conductivity of mixture in kJ/(m K s)
* Diffusivity of components in cm2/s

References:

* National Institute of Standards and Technology, NIST Chemistry WebBook, (n.d.). (accessed March 10, 2018).
* R.H. Perry, D.W. Green, Perry’s Chemical Engineering Handbook, 1997, McGraw-Hill, n.d.
* Poling, B.E., Prausnitz, J.M. and O’connell, J.P., 2001. Properties of gases and liquids. McGraw-Hill Education.