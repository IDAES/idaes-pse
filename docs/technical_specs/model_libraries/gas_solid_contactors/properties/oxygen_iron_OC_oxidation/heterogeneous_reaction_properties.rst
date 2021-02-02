Heterogeneous reaction properties
=================================

.. index::
   pair: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.hetero_reactions;HeteroReactionThermoParameterBlock
   pair: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.hetero_reactions;HeteroReactionThermoParameterData
   pair: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.hetero_reactions;HeteroReactionThermoStateBlock
   pair: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.hetero_reactions;HeteroReactionThermoStateBlockData

.. currentmodule:: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.hetero_reactions

Property package for the reaction of air (oxygen) with an iron-based OC. 
The  gas components modeled are O2, N2, CO2, H2O.
The solid components modeled are Fe2O3, Fe3O4, Al2O3.
Reaction scheme modeled is O2 + 4Fe3O4 => 6Fe2O3.

Main parameters:

* gas constant in kJ/(mol K),
* standard heat of reaction in kJ/mol,
* grain radius in m,
* molar density of oxygen carrier particle in mol/m3
* available volume for reaction in m3/(m3 OC),
* activation energy of reaction in kJ/mol,
* reaction order of reaction in (dimensionless)
* pre-exponential factor of reaction in :math:`mol^{1-n}m^{3n - 2}s^{-1}` where n is reaction order

The main methods supported are:

* reaction rate constant in :math:`mol^{1-n}m^{3n - 2}s^{-1}` where n is reaction order
* oxygen carrier conversion in (dimensionless)
* reaction rate in mol_rxn/(m3 s)

References:

* A. Abad, J. Adánez, F. García-Labiano, L.F. de Diego, P. Gayán, J. Celaya, Mapping of the range of operational conditions for Cu-, Fe-, and Ni-based oxygen carriers in chemical-looping combustion, Chem. Eng. Sci. 62 (2007) 533–549.