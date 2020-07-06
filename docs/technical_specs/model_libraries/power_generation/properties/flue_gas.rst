Flue Gas Property Package
=========================

.. index::
   pair: idaes.power_generation.properties.IdealProp_FlueGas;FlueGasParameterBlock
   pair: idaes.power_generation.properties.IdealProp_FlueGas;FlueGasParameterData
   pair: idaes.power_generation.properties.IdealProp_FlueGas;FlueGasStateBlock
   pair: idaes.power_generation.properties.IdealProp_FlueGas;FlueGasStateBlockData

.. currentmodule:: idaes.power_generation.properties.IdealProp_FlueGas

A flue gas property package has been developed to provide properties of combustion gases and air. 
The ideal gas property package includes the main components in flue gas: O2, N2, NO, CO2, H2O, SO2

Main parameters:

* molecular weight in kg/kg-mol indexed by component list,
* reference pressure & temperature in Pa and Kelvin,
* critical pressure and temperature in Pa and Kelvin indexed by component list, 
* gas constant in J/(mol K),
* constants for specific heat capacity in J/(mol K) indexed by component list and parameter A to H,
* vapor pressure coefficients (Antoine Eq.) P in Bar and T in K indexed by component list and parameters A to C, 

Source: NIST webbook (last update: 01/08/2020)

The main methods supported are: 

* heat capacity in J/(mol K),
* enthalpy in J/mol,
* entropy in J/(mol K),
* volumetric flowrate m3/s,
* viscosity of mixture in kg/(m s),
* thermal conductivity mixture in J / (m K s),
* molar density m3/mol,
* reduced pressure and temperature (unitless),
