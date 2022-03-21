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

**Flow basis**: Molar

**Units**: SI units

**State Variables**: 

The state block supports the following state variables:

* Component molar flowrate in mol/s,
* Pressure in bar,
* Component mole fraction in (dimensionless),
* Temperature in K

**Lists**:

* Component list - [CH4, CO2, H2O]
* Shomate parameter list - [1 to 8]
* Viscosity parameter list - [1 to 4]
* Thermal conductivity parameter list - [1 to 4]

**Parameters**:

.. csv-table::
   :header: "Parameter Name", "Symbol", "Description", "Units", "Reference"

   "``mw_comp``", ":math:`mw_{j}`", "Molecular weights of gas components indexed by component list", "kg/mol", "[1]"
   "``enth_mol_form_comp``", ":math:`H_{form,j}`", "Component molar heats of formation indexed by component list", "J/mol", "[1]"
   "``cp_param``", ":math:`CP_{j,i}`", "Heat capacity parameters indexed by component list and shomate list",  , "[1]"
   "``visc_d_param``", ":math:`\mu_{param,j,i}`", "Viscosity parameters indexed by component list and viscosity list",  , "[2]"
   "``therm_cond_param``", ":math:`k_{param,j,i}`", "Thermal conductivity parameters indexed by component list and thermal conductivity list",  , "[2]"   
   "``diff_vol_param``", ":math:`V_{param,j}`", "Diffusion volume parameters indexed by component list",  , "[3]"
   "``gas_const``", ":math:`R`", "Gas constant", "kJ/mol.K",     

**Variables**:

.. csv-table::
   :header: "Variable Name", "Symbol", "Description", "Units"

   "``flow_mol``", ":math:`F_{mass}`", "Component molar flowrate", "mol/s"
   "``pressure``", ":math:`P`","Pressure", "bar"
   "``mole_frac_comp``", ":math:`y_{j}`","Component mole fractions indexed by component list", "None"
   "``temperature``", ":math:`T`","Temperature", "K" 
   "``mw``", ":math:`mw`", "Molecular weight of gas mixture", "kg/mol"
   "``dens_mol``", ":math:`C_{g}`", "Molar density/concentration", "mol/m3"
   "``dens_mol_comp``", ":math:`C_{g,j}`", "Component molar concentration indexed by component list", "mol/m3"
   "``dens_mass``", ":math:`\rho_{mass}`", "Mass density", "kg/m3"
   "``visc_d``", ":math:`\mu_{vap}`", "Mixture dynamic viscosity", "kg/m.s"
   "``diffusion_comp``", ":math:`D_{vap,j}`", "Component diffusion in a gas mixture indexed by component list", "cm2/s"
   "``therm_cond``", ":math:`k_{vap}`", "Thermal conductivity of gas", "kJ/m.K.s"
   "``cp_mol_comp``", ":math:`c_{p,mol,j}`", "Pure component molar heat capacities indexed by component list", "J/mol.K"
   "``cp_mol``", ":math:`c_{p,mol}`", "Mixture heat capacity, mole-basis", "J/mol.K"
   "``cp_mass``", ":math:`c_{p,mass}`", "Mixture heat capacity, mass-basis", "J/kg.K"
   "``enth_mol_comp``", ":math:`H_{mol,j}`", "Pure component enthalpies indexed by component list", "J/mol"
   "``enth_mol``", ":math:`H_{mol}`", "Molar enthalpy of gas mixture", "J/mol"
   
**Methods**:

Sum of component fractions:

.. math:: 1 = \sum_j{y_{j}}

Molecular weight of gas mixture:

.. math:: mw = \sum_j{y_{j}mw_{j}}

Molar density:

.. math:: C_{g} = 100 \frac{P}{RT}

Component molar density:

.. math:: C_{g,j} = y_{j}C_{g}

Mass density:

.. math:: \rho_{mass} = mwC_{g}

Mixture dynamic viscosity, see reference [2] for parameters:

.. math:: \mu_{vap} = \sum_i{\frac{y_{i}\mu_{i}}{\sum_j{y_{j}{\left(\frac{mw_{j}}{mw_{i}}\right)}^{0.5}}}}
.. math:: \mu_{i} = \frac{\mu_{param,j,1}T^{\mu_{param,j,2}}}{1 + \frac{\mu_{param,j,3}}{T} + \frac{\mu_{param,j,4}}{T^2}}

Thermal conductivity, see reference [2] for parameters:

.. math:: k_{vap} = \sum_i{\frac{y_{i}k_{i}}{\sum_j{y_{j}A_{j,i}}}}
.. math:: k_{i} = \frac{k_{param,j,1}T^{k_{param,j,2}}}{1 + \frac{k_{param,j,3}}{T} + \frac{k_{param,j,4}}{T^2}}
.. math:: A_{j,i} = \frac{{\left(1 +  {\left( \frac{k_{j}}{k_{i}} \right)}^{0.5}  {\left( \frac{mw_{j}}{mw_{i}} \right)}^{0.25}  \right)}^{2}} {{{8{\left(1 + {\left(\frac{mw_{j}}{mw_{i}}\right)}  \right)}}}^{0.5}}

Diffusion of component in a multicomponent gas mixture, see reference [3] for parameters:

.. math:: D_{vap,j} = \frac{1 - y_{j}} {\sum_{j,j ≠ i}{\frac{y_{i}}{D_{j,i}}}}
.. math:: D_{j,i} = \frac{{0.00143} T^{1.75} {\left( \frac{1}{mw_{j}} + \frac{1}{mw_{i}} \right)}^{0.5} } {{P {\left( {V_{param,j}^{\frac{1}{3}}} {V_{param,i}^{\frac{1}{3}}}\right)} }^{2}} 

Molar heat capacity of component, see reference [1]:

.. math:: c_{p,mol,j} = CP_{j,1} + CP_{j,2}\bar{T} + CP_{j,3}\bar{T}^2 + CP_{j,4}\bar{T}^3 + \frac{CP_{j,5}}{\bar{T}^2}
.. math:: \bar{T} = 10^{-3}T

Molar heat capacity of gas mixture:

.. math:: c_{p,mol} = \sum_j{c_{p,mol,j}y_{j}}

Mass heat capacity of gas mixture:

.. math:: c_{p,mass} = \frac{c_{p,mol}}{mw}

Molar enthalpy of component, see reference [1]:

.. math:: H_{mol,j} = P_{j,1} + P_{j,2}\bar{T} + P_{j,3}\bar{T}^2 + P_{j,4}\bar{T}^3 + \frac{P_{j,5}}{\bar{T}^2}
.. math:: \bar{T} = 10^{-3}T

Molar enthalpy of gas mixture:

.. math:: H_{mole} = \sum_j{H_{mol,j}y_{j}}


**References**:

1. National Institute of Standards and Technology, NIST Chemistry WebBook, (n.d.). (accessed March 10, 2018).
2. R.H. Perry, D.W. Green, Perry’s Chemical Engineering Handbook, 1997, McGraw-Hill, n.d.
3. Poling, B.E., Prausnitz, J.M. and O’connell, J.P., 2001. Properties of gases and liquids. McGraw-Hill Education.