Solid properties
================

.. index::
   pair: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.solid_properties;SolidPhaseParameterBlock
   pair: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.solid_properties;SolidPhaseParameterData
   pair: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.solid_properties;SolidPhaseStateBlock
   pair: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.solid_properties;SolidPhaseStateBlockData

.. currentmodule:: idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.solid_properties

This property package provides the solid phase properties for the chemical looping combustion of an iron-based oxygen carrier.
The components modeled are Fe2O3, Fe3O4, and Al2O3.

**Flow basis**: Mass

**Units**: SI units

**State Variables**: 

The state block supports the following state variables:

* Component mass flowrate in kg/s,
* Particle porosity in (dimensionless),
* Component mass fraction in (dimensionless),
* Temperature in K

**Lists**:

* Component list - [Fe2O3, Fe3O4, Al2O3]
* Shomate parameter list - [1 to 8]


**Parameters**:

.. csv-table::
   :header: "Parameter Name", "Symbol", "Description", "Units", "Reference"

   "``mw_comp``", ":math:`mw_{j}`", "Molecular weights of solid components indexed by component list", "kg/mol", "[1]"
   "``dens_mass_comp_skeletal``", ":math:`\rho_{skeletal,j}`", "Skeletal density of solid components indexed by component list", "kg/m3", "[1]"
   "``cp_param``", ":math:`CP_{j,i}`", "Heat capacity parameters indexed by component list and shomate list", "J/ mol.K", "[1]"
   "``enth_mol_form_comp``", ":math:`H_{form,j}`", "Component molar heats of formation indexed by component list", "J/mol", "[1]"
   "``particle_dia``", ":math:`d_{p}`", "Diameter of solid particles", "m", "[2]"
   "``velocity_mf``", ":math:`v_{mf}`", "Velocity at minimum fluidization", "m/s", "[2]"
   "``voidage_mf``", ":math:`\varepsilon_{mf}`", "Voidage at minimum fluidization", "None", "[2]"
   "``therm_cond_sol``", ":math:`k_{sol}`", "Thermal conductivity of solid particles", "kJ/m.K.s", "[2]"
   

**Variables**:

.. csv-table::
   :header: "Variable Name", "Symbol", "Description", "Units"

   "``flow_mass``", ":math:`F_{mass}`", "Component mass flowrate", "kg/s"
   "``particle_porosity``", ":math:`\phi`","Porosity of oxygen carrier", "None"
   "``mass_frac_comp``", ":math:`x_{j}`","Component mass fractions indexed by component list", "None"
   "``temperature``", ":math:`T`","Temperature", "K" 
   "``dens_mass_skeletal``", ":math:`\rho_{skeletal}`", "Skeletal density of oxygen carrier", "kg/m3"
   "``dens_mass_particle``", ":math:`\rho_{particle}`", "Particle density of oxygen carrier", "kg/m3"
   "``cp_mol_comp``", ":math:`c_{p,mol,j}`", "Pure component solid heat capacities indexed by component list", "J/mol.K"
   "``cp_mass``", ":math:`c_{p,mass}`", "Mixture heat capacity, mass-basis", "J/kg.K"
   "``enth_mol_comp``", ":math:`H_{mol,j}`", "Pure component enthalpies indexed by component list", "J/mol"
   "``enth_mass``", ":math:`H_{mass}`", "Mixture specific enthalpy", "J/kg"


**Methods**:

Sum of component fractions:

.. math:: 1 = \sum_j{x_{j}}

Skeletal density of oxygen carrier:

.. math:: \rho_{skeletal} = \frac{1} {\sum_j{\frac{x_{j}}{\rho_{skeletal,j}}}}

Particle density of oxygen carrier:

.. math:: \rho_{particle} = {\left(1 - \phi \right)} \rho_{skeletal}

Molar heat capacity of component, see reference [1]:

.. math:: c_{p,mol,j} = CP_{j,1} + CP_{j,2}\bar{T} + CP_{j,3}\bar{T}^2 + CP_{j,4}\bar{T}^3 + \frac{CP_{j,5}}{\bar{T}^2}
.. math:: \bar{T} = 10^{-3}T

Mass heat capacity of oxygen carrier:

.. math:: c_{p,mass} = \sum_j{\frac{c_{p,mol,j}x_{j}}{mw_{j}}}

Molar enthalpy of component, see reference [1]:

.. math:: H_{mol,j} = P_{j,1} + P_{j,2}\bar{T} + P_{j,3}\bar{T}^2 + P_{j,4}\bar{T}^3 + \frac{P_{j,5}}{\bar{T}^2}
.. math:: \bar{T} = 10^{-3}T

Mass enthalpy of oxygen carrier:

.. math:: H_{mass} = \sum_j{\frac{H_{mol,j}x_{j}}{mw_{j}}}


**References**:

1. National Institute of Standards and Technology, NIST Chemistry WebBook, (n.d.). (accessed March 10, 2018).
2. Stevens R., Newby R., Shah V., Kuehn N., Keairns D., Guidance for NETL’s oxycombustion R&D program: Chemical looping combustion reference plant designs and sensitivity studies, 2014.