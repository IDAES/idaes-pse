Heterogeneous reaction properties
=================================

.. index::
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.hetero_reactions;HeteroReactionParameterBlock
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.hetero_reactions;HeteroReactionParameterData
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.hetero_reactions;HeteroReactionStateBlock
   pair: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.hetero_reactions;HeteroReactionStateBlockData

.. currentmodule:: idaes.gas_solid_contactors.properties.methane_iron_OC_reduction.hetero_reactions

Property package for the reaction of CH4 with an iron-based OC. More details of this reaction scheme can be found in reference [1].
The gas components modeled are CH4, CO2, H2O.
The solid components modeled are Fe2O3, Fe3O4, Al2O3.
Reaction scheme modeled is CH4 + 12Fe2O3 => 8Fe3O4 + CO2 + 2H2O

**Rate basis**: Molar

**Units**: SI units

**Lists**:

* Component list - [CH4, CO2, H2O, Fe2O3, Fe3O4, Al2O3]
* Reaction list - [R1]

**Parameters**:

.. csv-table::
   :header: "Parameter Name", "Symbol", "Description", "Units", "Reference"

   "``dh_rxn``", ":math:`H_{rxn}`", "Heat of reaction", "kJ/mol",
   "``grain_radius``", ":math:`r_{g}`", "Representative particle grain radius within oxygen carrier particle", "m", "[1]"
   "``dens_mol_sol``", ":math:`\rho_{sol,mol}`", "Molar density of oxygen carrier particle", "mol/m3", "[1]"
   "``a_vol``", ":math:`a_{vol}`", "Available reaction volume per volume of oxygen carrier", "None", "[1]"
   "``k0_rxn``", ":math:`k_{0}`", "Pre-exponential factor", ":math:`\frac{mol^{1-n_{i}} m^{3n_{i}-2}}{s}`", "[1]"
   "``energy_activation``", ":math:`E_{A}`", "'Activation energy", "kJ/mol", "[1]"   
   "``rxn_order``", ":math:`n_{i}`", "Reaction order indexed by reaction list", "None", "[1]"
   "``gas_const``", ":math:`R`", "Gas constant", "kJ/mol.K",
   "``rate_reaction_stoichiometry``", ":math:`b_{i,j}`", "Reaction Stoichiometry indexed by reaction list and component list", "kJ/mol.K", "[1]"      
   "``mw_comp``", ":math:`mw_{j}`", "Molecular weights of components indexed by component list", "kg/mol", "[1]"   
   
**Variables**:

.. csv-table::
   :header: "Variable Name", "Symbol", "Description", "Units"

   "``k_rxn``", ":math:`k`", "Rate constant", ":math:`\frac{mol^{1-n_{i}} m^{3n_{i}-2}}{s}`"
   "``OC_conv``", ":math:`X`","Fraction of oxygen carrier converted", "None"
   "``reaction_rate``", ":math:`rate_{rxn}`","Reaction rate", "mol_rxn/m3.s"
   "``temperature``", ":math:`T`","Temperature", "K"
   "``dens_mol_comp``", ":math:`C_{g,j}`", "Component molar concentration indexed by component list", "mol/m3"
   "``mass_frac_comp``", ":math:`x_{j}`","Component mass fractions indexed by component list", "None"
   "``particle_porosity``", ":math:`\phi`","Porosity of oxygen carrier", "None"   
   "``dens_mass_skeletal``", ":math:`\rho_{skeletal}`", "Skeletal density of oxygen carrier", "kg/m3"   
      
**Methods**:

Rate constant:

.. math:: k = k_{0}exp{\left( \frac{-E}{R T} \right)}

Reaction rate:

.. math:: rate_{rxn} = x_{Fe2O3}{\left( 1 - \phi \right)}\rho_{skeletal} \frac{a_{vol}}{mw_{Fe2O3}} \frac{{3} {k} b_{i,Fe2O3} {C_{g, CH4}}^{n_{i}}}{\rho_{sol,mol} r_{g}} {\left(1 - X\right)}^{\frac{2}{3}}

**References**:

1. A. Abad, J. Adánez, F. García-Labiano, L.F. de Diego, P. Gayán, J. Celaya, Mapping of the range of operational conditions for Cu-, Fe-, and Ni-based oxygen carriers in chemical-looping combustion, Chem. Eng. Sci. 62 (2007) 533–549.