Bubbling Fluidized Bed Reactor
========================================

The IDAES Bubbling Fluidized Bed Reactor (BFBR) model represents a unit operation where two material streams, 
a solid phase and a gas phase, pass through a linear vessel while undergoing chemical reaction(s). 
The BFBR model is represented as a 1-D axially discretized model with two phases (gas and solid), 
and two regions (bubble and emulsion). The model captures the gas-solid interaction between both phases and regions 
through reaction, mass and heat transfer.

**Assumptions:**

* Cloud-wake region effects are negligigble and are not modelled.
* Gas emulsion is at minimum fluidization conditions.
* Gas feeds into emulsion region before the excess enters into the bubble region.
* Gas and solids are well mixed in the radial direction but vary axially.

**Requirements:**

* Property package contains temperature and pressure variables.
* Property package contains minimum fluidization velocity and voidage parameters.

**The BFBR model equations are derived from:**

* A. Lee, D.C. Miller. A one-dimensional (1-D) three-region model for a bubbling fluidized-bed adsorber, Ind. Eng. Chem. Res. 52 (2013) 469–484.

Degrees of Freedom
------------------

BFBRs generally have at least 3 (or more) degrees of freedom, consisting of design and operating variables. The design variables of reactor length, diameter and number of orifices in the distributor are typically the minimum variables to be fixed.

Model Structure
---------------

The core BFBR unit model consists of two inlet ports (named gas_inlet and solid_inlet),
two outlet ports (named gas_outlet and solid_outlet), and three ControlVolume1DBlock 
Blocks (named bubble_region, gas_emulsion_region and solid_emulsion_region).


Construction Arguments
----------------------

The IDAES BFBR model has construction arguments specific to the whole unit and to the individual regions.

**Arguments that are applicable to the BFBR unit as a whole are:**

* finite_elements - sets the number of finite elements to use when discretizing the spatial domains (default = 10).
* length_domain_set - sets the list of point to use to initialize a new ContinuousSet (default = [0.0, 1.0]).
* transformation_method - sets the discretization method to use by the Pyomo TransformationFactory 
  to transform the spatial domain (default = dae.finite_difference):
  
  - dae.finite_difference - finite difference method.
  - dae.collocation - orthogonal collocation method.
  
* transformation_scheme - sets the scheme to use when transforming a domain. 
  Selected schemes should be compatible with the transformation_method chosen (default = None):
  
  - None - defaults to "BACKWARD" for finite difference transformation method and to "LAGRANGE-RADAU" for collocation transformation method
  - BACKWARD - use a finite difference transformation method.
  - FORWARD - use a finite difference transformation method.
  - LAGRANGE-RADAU - use a collocation transformation method.   
  
* collocation_points - sets the number of collocation points to use when discretizing the spatial domains (default = 3, collocation methods only).
* flow_type - indicates the flow arrangement within the unit to be modeled. Options are:

  - 'co-current' - (default) gas and solid streams both flow in the same direction (from x=0 to x=1)
  - 'counter-current' - gas and solid streams flow in opposite directions (gas from x=0 to x=1 and solid from x=1 to x=0).

* material_balance_type - indicates what type of energy balance should be constructed (default = MaterialBalanceType.componentTotal).

  - MaterialBalanceType.componentTotal - use total component balances.
  - MaterialBalanceType.total - use total material balance.
    
* energy_balance_type - indicates what type of energy balance should be constructed (default = EnergyBalanceType.enthalpyTotal).

  - EnergyBalanceType.none - excludes energy balances.
  - EnergyBalanceType.enthalpyTotal - single enthalpy balance for material.

* momentum_balance_type - indicates what type of momentum balance should be constructed (default = MomentumBalanceType.pressureTotal).

  - MomentumBalanceType.none - exclude momentum balances.
  - MomentumBalanceType.pressureTotal - single pressure balance for material.

* has_pressure_change - indicates whether terms for pressure change should be constructed (default = True).

  - True - include pressure change terms.
  - False - exclude pressure change terms.

**Arguments that are applicable to the gas phase:**

* property_package - property package to use when constructing bubble region Property Blocks (default = 'use_parent_value'). 
  This is provided as a Physical Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the bubble region Property Blocks when they are created (default = 'use_parent_value').
* reaction_package - reaction package to use when constructing bubble region Reaction Blocks (default = None). 
  This is provided as a Reaction Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* reaction_package_args - set of arguments to be passed to the bubble region Reaction Blocks when they are created (default = None).
* has_equilibrium_reactions - sets flag to indicate if terms of equilibrium controlled reactions should be constructed (default = False).

**Arguments that are applicable to the solid phase:**

* property_package - property package to use when constructing bubble region Property Blocks (default = 'use_parent_value'). 
  This is provided as a Physical Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the bubble region Property Blocks when they are created (default = 'use_parent_value').
* reaction_package - reaction package to use when constructing bubble region Reaction Blocks (default = None). 
  This is provided as a Reaction Parameter Block by the Flowsheet when creating the model. 
  If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* reaction_package_args - set of arguments to be passed to the bubble region Reaction Blocks when they are created (default = None).
* has_equilibrium_reactions - sets flag to indicate if terms of equilibrium controlled reactions should be constructed (default = False).

Additionally, BFBR units have the following construction arguments which are passed to all the ControlVolume1DBlock Blocks 
and are always specified to their default values.

========================= ==========================
Argument                  Default Value
========================= ==========================
dynamic                   useDefault
has_holdup                useDefault
========================= ==========================

Constraints
-----------

**Geometric Constraints**

Area of orifice:

.. math:: A_{or} = \frac{1}{n_{or}}

Bed cross-sectional area:

.. math:: A_{bed} = \pi \frac{D_{bed}^{2}}{4}

Area of bubble region:

.. math:: A_{b,t,x} = \delta_{t,x} A_{bed}

Area of gas emulsion region:

.. math:: A_{ge,t,x} = \delta_{e,t,x} \varepsilon_{e,t,x} A_{bed}

Area of solid emulsion region:

.. math:: A_{se,t,x} = \delta_{e,t,x} {\left(1 - \varepsilon_{e,t,x} \right)} A_{bed}

Length of bubble region:

.. math:: L_{b} = L_{bed}

Length of gas emulsion region:

.. math:: L_{ge} = L_{bed}

Length of solid emulsion region:

.. math:: L_{se} = L_{bed}

**Hydrodynamic Constraints**

Emulsion region volume fraction:

.. math:: \delta_{e,t,x} = 1 -\delta_{t,x}

Average cross-sectional voidage:

.. math:: \varepsilon_{t,x} = 1 - \left(1 - \varepsilon_{e,t,x} \right) \left(1 - \delta_{t,x} \right)

Emulsion region voidage:

.. math:: \varepsilon_{e,t,x} = \varepsilon_{mf,se}

Bubble growth coefficient:

.. math:: \gamma_{t,x} = \frac{0.0256}{v_{mf,se}} {\left(\frac{D_{bed}}{g} \right)}^{0.5}

Maximum bubble diameter:

.. math:: d_{bm,t,x}^{5}g = 2.59^{5} {\left([v_{g,t,x} - v_{ge,t,x}] A_{bed} \right)}^{2}

Bubble diameter (gas inlet, x = 0):

.. math:: d_{b,t,x} = 1.38g^{-0.2} {\left([v_{g,t,x} - v_{ge,t,x}]A_{or} \right)}^{0.4}

Bubble diameter (x > 0):

.. math:: \frac{dd_{b,t,x}}{ dx } = \frac{0.3}{D_{bed}} L_{bed} {\left(d_{bm,t,x} - d_{b,t,x} - \gamma_{t,x}{\left(D_{bed} d_{b,t,x}\right)}^{0.5}\right)}

Bubble rise velocity:

.. math:: v_{br,t,x}^{2} = 0.711^{2} g d_{b,t,x}

Bubble velocity:

.. math:: v_{b,t,x} = v_{g,t,x} - v_{mf,se} + v_{br,t,x}

Average gas density:

.. math:: C_{g,t,x} = \frac{F_{mol,b,t,x} C_{b,total,t,x} + F_{mol,ge,t,x} C_{ge,total,t,x}} {F_{mol,b,t,x} + F_{mol,ge,t,x}}

Superficial gas velocity:

.. math:: v_{g,t,x} = \frac{F_{mol,b,t,x} + F_{mol,ge,t,x}} {A_{bed} C_{g,t,x}}
                  
Bubble volume fraction:

.. math:: \delta_{t,x} v_{b,t,x} =  v_{ge,t,x}\delta_{e,t,x} - v_{g,t,x}

Gas emulsion pressure drop:

    if 'has_pressure_change' is 'True':

    .. math:: \Delta P_{ge,t,x} = - g {\left(1 - \varepsilon_{e,t,x} \right)} \rho_{mass,se,t,x}

    elif 'has_pressure_change' is 'False':
    
    .. math:: P_{ge,t,x} = P_{g,t,inlet}

**Mass Transfer Constraints**

Bubble to emulsion gas mass transfer coefficient:

.. math:: K_{be,t,x,j} d_{b,t,x}^{1.25} = 5.94 v_{mf,se} d_{b,t,x}^{0.25} + 5.85 D_{vap,ge,t,x,j}^{0.5} g^{0.25}

Bulk gas mass transfer:

    if :math:`C_{ge,total,t,x}` > :math:`C_{b,total,t,x}`:

    .. math:: K_{gbulkc,t,x,j} = 6 K_{d} \delta_{t,x} A_{bed} {\left(C_{ge,total,t,x} - C_{b,total,t,x} \right)} d_{b,t,x} y_{ge,t,x,j}

    else:

    .. math:: K_{gbulkc,t,x,j} = 6 K_{d} \delta_{t,x} A_{bed} {\left(C_{ge,total,t,x} - C_{b,total,t,x} \right)} d_{b,t,x} y_{b,t,x,j}
    
**Heat Transfer Constraints**

Bubble to emulsion gas heat transfer coefficient:

.. math:: H_{be,t,x,j} d_{b,t,x}^{1.25} = 4.5 v_{mf,se} c_{p\_vap,b,t,x} C_{b,total,t,x} d_{b,t,x}^{0.25} + 5.85 {\left(k_{vap,b,t,x}  C_{b,total,t,x} c_{p\_vap,b,t,x} \right)}^{0.5} g^{0.25}

Convective heat transfer coefficient:

.. math:: h_{tc,t,x} d_{p,se} = 0.03 k_{vap,e,t,x} {\left(v_{ge,t,x} d_{p,se} \frac{C_{ge,total,t,x}}{\mu_{vap,ge,t,x}} \right)}^{1.3}

Emulsion region gas-solids convective heat transfer:

.. math:: h_{t\_gs,t,x} d_{p,se} = 6 \delta_{e,t,x} {\left(1 - \varepsilon_{e,t,x} \right)} h_{tc,t,x} {\left(T_{ge,t,x} - T_{se,t,x} \right)}

Bulk gas heat transfer:

    if :math:`C_{ge,total,t,x}` > :math:`C_{b,total,t,x}`:

    .. math:: H_{gbulk,t,x} =  K_{d} \delta_{t,x} A_{bed} {\left(C_{ge,total,t,x} - C_{b,total,t,x} \right)} d_{b,t,x} H_{ge,t,x}

    else:

    .. math:: H_{gbulk,t,x} =  K_{d} \delta_{t,x} A_{bed} {\left(C_{ge,total,t,x} - C_{b,total,t,x} \right)} d_{b,t,x} H_{b,t,x}

**Mass and heat transfer terms in control volumes**

Bubble mass transfer '(p=vap)':

.. math:: M_{tr,b,t,x,p,j} = K_{gbulkc,t,x,j} - A_{b,t,x} K_{be,t,x,j} {\left(C_{b,total,t,x} - C_{ge,total,t,x} \right)}

Gas emulsion mass transfer '(p=vap)':

.. math:: M_{tr,ge,t,x,p,j} = - K_{gbulkc,t,x,j} + A_{b,t,x} K_{be,t,x,j} {\left(C_{b,total,t,x} - C_{ge,total,t,x} \right)} + r_{hetero,ge,t,x,j}

*if 'energy_balance_type' is not 'EnergyBalanceType.none':*

    Bubble heat transfer:

    .. math:: H_{tr, b, t,x} = H_{gbulk,t,x} - A_{b,t,x} H_{be,t,x,j} {\left(T_{b,t,x} - T_{ge,t,x} \right)}

    Gas emulsion heat transfer:

    .. math:: H_{tr, ge, t,x} = - H_{gbulk,t,x} + A_{b,t,x} H_{be,t,x,j} {\left(T_{b,t,x} - T_{ge,t,x} \right)} - h_{t\_gs,t,x} A_{bed}

    Solid emulsion heat transfer:

    .. math:: H_{tr, se, t,x} =  h_{t\_gs,t,x} A_{bed}

**Reaction constraints**

*if 'homogeneous reaction package' is not 'None':*
 
    Bubble rate reaction extent:

    .. math:: r_{ext,b,t,x,r} = A_{b,t,x} r_{b,t,x,r}

    Gas emulsion rate reaction extent:

    .. math:: r_{ext,ge,t,x,r} = A_{ge,t,x} r_{ge,t,x,r}

*if 'heterogeneous reaction package' is not 'None':*

    Solid emulsion rate reaction extent:

    .. math:: r_{ext,se,t,x,r} = A_{se,t,x} r_{se,t,x,r}

    Gas emulsion heterogeneous rate reaction extent:

    .. math:: r_{hetero,ge,t,x,j} = A_{se,t,x} \sum_{r}^{reactions} {s_{se,j,r} r_{se,t,x,r}}

**Flowrate constraints**

Bubble gas flowrate:

.. math:: F_{mol,b,t,x} = A_{bed} \delta_{t,x} v_{b,t,x} C_{b,total,t,x}

Emulsion gas flowrate:

.. math:: F_{mol,ge,t,x} = A_{bed} v_{ge,t,x} C_{ge,total,t,x}

**Inlet boundary conditions**

*if 'has_pressure_change' is 'True':*

    Gas emulsion pressure at inlet:

    .. math:: P_{ge,t,0} = P_{g,t,inlet} - \Delta P_{or}

Total gas balance at inlet:

.. math:: F_{mol,b,t,0} + F_{mol,ge,t,0} = F_{mol,g,t,inlet}

Solid particle porosity at inlet:

.. math:: \phi_{se,t,0} = \phi_{t,inlet}

Gas emulsion velocity at inlet:

.. math:: v_{ge,t,0} = v_{mf,se} 

Bubble mole fraction at inlet:

.. math:: y_{b,t,0,j} = y_{g,t,inlet,j}

Gas emulsion mole fraction at inlet:

.. math:: y_{ge,t,0,j} = y_{g,t,inlet,j}

Solid emulsion mass flow at inlet:

    *if 'flow_type' is 'co_current' x = 0 else if 'flow_type' is 'counter_current' x = 1:*

    .. math:: F_{mass,se,t,x} = F_{mass,s,t,inlet}

Solid emulsion mass fraction at inlet:

    *if 'flow_type' is 'co_current' x = 0 else if 'flow_type' is 'counter_current' x = 1:*

    .. math:: x_{se,t,x} = x_{s,t,inlet} 

*if 'energy_balance_type' is not 'EnergyBalanceType.none':*

    Gas inlet energy balance:

    .. math:: H_{b,t,0} + H_{ge,t,0} = H_{g,t,inlet}

    Gas emulsion temperature at inlet:

    .. math:: T_{ge,t,0} = T_{g,t,inlet}

*elif 'energy_balance_type' is 'EnergyBalanceType.none':*

    Isothermal bubble region: 

    .. math:: T_{b,t,x} = T_{g,t,inlet}
    
    Isothermal gas emulsion region:    
    
    .. math:: T_{ge,t,x} = T_{g,t,inlet}

    Isothermal solid emulsion region:      
    
    .. math:: T_{se,t,x} = T_{s,t,inlet}         

*if 'flow_type' is 'co_current' x = 0 else if 'flow_type' is 'counter_current' x = 1:*

    Solid inlet energy balance:

    .. math:: H_{se,t,x} = H_{s,t,inlet} 

**Outlet boundary conditions**

Gas emulsion pressure at outlet:

.. math:: P_{g,t,outlet} = P_{ge,t,1} 

Total gas balance at outlet:

.. math:: F_{mol,g,t,outlet} = F_{mol,b,t,1} + F_{mol,ge,t,1}

Solid outlet material balance:

    *if 'flow_type' is 'co_current' x = 1 else if 'flow_type' is 'counter_current' x = 0:*

    .. math:: F_{mass,s,t,outlet} = F_{mass,se,t,x}
    
Solid particle porosity at outlet:

.. math:: \phi_{t,outlet} = \phi_{se,t,1}

*if 'energy_balance_type' is not 'EnergyBalanceType.none':*

    Gas outlet energy balance:

    .. math:: H_{g,t,outlet} = H_{b,t,1} + H_{ge,t,1}

    Solid outlet energy balance:

        *if 'flow_type' is 'co_current' x = 1 else if 'flow_type' is 'counter_current' x = 0:*

        .. math:: H_{s,t,outlet} = H_{se,t,x}
        
*elif 'energy_balance_type' is 'EnergyBalanceType.none':*

    Gas outlet energy balance:

    .. math:: T_{g,t,outlet} = T_{ge,t,1}    

    Solid outlet energy balance:

        *if 'flow_type' is 'co_current' x = 1 else if 'flow_type' is 'counter_current' x = 0:*

        .. math:: T_{s,t,outlet} = T_{se,t,x}
        
Variables
---------

List of variables in the BFBR model:

=========================== ===================================================== ====================================================
Variable                    Name                                                  Notes
=========================== ===================================================== ====================================================
:math:`L_{bed}`             bed_height                                            Bed height
:math:`D_{bed}`             bed_diameter                                          Reactor diameter
:math:`A_{bed}`             bed_area                                              Reactor cross-sectional area
:math:`A_{or}`              area_orifice                                          Distributor plate area per orifice
:math:`n_{or}`              number_orifice                                        Number of distributor plate orifices per area
:math:`\delta_{t,x}`        delta                                                 Volume fraction occupied by bubble region
:math:`\delta_{e,t,x}`      delta_e                                               Volume fraction occupied by emulsion region                           
:math:`\varepsilon_{t,x}`   voidage_average                                       Cross-sectional average voidage
:math:`\varepsilon_{e,t,x}` voidage_emulsion                                      Emulsion region voidage fraction 
:math:`\phi_{se,t,x}`       solid_emulsion_region.particle_porosity               Particle porosity of solid   
:math:`\gamma_{t,x}`        bubble_growth_coeff                                   Bubble growth coefficient
:math:`d_{bm,t,x}`          bubble_diameter_max                                   Maximum theoretical bubble diameter 
:math:`d_{b,t,x}`           bubble_diameter                                       Average bubble diameter
:math:`v_{g,t,x}`           velocity_superficial_gas                              Gas superficial velocity  
:math:`v_{ge,t,x}`          velocity_emulsion_gas                                 Emulsion region superficial gas velocity  
:math:`v_{br,t,x}`          velocity_bubble_rise                                  Bubble rise velocity
:math:`v_{b,t,x}`           velocity_bubble                                       Average bubble diameter  
:math:`K_{be,t,x,j}`        Kbe                                                   Bubble to emulsion gas mass transfer coefficient
:math:`K_{gbulkc,t,x,j}`    Kgbulk_c                                              Gas phase component bulk transfer rate
:math:`H_{be,t,x,j}`        Hbe                                                   Bubble to emulsion gas heat transfer coefficient
:math:`c_{p\_vap,b,t,x}`    cp_mol                                                Mixture mole heat capacity
:math:`h_{tc,t,x}`          htc_conv                                              Gas to solid convective heat transfer coefficient
:math:`\mu_{vap,ge,t,x}`    visc_d                                                Mixture dynamic viscosity
:math:`h_{t\_gs,t,x}`       ht_conv                                               Gas to solid convective enthalpy transfer
:math:`H_{gbulk,t,x}`       Hgbulk                                                Bulk gas heat transfer between bubble and emulsion
:math:`r_{hetero,ge,t,x,j}` gas_emulsion_hetero_rxn                               Gas emulsion heterogeneous rate reaction generation  
:math:`L_{b}`               bubble_region.length                                  
:math:`L_{ge}`              gas_emulsion_region.length                            
:math:`L_{se}`              solid_emulsion_region.length                          
:math:`A_{b,t,x}`           bubble_region.area                          
:math:`A_{ge,t,x}`          gas_emulsion_region.area                          
:math:`A_{se,t,x}`          solid_emulsion_region.area                          
:math:`\Delta P_{ge,t,x}`   gas_emulsion_region.deltaP                            pressure drop across gas emulsion region
:math:`\rho_{mass,se,t,x}`  solid_emulsion_region.properties.dens_mass_particle   solid particle mass density         
:math:`D_{vap,ge,t,x,j}`    gas_emulsion_region.properties.diffusion_comp         gas component diffusion in gas emulsion region       
:math:`C_{b,total,t,x}`     bubble_region.properties.dens_mol_vap                 gas mole density in the bubble region 
:math:`C_{g,t,x}`           average_gas_density                                   average gas density
:math:`C_{ge,total,t,x}`    gas_emulsion_region.properties.dens_mol_vap           gas mole density in the emulsion region       
:math:`M_{tr,b,t,x,p,j}`    bubble_region.mass_transfer_term              
:math:`M_{tr,ge,t,x,p,j}`   gas_emulsion_region.mass_transfer_term             
:math:`M_{tr,se,t,x,p,j}`   solid_emulsion_region.mass_transfer_term            
:math:`r_{ext,b,t,x,r}`     bubble_region.rate_reaction_extent          
:math:`r_{ext,ge,t,x,r}`    gas_emulsion_region.rate_reaction_extent          
:math:`r_{ext,se,t,x,r}`    solid_emulsion_region.rate_reaction_extent          
:math:`r_{b,t,x,r}`         bubble_region.reactions.reaction_rate                 
:math:`r_{ge,t,x,r}`        gas_emulsion_region.reactions.reaction_rate                 
:math:`r_{se,t,x,r}`        solid_emulsion_region.reactions.reaction_rate                  
:math:`k_{vap,b,t,x}`       bubble_region.properties.therm_cond                   bubble region thermal conductivity  
:math:`k_{vap,e,t,x}`       gas_emulsion_region.properties.therm_cond             gas emulsion region thermal conductivity        
:math:`T_{b,t,x}`           bubble_region.properties.temperature                   
:math:`T_{ge,t,x}`          gas_emulsion_region.properties.temperature                   
:math:`T_{se,t,x}`          solid_emulsion_region.properties.temperature                   
:math:`H_{tr, b, t,x}`      bubble_region.heat                                    bubble region heat transfer term
:math:`H_{tr, ge, t,x}`     gas_emulsion_region.heat                              gas emulsion region heat transfer term
:math:`H_{tr, se, t,x}`     solid_emulsion_region.heat                            solid emulsion region heat transfer term
:math:`F_{mol,b,t,x}`       bubble_region.properties.flow_mol                  
:math:`F_{mol,ge,t,x}`      gas_emulsion_region.properties.flow_mol            
:math:`y_{b,t,x,j}`         bubble_region.properties.mole_frac                  
:math:`y_{ge,t,x,j}`        gas_emulsion_region.properties.mole_frac            
:math:`x_{se,t,x,j}`        solid_emulsion_region.properties.mass_frac          
:math:`P_{ge,t,x}`          gas_emulsion_region.properties.pressure             
:math:`F_{mol,g,t,inlet}`   gas_inlet.flow_mol                      
:math:`y_{g,t,inlet,j}`     gas_inlet.mole_frac                     
:math:`P_{g,t,inlet}`       gas_inlet.pressure                     
:math:`T_{g,t,inlet}`       gas_inlet.temperature                   
:math:`H_{g,t,inlet}`       gas_inlet.enthalpy                      
:math:`F_{mass,s,t,inlet}`  solid_inlet.flow_mass
:math:`\phi_{t,inlet}`      solid_inlet.particle_porosity                            
:math:`x_{s,t,inlet}`       solid_inlet.mass_frac                    
:math:`T_{s,t,inlet}`       solid_inlet.temperature                   
:math:`H_{s,t,inlet}`       solid_inlet.enthalpy                      
:math:`F_{mass,se,t,x}`     solid_emulsion_region.properties.flow_mass                     
:math:`H_{b,t,x}`           bubble_region.properties.enthalpy                      
:math:`H_{ge,t,x}`          gas_emulsion_region.properties.enthalpy                      
:math:`H_{se,t,x}`          solid_emulsion_region.properties.enthalpy                      
:math:`F_{mol,g,t,outlet}`  gas_outlet.flow_mol                      
:math:`y_{g,t,outlet,j}`    gas_outlet.mole_frac                     
:math:`P_{g,t,outlet}`      gas_outlet.pressure                      
:math:`T_{g,t,outlet}`      gas_outlet.temperature                   
:math:`H_{g,t,outlet}`      gas_outlet.enthalpy                      
:math:`F_{mass,s,t,outlet}` solid_outlet.flow_mass
:math:`\phi_{t,outlet}`     solid_outlet.particle_porosity                            
:math:`x_{s,t,outlet}`      solid_outlet.mass_frac                     
:math:`T_{s,t,outlet}`      solid_outlet.temperature                   
:math:`H_{s,t,outlet}`      solid_outlet.mass_enthalpy                      
:math:`v_{mf,se}`           solid_emulsion_region.properties.velocity_mf          velocity at minimum fluidization
:math:`\varepsilon_{mf,se}` solid_emulsion_region.properties.voidage_mf           voidage at minimum fluidization
:math:`K_{d}`               Kd                                                    bulk gas permeation coefficient
:math:`d_{p,se}`            solid_emulsion_region.properties._params.particle_dia                  
:math:`\Delta P_{or}`       deltaP_orifice                                        Pressure drop across orifice
=========================== ===================================================== ====================================================

Parameters
----------

List of parameters in the BFBR model:

============================ ============================= =========================================================================
Parameter                    Name                          Notes
============================ ============================= =========================================================================
:math:`s_{se,j,r}`           rate_reaction_stoichiometry   Reference to solid_emulsion_region.reactions.rate_reaction_stoichiometry 
============================ ============================= =========================================================================

Subscripts
----------

List of subscripts in the BFBR model:

====================== ====================
Subscript              Name                 
====================== ====================
:math:`b`              bubble region
:math:`e`              emulsion region 
:math:`g`              gas phase
:math:`ge`             gas_emulsion       
:math:`j`              component
:math:`p`              phase
:math:`s`              solid phase          
:math:`se`             solid_emulsion  
:math:`r`              reaction
:math:`t`              time 
:math:`x`              length
====================== ====================
                         
Initialization
--------------

The initialization method for this model will save the current state of the model before commencing initialization and reloads it afterwards. 
The state of the model will be the same after initialization, only the initial guesses for unfixed variables will be changed. 

The model allows for the passing of a dictionary of values of the state variables of the gas and solid phases that can be used as initial 
guesses for the state variables throughout the time and spatial domains of the model. This is optional but recommended. 
A typical guess could be values of the gas and solid inlet port variables at time t=0. 

The model initialization proceeds through a sequential hierarchical method where the model equations are deactivated at the start of the 
initialization routine, and the complexity of the model is built up through activation and solution of various sub-model blocks and equations 
at each initialization step. At each step the model variables are updated to better guesses obtained from the model solution at that step.

The initialization routine proceeds as follows:

* Step 1. Initialize the thermo-physical and transport properties model blocks
* Step 2. Initialize geometric constraints
* Step 3. Initialize the hydrodynamic properties
* Step 4a. Initialize pressure drop and mass balances without reactions
* Step 4b. Initialize mass balances with reactions
* Step 5. Initialize energy balances 

BFBR Class
----------

.. module:: idaes.gas_solid_contactors.unit_models.bubbling_fluidized_bed

.. autoclass:: BubblingFluidizedBed
  :members:

BFBRData Class
--------------

.. autoclass:: BubblingFluidizedBedData
  :members:
