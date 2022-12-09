Power Plant Costing Using NETL Baseline Reports
===============================================

.. contents:: Contents 
    :depth: 4

Introduction
------------

.. note:: This document outlines an example power plant costing library to support IDAES power generation unit models. The power plant costing method is available for most of the unit operations in power plants (Boiler, Feed Water Heaters, Compressor, Turbine, Condenser, etc.). For a general overview of costing capabilities in IDAES, see the :ref:`IDAES Process Costing Framework<reference_guides/core/costing/costing_framework:IDAES Process Costing Framework>`.

A capital cost methodology is developed in this module, both bare and erected cost and total plant cost are calculated based on costing correlations.

The Power Plant Costing Library contains 5 main costing functions: `get_PP_costing`, `get_SCO2_unit_cost`, `get_ASU_cost`, `get_fixed_OM_costs`, and `get_variable_OM_costs`.
The first function (`get_PP_costing`) can be called to include cost correlations for equipment typically used in simulation of 7 technologies: 

1. Supercritical pulverized coal plants (SCPC),
2. Subcritical pulverized coal plants,
3. Two-stage IGCC,
4. Single-stage IGCC,
5. Single-stage dry-feed IGCC,
6. natural gas air-fired plant (NGCC),
7. Advanced ultra-supercritical PC (AUSC).

Similarly, `get_sCO2_unit_cost` can be called to include cost correlations for equipment in supercritical CO2 power cycle plants and the method `get_ASU_cost` calls costing for equipment in air separation units. The method `costing_initialization` exists to initialize pp, sCO2 and ASU costing blocks.

The methods `get_fixed_OM_costs` and `get_variable_OM_costs` calculate operating and maintenance costs for fixed (operating labor, maintenance labor, admin/support labor, property taxes and insurance, maintenance materials) and variable (fuel, consumable and waste disposal costs for electricity production indexed by time) costs, respectively. The methods `initialize_fixed_OM_costs` and `initialize_variable_OM_costs` exist to initialize fixed and variable O&M costing blocks. To build all operating and maintanence costs sequentially, users may call the `build_process_costs` method.

Several reporting methods exist (`report`, `display_total_plant_costs`, `display_bare_erected_costs`, `display_equipment_costs`, `get_total_TPC`, `display_flowsheet_cost`) as well as a check on supercritical CO2 bounds (`check_sCO2_costing_bounds`).

Finally, `custom_power_plant_currency_units` appends custom units for power plant scenarios to the existing registered Pyomo currency units dictionary.

Details are given for each method later in this documentation; however, there are many similarities between methods as described below:

Costing sub-blocks
^^^^^^^^^^^^^^^^^^

In general, when `get_PP_costing` or `get_sCO2_unit_cost` is called on an instance of a unit model, a new sub-block is created 
on that unit named `costing` (i.e. `flowsheet.unit.costing`). All variables and constraints related to costing will be 
constructed within this new block (see detailed documentation for each unit for details on these variables and constraints).


Capital Cost Stages
^^^^^^^^^^^^^^^^^^^

There are multiple stages of capital cost, the lowest stage is the equipment cost which only includes
the cost of manufacturing the equipment. The next stage is the bare erected cost (BEC) which includes
the equipment cost and the cost of material and labor for installation. The final stage is the total
plant cost (TPC) which includes the BEC plus the engineering fee, process contingency,
and project contingency, all of which are typically estimated as a percentage of BEC.

.. math:: bare\_erected\_cost = equipment\_cost*(1 + material\_cost + labor\_cost)

.. math:: total\_plant\_cost = bare\_erected\_cost*(1 + eng\_fee + process\_contingency + project\_contingency)

.. note:: The equations above assume the additional costs (eng_fee or process and project contingencies) are given as percentages of BEC and TPC.

All costing methods calculate the bare erected and total plant costs. The sCO2 library is currently the only one
that includes an equipment cost. 

Dollar Year scaling
^^^^^^^^^^^^^^^^^^^

The value of money decreases over time due to inflation and missed investment opportunity.
Thus, all costs must be normalized to the same dollar year to be compared on a consistent basis.
This is done using a CE index and the following formula:

.. math:: bare\_erected\_cost = bare\_erected\_cost_{base\_year}*(CE\_index/CE\_index_{base\_year})

In the costing functions this equation is built into the constraint for the lowest level capital cost in the selected method.

Table 1. Base years of costing modules

=========================== ====================== 
Module                      Base Year
=========================== ======================
Power Plant Costing         2018
sCO2 Costing                2017     
ASU                         2011

=========================== ====================== 

When a `costing` block is created 
on the flowsheet object (i.e. `flowsheet.costing`), the methods automatically build any global parameters relating to costing under this block. The most 
common of these paramters is the CE index parameter. The CE index will be set to the base year of the method called, set by the argument `CE_index_year` which is allowed in most methods in this module.

Power Plant Costing Module
--------------------------

A costing module has been developed based on the capital cost scaling methodology from 
NETL's Bituminous Baseline Report Rev 4 [1]. It provides costing correlations for common 
variants of pulverized coal (PC), integrated gassification combined cycle (IGCC), and 
natural gas combined cycle (NGCC) power generation technologies. Users should refer to 
reference [2] for details of the costing correlations, however, a summary is provided below.


The module breaks down the cost of a power plant into separate accounts for each system 
within the plant. The accounts are scaled based on a process parameter that determines
the size of the equipment needed. The cost of the account is computed based on the scaled parameter,
reference parameter, reference cost, and scaling
exponent determined by NETL in [1]. This equation is similar to a six tenth factor approach, 
however, the exponents have been trained using several vendor quotes.

.. math:: scaled\_cost = reference\_cost*(\frac{scaled\_param}{reference\_param})^\alpha

where:

* scaled_cost - the cost of the system in Million dollars
* reference_cost - the cost of the reference system in thousands of dollars
* scaled_param - the value of the system's process parameter
* reference_param - the value of the reference system's process parameter
* alpha - scaling exponent

.. note:: The capital cost scaling equation can be applied to any capital cost stage. In the power plant costing library it is applied to the bare erected cost, while in the sCO2 library it is applied to the equipment cost.

The Power Plant costing method has the following arguments:

* blk : an existing unit model or Pyomo Block
* cost_accounts : A list of accounts or a string containing the name of a pre-named account. If the input is a list all accounts must share the same process parameter. Pre-named accounts are listed below.
* scaled_param : The Pyomo Variable representing the accounts' scaled parameter
* tech : The technology to cost, as different technologies have different accounts
* ccs : which reference parameter to use, as some accounts are costed using two different reference parameters; defaults to "B", and "A" is also a valid option
* CE_index_year : Chemical Engineering Cost Index base year, defaults to 2018; calling the registered Pyomo currency units dictionary of plant cost index values will allow conversion between base years within the flowsheet
* additional_costing_params : option to add a costing parameter dictionary to supplement existing account data
* use_additional_costing_params : True/False flag whether IDAES should use new data when additonal account names match existing account names, or fail with a useful error message; defaults to False, meaning pre-installed accounts in BB_costing_params and generic_ccs_costing will never be bypassed or overwritten by duplicated additional accounts unless this flag is set to True on the unit costing block argument level

 1. Supercritical PC,
 2. Subcritical PC, 
 3. Two-stage, slurry-feed IGCC 
 4. Single-stage, slurry-feed IGCC
 5. Single-stage, dry-feed IGCC,
 6. Natural Gas Combined Cycle (NGCC), 
 7. Advanced Ultrasupercritical PC


Many accounts scale using the same process parameter. For convenience the user is allowed to enter accounts as a list instead
of having to cost each account individually. If the accounts in the list do not use the same process parameter an error will be raised.

It is recognized that many users will be unfamiliar with the accounts in the Bituminous Baseline.
For this reason the cost_accounts argument will also accept a string with the name of a pre-named
account. Pre-nammed accounts aggregate the relevant accounts for certain systems. The pre-named
accounts for each technology can be found in the tables below.

Table 2. Pre-named Accounts for PC technologies

=========================== ============================ ============================ ==========
Pre-named Account           Accounts Included            Process Parameter            Units      
=========================== ============================ ============================ ==========
Coal Handling               1.1, 1.2, 1.3, 1.4, 1.9a     Coal Feed Rate               lb/hr           
Sorbent Handling            1.5, 1.6, 1.7, 1.8, 1.9b     Limestone Feed Rate          lb/hr  
Coal Feed                   2.1, 2.2, 2.9a               Coal Feed Rate               lb/hr     
Sorbent Feed                2.5, 2.6, 2.9b               Limestone Feed Rate          lb/hr
Feedwater System            3.1, 3.3                     HP BFW Flow Rate             lb/hr 
PC Boiler                   4.9                          HP BFW Flow Rate             lb/hr
Steam Turbine               8.1                          Steam Turbine Power          kW
Condenser                   8.3                          Condenser Duty               MMBtu/hr
Cooling Tower               9.1                          Cooling Tower Duty           MMBtu/hr
Circulating Water System    9.2, 9.3, 9.4, 9.6, 9.7      Circulating Water Flow Rate  gpm
Ash Handling                10.6, 10.7, 10.9             Total Ash Flow               lb/hr
=========================== ============================ ============================ ==========

Table 3. Pre-named Accounts for IGCC technologies

=========================== ========================================= ============================ ==========
Pre-named Account           Accounts Included                         Process Parameter            Units      
=========================== ========================================= ============================ ==========
Coal Handling               1.1, 1.2, 1.3, 1.4, 1.9                   Coal Feed Rate               lb/hr           
Coal Feed                   2.1, 2.2, 2.9                             Coal Feed Rate               lb/hr     
Feedwater System            3.1, 3.3                                  HP BFW Flow Rate             lb/hr 
Gasifier                    4.1                                       Coal Feed Rate               lb/hr
Syngas Cooler               4.2                                       Syngas Cooler Duty           MMBtu/hr
ASU                         4.3a                                      Oxygen Production            tpd
ASU Oxidant Compression     4.3b                                      Main Air Compressor Power    kW
Combustion Turbine          6.1, 6.3                                  Syngas Flowrate              lb/hr
Syngas Expander             6.2                                       Syngas Flowrate              lb/hr
HRSG                        7.1, 7.2                                  HRSG Duty                    MMBtu/hr
Steam Turbine               8.1                                       Steam Turbine Power          MW
Condenser                   8.3                                       Condenser Duty               MMBtu/hr
Cooling Tower               9.1                                       Cooling Tower Duty           MMBtu/hr
Circulating Water System    9.2, 9.3, 9.4, 9.6, 9.7                   Circulating Water Flow Rate  gpm
Slag Handling               10.1, 10.2, 10.3, 10.6, 10.7, 10.8, 10.9  Slag Production              lb/hr
=========================== ========================================= ============================ ==========

Table 4. Pre-named Accounts for NGCC technologies

=========================== ============================ ============================ ==========
Pre-named Account           Accounts Included            Process Parameter            Units      
=========================== ============================ ============================ ==========
Feedwater System            3.1, 3.3                     HP BFW Flow Rate             lb/hr 
Combustion Turbine          6.1, 6.3                     Fuel Gas Flow                lb/hr  
HRSG                        7.1, 7.2                     HRSG Duty                    MMBtu/hr     
Steam Turbine               8.1                          Steam Turbine Power          kW
Condenser                   8.3                          Condenser Duty               MMBtu/hr
Cooling Tower               9.1                          Cooling Tower Duty           MMBtu/hr
Circulating Water System    9.2, 9.3, 9.4, 9.6, 9.7      Circulating Water Flow Rate  gpm
=========================== ============================ ============================ ==========

The library has a 7th technology of AUSC. These operate at higher temperatures that traditional 
PC plants. The library contains modified correlation for the PC boiler, steam turbine, and steam piping
to reflect the use of high temperature materials.

Table 5. Pre-named Accounts for AUSC technologies

=========================== ============================ ============================ ==========
Pre-named Account           Accounts Included            Process Parameter            Units      
=========================== ============================ ============================ ==========
PC Boiler                   4.9                          HP BFW Flow Rate             lb/hr 
Steam Turbine               8.1                          Steam Turbine Power          kW
Steam Piping                8.4                          HP BFW Flow Rate             lb/hr
=========================== ============================ ============================ ==========


A call to get_PP_costing creates two variables and two constraints for each account in the list.
The variables are bare_erected_cost and total_plant_cost. Both variables are indexed
by the account number in string format. The function makes a new block called self.costing where
all variables and parameters associated with costing are stored.

.. note:: The bare_erected_cost and total_plant_cost are in Million dollars.



Example
^^^^^^^
Below is an example of how to add cost correlations to a flowsheet including a heat exchanger using the IDAES power plant costing module:


.. code:: python

    from pyomo.environ import (ConcreteModel, SolverFactory)
    from idaes.core import FlowsheetBlock
    from idaes.models.unit_models.heat_exchanger import \
        (HeatExchanger, HeatExchangerFlowPattern)
    from idaes.models.properties import iapws95
    from idaes.models_extra.power_generation.costing.power_plant_capcost import \
        QGESSCostingData
    
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    
    m.fs.unit = HeatExchanger(
        shell={"property_package": m.fs.properties},
        tube={"property_package": m.fs.properties},
        flow_pattern=HeatExchangerFlowPattern.countercurrent)
    # set inputs
    m.fs.unit.shell_inlet.flow_mol[0].fix(100*pyunits.mol/pyunits.s)
    m.fs.unit.shell_inlet.enth_mol[0].fix(3500*pyunits.mol/pyunits.s)
    m.fs.unit.shell_inlet.pressure[0].fix(101325*pyunits.Pa)
    
    m.fs.unit.tube_inlet.flow_mol[0].fix(100*pyunits.mol/pyunits.s)
    m.fs.unit.tube_inlet.enth_mol[0].fix(4000*pyunits.mol/pyunits.s)
    m.fs.unit.tube_inlet.pressure[0].fix(101325*pyunits.Pa)
    
    m.fs.unit.area.fix(1000*pyunits.m**2)
    m.fs.unit.overall_heat_transfer_coefficient.fix(100*pyunits.W/pyunits.m**2/pyunits.K)
    
    m.fs.unit.initialize()
    
    m.fs.unit.duty_MMbtu = pyo.Var(
        m.fs.time,
        initialize=1e5,
        doc="Condenser duty in MMbtu/hr",
        units=pyunits.MBtu/pyunits.hr)
    
    @m.fs.unit.Constraint(m.fs.time)
    def duty_conversion(b, t):
        conv_fact = 3.412/1e6 
        return b.duty_MMbtu[t] == conv_fact*b.heat_duty[t]
    
    hx_accounts = ["6.7.ccs"]
    m.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": hx_accounts,
            "scaled_param": m.fs.unit.duty_MMBtu,
            "tech": 1,
            "ccs": "A",
            "CE_index_year": "2013",
        },
    )

    # initialize costing equations
    QGESSCostingData.costing_initialization(m.fs.costing)
    
    opt = SolverFactory('ipopt')
    opt.options = {'tol': 1e-6, 'max_iter': 50}
    results = opt.solve(m, tee=True)
    
    QGESS.display_flowsheet_cost(m.fs.costing)


Supercritical CO2 Costing Module
--------------------------------

The sCO2 costing function, besides including the capital cost and engineering of the equipment, it can cost penalty due to the high temperature and pressure equipment required for sCO2 PC plants.
The function has has five arguments, self, equipment, scaled_param, temp_C, and n_equip.

* self : an existing unit model or Pyomo Block
* equipment : The type of equipment to be costed, see table 6
* scaled_param : The Pyomo Variable representing the component's scaled parameter
* temp_C : The Pyomo Variable representing the hottest temperature of the piece of equiment being costed. Some pieces of equipment do not have a temperature associated with them, so the default argument is None. This variable must have units of Celsius (Pyomo label `pyunits.C`).
* n_equip : The number of pieces of equipment to be costed. The function will evenly divide the scaled parameter between the number passed.
* CE_index_year : Chemical Engineering Cost Index base year, defaults to 2018; calling the registered Pyomo currency units dictionary of plant cost index values will allow conversion between base years within the flowsheet
* custom_accounts : Additional accounts to cost

The equipment cost is calculated using the following two equations. A temperature correction factor account for advanced materials needed at high temperatures.

.. math:: equipment\_cost = reference\_cost*(scaled\_parameter)^\alpha * temperature\_factor

.. math:: temperature\_factor = 1 + c*(T - T_{bp}) + d*(T - T_{bp})^2 & : if T \geq T_{bp}\\ (if  T > 550, otherwise  temperature\_factor = 1)
    
.. math:: T_{bp} = 550 C

The bare erected and total plant costs are calculated as shown in the introduction.
There are currently no estimates for the total plant cost components, so bare erected cost will be the same as total plant cost for now.

Five variables and constraints are created within the costing block. Three are for the equipment, bare erected, and total plant costs. One is for the temperature correction factor.
The last one is for the scaled parameter divided by n_equip.

Table 6. sCO2 Costing Library Components

=========================== ================= ============== 
Component                   Scaling Parameter Units               
=========================== ================= ============== 
Coal-fired heaters          :math:`Q`         :math:`MW_{th}`
Natural gas-fired heaters   :math:`Q`         :math:`MW_{th}`
Recuperators                :math:`UA`        :math:`W/K`
Direct air coolers          :math:`UA`        :math:`W/K`
Radial turbines             :math:`W_{sh}`    :math:`MW_{sh}`
Axial turbines              :math:`W_{sh}`    :math:`MW_{sh}`
IG centrifugal compressors  :math:`W_{sh}`    :math:`MW_{sh}`       
Barrel type compressors     :math:`V_{in}`    :math:`m^3/s`     
Gearboxes                   :math:`W_{e}`     :math:`MW_{sh}`   
Generators                  :math:`W_{e}`     :math:`MW_{e}`   
Explosion proof motors      :math:`W_{e}`     :math:`MW_{e}`
Synchronous motors          :math:`W_{e}`     :math:`MW_{e}`
Open drip-proof motors      :math:`W_{e}`     :math:`MW_{e}`
=========================== ================= ==============

Other Costing Modules
---------------------

Air Separation Unit Costing Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ASU costing function calculates total plant cost in the exact same way as the get_PP_costing function.
get_ASU_cost takes two arguments: self, and scaled_param. 

* self - a Pyomo Block or unit model
* scaled_param - The scaled parameter. For the ASU it is the oxygen flowrate in units of tons per day.
* CE_index_year : Chemical Engineering Cost Index base year, defaults to 2018; calling the registered Pyomo currency units dictionary of plant cost index values will allow conversion between base years within the flowsheet

The bare erected and total plant costs are calculated as shown in the introduction.
There are currently no estimates for the total plant cost components, so bare erected cost will be the same as total plant cost for now.

Fixed Operating & Maintenance Costs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Fixed O&M costing function adds constraints to estimate the labor, maintenance, support, taxes and other fixed costs. The method takes the following arguments: 

* b : Pyomo concrete model or flowsheet block
* net_power: production rate of the plant in MW, if provided will enable additional production-related cost calculations but not required to use method
* nameplate_capacity : rated plant output in MW, defaults to 650
* capacity_factor: factor adjusting variable costs for plant capacity; defaults to 85%.
* labor_rate : hourly rate of plant operators in project dollar year, defaults to 38.50
* labor_burden : a percentage multiplier used to estimate non-salary labor expenses, defaults to 30
* operators_per_shift : average number of operators per shift, defaults to 6
* tech : int 1-7 representing the catagories in get_PP_costing, used to determine maintenance costs, defaults to 1
* fixed_TPC : The TPC in $MM that will be used to determine fixed O&M costs. If the value is None, the function will try to use the TPC calculated from the individual units.
* CE_index_year : Chemical Engineering Cost Index base year, defaults to 2018; calling the registered Pyomo currency units dictionary of plant cost index values will allow conversion between base years within the flowsheet

The following maintenance cost percentages are assumed:

========== ============================= =========================
Technology Maintenance Labor / TPC Split Maintenance Labor Percent
========== ============================= =========================
1          0.4                           0.016                    
2          0.4                           0.016                    
3          0.35                          0.03                     
4          0.35                          0.03                     
5          0.35                          0.03                     
6          0.4                           0.019                    
7          0.4                           0.016                    
========== ============================= =========================

When this method is called, the following equations are added to the flowsheet as constraints:

.. math:: annual\_operating\_labor\_cost = operators\_per\_shift * labor\_rate * (1 + 0.01 * labor\_burden) * 8760

.. math:: maintenance\_labor\_cost = TPC * maintenance\_labor\_TPC\_split * maintenance\_labor\_percent

.. math:: admin\_and\_support\_labor\_cost = 0.25 * (annual\_operating\_labor\_cost + maintenance\_labor\_cost)

.. math:: property\_taxes\_and\_insurance = 0.02 * TPC

.. math:: total\_fixed\_OM\_cost = annual\_operating\_labor\_cost + maintenance\_labor\_cost + \\ admin\_and\_support\_labor\_cost + property\_taxes\_and\_insurance + other\_fixed\_costs

.. math:: maintenance\_material\_cost = \frac {TPC * maintenance\_material\_TPC\_split * maintenance\_material\_percent}{0.85 * nameplate\_capacity * 8760}

where 8760 is the number of operating hours per year, 0.25, 0.02 and 0.85 are assumed cost ratio coefficients, and the variable `other_fixed_costs` exists to allow for unaccounted fixed costs (default = 0).

Variable Operating & Maintenance Costs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Variable O&M costing function adds constraints to calculate correlations associated with fuel, consumable and waste disposal costs. The function may be used to calculate variable costs of producing electricity in $/MWh. The method takes the following arguments: 

* b: pyomo flowsheet block
* resources : a list of strings for the resorces to be costed
* rates : a list of pyomo vars for resource consumption rates
* prices : a dict of resource prices to be added to the premade dictionary
* CE_index_year : Chemical Engineering Cost Index base year, defaults to 2018; calling the registered Pyomo currency units dictionary of plant cost index values will allow conversion between base years within the flowsheet
* capacity_factor: factor adjusting variable costs for plant capacity; defaults to 85%.

The following default prices are assumed:

============================ ======= ==========================
Stream                       Price   Units                     
============================ ======= ==========================
natural_gas                  4.42    USD (2018) per Million BTU
coal                         51.96   USD (2018) per ton        
water                        1.90E-3 USD (2018) per gallon     
water_treatment_chemicals    550     USD (2018) per ton        
ammonia                      300     USD (2018) per ton        
SCR_catalyst                 150     USD (2018) per cubic foot 
triethylene_glycol           6.80    USD (2018) per gallon     
SCR_catalyst_waste           2.50    USD (2018) per cubic foot 
triethylene_glycol_waste     0.35    USD (2018) per gallon     
amine_purification_unit      38      USD (2018) per ton        
thermal_reclaimer_unit_waste 38      USD (2018) per ton        
============================ ======= ==========================

When this method is called, the following equations are added to the flowsheet as constraints:

.. math:: variable\_operating\_costs_{t, r} = resource\_prices_r * resource\_rates_{r, t} * 365 * 0.85

.. math:: total\_variable\_OM\_costs_t = \Sigma_r (variable\_operating\_costs_{t, r}) + maintenance\_material\_cost + other\_variable\_costs_t

where variables are indexed by resource `r` and time `t`, 0.85 is an assumed price ratio coefficient, and the variable `other_variable_costs` exists to allow for unaccounted variable costs (default = 0).

Build Process Costs
^^^^^^^^^^^^^^^^^^^

Users may quickly build all required process costs by calling the global method `m.fs.costing.build_process_costs`, which internally determines which fixed and variable cost quantities should be considered and/or calculated. The method takes the following arguments:

* total_plant_cost: The TPC in $MM that will be used to determine fixed O&M, costs. If the value is None, the function will try to use the TPC calculated from the individual units. This quantity should be a Pyomo Var or Param that will contain the TPC value.
* nameplate_capacity: rated plant output in MW
* labor_rate: hourly rate of plant operators in project dollar year
* labor_burden: a percentage multiplier used to estimate non-salary labor expenses
* operators_per_shift: average number of operators per shift
* tech: int 1-7 representing the catagories in get_PP_costing, used to determine maintenance costs
* land_cost: Expression, Var or Param to calculate land costs
* net_power: actual plant output in MW, only required if calculating variable costs
* resources: list setting resources to cost
* rates: list setting flow rates of resources
* prices: list setting prices of resources
* fixed_OM: True/False flag for calculating fixed O&M costs
* variable_OM: True/False flag for calculating variable O&M costs
* fuel: string setting fuel type for fuel costs
* chemicals: string setting chemicals type for chemicals costs
* waste: string setting waste type for waste costs
* tonne_CO2_capture: Var or value to use for tonnes of CO2 capture in one year
* CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars

Utility Functions
-----------------

Initialize Costing
^^^^^^^^^^^^^^^^^^

The `initialize_fixed_OM_costs` will initialize all fixed cost variable and constraint in the costing block. The `initialize_variable_OM_costs` does the same.

The `costing_initialization` function will initialize all the variable within every costing block in the flowsheet.
It takes one argument, the flowsheet object. It should be called after all the calls to 'get costing' functions are 
completed.

The function iterates through the flowsheet looking for costing blocks and calculates variables from constraints.

Total Flowsheet Cost Constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For optimization, a constraint summing all the total plant costs is required.
Calling build_flowsheet_cost_constraint(m) creates a variable named m.fs.costing.flowsheet_cost 
and builds the required constraint at the flowsheet level.

.. note:: The costing libraries can be used for simulation or optimization. For simulation, costing constraints can be built and solved after the flowsheet has been solved. For optimization, the costing constraints will need to be solved with the flowsheet.


Display Total Flowsheet Cost
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calling display_flowsheet_cost(m.fs.costing) will print the value of m.fs.costing.flowsheet_cost.


Display Individual Costs
^^^^^^^^^^^^^^^^^^^^^^^^

There are three functions for displaying individual costs.

* display_total_plant_costs(m.fs.costing)
* display_bare_erected_costs(m.fs.costing)
* display_equipment_costs(m.fs.costing)

Each one prints out a list of the costed blocks and the cost level of the function chosen.
The functions should be called after solving the model.

Calling the `report(m.fs.costing)` method will print specific total costs (e.g. total TPC, cost of electricity, annualized cost) if they exist and are calculated by the methods.

Checking Bounds
^^^^^^^^^^^^^^^

Currently, only the sCO2 module has support for checking bounds.

All costing methods have a range, outside of which the correlations become inaccurate.
Calling check_sCO2_costing_bounds(m.fs.costing) will display which components are within the proper range
and which are outside it. It should be called after the model is solved.

Adding Custom Currency Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The method `custom_power_plant_currency_units` allows the addition of custom currency units, such as specific months (e.g. USD_2008_Nov for the ASU reference costs) to the existing currency units dictionary containing only aggregates for each cost year.

References
----------

1. DOE/NETL-2015/1723 Cost and Performance Baseline for Fossil Energy Plants Volume 1a: Bituminus Coal (PC) and Natural Gas to Electricity Revision 3 and 4
2. DOE/NETL-341/013113 Quality Guidelines for Energy System Studies Capital Cost Scaling Methodology
3. NETL_PUB_21490 Techno-economic Evaluation of Utility-Scale Power Plants Based on the Indirect sCO2 Brayton Cycle. Charles White, David Gray, John Plunkett, Walter Shelton, Nathan Weiland, Travis Shultz. September 25, 2017
4. SCO2 Power Cycle Component Cost Correlations from DOE Data Spanning Multiple Scales and Applications. Nathan Weiland, Blake Lance, Sandeep Pidaparti. Proceedings of ASME Turbo Expo 2019: Turbomachinery Technical Conference and Exposition GT2019. June 17-21, 2019, Phoenix, Arizona, USA