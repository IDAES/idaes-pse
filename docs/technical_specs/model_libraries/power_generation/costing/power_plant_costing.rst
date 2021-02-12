Power Plant Costing Library
===========================

.. contents:: Contents 
    :depth: 4

Introduction
------------

.. note:: The power plant costing method is available for most of the unit operations in power plants (Boiler, Feed Water Heaters, Compressor, Turbine, Condenser, etc.).

A capital cost methodology is developed in this module, both bare and erected cost and total plant cost are calculated based on costing correlations. 
The Power Plant Costing Library contains two main costing functions `get_PP_costing` and `get_SCO2_unit_cost`.
The first function (get_PP_costing) can be called to include cost correlations for equipment typically used in simulation of 7 technologies: 

1. Supercritical pulverized coal plants (SCPC),
2. Subcritical pulverized coal plants,
3. Two-stage IGCC,
4. Single-stage IGCC,
5. Single-stage dry-feed IGCC,
6. natural gas air-fired plant (NGCC),
7. Advanced ultra-supercritical PC (AUSC).

Similarly, `get_sCO2_unit_cost` can be called to include cost correlations for equipment in supercritical CO2 power cycle plants.

Details are given for each method later in this documentation, 
however, there are many similarities between methods as discribed below:

Costing sub-blocks
^^^^^^^^^^^^^^^^^^

In general, when `get_PP_costing or get_sCO2_unit_cost` is called on an instance of a unit model, a new sub-block is created 
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

Dollar year scaling
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

The first time a 'get costing' function is called for a unit operation within a flowsheet, an additional `costing` block is created 
on the flowsheet object (i.e. `flowsheet.costing`) in order to hold any global parameters relating to costing. The most 
common of these paramters is the CE index parameter. The CE index will be set to the base year of the method called.

.. note:: The global paramters are created when the first instance of `get_costing` is called and use the values provided there for initialization. Subsequent `get_costing` calls use the existing paramters, and do not change the initialized values. i.e. any "year" argument provided to a `get_costing` call after the first will be ignored.

To manually set the dollar year the user must call m.fs.get_costing(year=2019) before any calls to a 'get costing' function are made.



Power Plant Costing Module
--------------------------

A default costing module has been developed based on the capital cost scaling methodology from 
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

* sacaled_cost - the cost of the system in Million dollars
* reference_cost - the cost of the reference system in thousands of dollars
* scaled_param - the value of the system's process parameter
* reference_param - the value of the reference system's process parameter
* alpha - scaling exponent

.. note:: The capital cost scaling equation can be applied to any capital cost stage. In the power plant costing library it is applied to the bare erected cost, while in the sCO2 library it is applied to the equipment cost.

The Power Plant costing method has five arguments, self, cost_accounts, scaled_param, units, and tech.

* self : an existing unit model or Pyomo Block
* cost_accounts : A list of accounts or a string containing the name of a pre-named account. If the input is a list all accounts must share the same process parameter. Pre-named accounts are listed below.
* scaled_param : The Pyomo Variable representing the accounts' scaled parameter
* tech : The technology to cost, different technologies have different accounts.

 1. Supercritical PC,
 2. Subcritical PC, 
 3. Two-stage, slurry-feed IGCC 
 4. Single-stage, slurry-feed IGCC
 5. Single-stage, dry-feed IGCC,
 6. Natural Gas Combined Cycle (NGCC), 
 7. Advanced Ultrasupercritical PC

* units : The user must pass a string with the units the scaled_param is in. It serves as a check to make sure the costing method is being used correctly.

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
Below is a simple example of how to add cost correlations to a flowsheet including a heat exchanger using the default IDAES costing module.


.. code:: python

    from pyomo.environ import (ConcreteModel, SolverFactory)
    from idaes.core import FlowsheetBlock
    from idaes.generic_models.unit_models.heat_exchanger import \
        (HeatExchanger, HeatExchangerFlowPattern)
    from idaes.generic_models.properties import iapws95
    from idaes.power_generation.costing.power_plant_costing import \
        (get_PP_costing, initialize_costing, display_total_plant_costs,
         display_flowsheet_cost)
    
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    
    m.fs.unit = HeatExchanger(default={
                "shell": {"property_package": m.fs.properties},
                "tube": {"property_package": m.fs.properties},
                "flow_pattern": HeatExchangerFlowPattern.countercurrent})
    # set inputs
    m.fs.unit.shell_inlet.flow_mol[0].fix(100)     # mol/s
    m.fs.unit.shell_inlet.enth_mol[0].fix(3500)    # j/s
    m.fs.unit.shell_inlet.pressure[0].fix(101325)  # Pa 
    
    m.fs.unit.tube_inlet.flow_mol[0].fix(100)
    m.fs.unit.tube_inlet.enth_mol[0].fix(4000)
    m.fs.unit.tube_inlet.pressure[0].fix(101325.0)
    
    m.fs.unit.area.fix(1000)  # m2
    m.fs.unit.overall_heat_transfer_coefficient.fix(100)  # W/m2K
    
    m.fs.unit.initialize()
    
    m.fs.unit.duty_MMbtu = pyo.Var(
        m.fs.time,
        initialize=1e5,
        doc="Condenser duty in MMbtu/hr")
    
    @m.fs.unit.Constraint(m.fs.time)
    def duty_conversion(b, t):
        conv_fact = 3.412/1e6 
        return b.duty_MMbtu[t] == conv_fact*b.heat_duty[t]
    
    get_PP_costing(m.fs.unit, 'Condenser', m.fs.unit.duty_MMbtu, 'MMBtu/hr', 1)
    # initialize costing equations
    initialize_costing(fs)
    
    opt = SolverFactory('ipopt')
    opt.options = {'tol': 1e-6, 'max_iter': 50}
    results = opt.solve(m, tee=True)
    
    display_total_plant_costs(fs)
    display_flowsheet_cost(fs)


Supercritical CO2 Costing Module
--------------------------------

The sCO2 costing function, besides including the capital cost and engineering of the equipment, it can cost penalty due to the high temperature and pressure equipment required for sCO2 PC plants.
The function has has five arguments, self, equipment, scaled_param, temp_C, and n_equip.

* self : an existing unit model or Pyomo Block
* equipment : The type of equipment to be costed, see table 6
* scaled_param : The Pyomo Variable representing the component's scaled parameter
* temp_C : The Pyomo Variable representing the hottest temperature of the piece of equiment being costed. Some pieces of equipment do not have a temperature associated with them, so the default argument is None.
* n_equip : The number of pieces of equipment to be costed. The function will evenly divide the scaled parameter between the number passed.

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

Air Separation Unit
^^^^^^^^^^^^^^^^^^^

The ASU costing function calculates total plant cost in the exact same way as the get_PP_costing function.
get_ASU_cost takes two arguments: self, and scaled_param. 

* self - a Pyomo Block or unit model
* scaled_param - The scaled parameter. For the ASU it is the oxygen flowrate in units of tons per day.

Utility Functions
-----------------

Initialize Costing
^^^^^^^^^^^^^^^^^^

The costing_initialization function will initialize all the variable within every costing block in the flowsheet.
It takes one argument, the flowsheet object. It should be called after all the calls to 'get costing' functions are 
completed.

The function iterates through the flowsheet looking for costing blocks and calculates variables from constraints.

Total Flowsheet Cost Constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For optimization, a constraint summing all the total plant costs is required.
Calling build_flowsheet_cost_constraint(m) creates a variable named m.fs.flowsheet_cost 
and builds the required constraint at the flowsheet level.

.. note:: The costing libraries can be used for simulation or optimization. For simulation, costing constraints can be built and solved after the flowsheet has been solved. For optimization, the costing constraints will need to be solved with the flowsheet.


Display Total Flowsheet Cost
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calling display_flowsheet_cost(m) will print the value of m.fs.flowsheet_cost.


Display Individual Costs
^^^^^^^^^^^^^^^^^^^^^^^^

There are three functions for displaying individual costs.

* display_total_plant_costs(fs)
* display_bare_erected_costs(fs)
* display_equipment_costs(fs)

Each one prints out a list of the costed blocks and the cost level of the function chosen.
The functions should be called after solving the model.

Checking Bounds
^^^^^^^^^^^^^^^

Currently, only the sCO2 module has support for checking bounds.

All costing methods have a range, outside of which the correlations become inaccurate.
Calling check_sCO2_costing_bounds(fs) will display which components are within the proper range
and which are outside it. It should be called after the model is solved.



References
----------

1. DOE/NETL-2015/1723 Cost and Performance Baseline for Fossil Energy Plants Volume 1a: Bituminus Coal (PC) and Natural Gas to Electricity Revision 3 and 4
2. DOE/NETL-341/013113 Quality Guidelines for Energy System Studies Capital Cost Scaling Methodology
3. NETL_PUB_21490 Techno-economic Evaluation of Utility-Scale Power Plants Based on the Indirect sCO2 Brayton Cycle. Charles White, David Gray, John Plunkett, Walter Shelton, Nathan Weiland, Travis Shultz. September 25, 2017
4. SCO2 Power Cycle Component Cost Correlations from DOE Data Spanning Multiple Scales and Applications. Nathan Weiland, Blake Lance, Sandeep Pidaparti. Proceedings of ASME Turbo Expo 2019: Turbomachinery Technical Conference and Exposition GT2019. June 17-21, 2019, Phoenix, Arizona, USA
