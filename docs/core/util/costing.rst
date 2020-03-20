Unit Model Costing
==================

.. contents:: Contents 
    :depth: 4

A generic costing "get_costing():" method has been developed for all the unit models in the generic unit model library. 
Currently, the method builds pruchase cost correlations from Seider et al. [1], however, this method allows the users to select other correlations to estimate the cost (i.e. other source, proprietary models, vendor quotes, etc.)

[1] Process and Product Design Principles: Synthesis, Analysis, and Evaluation. Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons. Chapter 22. Cost Accounting and Capital Cost Estimation

The get_costing method creates a block named "costing" that includes the cost model (variables, parameters, expressions, and constraints). 
The main variables are purchase cost and base cost, and each unit model has a different basis to estimate the cost. 
Table 1 shows the main variables that will be created in the block "costing", while table 2 shows the cost basis for all the unit models supported by this method.

Table 1. Main Variables added to the unit block ("self.costing").

=========================== ====================== ============ =============================================================================
Variable                    Symbol                 Units        Index Sets  Doc
=========================== ====================== ============ =============================================================================
Purchase cost               :math:`purchase\_cost` dollars      Heat transferred from hot side to the cold side
Base cost                   :math:`base\_cost`     unitless     Base cost
Pressure design             :math:`FP`             psig         Pressure design correction factor (only for heat exchangers)
Pump size factor            :math:`S`              gpm/ft^0.5   Pump Size factor used to compute pump base cost (only for pumps)
=========================== ====================== ============ =============================================================================
    
Table 2. Cost basis for each unit model.

=========================== =========================  =========== ===============================================================================
Unit Model                  Basis                      Units       Description
=========================== =========================  =========== ===============================================================================
heat exchanger              :math:`area`               ft^2        heat exchanger cost is computed as a function of the heat exchanger area
pump                        :math:`fluid_{work}`       ft^3/s      Pump fluid work is used to compute the pump size factor and head of the pump
compressor                  :math:`mechanical_{work}`  hp          Compressors are costed based on the mechanical work
turbine                     :math:`mechanical_{work}`  hp          Turbine cost is computed as a function of the mechanical work
=========================== =========================  =========== ===============================================================================

Example:
The example below demonstrates how to build the cost model as part of the heat exchanger unit model. 
As it can be observed, after calling the method 'm.fs.unit.get_costing()' the purchase cost and base cost variables are added to the unit block. 
These variables can be displayed using the display method "m.fs.unit.costing.purchase_cost.display()" and "m.fs.unit.costing.base_cost.display()"

Note that after building the costing model, the costing block is now part of the unit, and initialization of these constraints before solving the flowsheet is recommended.
The degrees of freedom remain unchanged because we added 2 independent variables and two constraints. 
Units: it is important to highlight that the costing method interrogates the property package to determine the units of this model, 
if the user provided the correct units in the metadata dictionary (see property models for additional information), the model units will be converted to the right units. 
For example: area is in m^2, while the cost correlations for heat exchanger cost require units to be in ft^2.

.. code:: python

    from pyomo.environ import (ConcreteModel, SolverFactory)
    from pyomo.util.calc_var_value import calculate_variable_from_constraint
    from idaes.core import FlowsheetBlock
    from idaes.generic_models.unit_models.heat_exchanger import \
        (HeatExchanger, HeatExchangerFlowPattern)
    from idaes.generic_models.properties import iapws95
    from idaes.core.util.model_statistics import degrees_of_freedom
    
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
    
    print(degrees_of_freedom(m))
    
    m.fs.unit.initialize()
    m.fs.unit.get_costing(L_factor='12ft')
    calculate_variable_from_constraint(
                m.fs.unit.costing.base_cost,
                m.fs.unit.costing.base_cost_eq)
    
    calculate_variable_from_constraint(
                m.fs.unit.costing.purchase_cost,
                m.fs.unit.costing.cp_cost_eq)
    
    opt = SolverFactory('ipopt')
    opt.options = {'tol': 1e-6,
                   'halt_on_ampl_error': 'no',
                   'max_iter': 50}
    results = opt.solve(m, tee=True)
    m.fs.unit.costing.base_cost.display()
    m.fs.unit.costing.purchase_cost.display()

Heat Exchanger Cost
-------------------

The Heat exchanger costing method is based on the area of the heat exchanger. 
This method computes the purchase cost (self.costing.purchase_cost) for a shell and tube heat 
exchanger (Eq. 22.43). First, the model computes the base cost (self.costing.base_cost = CB, for 4 differnt types
of heat exchangers, such as floating head, fixed head, U-tube, and
Kettle vaporizer) for a CE base cost index of 500. Then, additional factors will be used to calculate the final purchase cost, 
including construction material factor (FM), pressure design factor (FP), and tube length correction factor (FL).

.. math:: self.costing.purchase\_cost = FP*FM_{MAT}*FL*CB*(CE_{index}/500) (Eq. 22.43)

where:

* FP - is the pressure design correction factor

* FM_Mat - is the construction material correction factor

* FL - is the tube length correction factor

* CB - is the base cost

* CE - index is a global parameter that includes cost indexes for years 1980-2019

The heat exchanger costing method has three arguments, hx_type = heat exchanger type, FM_Mat = construction material factor, and FL = tube lenght factor.

* hx_type : 'floating_head', 'fixed_head', 'U-tube'*, 'Kettle_vap'

* material factor (FM): 'stain_steel'*, 'carb_steel'

* tube length (FL): '8ft'*, '12ft', '16ft', '20ft'

where '*' corresponds to the default options, FL and FM_MAT are pyomo-mutable parameters fixed based on user selection.

The base cost is computed dependind on the heat exchanger type selected by the user:

.. math:: self.costing.base\_cost = \exp{(\alpha_{1} - \alpha_{2}*\log{area} + \alpha_{3}*(\log{area})^{2})}


Table 3. Base cost factors for heat exchanger type.

================= ================== ================== ==================
Tube Length (ft)  :math:`\alpha_{1}` :math:`\alpha_{2}` :math:`\alpha_{3}`
================= ================== ================== ==================
floating_head     11.9052            0.8709             0.09005 
fixed_head        11.2927            0.8228             0.09861
U-tube            11.3852            0.9186             0.09790
Kettle_vap        12.2052            0.8709             0.09005
================= ================== ================== ==================


Table 4. Tube-Length correction factor.

================= =====
Tube Length (ft)  FL
================= =====
8                 1.25
12                1.12
16                1.05
20                1.00
================= =====

Construction material correction factor (FM_Mat) can be computed with Eq. 22.44 (Seider et al.)

.. math:: FM_{Mat} = a + (\frac{area}{100})^{b}     (Eq. 22.44)


Table 5. Materials of construction factors.

================================== ====== ======
Materials of Construction
Shell / Tube                       a      b
================================== ====== ======
Carbon steel / carbon steel        0.00   0.00
Carbon steel / brass               1.08   0.05
Carbon steel / stainless steel     1.75   0.13
Carbon steel / Monel               2.1    0.13
Carbon steel / titanium            5.2    0.16
Carbon steel / Cr-Mo steel         1.55   0.05
Cr-Mo steel / Cr-Mo steel          1.7    0.07
Stainless steel / stainless steel  2.7    0.07
Monel / Monel                      3.3    0.08
Titanium / titanium                9.6    0.06
================================== ====== ======


Pressure Changer Cost
---------------------
The pressure changer unit model is more complicated, because the pressure changer model can be imported to represent a pump, turbine, compressor, or a simply pressure changer (fan, blower, etc.).
The get_costing(): method currently supports costing of pumps, turbines, and compressors. The method authomatically interrogates the flowsheet object to determine if the unit is being used as a pump, turbine, or compressor. 
The additional arguments are required to build correlations for different type of pumps or compressors. 

Turbine Cost Model
*******************

We determine if the pressure changer is a turbine by interrogating the pressure changer config argument "config.compressor == False".
Turbine cost is based on the mechanical work of unit (work_mechanical), this correlation has been obtained using the NETL Report (DOE/NETL 2015).

.. math:: self.costing.purchase\_cost = 580*(mechanical_{work})^{0.81}

DOE/NETL, 2015, report. Cost and performance Baseline for Fossil Energy Plants. Volume 1a: Bituminous Coal (PC) and Natural Gas to Electricity. Revision 3

Pump Cost Model
****************

We determine if the pressure change is a pump if " unit.config.compressor == True, and unit.config.Thermodynamic.assumption.name == 'pump' "
Three main pump types are supported in this method. i) Centrifugal pumps, 2) External gear pumps, 3) Reciprocating Plunger pumps. Purchase cost is computed depending on user's inputs.
The argument pump_type will determine which cost correlations to be build (pump_type = 'centrifugal', 'external_gear', 'reciprocating')


Centrifugal pump (pump_type = 'centrifugal')
+++++++++++++++++++++++++++++++++++++++++++++

The centrifugal cost has two main components, the cost of the pump and the cost of the motor. The pump cost is based on the fluid work (work_fluid), pump head, and size factor. 
Additional arguments are required:

* pump_type_factor = '1.4'

* pump_motor_type_factor = 'open', 'enclosed', 'explosion_proof'


Based on users inputs the get_costing method builds base_cost and purchase_cost for both the pump and the motor. 
The unit purchase cost is obtained by adding the motor and pump costs.

.. math:: self.costing.purchase\_cost = self.costing.pump\_purchase\_cost + self.costing.motor\_purchase\_cost

To compute the purchase cost of the centrifugal pump, first we obtain the pump size factor (S) with Eq. 22.13, then we obtain the base cost with Eq. 22.14.
Finally, the purchase cost of the pump is obtained in Eq. 22.15.

.. math:: S = QH^{0.5}

.. math:: self.costing.pump\_base\_cost = \exp{(9.7171 - 0.6019*\log{S} + 0.0519*(\log{S})^{2})}

.. math:: self.costing.pump\_purchase\_cost = F_{T}*FM_{MAT}*self.costing.pump\_base\_cost*(CE_{index}/500)

where:

* S is the pump size factor (self.costing.size_factor)

* Q is the volumetric flowrate in gpm (depending on the model this variable can be found as self.unit.properties_in.flow_vol)

* H is the head of the pump in ft (which is defined as :math:`H = \Delta P/\rho_{liq}`)

* FT is the pump type factor (users must wisely select this factor based on the pump size factor, pump head range, and maximum motor hp)

* FM_Mat is the material factor for the pump

Table 6. Pump Type factor (Table 22.20 in Seider et al.).

====== ========= ======== ========= ========== ==================== =================
Case   FT factor # stages Shaft rpm Case-split Pump Head range (ft) Maximum Motor Hp
====== ========= ======== ========= ========== ==================== =================
'1.1'  1.00      1        3600      VSC        50  - 900            75 
'1.2'  1.50      1        1800      VSC        50  - 3500           200 
'1.3'  1.70      1        3600      HSC        100 - 1500           150  
'1.4'  2.00      1        1800      HSC        250 - 5000           250 
'2.1'  2.70      2        3600      HSC        50  - 1100           250 
'2.2'  8.90      2+       3600      HSC        100 - 1500           1450 
====== ========= ======== ========= ========== ==================== =================

For more details on how to select the FT factor, please see Seider et al.

Table 7. Materials of construction factors for centrifugal pumps and external gear pumps.

================= ======
Material Factor   FM_MAT
================= ======
Cast iron         1.00
Ductile iron      1.15
Cast steel        1.35
Bronze            1.90
Stainless steel   2.00
Hastelloy C       2.95 
Monel             3.30
Nickel            3.50
Titanium          9.70
================= ======

Electric Motor:

A centrifugal pump is usually driven by an electric motor, the self.costing.motor_purchase_cost is calculated based on the power consumption.

.. math:: self.motor_purchase_cost = FT * self.costing.motor\_base\_cost * (CE_{index}/500)  (Eq. 22.20)

.. math:: self.costing.motor\_base\_cost = \exp{(5.8259 + 0.13141\log{PC} + 0.053255(\log{PC})^{2} + 0.028628(\log{PC})^{3} - 0.0035549(\log{PC})^{4})}  (Eq. 22.19)

.. math:: PC = \frac{P_{T}}{\eta_{P}\eta_{M}} = \frac{P_{B}}{\eta_{M}} = \frac{Q H \rho}{33000\eta_{P}\eta_{M}}    (Eq. 22.16)

.. math:: \eta_{P} = -0.316 + 0.24015*\log{Q} - 0.01199(\log{Q})^{2}    (Eq. 22.17)

.. math:: \eta_{M} = 0.80 + 0.0319\log{PB} - 0.00182(\log{PB})^{2}   (Eq. 22.18)

Efficiencies are valid for PB in the range of 1 to 1500Hp and Q in the range of 50 to 5000 gpm

where:
* FT is the motor type factor

* PC is the power consumption in hp

* PT is the theoretical efficiency

* Q is the volumetric flowrate in gpm
* H is the pump head in ft

* PB is the pump brake hp

* nP is the fractional efficiency of the pump

* nM is the fractional efficiency of the motor

* :math:`\rho` is the liquid density in lb/gal

Table 8. FT Factors in Eq.(22.20) and Ranges for electric motors.

======================================== ======= =======
Type Motor Enclosure                     3600rpm 1800rpm
======================================== ======= =======
Open, drip-proof enclosure, 1 to 700Hp   1.0     0.90
Totally enclosed, fan-cooled, 1 to 250Hp 1.4     1.3
Explosion-proof enclosure, 1 to 25Hp     1.8     1.7
======================================== ======= =======

External Gear Pumps
++++++++++++++++++++

External gear pumps are not as common as the contrifugal pump, and various methods can be used to correlate purchase cost. 
Here the purchase cost is computed as a function of the volumetric flowrate (Q) in gpm.

.. math:: self.costing.pump\_base\_cost = \exp{(7.6964 + 0.1986\log{Q} + 0.0291(\log{Q})^{2})}           (Eq 22.21)

.. math:: self.costing.pump\_purchase\_cost = FM_{MAT} * self.costing.pump\_base\_cost * (CE_{index}/500)  (Eq. 22.22)


Reciprocal Plunger Pumps (pump_type = 'reciprocal'):
++++++++++++++++++++++++++++++++++++++++++++++++++++

.. math:: self.costing.pump\_base\_cost = \exp{(7.8103 + 0.26986\log{PB} + 0.06718(\log{PB})^{2})} (Eq. 22.23)

.. math:: self.costing.pump\_purchase\_cost = FM_{MAT} * self.costing.pump\_base\_cost * (CE_{index}/500)  (Eq. 22.22)

Table 9. Materials of construction factors for reciprocating plunger pumps.

=============== ==========
Material        FM_MAT
=============== ==========
Ductile iron    1.00
Ni-Al-Bronze    1.15
Carbon steel    1.50
Stainless steel 2.20
=============== ==========


Mover
******
We determine if the pressure changer is a mover if "unit.compressor == True and unit.config.Thermodynamic.assumption.name not 'pump ". 
If the unit represents a "Mover", the user can select to cost it as a compressor, fan, or blower. 
Therefore, the user must set the mover_type = 'compressor' or 'fan' or 'blower' (uper/lower case sensitive)

Compressor Cost
++++++++++++++++
The compressor cost is based on the mechanical work of the unit. Additional arguments are required to estimate the cost, such as compressor type, driver mover type, and material factor (FM_MAT).

* compressor_type = 'centrifugal', 'reciprocating', 'screw'

* driver_mover_type = 'electrical_motor', 'steam_turbine', 'gas_turbine'

* FM_mat = 'carbon_steel', 'stain_steel', 'nickel_alloy'

.. math:: self.costing.purchase\_cost = F_{D} F_{M} self.costing.base\_cost

.. math:: self.costing.base\_cost = \exp{(\alpha_{1} + \alpha_{2}*\log{mechanical_{work}})}

where: 
FD is the driver mover type factor and FM is the construction material factor.

Table 10. Compressor type factors.

================= ================== ==================
Compressor type   :math:`\alpha_{1}` :math:`\alpha_{2}`
================= ================== ==================
Centrifugal       7.5800             0.80
Reciprocating     7.9661             0.80
Screw Compressor  8.1238             0.7243
================= ================== ==================


Table 11. Driver mover type (for compressors only).

=============== ==========
Mover type      FD
=============== ==========
Electric Mover  1.00
Steam Turbine   1.15
Gas Turbine     1.25
=============== ==========

Table 12. Material of construction factor (for compressors only).

=============== ==========
Material        FM
=============== ==========
Cast iron       1.00
Stainless steel 1.15
Nickel alloy    1.25
=============== ==========

Fan and Blower Costing (work in progress)
+++++++++++++++++++++++++++++++++++++++++