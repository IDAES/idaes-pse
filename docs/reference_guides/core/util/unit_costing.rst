Unit Model Costing
==================

The IDAES Process Modeling Framework includes support for incorporating costing of unit 
operations into a flowsheet to allow for calculation and optimization of process costs. 
Cost Correlations are implemented using unit costing sub-modules to allow users to easily develop 
and incorporate their own costing models.

.. contents:: Contents 
    :depth: 4

Introduction
------------

All unit models within the core IDAES model library include a `get_costing` method which can be called to include
cost correlations for an instance of that unit. The `get_costing` method for each unit takes a number of arguments used 
to specify the basis for costing each piece of equipment. Details are given for each unit model later in this documentation, 
however, all `get_costing` methods take the following two arguments:
 
* module - this argument specifies the costing module to use when constructing the constraints and associated variables. if not provided, this defaults to the standard IDAES costing module.
* year - this argument sets the year to which all costs should be normalized (CE index 2010 to 2019)
 
When `get_costing` is called on an instance of a unit model, a new sub-block is created 
on that unit named `costing` (i.e. `flowsheet.unit.costing`). All variables and constraints related to costing will be 
constructed within this new block (see detailed documentation for each unit for details on these variables and constraints).

In addition, the first time `get_costing` is called for a unit operation within a flowsheet, an additional `costing` block is created 
on the flowsheet object (i.e. `flowsheet.unit.costing`) in order to hold any global parameters relating to costing. The most 
common of these paramters is the cost normalization parameter based on the year selected by the user.

The unit costing module also contains an `initialize` method which can be used to estimate initial values for costing variables based on the current state of the associated unit model. This method can be called directly from the `unit_costing` module to initialize a specific costing block, or can be incorporated into a unit model initialization procedure. This method has been incorporated into the `initialize` method of all the models in the core unit model library.
Therefore, if get_costing() is called before unit.initialize(), the initialize method will deactivate the costing block, initialize the unit model as normal, and then activate the costing block and initialize costing block.

.. note:: The global paramters are created when the first instance of `get_costing` is called and use the values provided there for initialization. Subsequent `get_costing` calls use the existing paramters, and do not change the initialized values. i.e. any "year" argument provided to a `get_costing` call after the first will be ignored.

 
Table 1. Main Variables added to the unit block ("self.costing").

=========================== ============================ ============ =============================================================================
Variable                    Symbol                       Units        Notes
=========================== ============================ ============ =============================================================================
Purchase cost               :math:`purchase\_cost`       dollars      Purchase cost
Base cost per unit          :math:`base\_cost_per_unit`  unitless     Base cost per unit
Base cost                   :math:`base\_cost`           unitless     Base cost (base cost per unit * number of units)
Number of units             :math:`number\_of\_units`    unitless     Number of units to be costed (to take advantage of the economics of scale)
=========================== ============================ ============ =============================================================================

.. note:: number of units by default is fixed to 1 and the user must unfix this variable to optimize the number of units. Also, `number of units` can be built as a continuous variable or an integer variable. If latest, the user must provide an mip solver. Use the global costing argument for this purpose (integer_n_units=True or False).

Example
-------
Below is a simple example of how to add cost correlations to a flowsheet including a heat exchanger using the default IDAES costing module.


.. code:: python

    from pyomo.environ import (ConcreteModel, SolverFactory)
    from pyomo.util.calc_var_value import calculate_variable_from_constraint
    from idaes.core import FlowsheetBlock
    from idaes.models.unit_models.heat_exchanger import \
        (HeatExchanger, HeatExchangerFlowPattern)
    from idaes.models.properties import iapws95
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
    
    m.fs.unit.get_costing(module=costing, length_factor='12ft')

    m.fs.unit.initialize()

    opt = SolverFactory('ipopt')
    opt.options = {'tol': 1e-6, 'max_iter': 50}
    results = opt.solve(m, tee=True)

Units
-----

It is important to highlight that the costing method interrogates the property 
package to determine the units of this model, if the user provided the correct 
units in the metadata dictionary (see property models for additional information), 
the model units will be converted to the right units. 
For example: in this example area is in m^2, while the cost correlations for heat 
exchangers require units to be in ft^2. Therefore, the costing method will convert 
the units to ft^2. The use of Pyomo-unit conversion tools is under development.

IDAES Costing Module
--------------------

A default costing module has been developed primarily based on base cost and purchase cost correlations 
from the following reference with some exceptions (noted in the documentation as appropiate).

Process and Product Design Principles: Synthesis, Analysis, and Evaluation. Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons. Chapter 22. Cost Accounting and Capital Cost Estimation

Users should refer to the reference above for details of the costing correlations, however, a summary of this methods is provided below.
    
Table 2. Cost basis for each unit model.

=========================== =========================  ===========
Unit Model                  Basis                      Units      
=========================== =========================  ===========
heat exchanger              :math:`area`               ft^2       
pump                        :math:`fluid_{work}`       ft^3/s     
compressor                  :math:`mechanical_{work}`  hp         
turbine                     :math:`mechanical_{work}`  hp         
vessels                     :math:`D and L`            ft         
fired heaters               :math:`heat\_duty`         BTU/hr         
=========================== =========================  ===========


Heat Exchanger Cost
^^^^^^^^^^^^^^^^^^^

.. module:: idaes.core.util.unit_costing

The purchse cost is computed based on the base unit cost and three correction factors (Eq. 22.43 in Seider et al.). The base cost is computed depending on the heat exchanger type selected by the user:

.. math:: self.costing.purchase\_cost = pressure\_factor*material\_factor*L\_factor*self.costing.base\_cost*(CE_{index}/500)

.. math:: self.costing.base\_cost\_per\unit = \exp{(\alpha_{1} - \alpha_{2}*\log{area*hx\_os} + \alpha_{3}*(\log{area*hx\_os})^{2})}

.. math:: self.costing.base\_cost = self.costing.base\_cost\_per\unit * self.costing.number\_of\_units

.. math:: area  = self.area / self.costing.number\_of\_units

where:

* pressure_factor - is the pressure design correction factor
* material_factor - is the construction material correction factor
* length_factor - is the tube length correction factor
* CE_index - is a global parameter for Chemical Enginering cost index for years 2010-2019
* hx_os - heat exchanger oversize factor (default = 1)
* area is a reference object and (self.area is the model variable)

The heat exchanger costing method has three arguments, hx_type = heat exchanger type, FM_Mat = construction material factor, and FL = tube length factor.

* hx_type : 'floating_head', 'fixed_head', 'U-tube'\*, 'Kettle_vap'
* material factor (Mat_factor): 'stain_steel'\*, 'carb_steel'
* tube length (length_factor): '8ft', '12ft'\*, '16ft', '20ft'

where '*' corresponds to the default options, FL and FM_MAT are pyomo-mutable parameters fixed based on user selection.


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

.. math:: material\_factor = a + (\frac{area}{100})^{b}


Table 5. Materials of construction factors.

================================== ====== ======
Materials of Construction
Shell / Tube                       a      b
================================== ====== ======
carbon steel/carbon steel          0.00   0.00
carbon steel/brass                 1.08   0.05
carbon steel/stainless steel       1.75   0.13
carbon steel/monel                 2.1    0.13
carbon steel/titanium              5.2    0.16
carbon steel/Cr-Mo steel           1.55   0.05
Cr-Mo steel/Cr-Mo steel            1.7    0.07
stainless steel/stainless steel    2.7    0.07
monel/monel                        3.3    0.08
titanium/titanium                  9.6    0.06
================================== ====== ======

Note that `Mat_factor` argument should be provided a string, for example: Mat_factor:'carbon steel/carbon steel'.

Pressure Changer Cost
^^^^^^^^^^^^^^^^^^^^^

The costing of a pressure changer unit model is more complicated, because the pressure changer 
model can be imported into the flowsheet object representing a pump, turbine, compressor, or a 
simply pressure changer (fan, blower, etc.). The `get_costing` method currently supports costing of pumps, turbines, and compressors. The method authomatically interrogates the flowsheet object to determine if the unit is being used as a pump, turbine, or compressor. 

The `get_costing` method authomatically determines if the unit model is being used as a pump, 
turbine, or compressor based on the `compressor` and `thermodynamic_assumption` configuration 
arguments provided by the user where creating the unit model. A summary of the decision logic is shown below.


========== =========== =========================
Unit Type  compressor  thermodynamic_assumption
========== =========== =========================
Turbine    False       Any
Pump       True        pump
Mover      True        not pump
========== =========== =========================

Additionally, some unit types have different sub-types which can be costed appropiately. In these cases, 
an additional argument is provided to `get_costing` to identify the sub-type to use which is detailed below.

Turbine Cost Model
""""""""""""""""""
The turbine cost is based on the mechanical work of unit (work_mechanical), this correlation has been obtained using the NETL Report (DOE/NETL 2015).

.. math:: self.costing.purchase\_cost = 580*(mechanical_{work})^{0.81}

DOE/NETL, 2015, report. Cost and performance Baseline for Fossil Energy Plants. Volume 1a: Bituminous Coal (PC) and Natural Gas to Electricity. Revision 3

Pump Cost Model
""""""""""""""""

Three subtypes are supported for costing of pumps, which can be set using the "pump_type" argument.

1) Centrifugal pumps (pump_type='centrifugal')
2) External gear pumps (pump_type='external')
3) Reciprocating Plunger pumps (pump_type='reciprocating')


Centrifugal Pump
++++++++++++++++

The centrifugal pump cost has two main components, the cost of the pump and the cost of the motor. The pump cost is based on the fluid work (work_fluid), pump head, and size factor. 
Additional arguments are required:

* pump_type_factor = '1.4' (see Table 6)
* pump_motor_type_factor = 'open', 'enclosed', 'explosion_proof'


Based on user's inputs the get_costing method builds base_cost and purchase_cost for both the pump and the motor. 
The unit purchase cost is obtained by adding the motor and pump costs.

.. math:: self.costing.purchase\_cost = self.costing.pump\_purchase\_cost + self.costing.motor\_purchase\_cost

To compute the purchase cost of the centrifugal pump, first we obtain the pump size factor (S) with Eq. 22.13, then we obtain the base cost with Eq. 22.14.
Finally, the purchase cost of the pump is obtained in Eq. 22.15. (Seider et al.)

.. math:: S = QH^{0.5}

.. math:: self.costing.pump\_base\_cost\_per\unit = \exp{(9.7171 - 0.6019*\log{S} + 0.0519*(\log{S})^{2})}

.. math:: self.costing.pump\_purchase\_cost = F_{T}*material\_factor*self.costing.pump\_base\_cost*(CE_{index}/500)

.. math:: self.costing.base\_cost = self.costing.pump\_base\_cost\_per\unit * self.costing.number\_of\_units

.. math:: Q  = self.Q / self.costing.number\_of\_units

.. note:: the same number of units have been considered for pumps and the pump motor

where:

* S is the pump size factor (`self.costing.size_factor`)
* Q is the volumetric flowrate in gpm (depending on the model this variable can be found as self.unit.properties_in.flow_vol)
* H is the head of the pump in ft (`self.pump_head`; which is defined as :math:`H = \Delta P/\rho_{liq}`)
* FT is a parameter fixed based on the pump_type_factor argument (users must wisely select this factor based on the pump size factor, pump head range, and maximum motor hp)
* material_factor is the material factor for the pump

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
cast iron         1.00
ductile iron      1.15
cast steel        1.35
bronze            1.90
stainless steel   2.00
hastelloy C       2.95 
monel             3.30
nickel            3.50
titanium          9.70
================= ======

Electric Motor:

A centrifugal pump is usually driven by an electric motor, the `self.costing.motor_purchase_cost` is calculated based on the power consumption.

.. math:: self.motor_purchase_cost = FT * self.costing.motor\_base\_cost * (CE_{index}/500)  (Eq. 22.20)

.. math:: self.costing.motor\_base\_cost = self.costing.motor\_base\_cost\_per\unit * self.costing.number\_of\_units

.. math:: Q  = self.Q / self.costing.number\_of\_units

.. math:: self.costing.self.costing.motor\_base\_cost\_per\unit = \exp{(5.8259 + 0.13141\log{PC} + 0.053255(\log{PC})^{2} + 0.028628(\log{PC})^{3} - 0.0035549(\log{PC})^{4})}  (Eq. 22.19)

.. math:: PC = \frac{P_{T}}{\eta_{P}\eta_{M}} = \frac{P_{B}}{\eta_{M}} = \frac{Q H \rho}{33000\eta_{P}\eta_{M}}    (Eq. 22.16)

.. math:: \eta_{P} = -0.316 + 0.24015*\log{Q} - 0.01199(\log{Q})^{2}    (Eq. 22.17)

.. math:: \eta_{M} = 0.80 + 0.0319\log{PB} - 0.00182(\log{PB})^{2}   (Eq. 22.18)

Efficiencies are valid for PB in the range of 1 to 1500Hp and Q in the range of 50 to 5000 gpm

where:

* motor_FT is the motor type correction factor
* PC is the power consumption in hp (`self.power_consumption_hp`; coded as a pyomo expression)
* Q is the volumetric flowrate in gpm (`self.Q_gpm`)
* H is the pump head in ft (`self.pump_head`)
* PB is the pump brake hp (`self.work`)
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
+++++++++++++++++++

External gear pumps are not as common as the contrifugal pump, and various methods can be used to correlate base cost. Eq. 22.21 in Seider et al.
Here the purchase cost is computed as a function of the volumetric flowrate (Q) in gpm Eq. 22.22 in Seider et al.


.. math:: self.costing.pump\_purchase\_cost = material\_factor * self.costing.pump\_base\_cost * (CE_{index}/500)

.. math:: self.costing.pump\_base\_cost = self.costing.pump\_base\_cost\_per\unit * self.costing.number\_of\_units

. math:: self.costing.self.costing.pump\_base\_cost\_per\unit  = \exp{(7.6964 + 0.1986\log{Q} + 0.0291(\log{Q})^{2})}

.. math:: Q  = self.Q / self.costing.number\_of\_units

Reciprocating Plunger Pumps
+++++++++++++++++++++++++++

The cost correlation method used here is based on the brake horsepower (PB).

.. math:: self.costing.pump\_purchase\_cost = material\_factor * self.costing.pump\_base\_cost * (CE_{index}/500)  (Eq. 22.22)

.. math:: self.costing.pump\_base\_cost = self.costing.pump\_base\_cost\_per\unit * self.costing.number\_of\_units

.. math:: self.costing.pump\_base\_cost\_per\unit = \exp{(7.8103 + 0.26986\log{PB} + 0.06718(\log{PB})^{2})} (Eq. 22.23)

.. math:: PB = f(Q)

.. math:: Q  = self.Q / self.costing.number\_of\_units

Table 9. Materials of construction factors for reciprocating plunger pumps.

=============== ==========
Material        Mat_factor
=============== ==========
ductile iron    1.00
Ni-Al-Bronze    1.15
carbon steel    1.50
stainless steel 2.20
=============== ==========


Mover (Compressor, Fan, Blower)
"""""""""""""""""""""""""""""""

If the unit represents a "Mover", the user can select to cost it as a compressor, fan, or blower. 
Therefore, the user must set the "mover_type" argument.

* mover_type= 'compressor' or 'fan' or 'blower' (uper/lower case sensitive)

Compressor Cost
+++++++++++++++
The compressor cost is based on the mechanical work of the unit. 
Additional arguments are required to estimate the cost such as compressor type, 
driver mover type, and material factor (Mat_factor).

* compressor_type = 'centrifugal', 'reciprocating', 'screw'
* driver_mover_type = 'electrical_motor', 'steam_turbine', 'gas_turbine'
* Mat_factor = 'carbon_steel', 'stain_steel', 'nickel_alloy'

.. math:: self.costing.purchase\_cost = (CE_{index}/500)* F_{D} * material\_factor * self.costing.base\_cost

.. math:: self.costing.base\_cost = self.costing.base\_cost\_per\_unit * self.costing.number\_of\_units

.. math:: self.costing.base\_cost\_per\_unit = \exp{(\alpha_{1} + \alpha_{2}*\log{mechanical_{work}})}

.. math:: mechanical_{work} = self.mechanical_{work} / self.costing.number\_of\_units

where: 

* FD is the driver mover type factor and FM is the construction material factor.

Table 10. Compressor type factors.

================= ================== ==================
Compressor type   :math:`\alpha_{1}` :math:`\alpha_{2}`
================= ================== ==================
Centrifugal       7.5800             0.80
Reciprocating     7.9661             0.80
Screw Compressor  8.1238             0.7243
================= ================== ==================


Table 11. Driver mover type (for compressors only).

=============== ===============
Mover type      FD (mover_type)
=============== ===============
Electric Mover  1.00
Steam Turbine   1.15
Gas Turbine     1.25
=============== ===============

Table 12. Material of construction factor (for compressors only).

=============== ===========
Material        Mat_factor
=============== ===========
Cast iron       1.00
Stainless steel 1.15
Nickel alloy    1.25
=============== ===========

Fan Cost
++++++++
The fan cost is a function of the actual cubic feet per minute (Q) entering the fan.
Additional arguments are required to estimate the fan cost such as mover_type='fan', fan_head_factor,
fan_type, and material factor (Mat_factor).

* fan_type = 'centrifugal_backward', 'centrifugal_straight', 'vane_axial', 'tube_axial'
* fan_head_factor = see table 14
* Mat_factor = 'carbon_steel', 'fiberglass', 'stain_steel', 'nickel_alloy'

To select the correct fan type users must calculate the total head in inH2O and select the proper fan type from table 13.
Additionally, the user must select the head factor (head_factor) from table 14.

Table 13. Typical Operating Ranges of Fans

=========================== ================ =================
Fan type                    Flow rate (ACFM)  Total head inH2O
ACFM^a inH2O
=========================== ================ =================
Centrifugal backward curved  1000-100000      1-40
Centrifugal straight radial  1000-20000       1-30
Vane axial                   1000-800000      0.02-16
Tube axial                   2000-800000      0.00-10
=========================== ================ =================


Finally, the purchase cost of the fan is given by base cost, material factor, and fan head factor. While, the base cost is given as a function of the ACFM (Q).

.. math:: self.costing.purchase\_cost = (CE_{index}/500) * head\_factor * material\_factor * self.costing.base\_cost

.. math:: self.costing.base\_cost = self.costing.base\_cost\_per\_unit * self.costing.number\_of\_units

.. math:: self.costing.base\_cost\_per\_unit = \exp{(\alpha_{1} - \alpha_{2}*\log{Q} + \alpha_{3}*(\log{Q})^{2})}

.. math:: Q  = self.Q / self.costing.number\_of\_units


Table 14. Head Factor, FH, for fans

============= =========================== =========================== ========== ==========
Head (in H2O) Centrifugal backward curved Centrifugal straight radial Vane axial Tube Axial
============= =========================== =========================== ========== ========== 
5-8           1.15                        1.15                        1.15       1.15
9-15          1.30                        1.30                        1.30 
16-30         1.45                        1.45
31-40         1.55 
============= =========================== =========================== ========== ========== 

Table 15. Materials of construction factor 

================ ======
Material Factor  FM
================ ======
carbon_steel     1
fiberglass       1.8
stain_steel      2.5
nickel_alloy     5.0
================ ======

Blower Cost
+++++++++++

The blower cost is based on the brake horsepower, which can be calculated with the inlet volumetric flow rate and pressure (cfm and lbf/in^2, respectivelly).
Additional arguments are required to estimate the blower cost such as mover_type='blower', blower_type, and material of construction factor (Mat_factor).

* blower_type = 'centrifugal', 'rotary'
* Mat_factor = 'carbon_steel', 'aluminum', 'fiberglass', 'stain_steel', 'nickel_alloy'

where the material factors given in table 15 for the fans can be used. In addition, centrifugal blowers are available with cast aluminum blades with Mat_factor = 0.60.

The purchase cost is given by the material factor and base cost. While, the base cost is given by the power consumption in horsepower (Pc). 

.. math:: self.costing.purchase\_cost = material\_factor * self.costing.base\_cost

.. math:: self.costing.base\_cost = self.costing.base\_cost\_per\_unit * self.costing.number\_of\_units

Centrigugal turbo blower (valid from PC = 5 to 1000 Hp):

.. math:: self.costing.base\_cost\_per\_unit = \exp{(6.8929 + 0.7900*\log{Pc})}

Rotary straight-lobe blower (valid from PC = 1 to 1000 Hp):

.. math:: self.costing.base\_cost\_per\_unit = \exp{(7.59176 + 0.79320*\log{Pc} - 0.012900*(\log{Pc})^{2})}

.. math:: Pc = f(Q)

.. math:: Q = self.Q / self.costing.number\_of\_units


Fired Heater
^^^^^^^^^^^^
Indirect fired heaters, also called fired heaters, process heaters, and furnaces, are used to heat or vaporize process streams at elevated temperatures (beyond where steam is usually employed).
This method computes the purchase cost of the fired heater based on the heat duty, fuel used (fired_type), pressure design, and materials of construction (Mat_factor).

* fuel_type = 'fuel', 'reformer', 'pyrolysis', 'hot_water', 'salts', 'dowtherm_a', 'steam_boiler'
* Mat_factor = see table 16

Table 16. Materials of construction factor

=============== ======
Material Factor (FM)
=============== ======
carbon_steel    1
Cr-Mo_alloy     1.4
stain_steel     1.7
=============== ======

The pressure design factor is given by (where P is pressure in psig and it is valid between 500 to 3000 psig):

.. math:: self.pressure\_factor == 0.986 - 0.0035*(P/500.00) + 0.0175*(P/500.00)^{2}

The base cost changes depending on the fuel type:
fuel:

.. math:: self.costing.base\_cost\_per\_unit = \exp{(0.32325 + 0.766*\log{heat\_duty})}

reformer:

.. math:: self.costing.base\_cost\_per\_unit = 0.859*heat\_duty^{0.81}

pyrolysis:

.. math:: self.costing.base\_cost\_per\_unit = 0.650*heat\_duty^{0.81}

hot_water:

.. math:: self.costing.base\_cost\_per\_unit = \exp{(9.593- 0.3769*\log{heat\_duty} + 0.03434*(\log{heat\_duty})^{2})}

salts:

.. math:: self.costing.base\_cost\_per\_unit = 12.32*heat\_duty^{0.64}

dowtherm_a:

.. math:: self.costing.base\_cost\_per\_unit = 12.74*heat\_duty^{0.65}

steam_boiler:

.. math:: self.costing.base\_cost\_per\_unit = 0.367*heat\_duty^{0.77}

.. math:: self.costing.base\_cost = self.costing.base\_cost\_per\_unit * self.costing.number\_of\_units

Finally, the purchase cost is given by:

.. math:: self.purchase\_cost = (CE_{index}/500) * pressure\_design * material\_factor * base\_cost


Cost of Pressure Vessels and Towers for Distillation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pressure vessels cost is based on the weight of the vessel, the cost of platforms and ladders can be included, and the cost of internal packing or trays can be calculated as well. 
This method constructs by defaul the cost of pressure vessels with platforms and ladders, and trays cost can be calculated if trays=True. This method requires a few arguments to build the cost of vessel. 
We recommend using this method to cost reactors (CSTR or PFR), flash tanks, vessels, and distillation columns.

* alignment = 'horizontal', 'vertical'
* Mat_factor = 'carbon_steel'
* weight_limit = 'option1', 'option2' (option 1: 1000 to 920,000 lb, option 2: 9000 to 2.5M lb only for vertical vessels)
* L_D_range = 'option1', 'option2' (option 1: 3 < D < 21, 12 < L < 40; option 2: 3 < D < 24, 27 < L < 170; all in ft D: diameter, L: length) only for vertical vessels
* PL='True', 'False': to build platforms and ladders cost
* plates = 'True', 'False': to build tray cost for distillation columns
* tray_mat_factor = 'carbon_steel' see table 18
* tray_type = 'sieve'
* number_tray = 10
* ref_parameter_diameter=None
* ref_parameter_length=None


By adding reference parameter, the method can be constructed in any pyomo costing block.
Since the generic models do not include the variables required to cost these type of units, the user must create the blocks and variables.
For example: m.fs.unit = Block(), m.fs.unit.diameter = Var(), m.fs.unit.length = Var(). Then m.fs.unit.costing = pyo.Block() and call vessel_costing method = vessel_costing(m.fs.unit.costing, args).

Table 17. Materials of construction factor and material density

================== ====== ==========================
Material Factor    (FM)   methal density (lb/in^3)
================== ====== ==========================
carbon_steel       1      0.284
low_alloy_steel    1.2    0.271
stain_steel_304    1.7    0.270
stain_steel_316    2.1    0.276
carpenter_20CB-3   3.2    0.292
nickel_200         5.4    0.3216
monel_400          3.6    0.319
inconel_600        3.9    0.3071
incoloy_825        3.7    0.2903
titanium           7.7    0.1628
================== ====== ==========================

Vessel Cost
"""""""""""

The weight of the unit is calculated based on the methal density, length, Diameter, and shell thickness. `shel_thickness` is a parameter initialized to 1.25, 
however, the user must calculate the shell wall minimum thickness computd from the ASME pressure vessel code (tp) add the average vessel thickness, the necessary wall thickness (tE), and select the appropriate shell_thickness.

.. math:: self.weight == \pi * ((D*12) + self.shell\_thickness) * ((L*12)+(0.8*D*12))*self.shell\_thickness*self.material\_density

The base cost of the vessel is given by:
Horizontal vessels (option1: 1000 < W < 920,000 lb):

.. math:: self.costing.base\_cost\_per\_unit = \exp{(8.9552 - 0.2330*\log{weight} + 0.04333*(\log{weight})^{2})}

Vertical vessels (option1: 4200 < W < 1M lb):

.. math:: self.costing.base\_cost\_per\_unit = \exp{(8.9552 - 0.2330*\log{weight} + 0.04333*(\log{weight})^{2})}

Vertical vessels (option2: 9,000 < W < 2.5M lb):

.. math:: self.costing.base\_cost\_per\_unit = \exp{(7.2756 - 0.18255*\log{weight} + 0.02297*(\log{weight})^{2})}

.. math:: self.costing.base\_cost = self.costing.base\_cost\_per\_unit * self.costing.number\_of\_units

.. math:: weight = self.weight / self.costing.number\_of\_units

The vessel purchase cost is given by:

.. math:: self.vessel\_purchase\_cost = (CE_{index}/500) * material\_factor * self.base\_cost + (self.base\_cost\_platf\_ladders * self.costing.number\_of\_units)

note that if PL = 'False', the cost of platforms and ladders is not included.

The final purchase cost is given by:

.. math:: self.purchase\_cost = self.vessel\_purchase\_cost + (self.purchase\_cost\_trays * self.costing.number\_of\_units)

note that if plates='False', the cost of trays is not included.


Base Cost of Platforms and ladders
""""""""""""""""""""""""""""""""""
The cost of platforms and ladders is based on the diamter and length in ft.
Horizontal vessels (option1: 3 < D < 12 ft):

.. math:: self.base\_cost\_platf\_ladders = 20059*D^{0.20294}

Vertical vessels (option1: 3 < D < 12 ft and 12 < L  < 40 ft):

.. math:: self.base\_cost\_platf\_ladders = 361.8*D^{0.73960} * L^{0.70684}

Vertical vessels (option2: 3 < D < 24 ft and 27 < L  < 170 ft):

.. math:: self.base\_cost\_platf\_ladders = 300.9*D^{0.63316} * L^{0.80161}


Purchase Cost of Plates
"""""""""""""""""""""""

The cost of plates is based on the number or trays, the type of trays used, and materials of construction. 
Tray type factor (tray_factor) is 1.0 for sieve trays, 1.18 for valve trays (valve), and 1.87 for bubble cap trays (bubble_cap). The number of trays factor (number_tray_factor) is equal to 1 if the number of trays is greater than 20. 
However, if the number of trays is less than 20, the number_tray_factor is given by:

.. math:: self.number\_tray\_factor = \frac{2.25}{1.0414^{NT}}

The materials of construction factor is calculated using the following equation:

.. math:: \alpha_1 + \alpha_2 * D

where alphas for different materials of construction are given in table 18.

Table 18. Materials of construction factor

================== ====== =======
Material           alpha1 alpha2
================== ====== =======
carbon_steel       1      0
stain_steel_303    1.189  0.0577
stain_steel_316    1.401  0.0724
carpenter_20CB-3   1.525  0.0788  
monel_400          2.306  0.1120
================== ====== =======

The tray base cost is then calculated as:

.. math:: self.base\_cost\_trays = 468.00*\exp{(0.1739*D)}

The purchase cost of the trays is given by:

.. math:: self.purchase\_cost\_trays = (CE_{index}/500)* self.number\_trays * self.number\_tray\_factor * self.type\_tray\_factor * self.tray\_material\_factor * self.base\_cost\_trays
