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

.. note:: This is a work in progress, and costing is currently only implemented for pressure changers and heat exchanger.

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

.. note:: The global paramters are created when the first instance of `get_costing` is called and use the values provided there for initialization. Subsequent `get_costing` calls use the existing paramters, and do not change the initialized values. i.e. any "year" argument provided to a `get_costing` call after the first will be ignored.

 
Table 1. Main Variables added to the unit block ("self.costing").

=========================== ====================== ============ =============================================================================
Variable                    Symbol                 Units        Notes
=========================== ====================== ============ =============================================================================
Purchase cost               :math:`purchase\_cost` dollars      Purchase cost
Base cost                   :math:`base\_cost`     unitless     Base cost
=========================== ====================== ============ =============================================================================

Example
-------
Below is a simple example of how to add cost correlations to a flowsheet including a heat exchanger using the default IDAES costing module.


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
    
    m.fs.unit.initialize()
    m.fs.unit.get_costing(module=costing, L_factor='12ft')
    # initialize costing equations
    calculate_variable_from_constraint(
                m.fs.unit.costing.base_cost,
                m.fs.unit.costing.base_cost_eq)
    
    calculate_variable_from_constraint(
                m.fs.unit.costing.purchase_cost,
                m.fs.unit.costing.cp_cost_eq)
    
    opt = SolverFactory('ipopt')
    opt.options = {'tol': 1e-6, 'max_iter': 50}
    results = opt.solve(m, tee=True)

Units
-----

It is important to highlight that the costing method interrogates the property 
package to determine the unit of this model, if the user provided the correct 
units in the metadata dictionary (see property models for additional information), 
the model units will be converted to the right units. 
For example: in this example area is in m^2, while the cost correlations for heat 
exchangers require units to be in ft^2. Therefore, the costing method will convert 
the units to ft^2. The use of Pyomo-unit conversion tools is under development.

IDAES Costing Module
--------------------

A default costing module has been developed primarily based on purchase cost correlations 
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
=========================== =========================  ===========


Heat Exchanger Cost
^^^^^^^^^^^^^^^^^^^

.. module:: idaes.core.util.unit_costing

The purchse cost is computed based on the base unit cost and three correction factors. The base cost is computed depending on the heat exchanger type selected by the user:

.. math:: self.costing.purchase\_cost = FP*FM_{MAT}*FL*self.costing.base\_cost*(CE_{index}/500) (Eq. 22.43)

.. math:: self.costing.base\_cost = \exp{(\alpha_{1} - \alpha_{2}*\log{area} + \alpha_{3}*(\log{area})^{2})}

where:

* FP - is the pressure design correction factor

* FM_Mat - is the construction material correction factor

* FL - is the tube length correction factor

* CE - index is a global parameter that includes cost indexes for years 2010-2019

The heat exchanger costing method has three arguments, hx_type = heat exchanger type, FM_Mat = construction material factor, and FL = tube lenght factor.

* hx_type : 'floating_head', 'fixed_head', 'U-tube', 'Kettle_vap'

* material factor (FM): 'stain_steel', 'carb_steel'

* tube length (FL): '8ft', '12ft', '16ft', '20ft'

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

The centrifugal cost has two main components, the cost of the pump and the cost of the motor. The pump cost is based on the fluid work (work_fluid), pump head, and size factor. 
Additional arguments are required:

* pump_type_factor = '1.4' (see Table 6)

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

* S is the pump size factor (`self.costing.size_factor`)

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

A centrifugal pump is usually driven by an electric motor, the `self.costing.motor_purchase_cost` is calculated based on the power consumption.

.. math:: self.motor_purchase_cost = FT * self.costing.motor\_base\_cost * (CE_{index}/500)  (Eq. 22.20)

.. math:: self.costing.motor\_base\_cost = \exp{(5.8259 + 0.13141\log{PC} + 0.053255(\log{PC})^{2} + 0.028628(\log{PC})^{3} - 0.0035549(\log{PC})^{4})}  (Eq. 22.19)

.. math:: PC = \frac{P_{T}}{\eta_{P}\eta_{M}} = \frac{P_{B}}{\eta_{M}} = \frac{Q H \rho}{33000\eta_{P}\eta_{M}}    (Eq. 22.16)

.. math:: \eta_{P} = -0.316 + 0.24015*\log{Q} - 0.01199(\log{Q})^{2}    (Eq. 22.17)

.. math:: \eta_{M} = 0.80 + 0.0319\log{PB} - 0.00182(\log{PB})^{2}   (Eq. 22.18)

Efficiencies are valid for PB in the range of 1 to 1500Hp and Q in the range of 50 to 5000 gpm

where:

* FT is the motor type correction factor

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
+++++++++++++++++++

External gear pumps are not as common as the contrifugal pump, and various methods can be used to correlate purchase cost. 
Here the purchase cost is computed as a function of the volumetric flowrate (Q) in gpm.

.. math:: self.costing.pump\_base\_cost = \exp{(7.6964 + 0.1986\log{Q} + 0.0291(\log{Q})^{2})}           (Eq 22.21)

.. math:: self.costing.pump\_purchase\_cost = FM_{MAT} * self.costing.pump\_base\_cost * (CE_{index}/500)  (Eq. 22.22)


Reciprocating Plunger Pumps
+++++++++++++++++++++++++++

The cost correlation method used here is based on the brake horsepower (PB).

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
"""""

If the unit represents a "Mover", the user can select to cost it as a compressor, fan, or blower. 
Therefore, the user must set the "mover_type" argument.
* mover_type= 'compressor' or 'fan' or 'blower' (uper/lower case sensitive)

Compressor Cost
+++++++++++++++
The compressor cost is based on the mechanical work of the unit. 
Additional arguments are required to estimate the cost such as compressor type, 
driver mover type, and material factor (FM_MAT).

* compressor_type = 'centrifugal', 'reciprocating', 'screw'

* driver_mover_type = 'electrical_motor', 'steam_turbine', 'gas_turbine'

* FM_mat = 'carbon_steel', 'stain_steel', 'nickel_alloy'

.. math:: self.costing.purchase\_cost = F_{D} F_{M} self.costing.base\_cost

.. math:: self.costing.base\_cost = \exp{(\alpha_{1} + \alpha_{2}*\log{mechanical_{work}})}

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
