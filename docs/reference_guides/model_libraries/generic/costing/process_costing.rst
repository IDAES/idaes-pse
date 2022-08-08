Process Costing Library
=======================

.. contents:: Contents 
    :depth: 4

Introduction
------------

.. note:: The power plant costing method is available for most of the unit operations in the IDAES Unit Model Library (Compressor, CSTR, Flash, Heater, HeatExchanger, HeatExchangerNTU, PFR, PressureChanger, Pump, StoichiometricReactor, Turbine, and generic vessels).

The process costing module outlines capital cost methodology for generic IDAES unit operations based on costing correlations that take unit configuration, material, size and other factors for unit costing into consideration.


The framework includes factors for vessel shell thickness and dimensions, material properties for equipment shells and tubes, and tray types for distillation columns. The framework assumes a 2018 cost basis, and does not yet support heat exchangers requiring both shell and tube area (HX1D).

When an instance of a unit model is passed to a `UnitModelCostingBlock`, a new sub-block is created 
on that unit named `costing` (i.e. `flowsheet.unit.costing`). All variables and constraints related to costing will be 
constructed within this new block (see detailed documentation for each unit for details on these variables and constraints).

Note that a `FlowsheetCostingBlock` must already exist, and users must also pass the correct costing method from the SSLW module as an argument (see the example below for details on proper implementation).


Table 1. Main Variables added to the unit block ("self.costing").

=========================== ============================== ============ ===============================================================================
Variable                    Symbol                         Units        Notes
=========================== ============================== ============ ===============================================================================
Purchase cost               :math:`purchase\_cost`         dollars      Purchase cost
Base cost per unit          :math:`base\_cost\_per\_unit`  unitless     Base cost per unit
Base cost                   :math:`base\_cost`             unitless     Base cost (base cost per unit * number of units)
Number of units             :math:`number\_of\_units`      unitless     Number of units to be costed (to take advantage of the economics of scale)
=========================== ============================== ============ ===============================================================================

.. note:: number of units by default is fixed to 1 and the user must unfix this variable to optimize the number of units. Also, `number of units` can be built as a continuous variable or an integer variable. If the latter, the user must provide an mip solver. Use the global costing argument for this purpose (integer_n_units=True or False).

Example
^^^^^^^
Below is a simple example of how to add cost correlations to a flowsheet including a heat exchanger using the default IDAES costing module.


.. code:: python

    from pyomo.environ import ConcreteModel, SolverFactory
    from pyomo.util.calc_var_value import calculate_variable_from_constraint
    from idaes.core import FlowsheetBlock
    from idaes.models.unit_models.heat_exchanger import HeatExchanger, HeatExchangerFlowPattern
    from idaes.models.properties import iapws95
    from idaes.core.util.model_statistics import degrees_of_freedom

    from idaes.models.costing.SSLW import SSLWCosting, SSLWCostingData
    from idaes.core import UnitModelCostingBlock
    from idaes.models.costing.SSLW import HeaterMaterial, HeaterSource
    
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
    
    m.fs.costing = SSLWCosting()
    m.fs.H101.costing = UnitModelCostingBlock(
        default={
            "flowsheet_costing_block": m.fs.costing,
            "costing_method": SSLWCostingData.cost_fired_heater,
            "costing_method_arguments": {
                    "material_type": HeaterMaterial.CarbonSteel,
                    "heat_source": HeaterSource.Fuel
            }
        }
    )

    m.fs.unit.initialize()

    opt = SolverFactory('ipopt')
    opt.options = {'tol': 1e-6, 'max_iter': 50}
    results = opt.solve(m, tee=True)

Units
^^^^^

It is important to highlight that the costing method interrogates the property 
package to determine the units of this model, if the user provided the correct 
units in the metadata dictionary (see property models for additional information), 
the model units will be converted to the right units. 
For example: in this example area is in m^2, while the cost correlations for heat 
exchangers require units to be in ft^2. Therefore, the costing method will convert 
the units to ft^2. The use of Pyomo-unit conversion tools is under development.

Available IDAES Costing Module Methods
--------------------------------------

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
vessels                     :math:`D , L`              ft         
fired heaters               :math:`heat\_duty`         BTU/hr         
=========================== =========================  ===========

Unit Mapping
^^^^^^^^^^^^

As a convenience for users, common unit models are mapped with appropriate costing methods. Therefore, for quick costing with defaults, users may import the `unit_mapping` dictionary and pass the unit model type as below:

.. code:: python

    unit_mapping = {  
            CSTR: cost_vertical_vessel,  
            Compressor: cost_compressor,  
            Flash: cost_vertical_vessel,  
            Heater: cost_fired_heater,  
            HeatExchanger: cost_heat_exchanger,  
            HeatExchangerNTU: cost_heat_exchanger,  
            PFR: cost_horizontal_vessel,  
            PressureChanger: cost_pressure_changer,  
            Pump: cost_pump,  
            StoichiometricReactor: cost_horizontal_vessel,  
            Turbine: cost_turbine,  
        }

    costing_method = unit_mapping[type(blk)]

The mapping is inheritance aware, e.g. pumps will try costing with `cost_pump()` and fall back to `cost_pressure_changer()`.

For example, using default unit mapping the heat exchanger example would add costing as:

.. code:: python

    m.fs.H101.costing = UnitModelCostingBlock(
        default={
            "flowsheet_costing_block": m.fs.costing,
            "costing_method": SSLWCostingData.unit_mapping[HeatExchanger],  # passing the imported class object
            "costing_method_arguments": {
                    "material_type": HeaterMaterial.CarbonSteel,
                    "heat_source": HeaterSource.Fuel
            }
        }
    )

The unit mapping is particularly useful when defining costing for a group of unit models with similar costing arguments and variables, such as vessels, by only passing the model class and not a specific costing method. In the code below, a flowsheet containing two `Flash` vessels, a `CSTR` unit and a `StoichiometricReactor` may be costed compactly via a loop:

.. code:: python

    from idaes.models.costing.SSLW import VesselMaterial
    from idaes.core.util.constants import Constants
    from pyomo.environ import Var, Constraint, units as pyunits, Param, value
    from idaes.models.unit_models import StoichiometricReactor, Flash, CSTR

    # Map unit models to unit classes
    # Will pass to unit_mapping which calls costing methods based on unit class
    unit_class_mapping = {
        m.fs.R101: StoichiometricReactor,
        m.fs.F101: Flash,
        m.fs.F102: Flash,
        m.fs.R102: CSTR}

    # Costing for vessels - m.fs.R101, m.fs.F101, m.fs.F102, m.fs.R102

    # Loop over units
    for unit in [m.fs.R101, m.fs.F101, m.fs.F102, m.fs.R102]:
        # Get correct unit class for unit model
        unit_class = unit_class_mapping[unit]
        
        # Add dimension variables and constraint if they don't exist
        if not hasattr(unit, "diameter"):
            unit.diameter = Var(initialize=1, units=pyunits.m)
        if not hasattr(unit, "length"):
            unit.length = Var(initialize=1, units=pyunits.m)
        if hasattr(unit, "volume"):  # if volume exists, set diameter from volume
            unit.volume_eq = Constraint(expr=unit.volume[0] == unit.length * unit.diameter**2 * 0.25 * Constants.pi)
        else:  # fix diameter directly
            unit.diameter.fix(0.2214 * pyunits.m)
        # Either way, fix L/D to calculate L from D
        # chose constant L/D of 3 for this example - this could also be a dictionary indexed by vessel/unit model
        unit.L_over_D = Constraint(expr=unit.length == 3 * unit.diameter)
            
        # Define vessel costing
        unit.costing = UnitModelCostingBlock(
            default={
                "flowsheet_costing_block": unit.parent_block().costing,  # e.g. m.fs.R101.costing
                "costing_method": SSLWCostingData.unit_mapping[unit_class],  # e.g. cost_vertical_vessel()
                "costing_method_arguments": {
                    "material_type": VesselMaterial.CarbonSteel,
                    "shell_thickness": 1.25 * pyunits.inch,
                    
                }
            }
        )

Heat Exchanger Cost
^^^^^^^^^^^^^^^^^^^

The purchse cost is computed based on the base unit cost and three correction factors (Eq. 22.43 in Seider et al.). The base cost is computed depending on the heat exchanger type selected by the user:

.. math:: self.costing.purchase\_cost = pressure\_factor*material\_factor*L\_factor*self.costing.base\_cost*(CE_{index}/500)

.. math:: self.costing.base\_cost\_per\_unit = \exp{(\alpha_{1} - \alpha_{2}*\log{area*hx\_os} + \alpha_{3}*(\log{area*hx\_os})^{2})}

.. math:: self.costing.base\_cost = self.costing.base\_cost\_per\_unit * self.costing.number\_of\_units

.. math:: area  = self.area / self.costing.number\_of\_units

where:

* pressure_factor - is the pressure design correction factor
* material_factor - is the construction material correction factor
* length_factor - is the tube length correction factor
* CE_index - is a global parameter for Chemical Enginering cost index for years 2010-2019
* hx_os - heat exchanger oversize factor (default = 1)
* area is a reference object and (self.area is the model variable)

The heat exchanger costing method has three arguments, hx_type = heat exchanger type, material_type = construction material, and tube_lengh = tube length.

* hx_type : 'floating_head', 'fixed_head', 'Utube'\*, 'Kettle_vap'
* material_type: 'CarbonSteelCarbonSteel', 'CarbonSteelBrass', 'CarbonSteelStainlessSteel'\*, 'CarbonSteelMonel', 'CarbonSteelTitanium', 'CarbonSteelCrMoSteel', 'CrMoSteelCrMoSteel', 'MonelMonel', 'TitaniumTitanium',
* tube_length: '8ft', '12ft'\*, '16ft', '20ft'

where '*' corresponds to the default options. Additionally, users may pass an argument 'integer' (defaults to `True``) whether the number of units should be restricted to an integer or not.


Table 3. Base cost factors for heat exchanger type.

================= ================== ================== ==================
HX Type           :math:`\alpha_{1}` :math:`\alpha_{2}` :math:`\alpha_{3}`
================= ================== ================== ==================
floating_head     11.9052            0.8709             0.09005 
fixed_head        11.2927            0.8228             0.09861
Utube             11.3852            0.9186             0.09790
Kettle_vap        12.2052            0.8709             0.09005
================= ================== ================== ==================


Table 4. Tube-Length correction factor.

===================    =====
HX Tube Length (ft)    FL
===================    =====
8                      1.25
12                     1.12
16                     1.05
20                     1.00
===================    =====

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

Pressure Changer Cost
^^^^^^^^^^^^^^^^^^^^^

The costing of a pressure changer unit model is more complicated, because the pressure changer 
model can be imported into the flowsheet object representing a pump, turbine, compressor, or a 
simply pressure changer (fan, blower, etc.). The `cost_pressure_changer` method currently supports costing of pumps, turbines, and compressors. The method authomatically interrogates the flowsheet object to determine if the unit is being used as a pump, turbine, or compressor. 

The `cost_pressure_changer()` method authomatically determines if the unit model is being used as a pump, 
turbine, or compressor based on the `compressor` and `thermodynamic_assumption` configuration 
arguments provided by the user where creating the unit model. A summary of the decision logic is shown below.


========== =========== =================================
Unit Type  compressor  thermodynamic_assumption
========== =========== =================================
Turbine    False       Any
Pump       True        pump
Mover      True        not pump, isothermal or adiabatic
========== =========== =================================

Additionally, some unit types have different sub-types which can be costed appropiately. In these cases, 
an additional argument is provided to `cost_pressure_changer` to identify the sub-type to use which is detailed below.

Turbine Cost Model
""""""""""""""""""
The turbine cost is based on the mechanical work of unit (work_mechanical), this correlation has been obtained using the NETL Report (DOE/NETL 2015).

.. math:: self.costing.purchase\_cost = 580*(mechanical_{work})^{0.81}

DOE/NETL, 2015, report. Cost and performance Baseline for Fossil Energy Plants. Volume 1a: Bituminous Coal (PC) and Natural Gas to Electricity. Revision 3

Users may pass an argument 'integer' (defaults to True) whether the number of units should be restricted to an integer or not.

Pump Cost Model
""""""""""""""""

Three subtypes are supported for costing of pumps, which can be set using the "pump_type" argument.

1) Centrifugal pumps (pump_type='Centrifugal')\*
2) External gear pumps (pump_type='ExternalGear')
3) Reciprocating Plunger pumps (pump_type='Reciprocating')

Additionally, there are an array of additional pump options:

* material_type: 'CastIron', 'DuctileIron', 'CastSteel', 'Bronze', 'StainlessSteel'\*, 'HastelloyC', 'Monel', 'Nickel', 'Titanium', 'NiAlBronze', 'CarbonSteel'

In the lists above, \* marks the default options. Additionally, users may pass an argument 'integer' (defaults to True) whether the number of units should be restricted to an integer or not.

Centrifugal Pump
++++++++++++++++

The centrifugal pump cost has two main components, the cost of the pump and the cost of the motor. The pump cost is based on the fluid work (work_fluid), pump head, and size factor. 
Additional arguments are required:

* pump_type_factor = empirical factor (see table 6) of 1.1, 1.2, 1.3, 1.4\*, 2.1, 2.2
* pump_motor_type_factor = 'open'\*, 'enclosed', 'explosion_proof'

In the lists above, \* marks the default options.

Based on user's inputs the costing method builds base_cost and purchase_cost for both the pump and the motor. 
The unit purchase cost is obtained by adding the motor and pump costs.

.. math:: self.costing.purchase\_cost = self.costing.pump\_purchase\_cost + self.costing.motor\_purchase\_cost

To compute the purchase cost of the centrifugal pump, first we obtain the pump size factor (S) with Eq. 22.13, then we obtain the base cost with Eq. 22.14.
Finally, the purchase cost of the pump is obtained in Eq. 22.15. (Seider et al.)

.. math:: S = QH^{0.5}

.. math:: self.costing.pump\_base\_cost\_per\_unit = \exp{(9.7171 - 0.6019*\log{S} + 0.0519*(\log{S})^{2})}

.. math:: self.costing.pump\_purchase\_cost = F_{T}*material\_factor*self.costing.pump\_base\_cost*(CE_{index}/500)

.. math:: self.costing.base\_cost = self.costing.pump\_base\_cost\_per\_unit * self.costing.number\_of\_units

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

Electric Motor Pumps
++++++++++++++++++++

A centrifugal pump is usually driven by an electric motor, the `self.costing.motor_purchase_cost` is calculated based on the power consumption.

.. math:: self.motor_purchase_cost = FT * self.costing.motor\_base\_cost * (CE_{index}/500)  (Eq. 22.20)

.. math:: self.costing.motor\_base\_cost = self.costing.motor\_base\_cost\_per\_unit * self.costing.number\_of\_units

.. math:: Q  = self.Q / self.costing.number\_of\_units

.. math:: self.costing.self.costing.motor\_base\_cost\_per\_unit = \exp{(5.8259 + 0.13141\log{PC} + 0.053255(\log{PC})^{2} + 0.028628(\log{PC})^{3} - 0.0035549(\log{PC})^{4})}  (Eq. 22.19)

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

.. math:: self.costing.pump\_base\_cost = self.costing.pump\_base\_cost\_per\_unit * self.costing.number\_of\_units

.. math:: self.costing.self.costing.pump\_base\_cost\_per\_unit  = \exp{(7.6964 + 0.1986\log{Q} + 0.0291(\log{Q})^{2})}

.. math:: Q  = self.Q / self.costing.number\_of\_units

Reciprocating Plunger Pumps
+++++++++++++++++++++++++++

The cost correlation method used here is based on the brake horsepower (PB).

.. math:: self.costing.pump\_purchase\_cost = material\_factor * self.costing.pump\_base\_cost * (CE_{index}/500)  (Eq. 22.22)

.. math:: self.costing.pump\_base\_cost = self.costing.pump\_base\_cost\_per\_unit * self.costing.number\_of\_units

.. math:: self.costing.pump\_base\_cost\_per\_unit = \exp{(7.8103 + 0.26986\log{PB} + 0.06718(\log{PB})^{2})} (Eq. 22.23)

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

* mover_type= 'compressor'\* or 'fan' or 'blower' (upper/lower case sensitive)

In the list above, \* marks the default option.

Compressor Cost
+++++++++++++++
The compressor cost is based on the mechanical work of the unit. 
Additional arguments are required to estimate the cost such as compressor type, 
driver mover type, and materials of construction (material_type).

* compressor_type = 'Centrifugal'\*, 'Reciprocating', 'Screw'
* driver_mover_type = 'ElectricMotor'\*, 'SteamTurbine', 'GasTurbine'
* material_type = 'CarbonSteel', 'StainlessSteel'\*, 'NickelAlloy'

where '*' corresponds to the default options. Additionally, users may pass an argument 'integer' (defaults to True) whether the number of units should be restricted to an integer or not.

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
fan_type, and material_type.

* fan_type = 'CentrifugalBackward'\*, 'CentrifugalStraight', 'VaneAxial', 'TubeAxial'
* fan_head_factor = defaults to 1.45, see table 14
* material_type = 'CarbonSteel', 'Fiberglass', 'StainlessSteel'\*, 'NickelAlloy'

where '*' corresponds to the default options. Additionally, users may pass an argument 'integer' (defaults to True) whether the number of units should be restricted to an integer or not.

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
Additional arguments are required to estimate the blower cost such as mover_type='blower', blower_type, and materials of construction (material_type).

* blower_type = 'Centrifugal'\*, 'Rotary'
* material_type = 'CarbonSteel', 'Aluminum', 'Fiberglass', 'StainlessSteel'\*, 'NickelAlloy'

where '*' corresponds to the default options. Additionally, users may pass an argument 'integer' (defaults to True) whether the number of units should be restricted to an integer or not.

The material factors given in table 15 for the fans can be used. In addition, centrifugal blowers are available with cast aluminum blades with Mat_factor = 0.60.

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
This method computes the purchase cost of the fired heater based on the heat duty, fuel used (fired_type), pressure design, and materials of construction (material_type).

* fuel_type = 'Fuel'\*, 'Reformer', 'Pyrolysis', 'HotWater', 'Salts', 'DowthermA', 'SteamBoiler'
* material_type = 'CarbonSteel'\*, 'CrMoSteel', 'StainlessSteel'

where '*' corresponds to the default options. Additionally, users may pass an argument 'integer' (defaults to True) whether the number of units should be restricted to an integer or not.

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
This method constructs by defaul the cost of pressure vessels with platforms and ladders, and trays cost can be calculated if trays= `True`. This method requires a few arguments to build the cost of vessel. 
We recommend using this method to cost reactors (`StoichiometricReactor`, `CSTR`` or `PFR``), flash tanks, vessels, and distillation columns.

* vertical = `False`, `True`
* material_type = 'CarbonSteel'\*, 'LowAlloySteel', 'StainlessSteel304', 'StainlessSteel316', 'Carpenter20CB3', 'Nickel200', 'Money400', 'Inconel600', 'Incoloy825', 'Titanium'
* shell_thickness = material property  including pressure allowance, defaults to 1.25 inches (Pyomo units)
* weight_limit = 1, 2 (1: 1000 to 920,000 lb, 2: 9000 to 2.5M lb only for vertical vessels)
* aspect_ratio_range = 1, 2 (1: 3 < D < 21, 12 < L < 40; 2: 3 < D < 24, 27 < L < 170; all in ft D: diameter, L: length) only for vertical vessels
* include_platform_ladders= `True`, `False` : to build platforms and ladders cost
* vessel_diameter = component to use as diameter variable, defaults to unit.diameter
* vessel_length = component to use as length variable, defaults to unit.length
* number_of_units = integer or Pyomo component to use for number of parallel units to be costed, defaults to 1
* number_of_trays = component to use for number of distillation trays in vessel, defaults to `None`
* tray_material = 'CarbonSteel', 'StainlessSteel303', 'StainlessSteel316', 'Carpenter20CB3', 'Monel'
* tray_type = 'Sieve', 'Valve', 'BubbleCap'

where '*' corresponds to the default options.

By adding reference parameter, the method can be constructed in any pyomo costing block.
Since the generic models do not include the variables required to cost these type of units, the user must create the blocks and variables.

For example:

.. code:: python

    m.fs.unit = Block()
    m.fs.unit.diameter = Var()
    m.fs.unit.length = Var().
    m.fs.unit.costing = pyo.Block()
    vessel_costing_method = SSLWCostingData.costing_vessel(m.fs.unit.costing, args).

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
Horizontal vessels (1: 1000 < W < 920,000 lb):

.. math:: self.costing.base\_cost\_per\_unit = \exp{(8.9552 - 0.2330*\log{weight} + 0.04333*(\log{weight})^{2})}

Vertical vessels (1: 4200 < W < 1M lb):

.. math:: self.costing.base\_cost\_per\_unit = \exp{(8.9552 - 0.2330*\log{weight} + 0.04333*(\log{weight})^{2})}

Vertical vessels (2: 9,000 < W < 2.5M lb):

.. math:: self.costing.base\_cost\_per\_unit = \exp{(7.2756 - 0.18255*\log{weight} + 0.02297*(\log{weight})^{2})}

.. math:: self.costing.base\_cost = self.costing.base\_cost\_per\_unit * self.costing.number\_of\_units

.. math:: weight = self.weight / self.costing.number\_of\_units

The vessel purchase cost is given by:

.. math:: self.vessel\_purchase\_cost = (CE_{index}/500) * material\_factor * self.base\_cost + (self.base\_cost\_platf\_ladders * self.costing.number\_of\_units)

note that if PL = `False`, the cost of platforms and ladders is not included.

The final purchase cost is given by:

.. math:: self.purchase\_cost = self.vessel\_purchase\_cost + (self.purchase\_cost\_trays * self.costing.number\_of\_units)

note that if plates=`False`, the cost of trays is not included.


Base Cost of Platforms and ladders
""""""""""""""""""""""""""""""""""
The cost of platforms and ladders is based on the diamter and length in ft.
Horizontal vessels (1: 3 < D < 12 ft):

.. math:: self.base\_cost\_platf\_ladders = 20059*D^{0.20294}

Vertical vessels (1: 3 < D < 12 ft and 12 < L  < 40 ft):

.. math:: self.base\_cost\_platf\_ladders = 361.8*D^{0.73960} * L^{0.70684}

Vertical vessels (2: 3 < D < 24 ft and 27 < L  < 170 ft):

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

Module Classes
^^^^^^^^^^^^^^^^^

.. module:: idaes.models.costing.SSLW

.. autoclass:: SSLWCosting
    :members:

.. autoclass:: SSLWCostingData
    :members:
