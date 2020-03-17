.. index::
    pair: idaes.core.flowsheet_model;FlowsheetBlockData
    pair: idaes.core.flowsheet_model;FlowsheetBlock

Unit Model Costing
==================

.. contents:: Contents 
    :depth: 2

A generic costing (get_costing():) method has been developed for all the units in the IDAES unit model library. All cost correlations have been obtained from Seider et al.

Process and Product Design Principles: Synthesis, Analysis, and Evaluation
Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
Chapter 22. Cost Accounting and Capital Cost Estimation
22.2 Cost Indexes and Capital Investment

The get_costing method will create a block named "costing" that includes the cost correlations and variables. 
The main variables are purchase cost and base cost, each unit has a different basis to estimate the cost. 
    


Heat Exchanger Cost
-------------------

Heat exchanger costing method is based on the area of the heat exchanger. 
This method computes the purchase cost (CP) for a shell and tube heat 
exchanger (Eq. 22.43), the model computes the base cost (CB for 4 type
of heat exchangers, such as floating head, fixed head, U-tube, and
Kettle vaporizer), construction material factor (FM), pressure design
factor (FP), and tube length correction factor (FL), using CE base cost
index of 500.

Cp = FP*FM*FL*CB*(CE_index/500)

The heat exchanger costing method has three arguments, h_type, material factor, and tube lenght. 

hx_type : 'floating_head', 'fixed_head', 'U-tube'*, 'Kettle_vap'

material factor (FM): 'stain_steel'*, 'carb_steel'

tube length (FL): '8ft'*, '12ft', '16ft', '20ft'

* --> default options

Pressure Changer Cost
---------------------
The pressure changer unit model is more complicated, because the pressure changer model can be imported to represent a pump, turbine, compressor, or a simply pressure changer (fan, blower, etc.).
The get_costing(): method currently supports costing of pumps, turbines, and compressors. The method authomatically interrogates the flowsheet object to determine if the unit is being used as a pump, turbine, or compressor. 
The additional arguments are required to build correlations for different type of pumps or compressors. 

Turbine: 
Turbine cost is based on the mechanical work of unit (work_mechanical).
We determine if the pressure changer is a turbine if "config.compressor == False".

Pump:
The pump cost has two main components, the cost of the pump and the cost of the motor. The pump cost is based on the fluid work (work_fluid), pump head, and size factor. 
We determine if the pressure change is a pump if " unit.config.compressor == True, and unit.config.Thermodynamic.assumption.name == 'pump' "
Additional arguments are required:
pump_type = 'centrifugal'
pump_type_factor = '1.4'
pump_motor_type_factor = 'open'

Mover:
We determine if the pressure changer is a mover if "unit.compressor == True and unit.config.Thermodynamic.assumption.name not 'pump ". 
However, the pressure changer could be a compressor, fan, or blower. 
Therefore, the user must set the mover_type = 'compressor', 'fan', 'blower'

The compressor cost is based on the mechanical work of the unit. Additional arguments are required for the compressor type, driver mover type, and material factor (FM_mat).
* compressor_type = 'centrifugal', 'reciprocating', 'screw'
* driver_mover_type = 'electrical_motor', 'steam_turbine', 'gas_turbine'
* FD_mat = 'carbon_steel', 'stain_steel', 'nickel_alloy' 


Flowsheet Classes
-----------------

.. module:: idaes.core.flowsheet_model

.. autoclass:: FlowsheetBlockData
    :members:

.. autoclass:: FlowsheetBlock
    :members:
