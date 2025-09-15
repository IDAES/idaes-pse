# -*- coding: utf-8 -*-
"""
Economics Calculations
MEA Solvent Based Carbon Capture System integrated with NGCC Power Plant 
(point source of flue gas)
"""

import idaes
from math import pi
import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
    )
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import units as pyunits
from idaes.core.solvers import get_solver
import os
from pyomo.common.fileutils import this_file_dir
import csv

def get_ngcc_solvent_cost(self, args_ngcc, args_solvent, args_compr,
                          export_economic_results=True,
                          overwrite_economic_results=True):
# -------------------------------------------------------------------------------------------
    # Create a main Costing block
    self.costing = QGESSCosting()
# -------------------------------------------------------------------------------------------
    # Access the scaling parameter values from the arguments
    # These objects are aliases for corresponding variables in the IDAES model
    feedwater_flowrate = pyo.units.convert(args_ngcc[0], to_units=pyo.units.lb/pyo.units.hr)  # lb/hr
    fuelgas_value = pyo.units.convert(args_ngcc[1], to_units=pyo.units.lb/pyo.units.hr)  # lb/hr
    fluegas_value = pyo.units.convert(args_ngcc[2], to_units=pyo.units.ft**3/pyo.units.min)  # ft3/min
    stack_flow_gas = pyo.units.convert(args_ngcc[3], to_units=pyo.units.ft**3/pyo.units.min)  # ft3/min

    # Total amount of CO2 captured in lb/hr (all trains combined)
    CO2CaptureRate = pyo.units.convert(args_solvent[0], to_units=pyo.units.lb/pyo.units.hr)  # lb/hr
    # Number of trains for the capture system
    NoTrain = pyo.units.convert(args_solvent[1], to_units=pyo.units.dimensionless)
    # Height of Absorber # m
    Abs_ht = pyo.units.convert(args_solvent[2], to_units=pyo.units.m)
    # Diameter of Absorber # m
    Abs_d = pyo.units.convert(args_solvent[3], to_units=pyo.units.m)
    # Height of Stripper # m
    Strip_ht = pyo.units.convert(args_solvent[4], to_units=pyo.units.m)
    # Diameter of Stripper # m
    Strip_d = pyo.units.convert(args_solvent[5], to_units=pyo.units.m)
    # solvent Makeup  # kg/tonne CO2
    solvent_makeup = pyo.units.convert(args_solvent[6], to_units=pyo.units.kg/pyo.units.tonne)
    # CO2 Emissions without capture lb/hr
    CO2_emission_nocap = pyo.units.convert(args_solvent[7], to_units=pyo.units.lb/pyo.units.hr)
    # CO2 Emissions with capture lb/hr
    CO2_emission_cap = pyo.units.convert(args_solvent[8], to_units=pyo.units.lb/pyo.units.hr)
    # solvent initial fill # kg
    solvent_fill_init = pyo.units.convert(args_solvent[9], to_units=pyo.units.kg)

    # Stripper Reboiler Duty (MW)
    Qreb = pyo.units.convert(args_solvent[10], to_units=pyo.units.MW)
    # Lean rich heat exchanger area (m2)
    LRHEXArea = pyo.units.convert(args_solvent[11], to_units=pyo.units.m**2)
    # Absorber intercooler details (temperature in deg K, heat duty in W)
    # ICOOL1Tin = pyo.units.convert(args_solvent[12], to_units=pyo.units.K)
    # ICOOL1Tout = pyo.units.convert(args_solvent[13], to_units=pyo.units.K)
    # QICOOL1 = pyo.units.convert(args_solvent[14], to_units=pyo.units.W)

    # ICOOL2Tin = pyo.units.convert(args_solvent[15], to_units=pyo.units.K)
    # ICOOL2Tout = pyo.units.convert(args_solvent[16], to_units=pyo.units.K)
    # QICOOL2 = pyo.units.convert(args_solvent[17], to_units=pyo.units.W)

    # Lean solvent cooler details (temperature in deg K, heat duty in W)
    # LSCoolerTin = pyo.units.convert(args_solvent[18], to_units=pyo.units.K)
    # LSCoolerTout = pyo.units.convert(args_solvent[19], to_units=pyo.units.K)
    # QLSCooler = pyo.units.convert(args_solvent[20], to_units=pyo.units.W)

    # Stripper condenser details (temperature in deg K, heat duty in W)
    StripcondTin = pyo.units.convert(args_solvent[21], to_units=pyo.units.K)
    StripcondTout = pyo.units.convert(args_solvent[22], to_units=pyo.units.K)
    QStripcond = pyo.units.convert(args_solvent[23], to_units=pyo.units.W)

    # Stripper outlet from reboiler temperature (deg K)
    Rebtempout = pyo.units.convert(args_solvent[24], to_units=pyo.units.K)

    # Flue gas blower net power requirement (kW)
    FGBlowerLoad = pyo.units.convert(args_solvent[25], to_units=pyo.units.kW)

    # Rich solvent pump net power requirement (kW)
    RichSolvPumpLoad = pyo.units.convert(args_solvent[26], to_units=pyo.units.kW)

    # Lean solvent pump net power requirement (kW)
    LeanSolvPumpLoad = pyo.units.convert(args_solvent[27], to_units=pyo.units.kW)

    # Fraction of steam extracted from the IP/LP steam turbine crossover - NGCC
    steam_extraction_ngcc = pyo.units.convert(args_solvent[28], to_units=pyo.units.dimensionless) # MW/MW

    # Packing height of the wash section in absorber column (m)
    wash_section_packing_ht = pyo.units.convert(args_solvent[29], to_units=pyo.units.m)

    # DCC Column 1 packing height (m)
    DCC_column1_packing_ht = pyo.units.convert(args_solvent[30], to_units=pyo.units.m)

    # DCC Column 1 diameter (m)
    DCC_column1_diameter = pyo.units.convert(args_solvent[31], to_units=pyo.units.m)

    # DCC Cooler 1 duty (W)
    DCC_cooler1_duty = pyo.units.convert(args_solvent[32], to_units=pyo.units.W)

    # DCC Cooler 1 temperature stream in (deg K)
    DCC_cooler1_streamtempin = pyo.units.convert(args_solvent[33], to_units=pyo.units.K)

    # DCC Cooler 1 temperature stream out (deg K)
    DCC_cooler1_streamtempout = pyo.units.convert(args_solvent[34], to_units=pyo.units.K)

    # DCC Pump 1 Load (kW)
    DCC_pump1_load = pyo.units.convert(args_solvent[35], to_units=pyo.units.kW)

    # DCC Pump 1 Flowrate (m3/s)
    DCC_pump1_flowrate = pyo.units.convert(args_solvent[36], to_units=pyo.units.m**3/pyo.units.s)

    # Rich solvent flowrate (m3/s)
    rich_solvent_flowrate = pyo.units.convert(args_solvent[37], to_units=pyo.units.m**3/pyo.units.s)

    # Lean solvent flowrate (m3/s)
    lean_solvent_flowrate = pyo.units.convert(args_solvent[38], to_units=pyo.units.m**3/pyo.units.s)

    # Absorber off gas flowrate (kg/hr)
    absorber_off_gas_flowrate = pyo.units.convert(args_solvent[39], to_units=pyo.units.kg/pyo.units.hr)

    # Flue gas flowrate at blower discharge (m3/s)
    fg_flowrate_blowerout = pyo.units.convert(args_solvent[40], to_units=pyo.units.m**3/pyo.units.s)

    # Pumparound 1 flowrate (m3/s)
    # pumparound_1_flowrate = pyo.units.convert(args_solvent[41], to_units=pyo.units.m**3/pyo.units.s)

    # Pumparound 2 flowrate (m3/s)
    # pumparound_2_flowrate = pyo.units.convert(args_solvent[42], to_units=pyo.units.m**3/pyo.units.s)
# -------------------------------------------------------------------------------------------
    # Absorber wash water pump net power requirement (kW)
    # Assume it scales linearly with flue gas flowrate
    self.costing.FG_flowrate_per_train_ref = pyo.Param(initialize=930909, mutable=False,
                                          units=pyo.units.ft**3/pyo.units.min) # acfm
    self.costing.WWashLoadRef = pyo.Param(initialize=71, mutable=False,
                             units=pyo.units.kW) # kW
    self.costing.WWashLoad = pyo.Expression(
        expr=pyo.units.convert(
            (fluegas_value/NoTrain)*self.costing.WWashLoadRef/self.costing.FG_flowrate_per_train_ref,
            to_units=pyo.units.kW)
        )
    
    # Stripper reflux pump net power requirement (kW)
    # Assume it scales linearly with CO2 capture flowrate
    self.costing.CO2CaptureRate_ref = pyo.Param(initialize=290676, mutable=False,
                                   units=pyo.units.kg/pyo.units.hr) # kg/hr
    self.costing.Stripperrefluxload_ref = pyo.Param(initialize=21.277, mutable=False,
                                       units=pyo.units.kW) # kW
    self.costing.Stripperrefluxload = pyo.Expression(
        expr=pyo.units.convert(
            (CO2CaptureRate/NoTrain)*(self.costing.Stripperrefluxload_ref/self.costing.CO2CaptureRate_ref),
            to_units=pyo.units.kW)
        )
    
    # Lean solvent pump 2 (after the solvent makeup is added) net power requirement (kW)
    # Assume it scales linearly with lean solvent flowrate
    self.costing.LSflow_ref = pyo.Param(initialize=21166, mutable=False,
                           units=pyo.units.gal/pyo.units.min) # gpm
    self.costing.LSpumpload_ref = pyo.Param(initialize=1071.35, mutable=False,
                               units=pyo.units.kW) # kW
    self.costing.LSpumpload = pyo.Expression(
        expr=pyo.units.convert(
            (lean_solvent_flowrate)*self.costing.LSpumpload_ref/self.costing.LSflow_ref,
            to_units=pyo.units.kW)
        )
    
    # Wash water cooler heat transfer surface area (m2)
    # Assume it scales linearly with flue gas flowrate
    self.costing.wwcooler_htsa_ref = pyo.Param(initialize=20, mutable=False,
                                  units=pyo.units.m**2) # m2
    self.costing.wwcooler_htsa = pyo.Expression(
        expr=pyo.units.convert(
            (fluegas_value/NoTrain)*self.costing.wwcooler_htsa_ref/self.costing.FG_flowrate_per_train_ref,
            to_units=pyo.units.m**2)
        )
       
    # Auxilliary load from CO2 compressors # hp
    CompAuxLoad_Cap = pyo.units.convert(args_compr[0], to_units=pyo.units.hp)
    # Heat Duty from compressor intercoolers # converting to W, original file used Btu/hr
    CompIntercool_Cap = pyo.units.convert(args_compr[1], to_units=pyo.units.W)
    
    # Additional Calculations
    # Low pressure steam condenser duty - NGCC
    self.costing.condenser_duty = pyo.Expression(
        expr=((1301.66 - (1139.70*steam_extraction_ngcc))*pyo.units.MBtu/pyo.units.hr)
        ) # MMBtu/hr, units from surrogate model
    # HRSG duty - NGCC
    self.costing.HRSG_duty = pyo.Expression(
        expr=((558.83 - (58.484*steam_extraction_ngcc))*3.96832*pyo.units.MBtu/pyo.units.hr)
        )# MMBtu/hr, units from surrogate model
    # Low pressure steam turbine power before generator losses
    self.costing.LPST_power = pyo.Expression(
        expr=((119.564 - (106.134*steam_extraction_ngcc))*pyo.units.MW)
        ) # MWe, units from surrogate model
    # High pressure steam turbine power before generator losses
    self.costing.HPST_power = pyo.Param(initialize=54.54, mutable=False, units=pyo.units.MW) # MWe
    # Intermediate pressure steam turbine power before generator losses
    self.costing.IPST_power = pyo.Param(initialize=85.51, mutable=False, units=pyo.units.MW) # MWe
    # Steam turbine efficiency
    self.costing.turbine_eff = pyo.Param(initialize=0.987186, mutable=False, units=pyo.units.dimensionless)
    # Aggregate power from steam turbine after generator losses
    self.costing.ST_gross_power = pyo.Expression(
        expr=self.costing.turbine_eff*(
            pyo.units.convert(self.costing.LPST_power, to_units=pyo.units.kW) +
            pyo.units.convert(self.costing.IPST_power, to_units=pyo.units.kW) +
            pyo.units.convert(self.costing.HPST_power, to_units=pyo.units.kW))
        ) # kW

    # Power from combustion turbine after generator losses
    self.costing.CT_gross_power = pyo.Expression(
        expr=pyo.units.convert(
            (469.56*pyo.units.MW),
            to_units=pyo.units.kW)
        ) # kW
    # Gas turbine power
    self.costing.gt_power = pyo.Expression(expr=(self.costing.CT_gross_power)) # kW
    # Plant gross power   
    self.costing.plant_gross_power = pyo.Expression(
        expr=pyo.units.convert(self.costing.CT_gross_power +
                               self.costing.ST_gross_power,
                               to_units=pyo.units.kW)
        ) # kW
    # STG CTG output
    self.costing.ctg_op = pyo.Param(initialize=540000, mutable=False, units=pyo.units.kW) # kW, same for case B31A and case B31B
    self.costing.stg_op_ref = pyo.Param(initialize=200, mutable=False, units=pyo.units.MW) # MW, case B31B
    self.costing.ST_gross_power_ref = pyo.Param(initialize=213, mutable=False, units=pyo.units.MW) # MW, case B31B
    self.costing.stg_op = pyo.Expression(
        expr=pyo.units.convert(
            self.costing.ST_gross_power*(self.costing.stg_op_ref/self.costing.ST_gross_power_ref),
            to_units=pyo.units.kW)
        ) # kW
    self.costing.stg_ctg_op = pyo.Expression(
        expr=pyo.units.convert(
            self.costing.ctg_op + self.costing.stg_op,
            to_units=pyo.units.kW)
        ) # kW
    
    # Cooling tower duty calculation
    # Aggregate cooling duty requirement in CCS, CO2 compression
    # Negative sign is used to adjust the sign convention and get a positive value
    self.costing.cooling_duty_CCS_compr = pyo.Expression(
        expr=pyo.units.convert(
            -NoTrain*(DCC_cooler1_duty + QStripcond) + CompIntercool_Cap,  # removed QICOOL1, QICOOL2, QLSCooler
            to_units=pyo.units.W)
        ) # W
    # Auxilliary cooling system requirement
    self.costing.aux_cooling_load = pyo.Param(initialize=146, mutable=False, units=pyo.units.GJ/pyo.units.hr) # GJ/hr
    # Cooling tower duty
    self.costing.cooling_tower_duty = pyo.Expression(
        expr=(
            pyo.units.convert(self.costing.cooling_duty_CCS_compr, to_units=pyo.units.MBtu/pyo.units.hr) +
            pyo.units.convert(self.costing.aux_cooling_load, to_units=pyo.units.MBtu/pyo.units.hr) +
            pyo.units.convert(self.costing.condenser_duty, to_units=pyo.units.MBtu/pyo.units.hr)
            )
        ) # MMBtu/hr
        
    # Circulating water flowrate calculation
    # Difference between circulating water return and supply temperatures
    # (cooling tower design specification)
    self.costing.dTCW = pyo.Param(initialize=11, mutable=False, units=pyo.units.K) # K
    # Average cp of water between 16 deg C and 27 deg C
    self.costing.cpwater_avg = pyo.Param(initialize=4.18147, mutable=False, units=pyo.units.kJ/pyo.units.kg/pyo.units.K) # kJ/kg.K
    # Density of water
    self.costing.rhowater = pyo.Param(initialize=996.38, mutable=False, units=pyo.units.kg/pyo.units.m**3) # kg/m3
    # Circulating water flowrate
    self.costing.cir_water_flowrate = pyo.Expression(
        expr=pyo.units.convert(
            self.costing.cooling_tower_duty/(self.costing.dTCW*self.costing.cpwater_avg*self.costing.rhowater),
            to_units=pyo.units.gal/pyo.units.min)
        ) # gpm
        
    # Raw water withdrawal flowrate calculation
    # Considered approximately equal to cooling tower makeup requirement
    self.costing.evap_losses = pyo.Expression(expr=0.008*self.costing.cir_water_flowrate*self.costing.dTCW/5.5/pyo.units.K) # gpm
    self.costing.drift_losses = pyo.Expression(expr=0.00001*self.costing.cir_water_flowrate) # gpm
    self.costing.CC = pyo.Param(initialize=4, mutable=False, units=pyo.units.dimensionless)
    self.costing.blowdown_losses = pyo.Expression(expr=self.costing.evap_losses/(self.costing.CC - 1)) # gpm
    self.costing.raw_water_withdrawal = pyo.Expression(expr=self.costing.evap_losses + self.costing.drift_losses + self.costing.blowdown_losses) # gpm
    
    # Process water discharge
    self.costing.process_water_discharge = pyo.Expression(expr=0.35*self.costing.raw_water_withdrawal) # gpm
    
    # Auxilliary load calculation
    # Circulating water flowrate case B31B
    self.costing.circ_water_flow_ref = pyo.Param(initialize=223629.5473, mutable=False,
                                    units=pyo.units.gal/pyo.units.min) # gpm
    # Circulating water pumps power requirement case B31B
    self.costing.circ_water_pumps_load_ref = pyo.Param(initialize=4580, mutable=False, units=pyo.units.kW) # kWe
    # Circulating water pumps power requirement - scaled
    self.costing.circ_water_pumps_load = pyo.Expression(
        expr=self.costing.circ_water_pumps_load_ref*(self.costing.cir_water_flowrate/self.costing.circ_water_flow_ref)
        )
    # Combustion turbine auxilliary load
    self.costing.CT_aux_load = pyo.Param(initialize=1020, mutable=False, units=pyo.units.kW) # kWe
    # Condensate pumps auxilliary load
    self.costing.condensate_pumps_aux_load = pyo.Param(initialize=170, mutable=False, units=pyo.units.kW) # kWe
    # Cooling tower fans power requirement case B31B
    self.costing.cool_tower_fans_load_ref = pyo.Param(initialize=2370, mutable=False, units=pyo.units.kW) # kWe
    # Cooling tower fans power requirement - scaled
    self.costing.cool_tower_fans_load = pyo.Expression(
        expr=self.costing.cool_tower_fans_load_ref*(self.costing.cir_water_flowrate/self.costing.circ_water_flow_ref)
        )
    # CO2 capture system auxilliary load
    self.costing.AuxLoad_CCS = pyo.Expression(
        expr=pyo.units.convert(
            NoTrain*(FGBlowerLoad + self.costing.WWashLoad + RichSolvPumpLoad + LeanSolvPumpLoad
                     + DCC_pump1_load + self.costing.Stripperrefluxload + self.costing.LSpumpload),
            to_units=pyo.units.kW)
        ) # kW
    # CO2 compression system auxilliary load
    self.costing.AuxLoad_compr = pyo.Expression(
        expr=pyo.units.convert(CompAuxLoad_Cap, to_units=pyo.units.kW)
        ) # kW
    # Feedwater pumps load
    self.costing.feedwater_pumps_load = pyo.Param(initialize=4830, mutable=False, units=pyo.units.kW) # kWe
    # Groundwater pumps load
    # Groundwater pumps power requirement case B31B
    self.costing.groundwater_pumps_load_ref = pyo.Param(initialize=430, mutable=False, units=pyo.units.kW) # kWe
    # Raw water withdrawal case B31B
    self.costing.raw_water_withdrawal_ref = pyo.Param(initialize=4773, mutable=False, units=pyo.units.gal/pyo.units.min) # gpm
    self.costing.groundwater_pumps_load = pyo.Expression(
        expr=self.costing.groundwater_pumps_load_ref*self.costing.raw_water_withdrawal/self.costing.raw_water_withdrawal_ref
        ) # kWe
    # SCR load
    self.costing.SCR_load = pyo.Param(initialize=2, mutable=False, units=pyo.units.kW) # kWe
    # Steam turbine auxilliaries 
    self.costing.ST_aux_load = pyo.Param(initialize=200, mutable=False, units=pyo.units.kW) # kWe
    # Transformer losses
    self.costing.transformer_losses = pyo.Param(initialize=2200, mutable=False, units=pyo.units.kW) # kWe
    # Miscellanous load
    self.costing.misc_loads = pyo.Param(initialize=570, mutable=False, units=pyo.units.kW) # kWe
    # Total auxilliary load
    self.costing.aux_load = pyo.Expression(
        expr=pyo.units.convert(
            self.costing.circ_water_pumps_load + self.costing.CT_aux_load +
            self.costing.condensate_pumps_aux_load + self.costing.cool_tower_fans_load +
            self.costing.AuxLoad_CCS + self.costing.AuxLoad_compr +
            self.costing.feedwater_pumps_load + self.costing.groundwater_pumps_load +
            self.costing.SCR_load + self.costing.ST_aux_load +
            self.costing.transformer_losses + self.costing.misc_loads,
            to_units=pyo.units.kW)
        )  # kW
        
    # Calculate Plant Net Power
    self.costing.plant_net_power = pyo.Expression(
        expr=pyo.units.convert(
            self.costing.plant_gross_power - self.costing.aux_load,
            to_units=pyo.units.MW)
        ) # MW

# -------------------------------------------------------------------------------------------
    # ***** TPC Calculation for NGCC and CCS *****
    
    import os
    import json
    from pyomo.common.fileutils import this_file_dir

    directory = this_file_dir()  # access the path to the folder containing this script

    with open(os.path.join(directory, "FEEDEquipmentDetailsFinal.json"), "r") as file:
        FEEDEquipmentDetailsFinal = json.load(file)
        
    pyo.units.load_definitions_from_strings(["USD_2021 = 500/708 * USD_CE500"])
    
    CE_index_year = "2018"
    CE_index_units = getattr(pyo.units, "MUSD_" + CE_index_year)

    # Accounts with Feedwater Flow to HP section of HRSG, as the
    # reference/scaling parameter - Exhibit 5-15

    FW_accounts = ['3.1', '3.3', '8.4']

    self.b1 = pyo.Block()

    self.b1.feedwater_flowrate = pyo.Var(initialize=feedwater_flowrate,units=pyunits.lb/pyunits.hr)  # lb/hr
    self.b1.feedwater_flowrate.fix(feedwater_flowrate)

    QGESSCostingData.get_PP_costing(
            self.b1,
            FW_accounts,
            self.b1.feedwater_flowrate,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b1_nocap = pyo.Block()

    self.b1_nocap.feedwater_flowrate = pyo.Var(initialize=1085751,units=pyunits.lb/pyunits.hr)  # lb/hr
    self.b1_nocap.feedwater_flowrate.fix(1085751)
    
    QGESSCostingData.get_PP_costing(
            self.b1_nocap,
            FW_accounts,
            self.b1_nocap.feedwater_flowrate,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with Raw water withdrawal as the reference/scaling parameter
    # Exhibit 5-14
    RW_withdraw_accounts = ['3.2', '3.4', '3.5', '9.5', '14.6']

    self.b2 = pyo.Block()

    self.b2.raw_water_withdrawal = pyo.Var(initialize=pyo.value(self.costing.raw_water_withdrawal),units=pyunits.gal/pyunits.min)  # gpm
    self.b2.raw_water_withdrawal.fix(pyo.value(self.costing.raw_water_withdrawal))
    
    QGESSCostingData.get_PP_costing(
            self.b2,
            RW_withdraw_accounts,
            self.b2.raw_water_withdrawal,
            6,
            CE_index_year=CE_index_year,
        )
    
    # No capture case
    self.b2_nocap = pyo.Block()

    self.b2_nocap.raw_water_withdrawal = pyo.Var(initialize=2902,units=pyunits.gal/pyunits.min)  # gpm
    self.b2_nocap.raw_water_withdrawal.fix(2902)
    QGESSCostingData.get_PP_costing(
            self.b2_nocap,
            RW_withdraw_accounts,
            self.b2_nocap.raw_water_withdrawal,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with fuel gas flowrate as the reference/scaling parameter
    # Exhibit 5-15 stream 2, Exhibit 5-8
    FuelG_accounts = ['3.6', '3.9', '6.1', '6.3', '6.4']

    self.b3 = pyo.Block()

    self.b3.fg_flowrate = pyo.Var(initialize=fuelgas_value,units=pyunits.lb/pyunits.hr)  # lb/hr
    self.b3.fg_flowrate.fix(fuelgas_value)
    
    QGESSCostingData.get_PP_costing(
            self.b3,
            FuelG_accounts,
            self.b3.fg_flowrate,
            6,
            CE_index_year=CE_index_year,
        )
    
    # No capture case
    self.b3_nocap = pyo.Block()

    self.b3_nocap.fg_flowrate = pyo.Var(initialize=205630,units=pyunits.lb/pyunits.hr)  # lb/hr
    self.b3_nocap.fg_flowrate.fix(205630)
    
    QGESSCostingData.get_PP_costing(
            self.b3_nocap,
            FuelG_accounts,
            self.b3_nocap.fg_flowrate,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with process water discharge as the reference/scaling parameter
    # Exhibit 5-14
    PW_discharge_accounts = ['3.7']

    self.b4 = pyo.Block()

    self.b4.process_water_discharge = pyo.Var(initialize=pyo.value(self.costing.process_water_discharge),units=pyunits.gal/pyunits.min)  # gpm
    self.b4.process_water_discharge.fix(pyo.value(self.costing.process_water_discharge))

    QGESSCostingData.get_PP_costing(
            self.b4,
            PW_discharge_accounts,
            self.b4.process_water_discharge,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b4_nocap = pyo.Block()

    self.b4_nocap.process_water_discharge = pyo.Var(initialize=657,units=pyunits.gal/pyunits.min)  # gpm
    self.b4_nocap.process_water_discharge.fix(657)
    
    QGESSCostingData.get_PP_costing(
            self.b4_nocap,
            PW_discharge_accounts,
            self.b4_nocap.process_water_discharge,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with flue gas flowrate as the reference/scaling parameter
    # Exhibit 5-15 stream 3, Exhibit 5-8
    FG_accounts = ['7.6']

    self.b5 = pyo.Block()

    self.b5.fg_flowrate = pyo.Var(initialize=(8658430/60)/0.025,units=pyunits.ft**3/pyunits.min)  # ft3/min
    self.b5.fg_flowrate.fix((8658430/60)/0.025)
    
    QGESSCostingData.get_PP_costing(
            self.b5,
            FG_accounts,
            self.b5.fg_flowrate,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b5_nocap = pyo.Block()

    self.b5_nocap.fg_flowrate = pyo.Var(initialize=(8658430/60)/0.025,units=pyunits.ft**3/pyunits.min)  # ft3/min
    self.b5_nocap.fg_flowrate.fix((8658430/60)/0.025)
    
    QGESSCostingData.get_PP_costing(
            self.b5_nocap,
            FG_accounts,
            self.b5_nocap.fg_flowrate,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with combustion turbine gross power as the reference/scaling
    # parameter
    # Exhibit 5-9
    CT_grosspower_accounts = ['6.5']

    self.b6 = pyo.Block()

    self.b6.ct_gross_power = pyo.Var(initialize=pyo.value(self.costing.CT_gross_power),units=pyunits.kW)  # kW
    self.b6.ct_gross_power.fix(pyo.value(self.costing.CT_gross_power))

    QGESSCostingData.get_PP_costing(
            self.b6,
            CT_grosspower_accounts,
            self.b6.ct_gross_power,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b6_nocap = pyo.Block()

    self.b6_nocap.ct_gross_power = pyo.Var(initialize=477*1000,units=pyunits.kW)  # kW
    self.b6_nocap.ct_gross_power.fix(477*1000)
    
    QGESSCostingData.get_PP_costing(
            self.b6_nocap,
            CT_grosspower_accounts,
            self.b6_nocap.ct_gross_power,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with HRSG duty as the reference/scaling parameter
    # Exhibit 5-8, streams 3 and 4
    HRSG_duty_accounts = ['7.1', '7.2']

    self.b7 = pyo.Block()

    self.b7.hrsg_duty = pyo.Var(initialize=pyo.value(self.costing.HRSG_duty),units=pyunits.MBtu/pyunits.hr)  # MMBtu/hr
    self.b7.hrsg_duty.fix(pyo.value(self.costing.HRSG_duty))
    
    QGESSCostingData.get_PP_costing(
            self.b7,
            HRSG_duty_accounts,
            self.b7.hrsg_duty,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b7_nocap = pyo.Block()

    self.b7_nocap.hrsg_duty = pyo.Var(initialize=-(-538.1+277.1)*8658430/(10**6),units=pyunits.MBtu/pyunits.hr)  # MMBtu/hr
    self.b7_nocap.hrsg_duty.fix(-(-538.1+277.1)*8658430/(10**6))

    QGESSCostingData.get_PP_costing(
            self.b7_nocap,
            HRSG_duty_accounts,
            self.b7_nocap.hrsg_duty,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with gas flow to stack as the reference/scaling parameter
    # Exhibit 5-8, stream 4
    Stack_flow_gas_accounts = ['7.3', '7.4', '7.5']

    self.b8 = pyo.Block()

    self.b8.stack_flow_gas = pyo.Var(initialize=stack_flow_gas,units=pyunits.ft**3/pyunits.min)  # acfm
    self.b8.stack_flow_gas.fix(stack_flow_gas)
    
    QGESSCostingData.get_PP_costing(
            self.b8,
            Stack_flow_gas_accounts,
            self.b8.stack_flow_gas,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b8_nocap = pyo.Block()

    self.b8_nocap.stack_flow_gas = pyo.Var(initialize=(8658430/60)/0.061,units=pyunits.ft**3/pyunits.min)  # ft3/min
    self.b8_nocap.stack_flow_gas.fix((8658430/60)/0.061)

    QGESSCostingData.get_PP_costing(
            self.b8_nocap,
            Stack_flow_gas_accounts,
            self.b8_nocap.stack_flow_gas,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with steam turbine gross power as the reference/scaling
    # parameter
    # Exhibit 5-9
    Steam_turbine_gross_power_accounts = ['8.1', '8.2', '8.5', '14.3']

    self.b9 = pyo.Block()

    self.b9.st_gross_power = pyo.Var(initialize=pyo.value(self.costing.ST_gross_power),units=pyunits.kW)  # kW
    self.b9.st_gross_power.fix(pyo.value(self.costing.ST_gross_power))

    QGESSCostingData.get_PP_costing(
            self.b9,
            Steam_turbine_gross_power_accounts,
            self.b9.st_gross_power,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b9_nocap = pyo.Block()

    self.b9_nocap.st_gross_power = pyo.Var(initialize=263*1000,units=pyunits.kW)  # kW
    self.b9_nocap.st_gross_power.fix(263*1000)
    
    QGESSCostingData.get_PP_costing(
            self.b9_nocap,
            Steam_turbine_gross_power_accounts,
            self.b9_nocap.st_gross_power,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with condenser duty as the reference/scaling parameter
    # Exhibit 5-9
    Condenser_duty_accounts = ['8.3']

    self.b10 = pyo.Block()

    self.b10.cond_duty = pyo.Var(initialize=pyo.value(self.costing.condenser_duty),units=pyunits.MBtu/pyunits.hr)  # MMBtu/hr
    self.b10.cond_duty.fix(pyo.value(self.costing.condenser_duty))

    QGESSCostingData.get_PP_costing(
            self.b10,
            Condenser_duty_accounts,
            self.b10.cond_duty,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b10_nocap = pyo.Block()

    self.b10_nocap.cond_duty = pyo.Var(initialize=1332,units=pyunits.MBtu/pyunits.hr)  # MMBtu/hr
    self.b10_nocap.cond_duty.fix(1332)
    
    QGESSCostingData.get_PP_costing(
            self.b10_nocap,
            Condenser_duty_accounts,
            self.b10_nocap.cond_duty,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with cooling tower duty as the reference/scaling parameter
    # Exhibit 5-16
    Cooling_tower_accounts = ['9.1']

    self.b11 = pyo.Block()

    # Obtain cooling tower duty in MMBtu/hr (includes condenser, Acid gas
    # removal, and other cooling loads)

    self.b11.cool_tower_duty = pyo.Var(initialize=pyo.value(self.costing.cooling_tower_duty),units=pyunits.MBtu/pyunits.hr)
    # MMBtu/hr
    self.b11.cool_tower_duty.fix(pyo.value(self.costing.cooling_tower_duty))

    QGESSCostingData.get_PP_costing(
            self.b11,
            Cooling_tower_accounts,
            self.b11.cool_tower_duty,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b11_nocap = pyo.Block()

    # Obtain cooling tower duty in MMBtu/hr (includes condenser, Acid gas
    # removal, and other cooling loads)

    self.b11_nocap.cool_tower_duty = pyo.Var(initialize=1357,units=pyunits.MBtu/pyunits.hr)
    # MMBtu/hr
    self.b11_nocap.cool_tower_duty.fix(1357)
    
    QGESSCostingData.get_PP_costing(
            self.b11_nocap,
            Cooling_tower_accounts,
            self.b11_nocap.cool_tower_duty,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with circulating water flowrate as the reference/scaling
    # parameter

    Circ_water_accounts = ['9.2', '9.3', '9.4', '9.6', '9.7', '14.5']

    self.b12 = pyo.Block()

    self.b12.circ_water_flow = pyo.Var(initialize=pyo.value(self.costing.cir_water_flowrate),units=pyunits.gal/pyunits.min)  # gpm
    self.b12.circ_water_flow.fix(pyo.value(self.costing.cir_water_flowrate))

    QGESSCostingData.get_PP_costing(
            self.b12,
            Circ_water_accounts,
            self.b12.circ_water_flow,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b12_nocap = pyo.Block()

    self.b12_nocap.circ_water_flow = pyo.Var(initialize=135967.5012,units=pyunits.gal/pyunits.min)  # gpm
    self.b12_nocap.circ_water_flow.fix(135967.5012)
    
    QGESSCostingData.get_PP_costing(
            self.b12_nocap,
            Circ_water_accounts,
            self.b12_nocap.circ_water_flow,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with total plant gross power as the reference/scaling parameter
    # Exhibit 5-9

    plant_gross_power_accounts = ['11.1', '11.7', '11.9', '13.1', '13.2',
                                  '13.3', '14.4', '14.7', '14.8', '14.9',
                                  '14.10']

    self.b13 = pyo.Block()

    self.b13.gross_power = pyo.Var(initialize=pyo.value(self.costing.plant_gross_power),units=pyunits.kW)  # kW
    self.b13.gross_power.fix(pyo.value(self.costing.plant_gross_power))

    QGESSCostingData.get_PP_costing(
            self.b13,
            plant_gross_power_accounts,
            self.b13.gross_power,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b13_nocap = pyo.Block()

    self.b13_nocap.gross_power = pyo.Var(initialize=740000,units=pyunits.kW)  # kW
    self.b13_nocap.gross_power.fix(740000)
    
    QGESSCostingData.get_PP_costing(
            self.b13_nocap,
            plant_gross_power_accounts,
            self.b13_nocap.gross_power,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with auxilliary load as the reference/scaling parameter
    # Exhibit 5-9

    auxilliary_load_accounts = ['11.2', '11.3', '11.4', '11.5', '11.6', '12.1',
                                '12.2', '12.3', '12.4', '12.5', '12.6', '12.7',
                                '12.8', '12.9']

    self.b14 = pyo.Block()

    self.b14.auxilliary_load = pyo.Var(initialize=pyo.value(self.costing.aux_load), units=pyo.units.kW)  # kW
    self.b14.auxilliary_load.fix(pyo.value(self.costing.aux_load))

    QGESSCostingData.get_PP_costing(
            self.b14,
            auxilliary_load_accounts,
            self.b14.auxilliary_load,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b14_nocap = pyo.Block()

    self.b14_nocap.auxilliary_load = pyo.Var(initialize=13552,units=pyunits.kW)  # kW
    self.b14_nocap.auxilliary_load.fix(13552)

    QGESSCostingData.get_PP_costing(
            self.b14_nocap,
            auxilliary_load_accounts,
            self.b14_nocap.auxilliary_load,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with STG,CTG output as the reference/scaling parameter
    # Case B31A Account 9 - Pg 501 rev 4 baseline report

    stg_ctg_accounts = ['11.8']

    self.b15 = pyo.Block()

    self.b15.stg_ctg_output = pyo.Var(initialize=pyo.value(self.costing.stg_ctg_op),units=pyunits.kW)  # kW
    self.b15.stg_ctg_output.fix(pyo.value(self.costing.stg_ctg_op))

    QGESSCostingData.get_PP_costing(
            self.b15,
            stg_ctg_accounts,
            self.b15.stg_ctg_output,
            6,
            CE_index_year=CE_index_year,
        )

    # No capture case
    self.b15_nocap = pyo.Block()

    self.b15_nocap.stg_ctg_output = pyo.Var(initialize=830000,units=pyunits.kW)  # kW
    self.b15_nocap.stg_ctg_output.fix(830000)

    QGESSCostingData.get_PP_costing(
            self.b15_nocap,
            stg_ctg_accounts,
            self.b15_nocap.stg_ctg_output,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )

    # Accounts with gas turbine power as the reference/scaling parameter
    # Exhibit 5-9

    gasturbine_accounts = ['14.1']

    self.b16 = pyo.Block()

    self.b16.gas_turbine_power = pyo.Var(initialize=pyo.value(self.costing.gt_power),units=pyunits.kW)  # kW
    self.b16.gas_turbine_power.fix(pyo.value(self.costing.gt_power))

    QGESSCostingData.get_PP_costing(
            self.b16,
            gasturbine_accounts,
            self.b16.gas_turbine_power,
            6,
            CE_index_year=CE_index_year,
        )

    self.b16_nocap = pyo.Block()

    self.b16_nocap.gas_turbine_power = pyo.Var(initialize=477*1000,units=pyunits.kW)  # kW
    self.b16_nocap.gas_turbine_power.fix(477*1000)

    QGESSCostingData.get_PP_costing(
            self.b16_nocap,
            gasturbine_accounts,
            self.b16_nocap.gas_turbine_power,
            6,
            ccs="A",
            CE_index_year=CE_index_year,
        )
    
    # ***Accounts with carbon capture system units***
    
    # Flue Gas Blower
    fg_blower_accounts = ["1.3"]
    self.b17 = pyo.Block()
    
    # Obtain blower hp
    self.costing.blower_hp = pyo.Expression(
        expr=pyo.units.convert(FGBlowerLoad, to_units=pyo.units.hp)
        )  # hp

    self.b17.fg_blower_hp = pyo.Var(initialize=pyo.value(self.costing.blower_hp), units=pyunits.hp)
    self.b17.fg_blower_hp.fix(pyo.value(self.costing.blower_hp))
    
    QGESSCostingData.get_PP_costing(
        self.b17,
        fg_blower_accounts,
        self.b17.fg_blower_hp,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # DCC Column
    DCC_column_accounts = ["1.6"]
    self.b18 = pyo.Block()
    
    # Obtain flue gas flowrate per train (acfm)
    self.costing.FG_flowrate_per_train = pyo.Expression(
        expr=fluegas_value/NoTrain) # acfm

    self.b18.fg_flowrate = pyo.Var(initialize=pyo.value(self.costing.FG_flowrate_per_train), units=pyunits.ft**3/pyunits.min)
    self.b18.fg_flowrate.fix()
    
    QGESSCostingData.get_PP_costing(
        self.b18,
        DCC_column_accounts,
        self.b18.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # DCC Column Packing
    DCC_column_packing_accounts = ["1.9"]
    self.b19 = pyo.Block()

    self.b19.fg_flowrate = pyo.Var(initialize=pyo.value(self.costing.FG_flowrate_per_train), units=pyunits.ft**3/pyunits.min)
    self.b19.fg_flowrate.fix(pyo.value(self.costing.FG_flowrate_per_train))
    
    QGESSCostingData.get_PP_costing(
        self.b19,
        DCC_column_packing_accounts,
        self.b19.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # DCC Pump
    DCC_pump_accounts = ["1.13"]
    self.b20 = pyo.Block()

    self.b20.fg_flowrate = pyo.Var(initialize=pyo.value(self.costing.FG_flowrate_per_train), units=pyunits.ft**3/pyunits.min)
    self.b20.fg_flowrate.fix(pyo.value(self.costing.FG_flowrate_per_train))
    QGESSCostingData.get_PP_costing(
        self.b20,
        DCC_pump_accounts,
        self.b20.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # DCC Cooler    
    DCC_cooler_accounts = ["1.16"]
    self.b21 = pyo.Block()

    self.b21.fg_flowrate = pyo.Var(initialize=pyo.value(self.costing.FG_flowrate_per_train), units=pyunits.ft**3/pyunits.min)
    self.b21.fg_flowrate.fix(pyo.value(self.costing.FG_flowrate_per_train))
    QGESSCostingData.get_PP_costing(
        self.b21,
        DCC_cooler_accounts,
        self.b21.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Absorber
    # Column volume calculation
    # Obtain absorber volume in m3
    self.costing.absorber_volume = pyo.Expression(
        expr=pyo.units.convert(
            pi*(1.25*(Abs_ht+wash_section_packing_ht) + (Abs_ht+wash_section_packing_ht)*0.75/10)*Abs_d**2/4,
            to_units=pyo.units.m**3
            )
        ) # m3
    
    absorber_accounts = ["2.3"]
    self.b22 = pyo.Block() 

    self.b22.absorber_volume = pyo.Var(initialize=pyo.value(self.costing.absorber_volume), units=pyunits.m**3)
    self.b22.absorber_volume.fix(pyo.value(self.costing.absorber_volume))
    QGESSCostingData.get_PP_costing(
        self.b22,
        absorber_accounts,
        self.b22.absorber_volume,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Absorber packing
    # Packing volume calculation
    # Obtain absorber packing volume in m3
    self.costing.absorber_packing_volume = pyo.Expression(
        expr=pyo.units.convert(
            pi*(Abs_ht+wash_section_packing_ht)*(Abs_d**2)/4,
            to_units=pyo.units.m**3
            )
        ) # m3
    
    absorber_packing_accounts = ["3.3"]
    self.b23 = pyo.Block()  

    self.b23.absorber_packing_volume = pyo.Var(
        initialize=pyo.value(self.costing.absorber_packing_volume), units=pyunits.m**3
    )
    self.b23.absorber_packing_volume.fix(pyo.value(self.costing.absorber_packing_volume))
    QGESSCostingData.get_PP_costing(
        self.b23,
        absorber_packing_accounts,
        self.b23.absorber_packing_volume,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Absorber intercooler 1
    # Cooling water conditions
    self.costing.CWTin = pyo.Param(initialize=16 + 273.15, mutable=False, units=pyo.units.K) # 16 deg C, converted to K
    self.costing.CWTout = pyo.Param(initialize=27 + 273.15, mutable=False, units=pyo.units.K) # 27 deg C, converted to K
    # IC1 heat transfer surface area calculation
    # self.costing.U_ICArea = pyo.Param(initialize=850, mutable=False, units=pyo.units.W/pyo.units.m**2/pyo.units.K)  # W/m2.K
    # self.costing.LMTD_IC1 = pyo.Expression(
    #     expr=((ICOOL1Tin - self.costing.CWTout) - (ICOOL1Tout - self.costing.CWTin))/\
    #     pyo.log((ICOOL1Tin - self.costing.CWTout)/(ICOOL1Tout - self.costing.CWTin))
    #     )
    # Obtain absorber IC1 heat transfer surface area in m2
    # self.costing.absorber_ic1_htsa = pyo.Expression(
    #     expr=-QICOOL1/(self.costing.U_ICArea*self.costing.LMTD_IC1)
    #     )  # m2
    
    # absorber_ic1_accounts = ["4.3.1"]
    # self.b24 = pyo.Block()

    # self.b24.absorber_ic1_htsa = pyo.Var(
    #     initialize=pyo.value(self.costing.absorber_ic1_htsa), units=pyunits.m**2
    # )
    # self.b24.absorber_ic1_htsa.fix(pyo.value(self.costing.absorber_ic1_htsa))
    # QGESSCostingData.get_PP_costing(
    #     self.b24,
    #     absorber_ic1_accounts,
    #     self.b24.absorber_ic1_htsa,
    #     6,
    #     CE_index_year=CE_index_year,
    #     additional_costing_params=FEEDEquipmentDetailsFinal,
    #     use_additional_costing_params=True,
    # )
    
    # Absorber intercooler 2
    # IC2 heat transfer surface area calculation
    # self.costing.LMTD_IC2 = pyo.Expression(
    #     expr=((ICOOL2Tin - self.costing.CWTout) - (ICOOL2Tout - self.costing.CWTin))/\
    #     pyo.log((ICOOL2Tin - self.costing.CWTout)/(ICOOL2Tout - self.costing.CWTin))
    #     )
    # Obtain absorber IC2 heat transfer surface area in m2
    # self.costing.absorber_ic2_htsa = pyo.Expression(
    #     expr=-QICOOL2/(self.costing.U_ICArea*self.costing.LMTD_IC2)
    #     )  # m2
    
    # absorber_ic2_accounts = ["4.3.2"]
    # self.b25 = pyo.Block()

    # self.b25.absorber_ic2_htsa = pyo.Var(
    #     initialize=pyo.value(self.costing.absorber_ic2_htsa), units=pyunits.m**2
    # )
    # self.b25.absorber_ic2_htsa.fix(pyo.value(self.costing.absorber_ic2_htsa))
    # QGESSCostingData.get_PP_costing(
    #     self.b25,
    #     absorber_ic2_accounts,
    #     self.b25.absorber_ic2_htsa,
    #     6,
    #     CE_index_year=CE_index_year,
    #     additional_costing_params=FEEDEquipmentDetailsFinal,
    #     use_additional_costing_params=True,
    # )

    # Absorber intercooling pump 1
    # absorber_ic1pump_accounts = ["5.3.1"]
    # self.b26 = pyo.Block()
    
    # Obtain absorber IC1 pump design flow in gpm
    # Consider 10 % design margin
    # self.costing.absorber_ic1pump_design_flow = 1.1*pyo.units.convert(
    #     pumparound_1_flowrate, to_units=pyo.units.gal/pyo.units.min) # gpm

    # self.b26.absorber_ic1pump_design_flow = pyo.Var(
    #     initialize=pyo.value(self.costing.absorber_ic1pump_design_flow), units=pyunits.gal/pyunits.min
    # )
    # self.b26.absorber_ic1pump_design_flow.fix(pyo.value(self.costing.absorber_ic1pump_design_flow))
    # QGESSCostingData.get_PP_costing(
    #     self.b26,
    #     absorber_ic1pump_accounts,
    #     self.b26.absorber_ic1pump_design_flow,
    #     6,
    #     CE_index_year=CE_index_year,
    #     additional_costing_params=FEEDEquipmentDetailsFinal,
    #     use_additional_costing_params=True,
    # )
    
    # Absorber intercooling pump 2
    # absorber_ic2pump_accounts = ["5.3.2"]
    # self.b27 = pyo.Block()
    
    # Obtain absorber IC2 pump design flow in gpm
    # Consider 10 % design margin
    # self.costing.absorber_ic2pump_design_flow = 1.1*pyo.units.convert(
    #     pumparound_2_flowrate, to_units=pyo.units.gal/pyo.units.min) # gpm

    # self.b27.absorber_ic2pump_design_flow = pyo.Var(
    #     initialize=pyo.value(self.costing.absorber_ic2pump_design_flow), units=pyunits.gal/pyunits.min
    # )
    # self.b27.absorber_ic2pump_design_flow.fix(pyo.value(self.costing.absorber_ic2pump_design_flow))
    # QGESSCostingData.get_PP_costing(
    #     self.b27,
    #     absorber_ic2pump_accounts,
    #     self.b27.absorber_ic2pump_design_flow,
    #     6,
    #     CE_index_year=CE_index_year,
    #     additional_costing_params=FEEDEquipmentDetailsFinal,
    #     use_additional_costing_params=True,
    # )
    
    # WW Recirculation Pump
    ww_pump_accounts = ["6.3"]
    self.b28 = pyo.Block()
    
    # Assume L/G ratio in the wash column is equivalent to 1
    # Assume 1 kg water is equivalent to 1 litre
    # Obtain wash water rated flow in gpm, considering 10 % margin
    self.costing.ww_rated_flow = 1.1*pyo.units.convert(
        absorber_off_gas_flowrate*pyo.units.L/pyo.units.kg,
        to_units=pyo.units.gal/pyo.units.min)  # gpm

    self.b28.ww_rated_flow = pyo.Var(
        initialize=pyo.value(self.costing.ww_rated_flow), units=pyunits.gal/pyunits.min
    )
    self.b28.ww_rated_flow.fix(pyo.value(self.costing.ww_rated_flow))
    QGESSCostingData.get_PP_costing(
        self.b28,
        ww_pump_accounts,
        self.b28.ww_rated_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # WW Cooler
    ww_cooler_accounts = ["7.3"]
    self.b29 = pyo.Block()
    
    # Obtain heat transfer surface area in m2

    self.b29.ww_cooler_htsa = pyo.Var(
        initialize=pyo.value(self.costing.wwcooler_htsa), units=pyunits.m**2
    )
    self.b29.ww_cooler_htsa.fix(pyo.value(self.costing.wwcooler_htsa))
    QGESSCostingData.get_PP_costing(
        self.b29,
        ww_cooler_accounts,
        self.b29.ww_cooler_htsa,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Rich solvent pump
    rich_solvent_pump_account = ["8.3"]
    self.b30 = pyo.Block()
    
    # Consider 10 % margin
    # Obtain Rich Solvent Rated Flowrate
    self.costing.rs_rated_flow = 1.1*pyo.units.convert(
        rich_solvent_flowrate, to_units=pyo.units.gal/pyo.units.min) # gpm

    self.b30.rs_rated_flow = pyo.Var(
        initialize=pyo.value(self.costing.rs_rated_flow), units=pyunits.gal/pyunits.min
    )
    self.b30.rs_rated_flow.fix(pyo.value(self.costing.rs_rated_flow))
    QGESSCostingData.get_PP_costing(
        self.b30,
        rich_solvent_pump_account,
        self.b30.rs_rated_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Lean rich heat exchanger
    lean_rich_hex_accounts = ["9.3"]
    self.b31 = pyo.Block()
    
    # Obtain lean rich heat exchanger area in m2

    self.b31.lean_rich_hex_area = pyo.Var(
        initialize=LRHEXArea, units=pyunits.m**2
    )
    self.b31.lean_rich_hex_area.fix(LRHEXArea)
    QGESSCostingData.get_PP_costing(
        self.b31,
        lean_rich_hex_accounts,
        self.b31.lean_rich_hex_area,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Stripper
    stripper_accounts = ["10.3"]
    self.b32 = pyo.Block()
    
    # Obtain stripper volume in m3
    self.costing.stripper_volume = pyo.Expression(
        expr=pi*(1.25*Strip_ht + Strip_ht*0.75/10)*Strip_d**2/4) # m3

    self.b32.stripper_volume = pyo.Var(initialize=pyo.value(self.costing.stripper_volume), units=pyunits.m**3)
    self.b32.stripper_volume.fix(pyo.value(self.costing.stripper_volume))
    QGESSCostingData.get_PP_costing(
        self.b32,
        stripper_accounts,
        self.b32.stripper_volume,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Stripper packing
    stripper_packing_accounts = ["11.3"]
    self.b33 = pyo.Block()
    
    # Obtain stripper packing volume in m3
    self.costing.stripper_packing_volume = pyo.Expression(
        expr=pi*Strip_ht*(Strip_d**2)/4)  # m3

    self.b33.stripper_packing_volume = pyo.Var(
        initialize=pyo.value(self.costing.stripper_packing_volume), units=pyunits.m**3
    )
    self.b33.stripper_packing_volume.fix(pyo.value(self.costing.stripper_packing_volume))
    QGESSCostingData.get_PP_costing(
        self.b33,
        stripper_packing_accounts,
        self.b33.stripper_packing_volume,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Stripper condenser
    # Area calculation
    self.costing.U_condenser = pyo.Param(initialize=375, mutable=False, units=pyo.units.W/pyo.units.m**2/pyo.units.K)  # W/m2.K
    self.costing.LMTD_condenser = pyo.Expression(
        expr=((StripcondTin - self.costing.CWTout) - (StripcondTout - self.costing.CWTin))/\
        pyo.log((StripcondTin - self.costing.CWTout)/(StripcondTout - self.costing.CWTin))
        )
    self.costing.CondArea = pyo.Expression(
        expr=-QStripcond/(self.costing.U_condenser*self.costing.LMTD_condenser)
        )
    
    stripper_condenser_accounts = ["12.3"]
    self.b34 = pyo.Block()
    
    # Obtain stripper condenser area in m2

    self.b34.stripper_condenser_area = pyo.Var(
        initialize=pyo.value(self.costing.CondArea), units=pyunits.m**2
    )
    self.b34.stripper_condenser_area.fix(pyo.value(self.costing.CondArea))
    QGESSCostingData.get_PP_costing(
        self.b34,
        stripper_condenser_accounts,
        self.b34.stripper_condenser_area,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Stripper reboiler
    # Area calculation
    self.costing.U_reboiler = pyo.Param(initialize=1050, mutable=False, units=pyo.units.W/pyo.units.K)  # W/m2.K
    self.costing.TSTEAM = pyo.Param(initialize=308 + 273.15, mutable=False, units=pyo.units.K) # 308 deg C, converted to K
    self.costing.TCONDENSATE = pyo.Param(initialize=151 + 273.15, mutable=False, units=pyo.units.K) # 151 deg C, converted to K
    self.costing.dT1 = pyo.Expression(
        expr=self.costing.TSTEAM-Rebtempout
        )
    self.costing.dT2 = pyo.Expression(
        expr=self.costing.TCONDENSATE-Rebtempout
        )
    self.costing.LMTDZONE1 = pyo.Expression(
        expr=(self.costing.dT1-self.costing.dT2)/(pyo.log(self.costing.dT1/self.costing.dT2))
        )
    self.costing.LMTDZONE2 = pyo.Expression(
        expr=self.costing.TCONDENSATE-Rebtempout
        )
    self.costing.QZONE1 = pyo.Expression(
        expr=pyo.units.convert(Qreb, to_units=pyo.units.W)/(1+6.9)
        )
    self.costing.AREAZONE1 = pyo.Expression(
        expr=self.costing.QZONE1/(self.costing.U_reboiler*self.costing.LMTDZONE1)
        )
    self.costing.QZONE2 = pyo.Expression(
        expr=pyo.units.convert(Qreb, to_units=pyo.units.W)-self.costing.QZONE1
        )
    self.costing.AREAZONE2 = pyo.Expression(
        expr=self.costing.QZONE2/(self.costing.U_reboiler*self.costing.LMTDZONE2)
        )
    self.costing.RbrArea = pyo.Expression(
        expr=self.costing.AREAZONE1+self.costing.AREAZONE2
        )
    
    stripper_reboiler_accounts = ["13.3"]
    self.b35 = pyo.Block()
    
    # Obtain stripper reboiler area in m2

    self.b35.stripper_reboiler_area = pyo.Var(
        initialize=pyo.value(self.costing.RbrArea), units=pyunits.m**2
    )
    self.b35.stripper_reboiler_area.fix(pyo.value(self.costing.RbrArea))
    QGESSCostingData.get_PP_costing(
        self.b35,
        stripper_reboiler_accounts,
        self.b35.stripper_reboiler_area,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Stripper reflux drum
    stripper_reflux_drum_accounts = ["14.3"]
    self.b36 = pyo.Block()
    
    self.b36.co2cap = pyo.Var(
        initialize=pyo.value(
            pyo.units.convert(CO2CaptureRate, to_units=pyo.units.kg/pyo.units.hr)/NoTrain),
        units=pyo.units.kg/pyo.units.hr
        )
    
    self.b36.co2cap.fix(pyo.value(
        pyo.units.convert(CO2CaptureRate, to_units=pyo.units.kg/pyo.units.hr)/NoTrain)
        )
    QGESSCostingData.get_PP_costing(
        self.b36,
        stripper_reflux_drum_accounts,
        self.b36.co2cap,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Lean solvent pump
    lean_solvent_pump_account = ["15.3"]
    self.b37 = pyo.Block()
    # Obtain Lean Solvent Rated Flowrate
    self.costing.ls_rated_flow = 1.1*pyo.units.convert(
        lean_solvent_flowrate, to_units=pyo.units.gal/pyo.units.min)  # gpm

    self.b37.ls_rated_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_rated_flow), units=pyunits.gal/pyunits.min
    )
    self.b37.ls_rated_flow.fix(pyo.value(self.costing.ls_rated_flow))
    QGESSCostingData.get_PP_costing(
        self.b37,
        lean_solvent_pump_account,
        self.b37.ls_rated_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Lean solvent cooler
    # Area calculation
    # self.costing.U_LSCoolerArea = pyo.Param(initialize=850, mutable=False, units=pyo.units.W/pyo.units.m**2/pyo.units.K)  # W/m2.K
    # self.costing.LMTD_LSCooler = pyo.Expression(
    #     expr=((LSCoolerTin - self.costing.CWTout) - (LSCoolerTout - self.costing.CWTin))/\
    #     pyo.log((LSCoolerTin - self.costing.CWTout)/(LSCoolerTout - self.costing.CWTin))
    #     )
    # self.costing.LSCoolerArea = pyo.Expression(
    #     expr=-QLSCooler/(self.costing.U_LSCoolerArea*self.costing.LMTD_LSCooler)
    #     )
        
    # lean_solvent_cooler_accounts = ["16.3"]
    # self.b38 = pyo.Block()
    
    # Obtain lean solvent cooler area in m2

    # self.b38.lean_solvent_cooler_area = pyo.Var(
    #     initialize=pyo.value(self.costing.LSCoolerArea), units=pyunits.m**2
    # )
    # self.b38.lean_solvent_cooler_area.fix(pyo.value(self.costing.LSCoolerArea))
    # QGESSCostingData.get_PP_costing(
    #     self.b38,
    #     lean_solvent_cooler_accounts,
    #     self.b38.lean_solvent_cooler_area,
    #     6,
    #     CE_index_year=CE_index_year,
    #     additional_costing_params=FEEDEquipmentDetailsFinal,
    #     use_additional_costing_params=True,
    # )

    # Solvent stripper reclaimer
    solvent_stripper_reclaimer_account = ["17.3"]
    self.b39 = pyo.Block()
    
    # Obtain the MEA makeup flowrate in kg/hr
    self.costing.solvent_makeup_flowrate = pyo.units.convert(
        solvent_makeup*CO2CaptureRate, to_units=pyo.units.kg/pyo.units.hr) # kg/hr

    self.b39.mea_makeup_flow = pyo.Var(
        initialize=pyo.value(self.costing.solvent_makeup_flowrate), units=pyunits.kg / pyunits.hr
    )
    self.b39.mea_makeup_flow.fix(pyo.value(self.costing.solvent_makeup_flowrate))
    QGESSCostingData.get_PP_costing(
        self.b39,
        solvent_stripper_reclaimer_account,
        self.b39.mea_makeup_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

    # Solvent filtration
    solvent_filtration_account = ["18.3"]
    self.b40 = pyo.Block()
    
    self.costing.ls_flowrate = pyo.units.convert(
        lean_solvent_flowrate, to_units=pyo.units.gal/pyo.units.min) # gpm

    self.b40.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate), units=pyunits.gal / pyunits.min
    )
    self.b40.ls_flow.fix(pyo.value(self.costing.ls_flowrate))
    QGESSCostingData.get_PP_costing(
        self.b40,
        solvent_filtration_account,
        self.b40.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Solvent Storage Tank
    solvent_storage_tank_account = ["19.3"]
    self.b41 = pyo.Block()

    self.b41.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate*NoTrain), units=pyunits.gal / pyunits.min
    )
    self.b41.ls_flow.fix(pyo.value(self.costing.ls_flowrate*NoTrain))
    QGESSCostingData.get_PP_costing(
        self.b41,
        solvent_storage_tank_account,
        self.b41.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Stripper Reflux Pump
    Stripper_Reflux_Pump_account = ["50.3"]
    self.b42 = pyo.Block()

    self.b42.co2cap = pyo.Var(
        initialize=pyo.value(
            pyo.units.convert(
                CO2CaptureRate/NoTrain,
                to_units=pyo.units.kg/pyo.units.hr)
            ),
            units=pyunits.kg / pyunits.hr)

    self.b42.co2cap.fix(pyo.value(
        pyo.units.convert(
            CO2CaptureRate/NoTrain,
            to_units=pyo.units.kg/pyo.units.hr)
        )
        )
    QGESSCostingData.get_PP_costing(
        self.b42,
        Stripper_Reflux_Pump_account,
        self.b42.co2cap,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Lean solvent pump 2
    LS_Pump2_account = ["65.1"]
    self.b43 = pyo.Block()

    self.b43.co2cap = pyo.Var(
        initialize=pyo.value(
            pyo.units.convert(
                CO2CaptureRate/NoTrain,
                to_units=pyo.units.kg/pyo.units.hr)
            ),
            units=pyunits.kg / pyunits.hr)

    self.b43.co2cap.fix(pyo.value(
        pyo.units.convert(
            CO2CaptureRate/NoTrain,
            to_units=pyo.units.kg/pyo.units.hr)
        )
        )
    QGESSCostingData.get_PP_costing(
        self.b43,
        LS_Pump2_account,
        self.b43.co2cap,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Solvent Sump
    Solvent_Sump_account = ["51.1"]
    self.b44 = pyo.Block()

    self.b44.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate), units=pyunits.gal / pyunits.min
    )
    
    self.b44.ls_flow.fix(pyo.value(self.costing.ls_flowrate))
    QGESSCostingData.get_PP_costing(
        self.b44,
        Solvent_Sump_account,
        self.b44.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Solvent Sump Pump
    Solvent_Sump_pump_account = ["52.1"]
    self.b45 = pyo.Block()

    self.b45.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate), units=pyunits.gal / pyunits.min
    )
    
    self.b45.ls_flow.fix(pyo.value(self.costing.ls_flowrate))
    QGESSCostingData.get_PP_costing(
        self.b45,
        Solvent_Sump_pump_account,
        self.b45.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Solvent Sump Filter
    Solvent_Sump_Filter_account = ["53.1"]
    self.b46 = pyo.Block()

    self.b46.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate), units=pyunits.gal / pyunits.min
    )
    
    self.b46.ls_flow.fix(pyo.value(self.costing.ls_flowrate))
    QGESSCostingData.get_PP_costing(
        self.b46,
        Solvent_Sump_Filter_account,
        self.b46.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Solvent Sump Pit Pump
    Solvent_Sump_Pit_Pump_account = ["54.1"]
    self.b47 = pyo.Block()

    self.b47.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate), units=pyunits.gal / pyunits.min
    )
    
    self.b47.ls_flow.fix(pyo.value(self.costing.ls_flowrate))
    QGESSCostingData.get_PP_costing(
        self.b47,
        Solvent_Sump_Pit_Pump_account,
        self.b47.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Pre Scrubber Pump
    Pre_Scrubber_Pump_account = ["55.1"]
    self.b48 = pyo.Block()

    self.b48.fg_flowrate = pyo.Var(initialize=pyo.value(self.costing.FG_flowrate_per_train), units=pyunits.ft**3/pyunits.min)
    
    self.b48.fg_flowrate.fix(pyo.value(self.costing.FG_flowrate_per_train))
    QGESSCostingData.get_PP_costing(
        self.b48,
        Pre_Scrubber_Pump_account,
        self.b48.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Reboiler condensate pot
    Reboiler_condensate_pot_account = ["56.1"]
    self.b49 = pyo.Block()

    self.b49.rebduty = pyo.Var(initialize=Qreb, units=pyunits.MW)
    
    self.b49.rebduty.fix(Qreb)
    QGESSCostingData.get_PP_costing(
        self.b49,
        Reboiler_condensate_pot_account,
        self.b49.rebduty,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # DCC Water Filter
    DCC_Water_Filter_account = ["57.1"]
    self.b50 = pyo.Block()

    self.b50.fg_flowrate = pyo.Var(initialize=pyo.value(self.costing.FG_flowrate_per_train), units=pyunits.ft**3/pyunits.min)
    
    self.b50.fg_flowrate.fix(pyo.value(self.costing.FG_flowrate_per_train))
    QGESSCostingData.get_PP_costing(
        self.b50,
        DCC_Water_Filter_account,
        self.b50.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Corrosion Inhib Package
    Corrosion_Inhib_Package_account = ["58.1"]
    self.b51 = pyo.Block()

    self.b51.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate), units=pyunits.gal / pyunits.min
    )
    
    self.b51.ls_flow.fix(pyo.value(self.costing.ls_flowrate))
    QGESSCostingData.get_PP_costing(
        self.b51,
        Corrosion_Inhib_Package_account,
        self.b51.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Antifoam Feed Package
    Antifoam_Feed_Package_account = ["59.1"]
    self.b52 = pyo.Block()

    self.b52.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate), units=pyunits.gal / pyunits.min
    )
    
    self.b52.ls_flow.fix(pyo.value(self.costing.ls_flowrate))
    QGESSCostingData.get_PP_costing(
        self.b52,
        Antifoam_Feed_Package_account,
        self.b52.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # NaOH Makeup Pump
    NaOH_Makeup_Pump_account = ["60.1"]
    self.b53 = pyo.Block()

    self.b53.fg_flowrate = pyo.Var(initialize=fluegas_value, units=pyunits.ft**3/pyunits.min)
    
    self.b53.fg_flowrate.fix(fluegas_value)
    QGESSCostingData.get_PP_costing(
        self.b53,
        NaOH_Makeup_Pump_account,
        self.b53.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Solvent Makeup Pump
    Solvent_Makeup_Pump_account = ["61.1"]
    self.b54 = pyo.Block()

    self.b54.ls_flow = pyo.Var(
        initialize=pyo.value(fluegas_value*NoTrain), units=pyunits.gal / pyunits.min
    )
    
    self.b54.ls_flow.fix(pyo.value(fluegas_value*NoTrain))
    QGESSCostingData.get_PP_costing(
        self.b54,
        Solvent_Makeup_Pump_account,
        self.b54.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # NaOH Storage Tank
    NaOH_Storage_Tank_account = ["62.1"]
    self.b55 = pyo.Block()

    self.b55.fg_flowrate = pyo.Var(initialize=fluegas_value, units=pyunits.ft**3/pyunits.min)
    
    self.b55.fg_flowrate.fix(fluegas_value)
    QGESSCostingData.get_PP_costing(
        self.b55,
        NaOH_Storage_Tank_account,
        self.b55.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Additional Line Items
    
    # Station Service Equipment
    Station_Service_Equipment_account = ["22.1"]
    self.b56 = pyo.Block()

    self.b56.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b56.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b56,
        Station_Service_Equipment_account,
        self.b56.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Switchgear & Motor Control
    Switchgear_Motor_Control_account = ["23.1"]
    self.b57 = pyo.Block()

    self.b57.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b57.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b57,
        Switchgear_Motor_Control_account,
        self.b57.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Conduit & Cable Tray
    Conduit_Cable_Tray_account = ["24.1"]
    self.b58 = pyo.Block()

    self.b58.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b58.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b58,
        Conduit_Cable_Tray_account,
        self.b58.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Wire & Cable
    Wire_Cable_account = ["25.1"]
    self.b59 = pyo.Block()

    self.b59.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b59.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b59,
        Wire_Cable_account,
        self.b59.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Main Power Transformers
    Main_Power_Transformers_account = ["26.1"]
    self.b60 = pyo.Block()

    self.b60.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b60.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b60,
        Main_Power_Transformers_account,
        self.b60.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Electrical Foundations
    Electrical_Foundations_account = ["27.1"]
    self.b61 = pyo.Block()

    self.b61.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b61.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b61,
        Electrical_Foundations_account,
        self.b61.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Control Boards, Panels & Racks
    Control_Boards_Panels_Racks_account = ["28.1"]
    self.b62 = pyo.Block()

    self.b62.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b62.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b62,
        Control_Boards_Panels_Racks_account,
        self.b62.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Distributed Control System Equipment
    Distributed_Control_System_Equipment_account = ["29.1"]
    self.b63 = pyo.Block()

    self.b63.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b63.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b63,
        Distributed_Control_System_Equipment_account,
        self.b63.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Instrument Wiring & Tubing
    Instrument_Wiring_Tubing_account = ["30.1"]
    self.b64 = pyo.Block()

    self.b64.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b64.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b64,
        Instrument_Wiring_Tubing_account,
        self.b64.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Other Instrumentation & Controls Equipment
    Other_Instrumentation_Controls_Equipment_account = ["31.1"]
    self.b65 = pyo.Block()

    self.b65.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b65.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b65,
        Other_Instrumentation_Controls_Equipment_account,
        self.b65.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Other Buildings & Structures
    Other_Buildings_Structures_account = ["33.1"]
    self.b66 = pyo.Block()

    self.b66.auxload = pyo.Var(initialize=pyo.value(self.costing.AuxLoad_CCS), units=pyunits.kW)
    
    self.b66.auxload.fix(pyo.value(self.costing.AuxLoad_CCS))
    QGESSCostingData.get_PP_costing(
        self.b66,
        Other_Buildings_Structures_account,
        self.b66.auxload,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Foundations
    Foundations_accounts = ["34.1"]
    self.b67 = pyo.Block()

    self.b67.fg_flowrate = pyo.Var(initialize=fluegas_value, units=pyunits.ft**3/pyunits.min)
    self.b67.fg_flowrate.fix(fluegas_value)
    
    QGESSCostingData.get_PP_costing(
        self.b67,
        Foundations_accounts,
        self.b67.fg_flowrate,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )
    
    # Interconnecting Piping
    Interconnecting_Piping_account = ["35.1"]
    self.b68 = pyo.Block()

    self.b68.ls_flow = pyo.Var(
        initialize=pyo.value(self.costing.ls_flowrate*NoTrain), units=pyunits.gal / pyunits.min
    )
    self.b68.ls_flow.fix(pyo.value(self.costing.ls_flowrate*NoTrain))
    QGESSCostingData.get_PP_costing(
        self.b68,
        Interconnecting_Piping_account,
        self.b68.ls_flow,
        6,
        CE_index_year=CE_index_year,
        additional_costing_params=FEEDEquipmentDetailsFinal,
        use_additional_costing_params=True,
    )

# -------------------------------------------------------------------------------------------
    # Build cost constraints    
    self.costing.build_process_costs(
        net_power=None,
        fixed_OM=False,
        variable_OM=False,
        resources=None,
        rates=None,
        prices=None,
        fuel=None,
        CE_index_year=CE_index_year,
    )

    print("\nCosting blocks built, initializing...\n")
    # Initialize costing
    QGESSCostingData.costing_initialization(self.costing)
    
    # Get default solver for testing
    solver = get_solver()

    print("\nCosting equations initialized, solving...\n")
    # Solve the model
    results = solver.solve(self.costing, tee=True)
# -------------------------------------------------------------------------------------------    
    # Obtain the total plant costs for each account - with capture

    self.costing.NGCC_TPC = pyo.Expression(
        expr=
        (pyo.units.convert(
            sum(self.b1.total_plant_cost[ac] for ac in FW_accounts) +
            sum(self.b2.total_plant_cost[ac] for ac in RW_withdraw_accounts) +
            sum(self.b3.total_plant_cost[ac] for ac in FuelG_accounts) +
            sum(self.b4.total_plant_cost[ac] for ac in PW_discharge_accounts) +
            sum(self.b5.total_plant_cost[ac] for ac in FG_accounts) +
            sum(self.b6.total_plant_cost[ac] for ac in CT_grosspower_accounts) +
            sum(self.b7.total_plant_cost[ac] for ac in HRSG_duty_accounts) +
            sum(self.b8.total_plant_cost[ac] for ac in Stack_flow_gas_accounts) +
            sum(self.b9.total_plant_cost[ac] for ac in Steam_turbine_gross_power_accounts) +
            sum(self.b10.total_plant_cost[ac] for ac in Condenser_duty_accounts) +
            sum(self.b11.total_plant_cost[ac] for ac in Cooling_tower_accounts) +
            sum(self.b12.total_plant_cost[ac] for ac in Circ_water_accounts) +
            sum(self.b13.total_plant_cost[ac] for ac in plant_gross_power_accounts) +
            sum(self.b14.total_plant_cost[ac] for ac in auxilliary_load_accounts) +
            sum(self.b15.total_plant_cost[ac] for ac in stg_ctg_accounts) +
            sum(self.b16.total_plant_cost[ac] for ac in gasturbine_accounts),
            to_units=getattr(pyunits, "USD_" + CE_index_year)  # $
            )
            )
        )
    
    # CCS TPC components
    self.costing.FG_blower_cost = pyo.Expression(
        expr=sum(self.b17.total_plant_cost[ac] for ac in fg_blower_accounts)
        )

    self.costing.DCC_column_cost = pyo.Expression(
        expr=sum(self.b18.total_plant_cost[ac] for ac in DCC_column_accounts)
        )
    
    self.costing.DCC_packing_cost = pyo.Expression(
        expr=sum(self.b19.total_plant_cost[ac] for ac in DCC_column_packing_accounts)
        )
    
    self.costing.DCC_pump_cost = pyo.Expression(
        expr=sum(self.b20.total_plant_cost[ac] for ac in DCC_pump_accounts)
        )
    
    self.costing.DCC_cooler_cost = pyo.Expression(
        expr=sum(self.b21.total_plant_cost[ac] for ac in DCC_cooler_accounts)
        )
    
    self.costing.absorber_column_cost = pyo.Expression(
        expr=sum(self.b22.total_plant_cost[ac] for ac in absorber_accounts)
        )
    
    self.costing.absorber_packing_cost = pyo.Expression(
        expr=sum(self.b23.total_plant_cost[ac] for ac in absorber_packing_accounts)
        )
    
    # self.costing.absorber_ic1_cost = pyo.Expression(
    #     expr=sum(self.b24.total_plant_cost[ac] for ac in absorber_ic1_accounts)
    #     )
    
    # self.costing.absorber_ic2_cost = pyo.Expression(
    #     expr=sum(self.b25.total_plant_cost[ac] for ac in absorber_ic2_accounts)
    #     )
    
    # self.costing.absorber_ic1pump_cost = pyo.Expression(
    #     expr=sum(self.b26.total_plant_cost[ac] for ac in absorber_ic1pump_accounts)
    #     )
    
    # self.costing.absorber_ic2pump_cost = pyo.Expression(
    #     expr=sum(self.b27.total_plant_cost[ac] for ac in absorber_ic2pump_accounts)
    #     )
    
    self.costing.ww_pump_cost = pyo.Expression(
        expr=sum(self.b28.total_plant_cost[ac] for ac in ww_pump_accounts)
        )
    
    self.costing.ww_cooler_cost = pyo.Expression(
        expr=sum(self.b29.total_plant_cost[ac] for ac in ww_cooler_accounts)
        )
    
    self.costing.rich_solvent_pump_cost = pyo.Expression(
        expr=sum(self.b30.total_plant_cost[ac] for ac in rich_solvent_pump_account)
        )
    
    self.costing.lean_rich_hex_cost = pyo.Expression(
        expr=sum(self.b31.total_plant_cost[ac] for ac in lean_rich_hex_accounts)
        )
    
    self.costing.stripper_column_cost = pyo.Expression(
        expr=sum(self.b32.total_plant_cost[ac] for ac in stripper_accounts)
        )
    
    self.costing.stripper_packing_cost = pyo.Expression(
        expr=sum(self.b33.total_plant_cost[ac] for ac in stripper_packing_accounts)
        )
    
    self.costing.stripper_condenser_cost = pyo.Expression(
        expr=sum(self.b34.total_plant_cost[ac] for ac in stripper_condenser_accounts)
        )
    
    self.costing.stripper_reboiler_cost = pyo.Expression(
        expr=sum(self.b35.total_plant_cost[ac] for ac in stripper_reboiler_accounts)
        )
    
    self.costing.stripper_reflux_drum_cost = pyo.Expression(
        expr=sum(self.b36.total_plant_cost[ac] for ac in stripper_reflux_drum_accounts)
        )
    
    self.costing.lean_solvent_pump_cost = pyo.Expression(
        expr=sum(self.b37.total_plant_cost[ac] for ac in lean_solvent_pump_account)
        )
    
    # self.costing.lean_solvent_cooler_cost = pyo.Expression(
    #     expr=sum(self.b38.total_plant_cost[ac] for ac in lean_solvent_cooler_accounts)
    #     )
    
    self.costing.solvent_stripper_reclaimer_cost = pyo.Expression(
        expr=sum(self.b39.total_plant_cost[ac] for ac in solvent_stripper_reclaimer_account)
        )
    
    self.costing.solvent_filtration_cost = pyo.Expression(
        expr=sum(self.b40.total_plant_cost[ac] for ac in solvent_filtration_account)
        )
    
    self.costing.Stripper_Reflux_Pump_cost = pyo.Expression(
        expr=sum(self.b42.total_plant_cost[ac] for ac in Stripper_Reflux_Pump_account)
        )
    
    self.costing.LS_Pump2_cost = pyo.Expression(
        expr=sum(self.b43.total_plant_cost[ac] for ac in LS_Pump2_account)
        )
    
    self.costing.Solvent_Sump_cost = pyo.Expression(
        expr=sum(self.b44.total_plant_cost[ac] for ac in Solvent_Sump_account)
        )
    
    self.costing.Solvent_Sump_pump_cost = pyo.Expression(
        expr=sum(self.b45.total_plant_cost[ac] for ac in Solvent_Sump_pump_account)
        )
    
    self.costing.Solvent_Sump_Filter_cost = pyo.Expression(
        expr=sum(self.b46.total_plant_cost[ac] for ac in Solvent_Sump_Filter_account)
        )
    
    self.costing.Solvent_Sump_Pit_Pump_cost = pyo.Expression(
        expr=sum(self.b47.total_plant_cost[ac] for ac in Solvent_Sump_Pit_Pump_account)
        )
    
    self.costing.Pre_Scrubber_Pump_cost = pyo.Expression(
        expr=sum(self.b48.total_plant_cost[ac] for ac in Pre_Scrubber_Pump_account)
        )
    
    self.costing.Reboiler_condensate_pot_cost = pyo.Expression(
        expr=sum(self.b49.total_plant_cost[ac] for ac in Reboiler_condensate_pot_account)
        )
    
    self.costing.DCC_Water_Filter_cost = pyo.Expression(
        expr=sum(self.b50.total_plant_cost[ac] for ac in DCC_Water_Filter_account)
        )
    
    self.costing.Corrosion_Inhib_Package_cost = pyo.Expression(
        expr=sum(self.b51.total_plant_cost[ac] for ac in Corrosion_Inhib_Package_account)
        )
    
    self.costing.Antifoam_Feed_Package_cost = pyo.Expression(
        expr=sum(self.b52.total_plant_cost[ac] for ac in Antifoam_Feed_Package_account)
        )

    self.costing.CCS_TPC = pyo.Expression(
        expr=
        (pyo.units.convert(
            self.costing.FG_blower_cost +
            self.costing.DCC_column_cost +
            self.costing.DCC_packing_cost +
            self.costing.DCC_pump_cost +
            self.costing.DCC_cooler_cost +
            self.costing.absorber_column_cost +
            self.costing.absorber_packing_cost +
            # self.costing.absorber_ic1_cost +
            # self.costing.absorber_ic2_cost +
            # self.costing.absorber_ic1pump_cost +
            # self.costing.absorber_ic2pump_cost +
            self.costing.ww_pump_cost +
            self.costing.ww_cooler_cost +
            self.costing.rich_solvent_pump_cost +
            self.costing.lean_rich_hex_cost +
            self.costing.stripper_column_cost +
            self.costing.stripper_packing_cost +
            self.costing.stripper_condenser_cost +
            self.costing.stripper_reboiler_cost +
            self.costing.stripper_reflux_drum_cost +
            self.costing.lean_solvent_pump_cost +
            # self.costing.lean_solvent_cooler_cost +
            self.costing.solvent_stripper_reclaimer_cost +
            self.costing.solvent_filtration_cost +
            self.costing.Stripper_Reflux_Pump_cost +
            self.costing.LS_Pump2_cost +
            self.costing.Solvent_Sump_cost +
            self.costing.Solvent_Sump_pump_cost +
            self.costing.Solvent_Sump_Filter_cost +
            self.costing.Solvent_Sump_Pit_Pump_cost +
            self.costing.Pre_Scrubber_Pump_cost +
            self.costing.Reboiler_condensate_pot_cost +
            self.costing.DCC_Water_Filter_cost +
            self.costing.Corrosion_Inhib_Package_cost +
            self.costing.Antifoam_Feed_Package_cost,
            to_units=getattr(pyunits, "USD_" + CE_index_year)  # $
            )
            )
        )

    self.costing.solvent_storage_tank_cost = pyo.Expression(
        expr=sum(self.b41.total_plant_cost[ac] for ac in solvent_storage_tank_account)
        )
    self.costing.NaOH_Makeup_Pump_cost = pyo.Expression(
        expr=sum(self.b53.total_plant_cost[ac] for ac in NaOH_Makeup_Pump_account)
        )
    self.costing.Solvent_Makeup_Pump_cost = pyo.Expression(
        expr=sum(self.b54.total_plant_cost[ac] for ac in Solvent_Makeup_Pump_account)
        )
    self.costing.NaOH_Storage_Tank_cost = pyo.Expression(
        expr=sum(self.b55.total_plant_cost[ac] for ac in NaOH_Storage_Tank_account)
        )
    self.costing.Station_Service_Equipment_cost = pyo.Expression(
        expr=sum(self.b56.total_plant_cost[ac] for ac in Station_Service_Equipment_account)
        )
    self.costing.Switchgear_Motor_Control_cost = pyo.Expression(
        expr=sum(self.b57.total_plant_cost[ac] for ac in Switchgear_Motor_Control_account)
        )
    self.costing.Conduit_Cable_Tray_cost = pyo.Expression(
        expr=sum(self.b58.total_plant_cost[ac] for ac in Conduit_Cable_Tray_account)
        )
    self.costing.Wire_Cable_cost = pyo.Expression(
        expr=sum(self.b59.total_plant_cost[ac] for ac in Wire_Cable_account)
        )
    self.costing.Main_Power_Transformers_cost = pyo.Expression(
        expr=sum(self.b60.total_plant_cost[ac] for ac in Main_Power_Transformers_account)
        )
    self.costing.Electrical_Foundations_cost = pyo.Expression(
        expr=sum(self.b61.total_plant_cost[ac] for ac in Electrical_Foundations_account)
        )
    self.costing.Control_Boards_Panels_Racks_cost = pyo.Expression(
        expr=sum(self.b62.total_plant_cost[ac] for ac in Control_Boards_Panels_Racks_account)
        )
    self.costing.Distributed_Control_System_Equipment_cost = pyo.Expression(
        expr=sum(self.b63.total_plant_cost[ac] for ac in Distributed_Control_System_Equipment_account)
        )
    self.costing.Instrument_Wiring_Tubing_cost = pyo.Expression(
        expr=sum(self.b64.total_plant_cost[ac] for ac in Instrument_Wiring_Tubing_account)
        )
    self.costing.Other_Instrumentation_Controls_Equipment_cost = pyo.Expression(
        expr=sum(self.b65.total_plant_cost[ac] for ac in Other_Instrumentation_Controls_Equipment_account)
        )
    self.costing.Other_Buildings_Structures_cost = pyo.Expression(
        expr=sum(self.b66.total_plant_cost[ac] for ac in Other_Buildings_Structures_account)
        )
    self.costing.Foundations_cost = pyo.Expression(
        expr=sum(self.b67.total_plant_cost[ac] for ac in Foundations_accounts)
        )
    self.costing.Interconnecting_Piping_cost = pyo.Expression(
        expr=sum(self.b68.total_plant_cost[ac] for ac in Interconnecting_Piping_account)
        )
    
    self.costing.RemovalSystemEquipment_TPC = pyo.Expression(
        expr=(
            (NoTrain*self.costing.CCS_TPC) +
            pyo.units.convert(
                self.costing.solvent_storage_tank_cost +
                self.costing.NaOH_Makeup_Pump_cost +
                self.costing.Solvent_Makeup_Pump_cost +
                self.costing.NaOH_Storage_Tank_cost +
                self.costing.Foundations_cost +
                self.costing.Interconnecting_Piping_cost,
                to_units=getattr(pyunits, "USD_" + CE_index_year)
                )
            )
        ) # $

    # Site Improvements
    self.costing.CO2CapTPC_ref = pyo.Param(initialize=337408894.9, mutable=False,
                                           units=pyo.units.USD_2018) # 2018 $
    self.costing.Site_Improvements_BEC_ref = pyo.Param(initialize=4287211.38811096, mutable=False,
                                                       units=pyo.units.USD_2021) # 2021 $
    self.costing.Site_Improvements_BEC = pyo.Expression(
        expr=pyo.units.convert(self.costing.Site_Improvements_BEC_ref, to_units=pyo.units.USD_2018)*\
        (pyo.units.convert(self.costing.RemovalSystemEquipment_TPC, to_units=pyo.units.USD_2018)/self.costing.CO2CapTPC_ref)**0.2
        )  # $ 2018
    self.costing.Site_Improvements_EPCC = pyo.Expression(expr=0.175*self.costing.Site_Improvements_BEC)
    self.costing.Site_Improvements_ProcCont = pyo.Expression(expr=0.17*self.costing.Site_Improvements_BEC)
    self.costing.Site_Improvements_ProjCont = pyo.Expression(expr=0.175*(
        self.costing.Site_Improvements_BEC + self.costing.Site_Improvements_EPCC + self.costing.Site_Improvements_ProcCont))
    self.costing.Site_Improvements_cost = pyo.Expression(expr=(
        pyunits.convert(self.costing.Site_Improvements_BEC + self.costing.Site_Improvements_EPCC +
                        self.costing.Site_Improvements_ProcCont + self.costing.Site_Improvements_ProjCont,
                        to_units=pyo.units.MUSD_2018))) # MM USD 2018
    
    self.costing.RemovalSystem_TPC = pyo.Expression(
        expr=
            self.costing.RemovalSystemEquipment_TPC +
            pyo.units.convert(
                self.costing.Station_Service_Equipment_cost +
                self.costing.Switchgear_Motor_Control_cost +
                self.costing.Conduit_Cable_Tray_cost +
                self.costing.Wire_Cable_cost +
                self.costing.Main_Power_Transformers_cost +
                self.costing.Electrical_Foundations_cost +
                self.costing.Control_Boards_Panels_Racks_cost +
                self.costing.Distributed_Control_System_Equipment_cost +
                self.costing.Instrument_Wiring_Tubing_cost +
                self.costing.Other_Instrumentation_Controls_Equipment_cost +
                self.costing.Other_Buildings_Structures_cost + 
                self.costing.Site_Improvements_cost,
                to_units=getattr(pyunits, "USD_" + CE_index_year)
                ) # $
            )


    # Obtain the total plant costs for each account - no capture case

    self.costing.NGCC_TPC_nocap = pyo.Expression(
        expr=
        (pyo.units.convert(
            sum(self.b1_nocap.total_plant_cost[ac] for ac in FW_accounts) +
            sum(self.b2_nocap.total_plant_cost[ac] for ac in RW_withdraw_accounts) +
            sum(self.b3_nocap.total_plant_cost[ac] for ac in FuelG_accounts) +
            sum(self.b4_nocap.total_plant_cost[ac] for ac in PW_discharge_accounts) +
            sum(self.b5_nocap.total_plant_cost[ac] for ac in FG_accounts) +
            sum(self.b6_nocap.total_plant_cost[ac] for ac in CT_grosspower_accounts) +
            sum(self.b7_nocap.total_plant_cost[ac] for ac in HRSG_duty_accounts) +
            sum(self.b8_nocap.total_plant_cost[ac] for ac in Stack_flow_gas_accounts) +
            sum(self.b9_nocap.total_plant_cost[ac] for ac in Steam_turbine_gross_power_accounts) +
            sum(self.b10_nocap.total_plant_cost[ac] for ac in Condenser_duty_accounts) +
            sum(self.b11_nocap.total_plant_cost[ac] for ac in Cooling_tower_accounts) +
            sum(self.b12_nocap.total_plant_cost[ac] for ac in Circ_water_accounts) +
            sum(self.b13_nocap.total_plant_cost[ac] for ac in plant_gross_power_accounts) +
            sum(self.b14_nocap.total_plant_cost[ac] for ac in auxilliary_load_accounts) +
            sum(self.b15_nocap.total_plant_cost[ac] for ac in stg_ctg_accounts) +
            sum(self.b16_nocap.total_plant_cost[ac] for ac in gasturbine_accounts),
            to_units=getattr(pyunits, "USD_" + CE_index_year)  # $
            )
            )
        )
    
    # MEA solvent initial fill cost
    self.costing.solvent_cost = pyo.Param(
        initialize=2.09,
        units=getattr(pyunits, "USD_" + CE_index_year)/pyo.units.kg) # $/kg
    self.costing.RemovalSystem_Equip_Adjust = pyo.Expression(
        expr=solvent_fill_init*self.costing.solvent_cost
        )

    # References: (B31B) (Exhibit 5-31, Rev 4 baseline report),
    # Exhibit 3-32 QGESS Capital Cost Scaling
    # Account 5.1 - CO2 Removal System
    self.costing.CC5_1 = pyo.Expression(
        expr=self.costing.RemovalSystem_TPC + self.costing.RemovalSystem_Equip_Adjust
        )
    # Account 5.4 CO2 Compression & Drying
    self.costing.CompAuxLoad_Base = pyo.Param(initialize=22918.0318, mutable=False, units=pyo.units.hp)  # hp
    self.costing.CC5_4 = pyo.Expression(
        expr=59674000*pyo.units.USD_2018*(CompAuxLoad_Cap/self.costing.CompAuxLoad_Base)**0.41
        )
    # Account 5.5 CO2 Compressor Aftercooler - assume its cost scales with CO2 capture flowrate
    self.costing.CO2capflow_Base = pyo.Param(initialize=223619.8061, mutable=False, units=pyo.units.kg/pyo.units.hr)  # kg/hr
    self.costing.CC5_5 = pyo.Expression(
        expr=498000*pyo.units.USD_2018*(pyo.units.convert(
            CO2CaptureRate,
            to_units=pyo.units.kg/pyo.units.hr
            )/self.costing.CO2capflow_Base)**0.6
        )
    # Account 5.12 Gas Cleanup Foundations
    self.costing.Reference_Capture = pyo.Param(initialize=493482.814, mutable=False, units=pyo.units.lb/pyo.units.hr)  # lb/hr
    self.costing.CC5_12 = pyo.Expression(
        expr=1145000*pyo.units.USD_2018*(CO2CaptureRate/self.costing.Reference_Capture)**0.79
        )

    # No of trains are included in the removal system
    self.costing.FG_Cleanup_TPC = pyo.Expression(
        expr= self.costing.CC5_1 + self.costing.CC5_4 + self.costing.CC5_5 + self.costing.CC5_12
        )

    # TPC of the integrated system
    self.costing.TPC = pyo.Expression(expr=self.costing.NGCC_TPC + self.costing.FG_Cleanup_TPC)

    # O&M cost calculation for NGCC with capture
    # Build cost constraints
    resources = ["natural_gas", "water", "water_treatment_chemicals",
                 "solvent", "SCR_catalyst", "ammonia", "triethylene_glycol",
                 "SCR_catalyst_waste", "triethylene_glycol_waste",
                 "amine_purification_unit waste", "thermal_reclaimer_unit_waste"]
    chemicals=["water_treatment_chemicals", "solvent", "SCR_catalyst",
                 "ammonia", "triethylene_glycol", ]
    self.costing.natural_gas_rate = pyo.Var(
        self.time,initialize=pyo.value(
            pyo.units.convert(fuelgas_value*22500*pyo.units.BTU/pyo.units.lb,  # heating value of 22500 BTU/lb
                              to_units=pyo.units.MBTU/pyo.units.hr)
            ),
        units=pyo.units.MBTU/pyo.units.hr)
    self.costing.water_rate = pyo.Var(
        self.time,initialize=pyo.value(
            pyo.units.convert(self.costing.raw_water_withdrawal,
                              to_units=pyo.units.gal/pyo.units.min)
            ),
        units=pyo.units.gal/pyo.units.min)
    self.costing.water_treatment_chemicals_rate = pyo.Var(
        self.time,initialize=pyo.value(
            pyo.units.convert(self.costing.raw_water_withdrawal/(4773*pyo.units.gal/pyo.units.min) # 4773 gpm is the BB ref value
                              * 10.2*pyo.units.ton/pyo.units.d,  # 10.2 ton/day is the BB ref value
                             to_units=pyo.units.ton/pyo.units.d)
            ),
        units=pyo.units.ton/pyo.units.d)
    self.costing.solvent_rate = pyo.Var(
        self.time,initialize=pyo.value(
            pyo.units.convert(solvent_makeup*CO2CaptureRate,
                              to_units=pyo.units.kg/pyo.units.hr)
            ),
        units=pyo.units.kg/pyo.units.hr)
    self.costing.SCR_catalyst_rate = pyo.Var(
        self.time,initialize=3.10,
        units=pyo.units.ft**3/pyo.units.d)
    self.costing.ammonia_rate = pyo.Var(
        self.time,initialize=3.5,
        units=pyo.units.ton/pyo.units.d)
    self.costing.triethylene_glycol_rate = pyo.Var(
        self.time,initialize=pyo.value(
            pyo.units.convert(CO2CaptureRate/(493588*pyo.units.lb/pyo.units.hr)  # 493588 lb/hr is the BB ref value
                              *394*pyo.units.gal/pyo.units.d,  # 394 gal/day is the BB ref value
                              to_units=pyo.units.gal/pyo.units.d)
            ),
        units=pyo.units.gal/pyo.units.d)
    self.costing.SCR_catalyst_waste_rate = pyo.Var(
        self.time,initialize=3.10,
        units=pyo.units.ft**3/pyo.units.d)
    self.costing.triethylene_glycol_waste_rate = pyo.Var(
        self.time,initialize=pyo.value(
            pyo.units.convert(CO2CaptureRate/(493588*pyo.units.lb/pyo.units.hr)  # 493588 lb/hr is the BB ref value
                              *394*pyo.units.gal/pyo.units.d,  # 394 gal/day is the BB ref value
                              to_units=pyo.units.gal/pyo.units.d)
            ),
        units=pyo.units.gal/pyo.units.d)
    self.costing.amine_purification_unit_waste_rate = pyo.Var(
        self.time,initialize=6.11,
        units=pyo.units.ton/pyo.units.d)
    self.costing.thermal_reclaimer_unit_waste_rate = pyo.Var(
        self.time,initialize=0.543,
        units=pyo.units.ton/pyo.units.d)
    rates = [self.costing.natural_gas_rate,
             self.costing.water_rate,
             self.costing.water_treatment_chemicals_rate,
             self.costing.solvent_rate,
             self.costing.SCR_catalyst_rate,
             self.costing.ammonia_rate,
             self.costing.triethylene_glycol_rate,
             self.costing.SCR_catalyst_waste_rate,
             self.costing.triethylene_glycol_waste_rate,
             self.costing.amine_purification_unit_waste_rate,
             self.costing.thermal_reclaimer_unit_waste_rate,]
    for var in rates:
        var.fix()
    prices = {"solvent": 2.09 * pyunits.USD_2018 / pyunits.kg}
    
    self.costing.total_plant_cost = pyo.Var(
        initialize=pyo.value(pyo.units.convert(self.costing.TPC, to_units=CE_index_units)), units=CE_index_units
        )
    self.costing.total_plant_cost.fix()

    print("\nSolving NGCC with capture costing model with operating cost variables...\n")
    results = solver.solve(self.costing, tee=True)

    self.ngcccap = QGESSCosting()
    self.ngcccap.total_plant_cost=pyo.Var(initialize=pyo.value(self.costing.total_plant_cost), units=CE_index_units)
    self.ngcccap.total_plant_cost.fix()
    self.ngcccap.net_power=pyo.Var(self.time, initialize=pyo.value(self.costing.plant_net_power), units=pyo.units.MW)
    self.ngcccap.net_power.fix()
    self.ngcccap.plant_gross_power=pyo.Var(initialize=pyo.value(self.costing.plant_gross_power), units=pyo.units.MW)
    self.ngcccap.plant_gross_power.fix()
    self.ngcccap.tonneCO2cap=pyo.Var(initialize=pyo.value(
        pyo.units.convert(CO2CaptureRate*8760*pyo.units.hr/pyo.units.year, to_units=pyo.units.tonne/pyo.units.year)
        ),
        units=pyo.units.tonne/pyo.units.year)
    self.ngcccap.tonneCO2cap.fix()
    self.ngcccap.land_cost=pyo.Var(initialize=3000*100*1e-6, units=CE_index_units)
    self.ngcccap.land_cost.fix()
    self.ngcccap.transport_cost_per_tonne=pyo.Var(initialize=10*1e-6, units=CE_index_units/pyo.units.tonne)
    self.ngcccap.transport_cost_per_tonne.fix()

    self.ngcccap.build_process_costs(
        total_plant_cost=self.ngcccap.total_plant_cost,
        nameplate_capacity=self.ngcccap.plant_gross_power,
        labor_rate=38.5,
        labor_burden=30,
        operators_per_shift=6.3,
        tech=6,
        land_cost=self.ngcccap.land_cost,
        net_power=self.ngcccap.net_power,
        fixed_OM=True,
        variable_OM=True,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
        chemicals=["water_treatment_chemicals", "solvent", "SCR_catalyst",
                   "ammonia", "triethylene_glycol", ],
        chemicals_inventory = ["water_treatment_chemicals", "solvent", "ammonia", "triethylene_glycol"],
        waste=["SCR_catalyst_waste", "triethylene_glycol_waste",
               "amine_purification_unit waste", "thermal_reclaimer_unit_waste", ],
        transport_cost=self.ngcccap.transport_cost_per_tonne,
        tonne_CO2_capture=self.ngcccap.tonneCO2cap,
        CE_index_year=CE_index_year,
    )

    # ADJUSTING VARIABLE BOUNDS FOR CONVERGENCE, MAY NEED TO ADD MORE HERE IN THE FUTURE
    self.ngcccap.total_fixed_OM_cost.setub(1e5)

    # Initialize costing
    print("\nNGCC With Cap block built, initializing...\n")
    QGESSCostingData.initialize_fixed_OM_costs(self.ngcccap)
    QGESSCostingData.initialize_variable_OM_costs(self.ngcccap)
    
    # Get default solver for testing
    solver = get_solver()

    # Solve the model
    print("\nNGCC With Cap block initialized, solving...\n")
    results = solver.solve(self.ngcccap, tee=True)

    # O&M cost calculation for NGCC without capture
    # Build cost constraints
    resources_nocap = ["natural_gas", "water", "water_treatment_chemicals",
                 "SCR_catalyst", "ammonia", "SCR_catalyst_waste",]
    chemicals_nocap=["water_treatment_chemicals", "SCR_catalyst",
                   "ammonia", ]
    self.costing.natural_gas_rate_nocap = pyo.Var(
        self.time,initialize=205630*22500*1e-6,
        units=pyo.units.MBTU/pyo.units.hr)
    self.costing.water_rate_nocap = pyo.Var(
        self.time,initialize=2902,
        units=pyo.units.gal/pyo.units.min)
    self.costing.water_treatment_chemicals_rate_nocap = pyo.Var(
        self.time,initialize=6.22,
        units=pyo.units.ton/pyo.units.d)
    self.costing.SCR_catalyst_rate_nocap = pyo.Var(
        self.time,initialize=3.10,
        units=pyo.units.ft**3/pyo.units.d)
    self.costing.ammonia_rate_nocap = pyo.Var(
        self.time,initialize=3.5,
        units=pyo.units.ton/pyo.units.d)
    self.costing.SCR_catalyst_waste_rate_nocap = pyo.Var(
        self.time,initialize=3.10,
        units=pyo.units.ft**3/pyo.units.d)
    rates_nocap = [self.costing.natural_gas_rate_nocap,
                   self.costing.water_rate_nocap,
                   self.costing.water_treatment_chemicals_rate_nocap,
                   self.costing.SCR_catalyst_rate_nocap,
                   self.costing.ammonia_rate_nocap,
                   self.costing.SCR_catalyst_waste_rate_nocap]
    for var in rates_nocap:
        var.fix()

    self.costing.total_plant_cost_nocap = pyo.Var(
        initialize=pyo.value(pyo.units.convert(self.costing.NGCC_TPC_nocap, to_units=CE_index_units)), units=CE_index_units
        )
    self.costing.total_plant_cost_nocap.fix()

    print("\nSolving NGCC no capture costing model with operating cost variables...\n")
    results = solver.solve(self.costing, tee=True)
    
    self.ngccnocap = QGESSCosting()
    self.ngccnocap.total_plant_cost=pyo.Var(initialize=pyo.value(self.costing.total_plant_cost_nocap), units=CE_index_units)
    self.ngccnocap.total_plant_cost.fix()
    self.ngccnocap.net_power=pyo.Var(self.time, initialize=pyo.value(self.costing.plant_net_power), units=pyo.units.MW)
    self.ngccnocap.net_power.fix()
    self.ngccnocap.plant_gross_power=pyo.Var(initialize=pyo.value(self.costing.plant_gross_power), units=pyo.units.MW)
    self.ngccnocap.plant_gross_power.fix()
    self.ngccnocap.land_cost=pyo.Var(initialize=3000*100*1e-6, units=CE_index_units)
    self.ngccnocap.land_cost.fix()
    self.ngccnocap.netpower=pyo.Var(initialize=727, units=pyo.units.MW)

    self.ngccnocap.build_process_costs(
        total_plant_cost=self.ngccnocap.total_plant_cost,
        nameplate_capacity=740000/1000,
        labor_rate=38.5,
        labor_burden=30,
        operators_per_shift=5,
        tech=6,
        land_cost=self.ngccnocap.land_cost,
        net_power=self.ngccnocap.netpower,
        fixed_OM=True,
        variable_OM=True,
        resources=resources_nocap,
        rates=rates_nocap,
        prices=dict(),
        fuel="natural_gas",
        chemicals=["water_treatment_chemicals", "SCR_catalyst",
                   "ammonia", ],
        chemicals_inventory=["water_treatment_chemicals", "ammonia"],
        waste=["SCR_catalyst_waste",],
        CE_index_year=CE_index_year,
    )


    # Initialize costing
    print("\nNGCC No Cap block built, initializing...\n")
    QGESSCostingData.initialize_fixed_OM_costs(self.ngccnocap)
    QGESSCostingData.initialize_variable_OM_costs(self.ngccnocap)
    
    # Get default solver for testing
    solver = get_solver()

    # Solve the model
    print("\nNGCC No Cap block initialized, solving...\n")
    results = solver.solve(self.ngccnocap, tee=True)

    # LCOE Calculation
    LCOE_costcomp_units = getattr(pyunits, "USD_" + CE_index_year)
    LCOE_units = getattr(pyunits, "USD_" + CE_index_year)/(pyunits.MW * pyunits.hr)
    CO2Cost_units = getattr(pyunits, "USD_" + CE_index_year)/(pyunits.tonne)
    
    # LCOE Components (NGCC with capture)

    self.costing.cap_factor = pyo.Var(initialize=self.ngcccap.capacity_factor, units=pyo.units.get_units(self.ngcccap.capacity_factor))
    self.costing.capital_cost_orig = pyo.Var(initialize=self.ngcccap.annualized_cost, units=pyo.units.get_units(self.ngcccap.annualized_cost))
    self.costing.TOC_fuel_supply_twomonths = (self.ngcccap.fuel_cost_OC * (2/2.25) * pyunits.year)
    self.costing.chemicals_non_init_fill = ["water_treatment_chemicals","ammonia", "triethylene_glycol","solvent"]
    self.costing.TOC_chemicals_non_init_fill = sum(self.ngcccap.variable_operating_costs[0,i] for i in self.costing.chemicals_non_init_fill)*1*pyunits.year/2
    self.costing.capital_cost = self.costing.capital_cost_orig + 0.0707*1.093*(- self.costing.TOC_fuel_supply_twomonths - self.costing.TOC_chemicals_non_init_fill)
    self.costing.fixed_cost = pyo.Var(initialize=self.ngcccap.total_fixed_OM_cost, units=pyo.units.get_units(self.ngcccap.total_fixed_OM_cost))
    self.costing.variable_cost = self.ngcccap.total_variable_OM_cost[0] * self.costing.cap_factor

    self.costing.capturecost = self.costing.capital_cost + self.costing.fixed_cost + self.costing.variable_cost
    self.costing.tonnescaptured = pyo.Var(initialize=self.ngcccap.tonne_CO2_capture, units=pyo.units.get_units(self.ngcccap.tonne_CO2_capture))

    self.costing.fuel_cost = pyo.Var(
        initialize=pyo.value(self.ngcccap.fuel_cost_OC * 12/2.25),
        units=pyo.units.get_units(self.ngcccap.fuel_cost_OC))  # converting back from 2.25 months to 1 year
    self.costing.transport_cost = pyo.Var(initialize=self.ngcccap.transport_cost, units=pyo.units.get_units(self.ngcccap.transport_cost))

    self.costing.capital_lcoe = self.costing.capital_cost/(self.costing.cap_factor * self.costing.plant_net_power * 8760 * pyunits.hr)
    self.costing.capital_lcoe = pyo.Var(
        initialize=pyo.value(pyunits.convert(self.costing.capital_lcoe, LCOE_units)),
        units=LCOE_units)

    self.costing.fixed_lcoe = self.costing.fixed_cost/(self.costing.cap_factor * self.costing.plant_net_power * 8760 * pyunits.hr)
    self.costing.fixed_lcoe = pyo.Var(
        initialize=pyo.value(pyunits.convert(self.costing.fixed_lcoe, LCOE_units)),
        units=LCOE_units)

    self.costing.variable_lcoe = self.costing.variable_cost/(self.costing.cap_factor * self.costing.plant_net_power * 8760 * pyunits.hr)
    self.costing.variable_lcoe = pyo.Var(
        initialize=pyo.value(pyunits.convert( self.costing.variable_lcoe * 1 * pyunits.year, LCOE_units)),
        units=LCOE_units)

    # fuel cost is included in total variable costs; reporting it individually here
    self.costing.fuel_lcoe = self.costing.fuel_cost/(self.costing.plant_net_power * 8760 * pyunits.hr)
    self.costing.fuel_lcoe = pyo.Var(
        initialize=pyo.value(pyunits.convert( self.costing.fuel_lcoe * 1 * pyunits.year, LCOE_units)),
        units=LCOE_units)

    # reporting nonfuel variable cost explicitly
    self.costing.nonfuel_variable_lcoe = pyo.Var(
        initialize=pyo.value(self.costing.variable_lcoe - self.costing.fuel_lcoe),
        units=LCOE_units)

    self.costing.transport_lcoe = self.costing.transport_cost/(self.costing.cap_factor * self.costing.plant_net_power * 8760 * pyunits.hr)
    self.costing.transport_lcoe = pyo.Var(
        initialize=pyo.value(pyunits.convert( self.costing.transport_lcoe, LCOE_units)),
        units=LCOE_units)

    self.costing.LCOE = pyo.Var(
        initialize=pyo.value(self.costing.capital_lcoe + self.costing.fixed_lcoe + self.costing.variable_lcoe + self.costing.transport_lcoe),
        units=LCOE_units)
    
    # LCOE Components (NGCC without capture)

    self.costing.cap_factor_nocap = pyo.Var(initialize=self.ngccnocap.capacity_factor, units=pyo.units.get_units(self.ngccnocap.capacity_factor))
    self.costing.capital_cost_nocap_orig = pyo.Var(initialize=self.ngccnocap.annualized_cost, units=pyo.units.get_units(self.ngccnocap.annualized_cost))
    self.costing.TOC_fuel_supply_twomonths_nocap = (self.ngccnocap.fuel_cost_OC * (2/2.25) * pyunits.year)
    self.costing.chemicals_non_init_fill_nocap = ["water_treatment_chemicals","ammonia"]
    self.costing.TOC_chemicals_non_init_fill_nocap = sum(self.ngccnocap.variable_operating_costs[0,i] for i in self.costing.chemicals_non_init_fill_nocap)*1*pyunits.year/2
    self.costing.capital_cost_nocap = self.costing.capital_cost_nocap_orig + 0.0707*1.093*(- self.costing.TOC_fuel_supply_twomonths_nocap - self.costing.TOC_chemicals_non_init_fill_nocap)
    self.costing.fixed_cost_nocap = pyo.Var(initialize=self.ngccnocap.total_fixed_OM_cost, units=pyo.units.get_units(self.ngccnocap.total_fixed_OM_cost))
    self.costing.variable_cost_nocap = self.ngccnocap.total_variable_OM_cost[0] * self.costing.cap_factor_nocap

    self.costing.tonnescaptured_nocap = pyo.Var(initialize=0, units=pyo.units.tonne)

    self.costing.fuel_cost_nocap = pyo.Var(
        initialize=pyo.value(self.ngccnocap.fuel_cost_OC * 12/2.25),
        units=pyo.units.get_units(self.ngccnocap.fuel_cost_OC))  # converting back from 2.25 months to 1 year
    self.costing.plant_net_power_nocap = pyo.Var(initialize=727, units=pyo.units.MW)

    self.costing.capital_lcoe_nocap = self.costing.capital_cost_nocap/(self.costing.cap_factor * self.costing.plant_net_power_nocap * 8760 * pyunits.hr)
    self.costing.capital_lcoe_nocap = pyo.Var(
        initialize=pyo.value(pyunits.convert(self.costing.capital_lcoe_nocap, LCOE_units)),
        units=LCOE_units)

    self.costing.fixed_lcoe_nocap = self.costing.fixed_cost_nocap/(self.costing.cap_factor * self.costing.plant_net_power_nocap * 8760 * pyunits.hr)
    self.costing.fixed_lcoe_nocap = pyo.Var(
        initialize=pyo.value(pyunits.convert(self.costing.fixed_lcoe_nocap, LCOE_units)),
        units=LCOE_units)

    self.costing.variable_lcoe_nocap = self.costing.variable_cost_nocap/(self.costing.cap_factor * self.costing.plant_net_power_nocap * 8760 * pyunits.hr)
    self.costing.variable_lcoe_nocap = pyo.Var(
        initialize=pyo.value(pyunits.convert(self.costing.variable_lcoe_nocap * 1 * pyunits.year, LCOE_units)),
        units=LCOE_units)

    # fuel cost is included in total variable costs; reporting it individually here
    self.costing.fuel_lcoe_nocap = self.costing.fuel_cost_nocap/(self.costing.plant_net_power_nocap * 8760 * pyunits.hr)
    self.costing.fuel_lcoe_nocap = pyo.Var(
        initialize=pyo.value(pyunits.convert(self.costing.fuel_lcoe_nocap * 1 * pyunits.year, LCOE_units)),
        units=LCOE_units)

    # reporting nonfuel variable cost explicitly
    self.costing.nonfuel_variable_lcoe_nocap = pyo.Var(
        initialize=pyo.value(self.costing.variable_lcoe_nocap - self.costing.fuel_lcoe_nocap),
        units=LCOE_units)

    self.costing.transport_lcoe_nocap = pyo.Var(initialize=0, units=LCOE_units)

    self.costing.LCOE_nocap = pyo.Var(
        initialize=pyo.value(self.costing.capital_lcoe_nocap + self.costing.fixed_lcoe_nocap + self.costing.variable_lcoe_nocap + self.costing.transport_lcoe_nocap),
        units=LCOE_units)
    
    # $/tonne CO2
    # Cost of CO2 Capture
    self.costing.Cost_of_capture = pyo.Var(
        initialize=pyo.value(
            pyo.units.convert(
                (self.costing.LCOE-self.costing.transport_lcoe-self.costing.LCOE_nocap)
                *self.costing.plant_net_power/CO2CaptureRate,
            to_units=CO2Cost_units)
            ),
        units=CO2Cost_units)
    # Cost of CO2 Avoided
    self.costing.Cost_CO2_avoided = pyo.Var(
        initialize=pyo.value(
            pyo.units.convert(
                (self.costing.LCOE-self.costing.LCOE_nocap)
                /(
                    CO2_emission_nocap/self.costing.plant_net_power_nocap-
                    CO2_emission_cap/self.costing.plant_net_power),
                to_units=CO2Cost_units)
            ),
        units=CO2Cost_units)

    print("\nSolve a final time to calculate results variables...\n")
    results = solver.solve(self.costing, tee=True)

# -------------------------------------------------------------------------------------------
    if export_economic_results is True:
        DNE = "COMPONENT DOES NOT EXIST"
        CCS_cost_summary_header = ['FG Blower','DCC Column','DCC Packing','DCC Pump','DCC Cooler','Absorber Column',
    	                           'Absorber Packing','Absorber Intercooler 1','Absorber Intercooler 2','Absorber Intercooler Pump 1',
    							   'Absorber Intercooler Pump 2','WW Pump','WW Cooler','Rich Solvent Pump','Lean Rich Heat Exchanger',
    							   'Stripper Column','Stripper Packing','Stripper Condenser','Stripper Reboiler','Stripper Reflux Drum',
    							   'Lean Solvent Pump','Lean Solvent Cooler','Solvent Stripper Reclaimer','Solvent Filtration',
    							   'Stripper Reflux Pump','LS_Pump2','Solvent Sump','Solvent Sump Pump','Solvent Sump Filter','Solvent Sump Pit Pump',
    							   'Pre Scrubber Pump','Reboiler Condensate Pot','DCC Water Filter','Corrosion Inhibition','Antifoam Feed Package',
    							   'Solvent Storage Tank','NaOH Makeup','Solvent Makeup Pump','NaOH Storage Tank','Foundations','Interconnecting Piping',
    							   'Station Service Equipment','Switchgear Motor Control','Conduit Cable Tray','Wire Cable','Main Power Transformers',
    							   'Electrical Foundations','Control Boards and Panels Rack','Distributed Control System Equipment','Instrument Wiring Tubing',
    							   'Other Instrumentation Controls Equipment','Other Buildings and Structures','Site Improvements']
        CCS_cost_summary_objects = [self.costing.FG_blower_cost,self.costing.DCC_column_cost,self.costing.DCC_packing_cost,self.costing.DCC_pump_cost,self.costing.DCC_cooler_cost,
                                    self.costing.absorber_column_cost,self.costing.absorber_packing_cost,
                                    DNE,DNE,DNE,DNE,  # self.costing.absorber_ic1_cost,self.costing.absorber_ic2_cost,self.costing.absorber_ic1pump_cost,self.costing.absorber_ic2pump_cost,
                                    self.costing.ww_pump_cost,self.costing.ww_cooler_cost,
                                    self.costing.rich_solvent_pump_cost,self.costing.lean_rich_hex_cost,self.costing.stripper_column_cost,self.costing.stripper_packing_cost,
                                    self.costing.stripper_condenser_cost,self.costing.stripper_reboiler_cost,self.costing.stripper_reflux_drum_cost,self.costing.lean_solvent_pump_cost,
                                    DNE,  # self.costing.lean_solvent_cooler_cost,
                                    self.costing.solvent_stripper_reclaimer_cost,self.costing.solvent_filtration_cost,self.costing.Stripper_Reflux_Pump_cost,
                                    self.costing.LS_Pump2_cost,self.costing.Solvent_Sump_cost,self.costing.Solvent_Sump_pump_cost,self.costing.Solvent_Sump_Filter_cost,
                                    self.costing.Solvent_Sump_Pit_Pump_cost,self.costing.Pre_Scrubber_Pump_cost,self.costing.Reboiler_condensate_pot_cost,self.costing.DCC_Water_Filter_cost,
                                    self.costing.Corrosion_Inhib_Package_cost,self.costing.Antifoam_Feed_Package_cost,self.costing.solvent_storage_tank_cost,self.costing.NaOH_Makeup_Pump_cost,
                                    self.costing.Solvent_Makeup_Pump_cost,self.costing.NaOH_Storage_Tank_cost,self.costing.Foundations_cost,self.costing.Interconnecting_Piping_cost,
                                    self.costing.Station_Service_Equipment_cost,self.costing.Switchgear_Motor_Control_cost,self.costing.Conduit_Cable_Tray_cost,self.costing.Wire_Cable_cost,
                                    self.costing.Main_Power_Transformers_cost,self.costing.Electrical_Foundations_cost,self.costing.Control_Boards_Panels_Racks_cost,
                                    self.costing.Distributed_Control_System_Equipment_cost,self.costing.Instrument_Wiring_Tubing_cost,self.costing. Other_Instrumentation_Controls_Equipment_cost,
                                    self.costing.Other_Buildings_Structures_cost,self.costing.Site_Improvements_cost]
        CCS_cost_summary_data = [pyo.value(obj) for obj in CCS_cost_summary_objects]
        CCS_cost_summary_units = [pyo.units.get_units(obj) for obj in CCS_cost_summary_objects]
    
        if overwrite_economic_results is True or not os.path.exists(os.path.join(this_file_dir(),'CCS_Cost_Summary.csv')):
            with open(os.path.join(this_file_dir(),'CCS_Cost_Summary.csv'),'w',encoding='UTF8') as f:
                writer = csv.writer(f)
                writer.writerow(CCS_cost_summary_header)
                writer.writerow(CCS_cost_summary_data)
                writer.writerow(CCS_cost_summary_units)
    
    
        cost_summary_header = ['Capital_LCOE','Fixed_LCOE','Nonfuel Variable_LCOE','Fuel_LCOE', 'Variable_LCOE', 'Transport_LCOE','LCOE','Cost of Capture ($/tonne CO2)',
                               'Cost of CO2 Avoided ($/tonne CO2)','Removal System TPC','Removal System Equip Adjust','CO2 Compression & Drying',
                               'FG Cleanup TPC','NGCC TPC']
        cost_summary_obj = [self.costing.capital_lcoe,self.costing.fixed_lcoe,self.costing.nonfuel_variable_lcoe,self.costing.fuel_lcoe,self.costing.variable_lcoe,
                             self.costing.transport_lcoe,self.costing.LCOE,self.costing.Cost_of_capture,self.costing.Cost_CO2_avoided,
                             self.costing.RemovalSystem_TPC,self.costing.RemovalSystem_Equip_Adjust,self.costing.CC5_4,
                             self.costing.FG_Cleanup_TPC,self.costing.NGCC_TPC]
        cost_summary_data = [pyo.value(obj) for obj in cost_summary_obj]
        cost_summary_units = [pyo.units.get_units(obj) for obj in cost_summary_obj]
    
        if overwrite_economic_results is True or not os.path.exists(os.path.join(this_file_dir(),'Cost_Summary.csv')):
            with open(os.path.join(this_file_dir(), 'Cost_Summary.csv'),'w',encoding='UTF8') as f:
                writer = csv.writer(f)
                writer.writerow(cost_summary_header)
                writer.writerow(cost_summary_data)
                writer.writerow(cost_summary_units)
    
    		
        process_summary_header = ['Plant Net Power (MWe)','CO2 Emissions (kg/hr)','Plant Gross Power (MWe)','Auxiliary Load (MWe)']
        process_summary_obj = [pyo.units.convert(self.costing.plant_net_power, to_units=pyo.units.MW),
                               pyo.units.convert(CO2_emission_cap, to_units=pyo.units.kg/pyo.units.hr),
                               pyo.units.convert(self.costing.plant_gross_power, to_units=pyo.units.MW),
                               pyo.units.convert(self.costing.aux_load, to_units=pyo.units.MW)]
        process_summary_data = [pyo.value(obj) for obj in process_summary_obj]
        process_summary_units = [pyo.units.get_units(obj) for obj in process_summary_obj]
    
        if overwrite_economic_results is True or not os.path.exists(os.path.join(this_file_dir(),'Process_Summary.csv')):
            with open(os.path.join(this_file_dir(), 'Process_Summary.csv'),'w',encoding='UTF8') as f:
                writer = csv.writer(f)
                writer.writerow(process_summary_header)
                writer.writerow(process_summary_data)
                writer.writerow(process_summary_units)

    return
    # old return statement, for reference if more reporting is added later
    # return Cost_of_capture, Cost_CO2_avoided, LCOE, LCOE_nocap, capital_lcoe_nocap, capital_lcoe, fixed_lcoe_nocap, fixed_lcoe, variable_lcoe_nocap, variable_lcoe, fuel_lcoe_nocap, fuel_lcoe, transport_lcoe, RemovalSystem_TPC, plant_net_power, absorber_column_cost, absorber_packing_cost, stripper_column_cost, stripper_packing_cost, stripper_reboiler_cost, stripper_condenser_cost, lean_rich_hex_cost, FG_blower_cost, DCC_column_cost, DCC_packing_cost, DCC_pump_cost, DCC_cooler_cost, absorber_ic1_cost, absorber_ic2_cost, absorber_ic1pump_cost, absorber_ic2pump_cost, ww_pump_cost, ww_cooler_cost, rich_solvent_pump_cost, stripper_reflux_drum_cost, lean_solvent_pump_cost, lean_solvent_cooler_cost, solvent_stripper_reclaimer_cost, solvent_filtration_cost, solvent_storage_tank_cost
