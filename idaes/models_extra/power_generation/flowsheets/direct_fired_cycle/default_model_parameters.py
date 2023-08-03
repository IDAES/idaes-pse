

"""
This script contains parameter data

In the current formulation, the fixed O&M cost is calculated as a fraction of
CAPEX, and it includes: TBD. e.g., insurance, labor, etc

The variable O&M cost (non-fuel) is calculated as follows:
Non-fuel VOM includes: TBD.
Total VOM = 3339.92 + 18048.07 = 21387.99. Split the VOM for inidividual units
in the ratio of their CAPEX. 1128855 / (1128855 + 545522) = 0.6742.

The VOM of the DFC cycle = 0.6742 * 21387.99 = 14419.78286 
 ==> 14419.78286 / (0.85 * 8760) = 1.93658 [in $1000/hr]
     This number will be scaled linearly with normalized DFC's capacity

The VOM of the ASU unit = 0.3258 * 21387.99 = 6968.207
 ==> 6968.207 / (0.85 * 8760) = 0.93583 [in $1000/hr]
     This number will be scaled linearly with normalized ASU's capacity

The VOM of the NLU unit = 5404.96
 ==> 5404.96 / (0.85 * 8760) = 0.7259 [in $1000/hr]
     This number will be scaled linearly with normalized NLU's capacity
"""

# TODO: Contact Sandeep for reference. 
# HHV of NG = 22499.17034 btu/lb = 0.0496 MMBtu/kg
NG_HHV = 0.0496
HR_TO_SEC = 3600
LBS_TO_KG = 0.453592

# Parameter data for DFC
# Calculating the FOM as 3.157% of the CAPEX. The percentage value is obtained from
# 35,641.27 (FOM) / 1,128,855 (CAPEX). The FOM is also assumed 
# to vary linearly with the capacity.
DFC_PARAMS = {
    "dfc_capacity"  : 838.11322145,          # [MW] Net power output 
    "ng_flow"       : 28.87,                 # [kg/s] NG flow at full capacity 
    "o2_ng_ratio"   : 3.785632948,           # [-] Oxygen/NG flowratio
    "capex"         : 1128855,               # [$1000] CAPEX of DFC 
    "fom_factor"    : 0.03157,               # [-] Multiplier for FOM 
    "op_curve_coeff": [1.2845, -0.2845],     # Coefficients of performance curve
    "co2_emission"  : 3.06,                  # [kg/MWh] CO2 emissions 
    "vom"           : 1.93658,               # [$1000/hr] Non-fuel VOM 
}

# Parameter data for Monolithic ASU
MONO_ASU_PARAMS = {
    "asu_capacity"      : 109.2912232,       # [kg/s] Maximum O2 flow 
    "power_requirement" : 162.1878617,       # [MW] Power requirement at max capacity 
    "capex"             : 545522,            # [$1000] CAPEX of ASU 
    "fom_factor"        : 0.03157,           # [-] Multiplier for FOM
    "op_curve_coeff"    : [0.9625, 0.0375],  # Coefficients of performance curve
    "vom"               : 0.93583,           # [$1000/hr] Non-electricity VOM 
}

# Parameter data for ASU
ASU_PARAMS = {
    "asu_capacity"      : 27.3228058,        # [kg/s] Maximum O2 flow 
    "power_requirement" : 40.546965425,      # [MW] Power requirement at max capacity 
    "capex"             : 136380.5,          # [$1000] CAPEX of ASU 
    "fom_factor"        : 0.03157,           # [-] Multiplier for FOM 
    "op_curve_coeff"    : [0.9625, 0.0375],  # coefficients of performance curve
    "vom"               : 0.2339575,         # [$1000/hr] Non-electricity VOM 
}

# Parameter data for NLU
NLU_PARAMS = {
    "nlu_capacity"      : 109.291223201103,   # [kg/s] Maximum O2 flow 
    "power_requirement" : 143.8024,           # [MW] Power requirement at max capacity 
    "capex"             : 171189.6,           # [$1000] CAPEX of NLU 
    "fom_factor"        : 0.03157,            # [-] Multiplier for FOM [-]
    "vom"               : 0.7259,             # [$1000/hr] Non-electricity VOM 
}

# Parameter data for LOx tank
O2_TANK_PARAMS = {
    "capex"             : [0.98167214, 2779.90543786],
    "fom_factor"        : 0.03157,            # [-] Multiplier for FOM
    "min_holdup"        : 0.1,                # [-] Minimum holdup fraction
    "power_requirement" : 1.47054,            # [MW] Power requirement for storage
    "o2_flow"           : 112.0865459,        # [kg/s] O2 flowrate for 
    "tank_constraint"   : "periodic",         # "Periodic" or "initial_holdup"
}


# Cashflow parameters
CASHFLOW_PARAMS = {
    "plant_life"        : 30,                 # [-] Plant lifetime in years
    "discount_rate"     : 0.075,              # [-] Discount rate
    "tax_rate"          : 0.2,                # [-] Corporate tax rate
    "electricity_cost"  : 5,                  # [$/MWh] Excess penalty for purchasing electricity
}
