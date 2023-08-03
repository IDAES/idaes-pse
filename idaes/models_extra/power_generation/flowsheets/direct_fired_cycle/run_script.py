import os
import copy
from idaes.power_generation.flowsheets.direct_fired_cycle.npv_optimization import (
    npv_model_dfc_asu,
    npv_model_dfc_asu_with_lox,
    npv_model_dfc_asu_nlu,
)
import idaes.power_generation.flowsheets.direct_fired_cycle. \
    default_model_parameters as dmp

DFC_PARAMS = dmp.DFC_PARAMS
ASU_PARAMS = dmp.ASU_PARAMS
TANK_PARAMS = dmp.O2_TANK_PARAMS
NLU_PARAMS = dmp.NLU_PARAMS
CASHFLOW_PARAMS = dmp.CASHFLOW_PARAMS


# for tax in [100, 150]:
#     for region in ["CAISO", "ERCOT", "MISO-W", "NYISO", "PJM-W"]:
#         m, sol = npv_model_dfc_asu(location=region, carbon_tax=tax)
#         # m, sol = npv_model_dfc_asu_with_lox(location=region, carbon_tax=tax)
#         # m, sol = npv_model_dfc_asu_nlu(location=region, carbon_tax=tax)

# for tax in [100]:
#     for region in ["CAISO"]:
#         m, sol = npv_model_dfc_asu(location=region, carbon_tax=tax)

# for tax, region in [(100, "CAISO"), (150, "ERCOT"), (150, "PJM-W")]:
#     print("*" * 80)
#     print("Running ", region, tax)
#     print("*" * 80)
#     m, sol = npv_model_dfc_asu_nlu(location=region, carbon_tax=tax)


for tax, region in [(100, "CAISO")]:
    print("*" * 80)
    print("Running ", region, tax)
    print("*" * 80)

    dfc_params = copy.deepcopy(DFC_PARAMS)
    asu_params = copy.deepcopy(ASU_PARAMS)
    cost_params = copy.deepcopy(CASHFLOW_PARAMS)

    # Modify any default parameters if needed
    # dfc_params["capex"] /= 2

    cwd = os.getcwd()
    if not os.path.exists(cwd + "\\results"):
        os.mkdir(os.path.join(cwd, "results"))

    m, sol = npv_model_dfc_asu(
        location=region,
        carbon_tax=tax,
        dfc_params=dfc_params,
        asu_params=asu_params,
        cost_params=cost_params,
        folder="results/",
    )
