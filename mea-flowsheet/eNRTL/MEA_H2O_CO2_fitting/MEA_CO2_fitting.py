import sys, os, logging
sys.path.append('eNRTL')


import pandas as pd
import numpy as np
from scipy.linalg import svd, schur, eigh
import matplotlib.pyplot as plt

import pyomo.environ as pyo
from pyomo.common.fileutils import this_file_dir
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.contrib.interior_point.inverse_reduced_hessian import inv_reduced_hessian_barrier

import idaes.logger as idaeslog
from idaes.core import Solvent
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.initialization import solve_indexed_blocks
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.pure import NIST
from MEA_eNRTL import get_prop_dict, create_heat_capacity_no_inherent_rxns, initialize_inherent_reactions

def loss(x):
    return 0.5*x ** 2

def get_estimated_params(m):
    param_list = []
    for rxn_name in [
        # MEA bicarbonate parameters totally characterized from pKa data
        # of MEA and CO2 in H2O
        #"MEA_bicarbonate_formation", 
        "MEA_carbamate_formation_combo",
    ]:
        rxn_obj = getattr(m.params, "reaction_"+rxn_name)
        param_list.append(rxn_obj.k_eq_coeff_1)
        param_list.append(rxn_obj.k_eq_coeff_2)
    return param_list

def create_and_scale_params(m):
    rxn_combinations = {
        "MEA_bicarbonate_formation": {
            "bicarbonate_formation": 1,
            "MEA_protonation": 1
        },
        "MEA_carbamate_formation_combo": {
            "MEA_carbamate_formation": 1,
            "MEA_protonation": 1
        },
    }
    config = get_prop_dict(
        ["H2O", "MEA", "CO2"], 
        excluded_rxns=["H2O_autoionization", "carbonate_formation"],
        rxn_combinations=rxn_combinations
    )
    assert "H3O^+" not in config["components"]
    assert "OH^-" not in config["components"]
    params = m.params = GenericParameterBlock(**config)

    gsf = iscale.get_scaling_factor
    scaling_factor_flow_mol = 1/100
    params.set_default_scaling("enth_mol_phase", 3e-4)
    params.set_default_scaling("pressure", 1e-5)
    params.set_default_scaling("temperature", 1)
    params.set_default_scaling("flow_mol", scaling_factor_flow_mol)
    params.set_default_scaling("flow_mol_phase", scaling_factor_flow_mol)

    params.set_default_scaling("flow_mass_phase", scaling_factor_flow_mol / 18e-3)  # MW mixture ~= 24 g/Mol
    params.set_default_scaling("dens_mol_phase", 1 / 18000)
    params.set_default_scaling("visc_d_phase", 700)
    params.set_default_scaling("log_k_eq", 1)

    mole_frac_scaling_factors = {
        "H2O": 2,
        "MEA": 2,
        "CO2": 5
    }
    mole_frac_true_scaling_factors = {
        "MEA": 1e1,
        "H2O": 2,
        "HCO3^-": 5e3,
        # "H3O^+": 5e5,
        # "CO3^2-": 1e12,
        # "OH^-": 5e11,
        # "HCO3^-": 1e3,
        # "H3O^+": 1e11,
        # "CO3^2-": 1e3,
        # "OH^-": 1e3,
        "MEAH^+": 1e1,
        "MEACOO^-": 5e2,
        "CO2": 1e3

    }
    inherent_rxn_scaling_factors = {
        "MEA_bicarbonate_formation": 5e2,
        "MEA_carbamate_formation_combo": 5e2,
    }
    for comp, sf_x in mole_frac_scaling_factors.items():
        params.set_default_scaling("mole_frac_comp", sf_x, index=comp)
        params.set_default_scaling("mole_frac_phase_comp", sf_x, index=("Liq", comp))
        params.set_default_scaling(
            "flow_mol_phase_comp",
            sf_x * scaling_factor_flow_mol,
            index=("Liq", comp)
        )

    for comp, sf_x in mole_frac_true_scaling_factors.items():
        params.set_default_scaling("mole_frac_phase_comp_true", sf_x, index=("Liq", comp))
        params.set_default_scaling(
            "flow_mol_phase_comp_true",
            sf_x * scaling_factor_flow_mol,
            index=("Liq", comp)
        )
    for rxn, sf_xi in inherent_rxn_scaling_factors.items():
        params.set_default_scaling(
            "apparent_inherent_reaction_extent",
            scaling_factor_flow_mol*sf_xi,
            index=rxn
        )

    iscale.set_scaling_factor(m.params.Liq.alpha, 1) 
    iscale.set_scaling_factor(m.params.Liq.tau_A, 1) # Reminder that it's well-scaled
    iscale.set_scaling_factor(m.params.Liq.tau_B, 1/300)

    for rxn_name in inherent_rxn_scaling_factors.keys():
        rxn_obj = getattr(m.params, "reaction_"+rxn_name)
        iscale.set_scaling_factor(rxn_obj.k_eq_coeff_1, 1)
        iscale.set_scaling_factor(rxn_obj.k_eq_coeff_2, 1/300)
        iscale.set_scaling_factor(rxn_obj.k_eq_coeff_3, 1)
        iscale.set_scaling_factor(rxn_obj.k_eq_coeff_4, 300)


optarg = {
    # 'bound_push' : 1e-22,
    'nlp_scaling_method': 'user-scaling',
    'linear_solver': 'ma57',
    'OF_ma57_automatic_scaling': 'yes',
    'max_iter': 300,
    'tol': 1e-8,
    'constr_viol_tol': 1e-8,
    'halt_on_ampl_error': 'no',
    # 'mu_strategy': 'monotone',
}

if __name__ == "__main__":
    logging.getLogger('pyomo.repn.plugins.nl_writer').setLevel(logging.ERROR)
    init_outlevel = idaeslog.DEBUG

    data = os.sep.join([this_file_dir(), "data"])
    m = pyo.ConcreteModel()
    create_and_scale_params(m)

    obj_expr = 0

    # df_hillard_loading = pd.read_csv(os.sep.join([data, "hillard_PZ_loading_concatenated.csv"]), index_col=None)

    # n_data_hillard_loading = len(df_hillard_loading ["temperature"])
    # m.hillard_loading = m.params.build_state_block(range(n_data_hillard_loading), defined_state=True)

    # for i, row in df_hillard_loading.iterrows():
    #     molality=row["PZ_molality"]
    #     n_PZ = molality
    #     n_H2O = 1/0.01802
    #     # Convention in papers from Gary Rochelle's group to assume 1 mole PZ can dissolve 2 moles CO2
    #     n_CO2 = 2 * n_PZ * row["CO2_loading"]
    #     n_tot = n_PZ + n_H2O + n_CO2

    #     m.hillard_loading[i].flow_mol.fix(100)
    #     m.hillard_loading[i].mole_frac_comp["H2O"].fix(n_H2O/n_tot)
    #     m.hillard_loading[i].mole_frac_comp["PZ"].fix(n_PZ/n_tot)
    #     m.hillard_loading[i].mole_frac_comp["CO2"].fix(n_CO2/n_tot)
    #     # Got 45 psig as the average temperature from Cuillane and Rochelle (2005)
    #     # Assuming it's the same here, because I can't find Hillard's master thesis.
    #     # Total pressure doesn't matter at this stage
    #     m.hillard_loading[i].pressure.fix(310264+101300) 
    #     m.hillard_loading[i].temperature.fix(row["temperature"] + 273.15) # Temperature in C
    #     obj_expr +=  (
    #         # CO2 partial pressure
    #         loss(
    #             m.hillard_loading[i].log_fug_phase_comp["Liq","CO2"]
    #             - pyo.log(row["CO2_partial_pressure"] * 1e3)  # Pressure in kPa
    #         )
    #     )

    df_jou_loading = pd.read_csv(os.sep.join([data, "jou_et_al_1995_loading.csv"]), index_col=None)

    n_data_jou_loading = len(df_jou_loading ["temperature"])
    m.jou_loading = m.params.build_state_block(range(n_data_jou_loading), defined_state=True)

    # 30 weight percent MEA
    x_MEA_unloaded = 0.3/(61.084/18.02+0.3*(1-61.084/18.02))
    x_H2O_unloaded = 1 - x_MEA_unloaded

    for i, row in df_jou_loading.iterrows():
        n_MEA = 100*x_MEA_unloaded
        n_H2O = 100*x_H2O_unloaded
        # Two sources of measurement for CO2 loading, sometimes only one is present, and they can disagree.
        # Take the average now, probably better to work them into the objective function eventually
        if np.isnan(row["CO2_loading_GC"]):
            CO2_loading = row["CO2_loading_BaCO3"]
        elif np.isnan(row["CO2_loading_BaCO3"]):
            CO2_loading = row["CO2_loading_GC"]
        else:
            CO2_loading = (row["CO2_loading_BaCO3"] + row["CO2_loading_GC"])/2

        n_CO2 = n_MEA * CO2_loading
        n_tot = n_MEA + n_H2O + n_CO2

        m.jou_loading[i].flow_mol.fix(100)
        m.jou_loading[i].mole_frac_comp["H2O"].fix(n_H2O/n_tot)
        m.jou_loading[i].mole_frac_comp["MEA"].fix(n_MEA/n_tot)
        m.jou_loading[i].mole_frac_comp["CO2"].fix(n_CO2/n_tot)
        # It isn't entirely clear whether the pressure recorded is the total pressure or the total pressure
        # minus a MEA partial pressure predicted via Raoult's law. I'm taking it as the former
        m.jou_loading[i].pressure.fix(row["total_pressure"]*1e3) # Pressure in kPa 
        m.jou_loading[i].temperature.fix(row["temperature"] + 273.15) # Temperature in C
        obj_expr +=  (
            # CO2 partial pressure
            loss(
                m.jou_loading[i].log_fug_phase_comp["Liq","CO2"]
                - pyo.log(row["CO2_pressure"]*1e3)  # Pressure in kPa
            )
        )

    #################################################################################
    ##### Enthalpy of absorption
    #################################################################################

    # df_liu_dH_abs = pd.read_csv(os.sep.join([data, "liu_dH_abs.csv"]), index_col=None)
    # n_data_liu_dH_abs = len(df_liu_dH_abs)
    # # Parse dataframe into dictionaries
    # dict_liu_dH_abs = {}
    # n_test = 0
    # for T in df_liu_dH_abs["temperature"].unique():
    #     dict_liu_dH_abs[T] = {}
    #     # Relies on the fact that there are the same number of tests for each T
    #     for test in df_liu_dH_abs["test"].unique():
    #         n_test += 1
    #         dict_liu_dH_abs[T][test] = {
    #             "loading": np.array(
    #                 df_liu_dH_abs.loc[
    #                     (df_liu_dH_abs["temperature"] == T) & (df_liu_dH_abs["test"] == test)]["loading"]
    #                 ),
    #             "dH_abs": np.array(
    #                 df_liu_dH_abs.loc[
    #                     (df_liu_dH_abs["temperature"] == T) & (df_liu_dH_abs["test"] == test)]["dH_abs"]
    #                 ),
    #         }
        

    # m.liu_dH_abs = m.params.build_state_block(range(n_data_liu_dH_abs + n_test), defined_state=True)
    # idx_start = 0
    # for T in dict_liu_dH_abs.keys():
    #     for test, subdict in dict_liu_dH_abs[T].items():
    #         subdict["range"] = range(idx_start, idx_start + len(subdict["loading"] + 1))
    #         # Start out with a completely unloaded mixture
    #         w_PZ = 0.069
    #         w_H2O = 1-w_PZ
    #         x_PZ = pyo.value(w_PZ/(w_PZ+(m.params.PZ.mw/m.params.H2O.mw)*w_H2O))
    #         n_PZ = 100 * x_PZ
    #         n_H2O = 100 * (1-x_PZ)
    #         n_CO2 = 0.01 # Add a tiny amount of CO2 to avoid numerical issues
    #         n_tot = n_PZ + n_H2O + n_CO2
    #         m.liu_dH_abs[idx_start].flow_mol.fix(n_tot)
    #         m.liu_dH_abs[idx_start].mole_frac_comp["H2O"].fix(n_H2O/n_tot)
    #         m.liu_dH_abs[idx_start].mole_frac_comp["PZ"].fix(n_PZ/n_tot)
    #         m.liu_dH_abs[idx_start].mole_frac_comp["CO2"].fix(n_CO2/n_tot)
    #         # Not a lot of information about pressure
    #         m.liu_dH_abs[idx_start].pressure.fix(2e5) 
    #         m.liu_dH_abs[idx_start].temperature.fix(T + 273.15)

    #         for k, loading in enumerate(subdict["loading"]):
    #             n_CO2 = n_PZ*loading
    #             n_tot = n_PZ + n_H2O + n_CO2
    #             m.liu_dH_abs[idx_start+k+1]
    #             m.liu_dH_abs[idx_start+k+1].flow_mol.fix(n_tot)
    #             m.liu_dH_abs[idx_start+k+1].mole_frac_comp["H2O"].fix(n_H2O/n_tot)
    #             m.liu_dH_abs[idx_start+k+1].mole_frac_comp["PZ"].fix(n_PZ/n_tot)
    #             m.liu_dH_abs[idx_start+k+1].mole_frac_comp["CO2"].fix(n_CO2/n_tot)
    #             # Not a lot of information about pressure
    #             m.liu_dH_abs[idx_start+k+1].pressure.fix(2e5) 
    #             m.liu_dH_abs[idx_start+k+1].temperature.fix(T + 273.15)

    # df_hartono_dH_abs = pd.read_csv(os.sep.join([data, "hartono_table_15.csv"]), index_col=None)
    # n_data_hartono_dH_abs = len(df_hartono_dH_abs)
    # # Parse dataframe into dictionaries
    # dict_hartono_dH_abs = {}
    # n_test = 0
    # for T in df_hartono_dH_abs["temperature"].unique():
    #     n_test += 1
    #     dict_hartono_dH_abs[T] = {
    #         "loading": np.array(
    #             df_hartono_dH_abs.loc[df_hartono_dH_abs["temperature"] == T]["loading"]
    #             ),
    #         "loading_uncertainty": np.array(
    #             df_hartono_dH_abs.loc[df_hartono_dH_abs["temperature"] == T]["loading_uncertainty"]
    #             ),
    #         "dH_abs": np.array(
    #             df_hartono_dH_abs.loc[df_hartono_dH_abs["temperature"] == T]["dH_abs"]
    #             ),
    #         "dH_abs_uncertainty": np.array(
    #             df_hartono_dH_abs.loc[df_hartono_dH_abs["temperature"] == T]["dH_uncertainty"]
    #             ),
    #     }      

    # m.hartono_dH_blk = m.params.build_state_block(range(n_data_hartono_dH_abs + n_test), defined_state=True)
    # idx_start = 0
    # for T, subdict in dict_hartono_dH_abs.items():
    #     subdict["range"] = range(idx_start, idx_start + len(subdict["loading"] + 1))
    #     # Start out with a completely unloaded mixture
    #     molality_PZ = 2.1
    #     n_H2O = pyo.value(1/m.params.H2O.mw)
    #     n_CO2 = 0.01
    #     n_tot = molality_PZ + n_H2O + n_CO2
    #     blk = m.hartono_dH_blk[idx_start]
    #     blk.flow_mol.fix(n_tot)
    #     blk.mole_frac_comp["H2O"].fix(n_H2O/n_tot)
    #     blk.mole_frac_comp["PZ"].fix(molality_PZ/n_tot)
    #     blk.mole_frac_comp["CO2"].fix(n_CO2/n_tot)
    #     # Not a lot of information about pressure
    #     blk.pressure.fix(2e5) 
    #     blk.temperature.fix(T)
    #     enth_abs_list = []
    #     for k, loading in enumerate(subdict["loading"]):
    #         blk_old = blk
    #         blk = m.hartono_dH_blk[idx_start+k+1]
    #         n_CO2 = molality_PZ*loading
    #         n_tot = n_PZ + n_H2O + n_CO2
    #         blk
    #         blk.flow_mol.fix(n_tot)
    #         blk.mole_frac_comp["H2O"].fix(n_H2O/n_tot)
    #         blk.mole_frac_comp["PZ"].fix(molality_PZ/n_tot)
    #         blk.mole_frac_comp["CO2"].fix(n_CO2/n_tot)
    #         # Not a lot of information about pressure
    #         blk.pressure.fix(2e5) 
    #         blk.temperature.fix(T)
    #         CO2_obj = m.params.CO2
    #         dH_abs_expr = -( # Negative sign bc values reported are positive
    #             blk.energy_internal_mol_phase["Liq"] * blk.flow_mol
    #             - blk_old.energy_internal_mol_phase["Liq"] * blk_old.flow_mol
    #             -NIST.enth_mol_ig_comp.return_expression(
    #                 blk,
    #                 CO2_obj, 
    #                 blk.temperature
    #             ) * (blk.flow_mol - blk_old.flow_mol)
    #             # There should be an additional term for PV work to pressurize
    #             # the vapor phase with CO2, but we can't account for it because
    #             # we have limited info about the vapor phase and the pressurized
    #             # cylinder the CO2 was being discharged from
    #         ) / (blk.flow_mol - blk_old.flow_mol)
    #         enth_abs_list.append(dH_abs_expr)
    #         # if loading <= 0.75: # threshold bc of missing vapor phase enthalpy
    #         #     obj_expr += loss(
    #         #         (
    #         #             dH_abs_expr*1e-3 # data is in kJ
    #         #             -subdict["dH_abs"][k]
    #         #         ) / subdict["dH_abs_uncertainty"][k]
    #         #     )
    #     subdict["enth_abs_expr"] = enth_abs_list
    #     idx_start += len(subdict["loading"]) + 1

                

    iscale.calculate_scaling_factors(m)
    print(f"DOF: {degrees_of_freedom(m)}")
    # initialize_inherent_reactions(m.hillard_loading)
    # m.hillard_loading.initialize(
    #     hold_state=False,
    #     outlvl=init_outlevel,
    #     optarg=optarg
    # )

    initialize_inherent_reactions(m.jou_loading)
    m.jou_loading.initialize(
        hold_state=False,
        outlvl=init_outlevel,
        optarg=optarg
    )
    # initialize_inherent_reactions(m.hartono_dH_blk)
    # m.hartono_dH_blk.initialize(
    #     hold_state=False,
    #     outlvl=init_outlevel,
    #     optarg=optarg
    # )
    # Need to create heat capacity equations after initialization because
    # initialization method doesn't know about them
    # create_heat_capacity_no_inherent_rxns(m.state_chen_54_59)
    # assert degrees_of_freedom(m.state_chen_54_59) == 0
    # for i, row in df_chen_54_59.iterrows():
    #     obj_expr +=  (
    #         # H2O partial pressure
    #         loss(
    #             0.1*(m.state_chen_54_59[i].Liq_cp - float(row["cp"]))
    #         )
    #     )

    iscale.calculate_scaling_factors(m)

    estimated_vars = get_estimated_params(m)
    # Apparently the model reduces to the Margules model at alpha=0,
    # which corresponds to random mixing. Alpha values less than zero
    # are not thermodynamically meaningful, because you can't get more
    # random than fully random.
    # m.params.Liq.alpha["H2O", "PZ"].lb = 0
    # m.params.Liq.alpha["H2O", "PZ"].ub = 1

    for var in estimated_vars:
        obj_expr += 0.01 * loss(
            (var - var.value) * iscale.get_scaling_factor(var)
        )
    print("Test1")
    m.obj = pyo.Objective(expr=obj_expr)
    print("Test2")
    m_scaled = pyo.TransformationFactory('core.scale_model').create_using(m, rename=False)
    print("Test3")
    optarg.pop("nlp_scaling_method", None) # Scaled model doesn't need user scaling
    optarg["max_iter"] = 500
    solver = get_solver(
        "ipopt",
        options=optarg
    )
    # solve_indexed_blocks(solver, m.state_chen_54_59, tee=True)
    estimated_vars_scaled = get_estimated_params(m_scaled)
    for var in estimated_vars_scaled:
        var.unfix()
    print("Test4")
    res = solver.solve(m_scaled, tee=True)
    print("Test5")
    pyo.assert_optimal_termination(res)

    


    inv_red_hess = inv_reduced_hessian_barrier(m_scaled, estimated_vars_scaled, tee=True, solver_options=optarg)
    W, V = eigh(inv_red_hess[1])

    print(f"Largest eigenvalue of inverse reduced Hessian: {W[-1]}")
    print("\n" + "Variables with most uncertainty:")
    for i in np.where(abs(V[:, -1]) > 0.3)[0]:
        print(str(i) + ": " + estimated_vars_scaled[i].name)

    pyo.TransformationFactory('core.scale_model').propagate_solution(m_scaled,m)
    print("==================================================")
    print("========== Variables with uncertainty ============")
    print("==================================================")
    def gsf(var):
        return iscale.get_scaling_factor(var, default=1)
    for i, var in enumerate(estimated_vars):
        print(f"{var.name}: {var.value:.3e} +/- {pyo.sqrt(inv_red_hess[1][i][i])/gsf(var):.3e}")

    colors = [
        "firebrick",
        "royalblue",
        "forestgreen",
        "goldenrod",
        "magenta", 
        "orangered",
        "crimson",
        "darkcyan"
    ]
    # print ("\n Data source: Hillard")
    # P_CO2_rel_err = 0
    # n_data = 0
    # for molality in df_hillard_loading["PZ_molality"].unique():
    #     fig = plt.figure()
    #     ax = fig.subplots()
    #     df_hillard_molality =df_hillard_loading.loc[df_hillard_loading["PZ_molality"] ==  molality]
    #     for j, T in enumerate(df_hillard_molality["temperature"].unique()):
    #         x = []
    #         y_data =[]
    #         y_model = []
    #         for i, row in df_hillard_molality.loc[df_hillard_molality["temperature"] ==  T].iterrows():
    #             n_data += 1
    #             loading = float(row["CO2_loading"])
    #             x.append(loading)
    #             y_data.append(row["CO2_partial_pressure"])
    #             candidate_idx = np.where(df_hillard_loading["CO2_partial_pressure"].values == row["CO2_partial_pressure"])
    #             assert len(candidate_idx) == 1
    #             assert len(candidate_idx[0]) == 1
    #             y_model.append(
    #                 pyo.value(m.hillard_loading[candidate_idx[0][0]].fug_phase_comp["Liq","CO2"] / 1e3)  # Convert from Pa to kPa
    #             )

    #         for y_data_pt, y_model_pt in zip(y_data, y_model):
    #             P_CO2_rel_err += abs((y_data_pt-y_model_pt)/y_data_pt) 

    #         ax.semilogy(x, y_data, linestyle="none", marker="o", color=colors[j])
    #         ax.semilogy(x, y_model, linestyle="-", marker="none", color=colors[j], label=f"T = {T}" )

    #     ax.set_xlabel("CO2 Loading (Rochelle convention)")
    #     ax.set_ylabel("Total vapor pressure (kPa)")
    #     ax.set_title(f"Hillard Loading PZ molality = {molality}")
    #     ax.legend()
    #     fig.show()
    # print(f"Total number of data points: {n_data}")
    # print(f"Average relative error in Co2 pressure: {P_CO2_rel_err/n_data}")

    print ("\n Data source: Jou")
    P_CO2_rel_err = 0
    n_data = 0


    fig = plt.figure()
    ax = fig.subplots()

    for j, T in enumerate(df_jou_loading["temperature"].unique()):
        x = []
        y_data =[]
        y_model = []
        for i, row in df_jou_loading.loc[df_jou_loading["temperature"] ==  T].iterrows():
            n_data +=1
            if np.isnan(row["CO2_loading_GC"]):
                CO2_loading = row["CO2_loading_BaCO3"]
            elif np.isnan(row["CO2_loading_BaCO3"]):
                CO2_loading = row["CO2_loading_GC"]
            else:
                CO2_loading = (row["CO2_loading_BaCO3"] + row["CO2_loading_GC"])/2
            x.append(CO2_loading)
            y_data.append(row["CO2_pressure"])
            candidate_idx = np.where(df_jou_loading["CO2_pressure"].values == row["CO2_pressure"])
            assert len(candidate_idx) == 1
            assert len(candidate_idx[0]) == 1
            y_model.append(
                pyo.value(m.jou_loading[candidate_idx[0][0]].fug_phase_comp["Liq","CO2"])/1e3  # Convert from Pa to kPa
            )
        
        for y_data_pt, y_model_pt in zip(y_data, y_model):
            P_CO2_rel_err += abs((y_data_pt-y_model_pt)/y_data_pt) 

        ax.semilogy(x, y_data, linestyle="none", marker="o", color=colors[j])
        ax.semilogy(x, y_model, linestyle="-", marker="none", color=colors[j], label=f"T = {T}" )
    ax.set_xlabel("CO2 Loading")
    ax.set_ylabel("Total vapor pressure (Pa)")
    ax.set_title(f"Jou Loading")
    ax.legend()
    fig.show()

    print(f"Total number of data points: {n_data}")
    print(f"Average relative error in CO2 pressure: {P_CO2_rel_err/n_data}")

    # fig = plt.figure()
    # ax = fig.subplots()

    # for i, T in enumerate(dict_hartono_dH_abs.keys()):
    #     subdict = dict_hartono_dH_abs[T]
    #     x = []
    #     y_data =[]
    #     y_model = []
    #     for k, loading in enumerate(subdict["loading"]):
    #         x.append(loading)
    #         y_data.append(subdict["dH_abs"][k])
    #         y_model.append(1e-3*pyo.value(subdict["enth_abs_expr"][k]))
                
    #     ax.plot(x, y_data, linestyle="none", marker="o", color=colors[i])
    #     ax.plot(x, y_model, linestyle="-", marker="none", color=colors[i], label=f"T = {T}" )
    # ax.set_xlabel("CO2 Loading")
    # ax.set_ylabel("Enth_Abs (kJ/mol)")
    # ax.set_title("Hartono et al. Enthalpy of Absorption")
    # ax.legend()
    # fig.show()
    print("ok, boomer")