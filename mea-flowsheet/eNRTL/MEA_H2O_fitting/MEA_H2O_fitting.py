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
from MEA_eNRTL import get_prop_dict, create_heat_capacity_no_inherent_rxns

def loss(x):
    return 0.5*x ** 2

def get_estimated_params(m, henry, cp):
    param_list = [
        m.params.Liq.tau_A["H2O", "MEA"],
        m.params.Liq.tau_B["H2O", "MEA"],
        m.params.Liq.tau_A["MEA", "H2O"],
        m.params.Liq.tau_B["MEA", "H2O"],
        m.params.Liq.alpha["H2O", "MEA"],
    ]
    if henry:
        param_list += [
            m.params.PZ.henry_coeff_1,
            # m.params.PZ.henry_coeff_2,
            m.params.PZ.henry_coeff_3,
            m.params.PZ.henry_coeff_4
        ]
    if cp:
        param_list += [
            m.params.MEA.cp_mol_liq_comp_coeff_A,
            m.params.MEA.cp_mol_liq_comp_coeff_B
        ]
    return param_list

def create_and_scale_params(m, henry, cp):
    config = get_prop_dict(["H2O", "MEA"], excluded_rxns=["H2O_autoionization", "MEA_protonation"])

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
    }
    for comp, sf_x in mole_frac_scaling_factors.items():
        params.set_default_scaling("mole_frac_comp", sf_x, index=comp)
        params.set_default_scaling("mole_frac_phase_comp", sf_x, index=("Liq", comp))
        params.set_default_scaling(
            "flow_mol_phase_comp",
            sf_x * scaling_factor_flow_mol,
            index=("Liq", comp)
        )
        params.set_default_scaling("mole_frac_phase_comp_true", sf_x, index=("Liq", comp))
        params.set_default_scaling(
            "flow_mol_phase_comp_true",
            sf_x * scaling_factor_flow_mol,
            index=("Liq", comp)
        )

    # for comp, sf_x in mole_frac_true_scaling_factors.items():
    #     params.set_default_scaling("mole_frac_phase_comp_true", sf_x, index=("Liq", comp))
    #     params.set_default_scaling(
    #         "flow_mol_phase_comp_true",
    #         sf_x * scaling_factor_flow_mol,
    #         index=("Liq", comp)
    #     )

    # params.set_default_scaling("apparent_inherent_reaction_extent", scaling_factor_flow_mol*1e4, index="H2O_autoionization")
    # params.set_default_scaling("apparent_inherent_reaction_extent", scaling_factor_flow_mol*1e4, index="PZ_protonation")

    iscale.set_scaling_factor(m.params.Liq.alpha, 5) 
    iscale.set_scaling_factor(m.params.Liq.tau_A, 1) # Reminder that it's well-scaled
    iscale.set_scaling_factor(m.params.Liq.tau_B, 1/300)

    if henry:
        iscale.set_scaling_factor(m.params.MEA.henry_coeff_1, 1/300)
        iscale.set_scaling_factor(m.params.MEA.henry_coeff_2, 1) # Reminder that it's well-scaled
        iscale.set_scaling_factor(m.params.MEA.henry_coeff_3, 300)
        iscale.set_scaling_factor(m.params.MEA.henry_coeff_4, 1) # Reminder that it's well-scaled

    if cp:
        iscale.set_scaling_factor(m.params.MEA.cp_mol_liq_comp_coeff_A, 3e-2)
        iscale.set_scaling_factor(m.params.MEA.cp_mol_liq_comp_coeff_B, 3e-2)

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
    init_outlevel = idaeslog.INFO
    henry=False
    include_page_data = False
    fit_pure_MEA_cp = False

    optimize = False

    data = os.sep.join([this_file_dir(), "data"])
    m = pyo.ConcreteModel()
    create_and_scale_params(m, henry=henry, cp=fit_pure_MEA_cp)

    # Initial guess taken from Zhang et al., 2011
    m.params.Liq.tau_A["H2O", "MEA"].set_value(0.1559)
    m.params.Liq.tau_B["H2O", "MEA"].set_value(110.80)
    m.params.Liq.tau_A["MEA", "H2O"].set_value(1.5201)
    m.params.Liq.tau_B["MEA", "H2O"].set_value(-910.30)
    m.params.Liq.alpha["H2O", "MEA"].set_value(0.2)

    obj_expr = 0

    df_nath = pd.read_csv(os.sep.join([data, "nath_bender_1983_T_x_P.csv"]), index_col=None)

    n_data_nath = len(df_nath["temperature"])
    m.state_nath = m.params.build_state_block(range(n_data_nath), defined_state=True)

    for i, row in df_nath.iterrows():
        m.state_nath[i].flow_mol.fix(100)
        m.state_nath[i].mole_frac_comp["H2O"].fix(row["mole_frac_H2O"])
        m.state_nath[i].mole_frac_comp["MEA"].fix(1-row["mole_frac_H2O"])
        m.state_nath[i].pressure.fix(row["total_pressure"]*133.3) # Pressure in Torr
        m.state_nath[i].temperature.fix(row["temperature"] + 273.15) # Temperature in C
        obj_expr +=  (
            # Total pressure
            loss(
                pyo.log((m.state_nath[i].fug_phase_comp["Liq","MEA"] + m.state_nath[i].fug_phase_comp["Liq","H2O"])/pyo.units.Pa) 
                - pyo.log(row["total_pressure"]*133.3)
            )
        )

    df_cai = pd.read_csv(os.sep.join([data, "cai_xie_wu_1996_T_x_y.csv"]), index_col=None)

    n_data_cai = len(df_cai["temperature"])
    m.state_cai = m.params.build_state_block(range(n_data_cai), defined_state=True)

    for i, row in df_cai.iterrows():
        m.state_cai[i].flow_mol.fix(100)
        m.state_cai[i].mole_frac_comp["H2O"].fix(row["mole_frac_liq_H2O"])
        m.state_cai[i].mole_frac_comp["MEA"].fix(1-row["mole_frac_liq_H2O"])
        m.state_cai[i].pressure.fix(row["pressure"]*1e3) # Pressure in kPa
        # Note: boiling temperature is actually dependent variable, but fix it here for initialization
        m.state_cai[i].temperature.fix(row["temperature"]) # Temperature in K
        obj_expr +=  (
            # Vapor phase water mole fraction
            loss(
                3*(m.state_cai[i].fug_phase_comp["Liq","H2O"] / (row["pressure"] *1e3) - row["mole_frac_vap_H2O"])
            )
            # Boiling point
            + loss(
                1e-1*(m.state_cai[i].temperature - row["temperature"])
            )
        )

    df_touhara_excess_enthalpy = pd.read_csv(os.sep.join([data, "touhara_et_al_1982_excess_enthalpy.csv"]), index_col=None)

    n_data_touhara_excess_enthalpy = len(df_touhara_excess_enthalpy["mole_frac_MEA"])
    m.state_touhara_excess_enthalpy = m.params.build_state_block(range(n_data_touhara_excess_enthalpy), defined_state=True)

    for i, row in df_touhara_excess_enthalpy.iterrows():
        m.state_touhara_excess_enthalpy[i].flow_mol.fix(100)
        m.state_touhara_excess_enthalpy[i].mole_frac_comp["H2O"].fix(1-row["mole_frac_MEA"])
        m.state_touhara_excess_enthalpy[i].mole_frac_comp["MEA"].fix(row["mole_frac_MEA"])
        # Unclear what pressure at which excess enthalpy measurements occured, assumed to be 1 atm
        m.state_touhara_excess_enthalpy[i].pressure.fix(1.013e5)
        m.state_touhara_excess_enthalpy[i].temperature.fix(298.15)
        obj_expr +=  (
            loss(
                5e-4*(m.state_touhara_excess_enthalpy[i].enth_mol_phase_excess["Liq"] + row["negative_excess_enthalpy"])
            )
        )

    df_posey_excess_enthalpy = pd.read_csv(os.sep.join([data, "posey_1996_excess_enthalpy.csv"]), index_col=None)

    n_data_posey_excess_enthalpy = len(df_posey_excess_enthalpy["mole_frac_MEA"])
    m.state_posey_excess_enthalpy = m.params.build_state_block(range(n_data_posey_excess_enthalpy), defined_state=True)

    for i, row in df_posey_excess_enthalpy.iterrows():
        m.state_posey_excess_enthalpy[i].flow_mol.fix(100)
        m.state_posey_excess_enthalpy[i].mole_frac_comp["H2O"].fix(1-row["mole_frac_MEA"])
        m.state_posey_excess_enthalpy[i].mole_frac_comp["MEA"].fix(row["mole_frac_MEA"])
        # Unclear what pressure at which excess enthalpy measurements occured, assumed to be 1 atm
        m.state_posey_excess_enthalpy[i].pressure.fix(1.013e5)
        m.state_posey_excess_enthalpy[i].temperature.fix(row["temperature"] + 273.15) # Temperature in C
        obj_expr +=  (
            # Scale contribution bc of so few data points
            5*loss(
                5e-4*(m.state_posey_excess_enthalpy[i].enth_mol_phase_excess["Liq"] - row["excess_enthalpy"])
            )
        )

    # Excluding heat capacity data because of disagreement between sources
    if include_page_data:
        df_page_et_al  = pd.read_csv(os.sep.join([data, "page_huot_jolicoeur_1992_cp_density.csv"]), index_col=None)
        n_data_chen_54_59 = len(df_page_et_al["temperature"])
        m.state_page_et_al = m.params.build_state_block(range(n_data_chen_54_59), defined_state=True)

        page_interp_bounds = {key: [float("NaN"), float("NaN")] for key in df_page_et_al["temperature"].unique()}
        for i, row in df_page_et_al.iterrows():
            x_MEA = float(row["mole_frac_MEA"])
            if x_MEA < 1e-4:
                x_MEA = 1e-4
                page_interp_bounds[row["temperature"]][0] = row["heat_capacity_mass"]
            elif x_MEA > 1 - 1e-4:
                x_MEA = 1 - 1e-4
                page_interp_bounds[row["temperature"]][1] = row["heat_capacity_mass"]
            x_H2O = 1 - x_MEA
            m.state_page_et_al[i].flow_mol.fix(100)
            m.state_page_et_al[i].mole_frac_comp["H2O"].fix(x_H2O)
            m.state_page_et_al[i].mole_frac_comp["MEA"].fix(x_MEA)
            m.state_page_et_al[i].pressure.fix(1.013e5) # No pressure listed in experimental paper
            m.state_page_et_al[i].temperature.fix(float(row["temperature"])+273.15) # Temperature in C

    print(f"DOF: {degrees_of_freedom(m)}")

    m.state_nath.initialize(
        hold_state=False,
        outlvl=init_outlevel,
        optarg=optarg
    )
    m.state_cai.initialize(
        hold_state=False,
        outlvl=init_outlevel,
        optarg=optarg
    )
    m.state_touhara_excess_enthalpy.initialize(
        hold_state=False,
        outlvl=init_outlevel,
        optarg=optarg
    )
    m.state_posey_excess_enthalpy.initialize(
        hold_state=False,
        outlvl=init_outlevel,
        optarg=optarg
    )
    if include_page_data:
        m.state_page_et_al.initialize(
            hold_state=False,
            outlvl=init_outlevel,
            optarg=optarg
        )

    # Create VLE constraints for Posey T-x-y data
    for i, row in df_cai.iterrows():
        m.state_cai[i].temperature.unfix()
        @m.state_cai[i].Constraint()
        def vapor_pressure_eqn(b):
            return (
                b.fug_phase_comp["Liq","H2O"] + b.fug_phase_comp["Liq","MEA"]
                == row["pressure"]*1e3
            ) 
        iscale.constraint_scaling_transform(m.state_cai[i].vapor_pressure_eqn, 1e-5)
    # Need to create heat capacity equations after initialization because
    # initialization method doesn't know about them
    if include_page_data:
        create_heat_capacity_no_inherent_rxns(m.state_page_et_al)
        assert degrees_of_freedom(m.state_page_et_al) == 0

        excess_cp_list = []
        for i, row in df_page_et_al.iterrows():
            if not np.isnan(float(row["heat_capacity_mass"])):
                cp0 = page_interp_bounds[row["temperature"]][0] * 18.02
                cp1 = page_interp_bounds[row["temperature"]][1] * 61.08
                x_MEA = float(row["mole_frac_MEA"])
                excess_cp = float(row["heat_capacity_mass"])*((1-x_MEA)*18.02 + x_MEA*61.08) - float((1-x_MEA)*cp0 + x_MEA*cp1)
                excess_cp_list.append(excess_cp)
                obj_expr +=  (
                    # H2O partial pressure
                    loss(
                        0.1*(m.state_page_et_al[i].Liq_cp_phase_excess - excess_cp)
                    )
                )

    iscale.calculate_scaling_factors(m)

    estimated_vars = get_estimated_params(m, henry=henry, cp=fit_pure_MEA_cp)
    # Apparently the model reduces to the Margules model at alpha=0,
    # which corresponds to random mixing. Alpha values less than zero
    # are not thermodynamically meaningful, because you can't get more
    # random than fully random.
    m.params.Liq.alpha["H2O", "MEA"].lb = 0
    # m.params.Liq.alpha["H2O", "MEA"].ub = 1

    for var in estimated_vars:
        obj_expr += 1e-4 * loss(
            (var - var.value) * iscale.get_scaling_factor(var)
        )

    m.obj = pyo.Objective(expr=obj_expr)
    m_scaled = pyo.TransformationFactory('core.scale_model').create_using(m, rename=False)
    optarg.pop("nlp_scaling_method", None) # Scaled model doesn't need user scaling
    optarg["max_iter"] = 500
    solver = get_solver(
        "ipopt",
        options=optarg
    )
    # solve_indexed_blocks(solver, m.state_page_et_al, tee=True)
    estimated_vars_scaled = get_estimated_params(m_scaled, henry=henry, cp=fit_pure_MEA_cp)
    if optimize:
        for var in estimated_vars_scaled:
            var.unfix()
        res = solver.solve(m_scaled, tee=True)
        pyo.assert_optimal_termination(res)

    
        optarg["bound_push"] = 1e-6
        # optarg["mu_init"] = 1e-4

        # try:
        #     inv_red_hess = inv_reduced_hessian_barrier(m_scaled, estimated_vars_scaled, tee=True, solver_options=optarg)
        #     W, V = eigh(inv_red_hess[1])

        #     print(f"Largest eigenvalue of inverse reduced Hessian: {W[-1]}")
        #     print("\n" + "Variables with most uncertainty:")
        #     for i in np.where(abs(V[:, -1]) > 0.3)[0]:
        #         print(str(i) + ": " + estimated_vars_scaled[i].name)
        #     uncertainty = True
        # except (ValueError, RuntimeError):
        #     print("Variable close to bound or Hessian singular, cannot perform Hessian analysis")
        #     uncertainty = False
        uncertainty=False
    else:
        uncertainty = False
        res = solver.solve(m_scaled, tee=True)
        pyo.assert_optimal_termination(res)

    pyo.TransformationFactory('core.scale_model').propagate_solution(m_scaled,m)
    print("==================================================")
    print("========== Variables with uncertainty ============")
    print("==================================================")
    def gsf(var):
        return iscale.get_scaling_factor(var, default=1)
    for i, var in enumerate(estimated_vars):
        if uncertainty:
            print(f"{var.name}: {var.value:.3e} +/- {pyo.sqrt(inv_red_hess[1][i][i])/gsf(var):.3e}")
        else:
            print(f"{var.name}: {var.value:.3e}")
        

    fig = plt.figure()
    ax = fig.subplots()
    colors1 = ["firebrick", "royalblue", "forestgreen"]
    colors2 = ["gold", "mediumturquoise", "darkorchid"]
    print ("\n Data source: Nath")
    P_tot_rel_err = 0
    n_data = 0
    for T, color1, color2 in zip([60, 78, 91.7], colors1, colors2):
        x = []
        P_tot_data = []
        P_tot_model = []
        for i, row in df_nath.loc[df_nath["temperature"] ==  T].iterrows():
            n_data+=1
            x_H2O = float(row["mole_frac_H2O"])
            x.append(x_H2O)
            P_tot_data.append(row["total_pressure"])
            P_tot_model.append(
                pyo.value(m.state_nath[i].fug_phase_comp["Liq","MEA"] 
                 + m.state_nath[i].fug_phase_comp["Liq","H2O"]) / 133.3 # Convert from Pa to Torr
                )

        for y_data, y_model in zip(P_tot_data,P_tot_model):
            P_tot_rel_err += abs((y_data-y_model)/y_data)

        ax.plot(x, P_tot_data, linestyle="none", marker="o", color=color1)
        ax.plot(x, P_tot_model, linestyle="-", marker="none", color=color1, label=f"T = {T}" )

   
    print(f"Total number of data points: {n_data}")
    print(f"Average relative error in total pressure: {P_tot_rel_err/n_data}")
    ax.set_xlabel("Water mole fraction")
    ax.set_ylabel("Total vapor pressure (Torr)")
    ax.set_title("Nath MEA-H2O Volatility")
    ax.legend()

    fig.show()

    fig = plt.figure()
    ax = fig.subplots()
    colors1 = ["firebrick", "royalblue", "forestgreen"]
    colors2 = ["gold", "mediumturquoise", "darkorchid"]
    print ("\n Data source: Cai Xie Wu T-x-y")
    T_rel_err = 0
    y_H2O_rel_err = 0
    n_data = 0
    for P, color1 in zip([101.33, 66.66], colors1):
        x = []
        T_data = []
        T_model = []
        x_Vap_MEA_data = []
        x_Vap_MEA_model = []
        for i, row in df_cai.loc[df_cai["pressure"] ==  P].iterrows():
            n_data+=1
            x_MEA = float(1-row["mole_frac_liq_H2O"])
            x.append(x_MEA)
            T_data.append(row["temperature"])
            T_model.append(pyo.value(m.state_cai[i].temperature))
            x_Vap_MEA_data.append(float(1-row["mole_frac_vap_H2O"]))
            x_Vap_MEA_model.append(
                pyo.value(
                    m.state_cai[i].fug_phase_comp["Liq","MEA"]/(P*1e3)
                )
            )

        for y_data, y_model in zip(T_data, T_model):
            T_rel_err += abs((y_data-y_model)/y_data)
        for y_data, y_model in zip(x_Vap_MEA_data, x_Vap_MEA_model):
            if y_data > 1e-5:
                y_H2O_rel_err += abs((y_data-y_model)/y_data)

        ax.plot(x, T_data, linestyle="none", marker="o", color=color1)
        ax.plot(x, T_model, linestyle="-", marker="none", color=color1, label=f"Liquid P = {P} kPa" )

        ax.plot(x_Vap_MEA_data, T_data, linestyle="none", marker="x", color=color1)
        ax.plot(x_Vap_MEA_model, T_model, linestyle="--", marker="none", color=color1, label=f"Vapor P = {P} kPa" )

   
    print(f"Total number of data points: {n_data}")
    print(f"Average relative error in temperature: {T_rel_err/n_data}")
    print(f"Average relative error in vapor water mole fraction: {y_H2O_rel_err/n_data}")
    ax.set_xlabel("MEA mole fraction")
    ax.set_ylabel("Temperature")
    ax.set_title("Cai MEA-H2O T-x-y")
    ax.legend()

    ax.set_xlim([0, 1])
    ax.set_ylim([350, 450])

    fig.show()


    fig = plt.figure()
    ax = fig.subplots()
    print ("\n Data source: Touhara enthalpy")
    enth_rel_err = 0
    n_data = 0
    x = []
    excess_enthalpy_data = []
    excess_enthalpy_model = []
    for i, row in df_touhara_excess_enthalpy.iterrows():
        n_data+=1
        x_MEA = float(row["mole_frac_MEA"])
        x.append(x_MEA)
        excess_enthalpy_data.append(-row["negative_excess_enthalpy"])
        excess_enthalpy_model.append(
            pyo.value(
                m.state_touhara_excess_enthalpy[i].enth_mol_phase_excess["Liq"]
            )
        )

    for y_data, y_model in zip(excess_enthalpy_data ,excess_enthalpy_model):
        enth_rel_err += abs((y_data-y_model)/y_data)

    ax.plot(x, excess_enthalpy_data, linestyle="none", marker="o", color=colors1[0], label="Touhara Data")
    ax.plot(x, excess_enthalpy_model, linestyle="-", marker="none", color=colors1[0], label="T = 25 C" )
   
    print(f"Total number of data points: {n_data}")
    print(f"Average relative error in excess enthalpy: {enth_rel_err/n_data}")

    print ("\n Data source: Posey enthalpy")
    enth_rel_err = 0
    n_data = 0
    x = []
    excess_enthalpy_data = []
    excess_enthalpy_model = []
    for i, row in df_posey_excess_enthalpy.iterrows():
        n_data+=1
        x_MEA = float(row["mole_frac_MEA"])
        x.append(x_MEA)
        excess_enthalpy_data.append(row["excess_enthalpy"])
        excess_enthalpy_model.append(
            pyo.value(
                m.state_posey_excess_enthalpy[i].enth_mol_phase_excess["Liq"]
            )
        )

    for y_data, y_model in zip(excess_enthalpy_data ,excess_enthalpy_model):
        enth_rel_err += abs((y_data-y_model)/y_data)

    ax.plot(x, excess_enthalpy_data, linestyle="none", marker="x", color=colors1[1], label="Posey Data")
    ax.plot(x, excess_enthalpy_model, linestyle="-", marker="none", color=colors1[1], label="T = 70 C" )
   
    print(f"Total number of data points: {n_data}")
    print(f"Average relative error in excess enthalpy: {enth_rel_err/n_data}")

    ax.set_xlabel("MEA mole fraction")
    ax.set_ylabel("Excess enthalpy (J/mol)")
    ax.set_title("MEA-H2O Excess Enthalpy")
    ax.legend()
    ax.set_ylim([0, -3000])
    fig.show()

    if include_page_data:
        fig = plt.figure()
        ax = fig.subplots()

        fig2 = plt.figure()
        ax2 = fig2.subplots()

        colors = ["firebrick", "royalblue", "forestgreen", "gold", "darkorchid"]

        print ("\n Data source: Page et al.")
        cp_rel_err = 0
        n_data = 0
        for T, color in zip(df_page_et_al["temperature"].unique(), colors):
            x = []
            y_data =[]
            y_model = []
            y_model_excess = []
            y_interp = []
            for i, row in df_page_et_al.loc[df_page_et_al["temperature"] ==  T].iterrows():
                n_data += 1
                x.append(float(row["mole_frac_MEA"]))
                y_data.append(pyo.value(1e3*float(row["heat_capacity_mass"])*m.state_page_et_al[i].mw))
                y_model.append(pyo.value(m.state_page_et_al[i].Liq_cp))
                y_model_excess.append(pyo.value(m.state_page_et_al[i].Liq_cp_phase_excess))
                # y_ideal.append(pyo.value(m.state_page_et_al[i].Liq_cp-m.state_page_et_al[i].Liq_cp_phase_excess))
            
            for y_data_pt, y_model_pt in zip(y_data, y_model):
                if not np.isnan(y_data_pt):
                    cp_rel_err += abs((y_data_pt-y_model_pt)/y_data_pt) 

            for x_pt in x:
                y_interp.append((1-x_pt)*y_data[0]+x_pt*y_data[-1])

            ax.plot(x, y_data, linestyle="none", marker="o", color=color)
            ax.plot(x, y_model, linestyle="-", marker="none", color=color, label=f"Temperature = {T}" )
            ax.plot(x, y_interp, linestyle="--", marker="none", color=color)

            ax2.plot(x, y_model_excess, linestyle="-", marker="none", color=color, label=f"Temperature = {T}" )
            ax2.plot(x, np.array(y_data)-np.array(y_interp), linestyle="none", marker="o", color=color)


        print(f"Total number of data points: {n_data}")
        print(f"Average relative error in c_p: {cp_rel_err/n_data}")

        ax.set_xlabel("Mole fraction MEA")
        ax.set_ylabel("Heat capacity (J/mol K)")
        ax.set_title("Page et al. Heat Capacity")
        ax.legend()

        fig.show()

        ax2.set_xlabel("Mole fraction MEA")
        ax2.set_ylabel("Excess heat capacity (J/mol K)")
        ax2.set_title("Page et al. Excess Heat Capacity")
        ax2.legend()
        fig2.show()
    print("ok, boomer")