#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for Flash unit model based on the complementarity formulation.
Author: Vibhav Dabadghao
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    value,
    units,
    Var,
    units as pyunits,
    Objective,
    SolverFactory,
    log,
)
from pyomo.opt import SolverStatus, TerminationCondition, ProblemFormat
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.core.expr.visitor import identify_variables
from pyomo.common.collections import ComponentSet

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.flash import Flash, EnergySplittingType
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_property_ncp import (
        GenericParameterBlock)
from idaes.models.properties.modular_properties.examples.ASU_PR \
    import configuration

from idaes.models.properties import iapws95
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
solver.options['linear_solver'] = 'MA57'
solver.options['bound_push'] = 1e-06  # default = 1e-02
solver.options['mu_init'] = 1e-06  # default = 1e-01
solver.options["bound_relax_factor"] = 1e-05  # default = 1e-08
solver.options['tol'] = 1e-05  # default = 1e-08  


# -----------------------------------------------------------------------------
def scale_variable(var, sf):
    if var.is_indexed():
        for _, v in var.items():
            iscale.set_scaling_factor(v, sf)
    else:
        iscale.set_scaling_factor(var, sf)
    return None


def apply_variable_scaling(m):
    b_in = m.fs.flash.control_volume.properties_in[0.0]
    b_out = m.fs.flash.control_volume.properties_out[0.0]
    for b in [b_in, b_out]:
        scale_variable(b.flow_mol, 1.0)
        scale_variable(b_in.flow_mol_phase, 1.0)
        scale_variable(b.mole_frac_phase_comp, 1.0)
        scale_variable(b.mole_frac_comp, 1.0)
        scale_variable(b.s_Vap_Liq, 1.0)
        scale_variable(b.gp_Vap_Liq, 1.0)
        scale_variable(b.gn_Vap_Liq, 1.0)
    return None


def rescale_constraint(con, sf_additional):
    if con.is_indexed():
        for _, c in con.items():
            sf_old = iscale.get_constraint_transform_applied_scaling_factor(c)
            if sf_old is not None:
                sf_new = sf_old * sf_additional
            else:
                sf_new = sf_additional
            iscale.constraint_scaling_transform(c, sf_new, overwrite=True)
    else:
        sf_old = iscale.get_constraint_transform_applied_scaling_factor(con)
        if sf_old is not None:
            sf_new = sf_old * sf_additional
        else:
            sf_new = sf_additional
        iscale.constraint_scaling_transform(con, sf_new, overwrite=True)
    return None


def apply_constraint_scaling(m):
    # Constraints on control_volume block
    rescale_constraint(m.fs.flash.control_volume.material_balances, 1e+02)
    rescale_constraint(m.fs.flash.control_volume.pressure_balance, 1e+05)

    # Constraints on state blocks
    b_in = m.fs.flash.control_volume.properties_in[0.0]
    b_out = m.fs.flash.control_volume.properties_out[0.0]
    for b in [b_in, b_out]:
        rescale_constraint(b.log_mole_frac_phase_comp_eqn, 1e-03)
        rescale_constraint(b.sum_mole_frac, 1e-03)        
        rescale_constraint(b.total_flow_balance, 1e+02)
        rescale_constraint(b.phase_fraction_constraint, 1e+02)
        rescale_constraint(b.tbar_constraint_Vap_Liq, 1e+02)
        if b == b_in:
            rescale_constraint(b.component_flow_balances, 1e-02)
            rescale_constraint(b.equilibrium_constraint, 1e+00)
        if b == b_out:
            rescale_constraint(b.sum_mole_frac_out, 1e-03)

    return None


def test_single_stage_flash():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = GenericParameterBlock(default=configuration)
    m.fs.flash = Flash(default={"property_package": m.fs.properties})
    
    iscale.calculate_scaling_factors(m)
    apply_variable_scaling(m)
    apply_constraint_scaling(m)

    # Specify inputs
    components = ['argon', 'nitrogen', 'oxygen']
    P_atm = 101325.0  # [Pa]
    z = {'argon': 0.05, 'nitrogen': 0.60, 'oxygen': 0.35}
    flow_mol = 1.0
    temperature = 95
    duty = 0

    m.fs.flash.inlet.flow_mol.fix(flow_mol)  # 100
    m.fs.flash.inlet.temperature.fix(temperature)  # 300
    m.fs.flash.deltaP.fix(0)
    for comp in components:
        m.fs.flash.inlet.mole_frac_comp[0, comp].fix(z[comp])
    m.fs.flash.inlet.pressure.fix(P_atm)
    m.fs.flash.heat_duty.fix(duty)

    m.fs.flash.initialize(outlvl=idaeslog.INFO)    
    results = solver.solve(m.fs.flash, tee=True)
    m.fs.flash.report()
    
    duty = -1000
    print(f'''
          ------------------------
          Solving for Q = {duty}
          ------------------------''')
    m.fs.flash.heat_duty.fix(duty)
    m.fs.flash.initialize(outlvl=idaeslog.INFO)
    results = solver.solve(m.fs.flash, tee=True)
    
    Tcm = 136.560995169
    Pcm = 3953826.947276
    pressure_fraction = 0.4
    pressure = pressure_fraction * Pcm
    m.fs.flash.inlet.pressure.fix(pressure)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.initialize(outlvl=idaeslog.INFO)
    
    m.fs.flash.heat_duty.fix(duty)
    results = solver.solve(m.fs.flash, tee=True, load_solutions=True)
    
    return m


def test_PQ_sweep():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = GenericParameterBlock(default=configuration)
    m.fs.flash = Flash(default={"property_package": m.fs.properties})
    
    iscale.calculate_scaling_factors(m)
    apply_variable_scaling(m)
    apply_constraint_scaling(m)

    # Specify inputs
    components = ['argon', 'nitrogen', 'oxygen']
    P_atm = 101325.0  # [Pa]
    z = {'argon': 0.05, 'nitrogen': 0.60, 'oxygen': 0.35}
    flow_mol = 1.0
    temperature = 95
    duty = 0

    m.fs.flash.inlet.flow_mol.fix(flow_mol)  # 100
    m.fs.flash.inlet.temperature.fix(temperature)  # 300
    m.fs.flash.deltaP.fix(0)
    for comp in components:
        m.fs.flash.inlet.mole_frac_comp[0, comp].fix(z[comp])
    m.fs.flash.inlet.pressure.fix(P_atm)
    m.fs.flash.heat_duty.fix(duty)

    m.fs.flash.initialize(outlvl=idaeslog.INFO)    
    results = solver.solve(m.fs.flash, tee=True)
    m.fs.flash.report()
    
    duty = -1000
    print(f'Solving for Q = {duty}')
    m.fs.flash.heat_duty.fix(duty)
    m.fs.flash.initialize(outlvl=idaeslog.INFO)
    results = solver.solve(m.fs.flash, tee=True)
    
    Tcm = 136.560995169
    Pcm = 3953826.947276
    pressure_fraction = 0.8
    pressure = pressure_fraction * Pcm
    m.fs.flash.inlet.pressure.fix(pressure)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.initialize(outlvl=idaeslog.INFO)

    Q = np.linspace(-1000, 7000, num=50)
    vf = []
    temperature = []
    tbar = []
    for i, duty in enumerate(Q):
        m.fs.flash.inlet.pressure.fix(pressure)
        m.fs.flash.heat_duty.fix(duty)
        
        print("Simulating with Q = ", value(m.fs.flash.heat_duty[0]))
        
        solver.options['bound_push'] = 1e-06  # default = 1e-02
        solver.options['mu_init'] = 1e-06  # default = 1e-01
        results = solver.solve(m.fs.flash, tee=True, load_solutions=True)

        if (results.solver.termination_condition == TerminationCondition.infeasible or
                results.solver.termination_condition ==
                TerminationCondition.maxIterations):
            while solver.options['bound_push'] <= 0.01 and solver.options['mu_init'] <= 0.01:
                solver.options['bound_push'] *= 10
                results = solver.solve(m.fs.flash, tee=True, load_solutions=True)
                if (results.solver.termination_condition ==
                        TerminationCondition.optimal):
                    break
                if solver.options['bound_push'] == 0.1 and solver.options['mu_init'] != 0.1:
                    solver.options['bound_push'] = 1e-6
                    solver.options['mu_init'] *= 10
    
        if (results.solver.termination_condition == TerminationCondition.optimal
            or results.solver.status == SolverStatus.ok):
            vf.append(value(m.fs.flash.vap_outlet.flow_mol[0]))
            tbar.append(
                value(m.fs.flash.control_volume.properties_out[0.0].tbar[('Vap', 'Liq')]))
            temperature.append(
                value(m.fs.flash.control_volume.properties_out[0.0].temperature))
        else:
            print('Solve failed')
            break

    plt.figure()
    plt.plot(temperature, vf, label=f'{pressure_fraction} $P_c$')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Vapor fraction [-]')
    plt.title(f'Vapor fraction variation with temperature, PR-CEOS')
    plt.legend()
    plt.tight_layout()

    return m


def get_jacobian(m):
    m.obj = Objective(expr=0)
    nlp = PyomoNLP(m)
    solvars = nlp.get_pyomo_variables()  # len = 41
    sol = nlp.init_primals()
    nlp.set_primals(sol)
    jac = nlp.evaluate_jacobian().toarray()  # 41x41
    cond = np.linalg.cond(jac)  # condition number
    jac_conditioning = {idx: (min(row), max(row)) for idx, row in enumerate(jac)}
    solcons = nlp.get_pyomo_constraints()
    solcons = {idx: s.name for idx, s in enumerate(solcons)}
    # solcons_dict = {con: con.getname() for con in solcons.values()}
    
    var_dict = {v.parent_component().name: v.parent_component() for v in solvars}
    # incidence_graph_names = {solcons_dict[solcons[i]]: 
    #                     list(set([var.name for var in identify_variables(solcons[i].expr)]))
    #                     for i in range(len(solcons))}
    # incidence_graph = {}
    # for i in range(len(solcons)):
    #     known_vars = ComponentSet()
    #     for var in identify_variables(solcons[i].expr):
    #         if var not in known_vars:
    #             known_vars.add(var)
    #     incidence_graph[solcons_dict[solcons[i]]] = list(known_vars)
    
    return [nlp, jac, jac_conditioning, cond, solcons, solvars, var_dict]


def check_scaling(m):
    jac, nlp = iscale.get_jacobian(m, scaled=True)
    # print("Extreme Jacobian entries:")
    # for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, large=1E3, small=0):
    #     print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]")
    print("Badly scaled variables:")
    for i in iscale.extreme_jacobian_columns(
            jac=jac, nlp=nlp, large=1E3, small=5E-3):
        print(f"    {i[0]:.2e}, [{i[1]}]")
    print("\n\n"+"Badly scaled constraints:")
    for i in iscale.extreme_jacobian_rows(
            jac=jac, nlp=nlp, large=1E3, small=5E-3):
        print(f"    {i[0]:.2e}, [{i[1]}]")
    print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")
    
    # from idaes.core.util.model_diagnostics import DegeneracyHunter
    # if not hasattr(m,"obj"):
    #     m.obj = Objective(expr=0)
    # dh = DegeneracyHunter(m, solver=SolverFactory('cbc'))
    # dh.check_rank_equality_constraints(dense=True)
    # variables = dh.nlp.get_pyomo_variables()
    # #ds = dh.find_candidate_equations()
    # import numpy as np
    # for i in np.where(abs(dh.v[:,-1])>0.1)[0]: 
    #     print(str(i) + ": " + variables[i].name)
    # for i in np.where(abs(dh.u[:,-1])>0.1)[0]: 
    #     print(str(i) + ": " + dh.eq_con_list[i].name)


if __name__ == "__main__":
    # m = test_single_stage_flash()
    m = test_PQ_sweep()
    
    
