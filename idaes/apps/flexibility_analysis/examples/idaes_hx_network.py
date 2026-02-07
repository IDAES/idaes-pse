#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This is a flexibility analysis example using an IDAES model of a 
heat exchanger network. 
"""
import logging
import pyomo.environ as pe
from pyomo.contrib.fbbt.fbbt import fbbt
from pyomo.network import Arc
from pyomo.core.base.block import _BlockData
from pyomo.contrib.solver.util import get_objective
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.core import FlowsheetBlock
from idaes.models.unit_models.heater import Heater

from idaes.core.util.initialization import propagate_state
from idaes.core.base.control_volume_base import ControlVolumeBlockData
import idaes.apps.flexibility_analysis as flexibility
from idaes.apps.flexibility_analysis.var_utils import BoundsManager
from idaes.apps.flexibility_analysis.simplify import simplify_expr


logging.basicConfig(level=logging.INFO)


def create_model():
    """
    This function creates an IDAES model of a
    heat exchanger network for flexibility analysis.
    """
    m = pe.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase="Vap",
        activity_coeff_model="Ideal",
        state_vars="FTPz",
    )

    m.fs.heater1 = Heater(property_package=m.fs.properties)
    m.fs.cooler1 = Heater(property_package=m.fs.properties)
    m.fs.heater2 = Heater(property_package=m.fs.properties)
    m.fs.cooler2 = Heater(property_package=m.fs.properties)
    m.fs.heater3 = Heater(property_package=m.fs.properties)
    m.fs.cooler3 = Heater(property_package=m.fs.properties)
    m.fs.cooler4 = Heater(property_package=m.fs.properties)

    m.fs.s1 = Arc(source=m.fs.heater1.outlet, destination=m.fs.heater2.inlet)
    m.fs.s2 = Arc(source=m.fs.cooler1.outlet, destination=m.fs.cooler4.inlet)
    m.fs.s3 = Arc(source=m.fs.cooler2.outlet, destination=m.fs.cooler3.inlet)

    m.fs.duty_cons = pe.ConstraintList()
    m.fs.duty_cons.add(m.fs.heater1.heat_duty[0] == -m.fs.cooler1.heat_duty[0])
    m.fs.duty_cons.add(m.fs.heater2.heat_duty[0] == -m.fs.cooler2.heat_duty[0])
    m.fs.duty_cons.add(m.fs.heater3.heat_duty[0] == -m.fs.cooler3.heat_duty[0])
    m.fs.duty_cons.add(m.fs.cooler1.heat_duty[0] <= 0)
    m.fs.duty_cons.add(m.fs.cooler2.heat_duty[0] <= 0)
    m.fs.duty_cons.add(m.fs.cooler3.heat_duty[0] <= 0)
    m.fs.duty_cons.add(m.fs.cooler4.heat_duty[0] <= 0)

    approach_limit = 5
    m.fs.temp_approach = pe.ConstraintList()
    m.fs.temp_approach.add(
        m.fs.cooler1.inlet.temperature[0]
        >= m.fs.heater1.outlet.temperature[0] + approach_limit
    )
    m.fs.temp_approach.add(
        m.fs.cooler1.outlet.temperature[0]
        >= m.fs.heater1.inlet.temperature[0] + approach_limit
    )
    m.fs.temp_approach.add(
        m.fs.cooler2.inlet.temperature[0]
        >= m.fs.heater2.outlet.temperature[0] + approach_limit
    )
    m.fs.temp_approach.add(
        m.fs.cooler2.outlet.temperature[0]
        >= m.fs.heater2.inlet.temperature[0] + approach_limit
    )
    m.fs.temp_approach.add(
        m.fs.cooler3.inlet.temperature[0]
        >= m.fs.heater3.outlet.temperature[0] + approach_limit
    )
    m.fs.temp_approach.add(
        m.fs.cooler3.outlet.temperature[0]
        >= m.fs.heater3.inlet.temperature[0] + approach_limit
    )

    # specs
    m.fs.heater1.inlet.temperature[0].fix(488)
    m.fs.cooler1.inlet.temperature[0].fix(720)
    m.fs.cooler2.inlet.temperature[0].fix(683)
    m.fs.heater2.outlet.temperature[0].fix(663)
    m.fs.heater3.inlet.temperature[0].fix(400)
    m.fs.heater3.outlet.temperature[0].fix(550)
    m.fs.cooler4.outlet.temperature[0].fix(450)
    m.fs.cooler3_out_temp_spec = pe.Constraint(
        expr=m.fs.cooler3.outlet.temperature[0] <= 423
    )

    m.fs.heater1.inlet.flow_mol[0].fix(1)
    m.fs.cooler1.inlet.flow_mol[0].fix(1)
    m.fs.cooler2.inlet.flow_mol[0].fix(1)
    # m.fs.heater3.inlet.flow_mol[0].fix(1)
    m.fs.heater3.inlet.flow_mol[0].setlb(0.8)
    m.fs.heater3.inlet.flow_mol[0].setub(1.2)

    m.fs.heater1.inlet.pressure[0].fix(101325)
    m.fs.cooler1.inlet.pressure[0].fix(101325)
    m.fs.cooler2.inlet.pressure[0].fix(101325)
    m.fs.heater3.inlet.pressure[0].fix(101325)

    m.fs.heater1.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.cooler1.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.cooler2.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.heater3.inlet.mole_frac_comp[0, "benzene"].fix(0.5)

    m.fs.heater1.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.cooler1.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.cooler2.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.heater3.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

    m.obj = pe.Objective(expr=-m.fs.cooler4.heat_duty[0])

    pe.TransformationFactory("network.expand_arcs").apply_to(m)
    scale_model(m)

    nominal_values = pe.ComponentMap()
    # nominal_values[m.fs.cooler1.inlet.temperature[0]] = 7.20
    nominal_values[m.fs.heater1.inlet.temperature[0]] = 4.88
    nominal_values[m.fs.cooler2.inlet.temperature[0]] = 6.83
    # nominal_values[m.fs.heater3.inlet.temperature[0]] = 4.00
    nominal_values[m.fs.heater1.inlet.flow_mol[0]] = 1
    # nominal_values[m.fs.cooler1.inlet.flow_mol[0]] = 1
    nominal_values[m.fs.cooler2.inlet.flow_mol[0]] = 1
    # nominal_values[m.fs.heater3.inlet.flow_mol[0]] = 1

    param_bounds = pe.ComponentMap()
    # p = m.fs.cooler1.inlet.temperature[0]
    # param_bounds[p] = (nominal_values[p] - 0.10, nominal_values[p] + 0.10)
    p = m.fs.heater1.inlet.temperature[0]
    param_bounds[p] = (nominal_values[p] - 0.10, nominal_values[p] + 0.10)
    p = m.fs.cooler2.inlet.temperature[0]
    param_bounds[p] = (nominal_values[p] - 0.10, nominal_values[p] + 0.10)
    # p = m.fs.heater3.inlet.temperature[0]
    # param_bounds[p] = (nominal_values[p] - 0.10, nominal_values[p] + 0.10)
    p = m.fs.heater1.inlet.flow_mol[0]
    param_bounds[p] = (nominal_values[p] - 0.2, nominal_values[p] + 0.2)
    # p = m.fs.cooler1.inlet.flow_mol[0]
    # param_bounds[p] = (nominal_values[p] - 0.2, nominal_values[p] + 0.2)
    p = m.fs.cooler2.inlet.flow_mol[0]
    param_bounds[p] = (nominal_values[p] - 0.2, nominal_values[p] + 0.2)
    # p = m.fs.heater3.inlet.flow_mol[0]
    # param_bounds[p] = (nominal_values[p] - 0.2, nominal_values[p] + 0.2)

    for p in nominal_values.keys():
        assert p.is_fixed()
        p.unfix()

    for c in m.component_data_objects(pe.Constraint, active=True, descend_into=True):
        body = simplify_expr(c.body)
        c.set_value((c.lower, body, c.upper))

    for p in nominal_values.keys():
        assert not p.is_fixed()
        p.fix()

    return m, nominal_values, param_bounds


def initialize(m):
    """
    Initialize the model
    """
    m.fs.heater1.heat_duty[0].fix(0.1000)
    m.fs.heater1.initialize()
    m.fs.heater1.heat_duty[0].unfix()
    propagate_state(m.fs.s1)
    m.fs.cooler1.heat_duty[0].fix(-m.fs.heater1.heat_duty[0])
    m.fs.cooler1.initialize()
    m.fs.cooler1.heat_duty[0].unfix()
    m.fs.heater2.inlet.temperature[0].fix()
    m.fs.heater2.initialize()
    m.fs.heater2.inlet.temperature[0].unfix()
    m.fs.cooler2.heat_duty[0].fix(-m.fs.heater2.heat_duty[0].value)
    m.fs.cooler2.initialize()
    m.fs.cooler2.heat_duty[0].unfix()
    m.fs.heater3.initialize()
    propagate_state(m.fs.s3)
    m.fs.cooler3.inlet.temperature[0].fix()
    m.fs.cooler3.heat_duty[0].fix(-m.fs.heater3.heat_duty[0].value)
    m.fs.cooler3.initialize()
    m.fs.cooler3.heat_duty[0].unfix()
    m.fs.cooler3.inlet.temperature[0].unfix()
    propagate_state(m.fs.s2)
    m.fs.cooler4.inlet.temperature[0].fix()
    m.fs.cooler4.initialize()
    m.fs.cooler4.inlet.temperature[0].unfix()


def get_var_bounds(m: _BlockData, param_bounds):
    """
    Generate a map with valid variable bounds for
    any possible realization of the uncertain parameters
    """
    for p, (p_lb, p_ub) in param_bounds.items():
        p.unfix()
        p.setlb(p_lb)
        p.setub(p_ub)
    bounds_manager = BoundsManager(m)
    bounds_manager.save_bounds()
    for b in m.block_data_objects(active=True, descend_into=True):
        if isinstance(b, ControlVolumeBlockData):
            b.properties_in[0].flow_mol.setlb(0)
            b.properties_in[0].flow_mol.setub(5)
            b.properties_out[0].flow_mol.setlb(0)
            b.properties_out[0].flow_mol.setub(5)
            b.properties_in[0].flow_mol_phase["Vap"].setlb(0)
            b.properties_in[0].flow_mol_phase["Vap"].setub(5)
            b.properties_out[0].flow_mol_phase["Vap"].setlb(0)
            b.properties_out[0].flow_mol_phase["Vap"].setub(5)
            b.properties_in[0].mole_frac_comp["benzene"].setlb(0)
            b.properties_in[0].mole_frac_comp["benzene"].setub(1)
            b.properties_in[0].mole_frac_comp["toluene"].setlb(0)
            b.properties_in[0].mole_frac_comp["toluene"].setub(1)
            b.properties_out[0].mole_frac_comp["benzene"].setlb(0)
            b.properties_out[0].mole_frac_comp["benzene"].setub(1)
            b.properties_out[0].mole_frac_comp["toluene"].setlb(0)
            b.properties_out[0].mole_frac_comp["toluene"].setub(1)
            b.properties_in[0].mole_frac_phase_comp["Vap", "benzene"].setlb(0)
            b.properties_in[0].mole_frac_phase_comp["Vap", "benzene"].setub(1)
            b.properties_in[0].mole_frac_phase_comp["Vap", "toluene"].setlb(0)
            b.properties_in[0].mole_frac_phase_comp["Vap", "toluene"].setub(1)
            b.properties_out[0].mole_frac_phase_comp["Vap", "benzene"].setlb(0)
            b.properties_out[0].mole_frac_phase_comp["Vap", "benzene"].setub(1)
            b.properties_out[0].mole_frac_phase_comp["Vap", "toluene"].setlb(0)
            b.properties_out[0].mole_frac_phase_comp["Vap", "toluene"].setub(1)
            b.properties_in[0].pressure.setlb(1.01325)
            b.properties_in[0].pressure.setub(1.01325)
            b.properties_out[0].pressure.setlb(1.01325)
            b.properties_out[0].pressure.setub(1.01325)
            b.properties_in[0].temperature.setlb(3.00)
            b.properties_in[0].temperature.setub(10.00)
            b.properties_out[0].temperature.setlb(3.00)
            b.properties_out[0].temperature.setub(10.00)
            cons_to_fbbt = list()
            cons_to_fbbt.extend(b.properties_in[0].eq_enth_mol_phase_comp.values())
            cons_to_fbbt.extend(b.properties_out[0].eq_enth_mol_phase_comp.values())
            cons_to_fbbt.extend(b.properties_in[0].eq_enth_mol_phase.values())
            cons_to_fbbt.extend(b.properties_out[0].eq_enth_mol_phase.values())
            cons_to_fbbt.extend(b.enthalpy_balances.values())
            for c in cons_to_fbbt:
                assert c.active
                if c.active:
                    fbbt(c)
                else:
                    c.activate()
                    fbbt(c)
                    c.deactivate()
    res = pe.ComponentMap()
    for v in m.component_data_objects(pe.Var, descend_into=True):
        res[v] = (v.lb, v.ub)
    bounds_manager.pop_bounds()
    return res


def scale_model(m):
    """
    This function scales the heat exchanger network model
    prior to performing flexibility analysis

    Parameters
    ----------
    m: _BlockData
        The IDAES model of the heat exchanger network
    """
    m.scaling_factor = pe.Suffix(direction=pe.Suffix.EXPORT)
    for b in m.block_data_objects(active=True, descend_into=True):
        if isinstance(b, ControlVolumeBlockData):
            m.scaling_factor[b.heat[0]] = 1e-4
            m.scaling_factor[b.properties_in[0].pressure] = 1e-5
            m.scaling_factor[b.properties_out[0].pressure] = 1e-5
            m.scaling_factor[b.properties_in[0].temperature] = 1e-2
            m.scaling_factor[b.properties_out[0].temperature] = 1e-2
            m.scaling_factor[b.properties_in[0].enth_mol_phase["Vap"]] = 1e-4
            m.scaling_factor[b.properties_out[0].enth_mol_phase["Vap"]] = 1e-4
            m.scaling_factor[
                b.properties_in[0].enth_mol_phase_comp["Vap", "benzene"]
            ] = 1e-4
            m.scaling_factor[
                b.properties_out[0].enth_mol_phase_comp["Vap", "benzene"]
            ] = 1e-4
            m.scaling_factor[
                b.properties_in[0].enth_mol_phase_comp["Vap", "toluene"]
            ] = 1e-4
            m.scaling_factor[
                b.properties_out[0].enth_mol_phase_comp["Vap", "toluene"]
            ] = 1e-4
            for c in b.enthalpy_balances.values():
                m.scaling_factor[c] = 1e-4
            for c in b.pressure_balance.values():
                m.scaling_factor[c] = 1e-4
            for c in b.properties_in[0].eq_enth_mol_phase.values():
                m.scaling_factor[c] = 1e-4
            for c in b.properties_out[0].eq_enth_mol_phase.values():
                m.scaling_factor[c] = 1e-4
            for c in b.properties_in[0].eq_enth_mol_phase_comp.values():
                m.scaling_factor[c] = 1e-4
            for c in b.properties_out[0].eq_enth_mol_phase_comp.values():
                m.scaling_factor[c] = 1e-4
    for c in m.fs.duty_cons.values():
        m.scaling_factor[c] = 1e-4
    m.scaling_factor[m.fs.s1_expanded.pressure_equality[0]] = 1e-4
    m.scaling_factor[m.fs.s2_expanded.pressure_equality[0]] = 1e-4
    m.scaling_factor[m.fs.s3_expanded.pressure_equality[0]] = 1e-4
    m.scaling_factor[get_objective(m)] = 1e-4

    pe.TransformationFactory("core.scale_model").apply_to(m, rename=False)


def main(method: flexibility.FlexTestMethod):
    """
    Run the example

    Parameters
    ----------
    method: flexibility.FlexTestMethod
        The method to use for the flexibility test
    """
    m, nominal_values, param_bounds = create_model()
    initialize(m)
    var_bounds = get_var_bounds(m, param_bounds)
    config = flexibility.FlexTestConfig()
    config.feasibility_tol = 1e-6
    config.terminate_early = False
    config.method = method
    config.minlp_solver = pe.SolverFactory("scip")
    config.minlp_solver.options["limits/time"] = 300
    config.sampling_config.solver = pe.SolverFactory("ipopt")
    config.sampling_config.strategy = flexibility.SamplingStrategy.grid
    config.sampling_config.num_points = 3
    if method == flexibility.FlexTestMethod.linear_decision_rule:
        config.decision_rule_config = flexibility.LinearDRConfig()
        config.decision_rule_config.solver = pe.SolverFactory("scip")
    elif method == flexibility.FlexTestMethod.relu_decision_rule:
        config.decision_rule_config = flexibility.ReluDRConfig()
        config.decision_rule_config.n_layers = 1
        config.decision_rule_config.n_nodes_per_layer = 10
        config.decision_rule_config.epochs = 10000
        config.decision_rule_config.batch_size = 50
        config.decision_rule_config.scale_inputs = True
        config.decision_rule_config.scale_outputs = True
    results = flexibility.solve_flextest(
        m=m,
        uncertain_params=list(nominal_values.keys()),
        param_nominal_values=nominal_values,
        param_bounds=param_bounds,
        controls=[
            m.fs.cooler4.control_volume.heat[0],
            m.fs.heater3.control_volume.properties_in[0].flow_mol,
        ],
        valid_var_bounds=var_bounds,
        config=config,
    )
    print(results)
    return results


if __name__ == "__main__":
    main(flexibility.FlexTestMethod.sampling)
