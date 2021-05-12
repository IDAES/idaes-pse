##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for turbine multistage model.

Author: John Eslick
"""
import pytest

import pyomo.environ as pyo
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Heater
from idaes.power_generation.unit_models.helm import HelmTurbineMultistage
from idaes.generic_models.properties import iapws95
from idaes.core.util.model_statistics import (
    degrees_of_freedom, activated_equalities_generator)
import idaes.core.util.scaling as iscale
from idaes.core.util import get_solver

# Set up solver
solver = get_solver()


@pytest.mark.unit
def build_turbine_for_run_test():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    # roughly based on NETL baseline studies
    m.fs.turb = HelmTurbineMultistage(default={
        "property_package": m.fs.properties,
        "num_hp": 7,
        "num_ip": 14,
        "num_lp": 11,
        "hp_split_locations": [4,7],
        "ip_split_locations": [5, 14],
        "lp_split_locations": [4,7,9,11],
        "hp_disconnect": [7],
        "ip_split_num_outlets": {14:3}})

    # Add reheater
    m.fs.reheat = Heater(default={"property_package": m.fs.properties})
    m.fs.hp_to_reheat = Arc(source=m.fs.turb.hp_split[7].outlet_1,
                            destination=m.fs.reheat.inlet)
    m.fs.reheat_to_ip = Arc(source=m.fs.reheat.outlet,
                            destination=m.fs.turb.ip_stages[1].inlet)
    return m


@pytest.mark.component
def test_initialize():
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_turbine_for_run_test()
    turb = m.fs.turb

    # Set the inlet of the turbine
    p = 2.4233e7
    hin = pyo.value(iapws95.htpx(T=880*pyo.units.K, P=p*pyo.units.Pa))
    m.fs.turb.inlet_split.inlet.enth_mol[0].fix(hin)
    m.fs.turb.inlet_split.inlet.flow_mol[0].fix(26000)
    m.fs.turb.inlet_split.inlet.pressure[0].fix(p)

    # Set the inlet of the ip section, which is disconnected
    # here to insert reheater
    p = 7.802e+06
    hin = pyo.value(iapws95.htpx(T=880*pyo.units.K, P=p*pyo.units.Pa))
    m.fs.turb.ip_stages[1].inlet.enth_mol[0].value = hin
    m.fs.turb.ip_stages[1].inlet.flow_mol[0].value = 25220.0
    m.fs.turb.ip_stages[1].inlet.pressure[0].value = p

    for i, s in turb.hp_stages.items():
        s.ratioP[:] = 0.88
        s.efficiency_isentropic[:] = 0.9
    for i, s in turb.ip_stages.items():
        s.ratioP[:] = 0.85
        s.efficiency_isentropic[:] = 0.9
    for i, s in turb.lp_stages.items():
        s.ratioP[:] = 0.82
        s.efficiency_isentropic[:] = 0.9

    turb.hp_split[4].split_fraction[0,"outlet_2"].fix(0.03)
    turb.hp_split[7].split_fraction[0,"outlet_2"].fix(0.03)
    turb.ip_split[5].split_fraction[0,"outlet_2"].fix(0.04)
    turb.ip_split[14].split_fraction[0,"outlet_2"].fix(0.04)
    turb.ip_split[14].split_fraction[0,"outlet_3"].fix(0.15)
    turb.lp_split[4].split_fraction[0,"outlet_2"].fix(0.04)
    turb.lp_split[7].split_fraction[0,"outlet_2"].fix(0.04)
    turb.lp_split[9].split_fraction[0,"outlet_2"].fix(0.04)
    turb.lp_split[11].split_fraction[0,"outlet_2"].fix(0.04)

    # Congiure with reheater for a full test
    turb.ip_stages[1].inlet.fix()
    for i in turb.throttle_valve:
        turb.throttle_valve[i].Cv.fix()
        turb.throttle_valve[i].valve_opening.fix()
    turb.inlet_split.inlet.flow_mol.unfix()
    turb.inlet_mix.use_equal_pressure_constraint()

    iscale.calculate_scaling_factors(m)
    turb.initialize(outlvl=1)
    turb.ip_stages[1].inlet.unfix()


    for t in m.fs.time:
        m.fs.reheat.inlet.flow_mol[t].value = \
            pyo.value(turb.hp_split[7].outlet_1_state[t].flow_mol)
        m.fs.reheat.inlet.enth_mol[t].value = \
            pyo.value(turb.hp_split[7].outlet_1_state[t].enth_mol)
        m.fs.reheat.inlet.pressure[t].value = \
            pyo.value(turb.hp_split[7].outlet_1_state[t].pressure)
    m.fs.reheat.initialize(outlvl=4)
    def reheat_T_rule(b, t):
        return m.fs.reheat.control_volume.properties_out[t].temperature == 880
    m.fs.reheat.temperature_out_equation = pyo.Constraint(
            m.fs.reheat.flowsheet().config.time,
            rule=reheat_T_rule)

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    m.fs.turb.outlet_stage.control_volume.properties_out[0].pressure.fix()

    assert degrees_of_freedom(m)==0
    solver.solve(m, tee=True)

    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4

    return m


@pytest.mark.component
def test_initialize_calc_cf():
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_turbine_for_run_test()
    turb = m.fs.turb

    # Set the inlet of the turbine
    p = 2.4233e7
    hin = pyo.value(iapws95.htpx(T=880*pyo.units.K, P=p*pyo.units.Pa))
    m.fs.turb.inlet_split.inlet.enth_mol[0].fix(hin)
    m.fs.turb.inlet_split.inlet.flow_mol[0].fix(26000)
    m.fs.turb.inlet_split.inlet.pressure[0].fix(p)

    # Set the inlet of the ip section, which is disconnected
    # here to insert reheater
    p = 7.802e+06
    hin = pyo.value(iapws95.htpx(T=880*pyo.units.K, P=p*pyo.units.Pa))
    m.fs.turb.ip_stages[1].inlet.enth_mol[0].value = hin
    #m.fs.turb.ip_stages[1].inlet.flow_mol[0].value = 25220.0
    m.fs.turb.ip_stages[1].inlet.pressure[0].value = p

    for i, s in turb.hp_stages.items():
        s.ratioP[:] = 0.88
        s.efficiency_isentropic[:] = 0.9
    for i, s in turb.ip_stages.items():
        s.ratioP[:] = 0.85
        s.efficiency_isentropic[:] = 0.9
    for i, s in turb.lp_stages.items():
        s.ratioP[:] = 0.82
        s.efficiency_isentropic[:] = 0.9

    turb.hp_split[4].split_fraction[0,"outlet_2"].fix(0.03)
    turb.hp_split[7].split_fraction[0,"outlet_2"].fix(0.03)
    turb.ip_split[5].split_fraction[0,"outlet_2"].fix(0.04)
    turb.ip_split[14].split_fraction[0,"outlet_2"].fix(0.04)
    turb.ip_split[14].split_fraction[0,"outlet_3"].fix(0.15)
    turb.lp_split[4].split_fraction[0,"outlet_2"].fix(0.04)
    turb.lp_split[7].split_fraction[0,"outlet_2"].fix(0.04)
    turb.lp_split[9].split_fraction[0,"outlet_2"].fix(0.04)
    turb.lp_split[11].split_fraction[0,"outlet_2"].fix(0.04)

    # Congiure with reheater for a full test
    turb.ip_stages[1].inlet.fix()
    turb.inlet_split.inlet.flow_mol.unfix()
    turb.inlet_mix.use_equal_pressure_constraint()
    for i in m.fs.turb.inlet_stage:
        m.fs.turb.inlet_stage[i].ratioP[0] = 0.6
        turb.throttle_valve[i].Cv.fix()
        turb.throttle_valve[i].valve_opening.fix()
    iscale.calculate_scaling_factors(m)
    turb.initialize(outlvl=1, calculate_inlet_cf=True, calculate_outlet_cf=True)
    turb.ip_stages[1].inlet.unfix()

    for t in m.fs.time:
        m.fs.reheat.inlet.flow_mol[t].value = \
            pyo.value(turb.hp_split[7].outlet_1_state[t].flow_mol)
        m.fs.reheat.inlet.enth_mol[t].value = \
            pyo.value(turb.hp_split[7].outlet_1_state[t].enth_mol)
        m.fs.reheat.inlet.pressure[t].value = \
            pyo.value(turb.hp_split[7].outlet_1_state[t].pressure)
    m.fs.reheat.initialize(outlvl=4)
    def reheat_T_rule(b, t):
        return m.fs.reheat.control_volume.properties_out[t].temperature == 880
    m.fs.reheat.temperature_out_equation = pyo.Constraint(
            m.fs.reheat.flowsheet().config.time,
            rule=reheat_T_rule)

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    m.fs.turb.outlet_stage.control_volume.properties_out[0].pressure.fix()

    assert degrees_of_freedom(m)==0
    solver.solve(m, tee=True)

    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4

    assert pyo.value(m.fs.turb.inlet_split.inlet.flow_mol[0]) == pytest.approx(26000)

    return m
