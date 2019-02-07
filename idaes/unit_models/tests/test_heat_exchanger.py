##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for 0D heat exchanger models.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory, value

from idaes.core import FlowsheetBlock
from idaes.unit_models import Heater, HeatExchanger
from idaes.property_models import iapws95_ph
from idaes.ui.report import degrees_of_freedom

# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
else:
    solver = None

@pytest.fixture()
def build_heater():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95_ph.Iapws95ParameterBlock()
    m.fs.heater = Heater(default={"property_package": m.fs.properties})
    return m

@pytest.fixture()
def build_heat_exchanger():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95_ph.Iapws95ParameterBlock()
    m.fs.heat_exchanger = HeatExchanger(default={
        "side_1":{"property_package": m.fs.properties},
        "side_2":{"property_package": m.fs.properties}})
    return m

def test_build_heat_exchanger(build_heat_exchanger):
    m = build_heat_exchanger
    assert hasattr(m.fs.heat_exchanger, "inlet_1")
    assert hasattr(m.fs.heat_exchanger, "outlet_1")
    assert hasattr(m.fs.heat_exchanger, "inlet_2")
    assert hasattr(m.fs.heat_exchanger, "outlet_2")

def test_initialize_heat_exchanger(build_heat_exchanger):
    m = build_heat_exchanger
    init_state1 = {
        "flow_mol":100,
        "pressure":101325,
        "enth_mol":4000}
    init_state2 = {
        "flow_mol":100,
        "pressure":101325,
        "enth_mol":3500}

    prop_in_1 = m.fs.heat_exchanger.side_1.properties_in[0]
    prop_out_1 = m.fs.heat_exchanger.side_1.properties_out[0]
    prop_in_2 = m.fs.heat_exchanger.side_2.properties_in[0]
    prop_out_2 = m.fs.heat_exchanger.side_2.properties_out[0]
    prop_in_1.flow_mol.fix(100)
    prop_in_1.pressure.fix(101325)
    prop_in_1.enth_mol.fix(4000)
    prop_in_2.flow_mol.fix(100)
    prop_in_2.pressure.fix(101325)
    prop_in_2.enth_mol.fix(3000)

    m.fs.heat_exchanger.heat_duty.value = 10000
    m.fs.heat_exchanger.initialize(state_args_1=init_state1,
                                   state_args_2=init_state2,
                                   outlvl=5)
    solver.solve(m)
    assert degrees_of_freedom(m) == 0
    print(value(m.fs.heat_exchanger.delta_temperature[0]))
    print(value(m.fs.heat_exchanger.side_1.heat[0]))
    print(value(m.fs.heat_exchanger.side_2.heat[0]))
    assert abs(value(prop_in_1.temperature) - 326.1667075078748) <= 1e-4
    assert abs(value(prop_out_1.temperature) - 313.81921851031814) <= 1e-4
    assert abs(value(prop_in_2.temperature) -  312.88896252921734) <= 1e-4
    assert abs(value(prop_out_2.temperature) - 325.23704823703537) <= 1e-4
    assert abs(value(prop_in_1.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_out_1.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_in_1.phase_frac["Vap"]) - 0) <= 1e-6
    assert abs(value(prop_out_1.phase_frac["Vap"]) - 0) <= 1e-6

def test_build_heater(build_heater):
    m = build_heater
    assert hasattr(m.fs.heater, "inlet")
    assert hasattr(m.fs.heater, "outlet")
    assert len(m.fs.heater.inlet[0].vars) == 3
    assert len(m.fs.heater.outlet[0].vars) == 3

    for port in [m.fs.heater.inlet, m.fs.heater.outlet]:
        assert hasattr(port[0], "flow_mol")
        assert hasattr(port[0], "enth_mol")
        assert hasattr(port[0], "pressure")

@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize_heater(build_heater):
    m = build_heater
    m.fs.heater.inlet[:].enth_mol.fix(4000)
    m.fs.heater.inlet[:].flow_mol.fix(100)
    m.fs.heater.inlet[:].pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*20000)
    m.fs.heater.initialize()

    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert abs(value(prop_in.temperature) - 326.1667075078748) <= 1e-4
    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
    assert abs(value(prop_in.phase_frac["Vap"]) - 0) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6

@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_heater_q1(build_heater):
    m = build_heater
    m.fs.heater.inlet[:].enth_mol.fix(4000)
    m.fs.heater.inlet[:].flow_mol.fix(100)
    m.fs.heater.inlet[:].pressure.fix(101325)
    m.fs.heater.heat_duty[0].fix(100*20000)
    m.fs.heater.initialize()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert abs(value(prop_in.temperature) - 326.1667075078748) <= 1e-4
    assert abs(value(prop_out.temperature) - 373.12429584768876) <= 1e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_out.phase_frac["Liq"]) - 0.5953218682380845) <= 1e-6
    assert abs(value(prop_in.phase_frac["Vap"]) - 0) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0.40467813176191547) <= 1e-6
