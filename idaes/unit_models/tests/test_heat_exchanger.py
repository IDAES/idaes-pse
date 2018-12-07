##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for 0D heat exchanger models.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory, value

from idaes.core import FlowsheetBlock
from idaes.unit_models import Heater
from idaes.property_models import iapws95_ph
from idaes.ui.report import degrees_of_freedom

# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6, 'bound_push': 1e-8}
else:
    solver = None

@pytest.fixture()
def build_heater():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95_ph.Iapws95ParameterBlock()
    m.fs.heater = Heater(default={"property_package": m.fs.properties})
    return m

def test_build_heater():
    m = build_heater()
    assert hasattr(m.fs.heater, "inlet")
    assert hasattr(m.fs.heater, "outlet")
    assert len(m.fs.heater.inlet[0].vars) == 3
    assert len(m.fs.heater.outlet[0].vars) == 3

    for port in [m.fs.heater.inlet, m.fs.heater.outlet]:
        assert hasattr(port[0], "flow_mol")
        assert hasattr(port[0], "enth_mol")
        assert hasattr(port[0], "pressure")

@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize_heater():
    m = build_heater()
    init_state = {
        "flow_mol":100,
        "pressure":101325,
        "enth_mol":4000
    }
    m.fs.heater.initialize(state_args=init_state)

    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    assert abs(value(prop_in.temperature) - 326.166978534) <= 1e-4
    assert abs(value(prop_out.temperature) - 326.166978534) <= 1e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_out.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_in.phase_frac["Vap"]) - 0) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0) <= 1e-6

@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_heater_q1():
    m = build_heater()
    init_state = {
        "flow_mol":100,
        "pressure":101325,
        "enth_mol":4000
    }
    m.fs.heater.initialize(state_args=init_state)
    m.fs.heater.control_volume.heat[0].fix(init_state["flow_mol"]*1000)
    prop_in = m.fs.heater.control_volume.properties_in[0]
    prop_out = m.fs.heater.control_volume.properties_out[0]
    prop_in.enth_mol.fix()
    prop_in.flow_mol.fix()
    prop_in.pressure.fix()
    assert degrees_of_freedom(m) == 0
    solver.solve(m)
    assert abs(value(prop_in.temperature) - 326.166978534) <= 1e-4
    assert abs(value(prop_out.temperature) - 333.743257954399) <= 1e-4
    assert abs(value(prop_in.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_out.phase_frac["Liq"]) - 1) <= 1e-6
    assert abs(value(prop_in.phase_frac["Vap"]) - 0) <= 1e-6
    assert abs(value(prop_out.phase_frac["Vap"]) - 0) <= 1e-6
