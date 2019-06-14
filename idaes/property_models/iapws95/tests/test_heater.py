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
Tests for IAPWS using a heater block.  This is an easy way to tes the different
calculation method options.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory, value

from idaes.core import FlowsheetBlock
from idaes.unit_models import Heater, HeatExchanger
from idaes.property_models import iapws95
from idaes.core.util.model_statistics import degrees_of_freedom

iapws95.prop_available = iapws95_available()

# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
else:
    solver = None

@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_heater_q1():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock(
        default={""})
    m.fs.heater = Heater(default={"property_package": m.fs.properties})
    m.fs.heater.inlet.enth_mol.fix(4000)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.inlet.pressure.fix(101325)
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
