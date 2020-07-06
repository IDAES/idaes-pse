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
Make sure the supercritical steam cycle example solves.
"""

__author__ = "Miguel Zamarripa"

import pytest
import pyomo.environ as pyo
import idaes.power_generation.flowsheets.\
    supercritical_power_plant.boiler_subflowsheet_build as blr
import idaes.power_generation.flowsheets.\
    supercritical_power_plant.SCPC_full_plant as SCPC
from idaes.core.util.model_statistics import (degrees_of_freedom)
from idaes.generic_models.properties import iapws95

solver_available = pyo.SolverFactory('ipopt').available()
prop_available = iapws95.iapws95_available()


@pytest.mark.unit
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_init():
    m, solver = blr.main()
    # initialize each unit at the time
    blr.initialize(m)
    blr.unfix_inlets(m)
    # check that the model solved properly and has 0 degrees of freedom
    assert(degrees_of_freedom(m) == 0)


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_boiler():
    m, solver = blr.main()
    # initialize each unit at the time
    blr.initialize(m)
    # unfix inlets to build arcs at the flowsheet level
    blr.unfix_inlets(m)
    m.fs.ATMP1.outlet.enth_mol[0].fix(62710.01)
    m.fs.ATMP1.SprayWater.flow_mol[0].unfix()
    result = solver.solve(m, tee=False)
    assert result.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert result.solver.status == pyo.SolverStatus.ok
    assert pyo.value(m.fs.ECON.side_1.properties_out[0].temperature) == \
        pytest.approx(521.009, 1)


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_power_plan():
    # SCPC.main imports and solves the SCPC Power Plant Flowsheet
    m, results = SCPC.main()
    assert results.solver.status == pyo.SolverStatus.ok
