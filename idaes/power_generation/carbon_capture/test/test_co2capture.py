###############################################################################
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
###############################################################################

"""
Test for the CO2Pure Unit based on surrogates

"""
import pytest
# Import Pyomo libraries
import pyomo.environ as pyo
# Import IDAES core
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.core.util.testing import get_default_solver
from idaes.power_generation.carbon_capture.piperazine_surrogates.\
    co2_capture_system import CO2Capture
from idaes.core.util.model_statistics import degrees_of_freedom
# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()

# -----------------------------------------------------------------------------
__author__ = "M. Zamarripa"
__version__ = "1.0.0"


@pytest.fixture(scope="module")
def build_ccs():
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing(year='2018')

    # CO2Pure unit based on surrogates
    m.fs.CCS_system = CO2Capture()

    # Fix the inlet states for solving
    # NETL Baseline report Exhibit 5-22 B31B Case
    # Flue gas (Stream 4) is 138,406 kgmol/hr
    fluegas = 38446.11  # mol/s  --> 138406 kmol/hr *1000 mol/kmol 1hr/3600s
    m.fs.CCS_system.inlet.temperature[:].fix(303.15)  # K (30 C)
    m.fs.CCS_system.inlet.pressure[:].fix(101325)  # Pa (1 atm)
    m.fs.CCS_system.inlet.flow_mol_comp[:, 'CO2'].fix(fluegas*0.0408)
    m.fs.CCS_system.inlet.flow_mol_comp[:, 'O2'].fix(fluegas*0.12)
    m.fs.CCS_system.inlet.flow_mol_comp[:, 'Ar'].fix(fluegas*0.0089)
    m.fs.CCS_system.inlet.flow_mol_comp[:, 'H2O'].fix(fluegas*0.0875)
    m.fs.CCS_system.inlet.flow_mol_comp[:, 'N2'].fix(fluegas*0.7428)
    m.fs.CCS_system.CO2_capture_rate.fix(0.9)  # 90 % CO2 Capture
    m.fs.CCS_system.Pz_mol.fix(5)
    m.fs.CCS_system.lean_loading.fix(0.25)
    return m


@pytest.mark.unit
def test_basic_build(build_ccs):
    """Make a model and make sure it doesn't throw exception"""
    m = build_ccs
    assert degrees_of_freedom(m) == 0
    assert isinstance(m.fs.CCS_system.CO2_capture_rate, pyo.Var)
    assert isinstance(m.fs.CCS_system.SRD, pyo.Var)
    # Check unit config arguments


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_run_model(build_ccs):
    m = build_ccs
    # Initialize the unit
    m.fs.CCS_system.initialize(outlvl=idaeslog.INFO)

    # Solve the model
    results = solver.solve(m, tee=True)
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    assert degrees_of_freedom(m) == 0
    assert (pytest.approx(4.042, abs=1e-3) ==
            pyo.value(m.fs.CCS_system.SRD[0]))
    assert (pytest.approx(251.188, abs=1e-3) ==
            pyo.value(m.fs.CCS_system.reboiler_duty[0]))
