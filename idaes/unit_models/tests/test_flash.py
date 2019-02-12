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
Tests for Flash unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, SolverFactory, TerminationCondition,
                           SolverStatus, value)

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType, useDefault)
from idaes.unit_models.flash import Flash as FL
from idaes.property_models.BTX_ideal_VLE import IdealParameterBlock
from idaes.ui.report import degrees_of_freedom


# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None

# -----------------------------------------------------------------------------
# Create a flowsheet for test
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})

m.fs.properties = IdealParameterBlock()
m.fs.flash = FL(default={"property_package": m.fs.properties})


def test_build():
    assert len(m.fs.flash.config) == 10
    assert not m.fs.flash.config.has_holdup
    assert m.fs.flash.config.material_balance_type == \
        MaterialBalanceType.componentPhase
    assert m.fs.flash.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.flash.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.flash.config.has_phase_equilibrium
    assert m.fs.flash.config.has_heat_transfer
    assert m.fs.flash.config.has_pressure_change


def test_setInputs():
    m.fs.flash.inlet[0].vars["flow_mol"].fix(1)
    m.fs.flash.inlet[0].vars["temperature"].fix(368)
    m.fs.flash.inlet[0].vars["pressure"].fix(101325)
    m.fs.flash.inlet[0].vars["mole_frac"]["benzene"].fix(0.5)
    m.fs.flash.inlet[0].vars["mole_frac"]["toluene"].fix(0.5)

    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)

    assert degrees_of_freedom(m) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization():
    m.fs.flash.initialize()
    results = solver.solve(m, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert (pytest.approx(0.6038, abs=1e-3) ==
            m.fs.flash.liq_outlet[0].vars["flow_mol"].value)
    assert (pytest.approx(0.3961, abs=1e-3) ==
            m.fs.flash.vap_outlet[0].vars["flow_mol"].value)
    assert (pytest.approx(368, abs=1e-3) ==
            m.fs.flash.liq_outlet[0].vars["temperature"].value)
    assert (pytest.approx(101325, abs=1e-3) ==
            m.fs.flash.liq_outlet[0].vars["pressure"].value)

    assert (pytest.approx(0.4121, abs=1e-3) ==
            value(m.fs.flash.liq_outlet[0].vars["mole_frac"]["benzene"]))
    assert (pytest.approx(0.5878, abs=1e-3) ==
            value(m.fs.flash.liq_outlet[0].vars["mole_frac"]["toluene"]))
    assert (pytest.approx(0.6339, abs=1e-3) ==
            value(m.fs.flash.vap_outlet[0].vars["mole_frac"]["benzene"]))
    assert (pytest.approx(0.3660, abs=1e-3) ==
            value(m.fs.flash.vap_outlet[0].vars["mole_frac"]["toluene"]))
