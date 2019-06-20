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
                        MomentumBalanceType)
from idaes.unit_models.flash import Flash as FL
from idaes.property_models.ideal.BTX_ideal_VLE import BTXParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom


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

m.fs.properties = BTXParameterBlock(default={"valid_phase": ('Liq', 'Vap')})
m.fs.flash = FL(default={"property_package": m.fs.properties})


def test_build():
    assert len(m.fs.flash.config) == 9
    assert not m.fs.flash.config.has_holdup
    assert m.fs.flash.config.material_balance_type == \
        MaterialBalanceType.componentPhase
    assert m.fs.flash.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.flash.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.flash.config.has_heat_transfer
    assert m.fs.flash.config.has_pressure_change


def test_setInputs():
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac[0, "toluene"].fix(0.5)

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

    m.fs.flash.liq_outlet.flow_mol.display()
    assert (pytest.approx(0.645, abs=1e-3) ==
            value(m.fs.flash.liq_outlet.flow_mol[0]))
    assert (pytest.approx(0.355, abs=1e-3) ==
            value(m.fs.flash.vap_outlet.flow_mol[0]))
    assert (pytest.approx(368, abs=1e-3) ==
            value(m.fs.flash.liq_outlet.temperature[0]))
    assert (pytest.approx(101325, abs=1e-3) ==
            value(m.fs.flash.liq_outlet.pressure[0]))

    assert (pytest.approx(0.421, abs=1e-3) ==
            value(m.fs.flash.liq_outlet.mole_frac[0, "benzene"]))
    assert (pytest.approx(0.579, abs=1e-3) ==
            value(m.fs.flash.liq_outlet.mole_frac[0, "toluene"]))
    assert (pytest.approx(0.643, abs=1e-3) ==
            value(m.fs.flash.vap_outlet.mole_frac[0, "benzene"]))
    assert (pytest.approx(0.357, abs=1e-3) ==
            value(m.fs.flash.vap_outlet.mole_frac[0, "toluene"]))
#
#
#def test_report():
#    m.fs.flash.report()
