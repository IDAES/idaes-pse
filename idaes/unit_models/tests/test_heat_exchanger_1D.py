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
Tests for Heat Exchanger 1D unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, SolverFactory, TerminationCondition,
                           SolverStatus)

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType)
from idaes.unit_models.heat_exchanger_1D import HeatExchanger1D as HX1D

from idaes.property_models.examples.BFW_properties import PhysicalParameterBlock
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

m.fs.properties = PhysicalParameterBlock()
m.fs.HX1D = HX1D(default={"shell_side": {"property_package": m.fs.properties},
                          "tube_side": {"property_package": m.fs.properties}})


def test_build():
    # Check for default setting for config block attributes
    assert len(m.fs.HX1D.config) == 7
    assert m.fs.HX1D.config.flow_type == "co_current"
    assert m.fs.HX1D.config.has_wall_conduction == "none"

    assert len(m.fs.HX1D.config.shell_side) == 12
    assert not m.fs.HX1D.config.shell_side.has_holdup
    assert m.fs.HX1D.config.shell_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX1D.config.shell_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX1D.config.shell_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.HX1D.config.shell_side.has_heat_transfer
    assert not m.fs.HX1D.config.shell_side.has_pressure_change
    assert not m.fs.HX1D.config.shell_side.has_phase_equilibrium

    assert len(m.fs.HX1D.config.tube_side) == 12
    assert not m.fs.HX1D.config.tube_side.has_holdup
    assert m.fs.HX1D.config.tube_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX1D.config.tube_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX1D.config.tube_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.HX1D.config.tube_side.has_heat_transfer
    assert not m.fs.HX1D.config.tube_side.has_pressure_change
    assert not m.fs.HX1D.config.tube_side.has_phase_equilibrium




    # Check for inlets/outlets construction
    assert hasattr(m.fs.HX1D, "shell_inlet")
    assert hasattr(m.fs.HX1D, "shell_outlet")
    assert hasattr(m.fs.HX1D, "tube_inlet")
    assert hasattr(m.fs.HX1D, "tube_outlet")

def test_setInputs():
    """Set inlet and operating conditions."""
    # HX dimensions
    m.fs.HX1D.d_shell.fix(1.04)
    m.fs.HX1D.d_tube_outer.fix(0.01167)
    m.fs.HX1D.d_tube_inner.fix(0.01067)
    m.fs.HX1D.N_tubes.fix(1176)
    m.fs.HX1D.shell_length.fix(4.85)
    m.fs.HX1D.tube_length.fix(4.85)
    m.fs.HX1D.shell_heat_transfer_coefficient.fix(2000)
    m.fs.HX1D.tube_heat_transfer_coefficient.fix(51000)

    # Boundary conditions
    for t in m.fs.time:
        # Unit HX01
        # NETL baseline study
        m.fs.HX1D.shell_inlet[t].vars["flow_mol"].fix(2300)  # mol/s
        m.fs.HX1D.shell_inlet[t].vars["temperature"].fix(676)  # K
        m.fs.HX1D.shell_inlet[t].vars["pressure"].fix(7.38E6)  # Pa
        m.fs.HX1D.shell_inlet[t].vars["vapor_frac"].fix(1)

        m.fs.HX1D.tube_inlet[t].vars["flow_mol"].fix(26.6)  # mol/s
        m.fs.HX1D.tube_inlet[t].vars["temperature"].fix(529)  # K
        m.fs.HX1D.tube_inlet[t].vars["pressure"].fix(2.65E7)  # Pa
        m.fs.HX1D.tube_inlet[t].vars["vapor_frac"].fix(0)

    assert degrees_of_freedom(m) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization():
    m.fs.HX1D.initialize()
    results = solver.solve(m, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert (pytest.approx(2300, abs=1e-3) ==
            m.fs.HX1D.shell_outlet[0].vars["flow_mol"].value)
    assert (pytest.approx(559.220, abs=1e-3) ==
            m.fs.HX1D.shell_outlet[0].vars["temperature"].value)
    assert (pytest.approx(7.38E6, abs=1e-3) ==
            m.fs.HX1D.shell_outlet[0].vars["pressure"].value)

    assert (pytest.approx(26.6, abs=1e-3) ==
            m.fs.HX1D.tube_outlet[0].vars["flow_mol"].value)
    assert (pytest.approx(541.301, abs=1e-3) ==
            m.fs.HX1D.tube_outlet[0].vars["temperature"].value)
    assert (pytest.approx(2.65E7, abs=1e-3) ==
            m.fs.HX1D.tube_outlet[0].vars["pressure"].value)
