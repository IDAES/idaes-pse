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

from idaes.property_models.examples.BFW_properties import BFWParameterBlock
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

m.fs.properties = BFWParameterBlock()

m.fs.HX_co_current = HX1D(
    default={"shell_side": {"property_package": m.fs.properties},
             "tube_side": {"property_package": m.fs.properties},
             "flow_type": "co_current"})

m.fs.HX_counter_current = HX1D(
    default={"shell_side": {"property_package": m.fs.properties},
             "tube_side": {"property_package": m.fs.properties},
             "flow_type": "counter_current"})


def test_build():
    # Check build for co-current configuration
    assert len(m.fs.HX_co_current.config) == 8
    assert m.fs.HX_co_current.config.flow_type == "co_current"
    assert m.fs.HX_co_current.config.has_wall_conduction == "none"

    assert len(m.fs.HX_co_current.config.shell_side) == 12
    assert not m.fs.HX_co_current.config.shell_side.has_holdup
    assert m.fs.HX_co_current.config.shell_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX_co_current.config.shell_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX_co_current.config.shell_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.HX_co_current.config.shell_side.has_heat_transfer
    assert not m.fs.HX_co_current.config.shell_side.has_pressure_change
    assert not m.fs.HX_co_current.config.shell_side.has_phase_equilibrium

    assert len(m.fs.HX_co_current.config.tube_side) == 12
    assert not m.fs.HX_co_current.config.tube_side.has_holdup
    assert m.fs.HX_co_current.config.tube_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX_co_current.config.tube_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX_co_current.config.tube_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.HX_co_current.config.tube_side.has_heat_transfer
    assert not m.fs.HX_co_current.config.tube_side.has_pressure_change
    assert not m.fs.HX_co_current.config.tube_side.has_phase_equilibrium

    # Check for inlets/outlets construction for co-current configuration
    assert hasattr(m.fs.HX_co_current, "shell_inlet")
    assert hasattr(m.fs.HX_co_current, "shell_outlet")
    assert hasattr(m.fs.HX_co_current, "tube_inlet")
    assert hasattr(m.fs.HX_co_current, "tube_outlet")

    # Check build for counter-current configuration
    assert len(m.fs.HX_counter_current.config) == 8
    assert m.fs.HX_counter_current.config.flow_type == "counter_current"
    assert m.fs.HX_counter_current.config.has_wall_conduction == "none"

    assert len(m.fs.HX_counter_current.config.shell_side) == 12
    assert not m.fs.HX_counter_current.config.shell_side.has_holdup
    assert m.fs.HX_counter_current.config.shell_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX_counter_current.config.shell_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX_counter_current.config.shell_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.HX_counter_current.config.shell_side.has_heat_transfer
    assert not m.fs.HX_counter_current.config.shell_side.has_pressure_change
    assert not m.fs.HX_counter_current.config.shell_side.has_phase_equilibrium

    assert len(m.fs.HX_counter_current.config.tube_side) == 12
    assert not m.fs.HX_counter_current.config.tube_side.has_holdup
    assert m.fs.HX_counter_current.config.tube_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX_counter_current.config.tube_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX_counter_current.config.tube_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.HX_counter_current.config.tube_side.has_heat_transfer
    assert not m.fs.HX_counter_current.config.tube_side.has_pressure_change
    assert not m.fs.HX_counter_current.config.tube_side.has_phase_equilibrium

    # Check for inlets/outlets construction for counter_current configuration
    assert hasattr(m.fs.HX_counter_current, "shell_inlet")
    assert hasattr(m.fs.HX_counter_current, "shell_outlet")
    assert hasattr(m.fs.HX_counter_current, "tube_inlet")
    assert hasattr(m.fs.HX_counter_current, "tube_outlet")


def test_setInputs():
    """Set inlet and operating conditions for co-current heat exchanger."""
    # HX dimensions
    m.fs.HX_co_current.d_shell.fix(1.04)
    m.fs.HX_co_current.d_tube_outer.fix(0.01167)
    m.fs.HX_co_current.d_tube_inner.fix(0.01067)
    m.fs.HX_co_current.N_tubes.fix(1176)
    m.fs.HX_co_current.shell_length.fix(4.85)
    m.fs.HX_co_current.tube_length.fix(4.85)
    m.fs.HX_co_current.shell_heat_transfer_coefficient.fix(2000)
    m.fs.HX_co_current.tube_heat_transfer_coefficient.fix(51000)

    # Boundary conditions
    for t in m.fs.time:
        # Unit HX01
        # NETL baseline study
        m.fs.HX_co_current.shell_inlet.flow_mol[0].fix(2300)  # mol/s
        m.fs.HX_co_current.shell_inlet.temperature[0].fix(676)  # K
        m.fs.HX_co_current.shell_inlet.pressure[0].fix(7.38E6)  # Pa
        m.fs.HX_co_current.shell_inlet.vapor_frac[0].fix(1)

        m.fs.HX_co_current.tube_inlet.flow_mol[0].fix(26.6)  # mol/s
        m.fs.HX_co_current.tube_inlet.temperature[0].fix(529)  # K
        m.fs.HX_co_current.tube_inlet.pressure[0].fix(2.65E7)  # Pa
        m.fs.HX_co_current.tube_inlet.vapor_frac[0].fix(0)

    assert degrees_of_freedom(m.fs.HX_co_current) == 0

    """Set inlet and operating conditions for counter-current
       heat exchanger."""

    # HX dimensions
    m.fs.HX_counter_current.d_shell.fix(1.04)
    m.fs.HX_counter_current.d_tube_outer.fix(0.01167)
    m.fs.HX_counter_current.d_tube_inner.fix(0.01067)
    m.fs.HX_counter_current.N_tubes.fix(1176)
    m.fs.HX_counter_current.shell_length.fix(4.85)
    m.fs.HX_counter_current.tube_length.fix(4.85)
    m.fs.HX_counter_current.shell_heat_transfer_coefficient.fix(2000)
    m.fs.HX_counter_current.tube_heat_transfer_coefficient.fix(51000)

    # Boundary conditions
    for t in m.fs.time:
        # Unit HX01
        # NETL baseline study
        m.fs.HX_counter_current.shell_inlet.flow_mol[0].fix(2300)  # mol/s
        m.fs.HX_counter_current.shell_inlet.temperature[0].fix(676)  # K
        m.fs.HX_counter_current.shell_inlet.pressure[0].fix(7.38E6)  # Pa
        m.fs.HX_counter_current.shell_inlet.vapor_frac[0].fix(1)

        m.fs.HX_counter_current.tube_inlet.flow_mol[0].fix(26.6)  # mol/s
        m.fs.HX_counter_current.tube_inlet.temperature[0].fix(529)  # K
        m.fs.HX_counter_current.tube_inlet.pressure[0].fix(2.65E7)  # Pa
        m.fs.HX_counter_current.tube_inlet.vapor_frac[0].fix(0)

    assert degrees_of_freedom(m.fs.HX_counter_current) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization():
    """Test initialize and solve for co-current heat exchanger."""
    m.fs.HX_co_current.initialize()
    results = solver.solve(m, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert (pytest.approx(2300, abs=1e-3) ==
            m.fs.HX_co_current.shell_outlet.flow_mol[0].value)
    assert (pytest.approx(559.220, abs=1e-3) ==
            m.fs.HX_co_current.shell_outlet.temperature[0].value)
    assert (pytest.approx(7.38E6, abs=1e-3) ==
            m.fs.HX_co_current.shell_outlet.pressure[0].value)

    assert (pytest.approx(26.6, abs=1e-3) ==
            m.fs.HX_co_current.tube_outlet.flow_mol[0].value)
    assert (pytest.approx(541.301, abs=1e-3) ==
            m.fs.HX_co_current.tube_outlet.temperature[0].value)
    assert (pytest.approx(2.65E7, abs=1e-3) ==
            m.fs.HX_co_current.tube_outlet.pressure[0].value)

    """Test initialize and solve for counter-current heat exchanger."""
    m.fs.HX_counter_current.initialize()
    results = solver.solve(m, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert (pytest.approx(2300, abs=1e-3) ==
            m.fs.HX_counter_current.shell_outlet.flow_mol[0].value)
    assert (pytest.approx(552.311, abs=1e-3) ==
            m.fs.HX_counter_current.shell_outlet.temperature[0].value)
    assert (pytest.approx(7.38E6, abs=1e-3) ==
            m.fs.HX_counter_current.shell_outlet.pressure[0].value)

    assert (pytest.approx(26.6, abs=1e-3) ==
            m.fs.HX_counter_current.tube_outlet.flow_mol[0].value)
    assert (pytest.approx(542.869, abs=1e-3) ==
            m.fs.HX_counter_current.tube_outlet.temperature[0].value)
    assert (pytest.approx(2.65E7, abs=1e-3) ==
            m.fs.HX_counter_current.tube_outlet.pressure[0].value)
