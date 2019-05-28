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
from idaes.unit_models.heat_exchanger_1D import WallConductionType
from idaes.unit_models.heat_exchanger import HeatExchangerFlowPattern

from idaes.property_models.ideal.BTX_ideal_VLE import BTXParameterBlock
from idaes.ui.report import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError

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

m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                             'Liq'})

# Default options
m.fs.HX_co_current = HX1D(
    default={"shell_side": {"property_package": m.fs.properties},
             "tube_side": {"property_package": m.fs.properties},
             "flow_type": HeatExchangerFlowPattern.cocurrent})
# Default options
m.fs.HX_counter_current = HX1D(
    default={"shell_side": {"property_package": m.fs.properties},
             "tube_side": {"property_package": m.fs.properties},
             "flow_type": HeatExchangerFlowPattern.countercurrent})


def test_build():
    # Check build for co-current configuration
    assert len(m.fs.HX_co_current.config) == 8
    assert m.fs.HX_co_current.config.flow_type == \
        HeatExchangerFlowPattern.cocurrent
    assert m.fs.HX_co_current.config.has_wall_conduction == \
        WallConductionType.zero_dimensional

    assert len(m.fs.HX_co_current.config.shell_side) == 11
    assert not m.fs.HX_co_current.config.shell_side.has_holdup
    assert m.fs.HX_co_current.config.shell_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX_co_current.config.shell_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX_co_current.config.shell_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert not m.fs.HX_co_current.config.shell_side.has_pressure_change
    assert not m.fs.HX_co_current.config.shell_side.has_phase_equilibrium

    assert len(m.fs.HX_co_current.config.tube_side) == 11
    assert not m.fs.HX_co_current.config.tube_side.has_holdup
    assert m.fs.HX_co_current.config.tube_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX_co_current.config.tube_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX_co_current.config.tube_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert not m.fs.HX_co_current.config.tube_side.has_pressure_change
    assert not m.fs.HX_co_current.config.tube_side.has_phase_equilibrium

    # Check for inlets/outlets construction for co-current configuration
    assert hasattr(m.fs.HX_co_current, "shell_inlet")
    assert hasattr(m.fs.HX_co_current, "shell_outlet")
    assert hasattr(m.fs.HX_co_current, "tube_inlet")
    assert hasattr(m.fs.HX_co_current, "tube_outlet")

    # Check build for counter-current configuration
    assert len(m.fs.HX_counter_current.config) == 8
    assert m.fs.HX_counter_current.config.flow_type == \
        HeatExchangerFlowPattern.countercurrent
    assert m.fs.HX_counter_current.config.has_wall_conduction == \
        WallConductionType.zero_dimensional

    assert len(m.fs.HX_counter_current.config.shell_side) == 11
    assert not m.fs.HX_counter_current.config.shell_side.has_holdup
    assert m.fs.HX_counter_current.config.shell_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX_counter_current.config.shell_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX_counter_current.config.shell_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert not m.fs.HX_counter_current.config.shell_side.has_pressure_change
    assert not m.fs.HX_counter_current.config.shell_side.has_phase_equilibrium

    assert len(m.fs.HX_counter_current.config.tube_side) == 11
    assert not m.fs.HX_counter_current.config.tube_side.has_holdup
    assert m.fs.HX_counter_current.config.tube_side.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.HX_counter_current.config.tube_side.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.HX_counter_current.config.tube_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
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
    m.fs.HX_co_current.N_tubes.fix(10)
    m.fs.HX_co_current.shell_length.fix(4.85)
    m.fs.HX_co_current.tube_length.fix(4.85)
    m.fs.HX_co_current.shell_heat_transfer_coefficient.fix(2000)
    m.fs.HX_co_current.tube_heat_transfer_coefficient.fix(51000)

    # Boundary conditions
    for t in m.fs.time:
        # Unit HX01
        m.fs.HX_co_current.shell_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.HX_co_current.shell_inlet.temperature[0].fix(365)  # K
        m.fs.HX_co_current.shell_inlet.pressure[0].fix(101325)  # Pa
        m.fs.HX_co_current.shell_inlet.mole_frac[0, "benzene"].fix(0.5)
        m.fs.HX_co_current.shell_inlet.mole_frac[0, "toluene"].fix(0.5)

        m.fs.HX_co_current.tube_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.HX_co_current.tube_inlet.temperature[0].fix(300)  # K
        m.fs.HX_co_current.tube_inlet.pressure[0].fix(101325)  # Pa
        m.fs.HX_co_current.tube_inlet.mole_frac[0, "benzene"].fix(0.5)
        m.fs.HX_co_current.tube_inlet.mole_frac[0, "toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.HX_co_current) == 0

    """Set inlet and operating conditions for counter-current
       heat exchanger."""

    # HX dimensions
    m.fs.HX_counter_current.d_shell.fix(1.04)
    m.fs.HX_counter_current.d_tube_outer.fix(0.01167)
    m.fs.HX_counter_current.d_tube_inner.fix(0.01067)
    m.fs.HX_counter_current.N_tubes.fix(10)
    m.fs.HX_counter_current.shell_length.fix(4.85)
    m.fs.HX_counter_current.tube_length.fix(4.85)
    m.fs.HX_counter_current.shell_heat_transfer_coefficient.fix(2000)
    m.fs.HX_counter_current.tube_heat_transfer_coefficient.fix(51000)

    # Boundary conditions
    for t in m.fs.time:
        # Unit HX01
        m.fs.HX_counter_current.shell_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.HX_counter_current.shell_inlet.temperature[0].fix(365)  # K
        m.fs.HX_counter_current.shell_inlet.pressure[0].fix(101325)  # Pa
        m.fs.HX_counter_current.shell_inlet.mole_frac[0, "benzene"].fix(0.5)
        m.fs.HX_counter_current.shell_inlet.mole_frac[0, "toluene"].fix(0.5)

        m.fs.HX_counter_current.tube_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.HX_counter_current.tube_inlet.temperature[0].fix(300)  # K
        m.fs.HX_counter_current.tube_inlet.pressure[0].fix(101325)  # Pa
        m.fs.HX_counter_current.tube_inlet.mole_frac[0, "benzene"].fix(0.5)
        m.fs.HX_counter_current.tube_inlet.mole_frac[0, "toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.HX_counter_current) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialization():
    """Test initialize and solve for co-current heat exchanger."""
    m.fs.HX_co_current.initialize()
    results = solver.solve(m, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert (pytest.approx(5, abs=1e-3) ==
            m.fs.HX_co_current.shell_outlet.flow_mol[0].value)
    assert (pytest.approx(322.669, abs=1e-3) ==
            m.fs.HX_co_current.shell_outlet.temperature[0].value)
    assert (pytest.approx(101325, abs=1e-3) ==
            m.fs.HX_co_current.shell_outlet.pressure[0].value)

    assert (pytest.approx(1, abs=1e-3) ==
            m.fs.HX_co_current.tube_outlet.flow_mol[0].value)
    assert (pytest.approx(322.463, abs=1e-3) ==
            m.fs.HX_co_current.tube_outlet.temperature[0].value)
    assert (pytest.approx(101325, abs=1e-3) ==
            m.fs.HX_co_current.tube_outlet.pressure[0].value)

    # Check for energy conservation
    shell_side = 5 * (m.fs.HX_co_current.shell.properties[0, 0].
                      enth_mol_phase['Liq'].value - m.fs.HX_co_current.shell.
                      properties[0, 1].
                      enth_mol_phase['Liq'].value)
    tube_side = 1 * 10 * (m.fs.HX_co_current.tube.properties[0, 1].
                          enth_mol_phase['Liq'].value -
                          m.fs.HX_co_current.tube.properties[0, 0].
                          enth_mol_phase['Liq'].value)
    assert (shell_side - tube_side) <= 1e-6

    """Test initialize and solve for counter-current heat exchanger."""
    m.fs.HX_counter_current.initialize()
    results = solver.solve(m, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert (pytest.approx(5, abs=1e-3) ==
            m.fs.HX_counter_current.shell_outlet.flow_mol[0].value)
    assert (pytest.approx(304.292, abs=1e-3) ==
            m.fs.HX_counter_current.shell_outlet.temperature[0].value)
    assert (pytest.approx(101325, abs=1e-3) ==
            m.fs.HX_counter_current.shell_outlet.pressure[0].value)

    assert (pytest.approx(1, abs=1e-3) ==
            m.fs.HX_counter_current.tube_outlet.flow_mol[0].value)
    assert (pytest.approx(331.435, abs=1e-3) ==
            m.fs.HX_counter_current.tube_outlet.temperature[0].value)
    assert (pytest.approx(101325, abs=1e-3) ==
            m.fs.HX_counter_current.tube_outlet.pressure[0].value)

    # Check for energy conservation
    shell_side = 5 * (m.fs.HX_counter_current.shell.properties[0, 0].
                      enth_mol_phase['Liq'].value -
                      m.fs.HX_counter_current.shell.properties[0, 1].
                      enth_mol_phase['Liq'].value)
    tube_side = 1 * 10 * (m.fs.HX_counter_current.tube.properties[0, 0].
                          enth_mol_phase['Liq'].value -
                          m.fs.HX_counter_current.tube.properties[0, 1].
                          enth_mol_phase['Liq'].value)
    assert (shell_side - tube_side) <= 1e-6


# Test the custom discretisation options
m.fs1 = FlowsheetBlock(default={"dynamic": False})


def test_custom_build():
    with pytest.raises(ConfigurationError):
        m.fs1.HX_co_current = HX1D(
            default={"shell_side": {"property_package": m.fs.properties,
                                    "transformation_scheme": "BACKWARD"},
                     "tube_side": {"property_package": m.fs.properties,
                                   "transformation_scheme": "FORWARD"},
                     "flow_type": HeatExchangerFlowPattern.cocurrent})

    with pytest.raises(ConfigurationError):
        m.fs1.HX_counter_current = HX1D(
            default={"shell_side": {"property_package": m.fs.properties,
                                    "transformation_method":
                                    "dae.finite_difference"},
                     "tube_side": {"property_package": m.fs.properties,
                                   "transformation_method":
                                   "dae.collocation"},
                     "flow_type": HeatExchangerFlowPattern.countercurrent})
