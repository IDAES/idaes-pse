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
Tests for Heat Exchanger 1D unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel, TerminationCondition,
                           SolverStatus, value)
from pyomo.common.config import ConfigBlock

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType, useDefault)
from idaes.generic_models.unit_models.heat_exchanger_1D import HeatExchanger1D as HX1D
from idaes.generic_models.unit_models.heat_exchanger_1D import WallConductionType
from idaes.generic_models.unit_models.heat_exchanger import HeatExchangerFlowPattern

from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.generic_models.properties import iapws95
from idaes.generic_models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock)

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from idaes.core.util.testing import (get_default_solver,
                                     PhysicalParameterTestBlock,
                                     initialization_tester)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = HX1D(default={
            "shell_side": {"property_package": m.fs.properties},
            "tube_side": {"property_package": m.fs.properties}})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 8
    assert isinstance(m.fs.unit.config.shell_side, ConfigBlock)
    assert isinstance(m.fs.unit.config.tube_side, ConfigBlock)
    assert m.fs.unit.config.flow_type == HeatExchangerFlowPattern.cocurrent
    assert m.fs.unit.config.has_wall_conduction == \
        WallConductionType.zero_dimensional
    assert m.fs.unit.config.finite_elements == 20
    assert m.fs.unit.config.collocation_points == 5

    # Check shell side config arguments
    assert len(m.fs.unit.config.shell_side) == 11
    assert m.fs.unit.config.shell_side.dynamic == useDefault
    assert m.fs.unit.config.shell_side.has_holdup == useDefault
    assert m.fs.unit.config.shell_side.material_balance_type == \
        MaterialBalanceType.useDefault
    assert m.fs.unit.config.shell_side.energy_balance_type == \
        EnergyBalanceType.useDefault
    assert m.fs.unit.config.shell_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.shell_side.has_pressure_change
    assert not m.fs.unit.config.shell_side.has_phase_equilibrium
    assert m.fs.unit.config.shell_side.transformation_method == \
        'dae.finite_difference'
    assert m.fs.unit.config.shell_side.transformation_scheme == 'BACKWARD'

    # Check tube side config arguments
    assert len(m.fs.unit.config.tube_side) == 11
    assert m.fs.unit.config.tube_side.dynamic == useDefault
    assert m.fs.unit.config.tube_side.has_holdup == useDefault
    assert m.fs.unit.config.tube_side.material_balance_type == \
        MaterialBalanceType.useDefault
    assert m.fs.unit.config.tube_side.energy_balance_type == \
        EnergyBalanceType.useDefault
    assert m.fs.unit.config.tube_side.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.tube_side.has_pressure_change
    assert not m.fs.unit.config.tube_side.has_phase_equilibrium
    assert m.fs.unit.config.tube_side.transformation_method == \
        'dae.finite_difference'
    assert m.fs.unit.config.tube_side.transformation_scheme == 'BACKWARD'


@pytest.mark.unit
def test_config_validation():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = BTXParameterBlock(default={"valid_phase": 'Liq'})

    with pytest.raises(ConfigurationError):
        m.fs.HX_co_current = HX1D(
            default={"shell_side": {"property_package": m.fs.properties,
                                    "transformation_scheme": "BACKWARD"},
                     "tube_side": {"property_package": m.fs.properties,
                                   "transformation_scheme": "FORWARD"},
                     "flow_type": HeatExchangerFlowPattern.cocurrent})

    with pytest.raises(ConfigurationError):
        m.fs.HX_counter_current = HX1D(
            default={"shell_side": {"property_package": m.fs.properties,
                                    "transformation_method":
                                    "dae.finite_difference"},
                     "tube_side": {"property_package": m.fs.properties,
                                   "transformation_method":
                                   "dae.collocation"},
                     "flow_type": HeatExchangerFlowPattern.countercurrent})


# -----------------------------------------------------------------------------
class TestBTX_cocurrent(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = BTXParameterBlock(default={"valid_phase": 'Liq'})

        m.fs.unit = HX1D(default={
                "shell_side": {"property_package": m.fs.properties},
                "tube_side": {"property_package": m.fs.properties},
                "flow_type": HeatExchangerFlowPattern.cocurrent})

        m.fs.unit.d_shell.fix(1.04)
        m.fs.unit.d_tube_outer.fix(0.01167)
        m.fs.unit.d_tube_inner.fix(0.01067)
        m.fs.unit.N_tubes.fix(10)
        m.fs.unit.shell_length.fix(4.85)
        m.fs.unit.tube_length.fix(4.85)
        m.fs.unit.shell_heat_transfer_coefficient.fix(2000)
        m.fs.unit.tube_heat_transfer_coefficient.fix(51000)

        m.fs.unit.shell_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.shell_inlet.temperature[0].fix(365)  # K
        m.fs.unit.shell_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.shell_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.shell_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.tube_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.tube_inlet.temperature[0].fix(300)  # K
        m.fs.unit.tube_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.tube_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.tube_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        return m

    @pytest.mark.unit
    @pytest.mark.build
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "shell_inlet")
        assert len(btx.fs.unit.shell_inlet.vars) == 4
        assert hasattr(btx.fs.unit.shell_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.shell_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.shell_inlet, "temperature")
        assert hasattr(btx.fs.unit.shell_inlet, "pressure")

        assert hasattr(btx.fs.unit, "tube_inlet")
        assert len(btx.fs.unit.tube_inlet.vars) == 4
        assert hasattr(btx.fs.unit.tube_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.tube_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.tube_inlet, "temperature")
        assert hasattr(btx.fs.unit.tube_inlet, "pressure")

        assert hasattr(btx.fs.unit, "shell_outlet")
        assert len(btx.fs.unit.shell_outlet.vars) == 4
        assert hasattr(btx.fs.unit.shell_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.shell_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.shell_outlet, "temperature")
        assert hasattr(btx.fs.unit.shell_outlet, "pressure")

        assert hasattr(btx.fs.unit, "tube_outlet")
        assert len(btx.fs.unit.tube_outlet.vars) == 4
        assert hasattr(btx.fs.unit.tube_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.tube_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.tube_outlet, "temperature")
        assert hasattr(btx.fs.unit.tube_outlet, "pressure")

        assert hasattr(btx.fs.unit, "shell_area")
        assert hasattr(btx.fs.unit, "shell_length")
        assert hasattr(btx.fs.unit, "tube_area")
        assert hasattr(btx.fs.unit, "tube_length")
        assert hasattr(btx.fs.unit, "pi")
        assert hasattr(btx.fs.unit, "d_shell")
        assert hasattr(btx.fs.unit, "d_tube_outer")
        assert hasattr(btx.fs.unit, "d_tube_inner")
        assert hasattr(btx.fs.unit, "N_tubes")
        assert hasattr(btx.fs.unit, "shell_heat_transfer_coefficient")
        assert hasattr(btx.fs.unit, "tube_heat_transfer_coefficient")
        assert hasattr(btx.fs.unit, "temperature_wall")
        assert hasattr(btx.fs.unit, "shell_heat_transfer_eq")
        assert hasattr(btx.fs.unit, "tube_heat_transfer_eq")
        assert hasattr(btx.fs.unit, "wall_0D_model")
        assert hasattr(btx.fs.unit, "area_calc_tube")
        assert hasattr(btx.fs.unit, "area_calc_shell")

        assert number_variables(btx) == 911
        assert number_total_constraints(btx) == 845
        assert number_unused_variables(btx) == 8

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.component
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert (pytest.approx(5, abs=1e-3) ==
                value(btx.fs.unit.shell_outlet.flow_mol[0]))
        assert (pytest.approx(322.669, abs=1e-3) ==
                value(btx.fs.unit.shell_outlet.temperature[0]))
        assert (pytest.approx(101325, abs=1e-3) ==
                value(btx.fs.unit.shell_outlet.pressure[0]))

        assert (pytest.approx(1, abs=1e-3) ==
                value(btx.fs.unit.tube_outlet.flow_mol[0]))
        assert (pytest.approx(322.463, abs=1e-3) ==
                value(btx.fs.unit.tube_outlet.temperature[0]))
        assert (pytest.approx(101325, abs=1e-3) ==
                value(btx.fs.unit.tube_outlet.pressure[0]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert abs(value(btx.fs.unit.shell_inlet.flow_mol[0] -
                         btx.fs.unit.shell_outlet.flow_mol[0])) <= 1e-6
        assert abs(value(btx.fs.unit.tube_inlet.flow_mol[0] -
                         btx.fs.unit.tube_outlet.flow_mol[0])) <= 1e-6

        shell_side = value(
                btx.fs.unit.shell_outlet.flow_mol[0] *
                (btx.fs.unit.shell.properties[0, 0].enth_mol_phase['Liq'] -
                 btx.fs.unit.shell.properties[0, 1].enth_mol_phase['Liq']))
        tube_side = value(
                btx.fs.unit.tube_outlet.flow_mol[0]*btx.fs.unit.N_tubes *
                (btx.fs.unit.tube.properties[0, 1].enth_mol_phase['Liq'] -
                 btx.fs.unit.tube.properties[0, 0].enth_mol_phase['Liq']))
        assert abs(shell_side - tube_side) <= 1e-6

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, btx):
        btx.fs.unit.report()


# -----------------------------------------------------------------------------
class TestBTX_countercurrent(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = BTXParameterBlock(default={"valid_phase": 'Liq'})

        m.fs.unit = HX1D(default={
                "shell_side": {"property_package": m.fs.properties},
                "tube_side": {"property_package": m.fs.properties},
                "flow_type": HeatExchangerFlowPattern.countercurrent})

        m.fs.unit.d_shell.fix(1.04)
        m.fs.unit.d_tube_outer.fix(0.01167)
        m.fs.unit.d_tube_inner.fix(0.01067)
        m.fs.unit.N_tubes.fix(10)
        m.fs.unit.shell_length.fix(4.85)
        m.fs.unit.tube_length.fix(4.85)
        m.fs.unit.shell_heat_transfer_coefficient.fix(2000)
        m.fs.unit.tube_heat_transfer_coefficient.fix(51000)

        m.fs.unit.shell_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.shell_inlet.temperature[0].fix(365)  # K
        m.fs.unit.shell_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.shell_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.shell_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.tube_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.tube_inlet.temperature[0].fix(300)  # K
        m.fs.unit.tube_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.tube_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.tube_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        return m

    @pytest.mark.unit
    @pytest.mark.build
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "shell_inlet")
        assert len(btx.fs.unit.shell_inlet.vars) == 4
        assert hasattr(btx.fs.unit.shell_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.shell_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.shell_inlet, "temperature")
        assert hasattr(btx.fs.unit.shell_inlet, "pressure")

        assert hasattr(btx.fs.unit, "tube_inlet")
        assert len(btx.fs.unit.tube_inlet.vars) == 4
        assert hasattr(btx.fs.unit.tube_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.tube_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.tube_inlet, "temperature")
        assert hasattr(btx.fs.unit.tube_inlet, "pressure")

        assert hasattr(btx.fs.unit, "shell_outlet")
        assert len(btx.fs.unit.shell_outlet.vars) == 4
        assert hasattr(btx.fs.unit.shell_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.shell_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.shell_outlet, "temperature")
        assert hasattr(btx.fs.unit.shell_outlet, "pressure")

        assert hasattr(btx.fs.unit, "tube_outlet")
        assert len(btx.fs.unit.tube_outlet.vars) == 4
        assert hasattr(btx.fs.unit.tube_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.tube_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.tube_outlet, "temperature")
        assert hasattr(btx.fs.unit.tube_outlet, "pressure")

        assert hasattr(btx.fs.unit, "shell_area")
        assert hasattr(btx.fs.unit, "shell_length")
        assert hasattr(btx.fs.unit, "tube_area")
        assert hasattr(btx.fs.unit, "tube_length")
        assert hasattr(btx.fs.unit, "pi")
        assert hasattr(btx.fs.unit, "d_shell")
        assert hasattr(btx.fs.unit, "d_tube_outer")
        assert hasattr(btx.fs.unit, "d_tube_inner")
        assert hasattr(btx.fs.unit, "N_tubes")
        assert hasattr(btx.fs.unit, "shell_heat_transfer_coefficient")
        assert hasattr(btx.fs.unit, "tube_heat_transfer_coefficient")
        assert hasattr(btx.fs.unit, "temperature_wall")
        assert hasattr(btx.fs.unit, "shell_heat_transfer_eq")
        assert hasattr(btx.fs.unit, "tube_heat_transfer_eq")
        assert hasattr(btx.fs.unit, "wall_0D_model")
        assert hasattr(btx.fs.unit, "area_calc_tube")
        assert hasattr(btx.fs.unit, "area_calc_shell")

        assert number_variables(btx) == 911
        assert number_total_constraints(btx) == 845
        assert number_unused_variables(btx) == 8

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx):
        initialization_tester(
                btx,
                optarg={'tol': 1e-6},
                shell_state_args={"flow_mol": 5,
                                  "temperature": 304,
                                  "pressure": 101325},
                tube_state_args={"flow_mol": 1,
                                 "temperature": 331.5,
                                 "pressure": 101325})

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert (pytest.approx(5, abs=1e-3) ==
                value(btx.fs.unit.shell_outlet.flow_mol[0]))
        assert (pytest.approx(304.292, abs=1e-3) ==
                value(btx.fs.unit.shell_outlet.temperature[0]))
        assert (pytest.approx(101325, abs=1e-3) ==
                value(btx.fs.unit.shell_outlet.pressure[0]))

        assert (pytest.approx(1, abs=1e-3) ==
                value(btx.fs.unit.tube_outlet.flow_mol[0]))
        assert (pytest.approx(331.435, abs=1e-3) ==
                value(btx.fs.unit.tube_outlet.temperature[0]))
        assert (pytest.approx(101325, abs=1e-3) ==
                value(btx.fs.unit.tube_outlet.pressure[0]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert abs(value(btx.fs.unit.shell_inlet.flow_mol[0] -
                         btx.fs.unit.shell_outlet.flow_mol[0])) <= 1e-6
        assert abs(value(btx.fs.unit.tube_inlet.flow_mol[0] -
                         btx.fs.unit.tube_outlet.flow_mol[0])) <= 1e-6

        shell_side = value(
                btx.fs.unit.shell_outlet.flow_mol[0] *
                (btx.fs.unit.shell.properties[0, 0].enth_mol_phase['Liq'] -
                 btx.fs.unit.shell.properties[0, 1].enth_mol_phase['Liq']))
        tube_side = value(
                btx.fs.unit.tube_outlet.flow_mol[0]*btx.fs.unit.N_tubes *
                (btx.fs.unit.tube.properties[0, 0].enth_mol_phase['Liq'] -
                 btx.fs.unit.tube.properties[0, 1].enth_mol_phase['Liq']))
        assert abs(shell_side - tube_side) <= 1e-6

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, btx):
        btx.fs.unit.report()


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(),
                    reason="IAPWS not available")
class TestIAPWS_cocurrent(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = iapws95.Iapws95ParameterBlock(default={
                "phase_presentation": iapws95.PhaseType.LG})

        m.fs.unit = HX1D(default={
                "shell_side": {"property_package": m.fs.properties},
                "tube_side": {"property_package": m.fs.properties},
                "flow_type": HeatExchangerFlowPattern.cocurrent})

        m.fs.unit.d_shell.fix(1.04)
        m.fs.unit.d_tube_outer.fix(0.01167)
        m.fs.unit.d_tube_inner.fix(0.01067)
        m.fs.unit.N_tubes.fix(10)
        m.fs.unit.shell_length.fix(4.85)
        m.fs.unit.tube_length.fix(4.85)
        m.fs.unit.shell_heat_transfer_coefficient.fix(2000)
        m.fs.unit.tube_heat_transfer_coefficient.fix(51000)

        m.fs.unit.shell_inlet.flow_mol[0].fix(5)
        m.fs.unit.shell_inlet.enth_mol[0].fix(50000)
        m.fs.unit.shell_inlet.pressure[0].fix(101325)

        m.fs.unit.tube_inlet.flow_mol[0].fix(5)
        m.fs.unit.tube_inlet.enth_mol[0].fix(7000)
        m.fs.unit.tube_inlet.pressure[0].fix(101325)

        return m

    @pytest.mark.unit
    @pytest.mark.build
    def test_build(self, iapws):
        assert len(iapws.fs.unit.shell_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.shell_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.shell_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.shell_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "shell_outlet")
        assert len(iapws.fs.unit.shell_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.shell_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.shell_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.shell_outlet, "pressure")

        assert len(iapws.fs.unit.tube_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.tube_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.tube_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.tube_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "tube_outlet")
        assert len(iapws.fs.unit.tube_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.tube_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.tube_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.tube_outlet, "pressure")

        assert hasattr(iapws.fs.unit, "shell_area")
        assert hasattr(iapws.fs.unit, "shell_length")
        assert hasattr(iapws.fs.unit, "tube_area")
        assert hasattr(iapws.fs.unit, "tube_length")
        assert hasattr(iapws.fs.unit, "pi")
        assert hasattr(iapws.fs.unit, "d_shell")
        assert hasattr(iapws.fs.unit, "d_tube_outer")
        assert hasattr(iapws.fs.unit, "d_tube_inner")
        assert hasattr(iapws.fs.unit, "N_tubes")
        assert hasattr(iapws.fs.unit, "shell_heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "tube_heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "temperature_wall")
        assert hasattr(iapws.fs.unit, "shell_heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "tube_heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "wall_0D_model")
        assert hasattr(iapws.fs.unit, "area_calc_tube")
        assert hasattr(iapws.fs.unit, "area_calc_shell")

        assert number_variables(iapws) == 659
        assert number_total_constraints(iapws) == 595
        assert number_unused_variables(iapws) == 10

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.initialization
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.unit
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(5, abs=1e-5) == \
            value(iapws.fs.unit.shell_outlet.flow_mol[0])
        assert pytest.approx(5, abs=1e-5) == \
            value(iapws.fs.unit.tube_outlet.flow_mol[0])

        assert pytest.approx(46298, abs=4e0) == \
            value(iapws.fs.unit.shell_outlet.enth_mol[0])
        assert pytest.approx(7370, abs=1e0) == \
            value(iapws.fs.unit.tube_outlet.enth_mol[0])

        assert pytest.approx(101325, abs=1e2) == \
            value(iapws.fs.unit.shell_outlet.pressure[0])
        assert pytest.approx(101325, abs=1e2) == \
            value(iapws.fs.unit.tube_outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert abs(value(iapws.fs.unit.shell_inlet.flow_mol[0] -
                         iapws.fs.unit.shell_outlet.flow_mol[0])) <= 1e-6
        assert abs(value(iapws.fs.unit.tube_inlet.flow_mol[0] -
                         iapws.fs.unit.tube_outlet.flow_mol[0])) <= 1e-6

        shell_side = value(
                iapws.fs.unit.shell_outlet.flow_mol[0] *
                (iapws.fs.unit.shell_inlet.enth_mol[0] -
                 iapws.fs.unit.shell_outlet.enth_mol[0]))
        tube_side = value(
                iapws.fs.unit.tube_outlet.flow_mol[0]*iapws.fs.unit.N_tubes *
                (iapws.fs.unit.tube_inlet.enth_mol[0] -
                 iapws.fs.unit.tube_outlet.enth_mol[0]))
        assert abs(shell_side + tube_side) <= 1e-6

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iapws):
        iapws.fs.unit.report()


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(),
                    reason="IAPWS not available")
class TestIAPWS_countercurrent(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = iapws95.Iapws95ParameterBlock(default={
                "phase_presentation": iapws95.PhaseType.LG})

        m.fs.unit = HX1D(default={
                "shell_side": {"property_package": m.fs.properties},
                "tube_side": {"property_package": m.fs.properties},
                "flow_type": HeatExchangerFlowPattern.countercurrent})

        m.fs.unit.d_shell.fix(1.04)
        m.fs.unit.d_tube_outer.fix(0.01167)
        m.fs.unit.d_tube_inner.fix(0.01067)
        m.fs.unit.N_tubes.fix(10)
        m.fs.unit.shell_length.fix(4.85)
        m.fs.unit.tube_length.fix(4.85)
        m.fs.unit.shell_heat_transfer_coefficient.fix(2000)
        m.fs.unit.tube_heat_transfer_coefficient.fix(51000)

        m.fs.unit.shell_inlet.flow_mol[0].fix(5)
        m.fs.unit.shell_inlet.enth_mol[0].fix(50000)
        m.fs.unit.shell_inlet.pressure[0].fix(101325)

        m.fs.unit.tube_inlet.flow_mol[0].fix(5)
        m.fs.unit.tube_inlet.enth_mol[0].fix(7000)
        m.fs.unit.tube_inlet.pressure[0].fix(101325)

        return m

    @pytest.mark.unit
    @pytest.mark.build
    def test_build(self, iapws):
        assert len(iapws.fs.unit.shell_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.shell_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.shell_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.shell_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "shell_outlet")
        assert len(iapws.fs.unit.shell_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.shell_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.shell_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.shell_outlet, "pressure")

        assert len(iapws.fs.unit.tube_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.tube_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.tube_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.tube_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "tube_outlet")
        assert len(iapws.fs.unit.tube_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.tube_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.tube_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.tube_outlet, "pressure")

        assert hasattr(iapws.fs.unit, "shell_area")
        assert hasattr(iapws.fs.unit, "shell_length")
        assert hasattr(iapws.fs.unit, "tube_area")
        assert hasattr(iapws.fs.unit, "tube_length")
        assert hasattr(iapws.fs.unit, "pi")
        assert hasattr(iapws.fs.unit, "d_shell")
        assert hasattr(iapws.fs.unit, "d_tube_outer")
        assert hasattr(iapws.fs.unit, "d_tube_inner")
        assert hasattr(iapws.fs.unit, "N_tubes")
        assert hasattr(iapws.fs.unit, "shell_heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "tube_heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "temperature_wall")
        assert hasattr(iapws.fs.unit, "shell_heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "tube_heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "wall_0D_model")
        assert hasattr(iapws.fs.unit, "area_calc_tube")
        assert hasattr(iapws.fs.unit, "area_calc_shell")

        assert number_variables(iapws) == 659
        assert number_total_constraints(iapws) == 595
        assert number_unused_variables(iapws) == 10

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.initialization
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(5, abs=1e-5) == \
            value(iapws.fs.unit.shell_outlet.flow_mol[0])
        assert pytest.approx(5, abs=1e-5) == \
            value(iapws.fs.unit.tube_outlet.flow_mol[0])

        assert pytest.approx(45359, abs=1e0) == \
            value(iapws.fs.unit.shell_outlet.enth_mol[0])
        assert pytest.approx(7464, abs=1e0) == \
            value(iapws.fs.unit.tube_outlet.enth_mol[0])

        assert pytest.approx(101325, abs=1e2) == \
            value(iapws.fs.unit.shell_outlet.pressure[0])
        assert pytest.approx(101325, abs=1e2) == \
            value(iapws.fs.unit.tube_outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert abs(value(iapws.fs.unit.shell_inlet.flow_mol[0] -
                         iapws.fs.unit.shell_outlet.flow_mol[0])) <= 1e-6
        assert abs(value(iapws.fs.unit.tube_inlet.flow_mol[0] -
                         iapws.fs.unit.tube_outlet.flow_mol[0])) <= 1e-6

        shell_side = value(
                iapws.fs.unit.shell_outlet.flow_mol[0] *
                (iapws.fs.unit.shell_inlet.enth_mol[0] -
                 iapws.fs.unit.shell_outlet.enth_mol[0]))
        tube_side = value(
                iapws.fs.unit.tube_outlet.flow_mol[0]*iapws.fs.unit.N_tubes *
                (iapws.fs.unit.tube_inlet.enth_mol[0] -
                 iapws.fs.unit.tube_outlet.enth_mol[0]))
        assert abs(shell_side + tube_side) <= 1e-6

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iapws):
        iapws.fs.unit.report()


# -----------------------------------------------------------------------------
class TestSaponification_cocurrent(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = HX1D(default={
                "shell_side": {"property_package": m.fs.properties},
                "tube_side": {"property_package": m.fs.properties},
                "flow_type": HeatExchangerFlowPattern.cocurrent})

        m.fs.unit.d_shell.fix(1.04)
        m.fs.unit.d_tube_outer.fix(0.01167)
        m.fs.unit.d_tube_inner.fix(0.01067)
        m.fs.unit.N_tubes.fix(10)
        m.fs.unit.shell_length.fix(4.85)
        m.fs.unit.tube_length.fix(4.85)
        m.fs.unit.shell_heat_transfer_coefficient.fix(2000)
        m.fs.unit.tube_heat_transfer_coefficient.fix(51000)

        m.fs.unit.shell_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.shell_inlet.temperature[0].fix(320)
        m.fs.unit.shell_inlet.pressure[0].fix(101325)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.tube_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.tube_inlet.temperature[0].fix(300)
        m.fs.unit.tube_inlet.pressure[0].fix(101325)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):
        assert len(sapon.fs.unit.shell_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.shell_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.shell_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.shell_inlet, "temperature")
        assert hasattr(sapon.fs.unit.shell_inlet, "pressure")

        assert len(sapon.fs.unit.shell_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.shell_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.shell_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.shell_outlet, "temperature")
        assert hasattr(sapon.fs.unit.shell_outlet, "pressure")

        assert len(sapon.fs.unit.tube_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.tube_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.tube_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.tube_inlet, "temperature")
        assert hasattr(sapon.fs.unit.tube_inlet, "pressure")

        assert len(sapon.fs.unit.tube_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.tube_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.tube_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.tube_outlet, "temperature")
        assert hasattr(sapon.fs.unit.tube_outlet, "pressure")

        assert hasattr(sapon.fs.unit, "shell_area")
        assert hasattr(sapon.fs.unit, "shell_length")
        assert hasattr(sapon.fs.unit, "tube_area")
        assert hasattr(sapon.fs.unit, "tube_length")
        assert hasattr(sapon.fs.unit, "pi")
        assert hasattr(sapon.fs.unit, "d_shell")
        assert hasattr(sapon.fs.unit, "d_tube_outer")
        assert hasattr(sapon.fs.unit, "d_tube_inner")
        assert hasattr(sapon.fs.unit, "N_tubes")
        assert hasattr(sapon.fs.unit, "shell_heat_transfer_coefficient")
        assert hasattr(sapon.fs.unit, "tube_heat_transfer_coefficient")
        assert hasattr(sapon.fs.unit, "temperature_wall")
        assert hasattr(sapon.fs.unit, "shell_heat_transfer_eq")
        assert hasattr(sapon.fs.unit, "tube_heat_transfer_eq")
        assert hasattr(sapon.fs.unit, "wall_0D_model")
        assert hasattr(sapon.fs.unit, "area_calc_tube")
        assert hasattr(sapon.fs.unit, "area_calc_shell")

        assert number_variables(sapon) == 1037
        assert number_total_constraints(sapon) == 959
        assert number_unused_variables(sapon) == 14

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.initialization
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(1e-3, abs=1e-6) == \
            value(sapon.fs.unit.shell_outlet.flow_vol[0])
        assert pytest.approx(1e-3, abs=1e-6) == \
            value(sapon.fs.unit.tube_outlet.flow_vol[0])

        assert 55388.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "EthylAcetate"])
        assert 0.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "SodiumAcetate"])
        assert 0.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "Ethanol"])

        assert 55388.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "EthylAcetate"])
        assert 0.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "SodiumAcetate"])
        assert 0.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "Ethanol"])

        assert pytest.approx(309.4, abs=1e-1) == \
            value(sapon.fs.unit.shell_outlet.temperature[0])
        assert pytest.approx(301.1, abs=1e-1) == \
            value(sapon.fs.unit.tube_outlet.temperature[0])

        assert pytest.approx(101325, abs=1e2) == \
            value(sapon.fs.unit.shell_outlet.pressure[0])
        assert pytest.approx(101325, abs=1e2) == \
            value(sapon.fs.unit.tube_outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        shell_side = value(
                sapon.fs.unit.shell_outlet.flow_vol[0] *
                sapon.fs.properties.dens_mol*sapon.fs.properties.cp_mol *
                (sapon.fs.unit.shell_inlet.temperature[0] -
                 sapon.fs.unit.shell_outlet.temperature[0]))
        tube_side = value(
                sapon.fs.unit.tube_outlet.flow_vol[0]*sapon.fs.unit.N_tubes *
                sapon.fs.properties.dens_mol*sapon.fs.properties.cp_mol *
                (sapon.fs.unit.tube_inlet.temperature[0] -
                 sapon.fs.unit.tube_outlet.temperature[0]))
        assert abs(shell_side + tube_side) <= 1e-6

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, sapon):
        sapon.fs.unit.report()


# -----------------------------------------------------------------------------
class TestSaponification_countercurrent(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = HX1D(default={
                "shell_side": {"property_package": m.fs.properties},
                "tube_side": {"property_package": m.fs.properties},
                "flow_type": HeatExchangerFlowPattern.countercurrent})

        m.fs.unit.d_shell.fix(1.04)
        m.fs.unit.d_tube_outer.fix(0.01167)
        m.fs.unit.d_tube_inner.fix(0.01067)
        m.fs.unit.N_tubes.fix(10)
        m.fs.unit.shell_length.fix(4.85)
        m.fs.unit.tube_length.fix(4.85)
        m.fs.unit.shell_heat_transfer_coefficient.fix(2000)
        m.fs.unit.tube_heat_transfer_coefficient.fix(51000)

        m.fs.unit.shell_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.shell_inlet.temperature[0].fix(320)
        m.fs.unit.shell_inlet.pressure[0].fix(101325)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.shell_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.tube_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.tube_inlet.temperature[0].fix(300)
        m.fs.unit.tube_inlet.pressure[0].fix(101325)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.tube_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):
        assert len(sapon.fs.unit.shell_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.shell_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.shell_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.shell_inlet, "temperature")
        assert hasattr(sapon.fs.unit.shell_inlet, "pressure")

        assert len(sapon.fs.unit.shell_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.shell_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.shell_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.shell_outlet, "temperature")
        assert hasattr(sapon.fs.unit.shell_outlet, "pressure")

        assert len(sapon.fs.unit.tube_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.tube_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.tube_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.tube_inlet, "temperature")
        assert hasattr(sapon.fs.unit.tube_inlet, "pressure")

        assert len(sapon.fs.unit.tube_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.tube_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.tube_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.tube_outlet, "temperature")
        assert hasattr(sapon.fs.unit.tube_outlet, "pressure")

        assert hasattr(sapon.fs.unit, "shell_area")
        assert hasattr(sapon.fs.unit, "shell_length")
        assert hasattr(sapon.fs.unit, "tube_area")
        assert hasattr(sapon.fs.unit, "tube_length")
        assert hasattr(sapon.fs.unit, "pi")
        assert hasattr(sapon.fs.unit, "d_shell")
        assert hasattr(sapon.fs.unit, "d_tube_outer")
        assert hasattr(sapon.fs.unit, "d_tube_inner")
        assert hasattr(sapon.fs.unit, "N_tubes")
        assert hasattr(sapon.fs.unit, "shell_heat_transfer_coefficient")
        assert hasattr(sapon.fs.unit, "tube_heat_transfer_coefficient")
        assert hasattr(sapon.fs.unit, "temperature_wall")
        assert hasattr(sapon.fs.unit, "shell_heat_transfer_eq")
        assert hasattr(sapon.fs.unit, "tube_heat_transfer_eq")
        assert hasattr(sapon.fs.unit, "wall_0D_model")
        assert hasattr(sapon.fs.unit, "area_calc_tube")
        assert hasattr(sapon.fs.unit, "area_calc_shell")

        assert number_variables(sapon) == 1037
        assert number_total_constraints(sapon) == 959
        assert number_unused_variables(sapon) == 14

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(1e-3, abs=1e-6) == \
            value(sapon.fs.unit.shell_outlet.flow_vol[0])
        assert pytest.approx(1e-3, abs=1e-6) == \
            value(sapon.fs.unit.tube_outlet.flow_vol[0])

        assert 55388.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "EthylAcetate"])
        assert 0.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "SodiumAcetate"])
        assert 0.0 == value(
                sapon.fs.unit.shell_inlet.conc_mol_comp[0, "Ethanol"])

        assert 55388.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "EthylAcetate"])
        assert 0.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "SodiumAcetate"])
        assert 0.0 == value(
                sapon.fs.unit.tube_inlet.conc_mol_comp[0, "Ethanol"])

        assert pytest.approx(309.2, abs=1e-1) == \
            value(sapon.fs.unit.shell_outlet.temperature[0])
        assert pytest.approx(301.1, abs=1e-1) == \
            value(sapon.fs.unit.tube_outlet.temperature[0])

        assert pytest.approx(101325, abs=1e2) == \
            value(sapon.fs.unit.shell_outlet.pressure[0])
        assert pytest.approx(101325, abs=1e2) == \
            value(sapon.fs.unit.tube_outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        shell_side = value(
                sapon.fs.unit.shell_outlet.flow_vol[0] *
                sapon.fs.properties.dens_mol*sapon.fs.properties.cp_mol *
                (sapon.fs.unit.shell_inlet.temperature[0] -
                 sapon.fs.unit.shell_outlet.temperature[0]))
        tube_side = value(
                sapon.fs.unit.tube_outlet.flow_vol[0]*sapon.fs.unit.N_tubes *
                sapon.fs.properties.dens_mol*sapon.fs.properties.cp_mol *
                (sapon.fs.unit.tube_inlet.temperature[0] -
                 sapon.fs.unit.tube_outlet.temperature[0]))
        assert abs(shell_side + tube_side) <= 1e-6

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, sapon):
        sapon.fs.unit.report()
