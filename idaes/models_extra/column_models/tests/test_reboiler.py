#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for Reboiler unit model. Tests 2 sets of state vars using the
ideal property package (FTPz and FcTP).

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import check_optimal_termination, ConcreteModel, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models_extra.column_models import Reboiler
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@pytest.mark.unit
def test_config():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = Reboiler(property_package=m.fs.properties)

    assert len(m.fs.unit.config) == 9

    assert not m.fs.unit.config.has_boilup_ratio
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert hasattr(m.fs.unit, "heat_duty")


class TestBTXIdeal:
    @pytest.fixture(scope="class")
    def btx_ftpz(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal"
        )

        m.fs.unit = Reboiler(property_package=m.fs.properties, has_boilup_ratio=True)

        # Fix the reboiler variables
        m.fs.unit.boilup_ratio.fix(1)

        # Fix the inputs (typically this will be the outlet liquid from the
        # bottom tray)
        m.fs.unit.inlet.flow_mol.fix(1)
        m.fs.unit.inlet.temperature.fix(363)
        m.fs.unit.inlet.pressure.fix(101325)
        m.fs.unit.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        return m

    @pytest.fixture(scope="class")
    def btx_fctp(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FcTP"
        )

        m.fs.unit = Reboiler(property_package=m.fs.properties, has_boilup_ratio=True)

        # Fix the reboiler variables
        m.fs.unit.boilup_ratio.fix(1)

        # Fix the inputs (typically this will be the outlet liquid from the
        # bottom tray)
        m.fs.unit.inlet.flow_mol_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.flow_mol_comp[0, "toluene"].fix(0.5)
        m.fs.unit.inlet.temperature.fix(363)
        m.fs.unit.inlet.pressure.fix(101325)

        return m

    @pytest.mark.unit
    def test_build(self, btx_ftpz, btx_fctp):

        assert hasattr(btx_ftpz.fs.unit, "boilup_ratio")
        assert hasattr(btx_ftpz.fs.unit, "eq_boilup_ratio")

        assert hasattr(btx_ftpz.fs.unit, "inlet")
        assert hasattr(btx_ftpz.fs.unit.inlet, "flow_mol")
        assert hasattr(btx_ftpz.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(btx_ftpz.fs.unit.inlet, "temperature")
        assert hasattr(btx_ftpz.fs.unit.inlet, "pressure")

        assert hasattr(btx_ftpz.fs.unit, "bottoms")
        assert hasattr(btx_ftpz.fs.unit.bottoms, "flow_mol")
        assert hasattr(btx_ftpz.fs.unit.bottoms, "mole_frac_comp")
        assert hasattr(btx_ftpz.fs.unit.bottoms, "temperature")
        assert hasattr(btx_ftpz.fs.unit.bottoms, "pressure")

        assert hasattr(btx_ftpz.fs.unit, "vapor_reboil")
        assert hasattr(btx_ftpz.fs.unit.vapor_reboil, "flow_mol")
        assert hasattr(btx_ftpz.fs.unit.vapor_reboil, "mole_frac_comp")
        assert hasattr(btx_ftpz.fs.unit.vapor_reboil, "temperature")
        assert hasattr(btx_ftpz.fs.unit.vapor_reboil, "pressure")

        assert number_variables(btx_ftpz.fs.unit) == 48
        assert number_total_constraints(btx_ftpz.fs.unit) == 42
        assert number_unused_variables(btx_ftpz) == 0

        assert hasattr(btx_fctp.fs.unit, "boilup_ratio")
        assert hasattr(btx_fctp.fs.unit, "eq_boilup_ratio")

        assert hasattr(btx_fctp.fs.unit, "inlet")
        assert hasattr(btx_fctp.fs.unit.inlet, "flow_mol_comp")
        assert hasattr(btx_fctp.fs.unit.inlet, "temperature")
        assert hasattr(btx_fctp.fs.unit.inlet, "pressure")

        assert hasattr(btx_fctp.fs.unit, "bottoms")
        assert hasattr(btx_fctp.fs.unit.bottoms, "flow_mol_comp")
        assert hasattr(btx_fctp.fs.unit.bottoms, "temperature")
        assert hasattr(btx_fctp.fs.unit.bottoms, "pressure")

        assert hasattr(btx_fctp.fs.unit, "vapor_reboil")
        assert hasattr(btx_fctp.fs.unit.vapor_reboil, "flow_mol_comp")
        assert hasattr(btx_fctp.fs.unit.vapor_reboil, "temperature")
        assert hasattr(btx_fctp.fs.unit.vapor_reboil, "pressure")

        assert number_variables(btx_fctp.fs.unit) == 50
        assert number_total_constraints(btx_fctp.fs.unit) == 45
        assert number_unused_variables(btx_fctp) == 0

    @pytest.mark.unit
    def test_dof(self, btx_ftpz, btx_fctp):
        assert degrees_of_freedom(btx_ftpz) == 0
        assert degrees_of_freedom(btx_fctp) == 0

    @pytest.mark.component
    def test_units_FTPz(self, btx_ftpz, btx_fctp):
        assert_units_consistent(btx_ftpz)

    @pytest.mark.component
    def test_units_FcTP(self, btx_fctp):
        assert_units_consistent(btx_fctp)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx_ftpz, btx_fctp):
        initialization_tester(btx_ftpz)
        initialization_tester(btx_fctp)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx_ftpz, btx_fctp):
        results = solver.solve(btx_ftpz)

        # Check for optimal solution
        assert check_optimal_termination(results)

        results = solver.solve(btx_fctp)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx_ftpz, btx_fctp):
        # Bottoms port
        assert pytest.approx(0.5, abs=1e-3) == value(
            btx_ftpz.fs.unit.bottoms.flow_mol[0]
        )
        assert pytest.approx(0.3891, abs=1e-3) == value(
            btx_ftpz.fs.unit.bottoms.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.6109, abs=1e-3) == value(
            btx_ftpz.fs.unit.bottoms.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(368.728, abs=1e-2) == value(
            btx_ftpz.fs.unit.bottoms.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_ftpz.fs.unit.bottoms.pressure[0]
        )

        # Vapor reboil port
        assert pytest.approx(0.5, abs=1e-3) == value(
            btx_ftpz.fs.unit.vapor_reboil.flow_mol[0]
        )
        assert pytest.approx(0.6108, abs=1e-3) == value(
            btx_ftpz.fs.unit.vapor_reboil.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.3892, abs=1e-3) == value(
            btx_ftpz.fs.unit.vapor_reboil.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(368.728, abs=1e-2) == value(
            btx_ftpz.fs.unit.vapor_reboil.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_ftpz.fs.unit.vapor_reboil.pressure[0]
        )

        # Unit level
        assert pytest.approx(16926.5, rel=1e-4) == value(btx_ftpz.fs.unit.heat_duty[0])

        # Reboiler when using FcTP

        # Bottoms port
        assert pytest.approx(0.19455, abs=1e-3) == value(
            btx_fctp.fs.unit.bottoms.flow_mol_comp[0, "benzene"]
        )
        assert pytest.approx(0.30545, abs=1e-3) == value(
            btx_fctp.fs.unit.bottoms.flow_mol_comp[0, "toluene"]
        )
        assert pytest.approx(368.728, abs=1e-2) == value(
            btx_fctp.fs.unit.bottoms.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_fctp.fs.unit.bottoms.pressure[0]
        )

        # Vapor reboil port
        assert pytest.approx(0.3054, abs=1e-3) == value(
            btx_fctp.fs.unit.vapor_reboil.flow_mol_comp[0, "benzene"]
        )
        assert pytest.approx(0.1946, abs=1e-3) == value(
            btx_fctp.fs.unit.vapor_reboil.flow_mol_comp[0, "toluene"]
        )
        assert pytest.approx(368.728, abs=1e-2) == value(
            btx_fctp.fs.unit.vapor_reboil.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_fctp.fs.unit.vapor_reboil.pressure[0]
        )

        # Unit level
        assert pytest.approx(16926.5, rel=1e-4) == value(btx_fctp.fs.unit.heat_duty[0])

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx_ftpz, btx_fctp):
        assert (
            abs(
                value(
                    btx_ftpz.fs.unit.inlet.flow_mol[0]
                    - (
                        btx_ftpz.fs.unit.bottoms.flow_mol[0]
                        + btx_ftpz.fs.unit.vapor_reboil.flow_mol[0]
                    )
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    btx_fctp.fs.unit.inlet.flow_mol_comp[0, "benzene"]
                    + btx_fctp.fs.unit.inlet.flow_mol_comp[0, "toluene"]
                    - (
                        btx_fctp.fs.unit.bottoms.flow_mol_comp[0, "benzene"]
                        + btx_fctp.fs.unit.bottoms.flow_mol_comp[0, "toluene"]
                        + btx_fctp.fs.unit.vapor_reboil.flow_mol_comp[0, "benzene"]
                        + btx_fctp.fs.unit.vapor_reboil.flow_mol_comp[0, "toluene"]
                    )
                )
            )
            <= 1e-6
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, btx_ftpz, btx_fctp):
        btx_ftpz.fs.unit.report()
        btx_fctp.fs.unit.report()
