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
Tests for conventional tray unit model (no feed, no side draws).

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import check_optimal_termination, ConcreteModel, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.models_extra.column_models import Tray
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

    m.fs.unit = Tray(
        property_package=m.fs.properties,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    assert len(m.fs.unit.config) == 9

    assert not m.fs.unit.config.is_feed_tray
    assert not m.fs.unit.config.has_liquid_side_draw
    assert not m.fs.unit.config.has_vapor_side_draw


class TestBTXIdeal:
    @pytest.fixture(scope="class")
    def btx_ftpz(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal"
        )
        m.fs.unit = Tray(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        # Fix the tray inputs (FTPz)
        m.fs.unit.liq_in.flow_mol.fix(1)
        m.fs.unit.liq_in.temperature.fix(369)
        m.fs.unit.liq_in.pressure.fix(101325)
        m.fs.unit.liq_in.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.liq_in.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.vap_in.flow_mol.fix(1)
        m.fs.unit.vap_in.temperature.fix(372)
        m.fs.unit.vap_in.pressure.fix(101325)
        m.fs.unit.vap_in.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.vap_in.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.deltaP.fix(0)
        m.fs.unit.heat_duty.fix(0)

        return m

    @pytest.fixture(scope="class")
    def btx_fctp(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FcTP"
        )
        m.fs.unit = Tray(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        # Fix the tray inputs (FcTP)
        m.fs.unit.liq_in.flow_mol_comp[0, "benzene"].fix(0.5)
        m.fs.unit.liq_in.flow_mol_comp[0, "toluene"].fix(0.5)
        m.fs.unit.liq_in.temperature.fix(369)
        m.fs.unit.liq_in.pressure.fix(101325)

        m.fs.unit.vap_in.flow_mol_comp[0, "benzene"].fix(0.5)
        m.fs.unit.vap_in.flow_mol_comp[0, "toluene"].fix(0.5)
        m.fs.unit.vap_in.temperature.fix(372)
        m.fs.unit.vap_in.pressure.fix(101325)

        m.fs.unit.deltaP.fix(0)
        m.fs.unit.heat_duty.fix(0)

        return m

    @pytest.mark.unit
    def test_build(self, btx_ftpz, btx_fctp):
        # General build
        assert hasattr(btx_ftpz.fs.unit, "material_mixing_equations")
        assert hasattr(btx_ftpz.fs.unit, "enthalpy_mixing_equations")
        assert hasattr(btx_ftpz.fs.unit, "pressure_drop_equation")

        assert btx_ftpz.fs.unit.config.has_heat_transfer
        assert hasattr(btx_ftpz.fs.unit, "heat_duty")

        assert btx_ftpz.fs.unit.config.has_pressure_change
        assert hasattr(btx_ftpz.fs.unit, "deltaP")

        # State blocks
        assert hasattr(btx_ftpz.fs.unit, "properties_in_liq")
        assert hasattr(btx_ftpz.fs.unit, "properties_in_vap")
        assert hasattr(btx_ftpz.fs.unit, "properties_out")

        # Ports
        assert hasattr(btx_ftpz.fs.unit, "liq_in")
        assert hasattr(btx_ftpz.fs.unit.liq_in, "flow_mol")
        assert hasattr(btx_ftpz.fs.unit.liq_in, "mole_frac_comp")
        assert hasattr(btx_ftpz.fs.unit.liq_in, "temperature")
        assert hasattr(btx_ftpz.fs.unit.liq_in, "pressure")

        assert hasattr(btx_ftpz.fs.unit, "vap_in")
        assert hasattr(btx_ftpz.fs.unit.vap_in, "flow_mol")
        assert hasattr(btx_ftpz.fs.unit.vap_in, "mole_frac_comp")
        assert hasattr(btx_ftpz.fs.unit.vap_in, "temperature")
        assert hasattr(btx_ftpz.fs.unit.vap_in, "pressure")

        assert hasattr(btx_ftpz.fs.unit, "liq_out")
        assert hasattr(btx_ftpz.fs.unit.liq_out, "flow_mol")
        assert hasattr(btx_ftpz.fs.unit.liq_out, "mole_frac_comp")
        assert hasattr(btx_ftpz.fs.unit.liq_out, "temperature")
        assert hasattr(btx_ftpz.fs.unit.liq_out, "pressure")

        assert hasattr(btx_ftpz.fs.unit, "vap_out")
        assert hasattr(btx_ftpz.fs.unit.vap_out, "flow_mol")
        assert hasattr(btx_ftpz.fs.unit.vap_out, "mole_frac_comp")
        assert hasattr(btx_ftpz.fs.unit.vap_out, "temperature")
        assert hasattr(btx_ftpz.fs.unit.vap_out, "pressure")

        assert not hasattr(btx_ftpz.fs.unit, "feed")

        assert not hasattr(btx_ftpz.fs.unit, "liq_side_draw")
        assert not hasattr(btx_ftpz.fs.unit, "liq_side_sf")

        assert not hasattr(btx_ftpz.fs.unit, "vap_side_draw")
        assert not hasattr(btx_ftpz.fs.unit, "vap_side_sf")

        assert number_variables(btx_ftpz.fs.unit) == 71
        assert number_total_constraints(btx_ftpz.fs.unit) == 59
        assert number_unused_variables(btx_ftpz) == 0

        # General build
        assert hasattr(btx_fctp.fs.unit, "material_mixing_equations")
        assert hasattr(btx_fctp.fs.unit, "enthalpy_mixing_equations")
        assert hasattr(btx_fctp.fs.unit, "pressure_drop_equation")

        assert btx_fctp.fs.unit.config.has_heat_transfer
        assert hasattr(btx_fctp.fs.unit, "heat_duty")

        assert btx_fctp.fs.unit.config.has_pressure_change
        assert hasattr(btx_fctp.fs.unit, "deltaP")

        # State blocks
        assert hasattr(btx_fctp.fs.unit, "properties_in_liq")
        assert hasattr(btx_fctp.fs.unit, "properties_in_vap")
        assert hasattr(btx_fctp.fs.unit, "properties_out")

        # Ports
        assert hasattr(btx_fctp.fs.unit, "liq_in")
        assert hasattr(btx_fctp.fs.unit.liq_in, "flow_mol_comp")
        assert hasattr(btx_fctp.fs.unit.liq_in, "temperature")
        assert hasattr(btx_fctp.fs.unit.liq_in, "pressure")

        assert hasattr(btx_fctp.fs.unit, "vap_in")
        assert hasattr(btx_fctp.fs.unit.vap_in, "flow_mol_comp")
        assert hasattr(btx_fctp.fs.unit.vap_in, "temperature")
        assert hasattr(btx_fctp.fs.unit.vap_in, "pressure")

        assert hasattr(btx_fctp.fs.unit, "liq_out")
        assert hasattr(btx_fctp.fs.unit.liq_out, "flow_mol_comp")
        assert hasattr(btx_fctp.fs.unit.liq_out, "temperature")
        assert hasattr(btx_fctp.fs.unit.liq_out, "pressure")

        assert hasattr(btx_fctp.fs.unit, "vap_out")
        assert hasattr(btx_fctp.fs.unit.vap_out, "flow_mol_comp")
        assert hasattr(btx_fctp.fs.unit.vap_out, "temperature")
        assert hasattr(btx_fctp.fs.unit.vap_out, "pressure")

        assert not hasattr(btx_fctp.fs.unit, "feed")

        assert not hasattr(btx_fctp.fs.unit, "liq_side_draw")
        assert not hasattr(btx_fctp.fs.unit, "liq_side_sf")

        assert not hasattr(btx_fctp.fs.unit, "vap_side_draw")
        assert not hasattr(btx_fctp.fs.unit, "vap_side_sf")

        assert number_variables(btx_fctp.fs.unit) == 74
        assert number_total_constraints(btx_fctp.fs.unit) == 64
        assert number_unused_variables(btx_fctp) == 0

    @pytest.mark.unit
    def test_dof(self, btx_ftpz, btx_fctp):
        assert degrees_of_freedom(btx_ftpz.fs.unit) == 0
        assert degrees_of_freedom(btx_fctp.fs.unit) == 0

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

        # liq_out port
        assert pytest.approx(0.463370, abs=1e-3) == value(
            btx_ftpz.fs.unit.liq_out.flow_mol[0]
        )
        assert pytest.approx(0.33313, abs=1e-3) == value(
            btx_ftpz.fs.unit.liq_out.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.66686, abs=1e-3) == value(
            btx_ftpz.fs.unit.liq_out.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(370.567, abs=1e-3) == value(
            btx_ftpz.fs.unit.liq_out.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_ftpz.fs.unit.liq_out.pressure[0]
        )

        # vap_out port
        assert pytest.approx(1.5366, abs=1e-3) == value(
            btx_ftpz.fs.unit.vap_out.flow_mol[0]
        )
        assert pytest.approx(0.55031, abs=1e-3) == value(
            btx_ftpz.fs.unit.vap_out.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.44968, abs=1e-3) == value(
            btx_ftpz.fs.unit.vap_out.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(370.567, abs=1e-3) == value(
            btx_ftpz.fs.unit.vap_out.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_ftpz.fs.unit.vap_out.pressure[0]
        )

        # liq_out port
        assert pytest.approx(0.15436, abs=1e-3) == value(
            btx_fctp.fs.unit.liq_out.flow_mol_comp[0, "benzene"]
        )
        assert pytest.approx(0.30900, abs=1e-3) == value(
            btx_fctp.fs.unit.liq_out.flow_mol_comp[0, "toluene"]
        )
        assert pytest.approx(370.567, abs=1e-3) == value(
            btx_fctp.fs.unit.liq_out.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_fctp.fs.unit.liq_out.pressure[0]
        )

        # vap_out port
        assert pytest.approx(0.84560, abs=1e-3) == value(
            btx_fctp.fs.unit.vap_out.flow_mol_comp[0, "benzene"]
        )
        assert pytest.approx(0.69097, abs=1e-3) == value(
            btx_fctp.fs.unit.vap_out.flow_mol_comp[0, "toluene"]
        )
        assert pytest.approx(370.567, abs=1e-3) == value(
            btx_fctp.fs.unit.vap_out.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_fctp.fs.unit.vap_out.pressure[0]
        )
