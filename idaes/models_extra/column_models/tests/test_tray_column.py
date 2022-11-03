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
Tests for tray column unit model (single feed tray, no side draws).

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import check_optimal_termination, ConcreteModel, value
from pyomo.util.check_units import assert_units_consistent
import pyomo.common.unittest as unittest

from idaes.core import FlowsheetBlock
from idaes.models_extra.column_models import TrayColumn
from idaes.models_extra.column_models.condenser import CondenserType, TemperatureSpec
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.performance import PerformanceBaseClass

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.properties.modular_properties.examples.BT_ideal import configuration

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@pytest.mark.unit
def test_config():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = TrayColumn(
        number_of_trays=10,
        feed_tray_location=5,
        condenser_type=CondenserType.totalCondenser,
        condenser_temperature_spec=TemperatureSpec.customTemperature,
        property_package=m.fs.properties,
        has_heat_transfer=False,
        has_pressure_change=False,
    )

    assert len(m.fs.unit.config) == 12

    assert m.fs.unit.feed_tray.config.is_feed_tray

    assert hasattr(m.fs.unit, "condenser")
    assert hasattr(m.fs.unit, "reboiler")

    assert hasattr(m.fs.unit, "rectification_section")
    assert hasattr(m.fs.unit, "stripping_section")


def build_model_btx_ftpz():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal"
    )

    m.fs.unit = TrayColumn(
        number_of_trays=10,
        feed_tray_location=5,
        condenser_type=CondenserType.totalCondenser,
        condenser_temperature_spec=TemperatureSpec.atBubblePoint,
        property_package=m.fs.properties,
        has_heat_transfer=False,
        has_pressure_change=False,
    )

    # Inlet feed conditions
    m.fs.unit.feed.flow_mol.fix(40)
    m.fs.unit.feed.temperature.fix(368)
    m.fs.unit.feed.pressure.fix(101325)
    m.fs.unit.feed.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.unit.feed.mole_frac_comp[0, "toluene"].fix(0.5)

    # unit level inputs
    m.fs.unit.condenser.reflux_ratio.fix(1.4)
    m.fs.unit.condenser.condenser_pressure.fix(101325)

    m.fs.unit.reboiler.boilup_ratio.fix(1.3)

    return m


@pytest.mark.performance
class Test_TrayColumn_Performance(PerformanceBaseClass, unittest.TestCase):
    def build_model(self):
        return build_model_btx_ftpz()

    def initialize_model(self, model):
        model.fs.unit.initialize()


class TestBTXIdeal:
    @pytest.fixture(scope="class")
    def btx_ftpz(self):
        return build_model_btx_ftpz()

    @pytest.fixture(scope="class")
    def btx_fctp(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FcTP"
        )

        m.fs.unit = TrayColumn(
            number_of_trays=10,
            feed_tray_location=5,
            condenser_type=CondenserType.totalCondenser,
            condenser_temperature_spec=TemperatureSpec.atBubblePoint,
            property_package=m.fs.properties,
            has_heat_transfer=False,
            has_pressure_change=False,
        )

        # Inlet feed conditions
        m.fs.unit.feed.flow_mol_comp[0, "benzene"].fix(20)
        m.fs.unit.feed.flow_mol_comp[0, "toluene"].fix(20)
        m.fs.unit.feed.temperature.fix(368)
        m.fs.unit.feed.pressure.fix(101325)

        # unit level inputs
        m.fs.unit.condenser.reflux_ratio.fix(1.4)
        m.fs.unit.condenser.condenser_pressure.fix(101325)

        m.fs.unit.reboiler.boilup_ratio.fix(1.3)

        return m

    @pytest.mark.unit
    def test_build(self, btx_ftpz, btx_fctp):
        assert len(btx_ftpz.fs.unit.config) == 12

        assert btx_ftpz.fs.unit.feed_tray.config.is_feed_tray

        assert hasattr(btx_ftpz.fs.unit, "condenser")
        assert hasattr(btx_ftpz.fs.unit, "reboiler")

        assert hasattr(btx_ftpz.fs.unit, "rectification_section")
        assert hasattr(btx_ftpz.fs.unit, "stripping_section")

        assert len(btx_fctp.fs.unit.config) == 12

        assert btx_fctp.fs.unit.feed_tray.config.is_feed_tray

        assert hasattr(btx_fctp.fs.unit, "condenser")
        assert hasattr(btx_fctp.fs.unit, "reboiler")

        assert hasattr(btx_fctp.fs.unit, "rectification_section")
        assert hasattr(btx_fctp.fs.unit, "stripping_section")

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

        assert check_optimal_termination(results)

        results = solver.solve(btx_fctp)

        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx_ftpz, btx_fctp):

        # Distillate port - btx_ftpz
        assert pytest.approx(18.978, rel=1e-2) == value(
            btx_ftpz.fs.unit.condenser.distillate.flow_mol[0]
        )
        assert pytest.approx(0.8882, rel=1e-2) == value(
            btx_ftpz.fs.unit.condenser.distillate.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.1118, rel=1e-2) == value(
            btx_ftpz.fs.unit.condenser.distillate.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(355.642, abs=1e-2) == value(
            btx_ftpz.fs.unit.condenser.distillate.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_ftpz.fs.unit.condenser.distillate.pressure[0]
        )

        # Bottoms port - btx_ftpz
        assert pytest.approx(21.021, rel=1e-2) == value(
            btx_ftpz.fs.unit.reboiler.bottoms.flow_mol[0]
        )
        assert pytest.approx(0.1495, rel=1e-2) == value(
            btx_ftpz.fs.unit.reboiler.bottoms.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.8504, rel=1e-3) == value(
            btx_ftpz.fs.unit.reboiler.bottoms.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(377.337, abs=1e-2) == value(
            btx_ftpz.fs.unit.reboiler.bottoms.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_ftpz.fs.unit.reboiler.bottoms.pressure[0]
        )

        # Distillate port - btx_fctp
        assert pytest.approx(16.856, rel=1e-2) == value(
            btx_fctp.fs.unit.condenser.distillate.flow_mol_comp[0, "benzene"]
        )
        assert pytest.approx(2.121, rel=1e-2) == value(
            btx_fctp.fs.unit.condenser.distillate.flow_mol_comp[0, "toluene"]
        )
        assert pytest.approx(355.638, abs=1e-2) == value(
            btx_fctp.fs.unit.condenser.distillate.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_fctp.fs.unit.condenser.distillate.pressure[0]
        )

        # # Bottoms port - btx_fctp
        assert pytest.approx(3.179, rel=1e-2) == value(
            btx_fctp.fs.unit.reboiler.bottoms.flow_mol_comp[0, "benzene"]
        )
        assert pytest.approx(17.876, rel=1e-3) == value(
            btx_fctp.fs.unit.reboiler.bottoms.flow_mol_comp[0, "toluene"]
        )
        assert pytest.approx(377.28, abs=1e-2) == value(
            btx_fctp.fs.unit.reboiler.bottoms.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_fctp.fs.unit.reboiler.bottoms.pressure[0]
        )


class TestBTXIdealGeneric:
    @pytest.fixture(scope="class")
    def btx_ftpz_generic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = GenericParameterBlock(**configuration)

        m.fs.unit = TrayColumn(
            number_of_trays=10,
            feed_tray_location=5,
            condenser_type=CondenserType.totalCondenser,
            condenser_temperature_spec=TemperatureSpec.atBubblePoint,
            property_package=m.fs.properties,
            has_heat_transfer=False,
            has_pressure_change=False,
        )

        # Inlet feed conditions
        m.fs.unit.feed.flow_mol.fix(100)
        m.fs.unit.feed.temperature.fix(368)
        m.fs.unit.feed.pressure.fix(101325)
        m.fs.unit.feed.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.feed.mole_frac_comp[0, "toluene"].fix(0.5)

        # unit level inputs
        m.fs.unit.condenser.reflux_ratio.fix(1.4)
        m.fs.unit.condenser.condenser_pressure.fix(101325)

        m.fs.unit.reboiler.boilup_ratio.fix(1.3)

        iscale.calculate_scaling_factors(m.fs.unit)

        return m

    @pytest.mark.unit
    def test_build(self, btx_ftpz_generic):
        assert len(btx_ftpz_generic.fs.unit.config) == 12

        assert btx_ftpz_generic.fs.unit.feed_tray.config.is_feed_tray

        assert hasattr(btx_ftpz_generic.fs.unit, "condenser")
        assert hasattr(btx_ftpz_generic.fs.unit, "reboiler")

        assert hasattr(btx_ftpz_generic.fs.unit, "rectification_section")
        assert hasattr(btx_ftpz_generic.fs.unit, "stripping_section")

        assert len(btx_ftpz_generic.fs.unit.config) == 12

    @pytest.mark.unit
    def test_dof(self, btx_ftpz_generic):
        assert degrees_of_freedom(btx_ftpz_generic) == 0

    @pytest.mark.component
    def test_units_FTPz(self, btx_ftpz_generic):
        assert_units_consistent(btx_ftpz_generic)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx_ftpz_generic):
        initialization_tester(btx_ftpz_generic)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx_ftpz_generic):
        results = solver.solve(btx_ftpz_generic)

        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx_ftpz_generic):

        # Distillate port - btx_ftpz
        assert pytest.approx(49.28, rel=1e-2) == value(
            btx_ftpz_generic.fs.unit.condenser.distillate.flow_mol[0]
        )
        assert pytest.approx(0.8784, rel=1e-2) == value(
            btx_ftpz_generic.fs.unit.condenser.distillate.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.1215, rel=1e-2) == value(
            btx_ftpz_generic.fs.unit.condenser.distillate.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(355.853, abs=1e-2) == value(
            btx_ftpz_generic.fs.unit.condenser.distillate.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_ftpz_generic.fs.unit.condenser.distillate.pressure[0]
        )

        # Bottoms port - btx_ftpz
        assert pytest.approx(50.711, rel=1e-2) == value(
            btx_ftpz_generic.fs.unit.reboiler.bottoms.flow_mol[0]
        )
        assert pytest.approx(0.1322, rel=1e-2) == value(
            btx_ftpz_generic.fs.unit.reboiler.bottoms.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.8677, rel=1e-3) == value(
            btx_ftpz_generic.fs.unit.reboiler.bottoms.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(378.0433, abs=1e-2) == value(
            btx_ftpz_generic.fs.unit.reboiler.bottoms.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx_ftpz_generic.fs.unit.reboiler.bottoms.pressure[0]
        )
