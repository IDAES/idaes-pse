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
Tests for feed with flash.
Authors: Andrew Lee
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    value,
    units as pyunits,
)
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.models.unit_models.feed_flash import FeedFlash, FlashType
from idaes.models.properties import iapws95
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
from pyomo.util.check_units import assert_units_consistent


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = FeedFlash(property_package=m.fs.properties)

    # Check unit config arguments
    assert len(m.fs.unit.config) == 6

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.flash_type == FlashType.isothermal
    assert m.fs.unit.config.property_package is m.fs.properties


# -----------------------------------------------------------------------------
class TestBTXIdeal(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal"
        )

        m.fs.unit = FeedFlash(property_package=m.fs.properties)

        m.fs.unit.flow_mol.fix(1)
        m.fs.unit.temperature.fix(368)
        m.fs.unit.pressure.fix(101325)
        m.fs.unit.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.mole_frac_comp[0, "toluene"].fix(0.5)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "flow_mol")
        assert hasattr(btx.fs.unit, "mole_frac_comp")
        assert hasattr(btx.fs.unit, "temperature")
        assert hasattr(btx.fs.unit, "pressure")

        assert hasattr(btx.fs.unit, "outlet")
        assert len(btx.fs.unit.outlet.vars) == 4
        assert hasattr(btx.fs.unit.outlet, "flow_mol")
        assert hasattr(btx.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.outlet, "temperature")
        assert hasattr(btx.fs.unit.outlet, "pressure")

        assert hasattr(btx.fs.unit, "isothermal")

        assert number_variables(btx) == 34
        assert number_total_constraints(btx) == 29
        assert number_unused_variables(btx) == 0

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict is None

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "flow_mol": getattr(pyunits.pint_registry, "mole/second"),
                "mole_frac_comp benzene": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "mole_frac_comp toluene": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "temperature": getattr(pyunits.pint_registry, "kelvin"),
                "pressure": getattr(pyunits.pint_registry, "Pa"),
            },
            "Outlet": {
                "flow_mol": pytest.approx(1.0, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(298.15, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert pytest.approx(101325.0, abs=1e3) == value(btx.fs.unit.outlet.pressure[0])
        assert pytest.approx(368.00, abs=1e-0) == value(
            btx.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(1.0, abs=1e-2) == value(btx.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(0.396, abs=1e-3) == value(
            btx.fs.unit.control_volume.properties_out[0].flow_mol_phase["Vap"]
        )


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock(
            phase_presentation=iapws95.PhaseType.LG
        )

        m.fs.unit = FeedFlash(
            property_package=m.fs.properties, flash_type=FlashType.isenthalpic
        )

        m.fs.unit.flow_mol.fix(100)
        m.fs.unit.enth_mol.fix(24000)
        m.fs.unit.pressure.fix(101325)

        return m

    @pytest.mark.unit
    def test_build(self, iapws):
        assert hasattr(iapws.fs.unit, "flow_mol")
        assert hasattr(iapws.fs.unit, "enth_mol")
        assert hasattr(iapws.fs.unit, "pressure")

        assert hasattr(iapws.fs.unit, "outlet")
        assert len(iapws.fs.unit.outlet.vars) == 3
        assert hasattr(iapws.fs.unit.outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.outlet, "pressure")

        assert hasattr(iapws.fs.unit, "isenthalpic")

        assert number_variables(iapws) == 6
        assert number_total_constraints(iapws) == 3
        assert number_unused_variables(iapws) == 0

    @pytest.mark.component
    def test_units(self, iapws):
        assert_units_consistent(iapws)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict is None

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, iapws):
        stable = iapws.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Molar Flow": getattr(pyunits.pint_registry, "mole/second"),
                "Mass Flow": getattr(pyunits.pint_registry, "kg/second"),
                "T": getattr(pyunits.pint_registry, "K"),
                "P": getattr(pyunits.pint_registry, "Pa"),
                "Vapor Fraction": getattr(pyunits.pint_registry, "dimensionless"),
                "Molar Enthalpy": getattr(pyunits.pint_registry, "J/mole"),
            },
            "Outlet": {
                "Molar Flow": pytest.approx(1, rel=1e-4),
                "Mass Flow": pytest.approx(1.8015e-2, rel=1e-4),
                "T": pytest.approx(270.4877112932641, rel=1e-4),
                "P": pytest.approx(11032305.8275, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(0.01102138712926277, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

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
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(101325.0, abs=1e3) == value(
            iapws.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(24000, abs=1e3) == value(iapws.fs.unit.outlet.enth_mol[0])
        assert pytest.approx(100.0, abs=1e-2) == value(iapws.fs.unit.outlet.flow_mol[0])

        assert pytest.approx(373.12, abs=1e-2) == value(
            iapws.fs.unit.control_volume.properties_out[0].temperature
        )
        assert pytest.approx(0.5953, abs=1e-4) == value(
            iapws.fs.unit.control_volume.properties_out[0].phase_frac["Liq"]
        )
