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
Tests for Flash unit model based on the complementarity formulation.
Author: Vibhav Dabadghao
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    value,
    units,
    Var,
    units as pyunits,
    Objective,
    SolverFactory,
    log,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models.unit_models.flash import Flash, EnergySplittingType
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.examples.BT_PR import configuration

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.unit = Flash(property_package=m.fs.properties)

    # Check unit config arguments
    assert len(m.fs.unit.config) == 11

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.ideal_separation
    assert m.fs.unit.config.has_heat_transfer
    assert m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.energy_split_basis == EnergySplittingType.equal_temperature


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_calc_scale():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GenericParameterBlock(**configuration)
    m.fs.unit = Flash(property_package=m.fs.properties)
    iscale.calculate_scaling_factors(m)


# -----------------------------------------------------------------------------
class TestBTXPengRobinson(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = GenericParameterBlock(**configuration)
        m.fs.unit = Flash(property_package=m.fs.properties)

        m.fs.unit.inlet.flow_mol.fix(1)
        m.fs.unit.inlet.temperature.fix(368)
        m.fs.unit.inlet.pressure.fix(101325)
        m.fs.unit.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit.inlet, "flow_mol")
        assert hasattr(btx.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.inlet, "temperature")
        assert hasattr(btx.fs.unit.inlet, "pressure")

        assert hasattr(btx.fs.unit, "vap_outlet")
        assert len(btx.fs.unit.vap_outlet.vars) == 4
        assert hasattr(btx.fs.unit.vap_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.vap_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.vap_outlet, "temperature")
        assert hasattr(btx.fs.unit.vap_outlet, "pressure")

        assert hasattr(btx.fs.unit, "liq_outlet")
        assert len(btx.fs.unit.liq_outlet.vars) == 4
        assert hasattr(btx.fs.unit.liq_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.liq_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.liq_outlet, "temperature")
        assert hasattr(btx.fs.unit.liq_outlet, "pressure")

        assert hasattr(btx.fs.unit, "split")
        assert hasattr(btx.fs.unit, "heat_duty")
        assert hasattr(btx.fs.unit, "deltaP")

        assert number_variables(btx) == 84
        assert number_total_constraints(btx) == 43
        assert number_unused_variables(btx) == 14

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_equivalent(btx.fs.unit.heat_duty[0], units.W)
        assert_units_equivalent(btx.fs.unit.deltaP[0], units.Pa)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Heat Duty": btx.fs.unit.heat_duty[0],
                "Pressure Change": btx.fs.unit.deltaP[0],
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "flow_mol": getattr(pyunits.pint_registry, "mole/s"),
                "mole_frac_comp benzene": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "mole_frac_comp toluene": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "temperature": getattr(pyunits.pint_registry, "K"),
                "pressure": getattr(pyunits.pint_registry, "Pa"),
            },
            "Inlet": {
                "flow_mol": pytest.approx(1.00, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(368, rel=1e-4),
                "pressure": pytest.approx(101325, rel=1e-4),
            },
            "Vapor Outlet": {
                "flow_mol": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(298.15, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Liquid Outlet": {
                "flow_mol": pytest.approx(0.5, rel=1e-4),
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
        assert pytest.approx(0.603, abs=1e-3) == value(
            btx.fs.unit.liq_outlet.flow_mol[0]
        )
        assert pytest.approx(0.396, abs=1e-3) == value(
            btx.fs.unit.vap_outlet.flow_mol[0]
        )
        assert pytest.approx(368, abs=1e-3) == value(
            btx.fs.unit.liq_outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx.fs.unit.liq_outlet.pressure[0]
        )

        assert pytest.approx(0.415, abs=1e-3) == value(
            btx.fs.unit.liq_outlet.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.584, abs=1e-3) == value(
            btx.fs.unit.liq_outlet.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(0.628, abs=1e-3) == value(
            btx.fs.unit.vap_outlet.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.372, abs=1e-3) == value(
            btx.fs.unit.vap_outlet.mole_frac_comp[0, "toluene"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.inlet.flow_mol[0]
                    - (
                        btx.fs.unit.vap_outlet.flow_mol[0]
                        + btx.fs.unit.liq_outlet.flow_mol[0]
                    )
                )
            )
            <= 1e-6
        )
