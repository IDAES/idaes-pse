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
Tests for Flash unit model.
Author: Jaffer Ghouse
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    value,
    units,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models.unit_models.flash import Flash, EnergySplittingType
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties import iapws95
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

    m.fs.properties = PhysicalParameterTestBlock()

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
    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.unit = Flash(property_package=m.fs.properties)
    iscale.calculate_scaling_factors(m)


# -----------------------------------------------------------------------------
class TestBTXIdeal(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal"
        )

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

        assert number_variables(btx) == 48
        assert number_total_constraints(btx) == 41
        assert number_unused_variables(btx) == 0

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_consistent(btx)
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

        assert pytest.approx(0.412, abs=1e-3) == value(
            btx.fs.unit.liq_outlet.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.588, abs=1e-3) == value(
            btx.fs.unit.liq_outlet.mole_frac_comp[0, "toluene"]
        )
        assert pytest.approx(0.634, abs=1e-3) == value(
            btx.fs.unit.vap_outlet.mole_frac_comp[0, "benzene"]
        )
        assert pytest.approx(0.366, abs=1e-3) == value(
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

        m.fs.unit = Flash(
            property_package=m.fs.properties,
            ideal_separation=False,
            energy_split_basis=EnergySplittingType.enthalpy_split,
        )

        m.fs.unit.inlet.flow_mol.fix(100)
        m.fs.unit.inlet.enth_mol.fix(24000)
        m.fs.unit.inlet.pressure.fix(101325)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iapws):
        assert len(iapws.fs.unit.inlet.vars) == 3
        assert hasattr(iapws.fs.unit.inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.inlet, "pressure")

        assert hasattr(iapws.fs.unit, "vap_outlet")
        assert len(iapws.fs.unit.vap_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.vap_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.vap_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.vap_outlet, "pressure")

        assert hasattr(iapws.fs.unit, "liq_outlet")
        assert len(iapws.fs.unit.liq_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.liq_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.liq_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.liq_outlet, "pressure")

        assert hasattr(iapws.fs.unit, "split")
        assert hasattr(iapws.fs.unit, "heat_duty")
        assert hasattr(iapws.fs.unit, "deltaP")

        assert number_variables(iapws) == 18
        assert number_total_constraints(iapws) == 13
        assert number_unused_variables(iapws) == 0

    @pytest.mark.component
    def test_units(self, iapws):
        assert_units_consistent(iapws)
        assert_units_equivalent(iapws.fs.unit.heat_duty[0], units.W)
        assert_units_equivalent(iapws.fs.unit.deltaP[0], units.Pa)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Heat Duty": iapws.fs.unit.heat_duty[0],
                "Pressure Change": iapws.fs.unit.deltaP[0],
            }
        }

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
            "Inlet": {
                "Molar Flow": pytest.approx(100, rel=1e-4),
                "Mass Flow": pytest.approx(1.8015, rel=1e-4),
                "T": pytest.approx(373.12429584768876, rel=1e-4),
                "P": pytest.approx(101325.0, rel=1e-4),
                "Vapor Fraction": pytest.approx(0.40467852502291135, abs=1e-4),
                "Molar Enthalpy": pytest.approx(24000.0, rel=1e-4),
            },
            "Vapor Outlet": {
                "Molar Flow": pytest.approx(1, rel=1e-4),
                "Mass Flow": pytest.approx(1.8015e-2, rel=1e-4),
                "T": pytest.approx(270.4877112932641, rel=1e-4),
                "P": pytest.approx(11032305.8275, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(0.01102138712926277, rel=1e-4),
            },
            "Liquid Outlet": {
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
        t_in = value(iapws.fs.unit.control_volume.properties_in[0].temperature)
        assert pytest.approx(373.12, abs=1e-2) == t_in
        assert pytest.approx(t_in, abs=1e-2) == value(
            iapws.fs.unit.control_volume.properties_out[0].temperature
        )
        assert pytest.approx(t_in, abs=1e-2) == value(
            iapws.fs.unit.split.Vap_state[0].temperature
        )
        assert pytest.approx(t_in, abs=1e-2) == value(
            iapws.fs.unit.split.Liq_state[0].temperature
        )

        assert pytest.approx(101325.0, abs=1e3) == value(
            iapws.fs.unit.vap_outlet.pressure[0]
        )
        assert pytest.approx(48200, abs=1e3) == value(
            iapws.fs.unit.vap_outlet.enth_mol[0]
        )
        assert pytest.approx(40.467, abs=1e-3) == value(
            iapws.fs.unit.vap_outlet.flow_mol[0]
        )

        assert pytest.approx(101325.0, abs=1e3) == value(
            iapws.fs.unit.liq_outlet.pressure[0]
        )
        assert pytest.approx(7549.4, abs=1) == value(
            iapws.fs.unit.liq_outlet.enth_mol[0]
        )
        assert pytest.approx(59.532, abs=1e-3) == value(
            iapws.fs.unit.liq_outlet.flow_mol[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert (
            abs(
                value(
                    iapws.fs.unit.inlet.flow_mol[0]
                    - (
                        iapws.fs.unit.vap_outlet.flow_mol[0]
                        + iapws.fs.unit.liq_outlet.flow_mol[0]
                    )
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    iapws.fs.unit.inlet.flow_mol[0] * iapws.fs.unit.inlet.enth_mol[0]
                    - (
                        iapws.fs.unit.vap_outlet.flow_mol[0]
                        * iapws.fs.unit.vap_outlet.enth_mol[0]
                        + iapws.fs.unit.liq_outlet.flow_mol[0]
                        * iapws.fs.unit.liq_outlet.enth_mol[0]
                    )
                )
            )
            <= 1e-6
        )
