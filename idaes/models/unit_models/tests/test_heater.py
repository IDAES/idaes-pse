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
Tests for 0D heat exchanger models.

Author: John Eslick
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    value,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models.unit_models import Heater

from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties import iapws95
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.modular_properties.examples.BT_PR import configuration

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = Heater(property_package=m.fs.properties)

    # Check unit config arguments
    assert len(m.fs.unit.config) == 9
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_phase_equilibrium
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties


# -----------------------------------------------------------------------------
class TestBTX(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(valid_phase="Liq")

        m.fs.unit = Heater(property_package=m.fs.properties, has_pressure_change=True)

        m.fs.unit.inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.inlet.temperature[0].fix(365)  # K
        m.fs.unit.inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.heat_duty.fix(-5000)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "inlet")
        assert len(btx.fs.unit.inlet.vars) == 4
        assert hasattr(btx.fs.unit.inlet, "flow_mol")
        assert hasattr(btx.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.inlet, "temperature")
        assert hasattr(btx.fs.unit.inlet, "pressure")

        assert hasattr(btx.fs.unit, "outlet")
        assert len(btx.fs.unit.outlet.vars) == 4
        assert hasattr(btx.fs.unit.outlet, "flow_mol")
        assert hasattr(btx.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.outlet, "temperature")
        assert hasattr(btx.fs.unit.outlet, "pressure")

        assert hasattr(btx.fs.unit, "heat_duty")
        assert hasattr(btx.fs.unit, "deltaP")

        assert number_variables(btx) == 24
        assert number_total_constraints(btx) == 17
        assert number_unused_variables(btx) == 0

    @pytest.mark.integration
    def test_units(self, btx):
        assert_units_equivalent(btx.fs.unit.control_volume.heat, pyunits.J / pyunits.s)
        assert_units_equivalent(btx.fs.unit.heat_duty[0], pyunits.J / pyunits.s)
        assert_units_equivalent(btx.fs.unit.deltaP[0], pyunits.Pa)
        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

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
        assert pytest.approx(5, abs=1e-3) == value(btx.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(358.9, abs=1e-1) == value(
            btx.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(btx.fs.unit.outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(value(btx.fs.unit.inlet.flow_mol[0] - btx.fs.unit.outlet.flow_mol[0]))
            <= 1e-6
        )

        assert (
            abs(
                value(
                    btx.fs.unit.inlet.flow_mol[0]
                    * btx.fs.unit.control_volume.properties_in[0].enth_mol_phase["Liq"]
                    - btx.fs.unit.outlet.flow_mol[0]
                    * btx.fs.unit.control_volume.properties_out[0].enth_mol_phase["Liq"]
                    + btx.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {"vars": {"Heat Duty": btx.fs.unit.heat_duty[0]}}


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = Heater(property_package=m.fs.properties, has_pressure_change=True)
        m.fs.unit.deltaP.fix(0)
        m.fs.unit.inlet.flow_mol[0].fix(5)
        m.fs.unit.inlet.enth_mol[0].fix(50000)
        m.fs.unit.inlet.pressure[0].fix(101325)

        m.fs.unit.heat_duty.fix(10000)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iapws):
        assert len(iapws.fs.unit.inlet.vars) == 3
        assert hasattr(iapws.fs.unit.inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.inlet, "pressure")

        assert hasattr(iapws.fs.unit, "outlet")
        assert len(iapws.fs.unit.outlet.vars) == 3
        assert hasattr(iapws.fs.unit.outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.outlet, "pressure")

        assert hasattr(iapws.fs.unit, "heat_duty")

        assert number_variables(iapws) == 8
        assert number_total_constraints(iapws) == 3
        assert number_unused_variables(iapws) == 0

    @pytest.mark.integration
    def test_units(self, iapws):
        assert_units_equivalent(
            iapws.fs.unit.control_volume.heat, pyunits.J / pyunits.s
        )
        assert_units_equivalent(iapws.fs.unit.heat_duty[0], pyunits.J / pyunits.s)
        assert_units_equivalent(iapws.fs.unit.deltaP[0], pyunits.Pa)
        assert_units_consistent(iapws)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

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
        assert pytest.approx(5, abs=1e-5) == value(iapws.fs.unit.outlet.flow_mol[0])

        assert pytest.approx(52000, abs=1e0) == value(iapws.fs.unit.outlet.enth_mol[0])

        assert pytest.approx(101325, abs=1e2) == value(iapws.fs.unit.outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert (
            abs(
                value(
                    iapws.fs.unit.inlet.flow_mol[0] - iapws.fs.unit.outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    iapws.fs.unit.inlet.flow_mol[0] * iapws.fs.unit.inlet.enth_mol[0]
                    - iapws.fs.unit.outlet.flow_mol[0]
                    * iapws.fs.unit.outlet.enth_mol[0]
                    + iapws.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict == {"vars": {"Heat Duty": iapws.fs.unit.heat_duty[0]}}

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_verify(self, iapws):
        # Test the heater model against known test cases
        cases = {
            "F": (1000, 1000, 1000, 1000),  # mol/s
            "Tin": (300, 300, 300, 800),  # K
            "Pin": (1000, 1000, 1000, 1000),  # kPa
            "xin": (0, 0, 0, 1),  # vapor fraction
            "Tout": (400, 400, 453.028, 300),  # K
            "Pout": (1000, 100, 1000, 1000),  # kPa
            "xout": (0, 1, 0.228898, 0),  # vapor fraction
            "duty": (7566.19, 47145, 20000, -61684.4),  # kW
        }

        for i in [0, 1, 2, 3]:
            F = cases["F"][i]
            Tin = cases["Tin"][i]
            Pin = cases["Pin"][i] * 1000
            hin = iapws95.htpx(T=Tin * pyunits.K, P=Pin * pyunits.Pa)
            Tout = cases["Tout"][i]
            Pout = cases["Pout"][i] * 1000
            xout = cases["xout"][i]
            duty = cases["duty"][i] * 1000
            prop_in = iapws.fs.unit.control_volume.properties_in[0]
            prop_out = iapws.fs.unit.control_volume.properties_out[0]

            prop_in.flow_mol = F
            prop_in.enth_mol = hin
            prop_in.pressure = Pin
            iapws.fs.unit.deltaP[0] = Pout - Pin
            iapws.fs.unit.heat_duty[0] = duty
            results = solver.solve(iapws)

            # Check for optimal solution
            assert check_optimal_termination(results)

            assert Tin == pytest.approx(value(prop_in.temperature), rel=1e-3)
            assert Tout == pytest.approx(value(prop_out.temperature), rel=1e-3)
            assert Pout == pytest.approx(value(prop_out.pressure), rel=1e-3)
            assert xout == pytest.approx(value(prop_out.vapor_frac), rel=1e-3)


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = Heater(property_package=m.fs.properties)

        m.fs.unit.inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.inlet.temperature[0].fix(320)
        m.fs.unit.inlet.pressure[0].fix(101325)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.heat_duty.fix(1000)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):
        assert len(sapon.fs.unit.inlet.vars) == 4
        assert hasattr(sapon.fs.unit.inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.inlet, "temperature")
        assert hasattr(sapon.fs.unit.inlet, "pressure")

        assert len(sapon.fs.unit.outlet.vars) == 4
        assert hasattr(sapon.fs.unit.outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.outlet, "temperature")
        assert hasattr(sapon.fs.unit.outlet, "pressure")

        assert hasattr(sapon.fs.unit, "heat_duty")

        assert number_variables(sapon) == 17
        assert number_total_constraints(sapon) == 8
        assert number_unused_variables(sapon) == 0

    @pytest.mark.integration
    def test_units(self, sapon):
        assert_units_equivalent(
            sapon.fs.unit.control_volume.heat, pyunits.J / pyunits.s
        )
        assert_units_equivalent(sapon.fs.unit.heat_duty[0], pyunits.J / pyunits.s)
        assert_units_consistent(sapon)

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
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(1e-3, abs=1e-6) == value(sapon.fs.unit.outlet.flow_vol[0])

        assert 55388.0 == value(sapon.fs.unit.inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(sapon.fs.unit.inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(sapon.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"])
        assert 0.0 == value(sapon.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"])
        assert 0.0 == value(sapon.fs.unit.inlet.conc_mol_comp[0, "Ethanol"])

        assert pytest.approx(320.2, abs=1e-1) == value(
            sapon.fs.unit.outlet.temperature[0]
        )

        assert pytest.approx(101325, abs=1e2) == value(sapon.fs.unit.outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0] - sapon.fs.unit.outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    sapon.fs.unit.outlet.flow_vol[0]
                    * sapon.fs.properties.dens_mol
                    * sapon.fs.properties.cp_mol
                    * (
                        sapon.fs.unit.inlet.temperature[0]
                        - sapon.fs.unit.outlet.temperature[0]
                    )
                    + sapon.fs.unit.heat_duty[0]
                )
            )
            <= 1e-3
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {"vars": {"Heat Duty": sapon.fs.unit.heat_duty[0]}}


# -----------------------------------------------------------------------------
class TestBT_Generic(object):
    @pytest.fixture(scope="class")
    def btg(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = GenericParameterBlock(**configuration)

        m.fs.unit = Heater(property_package=m.fs.properties, has_pressure_change=True)

        m.fs.unit.inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.inlet.temperature[0].fix(365)  # K
        m.fs.unit.inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.heat_duty.fix(-5000)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btg):
        assert hasattr(btg.fs.unit, "inlet")
        assert len(btg.fs.unit.inlet.vars) == 4
        assert hasattr(btg.fs.unit.inlet, "flow_mol")
        assert hasattr(btg.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(btg.fs.unit.inlet, "temperature")
        assert hasattr(btg.fs.unit.inlet, "pressure")

        assert hasattr(btg.fs.unit, "outlet")
        assert len(btg.fs.unit.outlet.vars) == 4
        assert hasattr(btg.fs.unit.outlet, "flow_mol")
        assert hasattr(btg.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(btg.fs.unit.outlet, "temperature")
        assert hasattr(btg.fs.unit.outlet, "pressure")

        assert hasattr(btg.fs.unit, "heat_duty")
        assert hasattr(btg.fs.unit, "deltaP")

        assert number_variables(btg) == 94
        assert number_total_constraints(btg) == 57
        # Unused vars are density parameters
        assert number_unused_variables(btg) == 10

    @pytest.mark.integration
    def test_units(self, btg):
        assert_units_equivalent(btg.fs.unit.control_volume.heat, pyunits.J / pyunits.s)
        assert_units_equivalent(btg.fs.unit.heat_duty[0], pyunits.J / pyunits.s)
        assert_units_equivalent(btg.fs.unit.deltaP[0], pyunits.Pa)
        assert_units_consistent(btg)

    @pytest.mark.unit
    def test_dof(self, btg):
        assert degrees_of_freedom(btg) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btg):
        initialization_tester(btg)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btg):
        results = solver.solve(btg)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btg):
        assert pytest.approx(5, abs=1e-3) == value(btg.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(358.6, abs=1e-1) == value(
            btg.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(btg.fs.unit.outlet.pressure[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btg):
        assert (
            abs(value(btg.fs.unit.inlet.flow_mol[0] - btg.fs.unit.outlet.flow_mol[0]))
            <= 1e-6
        )

        assert (
            abs(
                value(
                    btg.fs.unit.inlet.flow_mol[0]
                    * btg.fs.unit.control_volume.properties_in[0].enth_mol
                    - btg.fs.unit.outlet.flow_mol[0]
                    * btg.fs.unit.control_volume.properties_out[0].enth_mol
                    + btg.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btg):
        perf_dict = btg.fs.unit._get_performance_contents()

        assert perf_dict == {"vars": {"Heat Duty": btg.fs.unit.heat_duty[0]}}
