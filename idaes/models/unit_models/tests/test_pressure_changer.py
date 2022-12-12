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
Tests for Pressure Changer unit model.

Author: Andrew Lee, Emmanuel Ogbe
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    units,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)

from idaes.models.unit_models.pressure_changer import (
    PressureChanger,
    PressureChangerData,
    Turbine,
    Compressor,
    Pump,
    ThermodynamicAssumption,
)

from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties import iapws95
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.util.exceptions import BalanceTypeNotSupportedError, InitializationError
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_ThermodynamicAssumption():
    assert len(ThermodynamicAssumption) == 4


class TestPressureChanger(object):
    @pytest.mark.unit
    def test_config(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = PressureChanger(property_package=m.fs.properties)

        # Check unit config arguments
        assert len(m.fs.unit.config) == 12

        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert not m.fs.unit.config.has_phase_equilibrium
        assert m.fs.unit.config.compressor
        assert (
            m.fs.unit.config.thermodynamic_assumption
            == ThermodynamicAssumption.isothermal
        )
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_dynamic_build(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_units=units.s)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = PressureChanger(property_package=m.fs.properties)
        iscale.calculate_scaling_factors(m)

        assert hasattr(m.fs.unit, "volume")

    @pytest.mark.unit
    def test_pump(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.pump,
        )
        iscale.calculate_scaling_factors(m)

        assert isinstance(m.fs.unit.fluid_work_calculation, Constraint)

    @pytest.mark.unit
    def test_adiabatic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.adiabatic,
        )
        iscale.calculate_scaling_factors(m)

    @pytest.mark.unit
    def test_isentropic_comp_phase_balances(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.isentropic,
            material_balance_type=MaterialBalanceType.componentPhase,
        )
        iscale.calculate_scaling_factors(m)

        assert isinstance(m.fs.unit.state_material_balances, Constraint)
        assert len(m.fs.unit.state_material_balances) == 4

    @pytest.mark.unit
    def test_isentropic_comp_total_balances(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.isentropic,
            material_balance_type=MaterialBalanceType.componentTotal,
        )
        iscale.calculate_scaling_factors(m)

        assert isinstance(m.fs.unit.state_material_balances, Constraint)
        assert len(m.fs.unit.state_material_balances) == 2

    @pytest.mark.unit
    def test_isentropic_total_balances(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        with pytest.raises(BalanceTypeNotSupportedError):
            m.fs.unit = PressureChanger(
                property_package=m.fs.properties,
                thermodynamic_assumption=ThermodynamicAssumption.isentropic,
                material_balance_type=MaterialBalanceType.total,
            )

    @pytest.mark.unit
    def test_isentropic_total_element_balances(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        with pytest.raises(BalanceTypeNotSupportedError):
            m.fs.unit = PressureChanger(
                property_package=m.fs.properties,
                thermodynamic_assumption=ThermodynamicAssumption.isentropic,
                material_balance_type=MaterialBalanceType.elementTotal,
            )

    @pytest.mark.unit
    def test_isentropic_material_balances_none(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        with pytest.raises(BalanceTypeNotSupportedError):
            m.fs.unit = PressureChanger(
                property_package=m.fs.properties,
                thermodynamic_assumption=ThermodynamicAssumption.isentropic,
                material_balance_type=MaterialBalanceType.none,
            )


# -----------------------------------------------------------------------------
class TestBTX_isothermal(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(valid_phase="Liq")

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.isothermal,
        )

        m.fs.unit.inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.inlet.temperature[0].fix(365)  # K
        m.fs.unit.inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.deltaP.fix(50000)
        iscale.calculate_scaling_factors(m)

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

        assert hasattr(btx.fs.unit, "work_mechanical")
        assert hasattr(btx.fs.unit, "deltaP")
        assert isinstance(btx.fs.unit.ratioP, Var)
        assert isinstance(btx.fs.unit.ratioP_calculation, Constraint)

        assert number_variables(btx) == 25
        assert number_total_constraints(btx) == 19
        assert number_unused_variables(btx) == 0

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_consistent(btx)
        assert_units_equivalent(btx.fs.unit.work_mechanical[0], units.W)
        assert_units_equivalent(btx.fs.unit.deltaP[0], units.Pa)

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
        assert pytest.approx(365, abs=1e-2) == value(btx.fs.unit.outlet.temperature[0])
        assert pytest.approx(151325, abs=1e2) == value(btx.fs.unit.outlet.pressure[0])
        assert pytest.approx(0, abs=1e-6) == value(btx.fs.unit.work_mechanical[0])

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
                    btx.fs.unit.outlet.flow_mol[0]
                    * (
                        btx.fs.unit.control_volume.properties_in[0].enth_mol_phase[
                            "Liq"
                        ]
                        - btx.fs.unit.control_volume.properties_out[0].enth_mol_phase[
                            "Liq"
                        ]
                    )
                )
            )
            <= 1e-6
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Mechanical Work": btx.fs.unit.work_mechanical[0],
                "Pressure Ratio": btx.fs.unit.ratioP[0],
                "Pressure Change": btx.fs.unit.deltaP[0],
            }
        }


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.isentropic,
            compressor=True,
        )

        m.fs.unit.inlet.flow_mol[0].fix(100)
        m.fs.unit.inlet.enth_mol[0].fix(4000)
        m.fs.unit.inlet.pressure[0].fix(101325)

        m.fs.unit.deltaP.fix(50000)
        m.fs.unit.efficiency_isentropic.fix(0.9)
        iscale.calculate_scaling_factors(m)

        return m

    @pytest.fixture(scope="class")
    def iapws_turb(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.isentropic,
            compressor=False,
        )
        iscale.calculate_scaling_factors(m)
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

        assert hasattr(iapws.fs.unit, "work_mechanical")
        assert hasattr(iapws.fs.unit, "deltaP")
        assert isinstance(iapws.fs.unit.ratioP, Var)
        assert isinstance(iapws.fs.unit.ratioP_calculation, Constraint)

        assert isinstance(iapws.fs.unit.efficiency_isentropic, Var)
        assert isinstance(iapws.fs.unit.work_isentropic, Var)

        assert hasattr(iapws.fs.unit, "properties_isentropic")
        assert isinstance(iapws.fs.unit.isentropic_pressure, Constraint)
        assert isinstance(iapws.fs.unit.state_material_balances, Constraint)
        assert isinstance(iapws.fs.unit.isentropic, Constraint)
        assert isinstance(iapws.fs.unit.isentropic_energy_balance, Constraint)
        assert isinstance(iapws.fs.unit.actual_work, Constraint)

        assert number_variables(iapws) == 14
        assert number_total_constraints(iapws) == 9
        assert number_unused_variables(iapws) == 0

    @pytest.mark.component
    def test_units(self, iapws):
        assert_units_consistent(iapws)
        assert_units_equivalent(iapws.fs.unit.work_mechanical[0], units.W)
        assert_units_equivalent(iapws.fs.unit.deltaP[0], units.Pa)

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
        # Check that outlet and isentropic pressure are equal
        assert pytest.approx(
            value(iapws.fs.unit.properties_isentropic[0].pressure), 1e-6
        ) == value(iapws.fs.unit.outlet.pressure[0])

        # Check that inlet and isentropic entropies are equal
        assert pytest.approx(
            value(iapws.fs.unit.properties_isentropic[0].entr_mol), 1e-6
        ) == value(iapws.fs.unit.control_volume.properties_in[0].entr_mol)

        assert pytest.approx(100, abs=1e-5) == value(iapws.fs.unit.outlet.flow_mol[0])

        assert pytest.approx(4002, abs=1e0) == value(iapws.fs.unit.outlet.enth_mol[0])

        assert pytest.approx(151325, abs=1e2) == value(iapws.fs.unit.outlet.pressure[0])

        assert pytest.approx(101.43796915073504, abs=1e-1) == value(
            iapws.fs.unit.work_mechanical[0]
        )

        assert pytest.approx(91.29417223566153, abs=1e-1) == value(
            iapws.fs.unit.work_isentropic[0]
        )

        # For verification, check outlet and isentropic temperatures
        assert pytest.approx(326.170, 1e-5) == value(
            iapws.fs.unit.control_volume.properties_out[0].temperature
        )

        assert pytest.approx(326.170, 1e-5) == value(
            iapws.fs.unit.properties_isentropic[0].temperature
        )

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
                    iapws.fs.unit.outlet.flow_mol[0]
                    * (
                        iapws.fs.unit.inlet.enth_mol[0]
                        - iapws.fs.unit.outlet.enth_mol[0]
                    )
                    + iapws.fs.unit.work_mechanical[0]
                )
            )
            <= 1e-6
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_verify(self, iapws_turb):
        iapws = iapws_turb
        # Verify the turbine results against 3 known test cases
        # Case Data (90% isentropic efficency)
        cases = {
            "F": (1000, 1000, 1000),  # mol/s
            "Tin": (500, 800, 400),  # K
            "Pin": (1000, 10000, 200),  # kPa
            "W": (-1224.64, -1911.55, -1010.64),  # kW
            "Tout": (463.435, 742.992, 382.442),  # K
            "Pout": (700, 7000, 140),  # kPa
            "xout": (1.0, 1.0, 0.9886),  # vapor fraction
            "Tisen": (460.149, 738.224, 382.442),
        }

        for i in [0, 1, 2]:
            F = cases["F"][i]
            Tin = cases["Tin"][i]
            Tout = cases["Tout"][i]
            Pin = cases["Pin"][i] * 1000
            Pout = cases["Pout"][i] * 1000
            hin = iapws95.htpx(T=Tin * units.K, P=Pin * units.Pa)
            W = cases["W"][i] * 1000
            Tis = cases["Tisen"][i]
            xout = cases["xout"][i]

            iapws.fs.unit.inlet.flow_mol[0].fix(F)
            iapws.fs.unit.inlet.enth_mol[0].fix(hin)
            iapws.fs.unit.inlet.pressure[0].fix(Pin)
            iapws.fs.unit.deltaP.fix(Pout - Pin)
            iapws.fs.unit.efficiency_isentropic.fix(0.9)
            iapws.fs.unit.initialize(optarg={"tol": 1e-6})
            results = solver.solve(iapws)
            # Check for optimal solution
            assert check_optimal_termination(results)

            Tout = pytest.approx(cases["Tout"][i], rel=1e-2)
            Pout = pytest.approx(cases["Pout"][i] * 1000, rel=1e-2)
            W = pytest.approx(cases["W"][i] * 1000, rel=1e-2)
            xout = pytest.approx(xout, rel=1e-2)
            prop_out = iapws.fs.unit.control_volume.properties_out[0]
            prop_in = iapws.fs.unit.control_volume.properties_in[0]
            prop_is = iapws.fs.unit.properties_isentropic[0]

            assert value(prop_in.temperature) == pytest.approx(Tin, rel=1e-3)
            assert value(prop_is.temperature) == pytest.approx(Tis, rel=1e-3)
            assert value(iapws.fs.unit.control_volume.work[0]) == W
            assert value(prop_out.pressure) == Pout
            assert value(prop_out.temperature) == Tout
            assert value(prop_out.vapor_frac) == xout

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Mechanical Work": iapws.fs.unit.work_mechanical[0],
                "Pressure Ratio": iapws.fs.unit.ratioP[0],
                "Pressure Change": iapws.fs.unit.deltaP[0],
                "Isentropic Efficiency": iapws.fs.unit.efficiency_isentropic[0],
            }
        }

    @pytest.mark.component
    def test_initialization_error(self, iapws):
        iapws.fs.unit.outlet.flow_mol[0].fix(500)

        with pytest.raises(InitializationError):
            iapws.fs.unit.initialize()


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.pump,
            compressor=False,
        )

        m.fs.unit.inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.inlet.temperature[0].fix(320)
        m.fs.unit.inlet.pressure[0].fix(101325)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.deltaP.fix(-20000)
        m.fs.unit.efficiency_pump.fix(0.9)
        iscale.calculate_scaling_factors(m)

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

        assert hasattr(sapon.fs.unit, "work_mechanical")
        assert hasattr(sapon.fs.unit, "deltaP")
        assert isinstance(sapon.fs.unit.ratioP, Var)
        assert isinstance(sapon.fs.unit.ratioP_calculation, Constraint)

        assert isinstance(sapon.fs.unit.efficiency_pump, Var)
        assert isinstance(sapon.fs.unit.work_fluid, Var)
        assert isinstance(sapon.fs.unit.fluid_work_calculation, Constraint)
        assert isinstance(sapon.fs.unit.actual_work, Constraint)

        assert number_variables(sapon) == 21
        assert number_total_constraints(sapon) == 11
        assert number_unused_variables(sapon) == 0

    @pytest.mark.component
    def test_units(self, sapon):
        assert_units_consistent(sapon)
        assert_units_equivalent(sapon.fs.unit.work_mechanical[0], units.W)
        assert_units_equivalent(sapon.fs.unit.deltaP[0], units.Pa)

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

        assert pytest.approx(55388.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]
        )

        assert pytest.approx(320.0, abs=1e-1) == value(
            sapon.fs.unit.outlet.temperature[0]
        )

        assert pytest.approx(81325, abs=1e2) == value(sapon.fs.unit.outlet.pressure[0])

        assert pytest.approx(-18.0, abs=1e-2) == value(sapon.fs.unit.work_mechanical[0])
        assert pytest.approx(-20.0, abs=1e-2) == value(sapon.fs.unit.work_fluid[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
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
                    + sapon.fs.unit.work_mechanical[0]
                )
            )
            <= 1e-4
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Mechanical Work": sapon.fs.unit.work_mechanical[0],
                "Pressure Ratio": sapon.fs.unit.ratioP[0],
                "Pressure Change": sapon.fs.unit.deltaP[0],
                "Efficiency": sapon.fs.unit.efficiency_pump[0],
            }
        }


class TestTurbine(object):
    @pytest.mark.unit
    def test_config(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = Turbine(property_package=m.fs.properties)
        iscale.calculate_scaling_factors(m)

        assert isinstance(m.fs.unit, PressureChangerData)
        # Check unit config arguments
        assert len(m.fs.unit.config) == 12

        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert not m.fs.unit.config.has_phase_equilibrium
        assert not m.fs.unit.config.compressor
        assert (
            m.fs.unit.config.thermodynamic_assumption
            == ThermodynamicAssumption.isentropic
        )
        assert m.fs.unit.config.property_package is m.fs.properties

        assert_units_consistent(m.fs.unit)


class TestCompressor(object):
    @pytest.mark.unit
    def test_config(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = Compressor(property_package=m.fs.properties)
        iscale.calculate_scaling_factors(m)

        assert isinstance(m.fs.unit, PressureChangerData)
        # Check unit config arguments
        assert len(m.fs.unit.config) == 12

        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert not m.fs.unit.config.has_phase_equilibrium
        assert m.fs.unit.config.compressor
        assert (
            m.fs.unit.config.thermodynamic_assumption
            == ThermodynamicAssumption.isentropic
        )
        assert m.fs.unit.config.property_package is m.fs.properties

        assert_units_consistent(m.fs.unit)


class TestPump(object):
    @pytest.mark.unit
    def test_config(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = Pump(property_package=m.fs.properties)
        iscale.calculate_scaling_factors(m)

        assert isinstance(m.fs.unit, PressureChangerData)
        # Check unit config arguments
        assert len(m.fs.unit.config) == 12

        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert not m.fs.unit.config.has_phase_equilibrium
        assert m.fs.unit.config.compressor
        assert m.fs.unit.config.thermodynamic_assumption == ThermodynamicAssumption.pump
        assert m.fs.unit.config.property_package is m.fs.properties

        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_pump_work_term_added_w_energybalancetype_none(self):
        # Check that work term is created when energy balance type none
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = Pump(
            property_package=m.fs.properties, energy_balance_type=EnergyBalanceType.none
        )

        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.none
        assert hasattr(m.fs.unit.control_volume, "work")
        assert hasattr(m.fs.unit, "work_mechanical")

    @pytest.mark.unit
    def test_pressure_changer_work_term_added_w_energybalancetype_none(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = PressureChanger(
            property_package=m.fs.properties,
            thermodynamic_assumption=ThermodynamicAssumption.pump,
            energy_balance_type=EnergyBalanceType.none,
        )
        assert hasattr(m.fs.unit.control_volume, "work")
        assert hasattr(m.fs.unit, "work_mechanical")


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
class TestTurbinePerformance(object):
    @pytest.mark.component
    def test_turbine_performance_way1(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = iapws95.Iapws95ParameterBlock()

        def perf_callback(b):
            unit_hd = units.J / units.kg
            unit_vflw = units.m**3 / units.s

            @b.Constraint(m.fs.time)
            def pc_isen_eff_eqn(b, t):
                prnt = b.parent_block()
                vflw = prnt.control_volume.properties_in[t].flow_vol
                return prnt.efficiency_isentropic[t] == 0.9 / 3.975 * vflw / unit_vflw

            @b.Constraint(m.fs.time)
            def pc_isen_head_eqn(b, t):
                prnt = b.parent_block()
                vflw = prnt.control_volume.properties_in[t].flow_vol
                return (
                    b.head_isentropic[t] / 1000
                    == -75530.8 / 3.975 / 1000 * vflw / unit_vflw * unit_hd
                )

        m.fs.unit = Turbine(
            property_package=m.fs.properties,
            support_isentropic_performance_curves=True,
            isentropic_performance_curves={"build_callback": perf_callback},
        )

        # set inputs
        m.fs.unit.inlet.flow_mol[0].fix(1000)  # mol/s
        Tin = 500  # K
        Pin = 1000000  # Pa
        Pout = 700000  # Pa
        hin = iapws95.htpx(Tin * units.K, Pin * units.Pa)
        m.fs.unit.inlet.enth_mol[0].fix(hin)
        m.fs.unit.inlet.pressure[0].fix(Pin)

        m.fs.unit.initialize()
        assert degrees_of_freedom(m) == 0
        assert_units_consistent(m.fs.unit)
        solver.solve(m, tee=True)

        assert value(m.fs.unit.efficiency_isentropic[0]) == pytest.approx(0.9, rel=1e-3)
        assert value(m.fs.unit.deltaP[0]) == pytest.approx(-3e5, rel=1e-3)

    @pytest.mark.component
    def test_turbine_performance_way2(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = iapws95.Iapws95ParameterBlock()
        m.fs.unit = Turbine(
            property_package=m.fs.properties, support_isentropic_performance_curves=True
        )
        unit_hd = units.J / units.kg
        unit_vflw = units.m**3 / units.s

        @m.fs.unit.performance_curve.Constraint(m.fs.time)
        def pc_isen_eff_eqn(b, t):
            prnt = b.parent_block()
            vflw = prnt.control_volume.properties_in[t].flow_vol
            return prnt.efficiency_isentropic[t] == 0.9 / 3.975 * vflw / unit_vflw

        @m.fs.unit.performance_curve.Constraint(m.fs.time)
        def pc_isen_head_eqn(b, t):
            prnt = b.parent_block()
            vflw = prnt.control_volume.properties_in[t].flow_vol
            return (
                b.head_isentropic[t] / 1000
                == -75530.8 / 3.975 / 1000 * vflw / unit_vflw * unit_hd
            )

        # set inputs
        m.fs.unit.inlet.flow_mol[0].fix(1000)  # mol/s
        Tin = 500  # K
        Pin = 1000000  # Pa
        Pout = 700000  # Pa
        hin = iapws95.htpx(Tin * units.K, Pin * units.Pa)
        m.fs.unit.inlet.enth_mol[0].fix(hin)
        m.fs.unit.inlet.pressure[0].fix(Pin)

        m.fs.unit.initialize()
        assert degrees_of_freedom(m) == 0
        assert_units_consistent(m.fs.unit)
        solver.solve(m, tee=True)

        assert value(m.fs.unit.efficiency_isentropic[0]) == pytest.approx(0.9, rel=1e-3)
        assert value(m.fs.unit.deltaP[0]) == pytest.approx(-3e5, rel=1e-3)
