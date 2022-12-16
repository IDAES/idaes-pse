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

import pytest

from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Objective,
    value,
    units as pyunits,
)

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from idaes.core.util.scaling import constraint_scaling_transform
from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.unit_models import GibbsReactor
from idaes.models.properties.modular_properties.eos.ceos import Cubic
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.logger as idaeslog

SOUT = idaeslog.INFO

# Set module level pyest marker
pytestmark = pytest.mark.cubic_root
prop_available = cubic_roots_available()

# Get default solver for testing
solver = get_solver()


class TestNaturalGasProps(object):
    @pytest.fixture()
    def m(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(
            **get_prop(
                components=[
                    "H2",
                    "CO",
                    "H2O",
                    "CO2",
                    "CH4",
                    "C2H6",
                    "C3H8",
                    "C4H10",
                    "N2",
                    "O2",
                    "Ar",
                ]
            )
        )
        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        return m

    @pytest.mark.integration
    def test_T_sweep(self, m):
        assert_units_consistent(m)

        m.fs.obj = Objective(expr=(m.fs.state[1].temperature - 510) ** 2)

        for logP in range(8, 13, 1):
            m.fs.state[1].flow_mol.fix(100)
            m.fs.state[1].mole_frac_comp["H2"].fix(0.1)
            m.fs.state[1].mole_frac_comp["CO"].fix(0.1)
            m.fs.state[1].mole_frac_comp["H2O"].fix(0.1)
            m.fs.state[1].mole_frac_comp["CO2"].fix(0.1)
            m.fs.state[1].mole_frac_comp["O2"].fix(0.1)
            m.fs.state[1].mole_frac_comp["N2"].fix(0.1)
            m.fs.state[1].mole_frac_comp["Ar"].fix(0.1)
            m.fs.state[1].mole_frac_comp["CH4"].fix(0.1)
            m.fs.state[1].mole_frac_comp["C2H6"].fix(0.1)
            m.fs.state[1].mole_frac_comp["C3H8"].fix(0.05)
            m.fs.state[1].mole_frac_comp["C4H10"].fix(0.05)
            m.fs.state[1].temperature.fix(300)
            m.fs.state[1].pressure.fix(10 ** (0.5 * logP))

            m.fs.state.initialize()

            m.fs.state[1].temperature.unfix()

            results = solver.solve(m, tee=True)

            assert check_optimal_termination(results)
            assert -93000 == pytest.approx(
                value(m.fs.state[1].enth_mol_phase["Vap"]), 1
            )
            assert 250 == pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1)

    @pytest.mark.integration
    def test_P_sweep(self, m):
        for T in range(300, 1000, 200):
            m.fs.state[1].flow_mol.fix(100)
            m.fs.state[1].mole_frac_comp["H2"].fix(0.1)
            m.fs.state[1].mole_frac_comp["CO"].fix(0.1)
            m.fs.state[1].mole_frac_comp["H2O"].fix(0.1)
            m.fs.state[1].mole_frac_comp["CO2"].fix(0.1)
            m.fs.state[1].mole_frac_comp["O2"].fix(0.1)
            m.fs.state[1].mole_frac_comp["N2"].fix(0.1)
            m.fs.state[1].mole_frac_comp["Ar"].fix(0.1)
            m.fs.state[1].mole_frac_comp["CH4"].fix(0.1)
            m.fs.state[1].mole_frac_comp["C2H6"].fix(0.1)
            m.fs.state[1].mole_frac_comp["C3H8"].fix(0.05)
            m.fs.state[1].mole_frac_comp["C4H10"].fix(0.05)
            m.fs.state[1].temperature.fix(T)
            m.fs.state[1].pressure.fix(1e5)

            m.fs.state.initialize()

            results = solver.solve(m)

            assert check_optimal_termination(results)

            while m.fs.state[1].pressure.value <= 1e6:
                m.fs.state[1].pressure.value = m.fs.state[1].pressure.value + 2e5
                results = solver.solve(m)

                assert check_optimal_termination(results)

                assert -70000 >= value(m.fs.state[1].enth_mol_phase["Vap"])
                assert -105000 <= value(m.fs.state[1].enth_mol_phase["Vap"])
                assert 185 <= value(m.fs.state[1].entr_mol_phase["Vap"])
                assert 265 >= value(m.fs.state[1].entr_mol_phase["Vap"])

    @pytest.mark.component
    def test_gibbs(self, m):
        m.fs.props = GenericParameterBlock(
            **get_prop(components=["H2", "CO", "H2O", "CO2", "O2", "N2", "Ar", "CH4"])
        )
        m.fs.reactor = GibbsReactor(
            dynamic=False,
            has_heat_transfer=True,
            has_pressure_change=False,
            property_package=m.fs.props,
        )

        m.fs.reactor.inlet.flow_mol.fix(8000)  # mol/s
        m.fs.reactor.inlet.temperature.fix(600)  # K
        m.fs.reactor.inlet.pressure.fix(1e6)  # Pa

        m.fs.reactor.outlet.temperature.fix(600)  # K

        m.fs.reactor.inlet.mole_frac_comp[0, "H2"].fix(0)
        m.fs.reactor.inlet.mole_frac_comp[0, "CO"].fix(0)
        m.fs.reactor.inlet.mole_frac_comp[0, "H2O"].fix(0.3)
        m.fs.reactor.inlet.mole_frac_comp[0, "CO2"].fix(0)
        m.fs.reactor.inlet.mole_frac_comp[0, "O2"].fix(0.2)
        m.fs.reactor.inlet.mole_frac_comp[0, "N2"].fix(0.05)
        m.fs.reactor.inlet.mole_frac_comp[0, "Ar"].fix(0.05)
        m.fs.reactor.inlet.mole_frac_comp[0, "CH4"].fix(0.4)

        constraint_scaling_transform(
            m.fs.reactor.control_volume.enthalpy_balances[0.0], 1e-6
        )
        m.fs.reactor.gibbs_scaling = 1e-6

        results = solver.solve(m, tee=True)

        assert check_optimal_termination(results)

        assert 0.017 == pytest.approx(
            value(m.fs.reactor.outlet.mole_frac_comp[0, "H2"]), 1e-1
        )
        assert 0.0001 == pytest.approx(
            value(m.fs.reactor.outlet.mole_frac_comp[0, "CO"]), 1
        )
        assert 0.487 == pytest.approx(
            value(m.fs.reactor.outlet.mole_frac_comp[0, "H2O"]), 1e-1
        )
        assert 0.103 == pytest.approx(
            value(m.fs.reactor.outlet.mole_frac_comp[0, "CO2"]), 1e-1
        )
        assert 0.293 == pytest.approx(
            value(m.fs.reactor.outlet.mole_frac_comp[0, "CH4"]), 1e-1
        )

        assert -634e6 == pytest.approx(value(m.fs.reactor.heat_duty[0]), 1e-3)


# Reference: CoolProp: http://www.coolprop.org/index.html
data = {
    "phase": "Vap",
    "pressure": 101325,
    "temperature": 311.8730783132979,
    "enth_mol_phase": 24188.517506218643,
    "heat_capacity_ratio_phase": 1.290590129342867,
    "cp_mol_phase": 37.726122878482165,
    "cv_mol_phase": 29.2316840341025,
    "isentropic_speed_sound_phase": 279.2959607849875,
    "mole_frac_phase_comp": {"CO2": 0.94, "H2O": 0.06},
}


class Test_CO2_H2O_Properties:
    @pytest.fixture(scope="class")
    def build_model(self):
        m = ConcreteModel()

        # Properties
        comp_props = get_prop(components=["CO2", "H2O"], phases=["Vap", "Liq"])

        # Parameters block
        m.params = GenericParameterBlock(**comp_props)

        m.props = m.params.build_state_block(
            [1], defined_state=True, parameters=m.params, has_phase_equilibrium=True
        )

        m.props[1].flow_mol.fix(100)
        m.props[1].pressure.fix(101325)
        m.props[1].mole_frac_comp["CO2"].fix(0.94)
        m.props[1].mole_frac_comp["H2O"].fix(0.06)
        m.props[1].temperature.fix(311.8730783132979)

        assert degrees_of_freedom(m) == 0

        m.props.initialize()
        results = get_solver(options={"bound_push": 1e-8}).solve(m)

        assert check_optimal_termination(results)

        return m

    @pytest.mark.unit
    def test_cp_mol_phase(self, build_model):
        m = build_model
        assert (
            str(pyunits.get_units(Cubic.cp_mol_phase(m.props[1], "Vap")))
            == "kg*m**2/K/mol/s**2"
        )

        assert (
            pytest.approx(value(Cubic.cp_mol_phase(m.props[1], "Vap")), rel=0.1)
            == data["cp_mol_phase"]
        )

    @pytest.mark.unit
    def test_cv_mol_phase(self, build_model):
        m = build_model
        assert (
            str(pyunits.get_units(Cubic.cv_mol_phase(m.props[1], "Vap")))
            == "kg*m**2/K/mol/s**2"
        )

        assert (
            pytest.approx(value(Cubic.cv_mol_phase(m.props[1], "Vap")), rel=0.1)
            == data["cv_mol_phase"]
        )

    @pytest.mark.unit
    def test_heat_capacity_ratio_phase(self, build_model):
        m = build_model
        assert (
            str(pyunits.get_units(Cubic.heat_capacity_ratio_phase(m.props[1], "Vap")))
            == "None"
        )

        assert (
            pytest.approx(
                value(Cubic.heat_capacity_ratio_phase(m.props[1], "Vap")), rel=0.1
            )
            == data["heat_capacity_ratio_phase"]
        )

    @pytest.mark.unit
    def test_isentropic_speed_sound_phase(self, build_model):
        m = build_model
        assert (
            str(
                pyunits.get_units(Cubic.isentropic_speed_sound_phase(m.props[1], "Vap"))
            )
            == "m/s"
        )

        assert (
            pytest.approx(
                value(Cubic.isentropic_speed_sound_phase(m.props[1], "Vap")), rel=0.1
            )
            == data["isentropic_speed_sound_phase"]
        )
