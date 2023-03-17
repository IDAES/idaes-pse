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
Test in which viscosity and thermal conductivity are put to use in a heat transfer problem.
Chiefly tests that the properties are appearing on modular state blocks as expected.

Author: Douglas Allan, Jaffer Ghouse
"""
import pyomo.common.unittest as unittest

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    value,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import pyomo.environ as pyo

from idaes.core import (
    FlowsheetBlock,
)
from idaes.models.unit_models.shell_and_tube_1d import ShellAndTube1D as HX1D
from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)
from idaes.models.properties.modular_properties.pure.ChungPure import ChungViscosityPure


from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

from idaes.core.util.performance import PerformanceBaseClass

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


def build_model(eos, visc_d_phase_comp=None):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    comp_set = {"N2", "O2", "Ar", "H2O", "CO2"}
    config_dict = get_prop(comp_set, {"Vap"}, eos=eos)

    if visc_d_phase_comp is not None:
        # Wolfram Alpha
        dens_mass_crit_comp = {
            "N2": 313.3,
            "O2": 436.1,
            "Ar": 535.6,
            "H2O": 322.0,
            "CO2": 467.6,
        }
        for comp in comp_set:
            config_dict["components"][comp]["visc_d_phase_comp"][
                "Vap"
            ] = visc_d_phase_comp
            config_dict["components"][comp]["parameter_data"]["dipole_moment"] = (
                0.0,
                pyunits.debye,
            )
            config_dict["components"][comp]["parameter_data"][
                "association_factor_chung"
            ] = 0.0
            config_dict["components"][comp]["parameter_data"]["dens_mol_crit"] = (
                dens_mass_crit_comp[comp]
                / config_dict["components"][comp]["parameter_data"]["mw"][0],
                pyunits.mol / pyunits.m**3,
            )
        # From WolframAlpha
        config_dict["components"]["H2O"]["parameter_data"]["dipole_moment"] = (
            1.8546,
            pyunits.debye,
        )
        # From Table 9-1 Properties of Gases and Liquids 6th Ed.
        config_dict["components"]["H2O"]["parameter_data"][
            "association_factor_chung"
        ] = 0.076

    m.fs.properties = GenericParameterBlock(
        **config_dict,
        doc="Air property parameters",
    )
    m.fs.properties.set_default_scaling("enth_mol_phase", 1e-1)
    m.fs.properties.set_default_scaling("pressure", 1e-5)
    m.fs.properties.set_default_scaling("temperature", 1e-2)
    m.fs.properties.set_default_scaling("flow_mol", 1e-1)
    m.fs.properties.set_default_scaling("flow_mol_phase", 1e-1)
    _mf_scale = {"Ar": 100, "O2": 10, "N2": 10, "H2O": 100, "CO2": 1000}
    for comp, s in _mf_scale.items():
        m.fs.properties.set_default_scaling("mole_frac_comp", s, index=comp)
        m.fs.properties.set_default_scaling(
            "mole_frac_phase_comp", s, index=("Vap", comp)
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", s * 1e-1, index=("Vap", comp)
        )

    m.fs.unit = HX1D(
        hot_side={"property_package": m.fs.properties},
        cold_side={"property_package": m.fs.properties},
        hot_side_name="Shell",
        cold_side_name="Tube",
        flow_type=HeatExchangerFlowPattern.countercurrent,
    )
    m.fs.unit.tube_reynolds_number = pyo.Var(
        m.fs.time,
        m.fs.unit.cold_side.length_domain,
        units=pyo.units.dimensionless,
        initialize=4000,
    )

    @m.fs.unit.Constraint(m.fs.time, m.fs.unit.cold_side.length_domain)
    def tube_reynolds_number_eqn(b, t, x):
        return (
            b.tube_reynolds_number[t, x]
            == b.cold_side.properties[t, x].flow_mass_phase["Vap"]
            * b.tube_inner_diameter
            / b.cold_side.properties[t, x].visc_d_phase["Vap"]
            / b.cold_side.area
        )

    @m.fs.unit.Constraint(m.fs.time, m.fs.unit.cold_side.length_domain)
    def tube_heat_transfer_coeff_eqn(b, t, x):
        return b.cold_side_heat_transfer_coefficient[t, x] == (
            0.027
            * b.cold_side.properties[t, x].therm_cond_phase["Vap"]
            / b.tube_inner_diameter
            * b.tube_reynolds_number[t, x] ** 0.8
            * b.cold_side.properties[t, x].prandtl_number_phase["Vap"] ** (1 / 3)
        )

    m.fs.unit.shell_reynolds_number = pyo.Var(
        m.fs.time,
        m.fs.unit.hot_side.length_domain,
        units=pyo.units.dimensionless,
        initialize=4000,
    )

    @m.fs.unit.Constraint(m.fs.time, m.fs.unit.hot_side.length_domain)
    def shell_reynolds_number_eqn(b, t, x):
        return (
            b.shell_reynolds_number[t, x]
            == b.hot_side.properties[t, x].flow_mass_phase["Vap"]
            * b.tube_outer_diameter
            / b.hot_side.area
            / b.hot_side.properties[t, x].visc_d_phase["Vap"]
        )

    @m.fs.unit.Constraint(m.fs.time, m.fs.unit.hot_side.length_domain)
    def shell_heat_transfer_coeff_eqn(b, t, x):
        return b.hot_side_heat_transfer_coefficient[t, x] == (
            0.6
            * 0.254
            * b.hot_side.properties[t, x].therm_cond_phase["Vap"]
            / b.tube_outer_diameter
            * b.shell_reynolds_number[t, x] ** 0.632
            * b.hot_side.properties[t, x].prandtl_number_phase["Vap"] ** (1 / 3)
        )

    m.fs.unit.length.fix(4.85)

    m.fs.unit.shell_diameter.fix(0.4)
    m.fs.unit.tube_outer_diameter.fix(0.01167)
    m.fs.unit.tube_inner_diameter.fix(0.01067)
    m.fs.unit.number_of_tubes.fix(100)
    m.fs.unit.length.fix(4.85)

    m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
    m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
    m.fs.unit.hot_side_inlet.pressure[0].fix(101325)  # Pa
    m.fs.unit.hot_side_inlet.mole_frac_comp[0, "O2"].fix(0.2074)
    m.fs.unit.hot_side_inlet.mole_frac_comp[0, "H2O"].fix(0.0099)
    m.fs.unit.hot_side_inlet.mole_frac_comp[0, "CO2"].fix(0.0003)
    m.fs.unit.hot_side_inlet.mole_frac_comp[0, "N2"].fix(0.7732)
    m.fs.unit.hot_side_inlet.mole_frac_comp[0, "Ar"].fix(0.0092)

    m.fs.unit.cold_side_inlet.flow_mol[0].fix(10)  # mol/s
    m.fs.unit.cold_side_inlet.temperature[0].fix(300)  # K
    m.fs.unit.cold_side_inlet.pressure[0].fix(101325)  # Pa
    m.fs.unit.cold_side_inlet.mole_frac_comp[0, "O2"].fix(0.2074)
    m.fs.unit.cold_side_inlet.mole_frac_comp[0, "H2O"].fix(0.0099)
    m.fs.unit.cold_side_inlet.mole_frac_comp[0, "CO2"].fix(0.0003)
    m.fs.unit.cold_side_inlet.mole_frac_comp[0, "N2"].fix(0.7732)
    m.fs.unit.cold_side_inlet.mole_frac_comp[0, "Ar"].fix(0.0092)

    for t in m.fs.time:
        for x in m.fs.unit.cold_side.length_domain:
            iscale.set_scaling_factor(m.fs.unit.cold_side.heat[t, x], 1e-3)
            iscale.set_scaling_factor(m.fs.unit.tube_reynolds_number[t, x], 1e-4)
            iscale.constraint_scaling_transform(
                m.fs.unit.tube_reynolds_number_eqn[t, x], 1e-4
            )
            iscale.set_scaling_factor(
                m.fs.unit.cold_side_heat_transfer_coefficient[t, x], 1e-2
            )
            iscale.constraint_scaling_transform(
                m.fs.unit.tube_heat_transfer_coeff_eqn[t, x], 1e-2
            )
        for x in m.fs.unit.hot_side.length_domain:
            iscale.set_scaling_factor(m.fs.unit.hot_side.heat[t, x], 1e-3)
            iscale.set_scaling_factor(m.fs.unit.shell_reynolds_number[t, x], 1e-3)
            iscale.constraint_scaling_transform(
                m.fs.unit.shell_reynolds_number_eqn[t, x], 1e-3
            )
            iscale.set_scaling_factor(
                m.fs.unit.hot_side_heat_transfer_coefficient[t, x], 1e-2
            )
            iscale.constraint_scaling_transform(
                m.fs.unit.shell_heat_transfer_coeff_eqn[t, x], 1e-2
            )

    iscale.calculate_scaling_factors(m)

    return m


def initialize_model(m):
    m.fs.unit.cold_side_heat_transfer_coefficient.fix(51000)
    m.fs.unit.tube_reynolds_number.fix()
    m.fs.unit.tube_reynolds_number_eqn.deactivate()
    m.fs.unit.tube_heat_transfer_coeff_eqn.deactivate()
    m.fs.unit.hot_side_heat_transfer_coefficient.fix(2000)
    m.fs.unit.shell_reynolds_number.fix()
    m.fs.unit.shell_reynolds_number_eqn.deactivate()
    m.fs.unit.shell_heat_transfer_coeff_eqn.deactivate()

    m.fs.unit.initialize_build()

    m.fs.unit.cold_side_heat_transfer_coefficient.unfix()
    m.fs.unit.tube_reynolds_number.unfix()
    m.fs.unit.tube_reynolds_number_eqn.activate()
    m.fs.unit.tube_heat_transfer_coeff_eqn.activate()
    m.fs.unit.hot_side_heat_transfer_coefficient.unfix()
    m.fs.unit.shell_reynolds_number.unfix()
    m.fs.unit.shell_reynolds_number_eqn.activate()
    m.fs.unit.shell_heat_transfer_coeff_eqn.activate()
    for t in m.fs.time:
        for x in m.fs.unit.cold_side.length_domain:
            calculate_variable_from_constraint(
                m.fs.unit.tube_reynolds_number[t, x],
                m.fs.unit.tube_reynolds_number_eqn[t, x],
            )
            calculate_variable_from_constraint(
                m.fs.unit.cold_side_heat_transfer_coefficient[t, x],
                m.fs.unit.tube_heat_transfer_coeff_eqn[t, x],
            )
        for x in m.fs.unit.hot_side.length_domain:
            calculate_variable_from_constraint(
                m.fs.unit.shell_reynolds_number[t, x],
                m.fs.unit.shell_reynolds_number_eqn[t, x],
            )
            calculate_variable_from_constraint(
                m.fs.unit.hot_side_heat_transfer_coefficient[t, x],
                m.fs.unit.shell_heat_transfer_coeff_eqn[t, x],
            )
    return solver.solve(m, tee=True)


class Test_transport_properties_ideal(object):
    @pytest.fixture(scope="class")
    def hx(self):
        return build_model(eos=EosType.IDEAL)

    @pytest.mark.integration
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, hx):
        results = initialize_model(hx)
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_solution(self, hx):
        hx.fs.unit.temperature_wall.display()
        assert pytest.approx(5, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(313.390, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(10, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(325.767, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_conservation(self, hx):
        assert (
            abs(
                value(
                    hx.fs.unit.hot_side_inlet.flow_mol[0]
                    - hx.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    hx.fs.unit.cold_side_inlet.flow_mol[0]
                    - hx.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                hx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Vap"]
                - hx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Vap"]
            )
        )
        cold_side = value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                hx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Vap"]
                - hx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Vap"]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


class Test_transport_properties_ideal_chung(object):
    @pytest.fixture(scope="class")
    def hx(self):
        return build_model(eos=EosType.IDEAL, visc_d_phase_comp=ChungViscosityPure)

    @pytest.mark.integration
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, hx):
        results = initialize_model(hx)
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_solution(self, hx):
        hx.fs.unit.temperature_wall.display()
        assert pytest.approx(5, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(313.338, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(10, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(325.793, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_conservation(self, hx):
        assert (
            abs(
                value(
                    hx.fs.unit.hot_side_inlet.flow_mol[0]
                    - hx.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    hx.fs.unit.cold_side_inlet.flow_mol[0]
                    - hx.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                hx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Vap"]
                - hx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Vap"]
            )
        )
        cold_side = value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                hx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Vap"]
                - hx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Vap"]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


# -----------------------------------------------------------------------------
class Test_transport_properties_PR(object):
    @pytest.fixture(scope="class")
    def hx(self):
        return build_model(eos=EosType.PR)

    @pytest.mark.integration
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, hx):
        results = initialize_model(hx)
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_solution(self, hx):
        hx.fs.unit.temperature_wall.display()
        assert pytest.approx(5, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(313.405, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(10, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(325.752, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_conservation(self, hx):
        assert (
            abs(
                value(
                    hx.fs.unit.hot_side_inlet.flow_mol[0]
                    - hx.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    hx.fs.unit.cold_side_inlet.flow_mol[0]
                    - hx.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                hx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Vap"]
                - hx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Vap"]
            )
        )
        cold_side = value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                hx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Vap"]
                - hx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Vap"]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


@pytest.mark.performance
class TestCubicTransportPerformance(PerformanceBaseClass, unittest.TestCase):
    def build_model(self):
        return build_model(EosType.PR)

    def initialize_model(self, model):
        initialize_model(model)


if __name__ == "__main__":
    m = build_model(EosType.IDEAL, visc_d_phase_comp=ChungViscosityPure)
    initialize_model(m)
