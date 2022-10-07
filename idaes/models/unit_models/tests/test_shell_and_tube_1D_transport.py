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
Tests for Shell and Tube 1D unit model.

Author: Douglas Allan, Jaffer Ghouse
"""
import pytest
from io import StringIO

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    value,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

import idaes
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    useDefault,
)
from idaes.models.unit_models.shell_and_tube_1d import ShellAndTube1D as HX1D
from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop, EosType
from idaes.models.properties import iapws95

from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

# Imports to assemble BT-PR with different units
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
import idaes.models.properties.modular_properties.pure.RPP4 as RPP

from idaes.logger import DEBUG

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
class TestBTX_cocurrent(object):
    @pytest.fixture(scope="class")
    def hx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = GenericParameterBlock(
            **get_prop({"N2", "O2", "Ar", "H2O", "CO2"}, {"Vap"}, eos=EosType.IDEAL),
            doc="Air property parameters",
        )
        m.fs.properties.set_default_scaling("enth_mol_phase", 1e-1)
        m.fs.properties.set_default_scaling("pressure", 1e-5)
        m.fs.properties.set_default_scaling("temperature", 1e-2)
        m.fs.properties.set_default_scaling("flow_mol", 1E-1)
        m.fs.properties.set_default_scaling("flow_mol_phase", 1e-1)
        _mf_scale = {
            "Ar": 100,
            "O2": 10,
            "N2": 10,
            "H2O": 100,
            "CO2": 1000}
        for comp, s in _mf_scale.items():
            m.fs.properties.set_default_scaling("mole_frac_comp", s, index=comp)
            m.fs.properties.set_default_scaling("mole_frac_phase_comp", s, index=("Vap", comp))
            m.fs.properties.set_default_scaling("flow_mol_phase_comp", s * 1e-1, index=("Vap", comp))

        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            hot_side_name="Shell",
            cold_side_name="Tube",
            flow_type=HeatExchangerFlowPattern.cocurrent,
        )
        m.fs.unit.tube_reynolds_number = pyo.Var(
            m.fs.time,
            m.fs.unit.cold_side.length_domain,
            units=pyo.units.dimensionless,
            initialize=4000
        )

        @m.fs.unit.Constraint(m.fs.time, m.fs.unit.cold_side.length_domain)
        def tube_reynolds_number_eqn(b, t, x):
            return (
                    b.tube_reynolds_number[t, x]
                    == b.cold_side.properties[t, x].flow_mass_phase["Vap"] * b.tube_inner_diameter
                    / b.cold_side.properties[t, x].visc_d_phase["Vap"] / b.cold_side.area
            )

        @m.fs.unit.Constraint(m.fs.time, m.fs.unit.cold_side.length_domain)
        def tube_heat_transfer_coeff_eqn(b, t, x):
            return b.cold_side_heat_transfer_coefficient[t, x] == (
                    0.027 * b.cold_side.properties[t, x].therm_cond_phase["Vap"] / b.tube_inner_diameter
                    * b.tube_reynolds_number[t, x] ** 0.8 * b.cold_side.properties[t, x].number_prandtl_phase[
                        "Vap"] ** (1 / 3)
            )

        m.fs.unit.shell_reynolds_number = pyo.Var(
            m.fs.time,
            m.fs.unit.hot_side.length_domain,
            units=pyo.units.dimensionless,
            initialize=4000
        )

        @m.fs.unit.Constraint(m.fs.time, m.fs.unit.hot_side.length_domain)
        def shell_reynolds_number_eqn(b, t, x):
            return (
                    b.shell_reynolds_number[t, x]
                    == b.hot_side.properties[t, x].flow_mass_phase["Vap"]
                    * b.tube_outer_diameter / b.hot_side.area
                    / b.hot_side.properties[t, x].visc_d_phase["Vap"]
            )

        @m.fs.unit.Constraint(m.fs.time, m.fs.unit.hot_side.length_domain)
        def shell_heat_transfer_coeff_eqn(b, t, x):
            return b.hot_side_heat_transfer_coefficient[t, x] == (
                    0.6 * 0.254 * b.hot_side.properties[t, x].therm_cond_phase["Vap"] / b.tube_outer_diameter
                    * b.shell_reynolds_number[t, x] ** 0.632 * b.hot_side.properties[t, x].number_prandtl_phase[
                        "Vap"] ** (
                            1 / 3)
            )
        m.fs.unit.length.fix(4.85)

        m.fs.unit.shell_diameter.fix(0.4)
        m.fs.unit.tube_outer_diameter.fix(0.01167)
        m.fs.unit.tube_inner_diameter.fix(0.01067)
        m.fs.unit.number_of_tubes.fix(100)
        m.fs.unit.length.fix(4.85)
        #m.fs.unit.hot_side_heat_transfer_coefficient.fix(2000)
        #m.fs.unit.cold_side_heat_transfer_coefficient.fix(51000)

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
                iscale.constraint_scaling_transform(m.fs.unit.tube_reynolds_number_eqn[t, x], 1e-4)
                iscale.set_scaling_factor(m.fs.unit.cold_side_heat_transfer_coefficient[t, x], 1e-2)
                iscale.constraint_scaling_transform(m.fs.unit.tube_heat_transfer_coeff_eqn[t, x], 1e-2)
            for x in m.fs.unit.hot_side.length_domain:
                iscale.set_scaling_factor(m.fs.unit.hot_side.heat[t, x], 1e-3)
                iscale.set_scaling_factor(m.fs.unit.shell_reynolds_number[t, x], 1e-3)
                iscale.constraint_scaling_transform(m.fs.unit.shell_reynolds_number_eqn[t, x], 1e-3)
                iscale.set_scaling_factor(m.fs.unit.hot_side_heat_transfer_coefficient[t, x], 1e-2)
                iscale.constraint_scaling_transform(m.fs.unit.shell_heat_transfer_coeff_eqn[t, x], 1e-2)

        iscale.calculate_scaling_factors(m)

        return m

    @pytest.mark.component
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, hx):
        hx.fs.unit.cold_side_heat_transfer_coefficient.fix(51000)
        hx.fs.unit.tube_reynolds_number.fix()
        hx.fs.unit.tube_reynolds_number_eqn.deactivate()
        hx.fs.unit.tube_heat_transfer_coeff_eqn.deactivate()
        hx.fs.unit.hot_side_heat_transfer_coefficient.fix(2000)
        hx.fs.unit.shell_reynolds_number.fix()
        hx.fs.unit.shell_reynolds_number_eqn.deactivate()
        hx.fs.unit.shell_heat_transfer_coeff_eqn.deactivate()

        hx.fs.unit.initialize_build(outlvl=DEBUG)

        hx.fs.unit.cold_side_heat_transfer_coefficient.unfix()
        hx.fs.unit.tube_reynolds_number.unfix()
        hx.fs.unit.tube_reynolds_number_eqn.activate()
        hx.fs.unit.tube_heat_transfer_coeff_eqn.activate()
        hx.fs.unit.hot_side_heat_transfer_coefficient.unfix()
        hx.fs.unit.shell_reynolds_number.unfix()
        hx.fs.unit.shell_reynolds_number_eqn.activate()
        hx.fs.unit.shell_heat_transfer_coeff_eqn.activate()
        for t in hx.fs.time:
            for x in hx.fs.unit.cold_side.length_domain:
                calculate_variable_from_constraint(
                    hx.fs.unit.tube_reynolds_number[t,x],
                    hx.fs.unit.tube_reynolds_number_eqn[t,x]
                )
                calculate_variable_from_constraint(
                    hx.fs.unit.cold_side_heat_transfer_coefficient[t,x],
                    hx.fs.unit.tube_heat_transfer_coeff_eqn[t,x]
                )
            for x in hx.fs.unit.hot_side.length_domain:
                calculate_variable_from_constraint(
                    hx.fs.unit.shell_reynolds_number[t,x],
                    hx.fs.unit.shell_reynolds_number_eqn[t,x]
                )
                calculate_variable_from_constraint(
                    hx.fs.unit.hot_side_heat_transfer_coefficient[t,x],
                    hx.fs.unit.shell_heat_transfer_coeff_eqn[t,x]
                )
        results = solver.solve(hx, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)
        #initialization_tester(hx)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, hx):
        results = solver.solve(hx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, hx):
        hx.fs.unit.temperature_wall.display()
        assert pytest.approx(5, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(322.669, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(10, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(322.463, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
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
                hx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Liq"]
                - hx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Liq"]
            )
        )
        cold_side = value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                hx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Liq"]
                - hx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Liq"]
            )
        )
        assert abs(hot_side - cold_side) <= 1e-6


# -----------------------------------------------------------------------------
class Testhx_countercurrent(object):
    @pytest.fixture(scope="class")
    def hx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = hxParameterBlock(valid_phase="Liq")

        m.fs.unit = HX1D(
            shell_is_hot=False,
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_type=HeatExchangerFlowPattern.countercurrent,
        )

        m.fs.unit.length.fix(4.85)

        m.fs.unit.shell_diameter.fix(1.04)
        m.fs.unit.tube_outer_diameter.fix(0.01167)
        m.fs.unit.tube_inner_diameter.fix(0.01067)
        m.fs.unit.number_of_tubes.fix(10)
        m.fs.unit.length.fix(4.85)
        m.fs.unit.hot_side_heat_transfer_coefficient.fix(2000)
        m.fs.unit.cold_side_heat_transfer_coefficient.fix(51000)

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(10)  # mol/s
        m.fs.unit.cold_side_inlet.temperature[0].fix(300)  # K
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        iscale.calculate_scaling_factors(m.fs.unit)

        return m

    @pytest.mark.unit
    def test_build(self, hx):
        assert hx.fs.unit.Shell is hx.fs.unit.cold_side
        assert hx.fs.unit.Tube is hx.fs.unit.hot_side

        assert hasattr(hx.fs.unit, "Shell_inlet")
        assert len(hx.fs.unit.Shell_inlet.vars) == 4
        assert hasattr(hx.fs.unit.Shell_inlet, "flow_mol")
        assert hasattr(hx.fs.unit.Shell_inlet, "mole_frac_comp")
        assert hasattr(hx.fs.unit.Shell_inlet, "temperature")
        assert hasattr(hx.fs.unit.Shell_inlet, "pressure")

        assert hasattr(hx.fs.unit, "Tube_inlet")
        assert len(hx.fs.unit.Tube_inlet.vars) == 4
        assert hasattr(hx.fs.unit.Tube_inlet, "flow_mol")
        assert hasattr(hx.fs.unit.Tube_inlet, "mole_frac_comp")
        assert hasattr(hx.fs.unit.Tube_inlet, "temperature")
        assert hasattr(hx.fs.unit.Tube_inlet, "pressure")

        assert hasattr(hx.fs.unit, "Shell_outlet")
        assert len(hx.fs.unit.Shell_outlet.vars) == 4
        assert hasattr(hx.fs.unit.Shell_outlet, "flow_mol")
        assert hasattr(hx.fs.unit.Shell_outlet, "mole_frac_comp")
        assert hasattr(hx.fs.unit.Shell_outlet, "temperature")
        assert hasattr(hx.fs.unit.Shell_outlet, "pressure")

        assert hasattr(hx.fs.unit, "Tube_outlet")
        assert len(hx.fs.unit.Tube_outlet.vars) == 4
        assert hasattr(hx.fs.unit.Tube_outlet, "flow_mol")
        assert hasattr(hx.fs.unit.Tube_outlet, "mole_frac_comp")
        assert hasattr(hx.fs.unit.Tube_outlet, "temperature")
        assert hasattr(hx.fs.unit.Tube_outlet, "pressure")

        assert hasattr(hx.fs.unit, "shell_diameter")
        assert hasattr(hx.fs.unit, "tube_inner_diameter")
        assert hasattr(hx.fs.unit, "tube_outer_diameter")
        assert hasattr(hx.fs.unit, "number_of_tubes")
        assert hasattr(hx.fs.unit, "length")
        assert hasattr(hx.fs.unit, "tube_side_xsec_area_calc")
        assert hasattr(hx.fs.unit, "shell_side_xsec_area_calc")

        assert hasattr(hx.fs.unit, "hot_side_heat_transfer_coefficient")
        assert hasattr(hx.fs.unit, "cold_side_heat_transfer_coefficient")
        assert hasattr(hx.fs.unit, "temperature_wall")
        assert hasattr(hx.fs.unit, "hot_side_heat_transfer_eq")
        assert hasattr(hx.fs.unit, "cold_side_heat_transfer_eq")
        assert hasattr(hx.fs.unit, "heat_conservation")

        assert number_variables(hx) == 869
        assert number_total_constraints(hx) == 804
        assert number_unused_variables(hx) == 8

    @pytest.mark.integration
    def test_units(self, hx):
        assert_units_equivalent(hx.fs.unit.length, pyunits.m)
        assert_units_equivalent(hx.fs.unit.shell_diameter, pyunits.m)
        assert_units_equivalent(hx.fs.unit.tube_inner_diameter, pyunits.m)
        assert_units_equivalent(hx.fs.unit.tube_outer_diameter, pyunits.m)
        assert_units_equivalent(hx.fs.unit.number_of_tubes, pyunits.dimensionless)

        assert_units_equivalent(
            hx.fs.unit.hot_side_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(
            hx.fs.unit.cold_side_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(hx.fs.unit.temperature_wall, pyunits.K)

        assert_units_consistent(hx)

    @pytest.mark.unit
    def test_dof(self, hx):
        assert degrees_of_freedom(hx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, hx):
        perf_dict = hx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Length": hx.fs.unit.hot_side.length,
                "Shell Diameter": hx.fs.unit.shell_diameter,
                "Tube Inner Diameter": hx.fs.unit.tube_inner_diameter,
                "Tube Outer Diameter": hx.fs.unit.tube_outer_diameter,
                "Number of Tubes": hx.fs.unit.number_of_tubes,
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, hx):
        stable = hx.fs.unit._get_stream_table_contents()

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
            "Tube Inlet": {
                "flow_mol": pytest.approx(5.0, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(365, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Tube Outlet": {
                "flow_mol": pytest.approx(1, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(298.15, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Shell Inlet": {
                "flow_mol": pytest.approx(10.0, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(300, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Shell Outlet": {
                "flow_mol": pytest.approx(1, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(298.15, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, hx):
        initialization_tester(hx)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, hx):
        results = solver.solve(hx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, hx):
        assert pytest.approx(5, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(304.292, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(10, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(331.436, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
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
                hx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Liq"]
                - hx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Liq"]
            )
        )
        cold_side = value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                hx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Liq"]
                - hx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Liq"]
            )
        )
        assert abs(hot_side - cold_side) <= 1e-6


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS_cocurrent(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock(
            phase_presentation=iapws95.PhaseType.LG
        )

        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_type=HeatExchangerFlowPattern.cocurrent,
        )

        m.fs.unit.length.fix(4.85)

        m.fs.unit.shell_diameter.fix(1.04)
        m.fs.unit.tube_outer_diameter.fix(0.01167)
        m.fs.unit.tube_inner_diameter.fix(0.01067)
        m.fs.unit.number_of_tubes.fix(10)
        m.fs.unit.length.fix(4.85)
        m.fs.unit.hot_side_heat_transfer_coefficient.fix(2000)
        m.fs.unit.cold_side_heat_transfer_coefficient.fix(51000)

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)
        m.fs.unit.hot_side_inlet.enth_mol[0].fix(50000)
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(5)
        m.fs.unit.cold_side_inlet.enth_mol[0].fix(7000)
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)

        return m

    @pytest.mark.unit
    def test_build(self, iapws):
        assert len(iapws.fs.unit.hot_side_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.hot_side_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "hot_side_outlet")
        assert len(iapws.fs.unit.hot_side_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.hot_side_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.hot_side_outlet, "pressure")

        assert len(iapws.fs.unit.cold_side_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.cold_side_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "cold_side_outlet")
        assert len(iapws.fs.unit.cold_side_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.cold_side_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.cold_side_outlet, "pressure")

        assert hasattr(iapws.fs.unit, "shell_diameter")
        assert hasattr(iapws.fs.unit, "tube_inner_diameter")
        assert hasattr(iapws.fs.unit, "tube_outer_diameter")
        assert hasattr(iapws.fs.unit, "number_of_tubes")
        assert hasattr(iapws.fs.unit, "length")
        assert hasattr(iapws.fs.unit, "tube_side_xsec_area_calc")
        assert hasattr(iapws.fs.unit, "shell_side_xsec_area_calc")

        assert hasattr(iapws.fs.unit, "hot_side_heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "cold_side_heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "temperature_wall")
        assert hasattr(iapws.fs.unit, "hot_side_heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "cold_side_heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "heat_conservation")

        assert number_variables(iapws) == 617
        assert number_total_constraints(iapws) == 554
        assert number_unused_variables(iapws) == 10

    @pytest.mark.integration
    def test_units(self, iapws):
        assert_units_equivalent(iapws.fs.unit.length, pyunits.m)
        assert_units_equivalent(iapws.fs.unit.shell_diameter, pyunits.m)
        assert_units_equivalent(iapws.fs.unit.tube_inner_diameter, pyunits.m)
        assert_units_equivalent(iapws.fs.unit.tube_outer_diameter, pyunits.m)
        assert_units_equivalent(iapws.fs.unit.number_of_tubes, pyunits.dimensionless)

        assert_units_equivalent(
            iapws.fs.unit.hot_side_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(
            iapws.fs.unit.cold_side_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(iapws.fs.unit.temperature_wall, pyunits.K)

        assert_units_consistent(iapws)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Length": iapws.fs.unit.hot_side.length,
                "Shell Diameter": iapws.fs.unit.shell_diameter,
                "Tube Inner Diameter": iapws.fs.unit.tube_inner_diameter,
                "Tube Outer Diameter": iapws.fs.unit.tube_outer_diameter,
                "Number of Tubes": iapws.fs.unit.number_of_tubes,
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, iapws):
        stable = iapws.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Molar Flow (mol/s)": getattr(pyunits.pint_registry, "mole/second"),
                "Mass Flow (kg/s)": getattr(pyunits.pint_registry, "kg/second"),
                "T (K)": getattr(pyunits.pint_registry, "K"),
                "P (Pa)": getattr(pyunits.pint_registry, "Pa"),
                "Vapor Fraction": getattr(pyunits.pint_registry, "dimensionless"),
                "Molar Enthalpy (J/mol) Vap": getattr(pyunits.pint_registry, "J/mole"),
                "Molar Enthalpy (J/mol) Liq": getattr(pyunits.pint_registry, "J/mole"),
            },
            "Shell Inlet": {
                "Molar Flow (mol/s)": pytest.approx(5, rel=1e-4),
                "Mass Flow (kg/s)": pytest.approx(0.090076, rel=1e-4),
                "T (K)": pytest.approx(422.6, rel=1e-4),
                "P (Pa)": pytest.approx(101325, rel=1e-4),
                "Vapor Fraction": pytest.approx(1, abs=1e-4),
                "Molar Enthalpy (J/mol) Vap": pytest.approx(50000, rel=1e-4),
                "Molar Enthalpy (J/mol) Liq": pytest.approx(11342, rel=1e-4),
            },
            "Shell Outlet": {
                "Molar Flow (mol/s)": pytest.approx(1, rel=1e-4),
                "Mass Flow (kg/s)": pytest.approx(1.8015e-2, rel=1e-4),
                "T (K)": pytest.approx(286.34, rel=1e-4),
                "P (Pa)": pytest.approx(1e5, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy (J/mol) Vap": pytest.approx(2168.6, rel=1e-4),
                "Molar Enthalpy (J/mol) Liq": pytest.approx(1000, rel=1e-4),
            },
            "Tube Inlet": {
                "Molar Flow (mol/s)": pytest.approx(5, rel=1e-4),
                "Mass Flow (kg/s)": pytest.approx(0.090076, rel=1e-4),
                "T (K)": pytest.approx(365.88, rel=1e-4),
                "P (Pa)": pytest.approx(101325, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy (J/mol) Vap": pytest.approx(47926, rel=1e-4),
                "Molar Enthalpy (J/mol) Liq": pytest.approx(7000, rel=1e-4),
            },
            "Tube Outlet": {
                "Molar Flow (mol/s)": pytest.approx(1, rel=1e-4),
                "Mass Flow (kg/s)": pytest.approx(1.8015e-2, rel=1e-4),
                "T (K)": pytest.approx(286.34, rel=1e-4),
                "P (Pa)": pytest.approx(1e5, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy (J/mol) Vap": pytest.approx(2168.6, rel=1e-4),
                "Molar Enthalpy (J/mol) Liq": pytest.approx(1000, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.unit
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(5, rel=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(5, rel=1e-5) == value(
            iapws.fs.unit.cold_side_outlet.flow_mol[0]
        )

        assert pytest.approx(48200.5, rel=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.enth_mol[0]
        )
        assert pytest.approx(8799.463, rel=1e-5) == value(
            iapws.fs.unit.cold_side_outlet.enth_mol[0]
        )

        assert pytest.approx(101325, rel=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.pressure[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            iapws.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert (
            abs(
                value(
                    iapws.fs.unit.hot_side_inlet.flow_mol[0]
                    - iapws.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    iapws.fs.unit.cold_side_inlet.flow_mol[0]
                    - iapws.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            iapws.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                iapws.fs.unit.hot_side_inlet.enth_mol[0]
                - iapws.fs.unit.hot_side_outlet.enth_mol[0]
            )
        )
        cold_side = value(
            iapws.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                iapws.fs.unit.cold_side_inlet.enth_mol[0]
                - iapws.fs.unit.cold_side_outlet.enth_mol[0]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS_countercurrent(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock(
            phase_presentation=iapws95.PhaseType.LG
        )

        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_type=HeatExchangerFlowPattern.countercurrent,
        )

        m.fs.unit.length.fix(4.85)

        m.fs.unit.shell_diameter.fix(1.04)
        m.fs.unit.tube_outer_diameter.fix(0.01167)
        m.fs.unit.tube_inner_diameter.fix(0.01067)
        m.fs.unit.number_of_tubes.fix(10)
        m.fs.unit.length.fix(4.85)
        m.fs.unit.hot_side_heat_transfer_coefficient.fix(2000)
        m.fs.unit.cold_side_heat_transfer_coefficient.fix(51000)

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)
        m.fs.unit.hot_side_inlet.enth_mol[0].fix(50000)
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(5)
        m.fs.unit.cold_side_inlet.enth_mol[0].fix(7000)
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)

        return m

    @pytest.mark.unit
    def test_build(self, iapws):
        assert len(iapws.fs.unit.hot_side_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.hot_side_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "hot_side_outlet")
        assert len(iapws.fs.unit.hot_side_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.hot_side_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.hot_side_outlet, "pressure")

        assert len(iapws.fs.unit.cold_side_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.cold_side_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "cold_side_outlet")
        assert len(iapws.fs.unit.cold_side_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.cold_side_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.cold_side_outlet, "pressure")

        assert hasattr(iapws.fs.unit, "shell_diameter")
        assert hasattr(iapws.fs.unit, "tube_inner_diameter")
        assert hasattr(iapws.fs.unit, "tube_outer_diameter")
        assert hasattr(iapws.fs.unit, "number_of_tubes")
        assert hasattr(iapws.fs.unit, "length")
        assert hasattr(iapws.fs.unit, "tube_side_xsec_area_calc")
        assert hasattr(iapws.fs.unit, "shell_side_xsec_area_calc")

        assert hasattr(iapws.fs.unit, "hot_side_heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "cold_side_heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "temperature_wall")
        assert hasattr(iapws.fs.unit, "hot_side_heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "cold_side_heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "heat_conservation")

        assert number_variables(iapws) == 617
        assert number_total_constraints(iapws) == 554
        assert number_unused_variables(iapws) == 10

    @pytest.mark.integration
    def test_units(self, iapws):
        assert_units_equivalent(iapws.fs.unit.length, pyunits.m)
        assert_units_equivalent(iapws.fs.unit.shell_diameter, pyunits.m)
        assert_units_equivalent(iapws.fs.unit.tube_inner_diameter, pyunits.m)
        assert_units_equivalent(iapws.fs.unit.tube_outer_diameter, pyunits.m)
        assert_units_equivalent(iapws.fs.unit.number_of_tubes, pyunits.dimensionless)

        assert_units_equivalent(
            iapws.fs.unit.hot_side_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(
            iapws.fs.unit.cold_side_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(iapws.fs.unit.temperature_wall, pyunits.K)

        assert_units_consistent(iapws)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Length": iapws.fs.unit.hot_side.length,
                "Shell Diameter": iapws.fs.unit.shell_diameter,
                "Tube Inner Diameter": iapws.fs.unit.tube_inner_diameter,
                "Tube Outer Diameter": iapws.fs.unit.tube_outer_diameter,
                "Number of Tubes": iapws.fs.unit.number_of_tubes,
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, iapws):
        stable = iapws.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Molar Flow (mol/s)": getattr(pyunits.pint_registry, "mole/second"),
                "Mass Flow (kg/s)": getattr(pyunits.pint_registry, "kg/second"),
                "T (K)": getattr(pyunits.pint_registry, "K"),
                "P (Pa)": getattr(pyunits.pint_registry, "Pa"),
                "Vapor Fraction": getattr(pyunits.pint_registry, "dimensionless"),
                "Molar Enthalpy (J/mol) Vap": getattr(pyunits.pint_registry, "J/mole"),
                "Molar Enthalpy (J/mol) Liq": getattr(pyunits.pint_registry, "J/mole"),
            },
            "Shell Inlet": {
                "Molar Flow (mol/s)": pytest.approx(5, rel=1e-4),
                "Mass Flow (kg/s)": pytest.approx(0.090076, rel=1e-4),
                "T (K)": pytest.approx(422.6, rel=1e-4),
                "P (Pa)": pytest.approx(101325, rel=1e-4),
                "Vapor Fraction": pytest.approx(1, abs=1e-4),
                "Molar Enthalpy (J/mol) Vap": pytest.approx(50000, rel=1e-4),
                "Molar Enthalpy (J/mol) Liq": pytest.approx(11342, rel=1e-4),
            },
            "Shell Outlet": {
                "Molar Flow (mol/s)": pytest.approx(1, rel=1e-4),
                "Mass Flow (kg/s)": pytest.approx(1.8015e-2, rel=1e-4),
                "T (K)": pytest.approx(286.34, rel=1e-4),
                "P (Pa)": pytest.approx(1e5, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy (J/mol) Vap": pytest.approx(2168.6, rel=1e-4),
                "Molar Enthalpy (J/mol) Liq": pytest.approx(1000, rel=1e-4),
            },
            "Tube Inlet": {
                "Molar Flow (mol/s)": pytest.approx(5, rel=1e-4),
                "Mass Flow (kg/s)": pytest.approx(0.090076, rel=1e-4),
                "T (K)": pytest.approx(365.88, rel=1e-4),
                "P (Pa)": pytest.approx(101325, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy (J/mol) Vap": pytest.approx(47926, rel=1e-4),
                "Molar Enthalpy (J/mol) Liq": pytest.approx(7000, rel=1e-4),
            },
            "Tube Outlet": {
                "Molar Flow (mol/s)": pytest.approx(1, rel=1e-4),
                "Mass Flow (kg/s)": pytest.approx(1.8015e-2, rel=1e-4),
                "T (K)": pytest.approx(286.34, rel=1e-4),
                "P (Pa)": pytest.approx(1e5, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy (J/mol) Vap": pytest.approx(2168.6, rel=1e-4),
                "Molar Enthalpy (J/mol) Liq": pytest.approx(1000, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(5, rel=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(5, rel=1e-5) == value(
            iapws.fs.unit.cold_side_outlet.flow_mol[0]
        )

        assert pytest.approx(47654.1, rel=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.enth_mol[0]
        )
        assert pytest.approx(9345.86, rel=1e-5) == value(
            iapws.fs.unit.cold_side_outlet.enth_mol[0]
        )

        assert pytest.approx(101325, rel=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.pressure[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            iapws.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert (
            abs(
                value(
                    iapws.fs.unit.hot_side_inlet.flow_mol[0]
                    - iapws.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    iapws.fs.unit.cold_side_inlet.flow_mol[0]
                    - iapws.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            iapws.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                iapws.fs.unit.hot_side_inlet.enth_mol[0]
                - iapws.fs.unit.hot_side_outlet.enth_mol[0]
            )
        )
        cold_side = value(
            iapws.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                iapws.fs.unit.cold_side_inlet.enth_mol[0]
                - iapws.fs.unit.cold_side_outlet.enth_mol[0]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


# -----------------------------------------------------------------------------
class TestBT_Generic_cocurrent(object):
    @pytest.fixture(scope="class")
    def hx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # As we lack other example prop packs with units, take the generic
        # BT-PR package and change the base units
        configuration2 = {
            # Specifying components
            "components": {
                "benzene": {
                    "type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "pressure_sat_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (78.1136e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                        "temperature_crit": (562.2, pyunits.K),  # [1]
                        "omega": 0.212,  # [1]
                        "cp_mol_ig_comp_coeff": {
                            "A": (-3.392e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                            "B": (4.739e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (-3.017e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (7.130e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "enth_mol_form_vap_comp_ref": (
                            82.9e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            -269,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),  # [3]
                        "pressure_sat_comp_coeff": {
                            "A": (-6.98273, None),  # [1]
                            "B": (1.33213, None),
                            "C": (-2.62863, None),
                            "D": (-3.33399, None),
                        },
                    },
                },
                "toluene": {
                    "type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "pressure_sat_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (92.1405e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (41e5, pyunits.Pa),  # [1]
                        "temperature_crit": (591.8, pyunits.K),  # [1]
                        "omega": 0.263,  # [1]
                        "cp_mol_ig_comp_coeff": {
                            "A": (-2.435e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                            "B": (5.125e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (-2.765e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (4.911e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "enth_mol_form_vap_comp_ref": (
                            50.1e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            -321,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),  # [3]
                        "pressure_sat_comp_coeff": {
                            "A": (-7.28607, None),  # [1]
                            "B": (1.38091, None),
                            "C": (-2.83433, None),
                            "D": (-2.79168, None),
                        },
                    },
                },
            },
            # Specifying phases
            "phases": {
                "Liq": {
                    "type": LiquidPhase,
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                },
                "Vap": {
                    "type": VaporPhase,
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                },
            },
            # Set base units of measurement
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.t,
                "amount": pyunits.mol,
                "temperature": pyunits.degR,
            },
            # Specifying state definition
            "state_definition": FTPx,
            "state_bounds": {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 500, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            "pressure_ref": (101325, pyunits.Pa),
            "temperature_ref": (298.15, pyunits.K),
            # Defining phase equilibria
            "phases_in_equilibrium": [("Vap", "Liq")],
            "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
            "bubble_dew_method": LogBubbleDew,
            "parameter_data": {
                "PR_kappa": {
                    ("benzene", "benzene"): 0.000,
                    ("benzene", "toluene"): 0.000,
                    ("toluene", "benzene"): 0.000,
                    ("toluene", "toluene"): 0.000,
                }
            },
        }

        m.fs.properties = GenericParameterBlock(**configuration)
        m.fs.properties2 = GenericParameterBlock(**configuration2)

        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties2},
            flow_type=HeatExchangerFlowPattern.cocurrent,
        )

        m.fs.unit.shell_diameter.fix(1.04)
        m.fs.unit.tube_outer_diameter.fix(0.01167)
        m.fs.unit.tube_inner_diameter.fix(0.01067)
        m.fs.unit.number_of_tubes.fix(10)
        m.fs.unit.length.fix(4.85)
        m.fs.unit.hot_side_heat_transfer_coefficient.fix(2000)
        m.fs.unit.cold_side_heat_transfer_coefficient.fix(51000)

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.cold_side_inlet.temperature[0].fix(540)  # degR
        m.fs.unit.cold_side_inlet.pressure[0].fix(101.325)  # kPa
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        return m

    @pytest.mark.component
    def test_build(self, hx):
        assert hasattr(hx.fs.unit, "hot_side_inlet")
        assert len(hx.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(hx.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(hx.fs.unit.hot_side_inlet, "mole_frac_comp")
        assert hasattr(hx.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(hx.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(hx.fs.unit, "cold_side_inlet")
        assert len(hx.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(hx.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(hx.fs.unit.cold_side_inlet, "mole_frac_comp")
        assert hasattr(hx.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(hx.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(hx.fs.unit, "hot_side_outlet")
        assert len(hx.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(hx.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(hx.fs.unit.hot_side_outlet, "mole_frac_comp")
        assert hasattr(hx.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(hx.fs.unit.hot_side_outlet, "pressure")

        assert hasattr(hx.fs.unit, "cold_side_outlet")
        assert len(hx.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(hx.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(hx.fs.unit.cold_side_outlet, "mole_frac_comp")
        assert hasattr(hx.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(hx.fs.unit.cold_side_outlet, "pressure")

        assert hasattr(hx.fs.unit, "shell_diameter")
        assert hasattr(hx.fs.unit, "tube_inner_diameter")
        assert hasattr(hx.fs.unit, "tube_outer_diameter")
        assert hasattr(hx.fs.unit, "number_of_tubes")
        assert hasattr(hx.fs.unit, "length")
        assert hasattr(hx.fs.unit, "tube_side_xsec_area_calc")
        assert hasattr(hx.fs.unit, "shell_side_xsec_area_calc")

        assert hasattr(hx.fs.unit, "hot_side_heat_transfer_coefficient")
        assert hasattr(hx.fs.unit, "cold_side_heat_transfer_coefficient")
        assert hasattr(hx.fs.unit, "temperature_wall")
        assert hasattr(hx.fs.unit, "hot_side_heat_transfer_eq")
        assert hasattr(hx.fs.unit, "cold_side_heat_transfer_eq")
        assert hasattr(hx.fs.unit, "heat_conservation")

        assert number_variables(hx) == 2021
        assert number_total_constraints(hx) == 1890
        assert number_unused_variables(hx) == 34

    @pytest.mark.integration
    def test_units(self, hx):
        assert_units_equivalent(hx.fs.unit.length, pyunits.m)
        assert_units_equivalent(hx.fs.unit.shell_diameter, pyunits.m)
        assert_units_equivalent(hx.fs.unit.tube_inner_diameter, pyunits.m)
        assert_units_equivalent(hx.fs.unit.tube_outer_diameter, pyunits.m)
        assert_units_equivalent(hx.fs.unit.number_of_tubes, pyunits.dimensionless)

        assert_units_equivalent(
            hx.fs.unit.hot_side_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(
            hx.fs.unit.cold_side_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(hx.fs.unit.temperature_wall, pyunits.K)

        assert_units_consistent(hx)

    @pytest.mark.component
    def test_dof(self, hx):
        assert degrees_of_freedom(hx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, hx):
        perf_dict = hx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Length": hx.fs.unit.hot_side.length,
                "Shell Diameter": hx.fs.unit.shell_diameter,
                "Tube Inner Diameter": hx.fs.unit.tube_inner_diameter,
                "Tube Outer Diameter": hx.fs.unit.tube_outer_diameter,
                "Number of Tubes": hx.fs.unit.number_of_tubes,
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, hx):
        stable = hx.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Total Molar Flowrate": getattr(pyunits.pint_registry, "mole/second"),
                "Total Mole Fraction benzene": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "Total Mole Fraction toluene": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "Temperature": getattr(pyunits.pint_registry, "kelvin"),
                "Pressure": getattr(pyunits.pint_registry, "Pa"),
            },
            "Shell Inlet": {
                "Total Molar Flowrate": pytest.approx(5.0, rel=1e-4),
                "Total Mole Fraction benzene": pytest.approx(0.5, rel=1e-4),
                "Total Mole Fraction toluene": pytest.approx(0.5, rel=1e-4),
                "Temperature": pytest.approx(365, rel=1e-4),
                "Pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Shell Outlet": {
                "Total Molar Flowrate": pytest.approx(100.0, rel=1e-4),
                "Total Mole Fraction benzene": pytest.approx(0.5, rel=1e-4),
                "Total Mole Fraction toluene": pytest.approx(0.5, rel=1e-4),
                "Temperature": pytest.approx(300, rel=1e-4),
                "Pressure": pytest.approx(1e5, rel=1e-4),
            },
            "Tube Inlet": {
                "Total Molar Flowrate": pytest.approx(1.0, rel=1e-4),
                "Total Mole Fraction benzene": pytest.approx(0.5, rel=1e-4),
                "Total Mole Fraction toluene": pytest.approx(0.5, rel=1e-4),
                "Temperature": pytest.approx(300, rel=1e-4),
                "Pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Tube Outlet": {
                "Total Molar Flowrate": pytest.approx(100.0, rel=1e-4),
                "Total Mole Fraction benzene": pytest.approx(0.5, rel=1e-4),
                "Total Mole Fraction toluene": pytest.approx(0.5, rel=1e-4),
                "Temperature": pytest.approx(300, rel=1e-4),
                "Pressure": pytest.approx(1e5, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_initialize(self, hx):
        initialization_tester(hx)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_solve(self, hx):
        results = solver.solve(hx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_solution(self, hx):
        # Note hot side in K and cold side in degR
        assert pytest.approx(5, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(354.878, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            hx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(1, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(638.780, rel=1e-5) == value(
            hx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101.325, rel=1e-5) == value(
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
                hx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Liq"]
                - hx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Liq"]
            )
        )
        cold_side = value(
            pyunits.convert(
                hx.fs.unit.cold_side_outlet.flow_mol[0]
                * (
                    hx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Liq"]
                    - hx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Liq"]
                ),
                to_units=pyunits.W,
            )
        )
        assert abs((hot_side - cold_side) / hot_side) <= 3e-4

    @pytest.mark.component
    def test_initialization_error(self, hx):
        hx.fs.unit.hot_side_outlet.flow_mol[0].fix(20)

        with idaes.temporary_config_ctx():
            with pytest.raises(InitializationError):
                hx.fs.unit.initialize(optarg={"max_iter": 1})
