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
Tests for Heat Exchanger 1D unit model.

Author: Jaffer Ghouse
"""
import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    value,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
import pyomo.common.unittest as unittest

import idaes
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    useDefault,
)
from idaes.models.unit_models.heat_exchanger_1D import HeatExchanger1D as HX1D
from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.examples.BT_PR import configuration
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties import iapws95
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)

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
from idaes.core.util.performance import PerformanceBaseClass

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

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_bad_option():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(KeyError):
        m.fs.unit = HX1D(**{"I'm a bad option": "hot"})


@pytest.mark.unit
def test_bad_option2():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="cold_side_name cannot be 'hot_side'."
    ):
        m.fs.unit = HX1D(cold_side_name="hot_side")


@pytest.mark.unit
def test_bad_option3():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="hot_side_name cannot be 'cold_side'."
    ):
        m.fs.unit = HX1D(hot_side_name="cold_side")


@pytest.mark.unit
def test_bad_option4():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="cold_side_name cannot be 'cold_side'."
    ):
        m.fs.unit = HX1D(cold_side_name="cold_side")


@pytest.mark.unit
def test_bad_option5():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(ConfigurationError, match="hot_side_name cannot be 'hot_side'."):
        m.fs.unit = HX1D(hot_side_name="hot_side")


@pytest.mark.unit
def test_hot_and_cold_names_same():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        NameError,
        match="HeatExchanger hot and cold side cannot have the same name 'shell'.",
    ):
        m.fs.unit = HX1D(hot_side_name="shell", cold_side_name="shell")


@pytest.mark.unit
def test_hot_side_name_clash():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    with pytest.raises(
        ValueError,
        match="fs.unit could not assign hot side alias "
        "build as an attribute of that name already "
        "exists.",
    ):
        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            hot_side_name="build",
        )


@pytest.mark.unit
def test_cold_side_name_clash():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    with pytest.raises(
        ValueError,
        match="fs.unit could not assign cold side alias "
        "build as an attribute of that name already "
        "exists.",
    ):
        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            cold_side_name="build",
        )


@pytest.mark.unit
def test_user_names():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = HX1D(
        hot_side_name="shell",
        cold_side_name="tube",
        shell={"property_package": m.fs.properties},
        tube={"property_package": m.fs.properties},
    )

    assert m.fs.unit.config.hot_side.property_package is m.fs.properties
    assert m.fs.unit.config.cold_side.property_package is m.fs.properties

    assert m.fs.unit.shell is m.fs.unit.hot_side
    assert m.fs.unit.tube is m.fs.unit.cold_side

    assert m.fs.unit.shell_inlet is m.fs.unit.hot_side_inlet
    assert m.fs.unit.tube_inlet is m.fs.unit.cold_side_inlet
    assert m.fs.unit.shell_outlet is m.fs.unit.hot_side_outlet
    assert m.fs.unit.tube_outlet is m.fs.unit.cold_side_outlet


@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = HX1D(
        hot_side={"property_package": m.fs.properties},
        cold_side={"property_package": m.fs.properties},
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 9
    assert isinstance(m.fs.unit.config.hot_side, ConfigBlock)
    assert isinstance(m.fs.unit.config.cold_side, ConfigBlock)
    assert m.fs.unit.config.flow_type == HeatExchangerFlowPattern.cocurrent
    assert m.fs.unit.config.finite_elements == 20
    assert m.fs.unit.config.collocation_points == 5

    assert m.fs.unit.config.hot_side_name == None
    assert m.fs.unit.config.cold_side_name == None

    # Check hot side config arguments
    assert len(m.fs.unit.config.hot_side) == 11
    assert m.fs.unit.config.hot_side.dynamic == useDefault
    assert m.fs.unit.config.hot_side.has_holdup == useDefault
    assert (
        m.fs.unit.config.hot_side.material_balance_type
        == MaterialBalanceType.useDefault
    )
    assert m.fs.unit.config.hot_side.energy_balance_type == EnergyBalanceType.useDefault
    assert (
        m.fs.unit.config.hot_side.momentum_balance_type
        == MomentumBalanceType.pressureTotal
    )
    assert not m.fs.unit.config.hot_side.has_pressure_change
    assert not m.fs.unit.config.hot_side.has_phase_equilibrium
    assert m.fs.unit.config.hot_side.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.hot_side.transformation_scheme == "BACKWARD"

    # Check cold side config arguments
    assert len(m.fs.unit.config.cold_side) == 11
    assert m.fs.unit.config.cold_side.dynamic == useDefault
    assert m.fs.unit.config.cold_side.has_holdup == useDefault
    assert (
        m.fs.unit.config.cold_side.material_balance_type
        == MaterialBalanceType.useDefault
    )
    assert (
        m.fs.unit.config.cold_side.energy_balance_type == EnergyBalanceType.useDefault
    )
    assert (
        m.fs.unit.config.cold_side.momentum_balance_type
        == MomentumBalanceType.pressureTotal
    )
    assert not m.fs.unit.config.cold_side.has_pressure_change
    assert not m.fs.unit.config.cold_side.has_phase_equilibrium
    assert m.fs.unit.config.cold_side.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.cold_side.transformation_scheme == "BACKWARD"


@pytest.mark.unit
def test_config_validation():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = BTXParameterBlock(valid_phase="Liq")

    with pytest.raises(ConfigurationError):
        m.fs.HX_co_current = HX1D(
            hot_side={
                "property_package": m.fs.properties,
                "transformation_scheme": "BACKWARD",
            },
            cold_side={
                "property_package": m.fs.properties,
                "transformation_scheme": "FORWARD",
            },
            flow_type=HeatExchangerFlowPattern.cocurrent,
        )

    with pytest.raises(ConfigurationError):
        m.fs.HX_counter_current = HX1D(
            hot_side={
                "property_package": m.fs.properties,
                "transformation_method": "dae.finite_difference",
            },
            cold_side={
                "property_package": m.fs.properties,
                "transformation_method": "dae.collocation",
            },
            flow_type=HeatExchangerFlowPattern.countercurrent,
        )


# -----------------------------------------------------------------------------
class TestBTX_cocurrent(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(valid_phase="Liq")

        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            hot_side_name="Shell",
            cold_side_name="Tube",
            flow_type=HeatExchangerFlowPattern.cocurrent,
        )

        m.fs.unit.area.fix(0.5)
        m.fs.unit.length.fix(4.85)
        m.fs.unit.heat_transfer_coefficient.fix(500)

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.cold_side_inlet.temperature[0].fix(300)  # K
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        iscale.calculate_scaling_factors(m)

        return m

    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "Shell_inlet")
        assert len(btx.fs.unit.Shell_inlet.vars) == 4
        assert hasattr(btx.fs.unit.Shell_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.Shell_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.Shell_inlet, "temperature")
        assert hasattr(btx.fs.unit.Shell_inlet, "pressure")

        assert hasattr(btx.fs.unit, "Tube_inlet")
        assert len(btx.fs.unit.Tube_inlet.vars) == 4
        assert hasattr(btx.fs.unit.Tube_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.Tube_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.Tube_inlet, "temperature")
        assert hasattr(btx.fs.unit.Tube_inlet, "pressure")

        assert hasattr(btx.fs.unit, "hot_side_outlet")
        assert len(btx.fs.unit.Shell_outlet.vars) == 4
        assert hasattr(btx.fs.unit.Shell_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.Shell_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.Shell_outlet, "temperature")
        assert hasattr(btx.fs.unit.Shell_outlet, "pressure")

        assert hasattr(btx.fs.unit, "Tube_outlet")
        assert len(btx.fs.unit.Tube_outlet.vars) == 4
        assert hasattr(btx.fs.unit.Tube_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.Tube_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.Tube_outlet, "temperature")
        assert hasattr(btx.fs.unit.Tube_outlet, "pressure")

        assert hasattr(btx.fs.unit, "area")
        assert hasattr(btx.fs.unit, "length")
        assert hasattr(btx.fs.unit, "heat_transfer_coefficient")
        assert hasattr(btx.fs.unit, "heat_transfer_eq")
        assert hasattr(btx.fs.unit, "heat_conservation")

        assert number_variables(btx) == 824
        assert number_total_constraints(btx) == 781
        assert number_unused_variables(btx) == 10

    @pytest.mark.integration
    def test_units(self, btx):
        assert_units_equivalent(btx.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(btx.fs.unit.length, pyunits.m)

        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Area": btx.fs.unit.area,
                "Length": btx.fs.unit.hot_side.length,
            }
        }

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
            "Shell Inlet": {
                "flow_mol": pytest.approx(5.0, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(365, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Shell Outlet": {
                "flow_mol": pytest.approx(1, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(298.15, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Tube Inlet": {
                "flow_mol": pytest.approx(1.0, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(300, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Tube Outlet": {
                "flow_mol": pytest.approx(1, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(298.15, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.component
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert pytest.approx(5, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(356.414, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(1, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(346.06, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.hot_side_inlet.flow_mol[0]
                    - btx.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    btx.fs.unit.cold_side_inlet.flow_mol[0]
                    - btx.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                btx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Liq"]
                - btx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Liq"]
            )
        )
        cold_side = value(
            btx.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                btx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Liq"]
                - btx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Liq"]
            )
        )
        assert abs(hot_side - cold_side) <= 1e-6


# -----------------------------------------------------------------------------
class TestBTX_countercurrent(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(valid_phase="Liq")

        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_type=HeatExchangerFlowPattern.countercurrent,
        )

        m.fs.unit.length.fix(4.85)
        m.fs.unit.area.fix(0.5)
        m.fs.unit.heat_transfer_coefficient.fix(500)

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.cold_side_inlet.temperature[0].fix(300)  # K
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        iscale.calculate_scaling_factors(m.fs.unit)

        return m

    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "hot_side_inlet")
        assert len(btx.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_side_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(btx.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_side_inlet")
        assert len(btx.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_side_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(btx.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(btx.fs.unit, "hot_side_outlet")
        assert len(btx.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_side_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(btx.fs.unit.hot_side_outlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_side_outlet")
        assert len(btx.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_side_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(btx.fs.unit.cold_side_outlet, "pressure")

        assert hasattr(btx.fs.unit, "area")
        assert hasattr(btx.fs.unit, "length")
        assert hasattr(btx.fs.unit, "heat_transfer_coefficient")
        assert hasattr(btx.fs.unit, "heat_transfer_eq")
        assert hasattr(btx.fs.unit, "heat_conservation")

        assert number_variables(btx) == 824
        assert number_total_constraints(btx) == 781
        assert number_unused_variables(btx) == 10

    @pytest.mark.integration
    def test_units(self, btx):
        assert_units_equivalent(btx.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(btx.fs.unit.length, pyunits.m)
        assert_units_equivalent(
            btx.fs.unit.heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )

        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Area": btx.fs.unit.area,
                "Length": btx.fs.unit.hot_side.length,
            }
        }

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
            "Hot Side Inlet": {
                "flow_mol": pytest.approx(5.0, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(365, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Hot Side Outlet": {
                "flow_mol": pytest.approx(1, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(298.15, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Cold Side Inlet": {
                "flow_mol": pytest.approx(1.0, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(300, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Cold Side Outlet": {
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
    def test_initialize(self, btx):
        initialization_tester(
            btx,
            optarg={"tol": 1e-6},
            hot_side_state_args={"flow_mol": 5, "temperature": 304, "pressure": 101325},
            cold_side_state_args={
                "flow_mol": 1,
                "temperature": 331.5,
                "pressure": 101325,
            },
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert pytest.approx(5, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(355.505, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(1, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(350.67, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.hot_side_inlet.flow_mol[0]
                    - btx.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    btx.fs.unit.cold_side_inlet.flow_mol[0]
                    - btx.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                btx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Liq"]
                - btx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Liq"]
            )
        )
        cold_side = value(
            btx.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                btx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Liq"]
                - btx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Liq"]
            )
        )
        assert abs(hot_side - cold_side) <= 1e-6


# -----------------------------------------------------------------------------
def build_model():
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
    m.fs.unit.area.fix(0.5)
    m.fs.unit.heat_transfer_coefficient.fix(500)

    m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)
    m.fs.unit.hot_side_inlet.enth_mol[0].fix(50000)
    m.fs.unit.hot_side_inlet.pressure[0].fix(101325)

    m.fs.unit.cold_side_inlet.flow_mol[0].fix(5)
    m.fs.unit.cold_side_inlet.enth_mol[0].fix(7000)
    m.fs.unit.cold_side_inlet.pressure[0].fix(101325)

    return m


@pytest.mark.performance
class Test_HX1D_Performance(PerformanceBaseClass, unittest.TestCase):
    def build_model(self):
        return build_model()

    def initialize_model(self, model):
        model.fs.unit.initialize()


@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS_cocurrent(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        return build_model()

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

        assert hasattr(iapws.fs.unit, "area")
        assert hasattr(iapws.fs.unit, "length")
        assert hasattr(iapws.fs.unit, "heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "heat_conservation")

        assert number_variables(iapws) == 572
        assert number_total_constraints(iapws) == 531
        assert number_unused_variables(iapws) == 12

    @pytest.mark.integration
    def test_units(self, iapws):
        assert_units_equivalent(iapws.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(iapws.fs.unit.length, pyunits.m)
        assert_units_equivalent(
            iapws.fs.unit.heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )

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
                "Area": iapws.fs.unit.area,
                "Length": iapws.fs.unit.hot_side.length,
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
            "Hot Side Inlet": {
                "Molar Flow": pytest.approx(5, rel=1e-4),
                "Mass Flow": pytest.approx(0.090076, rel=1e-4),
                "T": pytest.approx(422.6, rel=1e-4),
                "P": pytest.approx(101325, rel=1e-4),
                "Vapor Fraction": pytest.approx(1, abs=1e-4),
                "Molar Enthalpy": pytest.approx(50000, rel=1e-4),
            },
            "Hot Side Outlet": {
                "Molar Flow": pytest.approx(1, rel=1e-4),
                "Mass Flow": pytest.approx(1.8015e-2, rel=1e-4),
                "T": pytest.approx(270.4877112932641, rel=1e-4),
                "P": pytest.approx(11032305.8275, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(0.01102138712926277, rel=1e-4),
            },
            "Cold Side Inlet": {
                "Molar Flow": pytest.approx(5, rel=1e-4),
                "Mass Flow": pytest.approx(0.090076, rel=1e-4),
                "T": pytest.approx(365.88, rel=1e-4),
                "P": pytest.approx(101325, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(7000.0, rel=1e-4),
            },
            "Cold Side Outlet": {
                "Molar Flow": pytest.approx(1, rel=1e-4),
                "Mass Flow": pytest.approx(1.8015e-2, rel=1e-4),
                "T": pytest.approx(270.4877112932641, rel=1e-4),
                "P": pytest.approx(11032305.8275, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(0.01102138712926277, rel=1e-4),
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
        assert_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(5, rel=1e-4) == value(
            iapws.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(5, rel=1e-4) == value(
            iapws.fs.unit.cold_side_outlet.flow_mol[0]
        )

        assert pytest.approx(48673.2, rel=1e-4) == value(
            iapws.fs.unit.hot_side_outlet.enth_mol[0]
        )
        assert pytest.approx(8326.77, rel=1e-4) == value(
            iapws.fs.unit.cold_side_outlet.enth_mol[0]
        )

        assert pytest.approx(101325, rel=1e-4) == value(
            iapws.fs.unit.hot_side_outlet.pressure[0]
        )
        assert pytest.approx(101325, rel=1e-4) == value(
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


# # -----------------------------------------------------------------------------
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
        m.fs.unit.area.fix(0.5)
        m.fs.unit.heat_transfer_coefficient.fix(500)

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

        assert hasattr(iapws.fs.unit, "area")
        assert hasattr(iapws.fs.unit, "length")
        assert hasattr(iapws.fs.unit, "heat_transfer_coefficient")
        assert hasattr(iapws.fs.unit, "heat_transfer_eq")
        assert hasattr(iapws.fs.unit, "heat_conservation")

        assert number_variables(iapws) == 572
        assert number_total_constraints(iapws) == 531
        assert number_unused_variables(iapws) == 12

    @pytest.mark.integration
    def test_units(self, iapws):
        assert_units_equivalent(iapws.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(iapws.fs.unit.length, pyunits.m)
        assert_units_equivalent(
            iapws.fs.unit.heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )

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
                "Area": iapws.fs.unit.area,
                "Length": iapws.fs.unit.hot_side.length,
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
            "Hot Side Inlet": {
                "Molar Flow": pytest.approx(5, rel=1e-4),
                "Mass Flow": pytest.approx(0.090076, rel=1e-4),
                "T": pytest.approx(422.6, rel=1e-4),
                "P": pytest.approx(101325, rel=1e-4),
                "Vapor Fraction": pytest.approx(1, abs=1e-4),
                "Molar Enthalpy": pytest.approx(50000, rel=1e-4),
            },
            "Hot Side Outlet": {
                "Molar Flow": pytest.approx(1, rel=1e-4),
                "Mass Flow": pytest.approx(1.8015e-2, rel=1e-4),
                "T": pytest.approx(270.4877112932641, rel=1e-4),
                "P": pytest.approx(11032305.8275, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(0.01102138712926277, rel=1e-4),
            },
            "Cold Side Inlet": {
                "Molar Flow": pytest.approx(5, rel=1e-4),
                "Mass Flow": pytest.approx(0.090076, rel=1e-4),
                "T": pytest.approx(365.88, rel=1e-4),
                "P": pytest.approx(101325, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(7000.0, rel=1e-4),
            },
            "Cold Side Outlet": {
                "Molar Flow": pytest.approx(1, rel=1e-4),
                "Mass Flow": pytest.approx(1.8015e-2, rel=1e-4),
                "T": pytest.approx(270.4877112932641, rel=1e-4),
                "P": pytest.approx(11032305.8275, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(0.01102138712926277, rel=1e-4),
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
        assert_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(5, rel=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(5, rel=1e-5) == value(
            iapws.fs.unit.cold_side_outlet.flow_mol[0]
        )

        assert pytest.approx(48599.9, rel=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.enth_mol[0]
        )
        assert pytest.approx(8400.11, rel=1e-5) == value(
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


# # -----------------------------------------------------------------------------
class TestSaponification_cocurrent(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_type=HeatExchangerFlowPattern.cocurrent,
        )

        m.fs.unit.length.fix(4.85)
        m.fs.unit.area.fix(0.5)
        m.fs.unit.heat_transfer_coefficient.fix(500)

        m.fs.unit.hot_side_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.hot_side_inlet.temperature[0].fix(320)
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.cold_side_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.cold_side_inlet.temperature[0].fix(300)
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        return m

    @pytest.mark.unit
    def test_build(self, sapon):
        assert len(sapon.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.hot_side_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "pressure")

        assert len(sapon.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.hot_side_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "pressure")

        assert len(sapon.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.cold_side_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "pressure")

        assert len(sapon.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.cold_side_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "pressure")

        assert hasattr(sapon.fs.unit, "area")
        assert hasattr(sapon.fs.unit, "length")
        assert hasattr(sapon.fs.unit, "heat_transfer_coefficient")
        assert hasattr(sapon.fs.unit, "heat_transfer_eq")
        assert hasattr(sapon.fs.unit, "heat_conservation")

        assert number_variables(sapon) == 950
        assert number_total_constraints(sapon) == 895
        assert number_unused_variables(sapon) == 16

    @pytest.mark.integration
    def test_units(self, sapon):
        assert_units_equivalent(sapon.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(sapon.fs.unit.length, pyunits.m)
        assert_units_equivalent(
            sapon.fs.unit.heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )

        assert_units_consistent(sapon)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Area": sapon.fs.unit.area,
                "Length": sapon.fs.unit.hot_side.length,
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, sapon):
        stable = sapon.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Volumetric Flowrate": getattr(pyunits.pint_registry, "m**3/second"),
                "Molar Concentration H2O": getattr(pyunits.pint_registry, "mole/m**3"),
                "Molar Concentration NaOH": getattr(pyunits.pint_registry, "mole/m**3"),
                "Molar Concentration EthylAcetate": getattr(
                    pyunits.pint_registry, "mole/m**3"
                ),
                "Molar Concentration SodiumAcetate": getattr(
                    pyunits.pint_registry, "mole/m**3"
                ),
                "Molar Concentration Ethanol": getattr(
                    pyunits.pint_registry, "mole/m**3"
                ),
                "Temperature": getattr(pyunits.pint_registry, "K"),
                "Pressure": getattr(pyunits.pint_registry, "Pa"),
            },
            "Hot Side Inlet": {
                "Volumetric Flowrate": pytest.approx(1e-3, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(55388, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(0, abs=1e-4),
                "Molar Concentration Ethanol": pytest.approx(0, abs=1e-4),
                "Temperature": pytest.approx(320, rel=1e-4),
                "Pressure": pytest.approx(1.0132e05, rel=1e-4),
            },
            "Hot Side Outlet": {
                "Volumetric Flowrate": pytest.approx(1.00, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(100.00, rel=1e-4),
                "Temperature": pytest.approx(298.15, rel=1e-4),
                "Pressure": pytest.approx(1.0132e05, rel=1e-4),
            },
            "Cold Side Inlet": {
                "Volumetric Flowrate": pytest.approx(1e-3, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(55388, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(0, abs=1e-4),
                "Molar Concentration Ethanol": pytest.approx(0, abs=1e-4),
                "Temperature": pytest.approx(300, rel=1e-4),
                "Pressure": pytest.approx(1.0132e05, rel=1e-4),
            },
            "Cold Side Outlet": {
                "Volumetric Flowrate": pytest.approx(1.00, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(100.00, rel=1e-4),
                "Temperature": pytest.approx(298.15, rel=1e-4),
                "Pressure": pytest.approx(1.0132e05, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(1e-3, rel=1e-5) == value(
            sapon.fs.unit.hot_side_outlet.flow_vol[0]
        )
        assert pytest.approx(1e-3, rel=1e-5) == value(
            sapon.fs.unit.cold_side_outlet.flow_vol[0]
        )

        assert 55388.0 == value(sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(
            sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert 0.0 == value(
            sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert 0.0 == value(sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "Ethanol"])

        assert 55388.0 == value(sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(
            sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert 0.0 == value(
            sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert 0.0 == value(sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "Ethanol"])

        assert pytest.approx(318.873, rel=1e-5) == value(
            sapon.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(301.126, rel=1e-5) == value(
            sapon.fs.unit.cold_side_outlet.temperature[0]
        )

        assert pytest.approx(101325, rel=1e-5) == value(
            sapon.fs.unit.hot_side_outlet.pressure[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            sapon.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        hot_side = value(
            sapon.fs.unit.hot_side_outlet.flow_vol[0]
            * sapon.fs.properties.dens_mol
            * sapon.fs.properties.cp_mol
            * (
                sapon.fs.unit.hot_side_inlet.temperature[0]
                - sapon.fs.unit.hot_side_outlet.temperature[0]
            )
        )
        cold_side = value(
            sapon.fs.unit.cold_side_outlet.flow_vol[0]
            * sapon.fs.properties.dens_mol
            * sapon.fs.properties.cp_mol
            * (
                sapon.fs.unit.cold_side_inlet.temperature[0]
                - sapon.fs.unit.cold_side_outlet.temperature[0]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


# # -----------------------------------------------------------------------------
class TestSaponification_countercurrent(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = HX1D(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_type=HeatExchangerFlowPattern.countercurrent,
        )

        m.fs.unit.length.fix(4.85)
        m.fs.unit.area.fix(0.5)
        m.fs.unit.heat_transfer_coefficient.fix(500)

        m.fs.unit.hot_side_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.hot_side_inlet.temperature[0].fix(320)
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.cold_side_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.cold_side_inlet.temperature[0].fix(300)
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        return m

    @pytest.mark.unit
    def test_build(self, sapon):
        assert len(sapon.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.hot_side_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "pressure")

        assert len(sapon.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.hot_side_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "pressure")

        assert len(sapon.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.cold_side_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "pressure")

        assert len(sapon.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.cold_side_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "pressure")

        assert hasattr(sapon.fs.unit, "area")
        assert hasattr(sapon.fs.unit, "length")
        assert hasattr(sapon.fs.unit, "heat_transfer_coefficient")
        assert hasattr(sapon.fs.unit, "heat_transfer_eq")
        assert hasattr(sapon.fs.unit, "heat_conservation")

        assert number_variables(sapon) == 950
        assert number_total_constraints(sapon) == 895
        assert number_unused_variables(sapon) == 16

    @pytest.mark.integration
    def test_units(self, sapon):
        assert_units_equivalent(sapon.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(sapon.fs.unit.length, pyunits.m)
        assert_units_equivalent(
            sapon.fs.unit.heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )

        assert_units_consistent(sapon)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Area": sapon.fs.unit.area,
                "Length": sapon.fs.unit.hot_side.length,
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, sapon):
        stable = sapon.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Volumetric Flowrate": getattr(pyunits.pint_registry, "m**3/second"),
                "Molar Concentration H2O": getattr(pyunits.pint_registry, "mole/m**3"),
                "Molar Concentration NaOH": getattr(pyunits.pint_registry, "mole/m**3"),
                "Molar Concentration EthylAcetate": getattr(
                    pyunits.pint_registry, "mole/m**3"
                ),
                "Molar Concentration SodiumAcetate": getattr(
                    pyunits.pint_registry, "mole/m**3"
                ),
                "Molar Concentration Ethanol": getattr(
                    pyunits.pint_registry, "mole/m**3"
                ),
                "Temperature": getattr(pyunits.pint_registry, "K"),
                "Pressure": getattr(pyunits.pint_registry, "Pa"),
            },
            "Hot Side Inlet": {
                "Volumetric Flowrate": pytest.approx(1e-3, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(55388, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(0, abs=1e-4),
                "Molar Concentration Ethanol": pytest.approx(0, abs=1e-4),
                "Temperature": pytest.approx(320, rel=1e-4),
                "Pressure": pytest.approx(1.0132e05, rel=1e-4),
            },
            "Hot Side Outlet": {
                "Volumetric Flowrate": pytest.approx(1.00, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(100.00, rel=1e-4),
                "Temperature": pytest.approx(298.15, rel=1e-4),
                "Pressure": pytest.approx(1.0132e05, rel=1e-4),
            },
            "Cold Side Inlet": {
                "Volumetric Flowrate": pytest.approx(1e-3, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(55388, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(0, abs=1e-4),
                "Molar Concentration Ethanol": pytest.approx(0, abs=1e-4),
                "Temperature": pytest.approx(300, rel=1e-4),
                "Pressure": pytest.approx(1.0132e05, rel=1e-4),
            },
            "Cold Side Outlet": {
                "Volumetric Flowrate": pytest.approx(1.00, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(100.00, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(100.00, rel=1e-4),
                "Temperature": pytest.approx(298.15, rel=1e-4),
                "Pressure": pytest.approx(1.0132e05, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(1e-3, rel=1e-5) == value(
            sapon.fs.unit.hot_side_outlet.flow_vol[0]
        )
        assert pytest.approx(1e-3, rel=1e-5) == value(
            sapon.fs.unit.cold_side_outlet.flow_vol[0]
        )

        assert 55388.0 == value(sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(
            sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert 0.0 == value(
            sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert 0.0 == value(sapon.fs.unit.hot_side_inlet.conc_mol_comp[0, "Ethanol"])

        assert 55388.0 == value(sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "H2O"])
        assert 100.0 == value(sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "NaOH"])
        assert 100.0 == value(
            sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert 0.0 == value(
            sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert 0.0 == value(sapon.fs.unit.cold_side_inlet.conc_mol_comp[0, "Ethanol"])

        assert pytest.approx(318.869, rel=1e-5) == value(
            sapon.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(301.131, rel=1e-5) == value(
            sapon.fs.unit.cold_side_outlet.temperature[0]
        )

        assert pytest.approx(101325, rel=1e-5) == value(
            sapon.fs.unit.hot_side_outlet.pressure[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            sapon.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        hot_side = value(
            sapon.fs.unit.hot_side_outlet.flow_vol[0]
            * sapon.fs.properties.dens_mol
            * sapon.fs.properties.cp_mol
            * (
                sapon.fs.unit.hot_side_inlet.temperature[0]
                - sapon.fs.unit.hot_side_outlet.temperature[0]
            )
        )
        cold_side = value(
            sapon.fs.unit.cold_side_outlet.flow_vol[0]
            * sapon.fs.properties.dens_mol
            * sapon.fs.properties.cp_mol
            * (
                sapon.fs.unit.cold_side_inlet.temperature[0]
                - sapon.fs.unit.cold_side_outlet.temperature[0]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


# # -----------------------------------------------------------------------------
class TestBT_Generic_cocurrent(object):
    @pytest.fixture(scope="class")
    def btx(self):
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

        m.fs.unit.length.fix(4.85)
        m.fs.unit.area.fix(0.5)
        m.fs.unit.heat_transfer_coefficient.fix(500)

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
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "hot_side_inlet")
        assert len(btx.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_side_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(btx.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_side_inlet")
        assert len(btx.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_side_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(btx.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(btx.fs.unit, "hot_side_outlet")
        assert len(btx.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_side_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(btx.fs.unit.hot_side_outlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_side_outlet")
        assert len(btx.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_side_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(btx.fs.unit.cold_side_outlet, "pressure")

        assert hasattr(btx.fs.unit, "area")
        assert hasattr(btx.fs.unit, "length")
        assert hasattr(btx.fs.unit, "heat_transfer_coefficient")
        assert hasattr(btx.fs.unit, "heat_transfer_eq")
        assert hasattr(btx.fs.unit, "heat_conservation")

        assert number_variables(btx) == 1976
        assert number_total_constraints(btx) == 1867
        assert number_unused_variables(btx) == 36

    @pytest.mark.integration
    def test_units(self, btx):
        assert_units_equivalent(btx.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(btx.fs.unit.length, pyunits.m)
        assert_units_equivalent(
            btx.fs.unit.heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )

        assert_units_consistent(btx)

    @pytest.mark.component
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Area": btx.fs.unit.area,
                "Length": btx.fs.unit.hot_side.length,
            }
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

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
            "Hot Side Inlet": {
                "Total Molar Flowrate": pytest.approx(5.0, rel=1e-4),
                "Total Mole Fraction benzene": pytest.approx(0.5, rel=1e-4),
                "Total Mole Fraction toluene": pytest.approx(0.5, rel=1e-4),
                "Temperature": pytest.approx(365, rel=1e-4),
                "Pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Hot Side Outlet": {
                "Total Molar Flowrate": pytest.approx(100.0, rel=1e-4),
                "Total Mole Fraction benzene": pytest.approx(0.5, rel=1e-4),
                "Total Mole Fraction toluene": pytest.approx(0.5, rel=1e-4),
                "Temperature": pytest.approx(300, rel=1e-4),
                "Pressure": pytest.approx(1e5, rel=1e-4),
            },
            "Cold Side Inlet": {
                "Total Molar Flowrate": pytest.approx(1.0, rel=1e-4),
                "Total Mole Fraction benzene": pytest.approx(0.5, rel=1e-4),
                "Total Mole Fraction toluene": pytest.approx(0.5, rel=1e-4),
                "Temperature": pytest.approx(300, rel=1e-4),
                "Pressure": pytest.approx(101325.0, rel=1e-4),
            },
            "Cold Side Outlet": {
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
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_solution(self, btx):
        # Note hot side in K and cold side in degR
        assert pytest.approx(5, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(356.390, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            btx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(1, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(625.015, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101.325, rel=1e-5) == value(
            btx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.hot_side_inlet.flow_mol[0]
                    - btx.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    btx.fs.unit.cold_side_inlet.flow_mol[0]
                    - btx.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                btx.fs.unit.hot_side.properties[0, 0].enth_mol_phase["Liq"]
                - btx.fs.unit.hot_side.properties[0, 1].enth_mol_phase["Liq"]
            )
        )
        cold_side = value(
            pyunits.convert(
                btx.fs.unit.cold_side_outlet.flow_mol[0]
                * (
                    btx.fs.unit.cold_side.properties[0, 1].enth_mol_phase["Liq"]
                    - btx.fs.unit.cold_side.properties[0, 0].enth_mol_phase["Liq"]
                ),
                to_units=pyunits.W,
            )
        )
        assert abs((hot_side - cold_side) / hot_side) <= 3e-4

    @pytest.mark.component
    def test_initialization_error(self, btx):
        btx.fs.unit.hot_side_outlet.flow_mol[0].fix(20)

        with idaes.temporary_config_ctx():
            with pytest.raises(InitializationError):
                btx.fs.unit.initialize(optarg={"max_iter": 1})
