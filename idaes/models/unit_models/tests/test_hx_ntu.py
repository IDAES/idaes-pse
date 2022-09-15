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
Tests for Plate Heat Exchanger unit model.
Author: Akula Paul, Andrew Lee
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    Param,
    units as pyunits,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.heat_exchanger_ntu import HeatExchangerNTU as HXNTU

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.column_models.properties.MEA_solvent import (
    configuration as aqueous_mea,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.testing import initialization_tester, PhysicalParameterTestBlock
from idaes.core.util.exceptions import ConfigurationError, InitializationError


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_bad_option():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(KeyError):
        m.fs.unit = HXNTU(**{"I'm a bad option": "hot"})


@pytest.mark.unit
def test_bad_option2():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="cold_side_name cannot be 'hot_side'."
    ):
        m.fs.unit = HXNTU(cold_side_name="hot_side")


@pytest.mark.unit
def test_bad_option3():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="hot_side_name cannot be 'cold_side'."
    ):
        m.fs.unit = HXNTU(hot_side_name="cold_side")


@pytest.mark.unit
def test_bad_option4():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="cold_side_name cannot be 'cold_side'."
    ):
        m.fs.unit = HXNTU(cold_side_name="cold_side")


@pytest.mark.unit
def test_bad_option5():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(ConfigurationError, match="hot_side_name cannot be 'hot_side'."):
        m.fs.unit = HXNTU(hot_side_name="hot_side")


@pytest.mark.unit
def test_hot_and_cold_names_same():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        NameError,
        match="HeatExchanger hot and cold side cannot have the same name 'shell'.",
    ):
        m.fs.unit = HXNTU(hot_side_name="shell", cold_side_name="shell")


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
        m.fs.unit = HXNTU(
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
        m.fs.unit = HXNTU(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            cold_side_name="build",
        )


@pytest.mark.unit
def test_user_names():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = HXNTU(
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


# -----------------------------------------------------------------------------
class TestHXNTU(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.hotside_properties = GenericParameterBlock(**aqueous_mea)
        m.fs.coldside_properties = GenericParameterBlock(**aqueous_mea)

        m.fs.unit = HXNTU(
            hot_side={
                "property_package": m.fs.hotside_properties,
                "has_pressure_change": True,
            },
            cold_side={
                "property_package": m.fs.coldside_properties,
                "has_pressure_change": True,
            },
        )

        # Hot fluid
        m.fs.unit.hot_side_inlet.flow_mol[0].fix(60.54879)
        m.fs.unit.hot_side_inlet.temperature[0].fix(392.23)
        m.fs.unit.hot_side_inlet.pressure[0].fix(202650)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "CO2"].fix(0.0158)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "H2O"].fix(0.8747)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "MEA"].fix(0.1095)

        # Cold fluid
        m.fs.unit.cold_side_inlet.flow_mol[0].fix(63.01910)
        m.fs.unit.cold_side_inlet.temperature[0].fix(326.36)
        m.fs.unit.cold_side_inlet.pressure[0].fix(202650)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "CO2"].fix(0.0414)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "H2O"].fix(0.8509)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "MEA"].fix(0.1077)

        # Unit design variables
        m.fs.unit.area.fix(100)
        m.fs.unit.heat_transfer_coefficient.fix(200)
        m.fs.unit.effectiveness.fix(0.7)

        m.fs.unit.hot_side.deltaP.fix(-2000)
        m.fs.unit.cold_side.deltaP.fix(-2000)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):

        assert hasattr(model.fs.unit, "hot_side_inlet")
        assert len(model.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(model.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(model.fs.unit.hot_side_inlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(model.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(model.fs.unit, "hot_side_outlet")
        assert len(model.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(model.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(model.fs.unit.hot_side_outlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(model.fs.unit.hot_side_outlet, "pressure")

        assert hasattr(model.fs.unit, "cold_side_inlet")
        assert len(model.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(model.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(model.fs.unit.cold_side_inlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(model.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(model.fs.unit, "cold_side_outlet")
        assert len(model.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(model.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(model.fs.unit.cold_side_outlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(model.fs.unit.cold_side_outlet, "pressure")

        assert hasattr(model.fs.unit, "heat_duty")
        assert hasattr(model.fs.unit.cold_side, "heat")
        assert hasattr(model.fs.unit.hot_side, "heat")

        assert hasattr(model.fs.unit.cold_side, "deltaP")
        assert hasattr(model.fs.unit.hot_side, "deltaP")

        assert isinstance(model.fs.unit.area, Var)
        assert isinstance(model.fs.unit.heat_transfer_coefficient, Var)
        assert isinstance(model.fs.unit.effectiveness, Var)

        assert isinstance(model.fs.unit.eps_cmin, Param)
        assert value(model.fs.unit.eps_cmin) == 1e-3

        assert isinstance(model.fs.unit.Cmin, Expression)
        assert isinstance(model.fs.unit.Cmax, Expression)
        assert isinstance(model.fs.unit.Cratio, Expression)
        assert isinstance(model.fs.unit.NTU, Expression)

        assert isinstance(model.fs.unit.energy_balance_constraint, Constraint)
        assert isinstance(model.fs.unit.heat_duty_constraint, Constraint)

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

        assert_units_equivalent(model.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(
            model.fs.unit.heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(model.fs.unit.effectiveness[0], pyunits.dimensionless)
        assert_units_equivalent(model.fs.unit.NTU[0], pyunits.dimensionless)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, model):
        perf_dict = model.fs.unit._get_performance_contents()

        assert perf_dict is None

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, model):
        stable = model.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Total Molar Flowrate": getattr(pyunits.pint_registry, "mole/second"),
                "Total Mole Fraction H2O": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "Total Mole Fraction MEA": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "Total Mole Fraction CO2": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "Temperature": getattr(pyunits.pint_registry, "kelvin"),
                "Pressure": getattr(pyunits.pint_registry, "Pa"),
            },
            "Hot Inlet": {
                "Total Molar Flowrate": pytest.approx(60.549, rel=1e-4),
                "Total Mole Fraction H2O": pytest.approx(0.87470, rel=1e-4),
                "Total Mole Fraction MEA": pytest.approx(0.10950, rel=1e-4),
                "Total Mole Fraction CO2": pytest.approx(0.015800, rel=1e-4),
                "Temperature": pytest.approx(392.23, rel=1e-4),
                "Pressure": pytest.approx(2.0265e05, rel=1e-4),
            },
            "Hot Outlet": {
                "Total Molar Flowrate": pytest.approx(1, rel=1e-4),
                "Total Mole Fraction H2O": pytest.approx(1 / 3, rel=1e-4),
                "Total Mole Fraction MEA": pytest.approx(1 / 3, rel=1e-4),
                "Total Mole Fraction CO2": pytest.approx(1 / 3, rel=1e-4),
                "Temperature": pytest.approx(298.15, rel=1e-4),
                "Pressure": pytest.approx(101325, rel=1e-4),
            },
            "Cold Inlet": {
                "Total Molar Flowrate": pytest.approx(63.019, rel=1e-4),
                "Total Mole Fraction H2O": pytest.approx(0.85090, rel=1e-4),
                "Total Mole Fraction MEA": pytest.approx(0.10770, rel=1e-4),
                "Total Mole Fraction CO2": pytest.approx(0.041400, rel=1e-4),
                "Temperature": pytest.approx(326.36, rel=1e-4),
                "Pressure": pytest.approx(2.0265e05, rel=1e-4),
            },
            "Cold Outlet": {
                "Total Molar Flowrate": pytest.approx(1, rel=1e-4),
                "Total Mole Fraction H2O": pytest.approx(1 / 3, rel=1e-4),
                "Total Mole Fraction MEA": pytest.approx(1 / 3, rel=1e-4),
                "Total Mole Fraction CO2": pytest.approx(1 / 3, rel=1e-4),
                "Temperature": pytest.approx(298.15, rel=1e-4),
                "Pressure": pytest.approx(101325, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(200650, rel=1e-5) == value(
            model.fs.unit.hot_side_outlet.pressure[0]
        )
        assert pytest.approx(200650, rel=1e-5) == value(
            model.fs.unit.cold_side_outlet.pressure[0]
        )

        assert pytest.approx(343.995, rel=1e-5) == value(
            model.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(374.333, rel=1e-5) == value(
            model.fs.unit.cold_side_outlet.temperature[0]
        )

        assert pytest.approx(0.015800, rel=1e-5) == value(
            model.fs.unit.hot_side_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.10950, rel=1e-5) == value(
            model.fs.unit.hot_side_outlet.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(0.87470, rel=1e-5) == value(
            model.fs.unit.hot_side_outlet.mole_frac_comp[0, "H2O"]
        )

        assert pytest.approx(0.041400, rel=1e-5) == value(
            model.fs.unit.cold_side_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.10770, rel=1e-5) == value(
            model.fs.unit.cold_side_outlet.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(0.85090, rel=1e-5) == value(
            model.fs.unit.cold_side_outlet.mole_frac_comp[0, "H2O"]
        )

        Cmin = value(
            model.fs.unit.hot_side.properties_in[0].flow_mol
            * model.fs.unit.hot_side.properties_in[0].cp_mol
        )
        assert value(model.fs.unit.Cmin[0]) == pytest.approx(Cmin, rel=1e-5)
        assert value(model.fs.unit.Cmax[0]) == pytest.approx(
            value(
                model.fs.unit.cold_side.properties_in[0].flow_mol
                * model.fs.unit.cold_side.properties_in[0].cp_mol
            ),
            rel=1e-5,
        )
        assert value(model.fs.unit.NTU[0]) == pytest.approx(
            value(
                model.fs.unit.area * model.fs.unit.heat_transfer_coefficient[0] / Cmin
            ),
            rel=1e-5,
        )

        assert pytest.approx(0.7 * Cmin * (392.23 - 326.36), rel=1e-5) == value(
            model.fs.unit.heat_duty[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        # Mass conservation test
        assert (
            abs(
                value(
                    model.fs.unit.hot_side_inlet.flow_mol[0]
                    - model.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        for j in model.fs.hotside_properties.apparent_species_set:
            assert (
                abs(
                    value(
                        model.fs.unit.hot_side_inlet.flow_mol[0]
                        * model.fs.unit.hot_side_inlet.mole_frac_comp[0, j]
                        - model.fs.unit.hot_side_outlet.flow_mol[0]
                        * model.fs.unit.hot_side_outlet.mole_frac_comp[0, j]
                    )
                )
                <= 1e-6
            )

        assert (
            abs(
                value(
                    model.fs.unit.cold_side_inlet.flow_mol[0]
                    - model.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        for j in model.fs.coldside_properties.apparent_species_set:
            assert (
                abs(
                    value(
                        model.fs.unit.cold_side_inlet.flow_mol[0]
                        * model.fs.unit.cold_side_inlet.mole_frac_comp[0, j]
                        - model.fs.unit.cold_side_outlet.flow_mol[0]
                        * model.fs.unit.cold_side_outlet.mole_frac_comp[0, j]
                    )
                )
                <= 1e-6
            )

        # Energy conservation test
        assert (
            abs(
                value(
                    model.fs.unit.hot_side.properties_in[0]._enthalpy_flow_term["Liq"]
                    + model.fs.unit.cold_side.properties_in[0]._enthalpy_flow_term[
                        "Liq"
                    ]
                    - model.fs.unit.hot_side.properties_out[0]._enthalpy_flow_term[
                        "Liq"
                    ]
                    - model.fs.unit.cold_side.properties_out[0]._enthalpy_flow_term[
                        "Liq"
                    ]
                )
            )
            <= 1e-6
        )

    @pytest.mark.component
    def test_initialization_error(self, model):
        model.fs.unit.hot_side_outlet.pressure[0].fix(1)

        with pytest.raises(InitializationError):
            model.fs.unit.initialize()
