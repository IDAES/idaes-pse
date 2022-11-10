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
Tests for costing package based on methods from:

    Process and Product Design Principles: Synthesis, Analysis, and
    Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment
"""
import pytest

from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    units as pyunits,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigValue

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import (
    Compressor,
    CSTR,
    Flash,
    Heater,
    HeatExchanger,
    HeatExchangerNTU,
    PFR,
    PressureChanger,
    Pump,
    StoichiometricReactor,
    Turbine,
)
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
)
from idaes.models.properties import iapws95

from idaes.models.costing.SSLW import (
    SSLWCosting,
    SSLWCostingData,
    HXMaterial,
    HXTubeLength,
    HXType,
    VesselMaterial,
    TrayType,
    TrayMaterial,
    HeaterMaterial,
    HeaterSource,
    PumpType,
    PumpMaterial,
    PumpMotorType,
    CompressorDriveType,
    CompressorMaterial,
    CompressorType,
    FanMaterial,
    FanType,
    BlowerType,
    BlowerMaterial,
)

import logging
from io import StringIO
from pyomo.common.log import LoggingIntercept

# Some more information about this module
__author__ = "Andrew Lee"


solver = get_solver()


@pytest.fixture
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock()

    m.fs.costing = SSLWCosting()

    # Add a placeholder to represent a unit model
    m.fs.unit = UnitModelBlock()

    return m


@pytest.mark.component
@pytest.mark.parametrize("material", HXMaterial)
@pytest.mark.parametrize("hxtype", HXType)
@pytest.mark.parametrize("tube_length", HXTubeLength)
def test_cost_heat_exchanger(model, material, hxtype, tube_length):
    model.fs.unit.area = Param(initialize=1000, units=pyunits.m**2)
    model.fs.unit.hot_side = Block()
    model.fs.unit.hot_side.properties_in = Block(model.fs.time)
    model.fs.unit.hot_side.properties_in[0].pressure = Param(
        initialize=2, units=pyunits.atm
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_heat_exchanger,
        costing_method_arguments={
            "hx_type": hxtype,
            "material_type": material,
            "tube_length": tube_length,
        },
    )

    assert isinstance(model.fs.unit.costing.base_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.number_of_units, Var)
    assert isinstance(model.fs.unit.costing.pressure_factor, Var)
    assert isinstance(model.fs.unit.costing.material_factor, Var)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.hx_material_eqn, Constraint)
    assert isinstance(model.fs.unit.costing.p_factor_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)

    # Test solution for one known case
    if (
        material == HXMaterial.StainlessSteelStainlessSteel
        and hxtype == HXType.Utube
        and tube_length == HXTubeLength.TwelveFoot
    ):
        assert pytest.approx(87704.6, 1e-5) == value(
            model.fs.unit.costing.base_cost_per_unit
        )
        assert pytest.approx(0.982982, 1e-5) == value(
            model.fs.unit.costing.pressure_factor
        )
        assert pytest.approx(4.08752, 1e-5) == value(
            model.fs.unit.costing.material_factor
        )

        assert pytest.approx(476063, 1e-5) == value(
            pyunits.convert(
                model.fs.unit.costing.capital_cost, to_units=pyunits.USD_2018
            )
        )


@pytest.mark.component
@pytest.mark.parametrize("material_type", VesselMaterial)
@pytest.mark.parametrize("weight_limit", [1, 2])
@pytest.mark.parametrize("aspect_ratio_range", [1, 2])
@pytest.mark.parametrize("include_pl", [True, False])
def test_cost_vessel(
    model, material_type, weight_limit, aspect_ratio_range, include_pl
):
    model.fs.unit.length = Param(initialize=0.00075, units=pyunits.m)
    model.fs.unit.diameter = Param(initialize=2, units=pyunits.m)

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_vessel,
        costing_method_arguments={
            "vertical": True,
            "material_type": material_type,
            "weight_limit": weight_limit,
            "aspect_ratio_range": aspect_ratio_range,
            "include_platforms_ladders": include_pl,
        },
    )

    assert isinstance(model.fs.unit.costing.shell_thickness, Param)
    assert isinstance(model.fs.unit.costing.material_factor, Param)
    assert isinstance(model.fs.unit.costing.material_density, Param)

    assert isinstance(model.fs.unit.costing.base_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.weight, Var)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.weight_eq, Constraint)

    # Platforms and ladders
    if include_pl:
        assert isinstance(model.fs.unit.costing.base_cost_platforms_ladders, Var)

        assert isinstance(model.fs.unit.costing.cost_platforms_ladders_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)

    # Test solution for one known case
    if (
        material_type == VesselMaterial.CarbonSteel
        and weight_limit == 1
        and aspect_ratio_range == 1
        and include_pl
    ):
        assert pytest.approx(35958, 1e-5) == value(
            pyunits.convert(
                model.fs.unit.costing.capital_cost, to_units=pyunits.USD_2018
            )
        )


@pytest.mark.component
@pytest.mark.parametrize("tray_material", TrayMaterial)
@pytest.mark.parametrize("tray_type", TrayType)
def test_cost_vessel_trays(model, tray_material, tray_type):
    model.fs.unit.length = Param(initialize=0.00075, units=pyunits.m)
    model.fs.unit.diameter = Param(initialize=2, units=pyunits.m)

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_vessel,
        costing_method_arguments={
            "vertical": True,
            "number_of_trays": 10,
            "tray_material": tray_material,
            "tray_type": tray_type,
        },
    )

    assert isinstance(model.fs.unit.costing.shell_thickness, Param)
    assert isinstance(model.fs.unit.costing.material_factor, Param)
    assert isinstance(model.fs.unit.costing.material_density, Param)

    assert isinstance(model.fs.unit.costing.base_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.weight, Var)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.weight_eq, Constraint)

    assert isinstance(model.fs.unit.costing.base_cost_platforms_ladders, Var)

    assert isinstance(model.fs.unit.costing.cost_platforms_ladders_eq, Constraint)

    assert isinstance(model.fs.unit.costing.tray_type_factor, Param)

    assert isinstance(model.fs.unit.costing.base_cost_trays, Var)
    assert isinstance(model.fs.unit.costing.tray_material_factor, Var)
    assert isinstance(model.fs.unit.costing.number_trays_factor, Var)
    assert isinstance(model.fs.unit.costing.base_cost_per_tray, Var)

    assert isinstance(model.fs.unit.costing.tray_material_factor_eq, Constraint)
    assert isinstance(model.fs.unit.costing.num_tray_factor_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.single_tray_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.tray_costing_constraint, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)


@pytest.mark.component
@pytest.mark.parametrize("material_type", VesselMaterial)
def test_cost_vessel_horizontal(model, material_type):
    model.fs.unit.length = Param(initialize=0.00075, units=pyunits.m)
    model.fs.unit.diameter = Param(initialize=2, units=pyunits.m)

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_vessel,
        costing_method_arguments={
            "vertical": False,
            "material_type": material_type,
            "weight_limit": 1,
            "include_platforms_ladders": False,
        },
    )

    assert isinstance(model.fs.unit.costing.shell_thickness, Param)
    assert isinstance(model.fs.unit.costing.material_factor, Param)
    assert isinstance(model.fs.unit.costing.material_density, Param)

    assert isinstance(model.fs.unit.costing.base_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.weight, Var)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.weight_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)


@pytest.mark.component
@pytest.mark.parametrize("material_type", HeaterMaterial)
@pytest.mark.parametrize("heat_source", HeaterSource)
def test_cost_fired_heater(model, material_type, heat_source):
    model.fs.unit.heat_duty = Param([0], initialize=1000, units=pyunits.kJ / pyunits.s)
    model.fs.unit.control_volume = Block()
    model.fs.unit.control_volume.properties_in = Block(model.fs.time)
    model.fs.unit.control_volume.properties_in[0].pressure = Param(
        initialize=2, units=pyunits.atm
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_fired_heater,
        costing_method_arguments={
            "material_type": material_type,
            "heat_source": heat_source,
        },
    )

    assert isinstance(model.fs.unit.costing.pressure_factor, Var)
    assert isinstance(model.fs.unit.costing.base_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.capital_cost, Var)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_cost_per_unit_eq, Constraint)
    assert isinstance(model.fs.unit.costing.pressure_factor_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)


@pytest.mark.component
def test_cost_turbine(model):
    model.fs.unit.config.declare("compressor", ConfigValue(default=False))
    model.fs.unit.work_mechanical = Param(
        [0], initialize=-1224638, units=pyunits.J / pyunits.s
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_turbine,
    )

    assert isinstance(model.fs.unit.costing.capital_cost, Var)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)

    assert pytest.approx(213199, 1e-5) == value(model.fs.unit.costing.capital_cost)


@pytest.mark.component
@pytest.mark.parametrize(
    "material_type",
    [
        PumpMaterial.CastIron,
        PumpMaterial.DuctileIron,
        PumpMaterial.CastSteel,
        PumpMaterial.Bronze,
        PumpMaterial.StainlessSteel,
        PumpMaterial.HastelloyC,
        PumpMaterial.Monel,
        PumpMaterial.Nickel,
        PumpMaterial.Titanium,
    ],
)
@pytest.mark.parametrize("pump_type_factor", [1.1, 1.2, 1.3, 1.4, 2.1, 2.2])
@pytest.mark.parametrize("motor_type", PumpMotorType)
def test_cost_pump_centrifugal(model, material_type, pump_type_factor, motor_type):
    model.fs.unit.config.declare("compressor", ConfigValue(default=True))
    model.fs.unit.work_mechanical = Param(
        [0], initialize=10143.79, units=pyunits.J / pyunits.s
    )
    model.fs.unit.deltaP = Param([0], initialize=50000, units=pyunits.Pa)
    model.fs.unit.control_volume = Block()
    model.fs.unit.control_volume.properties_in = Block(model.fs.time)
    model.fs.unit.control_volume.properties_in[0].dens_mass = Param(
        initialize=986.64, units=pyunits.kg / pyunits.m**3
    )
    model.fs.unit.control_volume.properties_in[0].flow_vol = Param(
        initialize=0.182592, units=pyunits.m**3 / pyunits.s
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_pump,
        costing_method_arguments={
            "pump_type": PumpType.Centrifugal,
            "material_type": material_type,
            "pump_type_factor": pump_type_factor,
            "motor_type": motor_type,
        },
    )

    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.number_of_units, Var)
    assert isinstance(model.fs.unit.costing.pump_head, Var)
    assert isinstance(model.fs.unit.costing.size_factor, Var)
    assert isinstance(model.fs.unit.costing.base_pump_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.pump_capital_cost, Var)
    assert isinstance(model.fs.unit.costing.base_motor_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.motor_capital_cost, Var)

    assert isinstance(model.fs.unit.costing.material_factor, Param)
    assert isinstance(model.fs.unit.costing.FT, Param)
    assert isinstance(model.fs.unit.costing.motor_FT, Param)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_pump_cost_per_unit_eq, Constraint)
    assert isinstance(model.fs.unit.costing.base_motor_cost_eq, Constraint)
    assert isinstance(model.fs.unit.costing.pump_capital_cost_eq, Constraint)
    assert isinstance(model.fs.unit.costing.motor_capital_cost_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)

    # Test solution for one known case
    if (
        material_type == PumpMaterial.Nickel
        and pump_type_factor == 1.4
        and motor_type == PumpMotorType.Enclosed
    ):
        assert pytest.approx(62868.3, 1e-5) == value(
            pyunits.convert(
                model.fs.unit.costing.capital_cost, to_units=pyunits.USD_2018
            )
        )


@pytest.mark.component
@pytest.mark.parametrize(
    "material_type",
    [
        PumpMaterial.CastIron,
        PumpMaterial.DuctileIron,
        PumpMaterial.CastSteel,
        PumpMaterial.Bronze,
        PumpMaterial.StainlessSteel,
        PumpMaterial.HastelloyC,
        PumpMaterial.Monel,
        PumpMaterial.Nickel,
        PumpMaterial.Titanium,
    ],
)
@pytest.mark.parametrize("motor_type", PumpMotorType)
def test_cost_pump_ExternalGear(model, material_type, motor_type):
    model.fs.unit.config.declare("compressor", ConfigValue(default=True))
    model.fs.unit.work_mechanical = Param(
        [0], initialize=10143.79, units=pyunits.J / pyunits.s
    )
    model.fs.unit.deltaP = Param([0], initialize=50000, units=pyunits.Pa)
    model.fs.unit.control_volume = Block()
    model.fs.unit.control_volume.properties_in = Block(model.fs.time)
    model.fs.unit.control_volume.properties_in[0].dens_mass = Param(
        initialize=986.64, units=pyunits.kg / pyunits.m**3
    )
    model.fs.unit.control_volume.properties_in[0].flow_vol = Param(
        initialize=0.182592, units=pyunits.m**3 / pyunits.s
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_pump,
        costing_method_arguments={
            "pump_type": PumpType.ExternalGear,
            "material_type": material_type,
            "motor_type": motor_type,
        },
    )

    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.number_of_units, Var)
    assert isinstance(model.fs.unit.costing.pump_head, Var)
    assert isinstance(model.fs.unit.costing.size_factor, Var)
    assert isinstance(model.fs.unit.costing.base_pump_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.pump_capital_cost, Var)
    assert isinstance(model.fs.unit.costing.base_motor_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.motor_capital_cost, Var)

    assert isinstance(model.fs.unit.costing.material_factor, Param)
    assert isinstance(model.fs.unit.costing.FT, Param)
    assert isinstance(model.fs.unit.costing.motor_FT, Param)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_pump_cost_per_unit_eq, Constraint)
    assert isinstance(model.fs.unit.costing.base_motor_cost_eq, Constraint)
    assert isinstance(model.fs.unit.costing.pump_capital_cost_eq, Constraint)
    assert isinstance(model.fs.unit.costing.motor_capital_cost_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)


@pytest.mark.component
@pytest.mark.parametrize(
    "material_type",
    [
        PumpMaterial.DuctileIron,
        PumpMaterial.StainlessSteel,
        PumpMaterial.NiAlBronze,
        PumpMaterial.CarbonSteel,
    ],
)
@pytest.mark.parametrize("motor_type", PumpMotorType)
def test_cost_pump_reciprocating(model, material_type, motor_type):
    model.fs.unit.config.declare("compressor", ConfigValue(default=True))
    model.fs.unit.work_mechanical = Param(
        [0], initialize=10143.79, units=pyunits.J / pyunits.s
    )
    model.fs.unit.deltaP = Param([0], initialize=50000, units=pyunits.Pa)
    model.fs.unit.control_volume = Block()
    model.fs.unit.control_volume.properties_in = Block(model.fs.time)
    model.fs.unit.control_volume.properties_in[0].dens_mass = Param(
        initialize=986.64, units=pyunits.kg / pyunits.m**3
    )
    model.fs.unit.control_volume.properties_in[0].flow_vol = Param(
        initialize=0.182592, units=pyunits.m**3 / pyunits.s
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_pump,
        costing_method_arguments={
            "pump_type": PumpType.Reciprocating,
            "material_type": material_type,
            "motor_type": motor_type,
        },
    )

    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.number_of_units, Var)
    assert isinstance(model.fs.unit.costing.pump_head, Var)
    assert isinstance(model.fs.unit.costing.size_factor, Var)
    assert isinstance(model.fs.unit.costing.base_pump_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.pump_capital_cost, Var)
    assert isinstance(model.fs.unit.costing.base_motor_cost_per_unit, Var)
    assert isinstance(model.fs.unit.costing.motor_capital_cost, Var)

    assert isinstance(model.fs.unit.costing.material_factor, Param)
    assert isinstance(model.fs.unit.costing.FT, Param)
    assert isinstance(model.fs.unit.costing.motor_FT, Param)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_pump_cost_per_unit_eq, Constraint)
    assert isinstance(model.fs.unit.costing.base_motor_cost_eq, Constraint)
    assert isinstance(model.fs.unit.costing.pump_capital_cost_eq, Constraint)
    assert isinstance(model.fs.unit.costing.motor_capital_cost_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)


@pytest.mark.component
@pytest.mark.parametrize("compressor_type", CompressorType)
@pytest.mark.parametrize("drive_type", CompressorDriveType)
@pytest.mark.parametrize("material_type", CompressorMaterial)
def test_cost_compressor(model, compressor_type, drive_type, material_type):
    model.fs.unit.config.declare("compressor", ConfigValue(default=True))
    model.fs.unit.config.declare(
        "thermodynamic_assumption",
        ConfigValue(default=ThermodynamicAssumption.isothermal),
    )
    model.fs.unit.config.thermodynamic_assumption = ThermodynamicAssumption.isentropic
    model.fs.unit.work_mechanical = Param(
        [0], initialize=101410.4, units=pyunits.J / pyunits.s
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_compressor,
        costing_method_arguments={
            "compressor_type": compressor_type,
            "drive_type": drive_type,
            "material_type": material_type,
        },
    )

    assert isinstance(model.fs.unit.costing.capital_cost, Var)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)

    if (
        compressor_type == CompressorType.Centrifugal
        and drive_type == CompressorDriveType.ElectricMotor
        and material_type == CompressorMaterial.StainlessSteel
    ):
        assert pytest.approx(300695, 1e-5) == value(
            pyunits.convert(
                model.fs.unit.costing.capital_cost, to_units=pyunits.USD_2018
            )
        )


@pytest.mark.component
@pytest.mark.parametrize("compressor_type", CompressorType)
@pytest.mark.parametrize("drive_type", CompressorDriveType)
@pytest.mark.parametrize("material_type", CompressorMaterial)
def test_not_compressor(model, compressor_type, drive_type, material_type):
    # Test exception for non-supported compressor flags
    model.fs.unit.config.declare("compressor", ConfigValue(default=False))
    model.fs.unit.config.declare(
        "thermodynamic_assumption",
        ConfigValue(default=ThermodynamicAssumption.isentropic),
    )
    model.fs.unit.work_mechanical = Param(
        [0], initialize=101410.4, units=pyunits.J / pyunits.s
    )

    expected_string = (
        "cost_compressor method is only appropriate for "
        "pressure changers with the compressor argument "
        "equal to True."
    )

    stream = StringIO()
    with LoggingIntercept(stream, "idaes", logging.WARNING):
        model.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method=SSLWCostingData.cost_compressor,
            costing_method_arguments={
                "compressor_type": compressor_type,
                "drive_type": drive_type,
                "material_type": material_type,
            },
        )

    assert expected_string in str(stream.getvalue())


@pytest.mark.component
@pytest.mark.parametrize("fan_type", FanType)
@pytest.mark.parametrize("material_type", FanMaterial)
def test_cost_fan(model, fan_type, material_type):
    model.fs.unit.config.declare("compressor", ConfigValue(default=True))
    model.fs.unit.control_volume = Block()
    model.fs.unit.control_volume.properties_in = Block(model.fs.time)
    model.fs.unit.control_volume.properties_in[0].flow_vol = Param(
        initialize=0.182592, units=pyunits.m**3 / pyunits.s
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_fan,
        costing_method_arguments={"fan_type": fan_type, "material_type": material_type},
    )

    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.number_of_units, Var)
    assert isinstance(model.fs.unit.costing.base_cost_per_unit, Var)

    assert isinstance(model.fs.unit.costing.head_factor, Param)
    assert isinstance(model.fs.unit.costing.material_factor, Param)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_cost_per_unit_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)

    # TODO: Test case for solution checking


@pytest.mark.component
@pytest.mark.parametrize("blower_type", BlowerType)
@pytest.mark.parametrize("material_type", BlowerMaterial)
def test_cost_blower(model, blower_type, material_type):
    model.fs.unit.config.declare("compressor", ConfigValue(default=True))
    model.fs.unit.work_mechanical = Param(
        [0], initialize=101410.4, units=pyunits.J / pyunits.s
    )

    model.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=model.fs.costing,
        costing_method=SSLWCostingData.cost_blower,
        costing_method_arguments={
            "blower_type": blower_type,
            "material_type": material_type,
        },
    )

    assert isinstance(model.fs.unit.costing.capital_cost, Var)
    assert isinstance(model.fs.unit.costing.number_of_units, Var)
    assert isinstance(model.fs.unit.costing.base_cost_per_unit, Var)

    assert isinstance(model.fs.unit.costing.material_factor, Param)

    assert isinstance(model.fs.unit.costing.capital_cost_constraint, Constraint)
    assert isinstance(model.fs.unit.costing.base_cost_per_unit_eq, Constraint)

    assert degrees_of_freedom(model) == 0
    assert_units_consistent(model.fs.unit.costing)

    res = solver.solve(model)

    assert check_optimal_termination(res)

    # TODO: Test case for solution checking


@pytest.mark.integration
class TestMapping:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock()

        m.fs.pparams = PhysicalParameterTestBlock()
        m.fs.rparams = ReactionParameterTestBlock(property_package=m.fs.pparams)

        m.fs.costing = SSLWCosting()

        return m

    def test_compressor(self, model):
        # Add examples of supported unit models and add costing
        model.fs.C101 = Compressor(property_package=model.fs.pparams)
        model.fs.C101.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert model.fs.C101.costing.drive_factor.value == 1
        assert model.fs.C101.costing.material_factor.value == 2.5

    def test_cstr(self, model):
        model.fs.R102 = CSTR(
            property_package=model.fs.pparams, reaction_package=model.fs.rparams
        )
        # Add length and diameter to CSTR
        model.fs.R102.length = 1 * pyunits.m
        model.fs.R102.diameter = 1 * pyunits.m
        model.fs.R102.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert model.fs.R102.costing.material_factor.value == 1
        assert model.fs.R102.costing.material_density.value == 0.284

    def test_flash(self, model):
        model.fs.F103 = Flash(property_package=model.fs.pparams)
        # Add length and diameter to Flash
        model.fs.F103.length = 1 * pyunits.m
        model.fs.F103.diameter = 1 * pyunits.m
        model.fs.F103.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert model.fs.F103.costing.material_factor.value == 1
        assert model.fs.F103.costing.material_density.value == 0.284

    def test_heater(self, model):
        # Add examples of supported unit models and add costing
        model.fs.unit = Heater(property_package=model.fs.pparams)
        model.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert model.fs.unit.costing.pressure_factor.value == 1.1
        assert model.fs.unit.costing.material_factor.value == 1

    def test_hx0D(self, model):
        # Add examples of supported unit models and add costing
        model.fs.unit = HeatExchanger(
            hot_side={"property_package": model.fs.pparams},
            cold_side={"property_package": model.fs.pparams},
        )
        model.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert model.fs.unit.costing.length_factor.value == 1.12

    # TODO : Test for HX1D once supported

    def test_hx_ntu(self):
        # Need a different property package here
        m = ConcreteModel()

        m.fs = FlowsheetBlock()

        m.fs.pparams = iapws95.Iapws95ParameterBlock()

        m.fs.costing = SSLWCosting()
        # Add examples of supported unit models and add costing
        m.fs.unit = HeatExchangerNTU(
            hot_side={"property_package": m.fs.pparams},
            cold_side={"property_package": m.fs.pparams},
        )
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        assert m.fs.unit.costing.length_factor.value == 1.12

    def test_pfr(self, model):
        model.fs.unit = PFR(
            property_package=model.fs.pparams, reaction_package=model.fs.rparams
        )
        # Add length and diameter to PFR
        model.fs.unit.length = 1 * pyunits.m
        model.fs.unit.diameter = 1 * pyunits.m
        model.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert model.fs.unit.costing.material_factor.value == 1
        assert model.fs.unit.costing.material_density.value == 0.284

    def test_pressure_changer(self, model):
        # Add examples of supported unit models and add costing
        model.fs.unit = PressureChanger(
            property_package=model.fs.pparams,
            thermodynamic_assumption=ThermodynamicAssumption.isentropic,
        )
        model.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert model.fs.unit.costing.drive_factor.value == 1
        assert model.fs.unit.costing.material_factor.value == 2.5

    def test_pump_compressor(self, model):
        # Test exception for non-supported compressor flags
        model.fs.unit = Compressor(
            property_package=model.fs.pparams,
            thermodynamic_assumption=ThermodynamicAssumption.pump,
        )

        expected_string = (
            "fs.unit - pressure changers with the pump "
            "assumption should use the cost_pump method."
        )

        stream = StringIO()
        with LoggingIntercept(stream, "idaes", logging.WARNING):
            model.fs.unit.costing = UnitModelCostingBlock(
                flowsheet_costing_block=model.fs.costing
            )

        assert expected_string in str(stream.getvalue())

    def test_isothermal_compressor(self, model):
        # Test exception for non-supported compressor flags
        model.fs.unit = Compressor(
            property_package=model.fs.pparams,
            thermodynamic_assumption=ThermodynamicAssumption.isothermal,
        )

        expected_string = (
            "fs.unit - pressure changers without isentropic "
            "assumption are too simple to be costed."
        )

        stream = StringIO()
        with LoggingIntercept(stream, "idaes", logging.WARNING):
            model.fs.unit.costing = UnitModelCostingBlock(
                flowsheet_costing_block=model.fs.costing
            )

        assert expected_string in str(stream.getvalue())

    def test_adiabatic_compressor(self, model):
        # Test exception for non-supported compressor flags
        model.fs.unit = Compressor(
            property_package=model.fs.pparams,
            thermodynamic_assumption=ThermodynamicAssumption.adiabatic,
        )

        expected_string = (
            "fs.unit - pressure changers without isentropic "
            "assumption are too simple to be costed."
        )

        stream = StringIO()
        with LoggingIntercept(stream, "idaes", logging.WARNING):
            model.fs.unit.costing = UnitModelCostingBlock(
                flowsheet_costing_block=model.fs.costing
            )

        assert expected_string in str(stream.getvalue())

    def test_pump(self):
        # Need a different property package here
        m = ConcreteModel()

        m.fs = FlowsheetBlock()

        m.fs.pparams = iapws95.Iapws95ParameterBlock()

        m.fs.costing = SSLWCosting()

        # Add examples of supported unit models and add costing
        m.fs.unit = Pump(property_package=m.fs.pparams)
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        assert hasattr(m.fs.unit.costing, "pump_head")

    def test_rstoich(self, model):
        model.fs.unit = StoichiometricReactor(
            property_package=model.fs.pparams, reaction_package=model.fs.rparams
        )
        # Add length and diameter to reactor
        model.fs.unit.length = 1 * pyunits.m
        model.fs.unit.diameter = 1 * pyunits.m
        model.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert model.fs.unit.costing.material_factor.value == 1
        assert model.fs.unit.costing.material_density.value == 0.284

    def test_turbine(self, model):
        # Add examples of supported unit models and add costing
        model.fs.unit = Turbine(property_package=model.fs.pparams)
        model.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )
        assert hasattr(model.fs.unit.costing, "capital_cost")
        assert not hasattr(model.fs.unit.costing, "material_factor")
