#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

# TODO:
# 1. Build out graph (go run Jupyter at TSA0D model last step)

"""
UI exports for 0D Fixed Bed TSA unit model.
"""
from pyomo.environ import ConcreteModel, SolverFactory, units
from pyomo.common.dependencies import attempt_import

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util as iutil
import idaes.logger as idaes_log

from idaes.models_extra.temperature_swing_adsorption import (
    FixedBedTSA0D,
    FixedBedTSA0DInitializer,
    Adsorbent,
    SteamCalculationType,
)

from idaes.models_extra.temperature_swing_adsorption.util import tsa_summary

# pylint: disable=possibly-used-before-assignment
_, flowsheet_processor_available = attempt_import("idaes_flowsheet_processor")
if flowsheet_processor_available:
    from idaes_flowsheet_processor import api

_log = idaes_log.getLogger(__name__)

model_name = "0D Fixed Bed TSA"
model_name_for_ui = "0D Fixed Bed TSA"
unit_name = "FixedBedTSA0D"


def export_to_ui() -> "api.FlowsheetInterface":
    """
    Export the model to UI, define the UI options and build the model.
    Return:
        fsi object.
    """
    fsi = api.FlowsheetInterface(
        name="0D Fixed Bed TSA",
        description=model_name_for_ui,  # The water tap UI use description as flowsheet name
        do_export=export,
        do_build=build,
        do_solve=solve,
        build_options={
            "adsorbent": {
                "name": "Adsorbent",
                "display_name": "Adsorbent material",
                "values_allowed": [
                    "zeolite_13x",
                    "mmen_mg_mof_74",
                    "polystyrene_amine",
                ],
                "value": "zeolite_13x",
                "category": unit_name,
            },
            "number_of_beds": {
                "name": "number_of_beds",
                "display_name": "Number of beds",
                "description": "Number of beds in fixed bed TSA system used to split the mole flow rate at feed",
                "min_val": 1,
                "max_val": 100000,
                "value": 1,
                "values_allowed": "int",
                "category": unit_name,
            },
            "transformation_method": {
                "name": "transformation_method",
                "display_name": "Transformation Method",
                "description": "Method to use for DAE transformation",
                "values_allowed": ["dae.finite_difference", "dae.collocation"],
                "value": "dae.collocation",
                "category": unit_name,
            },
            "transformation_scheme": {
                "name": "transformation_scheme",
                "display_name": "Transformation Scheme",
                "values_allowed": [
                    "useDefault",
                    "backward",
                    "forward",
                    "lagrangeRadau",
                ],
                "value": "lagrangeRadau",
                "category": unit_name,
            },
            # TODO
            # This part need to be conditional render based on transformation_method in UI
            # This parts need to conditional render based on transformation_method
            #  this is only show when dae.finite_difference is selected
            # "finite_elements": {
            #     "name": "finite_elements",
            #     "display_name": "Number of Finite Elements",
            #     "description": "Number of finite elements to use when discretizing time domain",
            #     "min_val": 0,
            #     "max_val": 10000,
            #     "value": 20,
            #     "values_allowed": "int",
            #     "category": unit_name,
            # },
            # TODO: this only show when dae.collocation is selected
            "collocation_points": {
                "name": "collocation_points",
                "display_name": "Collocation Points",
                "description": "Number of collocation points per finite element",
                "min_val": 0,
                "max_val": 10000,
                "value": 6,
                "values_allowed": "int",
                "category": unit_name,
            },
        },
    )
    return fsi


def build(build_options=None, **kwargs):
    """
    Build the model

    Return:
        the model object.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    unit_opt = {k: v for k, v in build_options.items() if v.category == unit_name}

    # build args for m.fs.tsa
    tsa_params = {
        "adsorbent": Adsorbent[unit_opt["adsorbent"].value],  # Use zeolite_13x
        "number_of_beds": unit_opt["number_of_beds"].value,  # Set to 1
        "transformation_method": unit_opt[
            "transformation_method"
        ].value,  # Use dae.collocation
        "transformation_scheme": unit_opt[
            "transformation_scheme"
        ].value,  # default Use lagrangeRadau
        "steam_calculation": SteamCalculationType.none,
    }

    # Conditional parameters collocation_points base on transformation_method is dae.collocation
    if (
        "collocation_points" in unit_opt
        and unit_opt["collocation_points"].value is not None
    ):
        tsa_params["collocation_points"] = unit_opt["collocation_points"].value

    # Conditional parameters finite_elements base on transformation_method is dae.finite_difference
    if "finite_elements" in unit_opt and unit_opt["finite_elements"].value is not None:
        tsa_params["finite_elements"] = unit_opt["finite_elements"].value

    m.fs.tsa = FixedBedTSA0D(**tsa_params)

    flue_gas = {
        "flow_mol_comp": {
            "H2O": 0.0,
            "CO2": 0.00960 * 0.12,
            "N2": 0.00960 * 0.88,
            "O2": 0.0,
        },
        "temperature": 300.0,
        "pressure": 1.0e5,
    }

    for i in m.fs.tsa.component_list:
        m.fs.tsa.inlet.flow_mol_comp[:, i].fix(flue_gas["flow_mol_comp"][i])
    m.fs.tsa.inlet.temperature.fix(flue_gas["temperature"])
    m.fs.tsa.inlet.pressure.fix(flue_gas["pressure"])

    m.fs.tsa.temperature_desorption.fix(430)
    m.fs.tsa.temperature_adsorption.fix(310)
    m.fs.tsa.temperature_heating.fix(440)
    m.fs.tsa.temperature_cooling.fix(300)
    m.fs.tsa.bed_diameter.fix(3 / 100)
    m.fs.tsa.bed_height.fix(1.2)

    DOF = degrees_of_freedom(m)
    _log.info(f"The DOF of the TSA unit is {DOF}")

    return m


def export(flowsheet=None, exports=None, build_options=None, **kwargs):
    """
    Define which variables to export to UI exports.

    Return:
        None.
    """
    fs = flowsheet

    for compound in "H2O", "CO2", "N2", "O2":
        exports.add(
            obj=fs.tsa.inlet.flow_mol_comp[0, compound],
            name=f"{compound} molar flow rate",
            ui_units=units.mole / units.s,
            display_units="mol/s",
            description=f"Inlet molar flow rate for {compound}",
            rounding=6,
            input_category="Feed",
            output_category="Feed",
            is_input=True,
            is_output=True,
        )

    # Feed
    exports.add(
        obj=fs.tsa.inlet.pressure[0],
        name="Pressure",
        ui_units=units.Pa,
        display_units="Pa",
        rounding=1,
        input_category="Feed",
        output_category="Feed",
        is_input=True,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.inlet.temperature[0],
        name="Temperature",
        ui_units=units.K,
        display_units="K",
        rounding=6,
        input_category="Feed",
        output_category="Feed",
        is_input=True,
        is_output=True,
    )

    # equipment parameters
    exports.add(
        obj=fs.tsa.number_beds,  # Number of beds
        name="Number of beds",
        ui_units=units.dimensionless,
        display_units="-",
        rounding=4,
        output_category="Equipment Parameters",
        is_input=False,
        is_output=True,
    )
    exports.add(
        obj=fs.tsa.bed_diameter,
        name="Bed diameter",  # column diameter from original model
        ui_units=units.m,
        display_units="m",
        rounding=6,
        input_category="Equipment Parameters",
        output_category="Equipment Parameters",
        is_input=True,
        is_output=True,
    )
    exports.add(
        obj=fs.tsa.bed_height,
        name="Bed length",  # column length from original model
        ui_units=units.m,
        display_units="m",
        rounding=6,
        input_category="Equipment Parameters",
        output_category="Equipment Parameters",
        is_input=True,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.temperature_adsorption,
        name="Adsorption temperature",
        ui_units=units.K,
        display_units="K",
        rounding=6,
        input_category="Operating Parameters",
        output_category="Operating Parameters",
        is_input=True,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.temperature_desorption,
        name="Desorption temperature",
        ui_units=units.K,
        display_units="K",
        rounding=6,
        input_category="Operating Parameters",
        output_category="Operating Parameters",
        is_input=True,
        is_output=True,
    )
    exports.add(
        obj=fs.tsa.temperature_heating,
        name="Heating temperature",
        ui_units=units.K,
        display_units="K",
        rounding=6,
        input_category="Operating Parameters",
        output_category="Operating Parameters",
        is_input=True,
        is_output=True,
    )
    exports.add(
        obj=fs.tsa.temperature_cooling,
        name="Cooling temperature",
        ui_units=units.K,
        display_units="K",
        rounding=6,
        input_category="Operating Parameters",
        output_category="Operating Parameters",
        is_input=True,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.bed_volume,
        name="Bed volume",  # column volume from original model
        ui_units=units.m**3,
        display_units="m3",
        rounding=8,
        output_category="Operating Parameters",
        is_input=False,
        is_output=True,
    )
    exports.add(
        obj=fs.tsa.mole_frac_in["CO2"] * 100,
        name="CO2 mole fraction at feed",
        ui_units=units.dimensionless,
        display_units="%",
        rounding=6,
        output_category="Operating Parameters",
        is_input=False,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.flow_mol_in_total,
        name="Feed flow rate",
        ui_units=units.mol / units.s,
        display_units="mol/s",
        rounding=6,
        input_category="Operating Parameters",
        output_category="Operating Parameters",
        is_input=False,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.velocity_in,
        name="Feed velocity",
        ui_units=units.m / units.s,
        display_units="m/s",
        rounding=5,
        output_category="Operating Parameters",
        is_input=False,
        is_output=True,
    )

    # Minimum fluidization velocity
    exports.add(
        obj=fs.tsa.velocity_mf,
        name="Minimum fluidization velocity",
        ui_units=units.m / units.s,
        display_units="m/s",
        rounding=4,
        input_category="Operating Parameters",
        output_category="Operating Parameters",
        is_input=False,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.cycle_time,  # Cycle time
        name="Cycle time",
        ui_units=units.h,
        display_units="h",
        rounding=6,
        input_category="Operating Parameters",
        output_category="Operating Parameters",
        is_input=False,
        is_output=True,
    )

    # Pressure and Physical Parameters
    exports.add(
        obj=fs.tsa.pressure_drop,  # Pressure drop
        name="Pressure drop",
        ui_units=units.Pa,
        display_units="Pa",
        rounding=1,
        input_category="Operating Parameters",
        output_category="Operating Parameters",
        is_input=False,
        is_output=True,
    )

    # Performance metrics
    exports.add(
        obj=fs.tsa.purity,  # Purity
        name="Purity",
        ui_units=units.dimensionless,
        display_units="-",
        rounding=5,
        input_category="Performance Metrics",
        output_category="Performance Metrics",
        is_input=False,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.recovery,  # Recovery
        name="Recovery",
        ui_units=units.dimensionless,
        display_units="-",
        rounding=5,
        input_category="Performance Metrics",
        output_category="Performance Metrics",
        is_input=False,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.specific_energy,
        name="Specific energy",
        ui_units=units.MJ / units.kg,
        display_units="MJ/kg CO2",
        rounding=4,
        input_category="Performance Metrics",
        output_category="Performance Metrics",
        is_input=False,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.productivity,
        name="Productivity",
        ui_units=units.kg / units.metric_ton / units.h,
        display_units="kg CO2/ton/h",
        rounding=3,
        input_category="Performance Metrics",
        output_category="Performance Metrics",
        is_input=False,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.heat_duty_total,
        name="Heat duty total",
        ui_units=units.MW,
        display_units="MW",
        rounding=8,
        output_category="Performance Metrics",
        is_input=False,
        is_output=True,
    )

    # CO2 Capture Performance
    exports.add(
        obj=fs.tsa.CO2_captured_bed_cycle,  # CO2 captured per cycle per bed
        name="CO2 captured in one cycle per bed",
        ui_units=units.kg,
        display_units="kg/cycle",
        rounding=6,
        output_category="Performance Metrics",
        is_input=False,
        is_output=True,
    )

    # CO2 Emissions
    exports.add(
        obj=fs.tsa.emissions_co2,
        name="Amount of CO2 to atmosphere",
        ui_units=units.mol / units.s,
        display_units="mol/s",
        rounding=6,
        output_category="Performance Metrics",
        is_input=False,
        is_output=True,
    )

    exports.add(
        obj=fs.tsa.emissions_co2_ppm,  # CO2 concentration in emissions
        name="Concentration of CO2 emitted to atmosphere",
        ui_units=units.dimensionless,
        display_units="ppm",
        rounding=6,
        output_category="Performance Metrics",
        is_input=False,
        is_output=True,
    )


def initialize(fs, **kwargs):
    """
    Initialize the model.

    Return:
        solver_options
    """

    solver_options = {
        "nlp_scaling_method": "user-scaling",
        "tol": 1e-6,
    }

    initializer = FixedBedTSA0DInitializer(
        output_level=idaes_log.INFO, solver_options=solver_options
    )

    initializer.initialize(fs.tsa)
    return solver_options


def solve(flowsheet=None, **kwargs):
    """
    Solve the model.

    Return:
        solver results.
    """
    fs = flowsheet

    iutil.scaling.calculate_scaling_factors(fs.tsa)

    # XXX: skip when do_initialize is in FlowsheetInterface ctor
    solver_options = initialize(fs)

    solver = SolverFactory("ipopt")
    solver.options = solver_options

    results = solver.solve(fs.model())

    tsa_summary(fs.tsa)

    return results
