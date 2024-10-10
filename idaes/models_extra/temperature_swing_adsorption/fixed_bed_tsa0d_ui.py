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

"""
UI exports for 0D Fixed Bed TSA unit model.
"""
from watertap.ui import fsapi as api
from pyomo.environ import ConcreteModel, units
from idaes.core.solvers import get_solver
from idaes.core import FlowsheetBlock
from idaes.models_extra.temperature_swing_adsorption import (
    FixedBedTSA0D,
    Adsorbent,
    FixedBedTSA0DInitializer,
    SteamCalculationType,
    TransformationScheme,
)

desc = """The Fixed Bed TSA0d model is implemented as a Temperature Swing Adsorption
(TSA) cycle. The model is an 0D equilibrium-based shortcut model. The model
assumes local adsorption equilibrium and takes into account heat transfer
mechanisms but neglects mass transfer resistances.""".replace(
    "\n", " "
)

unit_name = "FixedBedTSA0D"


def export_to_ui() -> api.FlowsheetInterface:
    fsi = api.FlowsheetInterface(
        name="0D Fixed Bed TSA",
        description=desc,
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
                "value": 120,
                "values_allowed": "int",
                "category": unit_name,
            },
        },
    )
    # TO-DO: The wrapper
    return fsi


def build(build_options=None, **kwargs):
    unit_opt = {k: v for k, v in build_options.items() if v.category == unit_name}
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.unit = FixedBedTSA0D(
        adsorbent=Adsorbent[unit_opt["adsorbent"].value],
        number_of_beds=unit_opt["number_of_beds"].value,
        transformation_method="dae.collocation",
        transformation_scheme=TransformationScheme.lagrangeRadau,
        finite_elements=20,
        collocation_points=6,
        compressor=False,
        steam_calculation=SteamCalculationType.none,
    )

    m.fs.unit.inlet.flow_mol_comp[0, "H2O"].fix(0)
    m.fs.unit.inlet.flow_mol_comp[0, "CO2"].fix(40)
    m.fs.unit.inlet.flow_mol_comp[0, "N2"].fix(99960)
    m.fs.unit.inlet.flow_mol_comp[0, "O2"].fix(0)
    m.fs.unit.inlet.temperature.fix(303.15)
    m.fs.unit.inlet.pressure.fix(100000)

    m.fs.unit.temperature_desorption.fix(470)
    m.fs.unit.temperature_adsorption.fix(310)
    m.fs.unit.temperature_heating.fix(500)
    m.fs.unit.temperature_cooling.fix(300)
    m.fs.unit.bed_diameter.fix(4)
    m.fs.unit.bed_height.fix(8)

    return m


def export(flowsheet=None, exports=None, build_options=None, **kwargs):
    fs = flowsheet
    # Feeds
    for compound in "H2O", "CO2", "N2", "O2":
        exports.add(
            obj=fs.unit.inlet.flow_mol_comp[0, compound],
            name=f"{compound} molar flow rate",
            ui_units=units.mole / units.s,
            display_units="mol/s",
            rounding=2,
            description=f"Inlet molar flow rate for {compound}",
            is_input=True,
            input_category="Feed",
            is_output=False,
        )
    exports.add(
        obj=fs.unit.inlet.temperature[0],
        name="Temperature",
        ui_units=units.K,
        display_units="K",
        is_input=True,
        input_category="Feed",
    )
    exports.add(
        obj=fs.unit.inlet.pressure[0],
        name="Pressure",
        ui_units=units.bar,
        display_units="bar",
        is_input=True,
        input_category="Feed",
    )
    # TODO: export all the variables
    # Parameters
    # Products


def initialize(fs, **kwargs):
    """Initialize the model."""
    initializer = FixedBedTSA0DInitializer()
    initializer.initialize(fs.unit, heating_time_guess=70000, cooling_time_guess=110000)


solver = get_solver()


def solve(flowsheet=None, **kwargs):
    """Solve the model."""
    fs = flowsheet
    initialize(fs)  # XXX: skip when do_initialize is in FlowsheetInterface ctor
    results = solver.solve(fs.model())
    return results
