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
This is an example of a subcritical boiler recirculation system. This model
does not represent any specific power plant, it shows a generic recirculation
system. This model is for demonstration and tutorial purposes only.

This example consist of a simulation of the internal recirculation system
of a subcritical coal fired boiler. The simulation inlet includes the
feed water flowrate and enthalpy, the unit dimensions/specifications,
and boiler load, which in this case is represented as fixed heat flux from
fire side.

The main outputs of the model are the feedwater pressure (to control the boiler
feedwater pump), and steam outlet to superheaters.

For a detailed description see Jupyter Notebook
authors: Boiler Subsystem Team (J. Ma, M. Zamarripa)
"""
import os

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import model_serializer as ms
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

# Import Unit Model Modules
from idaes.models.properties import iapws95

# Import IDAES standard unit model
import idaes.logger as idaeslog
from idaes.models_extra.power_generation.unit_models.drum import Drum
from idaes.models_extra.power_generation.unit_models.downcomer import Downcomer
from idaes.models_extra.power_generation.unit_models.waterwall_section import (
    WaterwallSection,
)
from idaes.models_extra.power_generation.properties.flue_gas_ideal import (
    FlueGasParameterBlock,
)

# import plotting libraries
import matplotlib.pyplot as plt

_log = idaeslog.getModelLogger(__name__)


def main(m=None):
    """Create concrete model, make the flowsheet object, fix some variables,
    and solve the problem.
    if m concrete model already exists this function only builds the unit
    models, fixes the inputs, and initializes the units.
    This is useful to import subsections of a flowsheet
    (i.e. boiler subsystem + steam cycle) if concrete model already exists
    we assume m.fs flowsheet, m.fs.prop_water and m.fs.prop_gas also exist"""
    if m is None:
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        # Add property packages to flowsheet library
        m.fs.prop_water = iapws95.Iapws95ParameterBlock()
        m.fs.prop_gas = FlueGasParameterBlock()

    create_model(m)
    set_inputs(m)
    initialize(m)
    return m


def create_model(m):
    """Create unit models"""

    m.fs.ww_zones = pyo.RangeSet(10)

    m.fs.drum = Drum(
        property_package=m.fs.prop_water,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
    )
    m.fs.downcomer = Downcomer(
        property_package=m.fs.prop_water, has_holdup=False, has_heat_transfer=True
    )
    m.fs.Waterwalls = WaterwallSection(
        m.fs.ww_zones,
        dynamic=False,
        has_holdup=False,
        property_package=m.fs.prop_water,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    m.fs.stream_drum_out = Arc(
        source=m.fs.drum.liquid_outlet, destination=m.fs.downcomer.inlet
    )
    m.fs.stream_dcmr_out = Arc(
        source=m.fs.downcomer.outlet, destination=m.fs.Waterwalls[1].inlet
    )

    def arc_rule(b, i):
        return {
            "source": m.fs.Waterwalls[i].outlet,
            "destination": m.fs.Waterwalls[i + 1].inlet,
        }

    m.arc = Arc(pyo.RangeSet(9), rule=arc_rule)

    m.fs.stream_ww14 = Arc(
        source=m.fs.Waterwalls[10].outlet, destination=m.fs.drum.water_steam_inlet
    )
    # pyomo arcs expand constraints for inlet = outlet ports
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)


def set_inputs(m):
    """fix variables for geometry and design data"""

    # drum
    m.fs.drum.drum_diameter.fix(1.2)
    m.fs.drum.drum_length.fix(15.3256)
    m.fs.drum.number_downcomers.fix(6)
    m.fs.drum.downcomer_diameter.fix(0.38)
    m.fs.drum.drum_level[:].fix(0.6)  # drum level, set as unit radius
    m.fs.drum.heat_duty[:].fix(0.0)  # assume no heat loss
    # downcomer inputs
    m.fs.downcomer.diameter.fix(0.3)
    m.fs.downcomer.number_downcomers.fix(6)  # number of downcomers
    m.fs.downcomer.height.fix(40.1)
    m.fs.downcomer.heat_duty[:].fix(0.0)

    # 10 waterwall sections
    for i in m.fs.ww_zones:
        m.fs.Waterwalls[i].tube_diameter.fix(0.047)
        m.fs.Waterwalls[i].tube_thickness.fix(0.00350)
        m.fs.Waterwalls[i].fin_thickness.fix(0.00455)
        m.fs.Waterwalls[i].slag_thickness[:].fix(0.001)
        m.fs.Waterwalls[i].fin_length.fix(0.0115)
        m.fs.Waterwalls[i].number_tubes.fix(610)
        m.fs.Waterwalls[i].fcorrection_dp.fix(1.2)

    # water wall section height (must be same as fire side model zones)
    m.fs.Waterwalls[1].height.fix(6.150)
    m.fs.Waterwalls[2].height.fix(3.150)
    m.fs.Waterwalls[3].height.fix(1.5)
    m.fs.Waterwalls[4].height.fix(1.450)
    m.fs.Waterwalls[5].height.fix(1.350)
    m.fs.Waterwalls[6].height.fix(1.250)
    m.fs.Waterwalls[7].height.fix(1.150)
    m.fs.Waterwalls[8].height.fix(1.350)
    m.fs.Waterwalls[9].height.fix(3.250)
    m.fs.Waterwalls[10].height.fix(3.450)

    # water wall section projected area
    m.fs.Waterwalls[1].projected_area.fix(320.0)
    m.fs.Waterwalls[2].projected_area.fix(150.3)
    m.fs.Waterwalls[3].projected_area.fix(70.8)
    m.fs.Waterwalls[4].projected_area.fix(70.0)
    m.fs.Waterwalls[5].projected_area.fix(58.6)
    m.fs.Waterwalls[6].projected_area.fix(58.6)
    m.fs.Waterwalls[7].projected_area.fix(50.1)
    m.fs.Waterwalls[8].projected_area.fix(65.6)
    m.fs.Waterwalls[9].projected_area.fix(145.6)
    m.fs.Waterwalls[10].projected_area.fix(165.5)


def initialize(
    m,
    outlvl=idaeslog.NOTSET,
    optarg={
        "tol": 1e-6,
        "max_iter": 40,
    },
):
    """Initialize unit models"""
    init_log = idaeslog.getInitLogger(m.name, outlvl, tag="flowsheet")
    solve_log = idaeslog.getSolveLogger(m.name, outlvl, tag="flowsheet")

    solver = get_solver()
    solver.options = optarg
    init_log.info_low("Starting initialization...")

    if not os.path.exists("subcritical_boiler_init.json.gz"):
        # 10 Waterwalls, initial guess for specific simulation
        m.fs.Waterwalls[1].heat_fireside[:].fix(2.3e7)
        m.fs.Waterwalls[2].heat_fireside[:].fix(1.5e7)
        m.fs.Waterwalls[3].heat_fireside[:].fix(6.8e6)
        m.fs.Waterwalls[4].heat_fireside[:].fix(1.2e7)
        m.fs.Waterwalls[5].heat_fireside[:].fix(1.2e7)
        m.fs.Waterwalls[6].heat_fireside[:].fix(1.2e7)
        m.fs.Waterwalls[7].heat_fireside[:].fix(1.0e7)
        m.fs.Waterwalls[8].heat_fireside[:].fix(9.9e6)
        m.fs.Waterwalls[9].heat_fireside[:].fix(2.1)
        m.fs.Waterwalls[10].heat_fireside[:].fix(2.0e7)

        state_args_water_steam = {
            "flow_mol": 199470.7831,  # mol/s
            "pressure": 10903981.9107,  # Pa
            "enth_mol": 26585.3483,
        }  # j/mol

        state_args_feedwater = {
            "flow_mol": 4630.6098,  # mol/s
            "pressure": 10903981.9107,  # Pa
            "enth_mol": 22723.907,
        }  # j/mol

        m.fs.drum.initialize(
            outlvl=outlvl,
            optarg=solver.options,
            state_args_water_steam=state_args_water_steam,
            state_args_feedwater=state_args_feedwater,
        )
        m.fs.downcomer.inlet.flow_mol[:].fix(m.fs.drum.liquid_outlet.flow_mol[0].value)
        m.fs.downcomer.inlet.pressure[:].fix(m.fs.drum.liquid_outlet.pressure[0].value)
        m.fs.downcomer.inlet.enth_mol[:].fix(m.fs.drum.liquid_outlet.enth_mol[0].value)

        m.fs.downcomer.initialize(
            state_args={
                "flow_mol": m.fs.drum.liquid_outlet.flow_mol[0].value,
                "pressure": m.fs.drum.liquid_outlet.pressure[0].value,
                "enth_mol": m.fs.drum.liquid_outlet.enth_mol[0].value,
            },
            outlvl=outlvl,
            optarg=solver.options,
        )

        m.fs.Waterwalls[1].inlet.flow_mol[:].fix(
            m.fs.downcomer.outlet.flow_mol[0].value
        )
        m.fs.Waterwalls[1].inlet.pressure[:].fix(
            m.fs.downcomer.outlet.pressure[0].value
        )
        m.fs.Waterwalls[1].inlet.enth_mol[:].fix(
            m.fs.downcomer.outlet.enth_mol[0].value
        )
        m.fs.Waterwalls[1].initialize(
            state_args={
                "flow_mol": m.fs.Waterwalls[1].inlet.flow_mol[0].value,
                "pressure": m.fs.Waterwalls[1].inlet.pressure[0].value,
                "enth_mol": m.fs.Waterwalls[1].inlet.enth_mol[0].value,
            },
            outlvl=6,
            optarg=solver.options,
        )

        for i in range(2, 11):
            m.fs.Waterwalls[i].inlet.flow_mol[:].fix(
                m.fs.Waterwalls[i - 1].outlet.flow_mol[0].value
            )
            m.fs.Waterwalls[i].inlet.pressure[:].fix(
                m.fs.Waterwalls[i - 1].outlet.pressure[0].value
            )
            m.fs.Waterwalls[i].inlet.enth_mol[:].fix(
                m.fs.Waterwalls[i - 1].outlet.enth_mol[0].value
            )
            m.fs.Waterwalls[i].initialize(
                state_args={
                    "flow_mol": m.fs.Waterwalls[i - 1].outlet.flow_mol[0].value,
                    "pressure": m.fs.Waterwalls[i - 1].outlet.pressure[0].value,
                    "enth_mol": m.fs.Waterwalls[i - 1].outlet.enth_mol[0].value,
                },
                outlvl=6,
                optarg=solver.options,
            )

        m.fs.drum.feedwater_inlet.flow_mol[:].fix()
        m.fs.drum.feedwater_inlet.pressure[:].unfix()
        m.fs.drum.feedwater_inlet.enth_mol[:].fix()

        m.fs.downcomer.inlet.flow_mol[:].unfix()
        m.fs.downcomer.inlet.pressure[:].unfix()
        m.fs.downcomer.inlet.enth_mol[:].unfix()
        print("solving full-space problem")
        for i in m.fs.ww_zones:
            m.fs.Waterwalls[i].inlet.flow_mol[:].unfix()
            m.fs.Waterwalls[i].inlet.pressure[:].unfix()
            m.fs.Waterwalls[i].inlet.enth_mol[:].unfix()

        df = degrees_of_freedom(m)
        if df != 0:
            raise ValueError("Check degrees of freedom: {}".format(df))
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            print("solving full-space problem")
            res = solver.solve(m, tee=slc.tee)
        init_log.info_low("Initialization Complete: {}".format(idaeslog.condition(res)))
        ms.to_json(m, fname="subcritical_boiler_init.json.gz")
    else:
        print("\n\nInitializing from json file")
        ms.from_json(m, fname="subcritical_boiler_init.json.gz")
    return m


def run_sensitivity():
    m = main()
    optarg = {"tol": 1e-6, "max_iter": 40}

    m.fs.drum.feedwater_inlet.flow_mol[:].fix()
    m.fs.drum.feedwater_inlet.pressure[:].unfix()
    m.fs.drum.feedwater_inlet.enth_mol[:].fix()
    # set scaling parameters
    for i in m.fs.ww_zones:
        iscale.set_scaling_factor(m.fs.Waterwalls[i].heat_flux_conv[0], 1e-5)
    iscale.calculate_scaling_factors(m)
    # solve flowsheet
    solver = get_solver()
    solver.options = optarg
    results = solver.solve(m, tee=True)
    print(results)

    print("\nDrum feedwater inlet")
    m.fs.drum.feedwater_inlet.display()
    print("\nWaterWall Outlet")
    m.fs.drum.water_steam_inlet.display()
    m.fs.Waterwalls[10].outlet.display()
    print("\ndrum flash inlet")
    m.fs.drum.flash.inlet.display()
    print("\ndrum liquid outlet")
    m.fs.drum.liquid_outlet.display()
    print("\ndowncomer inlet")
    m.fs.downcomer.inlet.display()
    print("\ndowncomer outlet")
    m.fs.downcomer.outlet.display()
    print("\nwaterwall inlet")
    m.fs.Waterwalls[1].inlet.display()
    print("\n vapor fraction zone 10")
    m.fs.Waterwalls[10].control_volume.properties_out[0].vapor_frac.display()

    # Since heat transfer indirectly calculates the recirculation flowrate
    # and drum level is fixed in this case, the feedwater flowrate
    # can be calculated as a function of the boiler load, by changing the
    # heat duty to waterwall.
    # This sensitivity analysis demonstrate that increasing the boiler load
    # coal flowrate (in this case heat duty), increases the feedwater flowrate
    # and the steam outlet
    vap_frac = []
    flow = []
    flow_it = [
        1.1,
        1.2,
        1.3,
        1.4,
        1.5,
        1.6,
        1.7,
        1.8,
        1.9,
        2.0,
        2.1,
        2.2,
        2.3,
        2.4,
        2.5,
        2.6,
        2.7,
        2.8,
        2.9,
        3.0,
    ]
    count = 0
    for i in flow_it:
        count = count + 1
        m.fs.Waterwalls[1].heat_fireside[:].fix(2.3e7 * i)
        m.fs.Waterwalls[2].heat_fireside[:].fix(1.5e7 * i)
        m.fs.Waterwalls[3].heat_fireside[:].fix(6.8e6 * i)
        m.fs.Waterwalls[4].heat_fireside[:].fix(1.2e7 * i)
        m.fs.Waterwalls[5].heat_fireside[:].fix(1.2e7 * i)
        m.fs.Waterwalls[6].heat_fireside[:].fix(1.2e7 * i)
        m.fs.Waterwalls[7].heat_fireside[:].fix(1.0e7 * i)
        m.fs.Waterwalls[8].heat_fireside[:].fix(9.9e6 * i)
        m.fs.Waterwalls[9].heat_fireside[:].fix(2.1e7 * i)
        m.fs.Waterwalls[10].heat_fireside[:].fix(2.0e7 * i)
        m.fs.drum.feedwater_inlet.flow_mol[:].unfix()
        m.fs.drum.feedwater_inlet.pressure[:].fix()
        m.fs.drum.feedwater_inlet.enth_mol[:].fix()
        results = solver.solve(m, tee=False)
        vap_frac.append(
            pyo.value(m.fs.Waterwalls[10].control_volume.properties_out[0].vapor_frac)
        )
        flow.append(pyo.value(m.fs.drum.feedwater_inlet.flow_mol[0]))

        if pyo.check_optimal_termination(results):
            print("iter " + str(count) + " = EXIT- Optimal Solution Found.")
        else:
            print("iter " + str(count) + " = infeasible")

    x = []
    for i in flow_it:
        x.append((i - 1) * 100)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x, flow, label=str("feedwater (steam) flowrate"))
    plt.legend()
    ax.set_xlabel("heat duty - % increment")
    ax.set_ylabel("Flowrate mol/s")
    plt.show()
    return m
