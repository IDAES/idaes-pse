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
This is an example supercritical pulverized coal (SCPC) power plant, including
steam cycle and boiler heat exchanger network model.

This simulation model consist of a ~595 MW gross coal fired power plant.
The dimensions and operating conditions used for this simulation do not
represent any specific coal-fired power plant.

This model is for demonstration and tutorial purposes only.
Before looking at the model, it may be useful to look
at the process flow diagram (PFD).

SCPC Power Plant

Inputs:
    Fresh Water (water make up)
    Throttle valve opening,
    BFW - boiler feed water (from Feed water heaters)
    Coal from pulverizers

Main Assumptions:
    Coal flowrate as a function of load, coal HHV is fixed and heat dutty
    splitt from fire side to water wall and platen superheater is fixed.

    Boiler heat exchanger network:
        Water Flow:
            Fresh water -> FWH's -> Economizer -> Water Wall -> Primary SH -> Platen SH -> Finishing Superheate -> HP Turbine -> Reheater -> IP Turbine
        Flue Gas Flow:
            Fire Ball -> Platen SH -> Finishing SH -> Reheater  -> o -> Economizer -> Air Preheater
                                                   -> Primary SH --^
        Steam Flow:
            Boiler -> HP Turbine -> Reheater -> IP Turbine
            HP, IP, and LP steam extractions to Feed Water Heaters


    Models used:
        - Mixers: Attemperator, Flue gas mix
        - Heater: Platen SH, Fire/Water side (simplified model),
                  Feed Water Heaters, Hot Tank, Condenser
        - BoilerHeatExchanger: Economizer, Primary SH, Finishing SH, Reheater
            + Shell and tube heat exchanger
                - tube side: Steam (side 1 holdup)
                - shell side: flue gas (side 2 holdup)
        - Steam Turbines
        - Pumps
    Property packages used:
        - IAPWS: Water/steam side
        - IDEAL GAS: Flue Gas side

"""

__author__ = "Miguel Zamarripa"

# Import Python libraries
import logging

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc

# IDAES Imports
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

_log = idaeslog.getModelLogger(__name__, logging.INFO)


def import_steam_cycle():
    # build concrete model
    # import steam cycle model and initialize flowsheet
    import idaes.models_extra.power_generation.flowsheets.supercritical_steam_cycle.supercritical_steam_cycle as steam_cycle

    m, solver = steam_cycle.main()
    return m, solver


def main():
    # import steam cycle and build concrete model
    m, solver = import_steam_cycle()
    print(degrees_of_freedom(m))
    # at this point we have a flowsheet with "steam cycle" that solves
    # correctly, with 0 degrees of freedom.

    # next step is to import and build the boiler heat exchanger network
    # importing the boiler heat exchanger network
    # from (boiler_subflowsheet_build.py)
    # this step appends all the boiler unit models into our model ("m")
    # model "m" has been created a few lines above
    import idaes.models_extra.power_generation.flowsheets.supercritical_power_plant.boiler_subflowsheet_build as blr

    # import the models (ECON, WW, PrSH, PlSH, FSH, Spliter, Mixer, Reheater)
    # see boiler_subflowhseet_build.py for a beter description
    blr.build_boiler(m.fs)
    # initialize boiler network models (one by one)
    blr.initialize(m)
    # at this point we have both flowsheets (steam cycle + boiler network)
    # in the same model/concrete object ("m"), however they are disconnected.
    # Here we want to solve them at the same time
    # this is a square problem (i.e. degrees of freedom = 0)
    print("solving square problem disconnected")
    results = solver.solve(m, tee=True)

    # at this point we want to connect the units in both flowsheets
    # Economizer inlet = Feed water heater 8 outlet (water)
    # HP inlet = Attemperator outlet (steam)
    # Reheater inlet (steam) = HP split 7 outlet (last stage of HP turbine)
    # IP inlet = Reheater outlet steam7
    blr.unfix_inlets(m)
    print("unfix inlet conditions, degreeso of freedom = " + str(degrees_of_freedom(m)))
    # user can save the initialization to a json file (uncomment next line)
    #    MS.to_json(m, fname = 'SCPC_full.json')
    #   later user can use the json file to initialize the model
    #   if this is the case comment out previous MS.to_json and uncomment next line
    #    MS.from_json(m, fname = 'SCPC_full.json')

    # deactivate constraints linking the FWH8 to HP turbine
    m.fs.boiler_pressure_drop.deactivate()
    m.fs.close_flow.deactivate()
    m.fs.turb.constraint_reheat_flow.deactivate()
    m.fs.turb.constraint_reheat_press.deactivate()
    m.fs.turb.constraint_reheat_temp.deactivate()
    m.fs.turb.inlet_split.inlet.enth_mol.unfix()
    m.fs.turb.inlet_split.inlet.pressure.unfix()
    # user can fix the boiler feed water pump pressure (uncomenting next line)
    #    m.fs.bfp.outlet.pressure[:].fix(26922222.222))

    m.fs.FHWtoECON = Arc(
        source=m.fs.fwh8.desuperheat.cold_side_outlet,
        destination=m.fs.ECON.cold_side_inlet,
    )

    m.fs.Att2HP = Arc(source=m.fs.ATMP1.outlet, destination=m.fs.turb.inlet_split.inlet)

    m.fs.HPout2RH = Arc(
        source=m.fs.turb.hp_split[7].outlet_1, destination=m.fs.RH.cold_side_inlet
    )

    m.fs.RHtoIP = Arc(
        source=m.fs.RH.cold_side_outlet, destination=m.fs.turb.ip_stages[1].inlet
    )

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # unfix boiler connections
    m.fs.ECON.cold_side_inlet.flow_mol.unfix()
    m.fs.ECON.cold_side_inlet.enth_mol[0].unfix()
    m.fs.ECON.cold_side_inlet.pressure[0].unfix()
    m.fs.RH.cold_side_inlet.flow_mol.unfix()
    m.fs.RH.cold_side_inlet.enth_mol[0].unfix()
    m.fs.RH.cold_side_inlet.pressure[0].unfix()
    m.fs.hotwell.makeup.flow_mol[:].setlb(-1.0)

    # if user has trouble with infeasible solutions, an easy test
    # is to deactivate the link to HP turbine
    # (m.fs.Att2HP_expanded "enth_mol and pressure" equalities)
    # and fix inlet pressure and enth_mol to turbine
    # (m.fs.turb.inlet_split.inlet)
    # (then double check the values from m.fs.ATMP1.outlet)
    #  m.fs.Att2HP_expanded.enth_mol_equality.deactivate()
    #  m.fs.Att2HP_expanded.pressure_equality.deactivate()
    m.fs.turb.inlet_split.inlet.pressure.fix(2.423e7)
    #    m.fs.turb.inlet_split.inlet.enth_mol.fix(62710.01)

    # finally, since we want to maintain High Pressure (HP) inlet temperature
    # constant (~866 K), we need to fix Attemperator enthalpy outlet
    # and unfix heat duty to Platen superheater, note that fixing enthalpy
    # to control temperature is only valid because pressure is also fixed
    m.fs.ATMP1.outlet.enth_mol[0].fix(62710.01)
    m.fs.PlSH.heat_duty[:].unfix()  # fix(5.5e7)
    #    m.fs.ATMP1.SprayWater.flow_mol[0].unfix()
    print("connecting flowsheets, degrees of freedom = " + str(degrees_of_freedom(m)))
    print("solving full plant model")
    solver.options = {
        "tol": 1e-6,
        "linear_solver": "ma27",
        "max_iter": 40,
    }
    # square problems tend to work better without bounds
    strip_bounds = pyo.TransformationFactory("contrib.strip_var_bounds")
    strip_bounds.apply_to(m, reversible=True)
    # this is the final solve with both flowsheets connected
    results = solver.solve(m, tee=True)
    return m, results


if __name__ == "__main__":
    m, results = main()
