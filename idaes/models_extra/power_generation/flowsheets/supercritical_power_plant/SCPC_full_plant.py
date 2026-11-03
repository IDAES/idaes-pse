#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
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
            Fresh water -> FWH's -> Economizer -> Water Wall -> Primary SH -> Platen SH
              -> Finishing Superheate -> HP Turbine -> Reheater -> IP Turbine
        Flue Gas Flow:
            Fire Ball -> Platen SH -> Finishing SH -> Reheater  -> o -> Economizer
            -> Air Preheater                        -> Primary SH --^
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


Author: Miguel Zamarripa
Last edited: 2024-06-12 by Francis Hanna
"""

# pylint: disable=missing-function-docstring

__author__ = "Miguel Zamarripa"


# Import Python libraries
import logging

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc

# IDAES Imports
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

# Import supporting IDAES flowsheets
# Steam cycle flowsheet
from idaes.models_extra.power_generation.flowsheets.supercritical_steam_cycle import (
    supercritical_steam_cycle as steam_cycle,
)

# Boiler heat exchanger network flowsheet
from idaes.models_extra.power_generation.flowsheets.supercritical_power_plant import (
    boiler_subflowsheet_build as blr,
)

# setup logger
_log = idaeslog.getModelLogger(__name__, logging.INFO)


def main():
    """
    This function builds the supercritical power plant model, including the steam cycle
    and boiler heat exchanger network.

    Logic
    -----
    1.  Build the steam cycle flowsheet using the setup_steam_cycle function.
        See: idaes/models_extra/power_generation/flowsheets/supercritical_steam_cycle.py
        for more details

    2.  Build the boiler heat exchanger network flowsheet using the build_boiler function.
        See the ``boiler_subflowsheet_build.py`` module for more details.
        This step involves appending the boiler unit models into the steam cycle model,
        resulting in a single model object containing both flowsheets.


    Notes
    -----
    1.  The initialized model may be saved to ``SCPC_full.json`` using
        ``MS.to_json()`` and restored in a later run using ``MS.from_json()``.

    2.  If the connected model produces an infeasible solution, the high-pressure
        turbine connection can be tested by deactivating the enthalpy and pressure
        equalities on ``m.fs.Att2HP_expanded`` and fixing the corresponding inlet
        conditions on ``m.fs.turb.inlet_split.inlet``. The fixed values should be
        checked against the outlet conditions from ``m.fs.ATMP1``.

    3.  To maintain the high-pressure turbine inlet temperature at approximately
        866 K, fix the attemperator outlet molar enthalpy and unfix the platen
        superheater heat duty. Using enthalpy to control temperature is valid here
        only when the corresponding pressure is also fixed.
    """
    # Build the steam cycle flowsheet
    _log.info("Building steam cycle flowsheet")
    m, solver = steam_cycle.main()
    # Assert that the degrees of freedom = 0, i.e. square problem
    _log.info(f"Degrees of freedom: {degrees_of_freedom(m)}")

    # Build the boiler heat exchanger network flowsheet
    _log.info("Building boiler heat exchanger network flowsheet")
    blr.build_boiler(m.fs)
    # initialize boiler network models
    _log.info("Initializing boiler heat exchanger network flowsheet")
    blr.initialize(m)

    # Model content: two disconnected square flowsheets
    # solve the disconnected flowsheets
    _log.info("Solving square problem disconnected")
    # this function is imported from the boiler_subflowsheet_build.py module
    # the scale_solve function is a wrapper for the solver.solve function,
    # which scales the model and solves it. Scaler logic - uses autosclaer to
    # scale the entire flowsheet without overwriting any predefined scaling factors
    results = blr.scale_solve(m)

    # Connect the two flowsheets - the main connection points include:
    # 1. Economizer inlet = Feed water heater 8 outlet (water)
    # 2. HP inlet = Attemperator outlet (steam)
    # 3. Reheater inlet (steam) = HP split 7 outlet (last stage of HP turbine)
    # 4. IP inlet = Reheater outlet steam7

    # Step 1 - unfix the inlet conditions of the boiler heat exchanger network flowsheet
    _log.info("Unfixing boiler subflowsheet inlet conditions")
    blr.unfix_inlets(m)
    _log.info(f"Degrees of freedom: {degrees_of_freedom(m)}")
    # Step 2 - deactivate constraints linking the FWH8 to HP turbine
    _log.info("Deactivating boiler heat exchanger network constraints")
    m.fs.boiler_pressure_drop.deactivate()
    m.fs.close_flow.deactivate()
    m.fs.turb.constraint_reheat_flow.deactivate()
    m.fs.turb.constraint_reheat_press.deactivate()
    m.fs.turb.constraint_reheat_temp.deactivate()
    m.fs.turb.inlet_split.inlet.enth_mol.unfix()
    m.fs.turb.inlet_split.inlet.pressure.unfix()
    # user can fix the boiler feed water pump pressure (uncommenting next line)
    #    m.fs.bfp.outlet.pressure[:].fix(26922222.222))

    # Step 3 - connect the two flowsheets using Arcs
    # FWH8 outlet to ECON inlet
    _log.info(
        "Creating new arcs to connect the steam cycle and boiler heat exchanger "
        "network flowsheets"
    )
    _log.info("Creating arc from FWH8 to ECON")
    m.fs.FHWtoECON = Arc(
        source=m.fs.fwh8.desuperheat.cold_side_outlet,
        destination=m.fs.ECON.cold_side_inlet,
    )
    # HP inlet to Attemperator outlet
    _log.info("Creating arc from Attemperator to HP")
    m.fs.Att2HP = Arc(source=m.fs.ATMP1.outlet, destination=m.fs.turb.inlet_split.inlet)
    # Reheater inlet to HP split 7 outlet
    _log.info("Creating arc from HP split 7 to Reheater")
    m.fs.HPout2RH = Arc(
        source=m.fs.turb.hp_split[7].outlet_1, destination=m.fs.RH.cold_side_inlet
    )
    # IP inlet to Reheater outlet
    _log.info("Creating arc from Reheater to IP")
    m.fs.RHtoIP = Arc(
        source=m.fs.RH.cold_side_outlet, destination=m.fs.turb.ip_stages[1].inlet
    )
    # expand the newly created arcs
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # Step 4: unfix boiler connections
    _log.info("Unfixing boiler connections")
    m.fs.ECON.cold_side_inlet.flow_mol.unfix()
    m.fs.ECON.cold_side_inlet.enth_mol[0].unfix()
    m.fs.ECON.cold_side_inlet.pressure[0].unfix()
    m.fs.RH.cold_side_inlet.flow_mol.unfix()
    m.fs.RH.cold_side_inlet.enth_mol[0].unfix()
    m.fs.RH.cold_side_inlet.pressure[0].unfix()
    m.fs.hotwell.makeup.flow_mol[:].setlb(-1.0)

    m.fs.turb.inlet_split.inlet.pressure.fix(2.423e7)
    #    m.fs.turb.inlet_split.inlet.enth_mol.fix(62710.01)

    # Adjust attemperator and platen superheater to maintain HP inlet temperature
    # Refer to Note (3)
    m.fs.ATMP1.outlet.enth_mol[0].fix(62710.01)
    m.fs.PlSH.heat_duty[:].unfix()  # fix(5.5e7)
    # m.fs.ATMP1.SprayWater.flow_mol[0].unfix()

    _log.info(
        "connecting flowsheets, degrees of freedom = " + str(degrees_of_freedom(m))
    )
    _log.info("solving full plant model")
    # setup solver options
    solver.options = {
        "tol": 1e-6,
        "linear_solver": "ma27",
        "max_iter": 40,
    }
    # square problems tend to work better without bounds
    strip_bounds = pyo.TransformationFactory("contrib.strip_var_bounds")
    strip_bounds.apply_to(m, reversible=True)
    # this is the final solve with both flowsheets connected
    _log.info("Solving the model with both flowsheets connected")
    results = blr.scale_solve(m)

    return m, results


if __name__ == "__main__":
    m, results = main()
