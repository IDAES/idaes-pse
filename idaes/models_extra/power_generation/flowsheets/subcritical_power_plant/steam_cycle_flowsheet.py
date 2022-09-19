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
Dynamic sub-flowsheet for a subcritical 300MWe steam cycle system
"""
# Import Python time library
import time
import matplotlib.pyplot as plt

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc

# Import IDAES
from idaes.core.util.initialization import propagate_state as _set_port
from idaes.models_extra.power_generation.unit_models.helm import (
    HelmTurbineStage as TurbineStage,
    HelmTurbineOutletStage as TurbineOutletStage,
    ValveFunctionType,
    HelmTurbineMultistage as TurbineMultistage,
    HelmMixer as Mixer,
    MomentumMixingType,
    HelmValve as SteamValve,
    HelmValve as WaterValve,
    HelmIsentropicCompressor as WaterPump,
    HelmSplitter as Separator,
    HelmNtuCondenser as Condenser,
)
from idaes.models_extra.power_generation.unit_models import (
    WaterTank,
    FWH0DDynamic as FWH0D,
)
from idaes.models.control.controller import (
    PIDController,
    ControllerType,
    ControllerMVBoundType,
)

from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.core.util.dyn_utils import copy_values_at_time, copy_non_time_indexed_values

__author__ = "J. Ma, M. Zamarripa, J. Eslick"


def add_unit_models(m):
    """Function to add unit models to the steam cycle sub-flowsheet"""
    fs = m.fs_main.fs_stc
    prop_water = m.fs_main.prop_water

    # Custom valve function for throttle valve
    def throttle_valve_function(blk):
        blk.Cv.fix(1)
        a = blk.vfa = pyo.Var(initialize=2.8904e-02, doc="Valve function parameter a")
        b = blk.vfb = pyo.Var(initialize=3.3497e-02, doc="Valve function parameter b")
        c = blk.vfc = pyo.Var(initialize=1.4514e-02, doc="Valve function parameter c")
        d = blk.vfd = pyo.Var(initialize=1.4533e-03, doc="Valve function parameter d")
        a.fix()
        b.fix()
        c.fix()
        d.fix()
        o = blk.valve_opening

        @blk.Expression(m.fs_main.time)
        def valve_function(bd, t):
            return a * o[t] ** 3 - b * o[t] ** 2 + c * o[t] - d

    # Unit model for multistage turbine including throttle valve
    fs.turb = TurbineMultistage(
        dynamic=False,
        property_package=prop_water,
        num_parallel_inlet_stages=1,
        throttle_valve_function=ValveFunctionType.custom,
        throttle_valve_function_callback=throttle_valve_function,
        num_hp=14,
        num_ip=9,
        num_lp=5,
        hp_split_locations=[14],
        ip_split_locations=[6, 9],
        lp_split_locations=[2, 4, 5],
        hp_disconnect=[14],
        hp_split_num_outlets={14: 2},
        ip_split_num_outlets={9: 3},
    )

    # Unit model for regulating valve of BFPT (boiler feed pump turbine)
    fs.bfp_turb_valve = SteamValve(dynamic=False, property_package=prop_water)

    # Unit model for main stage of BFPT
    fs.bfp_turb = TurbineStage(dynamic=False, property_package=prop_water)

    # Unit model for outlet stage of BFPT
    fs.bfp_turb_os = TurbineOutletStage(dynamic=False, property_package=prop_water)

    # Unit model for main condenser
    fs.condenser = Condenser(
        dynamic=False,
        shell={"has_pressure_change": False, "property_package": prop_water},
        tube={"has_pressure_change": False, "property_package": prop_water},
    )

    # Unit model for auxiliary condenser
    fs.aux_condenser = Condenser(
        dynamic=False,
        shell={"has_pressure_change": False, "property_package": prop_water},
        tube={"has_pressure_change": False, "property_package": prop_water},
    )

    # Unit model for condenser hotwell (hotwell tank modeled separately)
    # Modeled as a mixer of makeup, main, and auxiliary condenser water streams
    # Set momentum_mixing_type to none since aux_condenser outlet pressure
    # is usually not equal to the main condenser outlet pressure
    # We impose the constraints to let the mixed pressure equal to
    # the main condenser pressure and makeup water pressure
    fs.condenser_hotwell = Mixer(
        dynamic=False,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["main_condensate", "makeup", "aux_condensate"],
        property_package=prop_water,
    )

    # Unit model for water control valve between makeup tank and hotwell
    fs.makeup_valve = WaterValve(
        dynamic=False, has_holdup=False, phase="Liq", property_package=prop_water
    )

    # Unit model for hotwell tank with holdup for dynamic model
    # Modeled as a simple tank with constant cross section area and tank level
    fs.hotwell_tank = WaterTank(
        tank_type="simple_tank", has_holdup=True, property_package=prop_water
    )

    # Unit model for condensate pump
    fs.cond_pump = WaterPump(dynamic=False, property_package=prop_water)

    # Unit model for water control valve after hotwell tank
    # Used to control deaerator level
    fs.cond_valve = WaterValve(
        dynamic=False, has_holdup=False, phase="Liq", property_package=prop_water
    )

    # Unit model for feed water heater 1
    fs.fwh1 = FWH0D(
        has_desuperheat=False,
        has_drain_cooling=False,
        has_drain_mixer=True,
        condense={
            "tube": {"has_pressure_change": True},
            "shell": {"has_pressure_change": True},
            "has_holdup": True,
        },
        property_package=prop_water,
    )

    # Unit model for drain pump of FWH1
    fs.fwh1_drain_pump = WaterPump(dynamic=False, property_package=prop_water)

    # Unit model for mixer of FWH1 drain and condensate
    fs.fwh1_drain_return = Mixer(
        dynamic=False,
        inlet_list=["feedwater", "fwh1_drain"],
        property_package=prop_water,
        momentum_mixing_type=MomentumMixingType.equality,
    )

    # Unit model for feed water heater 2
    fs.fwh2 = FWH0D(
        has_desuperheat=False,
        has_drain_cooling=True,
        has_drain_mixer=True,
        condense={
            "tube": {"has_pressure_change": True},
            "shell": {"has_pressure_change": True},
            "has_holdup": True,
        },
        cooling={"dynamic": False, "has_holdup": False},
        property_package=prop_water,
    )

    # Unit model for water control valve between drain of fwh2 and fwh1
    fs.fwh2_valve = WaterValve(
        dynamic=False, has_holdup=False, phase="Liq", property_package=prop_water
    )

    # Unit model for feed water heater 3
    fs.fwh3 = FWH0D(
        has_desuperheat=False,
        has_drain_cooling=True,
        has_drain_mixer=False,
        condense={
            "tube": {"has_pressure_change": True},
            "shell": {"has_pressure_change": True},
            "has_holdup": True,
        },
        cooling={"dynamic": False, "has_holdup": False},
        property_package=prop_water,
    )

    # Unit model for control valve between drain of fwh3 and fwh2
    fs.fwh3_valve = WaterValve(
        dynamic=False, has_holdup=False, phase="Liq", property_package=prop_water
    )

    # Unit model for deaerator also known as fwh4
    # Modeled as mixer to mix extracted steam with condensate and drain
    # from fwh5
    # Using MomentumMixingType.equality for momentum_mixing_type
    # deaerator tank modeled separately for holdup in dyyamic model
    fs.fwh4_deair = Mixer(
        dynamic=False,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["steam", "drain", "feedwater"],
        property_package=prop_water,
    )

    # Unit model for deaerator water tank
    # Modeled as a horizontal cylindrical tank
    fs.da_tank = WaterTank(
        tank_type="horizontal_cylindrical_tank",
        has_holdup=True,
        property_package=prop_water,
    )

    # Unit model for electrical feedwater booster pump
    fs.booster = WaterPump(dynamic=False, property_package=prop_water)

    # Unit model for main boiler feed water pump driven by steam turbine
    fs.bfp = WaterPump(dynamic=False, property_package=prop_water)

    # Unit model for splitter for spray water stream for main attemperator
    fs.split_attemp = Separator(
        dynamic=False, property_package=prop_water, outlet_list=["FeedWater", "Spray"]
    )

    # Unit model for attemperator spray control valve
    fs.spray_valve = WaterValve(
        dynamic=False, has_holdup=False, phase="Liq", property_package=prop_water
    )

    # Unit model for feed water heater 5
    fs.fwh5 = FWH0D(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=True,
        condense={
            "tube": {"has_pressure_change": True},
            "shell": {"has_pressure_change": True},
            "has_holdup": True,
        },
        desuperheat={"dynamic": False},
        cooling={"dynamic": False, "has_holdup": False},
        property_package=prop_water,
    )

    # Unit model for water control valve drain of fwh5 and deaerator
    fs.fwh5_valve = WaterValve(
        dynamic=False, has_holdup=False, phase="Liq", property_package=prop_water
    )

    # Unit model for feed water heater 6
    fs.fwh6 = FWH0D(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=False,
        condense={
            "tube": {"has_pressure_change": True},
            "shell": {"has_pressure_change": True},
            "has_holdup": True,
        },
        desuperheat={"dynamic": False},
        cooling={"dynamic": False, "has_holdup": False},
        property_package=prop_water,
    )

    # Unit model for water control valve between drain of fwh6 and fwh5
    fs.fwh6_valve = WaterValve(
        dynamic=False, has_holdup=False, phase="Liq", property_package=prop_water
    )

    # Important process variables, declared and used in PID controllers
    # Variable for main steam temperature
    fs.temperature_main_steam = pyo.Var(fs.time, initialize=810)

    # Constraint to calculate main steam temperature
    @fs.Constraint(fs.time)
    def temperature_main_steam_eqn(b, t):
        return (
            b.temperature_main_steam[t]
            == b.turb.throttle_valve[1].control_volume.properties_in[t].temperature
        )

    # Variable for gross power output in MW
    fs.power_output = pyo.Var(fs.time, initialize=300, doc="gross power output in MW")

    # Constraint to calculate gross power output
    @fs.Constraint(fs.time)
    def power_output_eqn(b, t):
        return b.power_output[t] == -b.turb.power[t] / 1e6

    if m.dynamic is True:
        # Add PID controllers if the flowsheet model is a dynamic model
        # PI controller to control level of fwh2
        fs.fwh2_ctrl = PIDController(
            process_var=fs.fwh2.condense.level,
            manipulated_var=fs.fwh2_valve.valve_opening,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )

        # PI controller to control level of fwh3
        fs.fwh3_ctrl = PIDController(
            process_var=fs.fwh3.condense.level,
            manipulated_var=fs.fwh3_valve.valve_opening,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )

        # PI controller to control level of fwh5
        fs.fwh5_ctrl = PIDController(
            process_var=fs.fwh5.condense.level,
            manipulated_var=fs.fwh5_valve.valve_opening,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )

        # PI controller to control level of fwh6
        fs.fwh6_ctrl = PIDController(
            process_var=fs.fwh6.condense.level,
            manipulated_var=fs.fwh6_valve.valve_opening,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )

        # PI controller to control level of deaerator tank
        fs.da_ctrl = PIDController(
            process_var=fs.da_tank.tank_level,
            manipulated_var=fs.cond_valve.valve_opening,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )

        # PI controller to control level of hotwell tank
        fs.makeup_ctrl = PIDController(
            process_var=fs.hotwell_tank.tank_level,
            manipulated_var=fs.makeup_valve.valve_opening,
            type=ControllerType.PI,
            mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
            calculate_initial_integral=False,
        )

        # PID controller to control main steam temperature
        fs.spray_ctrl = PIDController(
            process_var=fs.temperature_main_steam,
            manipulated_var=fs.spray_valve.valve_opening,
            type=ControllerType.PID,
            mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
            calculate_initial_integral=False,
        )

    return m


def set_arcs_and_constraints(m):
    """Add arcs to connect streams on steam cycle sub-flowsheet"""
    fs = m.fs_main.fs_stc
    fs.S017 = Arc(
        source=fs.turb.outlet_stage.outlet, destination=fs.condenser.hot_side_inlet
    )
    fs.S053 = Arc(source=fs.bfp_turb.outlet, destination=fs.bfp_turb_os.inlet)
    fs.S046 = Arc(
        source=fs.bfp_turb_os.outlet, destination=fs.aux_condenser.hot_side_inlet
    )
    fs.S022 = Arc(
        source=fs.condenser.hot_side_outlet,
        destination=fs.condenser_hotwell.main_condensate,
    )
    fs.S047 = Arc(
        source=fs.aux_condenser.hot_side_outlet,
        destination=fs.condenser_hotwell.aux_condensate,
    )
    fs.S050b = Arc(
        source=fs.makeup_valve.outlet,
        destination=fs.condenser_hotwell.makeup,
    )
    fs.S031a = Arc(
        source=fs.fwh1.condense.cold_side_outlet,
        destination=fs.fwh1_drain_return.feedwater,
    )
    fs.S030 = Arc(
        source=fs.fwh1.condense.hot_side_outlet, destination=fs.fwh1_drain_pump.inlet
    )
    fs.S045 = Arc(
        source=fs.fwh1_drain_pump.outlet, destination=fs.fwh1_drain_return.fwh1_drain
    )
    fs.S024 = Arc(source=fs.condenser_hotwell.outlet, destination=fs.hotwell_tank.inlet)
    fs.S025 = Arc(source=fs.hotwell_tank.outlet, destination=fs.cond_pump.inlet)
    fs.S016 = Arc(
        source=fs.turb.lp_split[5].outlet_2, destination=fs.fwh1.drain_mix.steam
    )
    fs.S032 = Arc(
        source=fs.fwh2.cooling.hot_side_outlet, destination=fs.fwh2_valve.inlet
    )
    fs.S032b = Arc(source=fs.fwh2_valve.outlet, destination=fs.fwh1.drain_mix.drain)
    fs.S026 = Arc(source=fs.cond_pump.outlet, destination=fs.cond_valve.inlet)
    fs.S029 = Arc(
        source=fs.cond_valve.outlet, destination=fs.fwh1.condense.cold_side_inlet
    )
    fs.S031b = Arc(
        source=fs.fwh1_drain_return.outlet, destination=fs.fwh2.cooling.cold_side_inlet
    )
    fs.S015 = Arc(
        source=fs.turb.lp_split[4].outlet_2, destination=fs.fwh2.drain_mix.steam
    )
    fs.S034 = Arc(
        source=fs.fwh3.cooling.hot_side_outlet, destination=fs.fwh3_valve.inlet
    )
    fs.S034b = Arc(source=fs.fwh3_valve.outlet, destination=fs.fwh2.drain_mix.drain)
    fs.S033 = Arc(
        source=fs.fwh2.condense.cold_side_outlet,
        destination=fs.fwh3.cooling.cold_side_inlet,
    )
    fs.S014 = Arc(
        source=fs.turb.lp_split[2].outlet_2, destination=fs.fwh3.condense.hot_side_inlet
    )
    fs.S035 = Arc(
        source=fs.fwh3.condense.cold_side_outlet, destination=fs.fwh4_deair.feedwater
    )
    fs.S043 = Arc(source=fs.turb.ip_split[9].outlet_2, destination=fs.fwh4_deair.steam)
    fs.S011 = Arc(
        source=fs.turb.ip_split[9].outlet_3, destination=fs.bfp_turb_valve.inlet
    )
    fs.S052 = Arc(source=fs.bfp_turb_valve.outlet, destination=fs.bfp_turb.inlet)
    fs.S036 = Arc(source=fs.fwh4_deair.outlet, destination=fs.da_tank.inlet)
    fs.S036b = Arc(source=fs.da_tank.outlet, destination=fs.booster.inlet)
    fs.S038 = Arc(source=fs.booster.outlet, destination=fs.bfp.inlet)
    fs.S037 = Arc(source=fs.bfp.outlet, destination=fs.split_attemp.inlet)
    fs.S054 = Arc(source=fs.split_attemp.Spray, destination=fs.spray_valve.inlet)
    fs.S039 = Arc(
        source=fs.fwh5.cooling.hot_side_outlet, destination=fs.fwh5_valve.inlet
    )
    fs.S039b = Arc(source=fs.fwh5_valve.outlet, destination=fs.fwh4_deair.drain)
    fs.S051 = Arc(
        source=fs.split_attemp.FeedWater, destination=fs.fwh5.cooling.cold_side_inlet
    )
    fs.S010 = Arc(
        source=fs.turb.ip_split[6].outlet_2,
        destination=fs.fwh5.desuperheat.hot_side_inlet,
    )
    fs.S041 = Arc(
        source=fs.fwh6.cooling.hot_side_outlet, destination=fs.fwh6_valve.inlet
    )
    fs.S041b = Arc(source=fs.fwh6_valve.outlet, destination=fs.fwh5.drain_mix.drain)
    fs.S040 = Arc(
        source=fs.fwh5.desuperheat.cold_side_outlet,
        destination=fs.fwh6.cooling.cold_side_inlet,
    )
    fs.S006 = Arc(
        source=fs.turb.hp_split[14].outlet_2,
        destination=fs.fwh6.desuperheat.hot_side_inlet,
    )
    # Call Pyomo function to apply above arc connections
    pyo.TransformationFactory("network.expand_arcs").apply_to(fs)

    # Constraints on the steam cycle sub-flowsheet
    #
    # Constraint to set boiler feed water pump work equal to BFPT work
    # Note that there are to unit models for BFPT including an outlet stage
    @fs.Constraint(fs.time)
    def constraint_bfp_power(b, t):
        return 0 == (
            1e-6 * fs.bfp.control_volume.work[t]
            + 1e-6 * fs.bfp_turb.control_volume.work[t]
            + 1e-6 * fs.bfp_turb_os.control_volume.work[t]
        )

    # Constraint to set IP inlet steam flow equal to HP outlet steam flow
    @fs.turb.Constraint(fs.time)
    def constraint_reheat_flow(b, t):
        return b.ip_stages[1].inlet.flow_mol[t] == b.hp_split[14].outlet_1.flow_mol[t]

    # The mixer for the condenser hotwell is declared with 3 inlet streams
    # with mometum mixing type set to MomentumMixingType.none
    # Here we use two constraints to set the mixed state pressure to
    # the pressure of the auxiliary condenser pressure
    #
    # Constraint to set the pressure of makeup water to hotwell equal to
    # the pressure of auxiliary condenser, which is usually lower than
    # the pressure of main condenser in this model
    @fs.Constraint(fs.time)
    def makeup_water_pressure_constraint(b, t):
        return (
            b.condenser_hotwell.makeup_state[t].pressure * 1e-4
            == b.condenser_hotwell.aux_condensate_state[t].pressure * 1e-4
        )

    # Constrait to set the mixed state pressure equal to the pressure of
    # auxiliary condenser
    @fs.condenser_hotwell.Constraint(fs.time)
    def mixer_pressure_constraint(b, t):
        return (
            b.aux_condensate_state[t].pressure * 1e-4
            == b.mixed_state[t].pressure * 1e-4
        )

    # Constraint to set deaerator tank outlet enthalpy equal to
    # saturation enthalpy at inlet - 100 (sligtly sub-cooled)
    # This constraint determines the steam extraction flow rate and
    # is very important to avoid flash of deaerator tank when load is
    # ramping down, which will causes convergence issue if the flash happens
    @fs.Constraint(fs.time)
    def da_outlet_enthalpy_constraint(b, t):
        return 1e-3 * b.da_tank.outlet.enth_mol[t] == 1e-3 * (
            b.da_tank.control_volume.properties_in[t].enth_mol_sat_phase["Liq"] - 100
        )

    # Constraints to set the pressures of the two inlets (extracte steam
    # and drain from upper stream feed water heater) equal since the
    # drain_mix model of the feed water heater is declared with momentum
    # mixing type set to MomentumMixingType.none
    # Constraint for feed water heater 1
    @fs.Constraint(fs.time)
    def fwh1_drain_mixer_pressure_eqn(b, t):
        return (
            b.fwh1.drain_mix.drain.pressure[t] * 1e-4
            == b.fwh1.drain_mix.steam.pressure[t] * 1e-4
        )

    # Constraints to set the pressures of the two inlets (extracte steam
    # and drain from upper stream feed water heater) equal since the
    # drain_mix model of the feed water heater is declared with momentum
    # mixing type set to MomentumMixingType.none
    # Constraint for feed water heater 2
    @fs.Constraint(fs.time)
    def fwh2_drain_mixer_pressure_eqn(b, t):
        return (
            b.fwh2.drain_mix.drain.pressure[t] * 1e-4
            == b.fwh2.drain_mix.steam.pressure[t] * 1e-4
        )

    # Constraints to set the pressures of the two inlets (extracte steam
    # and drain from upper stream feed water heater) equal since the
    # drain_mix model of the feed water heater is declared with momentum
    # mixing type set to MomentumMixingType.none
    # Constraint for feed water heater 5
    @fs.Constraint(fs.time)
    def fwh5_drain_mixer_pressure_eqn(b, t):
        return (
            b.fwh5.drain_mix.drain.pressure[t] * 1e-5
            == b.fwh5.drain_mix.steam.pressure[t] * 1e-5
        )

    # Add a custom constraint for condensate pump using the pump curve,
    # the relationship between flow rate and pressure increase
    @fs.cond_pump.Constraint(fs.time)
    def cond_pump_curve_constraint(b, t):
        return 1e-5 * b.deltaP[t] == 1e-5 * (
            -5.08e-7 * b.inlet.flow_mol[t] ** 3
            + 4.09e-3 * b.inlet.flow_mol[t] ** 2
            - 56.6 * b.inlet.flow_mol[t]
            + 2e6
        )

    # Add a custom constraint for the booster pump using the pump curve,
    # the relationship between flow rate and pressure increase
    @fs.booster.Constraint(fs.time)
    def booster_pump_curve_constraint(b, t):
        return 1e-5 * b.deltaP[t] == 1e-5 * (
            -2.36e-7 * b.inlet.flow_mol[t] ** 3
            + 3.15e-3 * b.inlet.flow_mol[t] ** 2
            - 34.1 * b.inlet.flow_mol[t]
            + 1.1e6
        )

    # Constraint to set water flow to economizer slightly higher than
    # steam flow in such that there is makeup flow needed to offset blowdown
    # This constraint will be deactivated when the steam cycle sub-flowsheet
    # is combined with the boiler system sub-flowsheet
    @fs.Constraint(fs.time)
    def fw_flow_constraint(b, t):
        return (
            1e-3 * b.da_tank.outlet.flow_mol[t]
            == 1e-3 * b.turb.inlet_split.inlet.flow_mol[t] * 1.02
        )

    # The following expressions are used to calculate drain cooler
    # approach temperature (DCA)
    def rule_dca_no_cool(b, t):
        return (
            b.condense.hot_side.properties_out[t].temperature
            - b.condense.cold_side.properties_in[t].temperature
        )

    def rule_dca(b, t):
        return (
            b.cooling.hot_side.properties_out[t].temperature
            - b.cooling.cold_side.properties_in[t].temperature
        )

    fs.fwh1.dca = pyo.Expression(fs.time, rule=rule_dca_no_cool)
    fs.fwh2.dca = pyo.Expression(fs.time, rule=rule_dca)
    fs.fwh3.dca = pyo.Expression(fs.time, rule=rule_dca)
    fs.fwh5.dca = pyo.Expression(fs.time, rule=rule_dca)
    fs.fwh6.dca = pyo.Expression(fs.time, rule=rule_dca)

    return m


def set_inputs(m):
    """
    Set unit model inputs related to dimensions, pararameters,
    and fixed design and operating variables
    Some of fixed inputs are for initialization only and will be unfixed
    before the sub-flowsheet is solved
    """
    fs = m.fs_main.fs_stc

    # Multi-stage turbine inputs
    #
    # Throttle valve inputs
    fs.turb.throttle_valve[1].Cv.fix(0.7)
    fs.turb.throttle_valve[1].valve_opening.fix(0.8)

    # Inlet stage inputs
    # Inlet stage flow coefficient
    fs.turb.turbine_inlet_cf_fix(0.0015)
    fs.turb.inlet_stage[1].efficiency_mech.fix(0.99)
    fs.turb.inlet_stage[1].eff_nozzle.fix(0.93)
    fs.turb.inlet_stage[1].blade_reaction.fix(0.9)

    # Outlet stage inputs
    # Outlet stage flow coefficient
    fs.turb.turbine_outlet_cf_fix(0.0567)
    fs.turb.outlet_stage.efficiency_mech.fix(0.99)
    # Efficiency if steam is dry
    fs.turb.outlet_stage.eff_dry.fix(0.93)
    # Design volumetric flow rate
    fs.turb.outlet_stage.design_exhaust_flow_vol.fix(4900)

    # Set the efficencies and pressure ratios of stages
    # other than inlet and outlet stages
    for i, s in fs.turb.hp_stages.items():
        s.efficiency_mech.fix(0.99)
        s.ratioP.fix(0.930)
        s.efficiency_isentropic.fix(0.93)
    for i, s in fs.turb.ip_stages.items():
        s.efficiency_mech.fix(0.99)
        if i <= 6:
            s.ratioP.fix(0.87)
            s.efficiency_isentropic.fix(0.92)
        else:
            s.ratioP.fix(0.83)
            s.efficiency_isentropic.fix(0.96)
    for i, s in fs.turb.lp_stages.items():
        s.efficiency_mech.fix(0.99)
        if i <= 2:
            s.ratioP[:].fix(0.65)
            s.efficiency_isentropic[:].fix(0.865)
        elif i <= 4:
            s.ratioP[:].fix(0.70)
            s.efficiency_isentropic[:].fix(0.855)
        else:
            s.ratioP[:].fix(0.35)
            s.efficiency_isentropic[:].fix(0.835)

    # Set main steam inlet condition
    pin = 1.25e7
    # Steam enthalpy corresponding to 810 K at 12.5 MPa
    hin = 61958
    fs.turb.inlet_split.inlet.enth_mol[:].fix(hin)
    fs.turb.inlet_split.inlet.flow_mol[:].value = 11000
    fs.turb.inlet_split.inlet.flow_mol[:].unfix()
    fs.turb.inlet_split.inlet.pressure[:].fix(pin)
    # Set inlet mixer up for pressure driven flow
    fs.turb.inlet_mix.use_equal_pressure_constraint()

    # Set reheater steam inlet condition
    pin = 2.8e6
    # Steam enthalpy corresponding to 810 K at 2.8 MPa
    hin = 63754
    fs.turb.ip_stages[1].inlet.enth_mol[:].fix(hin)
    fs.turb.ip_stages[1].inlet.pressure[:].fix(pin)
    # Reheating steam flow is 90% of main steam flow as initial guess
    fs.turb.ip_stages[1].inlet.flow_mol[:].value = 11000 * 0.9

    # Fix split fractions for steam extraction initialization
    # They will be unfixed after initialization
    # Split to FWH6
    fs.turb.hp_split[14].split_fraction[:, "outlet_2"].fix(0.02)
    # Split to FWH5
    fs.turb.ip_split[6].split_fraction[:, "outlet_2"].fix(0.02)
    # Split to deaerator, use small value to avoid flash of deaerator
    fs.turb.ip_split[9].split_fraction[:, "outlet_2"].fix(0.01)
    # Split to BFPT
    fs.turb.ip_split[9].split_fraction[:, "outlet_3"].fix(0.04)
    # Split to FWH3
    fs.turb.lp_split[2].split_fraction[:, "outlet_2"].fix(0.02)
    # Split to FWH2
    fs.turb.lp_split[4].split_fraction[:, "outlet_2"].fix(0.02)
    # Split to FWH1
    fs.turb.lp_split[5].split_fraction[:, "outlet_2"].fix(0.02)

    # Set inputs for control valve of BFPT
    fs.bfp_turb_valve.Cv.fix(0.001475)
    fs.bfp_turb_valve.valve_opening.fix(0.85)

    # Set inputs for bfp turbine before outlet stage
    fs.bfp_turb.efficiency_mech.fix(0.99)
    fs.bfp_turb.efficiency_isentropic.fix(0.85)
    fs.bfp_turb.ratioP[:].fix(0.073)

    # Set inputs for bfp turbine outlet stage
    fs.bfp_turb_os.efficiency_mech.fix(0.99)
    fs.bfp_turb_os.eff_dry.fix(0.93)
    fs.bfp_turb_os.design_exhaust_flow_vol.fix(180)
    fs.bfp_turb_os.flow_coeff.fix(0.0035)

    # Set inputs for main condenser
    # Cooling water condition
    fs.condenser.cold_side_inlet.flow_mol.fix(300000)
    # Enthalpy at 24 C
    fs.condenser.cold_side_inlet.enth_mol.fix(1800)
    fs.condenser.cold_side_inlet.pressure.fix(500000)
    fs.condenser.area.fix(13000)
    fs.condenser.overall_heat_transfer_coefficient.fix(3100)

    # Set inputs for auxiliary condenser
    fs.aux_condenser.cold_side_inlet.flow_mol.fix(12000)
    fs.aux_condenser.cold_side_inlet.enth_mol.fix(1800)
    fs.aux_condenser.cold_side_inlet.pressure.fix(500000)
    fs.aux_condenser.area.fix(350)
    fs.aux_condenser.overall_heat_transfer_coefficient.fix(3100)

    # Set inputs for makeup water valve
    fs.makeup_valve.Cv.value = 2.1039
    fs.makeup_valve.valve_opening.fix(0.5)

    # Set makeup water control valve
    # Assume the inlet pressure is 1 bar
    fs.makeup_valve.inlet.flow_mol.fix(220)
    fs.makeup_valve.inlet.enth_mol.fix(1890)
    fs.makeup_valve.inlet.pressure.fix(1e5)

    # Set hotwell tank inputs
    fs.hotwell_tank.tank_cross_sect_area.fix(90)
    fs.hotwell_tank.tank_level.fix(0.7)

    # Set condensate pump inputs
    fs.cond_pump.efficiency_isentropic.fix(0.80)
    fs.cond_pump.deltaP[:].value = 1.4e6

    # Set condensate valve inputs, the valve is used to control da_tank level
    fs.cond_valve.Cv.value = 18.75
    fs.cond_valve.valve_opening.fix(0.85)

    # Set inputs for FWH1
    fs.fwh1.condense.area.fix(580)
    fs.fwh1.condense.overall_heat_transfer_coefficient.fix(950)
    fs.fwh1.condense.tube.deltaP[:].fix(0)
    # Inputs required for dynamic model
    fs.fwh1.condense.level.fix(0.275)
    fs.fwh1.condense.heater_diameter.fix(1.3)
    fs.fwh1.condense.vol_frac_shell.fix(0.75)
    fs.fwh1.condense.cond_sect_length.fix(8.0)
    fs.fwh1.condense.tube.volume.fix(1.5)

    # Set inputs for FWH1 drain pump
    fs.fwh1_drain_pump.deltaP[:].value = 7e5
    fs.fwh1_drain_pump.efficiency_isentropic.fix(0.8)

    # Set inputs for FWH2
    fs.fwh2.condense.area.fix(580)
    fs.fwh2.cooling.area.fix(70)
    fs.fwh2.condense.overall_heat_transfer_coefficient.fix(1100)
    fs.fwh2.cooling.overall_heat_transfer_coefficient.fix(675)
    fs.fwh2.condense.tube.deltaP[:].fix(0)
    # Inputs required for dynamic model
    fs.fwh2.condense.level.fix(0.275)
    fs.fwh2.condense.heater_diameter.fix(1.3)
    fs.fwh2.condense.vol_frac_shell.fix(0.7)
    fs.fwh2.condense.cond_sect_length.fix(7.5)
    fs.fwh2.condense.tube.volume.fix(1.5)

    # Set inputs for level control valve after FWH2
    fs.fwh2_valve.Cv.value = 6.878
    fs.fwh2_valve.valve_opening.fix(0.5)

    # Set inputs for FWH3
    fs.fwh3.condense.area.fix(650)
    fs.fwh3.cooling.area.fix(65)
    fs.fwh3.condense.overall_heat_transfer_coefficient.fix(1200)
    fs.fwh3.cooling.overall_heat_transfer_coefficient.fix(620)
    fs.fwh3.condense.tube.deltaP[:].fix(0)
    # Inputs required for dynamic model
    fs.fwh3.condense.level.fix(0.275)
    fs.fwh3.condense.heater_diameter.fix(1.3)
    fs.fwh3.condense.vol_frac_shell.fix(0.7)
    fs.fwh3.condense.cond_sect_length.fix(7.5)
    fs.fwh3.condense.tube.volume.fix(1.5)

    # Set inputs for level control valve after FWH3
    fs.fwh3_valve.Cv.value = 2.5358
    fs.fwh3_valve.valve_opening.fix(0.5)

    # Set inputs for deaerator tank
    fs.da_tank.tank_diameter.fix(4.0)
    fs.da_tank.tank_length.fix(20)
    fs.da_tank.tank_level.fix(2.75)

    # Set inputs for booster pump
    fs.booster.efficiency_isentropic.fix(0.8)
    fs.booster.outlet.pressure.fix(1.5e6)

    # Set inputs for main boiler feed water pump
    fs.bfp.efficiency_isentropic.fix(0.8)
    fs.bfp.outlet.pressure.fix(1.45e7)

    # Set input for spliter for main steam attemperator spray
    fs.split_attemp.split_fraction[:, "Spray"].fix(0.0007)

    # Set inputs for water spray control valve
    fs.spray_valve.Cv.value = 0.331
    fs.spray_valve.valve_opening.fix(0.25)

    # Set inputs for FWH5
    fs.fwh5.condense.area.fix(600)
    fs.fwh5.desuperheat.area.fix(85)
    fs.fwh5.cooling.area.fix(100)
    fs.fwh5.condense.overall_heat_transfer_coefficient.fix(1010)
    fs.fwh5.desuperheat.overall_heat_transfer_coefficient.fix(145)
    fs.fwh5.cooling.overall_heat_transfer_coefficient.fix(675)
    fs.fwh5.condense.tube.deltaP[:].fix(0)
    # Inputs reqired for dynamic model
    fs.fwh5.condense.level.fix(0.275)
    fs.fwh5.condense.heater_diameter.fix(1.4)
    fs.fwh5.condense.vol_frac_shell.fix(0.675)
    fs.fwh5.condense.cond_sect_length.fix(6.4)
    fs.fwh5.condense.tube.volume.fix(1.3)

    # Set inputs for level control valve after FWH5
    fs.fwh5_valve.Cv.value = 5.156
    fs.fwh5_valve.valve_opening.fix(0.5)

    # Set inputs for FWH6
    fs.fwh6.condense.area.fix(750)
    fs.fwh6.desuperheat.area.fix(125)
    fs.fwh6.cooling.area.fix(130.0)
    fs.fwh6.condense.overall_heat_transfer_coefficient.fix(1070)
    fs.fwh6.desuperheat.overall_heat_transfer_coefficient.fix(260)
    fs.fwh6.cooling.overall_heat_transfer_coefficient.fix(635)
    fs.fwh6.condense.tube.deltaP[:].fix(0)
    # Inputs required for dynamic model
    fs.fwh6.condense.level.fix(0.275)
    fs.fwh6.condense.heater_diameter.fix(1.35)
    fs.fwh6.condense.vol_frac_shell.fix(0.662)
    fs.fwh6.condense.cond_sect_length.fix(8.4)
    fs.fwh6.condense.tube.volume.fix(1.6)

    # Set inputs for level control valve after FWH6
    fs.fwh6_valve.Cv.value = 1.9986
    fs.fwh6_valve.valve_opening.fix(0.5)

    if m.dynamic is True:
        # Set PI or PID controller parameters
        fs.fwh2_ctrl.gain_p.fix(-1e-1)
        fs.fwh2_ctrl.gain_i.fix(-2e-2)
        fs.fwh2_ctrl.setpoint.fix(0.275)
        fs.fwh2_ctrl.mv_ref.fix(0.5)

        fs.fwh3_ctrl.gain_p.fix(-1e-1)
        fs.fwh3_ctrl.gain_i.fix(-2e-2)
        fs.fwh3_ctrl.setpoint.fix(0.275)
        fs.fwh3_ctrl.mv_ref.fix(0.5)

        fs.fwh5_ctrl.gain_p.fix(-1e-1)
        fs.fwh5_ctrl.gain_i.fix(-1e-1)
        fs.fwh5_ctrl.setpoint.fix(0.275)
        fs.fwh5_ctrl.mv_ref.fix(0.5)

        fs.fwh6_ctrl.gain_p.fix(-1e-1)
        fs.fwh6_ctrl.gain_i.fix(-1e-1)
        fs.fwh6_ctrl.setpoint.fix(0.275)
        fs.fwh6_ctrl.mv_ref.fix(0.5)

        fs.da_ctrl.gain_p.fix(1e-1)
        fs.da_ctrl.gain_i.fix(1e-1)
        fs.da_ctrl.setpoint.fix(2.75)
        fs.da_ctrl.mv_ref.fix(0.5)

        fs.makeup_ctrl.gain_p.fix(0.01)
        fs.makeup_ctrl.gain_i.fix(0.001)
        fs.makeup_ctrl.setpoint.fix(0.7)
        fs.makeup_ctrl.mv_ref.fix(0.35)

        fs.spray_ctrl.gain_p.fix(-5e-2)
        fs.spray_ctrl.gain_i.fix(-1e-4)
        fs.spray_ctrl.gain_d.fix(-1e-4)
        fs.spray_ctrl.setpoint.fix(810)
        fs.spray_ctrl.mv_ref.fix(0.25)
        # Currently we have to set this minimum opening to avoid converging to
        # negative spray water flow rate
        fs.spray_ctrl.mv_lb = 0.05

        # Set initial conditions for controller errors
        t0 = fs.time.first()
        fs.fwh2_ctrl.integral_of_error[t0].fix(0)
        fs.fwh3_ctrl.integral_of_error[t0].fix(0)
        fs.fwh5_ctrl.integral_of_error[t0].fix(0)
        fs.fwh6_ctrl.integral_of_error[t0].fix(0)
        fs.da_ctrl.integral_of_error[t0].fix(0)
        fs.makeup_ctrl.integral_of_error[t0].fix(0)
        fs.spray_ctrl.integral_of_error[t0].fix(0)
        fs.spray_ctrl.derivative_of_error[t0].fix(0)

    return m


def _add_heat_transfer_correlation(fs):
    """
    Add overall heat transfer coefficient calculation constraints
    based on feed water flow rate
    """
    _add_u_eq(fs.fwh6.condense)
    _add_u_eq(fs.fwh5.condense)
    _add_u_eq(fs.fwh3.condense)
    _add_u_eq(fs.fwh2.condense)
    _add_u_eq(fs.fwh1.condense)


def initialize(m):
    """Initialize unit operation models on the sub-flowsheet"""
    fs = m.fs_main.fs_stc
    start_time = time.time()
    outlvl = idaeslog.INFO_LOW
    _log = idaeslog.getLogger(fs.name, outlvl, tag="unit")
    solve_log = idaeslog.getSolveLogger(fs.name, outlvl, tag="unit")
    _log.info("Starting steam cycle initialization...")

    solver = get_solver()

    # Set initial condition for dynamic unit models
    t0 = 0
    if m.dynamic is True:
        t0 = fs.time.first()
        fs.hotwell_tank.set_initial_condition()
        fs.da_tank.set_initial_condition()
        fs.fwh1.set_initial_condition()
        fs.fwh2.set_initial_condition()
        fs.fwh3.set_initial_condition()
        fs.fwh5.set_initial_condition()
        fs.fwh6.set_initial_condition()

    # Initialize main turbine
    fs.turb.outlet_stage.control_volume.properties_out[:].pressure.fix(6000)
    if m.dynamic is False:
        fs.turb.initialize(
            outlvl=outlvl,
            optarg=solver.options,
            calculate_outlet_cf=False,
            calculate_inlet_cf=False,
        )

    # Initialize control valve of BFPT
    _set_port(fs.bfp_turb_valve.inlet, fs.turb.ip_split[9].outlet_3)
    if m.dynamic is False:
        fs.bfp_turb_valve.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize BFPT front stage
    _set_port(fs.bfp_turb.inlet, fs.bfp_turb_valve.outlet)
    if m.dynamic is False:
        fs.bfp_turb.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize BFPT outlet stage
    _set_port(fs.bfp_turb_os.inlet, fs.bfp_turb.outlet)
    fs.bfp_turb_os.control_volume.properties_out[:].pressure.fix(6000)
    if m.dynamic is False:
        fs.bfp_turb_os.initialize(
            calculate_cf=False, outlvl=outlvl, optarg=solver.options
        )
    fs.bfp_turb_os.control_volume.properties_out[:].pressure.unfix()

    # Initialize main condenser
    _set_port(fs.condenser.hot_side_inlet, fs.turb.outlet_stage.outlet)
    fs.turb.outlet_stage.control_volume.properties_out[:].pressure.unfix()
    if m.dynamic is False:
        fs.condenser.initialize(unfix="pressure")

    # Initialize auxiliary condenser
    _set_port(fs.aux_condenser.hot_side_inlet, fs.bfp_turb_os.outlet)
    if m.dynamic is False:
        fs.aux_condenser.initialize(unfix="pressure")

    # Initialize makeup valve
    if m.dynamic is False:
        fs.makeup_valve.Cv.fix()
        fs.makeup_valve.initialize()
        fs.makeup_valve.Cv.unfix()

    # Initialize hotwell mixer
    _set_port(fs.condenser_hotwell.main_condensate, fs.condenser.hot_side_outlet)
    _set_port(fs.condenser_hotwell.aux_condensate, fs.aux_condenser.hot_side_outlet)
    _set_port(fs.condenser_hotwell.makeup, fs.makeup_valve.outlet)
    if m.dynamic is False:
        fs.condenser_hotwell.initialize(outlvl=outlvl, optarg=solver.options)
    fs.condenser_hotwell.main_condensate.unfix()

    # Initialize hotwell tank
    _set_port(fs.hotwell_tank.inlet, fs.condenser_hotwell.outlet)
    if m.dynamic is False:
        fs.hotwell_tank.initialize()

    # Initialize condensate pump
    _set_port(fs.cond_pump.inlet, fs.hotwell_tank.outlet)
    if m.dynamic is False:
        fs.cond_pump.deltaP.fix()
        fs.cond_pump.cond_pump_curve_constraint.deactivate()
        fs.cond_pump.initialize(outlvl=outlvl, optarg=solver.options)
        fs.cond_pump.deltaP.unfix()
        fs.cond_pump.cond_pump_curve_constraint.activate()

    # Initialize condensate valve
    _set_port(fs.cond_valve.inlet, fs.cond_pump.outlet)
    if m.dynamic is False:
        fs.cond_valve.Cv.fix()
        fs.cond_valve.initialize(outlvl=outlvl, optarg=solver.options)
        fs.cond_valve.Cv.unfix()

    # Set some initial inlet values and initialize fwh1
    fs.fwh1.drain_mix.drain.flow_mol[:] = 1000
    fs.fwh1.drain_mix.drain.pressure[:] = (
        fs.turb.lp_split[5].outlet_2.pressure[t0].value
    )
    fs.fwh1.drain_mix.drain.enth_mol[:] = 6117
    _set_port(fs.fwh1.condense.cold_side_inlet, fs.cond_valve.outlet)
    _set_port(fs.fwh1.drain_mix.steam, fs.turb.lp_split[5].outlet_2)

    if m.dynamic is False:
        fs.fwh1.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize fwh1 drain pump
    _set_port(fs.fwh1_drain_pump.inlet, fs.fwh1.condense.hot_side_outlet)
    fs.fwh1_drain_pump.control_volume.properties_out[:].pressure.fix(
        fs.fwh1.condense.cold_side.properties_out[t0].pressure.value
    )
    if m.dynamic is False:
        fs.fwh1_drain_pump.initialize(outlvl=outlvl, optarg=solver.options)
    fs.fwh1_drain_pump.control_volume.properties_out[:].pressure.unfix()

    # Initialize mixer to add fwh1 drain to feedwater
    _set_port(fs.fwh1_drain_return.feedwater, fs.fwh1.condense.cold_side_outlet)
    _set_port(fs.fwh1_drain_return.fwh1_drain, fs.fwh1_drain_pump.outlet)
    if m.dynamic is False:
        fs.fwh1_drain_return.initialize(outlvl=outlvl, optarg=solver.options)
    fs.fwh1_drain_return.feedwater.unfix()
    fs.fwh1_drain_return.fwh1_drain.unfix()

    # Set some initial inlet values and initialize fwh2
    fs.fwh2.drain_mix.drain.flow_mol[:] = 685
    fs.fwh2.drain_mix.drain.pressure[:] = (
        fs.turb.lp_split[4].outlet_2.pressure[t0].value
    )
    fs.fwh2.drain_mix.drain.enth_mol[:] = 9100
    _set_port(fs.fwh2.cooling.cold_side_inlet, fs.fwh1.condense.cold_side_outlet)
    _set_port(fs.fwh2.drain_mix.steam, fs.turb.lp_split[4].outlet_2)
    if m.dynamic is False:
        fs.fwh2.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize fwh2_valve
    _set_port(fs.fwh2_valve.inlet, fs.fwh2.cooling.hot_side_outlet)
    if m.dynamic is False:
        # use a lower flow rate to avoid too low exit pressure
        fs.fwh2_valve.inlet.flow_mol[:].value = (
            fs.fwh2_valve.inlet.flow_mol[t0].value * 0.75
        )
        fs.fwh2_valve.Cv.fix()
        fs.fwh2_valve.initialize(outlvl=outlvl, optarg=solver.options)
        fs.fwh2_valve.Cv.unfix()

    # Set some initial inlet values and initialize fwh3
    _set_port(fs.fwh3.cooling.cold_side_inlet, fs.fwh2.condense.cold_side_outlet)
    _set_port(fs.fwh3.condense.hot_side_inlet, fs.turb.lp_split[2].outlet_2)
    if m.dynamic is False:
        fs.fwh3.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize fwh3_valve
    _set_port(fs.fwh3_valve.inlet, fs.fwh3.cooling.hot_side_outlet)
    if m.dynamic is False:
        fs.fwh3_valve.Cv.fix()
        fs.fwh3_valve.initialize(outlvl=outlvl, optarg=solver.options)
        fs.fwh3_valve.Cv.unfix()

    # Set some initial inlet values and initialize fwh4
    fs.fwh4_deair.drain.flow_mol[:] = 10
    fs.fwh4_deair.drain.pressure[:] = 1.3e6
    fs.fwh4_deair.drain.enth_mol[:] = 13630
    _set_port(fs.fwh4_deair.feedwater, fs.fwh3.condense.cold_side_outlet)
    _set_port(fs.fwh4_deair.steam, fs.turb.ip_split[9].outlet_2)
    if m.dynamic is False:
        fs.fwh4_deair.initialize(outlvl=outlvl, optarg=solver.options)
    fs.fwh4_deair.feedwater.unfix()
    fs.fwh4_deair.steam.unfix()
    fs.fwh4_deair.drain.unfix()

    # Initialize deaerator tank
    _set_port(fs.da_tank.inlet, fs.fwh4_deair.outlet)
    if m.dynamic is False:
        fs.da_tank.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize booster pump
    _set_port(fs.booster.inlet, fs.da_tank.outlet)
    if m.dynamic is False:
        fs.booster.booster_pump_curve_constraint.deactivate()
        fs.booster.initialize(outlvl=outlvl, optarg=solver.options)
        fs.booster.booster_pump_curve_constraint.activate()

    # Initialize main boiler feed water pump
    _set_port(fs.bfp.inlet, fs.booster.outlet)
    if m.dynamic is False:
        fs.bfp.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize split_attemp
    _set_port(fs.split_attemp.inlet, fs.bfp.outlet)
    if m.dynamic is False:
        fs.split_attemp.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize spray_valve
    _set_port(fs.spray_valve.inlet, fs.split_attemp.Spray)
    if m.dynamic is False:
        fs.spray_valve.Cv.fix()
        fs.spray_valve.initialize(outlvl=outlvl, optarg=solver.options)
        fs.spray_valve.Cv.unfix()

    # Initialize fwh5
    fs.fwh5.drain_mix.drain.flow_mol[:] = 50
    fs.fwh5.drain_mix.drain.pressure[:] = 3.5e6
    fs.fwh5.drain_mix.drain.enth_mol[:] = 15000
    _set_port(fs.fwh5.cooling.cold_side_inlet, fs.split_attemp.FeedWater)
    _set_port(fs.fwh5.desuperheat.hot_side_inlet, fs.turb.ip_split[6].outlet_2)
    if m.dynamic is False:
        fs.fwh5.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize fwh5_valve
    _set_port(fs.fwh5_valve.inlet, fs.fwh5.cooling.hot_side_outlet)
    if m.dynamic is False:
        fs.fwh5_valve.Cv.fix()
        fs.fwh5_valve.initialize(outlvl=outlvl, optarg=solver.options)
        fs.fwh5_valve.Cv.unfix()

    # Set some initial inlet values and initialize fwh3
    _set_port(fs.fwh6.cooling.cold_side_inlet, fs.fwh5.desuperheat.cold_side_outlet)
    _set_port(fs.fwh6.desuperheat.hot_side_inlet, fs.turb.hp_split[14].outlet_2)
    if m.dynamic is False:
        fs.fwh6.initialize(outlvl=outlvl, optarg=solver.options)

    # Initialize fwh6_valve
    _set_port(fs.fwh6_valve.inlet, fs.fwh6.cooling.hot_side_outlet)
    if m.dynamic is False:
        fs.fwh6_valve.Cv.fix()
        fs.fwh6_valve.initialize(outlvl=outlvl, optarg=solver.options)
        fs.fwh6_valve.Cv.unfix()

    # Unfix stream connections before solving entire sub-flowsheet
    fs.turb.lp_split[5].split_fraction[:, "outlet_2"].unfix()
    fs.turb.lp_split[4].split_fraction[:, "outlet_2"].unfix()
    fs.turb.lp_split[2].split_fraction[:, "outlet_2"].unfix()
    fs.turb.ip_split[6].split_fraction[:, "outlet_2"].unfix()
    fs.turb.hp_split[14].split_fraction[:, "outlet_2"].unfix()
    fs.turb.ip_split[9].split_fraction[:, "outlet_3"].unfix()
    fs.turb.ip_split[9].split_fraction[:, "outlet_2"].unfix()

    # Set the overall heat transfer coefficients for all FWHs
    # to their designed values
    fs.fwh1.condense.overall_heat_transfer_coefficient.fix(2800)
    fs.fwh2.condense.overall_heat_transfer_coefficient.fix(3250)
    fs.fwh2.cooling.overall_heat_transfer_coefficient.fix(2000)
    fs.fwh3.condense.overall_heat_transfer_coefficient.fix(3600)
    fs.fwh3.cooling.overall_heat_transfer_coefficient.fix(1850)
    fs.fwh5.condense.overall_heat_transfer_coefficient.fix(3050)
    fs.fwh5.desuperheat.overall_heat_transfer_coefficient.fix(450)
    fs.fwh5.cooling.overall_heat_transfer_coefficient.fix(2000)
    fs.fwh6.condense.overall_heat_transfer_coefficient.fix(3200)
    fs.fwh6.desuperheat.overall_heat_transfer_coefficient.fix(780)
    fs.fwh6.cooling.overall_heat_transfer_coefficient.fix(1900)

    # unfix booster pump outlet pressure since it is specified by pump curve
    fs.booster.outlet.pressure.unfix()

    # fix makeup stream inlet pressure and enthalpy but unfix the flow
    fs.makeup_valve.inlet.pressure.fix()
    fs.makeup_valve.inlet.enth_mol.fix()
    fs.makeup_valve.inlet.flow_mol.unfix()

    # unfix spray split fraction, fix back pressure, valve Cv and opening
    # let code calculate spray flow rate
    fs.split_attemp.split_fraction[:, "Spray"].unfix()
    fs.spray_valve.outlet.pressure.fix(1.1e7)
    fs.spray_valve.valve_opening.fix()
    fs.spray_valve.Cv.fix()

    # fix all other valves' Cv's and unfix their openings
    fs.fwh2_valve.Cv.fix()
    fs.fwh3_valve.Cv.fix()
    fs.fwh5_valve.Cv.fix()
    fs.fwh6_valve.Cv.fix()
    fs.cond_valve.Cv.fix()
    fs.makeup_valve.Cv.fix()
    fs.bfp_turb_valve.Cv.fix()
    fs.fwh2_valve.valve_opening.unfix()
    fs.fwh3_valve.valve_opening.unfix()
    fs.fwh5_valve.valve_opening.unfix()
    fs.fwh6_valve.valve_opening.unfix()
    fs.makeup_valve.valve_opening.unfix()
    fs.cond_valve.valve_opening.unfix()
    fs.bfp_turb_valve.valve_opening.unfix()

    if m.dynamic is False:
        # check degree of freedom
        dof = degrees_of_freedom(fs)
        _log.info("Degrees of freedom before solving= {}".format(dof))
        assert dof == 0
        # finally solve the sub-flowsheet
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(fs, tee=slc.tee)
            _log.info("Solving bounded problem: {}".format(idaeslog.condition(res)))
            _log.info(
                "Steam cycle initialized in {:.1f}s".format(time.time() - start_time)
            )

        _log.info(
            "main steam temp={}".format(
                pyo.value(
                    fs.turb.throttle_valve[1]
                    .control_volume.properties_in[0]
                    .temperature
                )
            )
        )
        _log.info(
            "main steam enth={}".format(
                pyo.value(fs.turb.inlet_split.inlet.enth_mol[0])
            )
        )
        _log.info(
            "main steam flow_mol={}".format(
                pyo.value(fs.turb.inlet_split.inlet.flow_mol[0])
            )
        )
        _log.info(
            "main steam pressure={}".format(
                pyo.value(fs.turb.inlet_split.inlet.pressure[0])
            )
        )
        _log.info(
            "bfp turb steam flow={}".format(pyo.value(fs.bfp_turb.inlet.flow_mol[0]))
        )
        _log.info("FW pressure={}".format(pyo.value(fs.bfp.outlet.pressure[0])))
        _log.info(
            "HP stage 1 inlet enth_mol={}".format(
                pyo.value(fs.turb.hp_stages[1].inlet.enth_mol[0])
            )
        )
        _log.info(
            "HP stage 1 inlet flow_mol={}".format(
                pyo.value(fs.turb.hp_stages[1].inlet.flow_mol[0])
            )
        )
        _log.info(
            "HP stage 1 inlet pressure={}".format(
                pyo.value(fs.turb.hp_stages[1].inlet.pressure[0])
            )
        )
        _log.info(
            "IP stage 1 inlet enth_mol={}".format(
                pyo.value(fs.turb.ip_stages[1].inlet.enth_mol[0])
            )
        )
        _log.info(
            "IP stage 1 inlet flow_mol={}".format(
                pyo.value(fs.turb.ip_stages[1].inlet.flow_mol[0])
            )
        )
        _log.info(
            "IP stage 1 inlet pressure={}".format(
                pyo.value(fs.turb.ip_stages[1].inlet.pressure[0])
            )
        )
        _log.info(
            "Outlet stage enth_mol={}".format(
                pyo.value(fs.turb.outlet_stage.outlet.enth_mol[0])
            )
        )
        _log.info(
            "Outlet stage flow_mol={}".format(
                pyo.value(fs.turb.outlet_stage.outlet.flow_mol[0])
            )
        )
        _log.info(
            "Outlet stage pressure={}".format(
                pyo.value(fs.turb.outlet_stage.outlet.pressure[0])
            )
        )
        _log.info(
            "Power output of main turbine={}".format(pyo.value(fs.power_output[0]))
        )
        _log.info(
            "Power output of bfp turbine={}".format(
                -1e-6
                * pyo.value(
                    fs.bfp_turb.control_volume.work[0]
                    + fs.bfp_turb_os.control_volume.work[0]
                )
            )
        )
        _log.info(
            "FWH6 outlet enth_mol={}".format(
                pyo.value(fs.fwh6.desuperheat.cold_side_outlet.enth_mol[0])
            )
        )
        _log.info(
            "FWH6 outlet flow_mol={}".format(
                pyo.value(fs.fwh6.desuperheat.cold_side_outlet.flow_mol[0])
            )
        )
        _log.info(
            "FWH6 outlet pressure={}".format(
                pyo.value(fs.fwh6.desuperheat.cold_side_outlet.pressure[0])
            )
        )
        _log.info(
            "Outlet stage ratioP={}".format(pyo.value(fs.turb.outlet_stage.ratioP[0]))
        )
        _log.info(
            "Outlet stage inlet P={}".format(
                pyo.value(fs.turb.outlet_stage.inlet.pressure[0])
            )
        )
        _log.info(
            "Outlet stage outlet P={}".format(
                pyo.value(fs.turb.outlet_stage.outlet.pressure[0])
            )
        )
        _log.info(
            "bfp turbine out T={}".format(
                pyo.value(fs.bfp_turb.control_volume.properties_out[0].temperature)
            )
        )
        _log.info(
            "bfp turbine out P={}".format(
                pyo.value(fs.bfp_turb.control_volume.properties_out[0].pressure)
            )
        )
        _log.info(
            "bfp turbine out H={}".format(
                pyo.value(fs.bfp_turb.control_volume.properties_out[0].enth_mol)
            )
        )
        _log.info(
            "bfp turbine out flow={}".format(
                pyo.value(fs.bfp_turb.control_volume.properties_out[0].flow_mol)
            )
        )
        _log.info(
            "ip_split9 fraction to DA={}".format(
                pyo.value(fs.turb.ip_split[9].split_fraction[0, "outlet_2"])
            )
        )
        _log.info(
            "ip_split9 fraction to bfpt={}".format(
                pyo.value(fs.turb.ip_split[9].split_fraction[0, "outlet_3"])
            )
        )
        _log.info("booster outlet pres={}".format(fs.booster.outlet.pressure[0].value))
        _log.info(
            "booster inlet temp={}".format(
                pyo.value(fs.booster.control_volume.properties_in[0].temperature)
            )
        )
        _log.info(
            "booster inlet sat temp={}".format(
                pyo.value(fs.booster.control_volume.properties_in[0].temperature_sat)
            )
        )
        # check pressure at mixer inlets
        _log.info(
            "hotwell main inlet pressure={}".format(
                fs.condenser_hotwell.main_condensate.pressure[0].value
            )
        )
        _log.info(
            "hotwell aux inlet pressure={}".format(
                fs.condenser_hotwell.aux_condensate.pressure[0].value
            )
        )
        _log.info(
            "hotwell makeup inlet pressure={}".format(
                fs.condenser_hotwell.makeup.pressure[0].value
            )
        )
        _log.info(
            "drain return feedwater inlet pressure={}".format(
                fs.fwh1_drain_return.feedwater.pressure[0].value
            )
        )
        _log.info(
            "drain return drain inlet pressure={}".format(
                fs.fwh1_drain_return.fwh1_drain.pressure[0].value
            )
        )
        _log.info(
            "DA feedwater inlet pressure={}".format(
                fs.fwh4_deair.feedwater.pressure[0].value
            )
        )
        _log.info(
            "DA drain inlet pressure={}".format(fs.fwh4_deair.drain.pressure[0].value)
        )
        _log.info(
            "DA steam inlet pressure={}".format(fs.fwh4_deair.steam.pressure[0].value)
        )
        _log.info(
            "DA outlet pressure={}".format(fs.fwh4_deair.outlet.pressure[0].value)
        )
        _log.info(
            "fwh2 drain pressure before valve={}".format(
                fs.fwh2_valve.inlet.pressure[0].value
            )
        )
        _log.info(
            "fwh2 drain pressure after valve={}".format(
                fs.fwh2_valve.outlet.pressure[0].value
            )
        )
        _log.info(
            "fwh1 shell pressure ={}".format(fs.fwh1.drain_mix.steam.pressure[0].value)
        )
        _log.info(
            "fwh3 drain pressure before valve={}".format(
                fs.fwh3_valve.inlet.pressure[0].value
            )
        )
        _log.info(
            "fwh3 drain pressure after valve={}".format(
                fs.fwh3_valve.outlet.pressure[0].value
            )
        )
        _log.info(
            "fwh2 shell pressure ={}".format(fs.fwh2.drain_mix.steam.pressure[0].value)
        )
        _log.info(
            "fwh5 drain pressure before valve={}".format(
                fs.fwh5_valve.inlet.pressure[0].value
            )
        )
        _log.info(
            "fwh5 drain pressure after valve={}".format(
                fs.fwh5_valve.outlet.pressure[0].value
            )
        )
        _log.info(
            "fwh6 drain pressure before valve={}".format(
                fs.fwh6_valve.inlet.pressure[0].value
            )
        )
        _log.info(
            "fwh6 drain pressure after valve={}".format(
                fs.fwh6_valve.outlet.pressure[0].value
            )
        )
        _log.info(
            "fwh5 shell pressure ={}".format(fs.fwh5.drain_mix.steam.pressure[0].value)
        )
        _log.info("Cv fwh2={}".format(fs.fwh2_valve.Cv.value))
        _log.info("Cv fwh3={}".format(fs.fwh3_valve.Cv.value))
        _log.info("Cv fwh5={}".format(fs.fwh5_valve.Cv.value))
        _log.info("Cv fwh6={}".format(fs.fwh6_valve.Cv.value))
        _log.info("Cv cond valve={}".format(fs.cond_valve.Cv.value))
        _log.info("Cv makeup valve={}".format(fs.makeup_valve.Cv.value))
        _log.info("Cv spray valve={}".format(fs.spray_valve.Cv.value))
        _log.info("Cv bfp turb valve={}".format(fs.bfp_turb_valve.Cv.value))
        _log.info("Cv throttle valve={}".format(fs.turb.throttle_valve[1].Cv.value))
        _log.info("valve opening fwh2={}".format(fs.fwh2_valve.valve_opening[0].value))
        _log.info("valve opening fwh3={}".format(fs.fwh3_valve.valve_opening[0].value))
        _log.info("valve opening fwh5={}".format(fs.fwh5_valve.valve_opening[0].value))
        _log.info("valve opening fwh6={}".format(fs.fwh6_valve.valve_opening[0].value))
        _log.info("valve opening cond={}".format(fs.cond_valve.valve_opening[0].value))
        _log.info(
            "valve opening makeup={}".format(fs.makeup_valve.valve_opening[0].value)
        )
        _log.info(
            "valve opening spray={}".format(fs.spray_valve.valve_opening[0].value)
        )
        _log.info(
            "valve opening bfp turb={}".format(fs.bfp_turb_valve.valve_opening[0].value)
        )
        _log.info(
            "valve opening throttle={}".format(
                fs.turb.throttle_valve[1].valve_opening[0].value
            )
        )
        _log.info("fwh2 level={}".format(fs.fwh2.condense.level[0].value))
        _log.info("fwh3 level={}".format(fs.fwh3.condense.level[0].value))
        _log.info("fwh5 level={}".format(fs.fwh5.condense.level[0].value))
        _log.info("fwh6 level={}".format(fs.fwh6.condense.level[0].value))
        _log.info("hotwell level={}".format(fs.hotwell_tank.tank_level[0].value))
        _log.info("da level={}".format(fs.da_tank.tank_level[0].value))
        _log.info("cond_pump deltaP={}".format(fs.cond_pump.deltaP[0].value))
        _log.info(
            "cond_pump work={}".format(pyo.value(fs.cond_pump.work_mechanical[0]))
        )
        _log.info("cond_valve deltaP={}".format(fs.cond_valve.deltaP[0].value))
        _log.info("makeup flow={}".format(fs.makeup_valve.outlet.flow_mol[0].value))
        _log.info("spray flow={}".format(fs.spray_valve.outlet.flow_mol[0].value))
        _log.info(
            "spray split fraction={}".format(
                fs.split_attemp.split_fraction[0, "Spray"].value
            )
        )
        _log.info(
            "outlet stage inlet temp={}".format(
                pyo.value(
                    fs.turb.outlet_stage.control_volume.properties_in[0].temperature
                )
            )
        )
        _log.info(
            "outlet stage inlet pres={}".format(
                pyo.value(fs.turb.outlet_stage.control_volume.properties_in[0].pressure)
            )
        )
        _log.info(
            "outlet stage inlet flow_mol={}".format(
                pyo.value(fs.turb.outlet_stage.control_volume.properties_in[0].flow_mol)
            )
        )
        _log.info(
            "outlet stage outlet flow_vol={}".format(
                pyo.value(
                    fs.turb.outlet_stage.control_volume.properties_out[0].flow_vol
                )
            )
        )
        _log.info(
            "outlet stage inlet vapor frac={}".format(
                pyo.value(
                    fs.turb.outlet_stage.control_volume.properties_in[0].vapor_frac
                )
            )
        )
        _log.info(
            "outlet state outlet vapor frac={}".format(
                pyo.value(
                    fs.turb.outlet_stage.control_volume.properties_out[0].vapor_frac
                )
            )
        )
        _log.info(
            "outlet stage efficiency={}".format(
                pyo.value(fs.turb.outlet_stage.efficiency_isentropic[0])
            )
        )
        _log.info(
            "inlet stage efficiency={}".format(
                pyo.value(fs.turb.inlet_stage[1].efficiency_isentropic[0])
            )
        )
        _log.info(
            "HP stage 1 efficiency={}".format(
                pyo.value(fs.turb.hp_stages[1].efficiency_isentropic[0])
            )
        )
        _log.info(
            "IP stage 1 efficiency={}".format(
                pyo.value(fs.turb.ip_stages[1].efficiency_isentropic[0])
            )
        )
        _log.info(
            "LP stage 1 efficiency={}".format(
                pyo.value(fs.turb.lp_stages[1].efficiency_isentropic[0])
            )
        )
        _log.info("Tel outlet stage={}".format(pyo.value(fs.turb.outlet_stage.tel[0])))
        _log.info(
            "delta_h ={}".format(
                pyo.value(fs.turb.outlet_stage.delta_enth_isentropic[0])
            )
        )
        _log.info(
            "outlet stage power ={}".format(
                -1e-6 * pyo.value(fs.turb.outlet_stage.power_shaft[0])
            )
        )
        _log.info(
            "inlet stage pressure ratio={}".format(
                pyo.value(fs.turb.inlet_stage[1].ratioP[0])
            )
        )
        _log.info(
            "bfp turb valve inlet pressure={}".format(
                fs.bfp_turb_valve.inlet.pressure[0].value
            )
        )
        _log.info(
            "bfp turb inlet pressure={}".format(fs.bfp_turb.inlet.pressure[0].value)
        )
        _log.info(
            "bfp turb outlet flow={}".format(fs.bfp_turb.outlet.flow_mol[0].value)
        )
        _log.info(
            "bfp turb outlet enth={}".format(fs.bfp_turb.outlet.enth_mol[0].value)
        )
        _log.info(
            "bfp turb outlet pres={}".format(fs.bfp_turb.outlet.pressure[0].value)
        )
        _log.info(
            "bfp turb work={}".format(-1e-6 * fs.bfp_turb.control_volume.work[0].value)
        )
        _log.info(
            "bfp turb efficiency={}".format(fs.bfp_turb.efficiency_isentropic[0].value)
        )
        _log.info(
            "bfp turb os outlet flow={}".format(fs.bfp_turb_os.outlet.flow_mol[0].value)
        )
        _log.info(
            "bfp turb os outlet enth={}".format(fs.bfp_turb_os.outlet.enth_mol[0].value)
        )
        _log.info(
            "bfp turb os outlet pres={}".format(fs.bfp_turb_os.outlet.pressure[0].value)
        )
        _log.info(
            "bfp turb os outlet vapor frac={}".format(
                pyo.value(fs.bfp_turb_os.control_volume.properties_out[0].vapor_frac)
            )
        )
        _log.info(
            "bfp turb os work={}".format(
                -1e-6 * fs.bfp_turb_os.control_volume.work[0].value
            )
        )
        _log.info(
            "bfp turb os efficiency={}".format(
                fs.bfp_turb_os.efficiency_isentropic[0].value
            )
        )
        _log.info("Tel turb os ={}".format(pyo.value(fs.bfp_turb_os.tel[0])))
        _log.info("bfp turb os flow coef={}".format(fs.bfp_turb_os.flow_coeff.value))
        _log.info(
            "bfp turb os out vol flow={}".format(
                pyo.value(fs.bfp_turb_os.control_volume.properties_out[0].flow_vol)
            )
        )

        # Since the constraint to calculate flow rate based on valve opening
        # and pressure drop is based on the square of flow rate and opening,
        # negative a valve opening is mathmatically valid.  Set valve openings
        # to physically valid positive numbers
        valve_open = fs.fwh2_valve.valve_opening[0].value
        if valve_open < 0:
            fs.fwh2_valve.valve_opening.value = -valve_open
        valve_open = fs.fwh3_valve.valve_opening[0].value
        if valve_open < 0:
            fs.fwh3_valve.valve_opening[:].value = -valve_open
        valve_open = fs.fwh5_valve.valve_opening[0].value
        if valve_open < 0:
            fs.fwh5_valve.valve_opening[:].value = -valve_open
        valve_open = fs.fwh6_valve.valve_opening[0].value
        if valve_open < 0:
            fs.fwh6_valve.valve_opening[:].value = -valve_open

        _log.info("Adding heat transfer coefficent correations...")
        _add_heat_transfer_correlation(fs)

    else:
        _add_heat_transfer_correlation(fs)

        fs.fwh2.condense.level.unfix()
        fs.fwh3.condense.level.unfix()
        fs.fwh5.condense.level.unfix()
        fs.fwh6.condense.level.unfix()
        fs.hotwell_tank.tank_level.unfix()
        fs.da_tank.tank_level.unfix()

        fs.fwh2_valve.valve_opening.unfix()
        fs.fwh3_valve.valve_opening.unfix()
        fs.fwh5_valve.valve_opening.unfix()
        fs.fwh6_valve.valve_opening.unfix()
        fs.makeup_valve.valve_opening.unfix()
        fs.cond_valve.valve_opening.unfix()
        fs.spray_valve.valve_opening.unfix()

    return m


def _add_u_eq(blk, uex=0.8):
    """Add heat transfer coefficent adjustment for feed water flow rate.
    This is based on knowing the heat transfer coefficent at a particular flow
    and assuming the heat transfer coefficent is porportial to feed water
    flow rate raised to certain power (typically 0.8)

    Args:
        blk: Heat exchanger block to add correlation to
        uex: Correlation parameter value (defalut 0.8)

    Returns:
        None
    """
    ti = blk.flowsheet().time
    blk.U0 = pyo.Var(ti)
    blk.f0 = pyo.Var(ti)
    blk.uex = pyo.Var(ti, initialize=uex)
    for t in ti:
        blk.U0[t].value = blk.overall_heat_transfer_coefficient[t].value
        blk.f0[t].value = blk.tube.properties_in[t].flow_mol.value
    blk.overall_heat_transfer_coefficient.unfix()
    blk.U0.fix()
    blk.uex.fix()
    blk.f0.fix()

    @blk.Constraint(ti)
    def U_eq(b, t):
        return (
            b.overall_heat_transfer_coefficient[t]
            == b.U0[t] * (b.tube.properties_in[t].flow_mol / b.f0[t]) ** b.uex[t]
        )


def set_scaling_factors(m):
    """Set scaling factors for variables and expressions. These are used for
    variable scaling and used by the framework to scale constraints.

    Args:
        m: plant model to set scaling factors for.

    Returns:
        None
    """

    # Set steam cycle scale factors
    fs = m.fs_main.fs_stc

    iscale.set_scaling_factor(fs.condenser.hot_side.heat, 1e-9)
    iscale.set_scaling_factor(fs.condenser.cold_side.heat, 1e-9)

    iscale.set_scaling_factor(fs.aux_condenser.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.aux_condenser.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.hotwell_tank.control_volume.energy_holdup, 1e-9)
    iscale.set_scaling_factor(fs.hotwell_tank.control_volume.material_holdup, 1e-6)

    iscale.set_scaling_factor(fs.fwh1.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh1.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh1.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh1.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh1.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh1.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.fwh2.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh2.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh2.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh2.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh2.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh2.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.fwh3.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh3.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh3.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh3.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh3.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh3.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.da_tank.control_volume.energy_holdup, 1e-10)
    iscale.set_scaling_factor(fs.da_tank.control_volume.material_holdup, 1e-6)

    iscale.set_scaling_factor(fs.fwh5.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh5.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh5.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh5.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh5.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh5.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.fwh6.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh6.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh6.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh6.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh6.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh6.condense.cold_side.heat, 1e-7)

    # scaling factor for control valves
    for t in m.fs_main.time:
        iscale.set_scaling_factor(
            fs.spray_valve.control_volume.properties_in[t].flow_mol, 0.05
        )
        iscale.set_scaling_factor(
            fs.makeup_valve.control_volume.properties_in[t].flow_mol, 0.01
        )
        iscale.set_scaling_factor(
            fs.fwh2_valve.control_volume.properties_in[t].flow_mol, 0.001
        )
        iscale.set_scaling_factor(
            fs.fwh3_valve.control_volume.properties_in[t].flow_mol, 0.001
        )
        iscale.set_scaling_factor(
            fs.fwh5_valve.control_volume.properties_in[t].flow_mol, 0.001
        )
        iscale.set_scaling_factor(
            fs.fwh6_valve.control_volume.properties_in[t].flow_mol, 0.001
        )
        iscale.set_scaling_factor(
            fs.bfp_turb_valve.control_volume.properties_in[t].flow_mol, 0.01
        )

    # Calculate calculated scaling factors
    iscale.calculate_scaling_factors(m)


def main_steady_state():
    m_ss = get_model(dynamic=False)
    solver = get_solver()

    dof = degrees_of_freedom(m_ss)
    print("dof of full model", dof)
    # solving dynamic model at steady-state
    print("solving dynamic model at steady-state...")
    solver.solve(m_ss, tee=True)
    return m_ss


def main_dynamic():
    m_ss = get_model(dynamic=False)
    m_dyn = get_model(dynamic=True)
    copy_non_time_indexed_values(
        m_dyn.fs_main, m_ss.fs_main, copy_fixed=True, outlvl=idaeslog.ERROR
    )
    for t in m_dyn.fs_main.time:
        copy_values_at_time(
            m_dyn.fs_main, m_ss.fs_main, t, 0.0, copy_fixed=True, outlvl=idaeslog.ERROR
        )
    t0 = 0
    # estimate integral error for the PI controller
    m_dyn.fs_main.fs_stc.fwh2_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.fwh2_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.fwh3_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.fwh3_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.fwh5_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.fwh5_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.fwh6_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.fwh6_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.makeup_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.makeup_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.da_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.cond_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.spray_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.spray_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error[:].value = pyo.value(
        m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error_ref[t0]
    )
    m_dyn.fs_main.fs_stc.makeup_ctrl.integral_of_error[:].value = pyo.value(
        m_dyn.fs_main.fs_stc.makeup_ctrl.integral_of_error_ref[t0]
    )

    m_dyn.fs_main.fs_stc.fwh2.condense.level[0].fix()
    m_dyn.fs_main.fs_stc.fwh3.condense.level[0].fix()
    m_dyn.fs_main.fs_stc.fwh5.condense.level[0].fix()
    m_dyn.fs_main.fs_stc.fwh6.condense.level[0].fix()
    m_dyn.fs_main.fs_stc.hotwell_tank.tank_level[0].fix()
    m_dyn.fs_main.fs_stc.da_tank.tank_level[0].fix()

    m_dyn.fs_main.fs_stc.spray_valve.valve_opening[0].fix()

    solver = get_solver()

    outlvl = idaeslog.DEBUG
    _log = idaeslog.getLogger(m_dyn.name, outlvl, tag="flowsheet")
    dof = degrees_of_freedom(m_dyn.fs_main)
    _log.debug("dof full model={}".format(dof))
    # solving dynamic model at steady-state
    print("solving dynamic model at steady-state...")
    solver.solve(m_dyn.fs_main, tee=True)

    _log.debug(
        "main steam enth={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.inlet_split.inlet.enth_mol[0])
        )
    )
    _log.debug(
        "main steam flow_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0])
        )
    )
    _log.debug(
        "main steam pressure={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.inlet_split.inlet.pressure[0])
        )
    )
    _log.debug(
        "throttle_valve Cv={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.throttle_valve[1].Cv)
        )
    )
    _log.debug(
        "throttle_valve opening={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.throttle_valve[1].valve_opening[0])
        )
    )
    _log.debug(
        "HP stage 1 inlet enth_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.hp_stages[1].inlet.enth_mol[0])
        )
    )
    _log.debug(
        "HP stage 1 inlet flow_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.hp_stages[1].inlet.flow_mol[0])
        )
    )
    _log.debug(
        "HP stage 1 inlet pressure={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.hp_stages[1].inlet.pressure[0])
        )
    )
    _log.debug(
        "IP stage 1 inlet enth_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.ip_stages[1].inlet.enth_mol[0])
        )
    )
    _log.debug(
        "IP stage 1 inlet flow_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.ip_stages[1].inlet.flow_mol[0])
        )
    )
    _log.debug(
        "IP stage 1 inlet pressure={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.ip_stages[1].inlet.pressure[0])
        )
    )
    _log.debug(
        "Outlet stage enth_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.outlet_stage.outlet.enth_mol[0])
        )
    )
    _log.debug(
        "Outlet stage flow_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.outlet_stage.outlet.flow_mol[0])
        )
    )
    _log.debug(
        "Outlet stage pressure={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.outlet_stage.outlet.pressure[0])
        )
    )
    _log.debug(
        "Power output of main turbine={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.power_output[0])
        )
    )
    _log.debug(
        "Power output of bfp turbine={}".format(
            pyo.value(
                m_dyn.fs_main.fs_stc.bfp_turb.control_volume.work[0]
                + m_dyn.fs_main.fs_stc.bfp_turb_os.control_volume.work[0]
            )
            * (-1e-6)
        )
    )
    _log.debug(
        "FWH6 outlet enth_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.fwh6.desuperheat.outlet_2.enth_mol[0])
        )
    )
    _log.debug(
        "FWH6 outlet flow_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.fwh6.desuperheat.outlet_2.flow_mol[0])
        )
    )
    _log.debug(
        "FWH6 outlet pressure={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.fwh6.desuperheat.outlet_2.pressure[0])
        )
    )
    _log.debug(
        "water makeup flow={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.condenser_hotwell.makeup.flow_mol[0])
        )
    )
    _log.debug(
        "spray flow={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.spray_valve.outlet.flow_mol[0])
        )
    )
    _log.debug("Cv fwh2={}".format(m_dyn.fs_main.fs_stc.fwh2_valve.Cv.value))
    _log.debug("Cv fwh3={}".format(m_dyn.fs_main.fs_stc.fwh3_valve.Cv.value))
    _log.debug("Cv fwh5={}".format(m_dyn.fs_main.fs_stc.fwh5_valve.Cv.value))
    _log.debug("Cv fwh6={}".format(m_dyn.fs_main.fs_stc.fwh6_valve.Cv.value))
    _log.debug("Cv cond_valve={}".format(m_dyn.fs_main.fs_stc.cond_valve.Cv.value))
    _log.debug("Cv makeup_valve={}".format(m_dyn.fs_main.fs_stc.makeup_valve.Cv.value))
    _log.debug("Cv spray_valve={}".format(m_dyn.fs_main.fs_stc.spray_valve.Cv.value))
    _log.debug(
        "valve opening fwh2={}".format(
            m_dyn.fs_main.fs_stc.fwh2_valve.valve_opening[0].value
        )
    )
    _log.debug(
        "valve opening fwh3={}".format(
            m_dyn.fs_main.fs_stc.fwh3_valve.valve_opening[0].value
        )
    )
    _log.debug(
        "valve opening fwh5={}".format(
            m_dyn.fs_main.fs_stc.fwh5_valve.valve_opening[0].value
        )
    )
    _log.debug(
        "valve opening fwh6={}".format(
            m_dyn.fs_main.fs_stc.fwh6_valve.valve_opening[0].value
        )
    )
    _log.debug(
        "valve opening cond_valve={}".format(
            m_dyn.fs_main.fs_stc.cond_valve.valve_opening[0].value
        )
    )
    _log.debug(
        "valve opening makeup_valve={}".format(
            m_dyn.fs_main.fs_stc.makeup_valve.valve_opening[0].value
        )
    )
    _log.debug(
        "valve opening spray={}".format(
            m_dyn.fs_main.fs_stc.spray_valve.valve_opening[0].value
        )
    )
    _log.debug(
        "fwh2 level={}".format(m_dyn.fs_main.fs_stc.fwh2.condense.level[0].value)
    )
    _log.debug(
        "fwh3 level={}".format(m_dyn.fs_main.fs_stc.fwh3.condense.level[0].value)
    )
    _log.debug(
        "fwh5 level={}".format(m_dyn.fs_main.fs_stc.fwh5.condense.level[0].value)
    )
    _log.debug(
        "fwh6 level={}".format(m_dyn.fs_main.fs_stc.fwh6.condense.level[0].value)
    )
    _log.debug(
        "hotwell tank level={}".format(
            m_dyn.fs_main.fs_stc.hotwell_tank.tank_level[0].value
        )
    )
    _log.debug(
        "da tank level={}".format(m_dyn.fs_main.fs_stc.da_tank.tank_level[0].value)
    )
    _log.debug(
        "makeup flow={}".format(
            m_dyn.fs_main.fs_stc.makeup_valve.outlet.flow_mol[0].value
        )
    )
    _log.debug(
        "spray flow={}".format(
            m_dyn.fs_main.fs_stc.spray_valve.outlet.flow_mol[0].value
        )
    )
    _log.debug(
        "split_attemp fraction={}".format(
            m_dyn.fs_main.fs_stc.split_attemp.split_fraction[0, "Spray"].value
        )
    )

    # impose step change for dynamic model
    for t in m_dyn.fs_main.time:
        if t >= 30:
            m_dyn.fs_main.fs_stc.turb.ip_stages[1].inlet.pressure[t].fix(
                m_dyn.fs_main.fs_stc.turb.ip_stages[1].inlet.pressure[0].value * 1.05
            )
        else:
            m_dyn.fs_main.fs_stc.turb.ip_stages[1].inlet.pressure[t].fix(
                m_dyn.fs_main.fs_stc.turb.ip_stages[1].inlet.pressure[0].value
            )
    dof = degrees_of_freedom(m_dyn.fs_main)
    _log.debug("dof of full model={}".format(dof))
    # solving dynamic model
    _log.debug("solving dynamic model...")
    solver.solve(m_dyn.fs_main, tee=True)
    _log.debug(
        "Power output of main turbine={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.power_output[60])
        )
    )
    _log.debug(
        "Power output of bfp turbine={}".format(
            pyo.value(
                m_dyn.fs_main.fs_stc.bfp_turb.control_volume.work[60]
                + m_dyn.fs_main.fs_stc.bfp_turb_os.control_volume.work[60]
            )
            * (-1e-6)
        )
    )
    _log.debug(
        "main steam enth={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.inlet_split.inlet.enth_mol[60])
        )
    )
    _log.debug(
        "main steam flow_mol={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[60])
        )
    )
    _log.debug(
        "main steam pressure={}".format(
            pyo.value(m_dyn.fs_main.fs_stc.turb.inlet_split.inlet.pressure[60])
        )
    )
    _log.debug(
        "fw pressure={}".format(pyo.value(m_dyn.fs_main.fs_stc.bfp.outlet.pressure[60]))
    )

    # prepare for plots
    time = []
    pres_ip = []
    power_gross = []
    da_level = []
    cond_valve_open = []
    for t in m_dyn.fs_main.time:
        time.append(t)
        pres_ip.append(m_dyn.fs_main.fs_stc.turb.ip_stages[1].inlet.pressure[t].value)
        power_gross.append(pyo.value(m_dyn.fs_main.fs_stc.power_output[t]))
        da_level.append(m_dyn.fs_main.fs_stc.da_tank.tank_level[t].value)
        cond_valve_open.append(m_dyn.fs_main.fs_stc.cond_valve.valve_opening[t].value)

    # plot the dynamic solution
    plt.figure(1)
    plt.plot(time, pres_ip)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("IP Inlet Pressure [Pa]")
    plt.show(block=False)

    plt.figure(2)
    plt.plot(time, power_gross)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Gross Power [MW]")
    plt.show(block=False)

    plt.figure(3)
    plt.plot(time, da_level)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("DA Tank Level [m]")
    plt.show(block=False)

    plt.figure(4)
    plt.plot(time, cond_valve_open)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Condensate Valve Opening")
    plt.show(block=True)

    return m_dyn


def get_model(dynamic=True):
    m = pyo.ConcreteModel()
    m.dynamic = dynamic
    if m.dynamic:
        m.fs_main = FlowsheetBlock(
            dynamic=True, time_set=[0, 60], time_units=pyo.units.s
        )
    else:
        m.fs_main = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs_main.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs_main.fs_stc = FlowsheetBlock(time_units=pyo.units.s)
    m = add_unit_models(m)
    if m.dynamic:
        m.discretizer = pyo.TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=2, wrt=m.fs_main.time, scheme="BACKWARD")
    m = set_arcs_and_constraints(m)
    m = set_inputs(m)
    set_scaling_factors(m)
    m = initialize(m)
    return m


if __name__ == "__main__":
    # This method builds and runs a steam cycle flowsheet, the flowsheet
    # includes the Turbine train, Condenser, Feed Water Heaters and Pumps,
    # fixed inlets are steam flowrates from the boiler (Main Steam and
    # Hot Reheat) and makeup of water.
    m = main_steady_state()

    # This method builds and runs a dynamic simulation of a subcritical steam
    # cycle flowsheet.
    # m = main_dynamic()
