##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

import re
import os

import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.common.fileutils import this_file_dir

import idaes
from idaes.core.solvers import use_idaes_solver_configuration_defaults
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.initialization as iinit
import idaes.power_generation.unit_models.helm as helm
import idaes.generic_models.unit_models as gum
from idaes.generic_models.properties import iapws95
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.solvers import use_idaes_solver_configuration_defaults
import idaes.core.util.tables as tables
import idaes.core.util as iutil


def tag_model(m):
    """
    Create to dictionaries. One is called tags with tags for keys and expressions
    for values. The second is called tag_format with tags for keys and a
    formatting string.  The formating string controls the number format and can
    be used to add additional text, to report units for example. The two
    dictionaries are returned, and also attached to the model.

    Args:
        m: (ConcreteModel) model to tag

    Returns:
        (dict) tags, (dict) tag_format
    """
    tags = {}  # dict of with tag keys and expressions for their values
    tag_format = {}  # format string for the tags

    def new_tag(name, expr, format):
        # funcion to keep it more compact
        tags[name] = expr
        tag_format[name] = format

    # Create a dict with Arc name keys and state block values
    stream_states = tables.stream_states_dict(
        tables.arcs_to_stream_dict(
            m.fs,
            additional={  # streams that are half in HRSG, and may not have arcs
                "t01": m.fs.steam_turbine.inlet_split.inlet,
                "t02": m.fs.steam_turbine.hp_stages[7].outlet,
                "t03": m.fs.steam_turbine.ip_stages[1].inlet,
                "t05": m.fs.steam_turbine_lp_mix.hrsg,
                "t11": m.fs.return_mix.outlet,
                "t12": m.fs.hotwell.makeup,
                "t13": m.fs.return_mix.reboiler,
                "t14": m.fs.return_mix.dryer,
                "t15": m.fs.return_mix.reclaimer,
                "t17": m.fs.steam_turbine_lp_split.reboiler,
                "cw01": m.fs.main_condenser.tube_inlet,
                "cw02": m.fs.main_condenser.tube_outlet,
            },
            descend_into=False,
        )
    )
    for i, s in stream_states.items():  # create the tags for stream quantities
        # Just catch the streams in the steam turbine, since there will probably
        # be a bunch of other stuff in the flowsheet
        if re.match(r"^t[0-9][0-9]$", i) or re.match(r"^cw[0-9][0-9]$", i):
            new_tag(f"{i}_F", expr=s.flow_mol / 1000, format="{:.4f} kmol/s")
            new_tag(f"{i}_P", expr=s.pressure / 1000, format="{:.3f} kPa")
            new_tag(f"{i}_T", expr=s.temperature, format="{:.2f} K")
            new_tag(f"{i}_H", expr=s.enth_mol / 1000, format="{:.2f} kJ/mol")
            new_tag(f"{i}_X", expr=s.vapor_frac * 100, format="{:.2f} %")

    if hasattr(m, "tags"):
        m.tags.update(tags)
        m.tag_format.update(tag_format)
    else:
        m.tags = tags
        m.tag_format = tag_format
    return tags, tag_format


def write_pfd_results(filename, tags, tag_format):
    """
    Write simulation results in a template PFD in svg format and save as
    filename.

    Args:
        filename: (str) file namd for output
        tags: (dict) tag keys and expression values
        tag_format: (dict) tag keys and format string values

    Returns:
        None
    """
    template = os.path.join(this_file_dir(), "steam_turbine_template.svg")
    with open(template, "r") as f:
        iutil.svg_tag(svg=f, tags=tags, outfile=filename, tag_format=tag_format)


def get_model(m=None, init=True):
    if m is None:
        m = pyo.ConcreteModel()
    if not hasattr(m, "fs"):
        m.fs = FlowsheetBlock(default={"dynamic": False})
    if not hasattr(m.fs, "prop_water"):
        m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    solver = pyo.SolverFactory("ipopt")
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    m.fs.steam_turbine = helm.HelmTurbineMultistage(
        default={
            "property_package": m.fs.prop_water,
            "num_parallel_inlet_stages": 4,
            "num_hp": 7,  # at full load ave P ratio about 0.8238 with inlet stage
            "num_ip": 10,  # at full load ave P ratio about 0.8264
            "num_lp": 11,  # at full load ave P ratio about 0.7194 with outlet stage
            "hp_disconnect": [7],  # disconected for reheater
            "ip_disconnect": [10],
        }
    )  # disconnected for HRSG LP steam mix
    m.fs.steam_turbine_lp_mix = helm.HelmMixer(
        default={
            "property_package": m.fs.prop_water,
            "momentum_mixing_type": helm.MomentumMixingType.none,
            "inlet_list": ["turbine", "hrsg"],
        }
    )
    # The mixer for LP steam from Turbine and HRSG is assumed to be at turbine P
    @m.fs.steam_turbine_lp_mix.Constraint(m.fs.time)
    def lp_mixer_pressure_constraint(b, t):
        return 1e-6 * b.turbine_state[t].pressure == 1e-6 * b.mixed_state[t].pressure

    # Steam splitter for capture
    m.fs.steam_turbine_lp_split = helm.HelmSplitter(
        default={
            "property_package": m.fs.prop_water,
            "outlet_list": ["turbine", "reboiler"],
        }
    )
    m.fs.dummy_reheat = gum.Heater(default={"property_package": m.fs.prop_water})
    m.fs.dummy_reheat.temperature_out = pyo.Var(m.fs.time, initialize=850)

    @m.fs.dummy_reheat.Constraint(m.fs.time)
    def temperature_eqn(b, t):
        return b.temperature_out[t] == b.control_volume.properties_out[t].temperature

    m.fs.dummy_reheat.temperature_out.fix(858)
    m.fs.main_condenser = helm.HelmNtuCondenser(
        default={
            "shell": {
                "has_pressure_change": False,
                "property_package": m.fs.prop_water,
            },
            "tube": {"has_pressure_change": False, "property_package": m.fs.prop_water},
        }
    )
    m.fs.hotwell = helm.HelmMixer(
        default={
            "momentum_mixing_type": helm.MomentumMixingType.none,
            "inlet_list": ["condensate", "makeup"],
            "property_package": m.fs.prop_water,
        }
    )
    # The hotwell is assumed to be at the same pressure as the condenser.
    @m.fs.hotwell.Constraint(m.fs.time)
    def hotwell_pressure_constraint(b, t):
        return 1e-6 * b.condensate_state[t].pressure == 1e-6 * b.mixed_state[t].pressure

    m.fs.cond_pump = helm.HelmIsentropicCompressor(
        default={"property_package": m.fs.prop_water}
    )
    # condensate return mixer
    m.fs.return_mix = helm.HelmMixer(
        default={
            "property_package": m.fs.prop_water,
            "momentum_mixing_type": helm.MomentumMixingType.none,
            "inlet_list": ["pump", "reboiler", "dryer", "reclaimer"],
        }
    )

    @m.fs.return_mix.Constraint(m.fs.time)
    def return_mixer_pressure_constraint(b, t):
        return 1e-6 * b.pump_state[t].pressure == 1e-6 * b.mixed_state[t].pressure

    # A few more variables and constraints
    m.fs.hp_steam_temperature = pyo.Var(m.fs.time, initialize=850)
    m.fs.hot_reheat_temperature = pyo.Var(m.fs.time, initialize=850)

    @m.fs.Constraint(m.fs.time)
    def main_steam_temperature_eqn(b, t):
        return (
            b.hp_steam_temperature[t]
            == b.steam_turbine.inlet_split.mixed_state[t].temperature
        )

    @m.fs.Constraint(m.fs.time)
    def reheat_steam_temperature_eqn(b, t):
        return (
            b.hot_reheat_temperature[t]
            == b.steam_turbine.ip_stages[1].control_volume.properties_out[t].temperature
        )

    # arcs
    m.fs.t02_dummy = Arc(
        source=m.fs.steam_turbine.hp_stages[7].outlet,
        destination=m.fs.dummy_reheat.inlet,
    )
    m.fs.t03_dummy = Arc(
        source=m.fs.dummy_reheat.outlet,
        destination=m.fs.steam_turbine.ip_stages[1].inlet,
    )
    m.fs.t04 = Arc(
        source=m.fs.steam_turbine.ip_stages[10].outlet,
        destination=m.fs.steam_turbine_lp_mix.turbine,
    )
    m.fs.t06 = Arc(
        source=m.fs.steam_turbine_lp_split.turbine,
        destination=m.fs.steam_turbine.lp_stages[1].inlet,
    )
    m.fs.t16 = Arc(
        source=m.fs.steam_turbine_lp_mix.outlet,
        destination=m.fs.steam_turbine_lp_split.inlet,
    )
    m.fs.t07 = Arc(
        source=m.fs.steam_turbine.outlet_stage.outlet,
        destination=m.fs.main_condenser.shell_inlet,
    )
    m.fs.t08 = Arc(
        source=m.fs.main_condenser.shell_outlet, destination=m.fs.hotwell.condensate
    )
    m.fs.t09 = Arc(source=m.fs.hotwell.outlet, destination=m.fs.cond_pump.inlet)
    m.fs.t10 = Arc(source=m.fs.cond_pump.outlet, destination=m.fs.return_mix.pump)
    expand_arcs.apply_to(m)
    # set some inputs and get ready to initialize
    # Set the inlet of the turbine
    p = 16.5e6
    hin = pyo.value(iapws95.htpx(T=858 * pyo.units.K, P=p * pyo.units.Pa))
    m.fs.steam_turbine.inlet_split.inlet.enth_mol[0].fix(hin)
    m.fs.steam_turbine.inlet_split.inlet.flow_mol[0].fix(6.5064e3)
    m.fs.steam_turbine.inlet_split.inlet.pressure[0].fix(p)
    # set conditions for disconnected IP section, will take flow from HP
    p = 3.5e6
    hin = pyo.value(iapws95.htpx(T=858 * pyo.units.K, P=p * pyo.units.Pa))
    m.fs.steam_turbine.ip_stages[1].inlet.flow_mol[0].value = 6.5064e3
    m.fs.steam_turbine.ip_stages[1].inlet.enth_mol[0].value = hin
    m.fs.steam_turbine.ip_stages[1].inlet.pressure[0].value = p
    # set conditions for disconnected LP section, will take flow from HP
    p = 0.6e6
    hin = pyo.value(iapws95.htpx(T=573 * pyo.units.K, P=p * pyo.units.Pa))
    m.fs.steam_turbine.lp_stages[1].inlet.flow_mol[0].value = 9.1e3
    m.fs.steam_turbine.lp_stages[1].inlet.enth_mol[0].value = hin
    m.fs.steam_turbine.lp_stages[1].inlet.pressure[0].value = p
    # use same conditions for lp steam from HRSG, but also specify flow
    hin = pyo.value(iapws95.htpx(T=601 * pyo.units.K, P=p * pyo.units.Pa))
    m.fs.steam_turbine_lp_mix.hrsg.enth_mol.fix(hin)
    m.fs.steam_turbine_lp_mix.hrsg.pressure.fix(p)
    m.fs.steam_turbine_lp_mix.hrsg.flow_mol.fix(3000)

    # m.fs.steam_turbine_lp_split.split_fraction[0, "reboiler"].fix(0.001)
    m.fs.steam_turbine_lp_split.reboiler.flow_mol.fix(4000.0)

    for i, s in m.fs.steam_turbine.hp_stages.items():
        s.ratioP[:] = 0.8238
        s.efficiency_isentropic[:] = 0.89
        iscale.set_scaling_factor(s.control_volume.work, 1e-6)
    for i, s in m.fs.steam_turbine.ip_stages.items():
        s.ratioP[:] = 0.8264
        s.efficiency_isentropic[:] = 0.89
        iscale.set_scaling_factor(s.control_volume.work, 1e-6)
    for i, s in m.fs.steam_turbine.lp_stages.items():
        s.ratioP[:] = 0.75
        s.efficiency_isentropic[:] = 0.89
        iscale.set_scaling_factor(s.control_volume.work, 1e-6)

    m.fs.steam_turbine.outlet_stage.design_exhaust_flow_vol.fix(2300)
    # Unfix the main steam flow for pressure driven flow
    m.fs.steam_turbine.inlet_split.inlet.flow_mol.unfix()
    m.fs.steam_turbine.outlet_stage.control_volume.properties_out[0].pressure.fix(
        0.01e6
    )
    # equal outlet pressure from the throttle valves
    m.fs.steam_turbine.inlet_mix.use_equal_pressure_constraint()
    # setup the parallel inlet stages
    for i, s in m.fs.steam_turbine.inlet_stage.items():
        iscale.set_scaling_factor(s.control_volume.work, 1e-6)
        s.ratioP[0] = 0.82
        m.fs.steam_turbine.throttle_valve[i].Cv.fix()
        m.fs.steam_turbine.throttle_valve[i].valve_opening.fix(0.85)
    iscale.set_scaling_factor(m.fs.steam_turbine.outlet_stage.control_volume.work, 1e-6)

    m.fs.main_condenser.tube_inlet.flow_mol.fix(2e5)
    m.fs.main_condenser.tube_inlet.enth_mol.fix(1900)
    m.fs.main_condenser.tube_inlet.pressure.fix(5e5)
    m.fs.main_condenser.area.fix(5000)
    m.fs.main_condenser.overall_heat_transfer_coefficient.fix(15000)

    m.fs.hotwell.makeup.flow_mol[:].fix(1)
    m.fs.hotwell.makeup.enth_mol.fix(2500)
    m.fs.hotwell.makeup.pressure.fix(101325)
    m.fs.cond_pump.efficiency_isentropic.fix(0.80)
    m.fs.cond_pump.deltaP.fix(2.5e6)

    p = 0.49e6
    hin = pyo.value(iapws95.htpx(T=424 * pyo.units.K, P=p * pyo.units.Pa))
    m.fs.return_mix.reboiler.enth_mol[0].fix(hin)
    m.fs.return_mix.reboiler.flow_mol[0].fix(4000.0)
    m.fs.return_mix.reboiler.pressure[0].fix(p)
    p = 2e6
    hin = pyo.value(iapws95.htpx(T=487 * pyo.units.K, P=p * pyo.units.Pa))
    m.fs.return_mix.reclaimer.enth_mol[0].fix(hin)
    m.fs.return_mix.reclaimer.flow_mol[0].fix(0.036)
    m.fs.return_mix.reclaimer.pressure[0].fix(p)
    p = 1.6e6
    hin = pyo.value(iapws95.htpx(T=476 * pyo.units.K, P=p * pyo.units.Pa))
    m.fs.return_mix.dryer.enth_mol[0].fix(hin)
    m.fs.return_mix.dryer.flow_mol[0].fix(0.0019)
    m.fs.return_mix.dryer.pressure[0].fix(p)

    iscale.set_scaling_factor(m.fs.dummy_reheat.control_volume.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.main_condenser.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.main_condenser.tube.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.cond_pump.control_volume.work, 1e-6)

    # scaling
    iscale.calculate_scaling_factors(m)
    # This initializtion setup will use the inlet stage pressure ratios to
    # calculate flow coefficients for the inlet stages and the steam flow to
    # calculate the flow coefficient for the outlet stage.
    # fix the IP and LP inlets since they are disconnected
    m.fs.steam_turbine.ip_stages[1].inlet.fix()
    m.fs.steam_turbine.lp_stages[1].inlet.fix()
    if init:
        m.fs.steam_turbine.initialize(
            outlvl=idaeslog.INFO,
            copy_disconneted_flow=False,
            calculate_inlet_cf=True,
            calculate_outlet_cf=True,
        )
    m.fs.steam_turbine.ip_stages[1].inlet.unfix()
    m.fs.steam_turbine.lp_stages[1].inlet.unfix()

    iinit.propagate_state(arc=m.fs.t02_dummy)
    if init:
        m.fs.dummy_reheat.initialize(outlvl=idaeslog.INFO)
    iinit.propagate_state(arc=m.fs.t03_dummy)

    iinit.propagate_state(arc=m.fs.t04)
    if init:
        m.fs.steam_turbine_lp_mix.initialize(outlvl=idaeslog.INFO)
    iinit.propagate_state(arc=m.fs.t06)

    iinit.propagate_state(arc=m.fs.t07)
    if init:
        m.fs.main_condenser.initialize(outlvl=idaeslog.INFO, unfix="pressure")

    iinit.propagate_state(arc=m.fs.t08)
    if init:
        m.fs.hotwell.initialize(outlvl=idaeslog.INFO)

    iinit.propagate_state(arc=m.fs.t09)
    if init:
        m.fs.cond_pump.initialize(outlvl=idaeslog.INFO)

    # Unfix the turbine outlet pressure, is determined by the condenser model
    m.fs.steam_turbine.outlet_stage.control_volume.properties_out[0].pressure.unfix()

    if init:
        print(
            f"Condenser pressure {pyo.value(m.fs.main_condenser.shell.properties_out[0].pressure)}"
        )
        print(f"DoF: {degrees_of_freedom(m)}")
        solver.solve(m, tee=True)

        print("Adjust flow coefficient for steam extraction")
        m.fs.steam_turbine.outlet_stage.flow_coeff.fix(0.10)
        solver.solve(m, tee=True)

    tag_model(m)
    return m, solver


if __name__ == "__main__":
    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt.options.nlp_scaling_method = "user-scaling"
    m, solver = get_model()
