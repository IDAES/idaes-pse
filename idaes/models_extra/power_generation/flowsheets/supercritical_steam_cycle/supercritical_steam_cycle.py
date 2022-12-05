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

"""This is an example supercritical pulverized coal (SCPC) power plant steam cycle
model.  This model doesn't represent any specific power plant, but it represents
what could be considered a typical SCPC plant, producing around 620 MW gross.
This model is for demonstration and tutorial purposes only. Before looking at the
model, it may be useful to look at the process flow diagram (PFD).
"""

__author__ = "John Eslick, Maojian Wang"

# Import Python libraries
from collections import OrderedDict
import os
import argparse
import logging

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc, Port
from pyomo.common.fileutils import this_file_dir

# IDAES Imports
from idaes.core import FlowsheetBlock  # Flowsheet class
from idaes.core.util import model_serializer as ms  # load/save model state
from idaes.core.util.tags import svg_tag, ModelTagGroup  # place numbers/text in an SVG
from idaes.models.properties import iapws95  # steam properties
from idaes.models_extra.power_generation.unit_models.helm import (
    HelmTurbineMultistage,
    HelmMixer,
    HelmIsentropicCompressor,
    HelmTurbineStage,
    HelmNtuCondenser as Condenser,
)
from idaes.models_extra.power_generation.unit_models import FWH0D
from idaes.models.unit_models import (  # basic IDAES unit models, and enum
    MomentumMixingType,  # Enum type for mixer pressure calculation selection
)
from idaes.core.util.initialization import (
    propagate_state as _set_port,
)  # for model intialization
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

# Callback used to construct heat exchangers with the Underwood approx. for LMTD
from idaes.models.unit_models.heat_exchanger import (
    delta_temperature_underwood_callback,
)
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale


_log = idaeslog.getModelLogger(__name__, logging.INFO)


def create_model():
    """Create the flowsheet and add unit models. Fixing model inputs is done
    in a separate function to try to keep this fairly clean and easy to follow.

    Args:
        None

    Returns:
        (ConcreteModel) Steam cycle model
    """
    ############################################################################
    #  Flowsheet and Properties                                                #
    ############################################################################
    m = pyo.ConcreteModel(name="Steam Cycle Model")
    m.fs = FlowsheetBlock(dynamic=False)  # Add steady state flowsheet

    # A physical property parameter block for IAPWS-95 with pressure and enthalpy
    # (PH) state variables.  Usually pressure and enthalpy state variables are
    # more robust especially when the phases are unknown.
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        phase_presentation=iapws95.PhaseType.MIX
    )

    # A physical property parameter block with temperature, pressure and vapor
    # fraction (TPx) state variables. There are a few instances where the vapor
    # fraction is known and the temperature and pressure state variables are
    # preferable.
    m.fs.prop_water_tpx = iapws95.Iapws95ParameterBlock(
        phase_presentation=iapws95.PhaseType.LG, state_vars=iapws95.StateVars.TPX
    )
    ############################################################################
    #  Turbine with fill-in reheat constraints                                 #
    ############################################################################
    # The TurbineMultistage class allows creation of the full turbine model by
    # providing several configuration options, including: throttle valves;
    # high-, intermediate-, and low-pressure sections; steam extractions; and
    # pressure driven flow.  See the IDAES documentation for details.
    m.fs.turb = HelmTurbineMultistage(
        property_package=m.fs.prop_water,
        num_parallel_inlet_stages=4,
        num_hp=7,
        num_ip=10,
        num_lp=11,
        hp_split_locations=[4, 7],
        ip_split_locations=[5, 10],
        lp_split_locations=[4, 8, 10, 11],
        hp_disconnect=[7],
        ip_split_num_outlets={10: 3},
    )
    # This model is only the steam cycle, and the reheater is part of the boiler.
    # To fill in the reheater gap, a few constraints for the flow, pressure drop,
    # and outlet temperature are added. A detailed boiler model can be coupled later.
    #
    # hp_split[7] is the splitter directly after the last HP stage.  The splitter
    # outlet "outlet_1" is always taken to be the main steam flow through the turbine.
    # When the turbine model was instantiated the stream from the HP section to the IP
    # section was omitted, so the reheater could be inserted.

    # The flow constraint sets flow from outlet_1 of the splitter equal to
    # flow into the IP turbine.
    @m.fs.turb.Constraint(m.fs.time)
    def constraint_reheat_flow(b, t):
        return b.ip_stages[1].inlet.flow_mol[t] == b.hp_split[7].outlet_1.flow_mol[t]

    # Create a variable for pressure change in the reheater (assuming
    # reheat_delta_p should be negative).
    m.fs.turb.reheat_delta_p = pyo.Var(m.fs.time, initialize=0, units=pyo.units.Pa)

    # Add a constraint to calculate the IP section inlet pressure based on the
    # pressure drop in the reheater and the outlet pressure of the HP section.
    @m.fs.turb.Constraint(m.fs.time)
    def constraint_reheat_press(b, t):
        return (
            b.ip_stages[1].inlet.pressure[t]
            == b.hp_split[7].outlet_1.pressure[t] + b.reheat_delta_p[t]
        )

    # Create a variable for reheat temperature and fix it to the desired reheater
    # outlet temperature
    m.fs.turb.reheat_out_T = pyo.Var(m.fs.time, initialize=866, units=pyo.units.K)

    # Create a constraint for the IP section inlet temperature.
    @m.fs.turb.Constraint(m.fs.time)
    def constraint_reheat_temp(b, t):
        return (
            b.ip_stages[1].control_volume.properties_in[t].temperature
            == b.reheat_out_T[t]
        )

    ############################################################################
    #  Add Condenser/hotwell/condensate pump                                   #
    ############################################################################
    # Add a mixer for all the streams coming into the condenser.  In this case the
    # main steam, and the boiler feed pump turbine outlet go to the condenser
    m.fs.condenser_mix = HelmMixer(
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["main", "bfpt"],
        property_package=m.fs.prop_water,
    )
    # The pressure in the mixer comes from the connection to the condenser.  All
    # the streams coming in and going out of the mixer are equal, but we created
    # the mixer with no calculation for the unit pressure. Here a constraint that
    # specifies that the mixer pressure is equal to the main steam pressure is
    # added.  There is also a constraint that specifies the that BFP turbine outlet
    # pressure is the same as the condenser pressure.  Combined with the stream
    # connections between units, these constraints effectively specify that the
    # mixer inlet and outlet streams all have the same pressure.
    @m.fs.condenser_mix.Constraint(m.fs.time)
    def mixer_pressure_constraint(b, t):
        return b.main_state[t].pressure == b.mixed_state[t].pressure

    # PYLINT-TODO the name "mixer_pressure_constraint" is reused below to define other constraint functions,
    # causing pylint to report function-redefined errors
    # this likely does not actually cause issues at runtime,
    # but it could be worth to check anyway if the pylint errors can be addressed
    # e.g. by giving a unique name to each of the affected functions
    # pylint: disable=function-redefined

    # The condenser model uses the physical property model with TPx state
    # variables, while the rest of the model uses PH state variables. To
    # translate between the two property calculations, an extra port is added to
    # the mixer which contains temperature, pressure, and vapor fraction
    # quantities.
    m.fs.condenser_mix._flow_mol_ref = pyo.Reference(
        m.fs.condenser_mix.mixed_state[:].flow_mol
    )
    m.fs.condenser_mix._temperature_ref = pyo.Reference(
        m.fs.condenser_mix.mixed_state[:].temperature
    )
    m.fs.condenser_mix._pressure_ref = pyo.Reference(
        m.fs.condenser_mix.mixed_state[:].pressure
    )
    m.fs.condenser_mix._vapor_frac_ref = pyo.Reference(
        m.fs.condenser_mix.mixed_state[:].vapor_frac
    )

    m.fs.condenser_mix.outlet_tpx = Port(
        initialize={
            "flow_mol": m.fs.condenser_mix._flow_mol_ref,
            "temperature": m.fs.condenser_mix._temperature_ref,
            "pressure": m.fs.condenser_mix._pressure_ref,
            "vapor_frac": m.fs.condenser_mix._vapor_frac_ref,
        }
    )

    # Add NTU condenser model
    m.fs.condenser = Condenser(
        dynamic=False,
        shell={"has_pressure_change": False, "property_package": m.fs.prop_water},
        tube={"has_pressure_change": False, "property_package": m.fs.prop_water},
    )

    # Add the condenser hotwell.  In steady state a mixer will work.  This is
    # where makeup water is added if needed.
    m.fs.hotwell = HelmMixer(
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["condensate", "makeup"],
        property_package=m.fs.prop_water,
    )

    # The hotwell is assumed to be at the same pressure as the condenser.
    @m.fs.hotwell.Constraint(m.fs.time)
    def mixer_pressure_constraint(b, t):
        return b.condensate_state[t].pressure == b.mixed_state[t].pressure

    # Condensate pump (Use compressor model, since it is more robust if vapor form)
    m.fs.cond_pump = HelmIsentropicCompressor(property_package=m.fs.prop_water)
    ############################################################################
    #  Add low pressure feedwater heaters                                      #
    ############################################################################
    # All the feedwater heater sections will be set to use the Underwood
    # approximation for LMTD, so create the fwh_config dict to make the config
    # slightly cleaner
    fwh_config = {"delta_temperature_callback": delta_temperature_underwood_callback}

    # The feedwater heater model allows feedwater heaters with a desuperheat,
    # condensing, and subcooling section to be added an a reasonably simple way.
    # See the IDAES documentation for more information of configuring feedwater
    # heaters
    m.fs.fwh1 = FWH0D(
        has_desuperheat=False,
        has_drain_cooling=False,
        has_drain_mixer=True,
        property_package=m.fs.prop_water,
        condense=fwh_config,
    )
    # pump for fwh1 condensate, to pump it ahead and mix with feedwater
    m.fs.fwh1_pump = HelmIsentropicCompressor(property_package=m.fs.prop_water)
    # Mix the FWH1 drain back into the feedwater
    m.fs.fwh1_return = HelmMixer(
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["feedwater", "fwh1_drain"],
        property_package=m.fs.prop_water,
    )

    # Set the mixer pressure to the feedwater pressure
    @m.fs.fwh1_return.Constraint(m.fs.time)
    def mixer_pressure_constraint(b, t):
        return b.feedwater_state[t].pressure == b.mixed_state[t].pressure

    # Add the rest of the low pressure feedwater heaters
    m.fs.fwh2 = FWH0D(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=True,
        property_package=m.fs.prop_water,
        desuperheat=fwh_config,
        cooling=fwh_config,
        condense=fwh_config,
    )
    m.fs.fwh3 = FWH0D(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=True,
        property_package=m.fs.prop_water,
        desuperheat=fwh_config,
        cooling=fwh_config,
        condense=fwh_config,
    )
    m.fs.fwh4 = FWH0D(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=False,
        property_package=m.fs.prop_water,
        desuperheat=fwh_config,
        cooling=fwh_config,
        condense=fwh_config,
    )
    ############################################################################
    #  Add deaerator and boiler feed pump (BFP)                                #
    ############################################################################
    # The deaerator is basically an open tank with multiple inlets.  For steady-
    # state, a mixer model is sufficient.
    m.fs.fwh5_da = HelmMixer(
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["steam", "drain", "feedwater"],
        property_package=m.fs.prop_water,
    )

    @m.fs.fwh5_da.Constraint(m.fs.time)
    def mixer_pressure_constraint(b, t):
        # Not sure about deaerator pressure, so assume same as feedwater inlet
        return b.feedwater_state[t].pressure == b.mixed_state[t].pressure

    # Add the boiler feed pump and boiler feed pump turbine
    m.fs.bfp = HelmIsentropicCompressor(property_package=m.fs.prop_water)

    m.fs.bfpt = HelmTurbineStage(property_package=m.fs.prop_water)

    # The boiler feed pump outlet pressure is the same as the condenser
    @m.fs.Constraint(m.fs.time)
    def constraint_out_pressure(b, t):
        return (
            b.bfpt.control_volume.properties_out[t].pressure
            == b.condenser.shell.properties_out[t].pressure
        )

    # Instead of specifying a fixed efficiency, specify that the steam is just
    # starting to condense at the outlet of the boiler feed pump turbine.  This
    # ensures approximately the right behavior in the turbine.  With a fixed
    # efficiency, depending on the conditions you can get odd things like steam
    # fully condensing in the turbine.
    @m.fs.Constraint(m.fs.time)
    def constraint_out_enthalpy(b, t):
        return (
            b.bfpt.control_volume.properties_out[t].enth_mol
            == b.bfpt.control_volume.properties_out[t].enth_mol_sat_phase["Vap"]
            - 200 * pyo.units.J / pyo.units.mol
        )

    # The boiler feed pump power is the same as the power generated by the
    # boiler feed pump turbine. This constraint determines the steam flow to the
    # BFP turbine. The turbine work is negative for power out, while pump work
    # is positive for power in.
    @m.fs.Constraint(m.fs.time)
    def constraint_bfp_power(b, t):
        return 0 == b.bfp.control_volume.work[t] + b.bfpt.control_volume.work[t]

    ############################################################################
    #  Add high pressure feedwater heaters                                     #
    ############################################################################
    m.fs.fwh6 = FWH0D(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=True,
        property_package=m.fs.prop_water,
        desuperheat=fwh_config,
        cooling=fwh_config,
        condense=fwh_config,
    )
    m.fs.fwh7 = FWH0D(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=True,
        property_package=m.fs.prop_water,
        desuperheat=fwh_config,
        cooling=fwh_config,
        condense=fwh_config,
    )
    m.fs.fwh8 = FWH0D(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=False,
        property_package=m.fs.prop_water,
        desuperheat=fwh_config,
        cooling=fwh_config,
        condense=fwh_config,
    )
    ############################################################################
    #  Additional Constraints/Expressions                                      #
    ############################################################################

    # Add a few constraints to allow a for complete plant results despite the
    # lack of a detailed boiler model.

    # Boiler pressure drop
    m.fs.boiler_pressure_drop_fraction = pyo.Var(
        m.fs.time,
        initialize=0.01,
        doc="Fraction of pressure lost from boiler feed pump and turbine inlet",
    )

    @m.fs.Constraint(m.fs.time)
    def boiler_pressure_drop(b, t):
        return (
            m.fs.bfp.control_volume.properties_out[t].pressure
            * (1 - b.boiler_pressure_drop_fraction[t])
            == m.fs.turb.inlet_split.mixed_state[t].pressure
        )

    # Again, since the boiler is missing, set the flow of steam into the turbine
    # equal to the flow of feedwater out of the last feedwater heater.
    @m.fs.Constraint(m.fs.time)
    def close_flow(b, t):
        return (
            m.fs.bfp.control_volume.properties_out[t].flow_mol
            == m.fs.turb.inlet_split.mixed_state[t].flow_mol
        )

    # Calculate the amount of heat that is added in the boiler, including the
    # reheater.
    @m.fs.Expression(m.fs.time)
    def boiler_heat(b, t):
        return (
            b.turb.inlet_split.mixed_state[t].enth_mol
            * b.turb.inlet_split.mixed_state[t].flow_mol
            - b.fwh8.desuperheat.cold_side.properties_out[t].enth_mol
            * b.fwh8.desuperheat.cold_side.properties_out[t].flow_mol
            + b.turb.ip_stages[1].control_volume.properties_in[t].enth_mol
            * b.turb.ip_stages[1].control_volume.properties_in[t].flow_mol
            - b.turb.hp_split[7].outlet_1.enth_mol[t]
            * b.turb.hp_split[7].outlet_1.flow_mol[t]
        )

    # Calculate the efficiency of the steam cycle.  This doesn't account for
    # heat loss in the boiler, so actual plant efficiency would be lower.
    @m.fs.Expression(m.fs.time)
    def steam_cycle_eff(b, t):
        return -100 * b.turb.power[t] / b.boiler_heat[t]

    ############################################################################
    ##  Create the stream Arcs                                                ##
    ############################################################################

    ############################################################################
    #  Connect turbine and condenser units                                     #
    ############################################################################
    m.fs.EXHST_MAIN = Arc(
        source=m.fs.turb.outlet_stage.outlet, destination=m.fs.condenser_mix.main
    )

    m.fs.condenser_mix_to_condenser = Arc(
        source=m.fs.condenser_mix.outlet, destination=m.fs.condenser.hot_side_inlet
    )

    m.fs.COND_01 = Arc(
        source=m.fs.condenser.hot_side_outlet, destination=m.fs.hotwell.condensate
    )

    m.fs.COND_02 = Arc(source=m.fs.hotwell.outlet, destination=m.fs.cond_pump.inlet)
    ############################################################################
    #  Low pressure FWHs                                                       #
    ############################################################################
    m.fs.EXTR_LP11 = Arc(
        source=m.fs.turb.lp_split[11].outlet_2, destination=m.fs.fwh1.drain_mix.steam
    )
    m.fs.COND_03 = Arc(
        source=m.fs.cond_pump.outlet, destination=m.fs.fwh1.condense.cold_side_inlet
    )
    m.fs.FWH1_DRN1 = Arc(
        source=m.fs.fwh1.condense.hot_side_outlet, destination=m.fs.fwh1_pump.inlet
    )
    m.fs.FWH1_DRN2 = Arc(
        source=m.fs.fwh1_pump.outlet, destination=m.fs.fwh1_return.fwh1_drain
    )
    m.fs.FW01A = Arc(
        source=m.fs.fwh1.condense.cold_side_outlet,
        destination=m.fs.fwh1_return.feedwater,
    )
    # fwh2
    m.fs.FW01B = Arc(
        source=m.fs.fwh1_return.outlet, destination=m.fs.fwh2.cooling.cold_side_inlet
    )
    m.fs.FWH2_DRN = Arc(
        source=m.fs.fwh2.cooling.hot_side_outlet, destination=m.fs.fwh1.drain_mix.drain
    )
    m.fs.EXTR_LP10 = Arc(
        source=m.fs.turb.lp_split[10].outlet_2,
        destination=m.fs.fwh2.desuperheat.hot_side_inlet,
    )
    # fwh3
    m.fs.FW02 = Arc(
        source=m.fs.fwh2.desuperheat.cold_side_outlet,
        destination=m.fs.fwh3.cooling.cold_side_inlet,
    )
    m.fs.FWH3_DRN = Arc(
        source=m.fs.fwh3.cooling.hot_side_outlet, destination=m.fs.fwh2.drain_mix.drain
    )
    m.fs.EXTR_LP8 = Arc(
        source=m.fs.turb.lp_split[8].outlet_2,
        destination=m.fs.fwh3.desuperheat.hot_side_inlet,
    )
    # fwh4
    m.fs.FW03 = Arc(
        source=m.fs.fwh3.desuperheat.cold_side_outlet,
        destination=m.fs.fwh4.cooling.cold_side_inlet,
    )
    m.fs.FWH4_DRN = Arc(
        source=m.fs.fwh4.cooling.hot_side_outlet, destination=m.fs.fwh3.drain_mix.drain
    )
    m.fs.EXTR_LP4 = Arc(
        source=m.fs.turb.lp_split[4].outlet_2,
        destination=m.fs.fwh4.desuperheat.hot_side_inlet,
    )
    ############################################################################
    #  FWH5 (Deaerator) and boiler feed pump (BFP)                             #
    ############################################################################
    m.fs.FW04 = Arc(
        source=m.fs.fwh4.desuperheat.cold_side_outlet,
        destination=m.fs.fwh5_da.feedwater,
    )
    m.fs.EXTR_IP10 = Arc(
        source=m.fs.turb.ip_split[10].outlet_2, destination=m.fs.fwh5_da.steam
    )
    m.fs.FW05A = Arc(source=m.fs.fwh5_da.outlet, destination=m.fs.bfp.inlet)
    m.fs.EXTR_BFPT_A = Arc(
        source=m.fs.turb.ip_split[10].outlet_3, destination=m.fs.bfpt.inlet
    )
    m.fs.EXHST_BFPT = Arc(source=m.fs.bfpt.outlet, destination=m.fs.condenser_mix.bfpt)
    ############################################################################
    #  High-pressure feedwater heaters                                         #
    ############################################################################
    # fwh6
    m.fs.FW05B = Arc(
        source=m.fs.bfp.outlet, destination=m.fs.fwh6.cooling.cold_side_inlet
    )
    m.fs.FWH6_DRN = Arc(
        source=m.fs.fwh6.cooling.hot_side_outlet, destination=m.fs.fwh5_da.drain
    )
    m.fs.EXTR_IP5 = Arc(
        source=m.fs.turb.ip_split[5].outlet_2,
        destination=m.fs.fwh6.desuperheat.hot_side_inlet,
    )
    # fwh7
    m.fs.FW06 = Arc(
        source=m.fs.fwh6.desuperheat.cold_side_outlet,
        destination=m.fs.fwh7.cooling.cold_side_inlet,
    )
    m.fs.FWH7_DRN = Arc(
        source=m.fs.fwh7.cooling.hot_side_outlet, destination=m.fs.fwh6.drain_mix.drain
    )
    m.fs.EXTR_HP7 = Arc(
        source=m.fs.turb.hp_split[7].outlet_2,
        destination=m.fs.fwh7.desuperheat.hot_side_inlet,
    )
    # fwh8
    m.fs.FW07 = Arc(
        source=m.fs.fwh7.desuperheat.cold_side_outlet,
        destination=m.fs.fwh8.cooling.cold_side_inlet,
    )
    m.fs.FWH8_DRN = Arc(
        source=m.fs.fwh8.cooling.hot_side_outlet, destination=m.fs.fwh7.drain_mix.drain
    )
    m.fs.EXTR_HP4 = Arc(
        source=m.fs.turb.hp_split[4].outlet_2,
        destination=m.fs.fwh8.desuperheat.hot_side_inlet,
    )

    ############################################################################
    # Turn the Arcs into constraints and return the model                      #
    ############################################################################
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    return m


def _stream_dict(m):
    """Adds _streams to m, which contains a dictionary of streams for display

    Args:
        m (ConcreteModel): A Pyomo model from create_model()

    Returns:
        None
    """

    # Put all the top-level streams in the stream dictionary
    m._streams = dict(
        [(c.getname(), c) for c in m.fs.component_objects(Arc, descend_into=False)]
    )

    # There are some additional streams we are interested in that are either
    # inlets or outlets, where there is no arc or arcs buried in unit models.
    # Next those are added specifically.
    m._streams.update(
        {
            "STEAM_MAIN": m.fs.turb.inlet_split.mixed_state,
            "THRTL1": m.fs.turb.inlet_stage[1].control_volume.properties_in,
            "THRTL2": m.fs.turb.inlet_stage[2].control_volume.properties_in,
            "THRTL3": m.fs.turb.inlet_stage[3].control_volume.properties_in,
            "THRTL4": m.fs.turb.inlet_stage[4].control_volume.properties_in,
            "RHT_COLD": m.fs.turb.hp_split[7].outlet_1_state,
            "RHT_HOT": m.fs.turb.ip_stages[1].control_volume.properties_in,
            "STEAM_LP": m.fs.turb.ip_split[10].outlet_1_state,
            "CW01": m.fs.condenser.tube.properties_in,
            "CW02": m.fs.condenser.tube.properties_out,
            "MAKEUP_01": m.fs.hotwell.makeup_state,
            "FW08": m.fs.fwh8.desuperheat.cold_side.properties_out,
        }
    )

    # Alphabetize to find streams easier in the tabular stream table view.
    m._streams = OrderedDict(sorted(m._streams.items()))


def set_model_input(m):
    """Fix some variables and set values. Generally get the model ready to run
    in simulation mode (0 degrees of freedom).

    Args:
        m (ConcreteModel): A Pyomo model from create_model()

    Returns:
        None
    """
    ############################################################################
    #  Turbine input                                                           #
    ############################################################################
    main_steam_pressure = 2.423e7
    m.fs.turb.turbine_inlet_cf_fix(4.6e-4)
    m.fs.turb.turbine_outlet_cf_fix(0.32)
    # Set the turbine steam inlet conditions and flow guess for init
    m.fs.turb.inlet_split.inlet.enth_mol.fix(62710)
    m.fs.turb.inlet_split.inlet.pressure.fix(main_steam_pressure)
    m.fs.boiler_pressure_drop_fraction.fix(0.1)
    m.fs.turb.inlet_split.inlet.flow_mol[:].value = 21000
    m.fs.turb.inlet_split.inlet.flow_mol.unfix()  # Pressure-driven
    m.fs.turb.inlet_mix.use_equal_pressure_constraint()  # Pressure-driven mix
    # Set throttle valve
    m.fs.turb.throttle_cv_fix(0.0009)
    m.fs.turb.throttle_valve[1].valve_opening.fix(0.9)
    m.fs.turb.throttle_valve[2].valve_opening.fix(0.9)
    m.fs.turb.throttle_valve[3].valve_opening.fix(0.9)
    m.fs.turb.throttle_valve[4].valve_opening.fix(0.9)
    # Set the efficiency and pressure ratios of stages other than inlet and outlet
    for i, s in m.fs.turb.hp_stages.items():
        s.ratioP.fix(0.80)
        s.efficiency_isentropic.fix(0.9)
    for i, s in m.fs.turb.ip_stages.items():
        s.ratioP.fix(0.78)
        s.efficiency_isentropic.fix(0.9)
    for i, s in m.fs.turb.lp_stages.items():
        s.ratioP[:].fix(0.76)
        s.efficiency_isentropic[:].fix(0.9)
    ############################################################################
    #  Reheater                                                                #
    ############################################################################
    m.fs.turb.reheat_delta_p.fix(0)
    m.fs.turb.reheat_out_T.fix(866)
    ############################################################################
    #  Condenser section inputs                                                #
    ############################################################################
    m.fs.condenser.cold_side_inlet.flow_mol.fix(2500000)
    m.fs.condenser.cold_side_inlet.enth_mol.fix(1700)
    m.fs.condenser.cold_side_inlet.pressure.fix(500000)
    m.fs.condenser.area.fix(13000)
    m.fs.condenser.overall_heat_transfer_coefficient.fix(15000)

    m.fs.hotwell.makeup.flow_mol[:].value = 1  # don't fix is calculated
    m.fs.hotwell.makeup.enth_mol.fix(2500)
    m.fs.hotwell.makeup.pressure.fix(101325)
    m.fs.cond_pump.efficiency_isentropic.fix(0.80)
    m.fs.cond_pump.deltaP.fix(1e6)
    ############################################################################
    #  Low-pressure FWH section inputs                                         #
    ############################################################################
    # fwh1
    # Heat transfer coefficent correlation constraints can be added to the
    # feedwater heaters, but to keep this example simple, they are fixed
    # constant values here.
    m.fs.fwh1.condense.area.fix(400)
    m.fs.fwh1.condense.overall_heat_transfer_coefficient.fix(2000)
    # fwh1 pump
    m.fs.fwh1_pump.efficiency_isentropic.fix(0.80)
    m.fs.fwh1_pump.deltaP.fix(1.2e6)  # need pressure higher than feedwater
    # fwh2
    m.fs.fwh2.condense.area.fix(150)
    m.fs.fwh2.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh2.desuperheat.area.fix(75)
    m.fs.fwh2.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh2.cooling.area.fix(75)
    m.fs.fwh2.cooling.overall_heat_transfer_coefficient.fix(300)
    # fwh3
    m.fs.fwh3.condense.area.fix(100)
    m.fs.fwh3.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh3.desuperheat.area.fix(50)
    m.fs.fwh3.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh3.cooling.area.fix(50)
    m.fs.fwh3.cooling.overall_heat_transfer_coefficient.fix(300)
    # fwh4
    m.fs.fwh4.condense.area.fix(100)
    m.fs.fwh4.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh4.desuperheat.area.fix(50)
    m.fs.fwh4.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh4.cooling.area.fix(50)
    m.fs.fwh4.cooling.overall_heat_transfer_coefficient.fix(300)
    ############################################################################
    #  Deaerator and boiler feed pump (BFP) input                              #
    ############################################################################
    m.fs.bfp.efficiency_isentropic.fix(0.80)
    m.fs.bfp.outlet.pressure[:].value = main_steam_pressure * 1.1  # guess
    m.fs.bfpt.efficiency_isentropic.value = 0.80  # don't fix, just initial guess
    ############################################################################
    #  High-pressure feedwater heater                                          #
    ############################################################################
    # fwh6
    m.fs.fwh6.condense.area.fix(300)
    m.fs.fwh6.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh6.desuperheat.area.fix(150)
    m.fs.fwh6.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh6.cooling.area.fix(150)
    m.fs.fwh6.cooling.overall_heat_transfer_coefficient.fix(300)
    # fwh7
    m.fs.fwh7.condense.area.fix(200)
    m.fs.fwh7.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh7.desuperheat.area.fix(100)
    m.fs.fwh7.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh7.cooling.area.fix(100)
    m.fs.fwh7.cooling.overall_heat_transfer_coefficient.fix(300)
    # fwh8
    m.fs.fwh8.condense.area.fix(200)
    m.fs.fwh8.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh8.desuperheat.area.fix(100)
    m.fs.fwh8.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh8.cooling.area.fix(100)
    m.fs.fwh8.cooling.overall_heat_transfer_coefficient.fix(300)
    # Now all the model input has been specified.


def initialize(m, fileinput=None, outlvl=idaeslog.NOTSET):
    """Initialize a mode from create_model(), set model inputs before
    initializing.

    Args:
        m (ConcreteModel): A Pyomo model from create_model()
        fileinput (str|None): File to load initialized model state from. If a
            file is supplied skip initialization routine. If None, initialize.

    Returns:
        solver: A Pyomo solver object, that can be used to solve the model.

    """
    init_log = idaeslog.getInitLogger(m.name, outlvl, tag="flowsheet")
    solve_log = idaeslog.getSolveLogger(m.name, outlvl, tag="flowsheet")

    # set scaling factors

    iscale.set_scaling_factor(m.fs.condenser.hot_side.heat, 1e-9)
    iscale.set_scaling_factor(m.fs.condenser.cold_side.heat, 1e-9)

    iscale.set_scaling_factor(m.fs.fwh1.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.fwh1.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(m.fs.fwh2.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.fwh2.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(m.fs.fwh3.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.fwh3.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(m.fs.fwh4.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.fwh4.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(m.fs.fwh6.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.fwh6.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(m.fs.fwh7.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.fwh7.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(m.fs.fwh8.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(m.fs.fwh8.condense.cold_side.heat, 1e-7)

    iscale.calculate_scaling_factors(m)

    solver = get_solver()

    if fileinput is not None:
        init_log.info("Loading initial values from file: {}".format(fileinput))
        ms.from_json(m, fname=fileinput)
        return solver

    init_log.info("Starting initialization")
    ############################################################################
    #  Initialize turbine                                                      #
    ############################################################################
    # Extraction rates are calculated from the feedwater heater models, so to
    # initialize the turbine fix some initial guesses. They get unfixed after
    # solving the turbine
    m.fs.turb.outlet_stage.control_volume.properties_out[:].pressure.fix(3500)
    m.fs.turb.lp_split[11].split_fraction[:, "outlet_2"].fix(0.04403)
    m.fs.turb.lp_split[10].split_fraction[:, "outlet_2"].fix(0.04025)
    m.fs.turb.lp_split[8].split_fraction[:, "outlet_2"].fix(0.04362)
    m.fs.turb.lp_split[4].split_fraction[:, "outlet_2"].fix(0.08025)
    m.fs.turb.ip_split[10].split_fraction[:, "outlet_2"].fix(0.045)
    m.fs.turb.ip_split[10].split_fraction[:, "outlet_3"].fix(0.04)
    m.fs.turb.ip_split[5].split_fraction[:, "outlet_2"].fix(0.05557)
    m.fs.turb.hp_split[7].split_fraction[:, "outlet_2"].fix(0.09741)
    m.fs.turb.hp_split[4].split_fraction[:, "outlet_2"].fix(0.0740)
    # Put in a rough initial guess for the IP section inlet, since it is
    # disconnected from the HP section for the reheater.
    ip1_pin = 5.35e6
    ip1_hin = pyo.value(iapws95.htpx(T=866 * pyo.units.K, P=ip1_pin * pyo.units.Pa))
    ip1_fin = pyo.value(m.fs.turb.inlet_split.inlet.flow_mol[0])
    m.fs.turb.ip_stages[1].inlet.enth_mol[:].value = ip1_hin
    m.fs.turb.ip_stages[1].inlet.flow_mol[:].value = ip1_fin
    m.fs.turb.ip_stages[1].inlet.pressure[:].value = ip1_pin
    # initialize turbine
    assert degrees_of_freedom(m.fs.turb) == 0
    m.fs.turb.initialize(outlvl=outlvl, optarg=solver.options)
    with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        res = solver.solve(m.fs.turb, tee=slc.tee)
    init_log.info("Full turbine solve complete: {}".format(idaeslog.condition(res)))
    # The turbine outlet pressure is determined by the condenser once hooked up
    m.fs.turb.outlet_stage.control_volume.properties_out[:].pressure.unfix()
    # Extraction rates are calculated once the feedwater heater models are
    # added so unfix the splits.
    m.fs.turb.lp_split[11].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.lp_split[10].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.lp_split[8].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.lp_split[4].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.ip_split[10].split_fraction[:, "outlet_3"].unfix()
    m.fs.turb.ip_split[5].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.hp_split[7].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.hp_split[4].split_fraction[:, "outlet_2"].unfix()
    # Initialize the boiler feed pump turbine.
    _set_port(m.fs.bfpt.inlet, m.fs.turb.ip_split[10].outlet_3)
    m.fs.bfpt.control_volume.properties_out[:].pressure.fix(10000)
    m.fs.bfpt.efficiency_isentropic.fix()
    m.fs.bfpt.initialize(outlvl=idaeslog.DEBUG, optarg=solver.options)
    m.fs.bfpt.control_volume.properties_out[:].pressure.unfix()
    m.fs.bfpt.efficiency_isentropic.unfix()
    ############################################################################
    #  Condenser section                                                       #
    ############################################################################
    # initialize condenser mixer
    _set_port(m.fs.condenser_mix.main, m.fs.turb.outlet_stage.outlet)
    _set_port(m.fs.condenser_mix.bfpt, m.fs.bfpt.outlet)
    m.fs.condenser_mix.initialize(outlvl=outlvl, optarg=solver.options)
    # initialize condenser hx

    _set_port(m.fs.condenser.hot_side_inlet, m.fs.condenser_mix.outlet)
    _set_port(m.fs.condenser.cold_side_outlet, m.fs.condenser.cold_side_inlet)
    m.fs.condenser.initialize(outlvl=outlvl, optarg=solver.options)

    # initialize hotwell
    _set_port(m.fs.hotwell.condensate, m.fs.condenser.hot_side_outlet)

    m.fs.hotwell.initialize(outlvl=outlvl, optarg=solver.options)
    m.fs.hotwell.condensate.unfix()
    # initialize condensate pump
    _set_port(m.fs.cond_pump.inlet, m.fs.hotwell.outlet)
    m.fs.cond_pump.initialize(outlvl=outlvl, optarg=solver.options)
    ############################################################################
    #  Low-pressure FWH section                                                #
    ############################################################################
    # fwh1
    m.fs.fwh1.drain_mix.drain.flow_mol[:] = 1000
    m.fs.fwh1.drain_mix.drain.pressure[:] = 1e5
    m.fs.fwh1.drain_mix.drain.enth_mol[:] = 6117
    _set_port(m.fs.fwh1.condense.cold_side_inlet, m.fs.cond_pump.outlet)
    _set_port(m.fs.fwh1.drain_mix.steam, m.fs.turb.lp_split[11].outlet_2)
    m.fs.fwh1.initialize(outlvl=outlvl, optarg=solver.options)
    # initialize fwh1 drain pump
    _set_port(m.fs.fwh1_pump.inlet, m.fs.fwh1.condense.hot_side_outlet)
    m.fs.fwh1_pump.initialize(outlvl=5, optarg=solver.options)
    # initialize mixer to add fwh1 drain to feedwater
    _set_port(m.fs.fwh1_return.feedwater, m.fs.fwh1.condense.cold_side_outlet)
    _set_port(m.fs.fwh1_return.fwh1_drain, m.fs.fwh1.condense.hot_side_outlet)
    m.fs.fwh1_return.initialize(outlvl=outlvl, optarg=solver.options)
    m.fs.fwh1_return.feedwater.unfix()
    m.fs.fwh1_return.fwh1_drain.unfix()
    # fwh2
    m.fs.fwh2.drain_mix.drain.flow_mol[:] = 100
    m.fs.fwh2.drain_mix.drain.pressure[:] = 1.5e5
    m.fs.fwh2.drain_mix.drain.enth_mol[:] = 7000
    _set_port(m.fs.fwh2.cooling.cold_side_inlet, m.fs.fwh1_return.outlet)
    _set_port(m.fs.fwh2.desuperheat.hot_side_inlet, m.fs.turb.lp_split[10].outlet_2)
    m.fs.fwh2.initialize(outlvl=outlvl, optarg=solver.options)
    # fwh3
    m.fs.fwh3.drain_mix.drain.flow_mol[:] = 100
    m.fs.fwh3.drain_mix.drain.pressure[:] = 2.5e5
    m.fs.fwh3.drain_mix.drain.enth_mol[:] = 8000
    _set_port(m.fs.fwh3.cooling.cold_side_inlet, m.fs.fwh2.desuperheat.cold_side_outlet)
    _set_port(m.fs.fwh3.desuperheat.hot_side_inlet, m.fs.turb.lp_split[8].outlet_2)
    m.fs.fwh3.initialize(outlvl=outlvl, optarg=solver.options)
    # fwh4
    _set_port(m.fs.fwh4.cooling.cold_side_inlet, m.fs.fwh3.desuperheat.cold_side_outlet)
    _set_port(m.fs.fwh4.desuperheat.hot_side_inlet, m.fs.turb.lp_split[4].outlet_2)
    m.fs.fwh4.initialize(outlvl=outlvl, optarg=solver.options)
    ############################################################################
    #  boiler feed pump and deaerator                                          #
    ############################################################################
    _set_port(m.fs.fwh5_da.feedwater, m.fs.fwh4.desuperheat.cold_side_outlet)
    _set_port(m.fs.fwh5_da.steam, m.fs.turb.ip_split[10].outlet_2)
    m.fs.fwh5_da.drain.flow_mol[:] = 2000
    m.fs.fwh5_da.drain.pressure[:] = 3e6
    m.fs.fwh5_da.drain.enth_mol[:] = 9000
    m.fs.fwh5_da.initialize(outlvl=outlvl, optarg=solver.options)
    _set_port(m.fs.bfp.inlet, m.fs.fwh5_da.outlet)
    m.fs.bfp.control_volume.properties_out[:].pressure.fix()
    m.fs.bfp.initialize(outlvl=outlvl, optarg=solver.options)
    m.fs.bfp.control_volume.properties_out[:].pressure.unfix()
    ############################################################################
    #  High-pressure feedwater heaters                                         #
    ############################################################################
    # fwh6
    m.fs.fwh6.drain_mix.drain.flow_mol[:] = 1000
    m.fs.fwh6.drain_mix.drain.pressure[:] = 1e7
    m.fs.fwh6.drain_mix.drain.enth_mol[:] = 9500
    _set_port(m.fs.fwh6.cooling.cold_side_inlet, m.fs.bfp.outlet)
    _set_port(m.fs.fwh6.desuperheat.hot_side_inlet, m.fs.turb.ip_split[5].outlet_2)
    m.fs.fwh6.initialize(outlvl=outlvl, optarg=solver.options)
    # fwh7
    m.fs.fwh7.drain_mix.drain.flow_mol[:] = 2000
    m.fs.fwh7.drain_mix.drain.pressure[:] = 1e7
    m.fs.fwh7.drain_mix.drain.enth_mol[:] = 9500
    _set_port(m.fs.fwh7.cooling.cold_side_inlet, m.fs.fwh6.desuperheat.cold_side_outlet)
    _set_port(m.fs.fwh7.desuperheat.hot_side_inlet, m.fs.turb.hp_split[7].outlet_2)
    m.fs.fwh7.initialize(outlvl=outlvl, optarg=solver.options)
    # fwh8
    _set_port(m.fs.fwh8.cooling.cold_side_inlet, m.fs.fwh7.desuperheat.cold_side_outlet)
    _set_port(m.fs.fwh8.desuperheat.hot_side_inlet, m.fs.turb.hp_split[4].outlet_2)
    m.fs.fwh8.initialize(outlvl=outlvl, optarg=solver.options)
    with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        res = solver.solve(m, tee=slc.tee)
    init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    return solver


def pfd_result(m, df, svg):
    """Insert model results into the PFD and return a new SVG string, which
    can be displayed, further edited, or saved to a file.

    Args:
        m (ConcreteModel): A steam cycle model
        df (Pandas DataFrame): Stream table
        svg (FILE*, str, bytes): Original svg as either a file-like object,
            a string, or a byte array.

    Returns:
        (str): SVG content.
    """
    tags = {}  # dict of tags and data to insert into SVG
    for i in df.index:  # Create entires for streams
        tags[i + "_F"] = df.loc[i, "Molar Flow"]
        tags[i + "_T"] = df.loc[i, "T"]
        tags[i + "_P"] = df.loc[i, "P"]
        tags[i + "_X"] = df.loc[i, "Vapor Fraction"]
    # Add some additional quntities from the model to report
    tags["gross_power"] = -pyo.value(m.fs.turb.power[0])
    tags["gross_power_mw"] = -pyo.value(m.fs.turb.power[0]) * 1e-6
    tags["steam_mass_flow"] = df.loc["STEAM_MAIN", "Mass Flow"]
    tags["sc_eff"] = pyo.value(m.fs.steam_cycle_eff[0])
    tags["boiler_heat"] = pyo.value(m.fs.boiler_heat[0]) * 1e-6
    tags["steam_pressure"] = df.loc["STEAM_MAIN", "P"] / 1000.0
    tags["cond_pressure"] = df.loc["EXHST_MAIN", "P"] / 1000.0
    tags["bfp_power"] = pyo.value(m.fs.bfp.work_mechanical[0])
    tags["bfp_eff"] = pyo.value(m.fs.bfp.efficiency_isentropic[0]) * 100
    tags["bfpt_power"] = pyo.value(m.fs.bfpt.work_mechanical[0])
    tags["bfpt_eff"] = pyo.value(m.fs.bfpt.efficiency_isentropic[0]) * 100

    tag_group = ModelTagGroup()
    for t, v in tags.items():
        tag_group.add(t, v, format_string="{:.3f}")
    if svg is None:
        fname = os.path.join(this_file_dir(), "supercritical_steam_cycle.svg")
        with open(fname, "r") as f:
            svg = f.read()
    s = svg_tag(tag_group=tag_group, svg=svg)
    return s


def main(initialize_from_file=None, store_initialization=None):
    """Create and initalize a model and solver

    Args:
        None

    Returns:
        A tuple of a model and solver
    """
    m = create_model()
    _stream_dict(m)
    set_model_input(m)
    if initialize_from_file is None:
        solver = initialize(m)
    else:
        solver = initialize(m, fileinput=initialize_from_file)
    solver.solve(m, tee=True)
    if store_initialization is not None:
        ms.to_json(m, fname=store_initialization)
    return m, solver


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--initialize_from_file",
        help="File from which to load initialized values. If specified, the "
        "initialization proceedure will be skipped.",
        default=None,
    )
    parser.add_argument(
        "--store_initialization",
        help="If specified, store" " the initialized model values, to reload later.",
        default=None,
    )
    args = parser.parse_args()
    m, solver = main(
        initialize_from_file=args.initialize_from_file,
        store_initialization=args.store_initialization,
    )
