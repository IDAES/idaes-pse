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
"""
This is an example natural gas fuel cell (NGFC) power plant without CCS.
It uses a reduced order model created by PNNL to calculate the performance
of the solid oxide fuel cell based on a set a flowsheet conditions. More
details can be found in the associated jupyter notebook.
"""

import os
from collections import OrderedDict

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.network import Arc
from pyomo.common.fileutils import this_file_dir
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.util.infeasible import log_infeasible_constraints

# IDAES Imports
from idaes.core import FlowsheetBlock
from idaes.core.util import copy_port_values
from idaes.core.util import model_serializer as ms
from idaes.core.util.misc import svg_tag
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import create_stream_table_dataframe

import idaes.core.util.scaling as iscale

from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)
from idaes.generic_models.properties.core.generic.generic_reaction import (
    GenericReactionParameterBlock)

from idaes.generic_models.unit_models import (
    Mixer,
    Heater,
    HeatExchanger,
    PressureChanger,
    GibbsReactor,
    StoichiometricReactor,
    Separator,
    Translator)
from idaes.generic_models.unit_models.heat_exchanger import \
    delta_temperature_underwood_callback
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
from idaes.generic_models.unit_models.separator import SplittingType
from idaes.generic_models.unit_models.mixer import MomentumMixingType

from idaes.power_generation.properties.natural_gas_PR import get_prop, get_rxn
from idaes.power_generation.properties.NGFC.ROM.SOFC_ROM import \
    build_SOFC_ROM, initialize_SOFC_ROM

import logging


def build_power_island(m):
    # create property packages - 3 property packages and 1 reaction
    NG_config = get_prop(
        components=['H2', 'CO', "H2O", 'CO2', 'CH4', "C2H6", "C3H8", "C4H10",
                    'N2', 'O2', 'Ar'])
    m.fs.NG_props = GenericParameterBlock(default=NG_config)

    syn_config = get_prop(
        components=["H2", "CO", "H2O", "CO2", "CH4", "N2", "O2", "Ar"])
    m.fs.syn_props = GenericParameterBlock(default=syn_config)

    air_config = get_prop(
        components=['H2O', 'CO2', 'N2', 'O2', 'Ar'])
    m.fs.air_props = GenericParameterBlock(default=air_config)

    m.fs.rxn_props = GenericReactionParameterBlock(
        default=get_rxn(
            m.fs.syn_props, reactions=["h2_cmb", "co_cmb", "ch4_cmb"]))

    # build anode side units
    m.fs.anode_mix = Mixer(
        default={"inlet_list": ["feed", "recycle"],
                 "property_package": m.fs.NG_props})

    m.fs.anode_hx = HeatExchanger(
        default={"delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.NG_props,
                          "has_pressure_change": True}})

    m.fs.prereformer = GibbsReactor(
        default={"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "inert_species": ["N2", "Ar", "O2"],
                 "property_package": m.fs.NG_props})

    m.fs.anode_translator = Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": m.fs.NG_props,
                 "outlet_property_package": m.fs.syn_props})

    # C2H6, C3H8, and C4H10 are not present in the system after the prereformer
    # the property package is changed to improve model robustness
    @m.fs.anode_translator.Constraint(m.fs.time)
    def anode_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.anode_translator.Constraint(m.fs.time)
    def anode_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.anode_translator.Constraint(m.fs.time)
    def anode_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.anode_translator.Constraint(m.fs.time, m.fs.syn_props.component_list)
    def anode_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    m.fs.fuel_cell_mix = Mixer(
        default={"inlet_list": ["fuel_inlet", "ion_inlet"],
                 "momentum_mixing_type": MomentumMixingType.none,
                 "property_package": m.fs.syn_props})

    # outlet pressure should be equal to fuel inlet pressure because oxygen is
    # traveling through a solid electrolyte
    @m.fs.fuel_cell_mix.Constraint(m.fs.time)
    def ion_mix_constraint(b, t):
        return b.outlet.pressure[t] == b.fuel_inlet.pressure[t]

    m.fs.anode = GibbsReactor(
        default={"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "inert_species": ["N2", "Ar"],
                 "property_package": m.fs.syn_props})

    m.fs.anode_recycle = Separator(
        default={"outlet_list": ["exhaust", "recycle"],
                 "property_package": m.fs.syn_props})

    m.fs.anode_blower = PressureChanger(
        default={'compressor': True,
                 'property_package': m.fs.syn_props,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic})

    m.fs.recycle_translator = Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": m.fs.syn_props,
                 "outlet_property_package": m.fs.NG_props})

    # recycle is fed back to a part of the system that contains higher
    # hydrocarbons, so property package must be translated
    @m.fs.recycle_translator.Constraint(m.fs.time)
    def recycle_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.recycle_translator.Constraint(m.fs.time)
    def recycle_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.recycle_translator.Constraint(m.fs.time)
    def recycle_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.recycle_translator.Constraint(m.fs.time,
                                        m.fs.syn_props.component_list)
    def recycle_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    m.fs.recycle_translator.outlet.mole_frac_comp[0, 'C2H6'].fix(0)
    m.fs.recycle_translator.outlet.mole_frac_comp[0, 'C3H8'].fix(0)
    m.fs.recycle_translator.outlet.mole_frac_comp[0, 'C4H10'].fix(0)

    # build cathode side units
    m.fs.air_blower = PressureChanger(
        default={'compressor': True,
                 'property_package': m.fs.air_props,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic})

    m.fs.cathode_hx = HeatExchanger(
        default={"delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.air_props,
                           "has_pressure_change": True},
                 "tube": {"property_package": m.fs.air_props,
                          "has_pressure_change": True}})

    m.fs.cathode_mix = Mixer(
        default={"inlet_list": ["feed", "recycle"],
                 "property_package": m.fs.air_props})

    # represents oxygen moving across SOFC electrolyte
    m.fs.cathode = Separator(
        default={"outlet_list": ["air_outlet", "ion_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.air_props})

    m.fs.cathode.split_fraction[0, 'ion_outlet', 'H2O'].fix(0)
    m.fs.cathode.split_fraction[0, 'ion_outlet', 'CO2'].fix(0)
    m.fs.cathode.split_fraction[0, 'ion_outlet', 'N2'].fix(0)
    m.fs.cathode.split_fraction[0, 'ion_outlet', 'Ar'].fix(0)

    m.fs.cathode_translator = Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": m.fs.air_props,
                 "outlet_property_package": m.fs.syn_props})

    # cathode side O2 crosses the electrolyte to anode side
    @m.fs.cathode_translator.Constraint(m.fs.time)
    def cathode_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.cathode_translator.Constraint(m.fs.time)
    def cathode_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.cathode_translator.Constraint(m.fs.time)
    def cathode_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    m.fs.cathode_translator.outlet.mole_frac_comp.fix(0)
    m.fs.cathode_translator.outlet.mole_frac_comp[0, 'O2'].fix(1)

    m.fs.cathode_heat = Heater(
        default={"has_pressure_change": True,
                 "property_package": m.fs.air_props})

    m.fs.cathode_recycle = Separator(
        default={"outlet_list": ["exhaust", "recycle"],
                 "property_package": m.fs.air_props})

    m.fs.cathode_blower = PressureChanger(
        default={'compressor': True,
                 'property_package': m.fs.air_props,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic})

    # build combustor and HRSG units
    m.fs.cathode_exhaust_split = Separator(
        default={"outlet_list": ["exhaust_outlet", "combustor_outlet"],
                 "property_package": m.fs.air_props})

    m.fs.cathode_expander = PressureChanger(
        default={"compressor": False,
                 "property_package": m.fs.air_props,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic})

    m.fs.cathode_HRSG = Heater(
        default={"has_pressure_change": True,
                 "property_package": m.fs.air_props})

    m.fs.cathode_exhaust_translator = Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": m.fs.air_props,
                 "outlet_property_package": m.fs.syn_props})

    # for the portion of the cathode exhaust entering the combustor
    @m.fs.cathode_exhaust_translator.Constraint(m.fs.time)
    def cathode_exhaust_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.cathode_exhaust_translator.Constraint(m.fs.time)
    def cathode_exhaust_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.cathode_exhaust_translator.Constraint(m.fs.time)
    def cathode_exhaust_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.cathode_exhaust_translator.Constraint(m.fs.time,
                                                m.fs.air_props.component_list)
    def cathode_exhaust_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    for j in m.fs.syn_props.component_list:
        if j not in m.fs.air_props.component_list:
            m.fs.cathode_exhaust_translator.outlet.mole_frac_comp[0, j].fix(0)

    m.fs.combustor_mix = Mixer(
        default={"inlet_list": ["anode_inlet", "cathode_inlet"],
                 "property_package": m.fs.syn_props})

    m.fs.combustor = StoichiometricReactor(
        default={"has_heat_of_reaction": False,
                 "has_heat_transfer": False,
                 "has_pressure_change": True,
                 "property_package": m.fs.syn_props,
                 "reaction_package": m.fs.rxn_props})

    m.fs.combustor_expander = PressureChanger(
        default={"compressor": False,
                 "property_package": m.fs.syn_props,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic})

    m.fs.anode_HRSG = Heater(
        default={"has_pressure_change": True,
                 "property_package": m.fs.syn_props})

    # build arcs to connect unit operations
    # arcs for anode side
    m.fs.ANODE_MIXER = Arc(
        source=m.fs.anode_mix.outlet, destination=m.fs.anode_hx.tube_inlet)

    m.fs.PREREF_IN = Arc(
        source=m.fs.anode_hx.tube_outlet, destination=m.fs.prereformer.inlet)

    m.fs.PREREF_OUT = Arc(
        source=m.fs.prereformer.outlet,
        destination=m.fs.anode_translator.inlet)

    m.fs.ANODE_TRANS_OUT = Arc(
        source=m.fs.anode_translator.outlet,
        destination=m.fs.fuel_cell_mix.fuel_inlet)

    m.fs.ANODE_IN = Arc(
        source=m.fs.fuel_cell_mix.outlet, destination=m.fs.anode.inlet)

    m.fs.ANODE = Arc(
        source=m.fs.anode.outlet, destination=m.fs.anode_recycle.inlet)

    m.fs.ANODE_RECYCLE_HX = Arc(
        source=m.fs.anode_recycle.exhaust,
        destination=m.fs.anode_hx.shell_inlet)

    m.fs.ANODE_RECYCLE = Arc(
        source=m.fs.anode_recycle.recycle, destination=m.fs.anode_blower.inlet)

    m.fs.ANODE_BLOWER = Arc(
        source=m.fs.anode_blower.outlet,
        destination=m.fs.recycle_translator.inlet)

    m.fs.REC_TRANS_OUT = Arc(
        source=m.fs.recycle_translator.outlet,
        destination=m.fs.anode_mix.recycle)

    # arcs for cathode side
    m.fs.AIR_BLOWER = Arc(
        source=m.fs.air_blower.outlet, destination=m.fs.cathode_hx.tube_inlet)

    m.fs.CATHODE_HX_COLD = Arc(
        source=m.fs.cathode_hx.tube_outlet, destination=m.fs.cathode_mix.feed)

    m.fs.CATHODE_MIXER = Arc(
        source=m.fs.cathode_mix.outlet, destination=m.fs.cathode.inlet)

    m.fs.O2_IONS = Arc(
        source=m.fs.cathode.ion_outlet,
        destination=m.fs.cathode_translator.inlet)

    m.fs.O2_IONS_TRANSLATE = Arc(
        source=m.fs.cathode_translator.outlet,
        destination=m.fs.fuel_cell_mix.ion_inlet)

    m.fs.CATHODE = Arc(
        source=m.fs.cathode.air_outlet, destination=m.fs.cathode_heat.inlet)

    m.fs.CATHODE_HEAT = Arc(
        source=m.fs.cathode_heat.outlet,
        destination=m.fs.cathode_recycle.inlet)

    m.fs.CATHODE_RECYCLE_HX = Arc(
        source=m.fs.cathode_recycle.exhaust,
        destination=m.fs.cathode_hx.shell_inlet)

    m.fs.CATHODE_RECYCLE = Arc(
        source=m.fs.cathode_recycle.recycle,
        destination=m.fs.cathode_blower.inlet)

    m.fs.CATHODE_BLOWER = Arc(
        source=m.fs.cathode_blower.outlet,
        destination=m.fs.cathode_mix.recycle)

    # arcs for combustor and HRSG
    m.fs.CATHODE_HX_HOT_OUT = Arc(
        source=m.fs.cathode_hx.shell_outlet,
        destination=m.fs.cathode_exhaust_split.inlet)

    m.fs.CATHODE_EXPANDER_IN = Arc(
        source=m.fs.cathode_exhaust_split.exhaust_outlet,
        destination=m.fs.cathode_expander.inlet)

    m.fs.CATHODE_HRSG_IN = Arc(
        source=m.fs.cathode_expander.outlet,
        destination=m.fs.cathode_HRSG.inlet)

    m.fs.CATH_EXH_TRANS_IN = Arc(
        source=m.fs.cathode_exhaust_split.combustor_outlet,
        destination=m.fs.cathode_exhaust_translator.inlet)

    m.fs.CATH_EXH_TRANS_OUT = Arc(
        source=m.fs.cathode_exhaust_translator.outlet,
        destination=m.fs.combustor_mix.cathode_inlet)

    m.fs.ANODE_HX_HOT_OUT = Arc(
        source=m.fs.anode_hx.shell_outlet,
        destination=m.fs.combustor_mix.anode_inlet)

    m.fs.COMBUST_IN = Arc(
        source=m.fs.combustor_mix.outlet,
        destination=m.fs.combustor.inlet)

    m.fs.COMBUST_EXPANDER_IN = Arc(
        source=m.fs.combustor.outlet,
        destination=m.fs.combustor_expander.inlet)

    m.fs.ANODE_HRSG_IN = Arc(
        source=m.fs.combustor_expander.outlet,
        destination=m.fs.anode_HRSG.inlet)

    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)


def set_power_island_inputs(m):
    ###################
    # anode side inputs
    ###################

    # syngas feed conditions
    m.fs.anode_mix.feed.flow_mol.fix(3639)  # mol/s
    m.fs.anode_mix.feed.temperature.fix(621.3)  # K
    m.fs.anode_mix.feed.pressure.fix(137895)  # Pa, equal to 20 psia
    m.fs.anode_mix.feed.mole_frac_comp[0, 'CH4'].fix(0.1787)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'C2H6'].fix(0.0061)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'C3H8'].fix(0.0013)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'C4H10'].fix(0.0008)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'CO'].fix(0.0909)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'CO2'].fix(0.0438)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'H2'].fix(0.2752)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'H2O'].fix(0.1118)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'N2'].fix(0.2879)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'O2'].fix(0.0000)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'Ar'].fix(0.0034)

    # anode heat exchanger
    m.fs.anode_hx.tube.deltaP.fix(-137.9)  # equal to -0.02 psi
    m.fs.anode_hx.shell.deltaP.fix(-137.9)  # equal to -0.02 psi
    m.fs.anode_hx.area.fix(12664)  # m2
    m.fs.anode_hx.overall_heat_transfer_coefficient.fix(80)  # W/m2K

    # prereformer and anode
    m.fs.prereformer.heat_duty.fix(0)
    m.fs.prereformer.deltaP.fix(-137.9)  # equal to -0.02 psi

    m.fs.anode.outlet.temperature.fix(979.93)  # K, does not stay fixed
    m.fs.anode.outlet.pressure.fix(137137)  # Pa, equal to 19.89 psia

    # anode recycle and blower
    m.fs.anode_recycle.split_fraction[0, 'recycle'].fix(0.445)

    m.fs.anode_blower.outlet.pressure.fix(137888)  # Pa, equal to 20 psia
    m.fs.anode_blower.efficiency_isentropic.fix(0.8)

    #####################
    # cathode side inputs
    #####################

    # air feed conditions
    m.fs.air_blower.inlet.flow_mol.fix(17915.5)  # mol/s
    m.fs.air_blower.inlet.temperature.fix(288.15)  # K, equal to 59 F
    m.fs.air_blower.inlet.pressure.fix(101353)  # Pa, equal to 14.7 psia
    m.fs.air_blower.inlet.mole_frac_comp[0, 'H2O'].fix(0.0104)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'CO2'].fix(0.0003)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'N2'].fix(0.7722)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'O2'].fix(0.2077)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'Ar'].fix(0.0094)

    # air blower
    m.fs.air_blower.outlet.pressure.fix(111006)  # equal to 16.1 psi
    m.fs.air_blower.efficiency_isentropic.fix(0.82)

    # cathode heat exchanger
    m.fs.cathode_hx.tube_outlet.pressure.fix(105490)  # equal to 15.3 psia
    m.fs.cathode_hx.shell.deltaP.fix(-1379)  # equal to -0.2 psi
    m.fs.cathode_hx.area.fix(17159)  # m2
    m.fs.cathode_hx.overall_heat_transfer_coefficient.fix(80)  # W/m2K

    # cathode
    m.fs.cathode.ion_outlet.flow_mol.fix(1670)

    m.fs.cathode_heat.outlet.temperature.fix(999.4)  # K
    m.fs.cathode_heat.outlet.pressure.fix(104111)  # Pa, equal to 15.1 psi

    # cathode recycle and blower
    m.fs.cathode_recycle.split_fraction[0, 'recycle'].fix(0.5)

    m.fs.cathode_blower.outlet.pressure.fix(105490)  # Pa, equal to 15.3 psi
    m.fs.cathode_blower.efficiency_isentropic.fix(0.8)

    #####################
    # HRSG section inputs
    #####################

    m.fs.cathode_exhaust_split.split_fraction[0, 'combustor_outlet'].fix(0.3)

    m.fs.cathode_expander.deltaP.fix(-6.9)  # Pa, equal to -0.001 psi
    m.fs.cathode_expander.efficiency_isentropic.fix(0.95)

    m.fs.cathode_HRSG.outlet.temperature.fix(405.7)  # K, equal to 270 F
    m.fs.cathode_HRSG.deltaP.fix(-1379)  # Pa, equal to -0.2 psi

    m.fs.combustor.deltaP.fix(-6895)  # equal to -1 psi
    m.fs.combustor.outlet.mole_frac_comp[0, "H2"].fix(0)
    m.fs.combustor.outlet.mole_frac_comp[0, "CO"].fix(0)
    m.fs.combustor.outlet.mole_frac_comp[0, "CH4"].fix(0)

    m.fs.combustor_expander.deltaP.fix(-6.9)  # Pa, equal to -0.001 psi
    m.fs.combustor_expander.efficiency_isentropic.fix(0.95)

    m.fs.anode_HRSG.outlet.temperature.fix(405.7)  # K, equal to 270 F
    m.fs.anode_HRSG.deltaP.fix(-1379)  # Pa, equal to -0.2 psi


def build_reformer(m):
    # build unit models
    m.fs.reformer_recuperator = HeatExchanger(
        default={"delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.NG_props},
                 "tube": {"property_package": m.fs.NG_props}})

    m.fs.NG_expander = PressureChanger(
        default={'compressor': False,
                 'property_package': m.fs.NG_props,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic})

    m.fs.reformer_bypass = Separator(
        default={"outlet_list": ["reformer_outlet", "bypass_outlet"],
                 "property_package": m.fs.NG_props})

    m.fs.air_compressor_s1 = PressureChanger(
        default={"compressor": True,
                 "property_package": m.fs.NG_props,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    # sets air flow to reformer based on NG flow
    @m.fs.Constraint()
    def reformer_air_rule(fs):
        ng_flow = m.fs.reformer_recuperator.tube_inlet.flow_mol[0]
        air_flow = m.fs.air_compressor_s1.inlet.flow_mol[0]
        IR = m.fs.reformer_bypass.split_fraction[0, 'bypass_outlet']
        return air_flow == 2.868*(1 - IR)*ng_flow

    m.fs.intercooler_s1 = Heater(
        default={"property_package": m.fs.NG_props,
                 "has_pressure_change": True})

    m.fs.air_compressor_s2 = PressureChanger(
        default={"compressor": True,
                 "property_package": m.fs.NG_props,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    m.fs.intercooler_s2 = Heater(
        default={"property_package": m.fs.NG_props,
                 "has_pressure_change": True})

    m.fs.reformer_mix = Mixer(
        default={"inlet_list": ["gas_inlet", "oxygen_inlet", "steam_inlet"],
                 "property_package": m.fs.NG_props})

    # sets steam flow to reformer based on NG flow
    @m.fs.Constraint()
    def reformer_steam_flow(fs):
        ng_flow = m.fs.reformer_recuperator.tube_inlet.flow_mol[0]
        steam_flow = m.fs.reformer_mix.steam_inlet.flow_mol[0]
        IR = m.fs.reformer_bypass.split_fraction[0, 'bypass_outlet']
        return steam_flow == (1 - IR)*ng_flow

    m.fs.reformer = GibbsReactor(
        default={"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "inert_species": ["N2", "Ar"],
                 "property_package": m.fs.NG_props})

    m.fs.bypass_rejoin = Mixer(
        default={"inlet_list": ["syngas_inlet", "bypass_inlet"],
                 "property_package": m.fs.NG_props})

    # build and connect arcs
    m.fs.TO_NG_EXP = Arc(
        source=m.fs.reformer_recuperator.tube_outlet,
        destination=m.fs.NG_expander.inlet)

    m.fs.NG_EXP_OUT = Arc(
        source=m.fs.NG_expander.outlet, destination=m.fs.reformer_bypass.inlet)

    m.fs.TO_REF = Arc(
        source=m.fs.reformer_bypass.reformer_outlet,
        destination=m.fs.reformer_mix.gas_inlet)

    m.fs.STAGE_1_OUT = Arc(
        source=m.fs.air_compressor_s1.outlet,
        destination=m.fs.intercooler_s1.inlet)

    m.fs.IC_1_OUT = Arc(
        source=m.fs.intercooler_s1.outlet,
        destination=m.fs.air_compressor_s2.inlet)

    m.fs.STAGE_2_OUT = Arc(
        source=m.fs.air_compressor_s2.outlet,
        destination=m.fs.intercooler_s2.inlet)

    m.fs.IC_2_OUT = Arc(
        source=m.fs.intercooler_s2.outlet,
        destination=m.fs.reformer_mix.oxygen_inlet)

    m.fs.REF_IN = Arc(
        source=m.fs.reformer_mix.outlet, destination=m.fs.reformer.inlet)

    m.fs.REF_OUT = Arc(
        source=m.fs.reformer.outlet,
        destination=m.fs.reformer_recuperator.shell_inlet)

    m.fs.REF_RECUP_OUT = Arc(
        source=m.fs.reformer_recuperator.shell_outlet,
        destination=m.fs.bypass_rejoin.syngas_inlet)

    m.fs.REF_BYPASS = Arc(
        source=m.fs.reformer_bypass.bypass_outlet,
        destination=m.fs.bypass_rejoin.bypass_inlet)

    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)


def set_reformer_inputs(m):
    # natural gas feed conditions
    m.fs.reformer_recuperator.tube_inlet.flow_mol.fix(1161.9)  # mol/s
    m.fs.reformer_recuperator.tube_inlet.temperature.fix(288.15)  # K
    m.fs.reformer_recuperator.tube_inlet.pressure.fix(
        3447379)  # Pa, equal to 500 psia
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'CH4'].fix(0.931)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'C2H6'].fix(0.032)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'C3H8'].fix(0.007)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'C4H10'].fix(0.004)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'CO'].fix(0)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'CO2'].fix(0.01)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'H2'].fix(0)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'H2O'].fix(0)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'N2'].fix(0.016)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'O2'].fix(0)
    m.fs.reformer_recuperator.tube_inlet.mole_frac_comp[0, 'Ar'].fix(0)

    # recuperator conditions
    m.fs.reformer_recuperator.tube_outlet.temperature.fix(1010.93)
    m.fs.reformer_recuperator.overall_heat_transfer_coefficient.fix(80)

    # natural gas expander conditions
    m.fs.NG_expander.outlet.pressure.fix(206843)  # equal to 30 psia
    m.fs.NG_expander.efficiency_isentropic.fix(0.9)

    # internal reformation percentage
    m.fs.reformer_bypass.split_fraction[0, 'bypass_outlet'].fix(0.6)

    # air to reformer
    m.fs.air_compressor_s1.inlet.flow_mol[0] == 1332.9  # mol/s
    m.fs.air_compressor_s1.inlet.temperature.fix(288.15)  # K
    m.fs.air_compressor_s1.inlet.pressure.fix(101353)  # Pa, equal to 14.7 psia
    m.fs.air_compressor_s1.inlet.mole_frac_comp.fix(0)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'CO2'].fix(0.0003)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'H2O'].fix(0.0104)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'N2'].fix(0.7722)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'O2'].fix(0.2077)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'Ar'].fix(0.0094)

    # air compressors and intercoolers
    m.fs.air_compressor_s1.outlet.pressure.fix(144790)  # Pa, equal to 21 psia
    m.fs.air_compressor_s1.efficiency_isentropic.fix(0.84)

    m.fs.intercooler_s1.outlet.temperature.fix(310.93)  # K, equal to 100 F
    m.fs.intercooler_s1.deltaP.fix(-3447)  # Pa, equal to -0.5 psi

    m.fs.air_compressor_s2.outlet.pressure.fix(206843)  # Pa, equal to 30 psia
    m.fs.air_compressor_s2.efficiency_isentropic.fix(0.84)

    m.fs.intercooler_s2.outlet.temperature.fix(310.93)  # K, equal to 100 F
    m.fs.intercooler_s2.deltaP.fix(-3447)  # Pa, equal to -0.5 psi

    # steam to reformer
    m.fs.reformer_mix.steam_inlet.flow_mol[0] == 464.77  # mol/s
    m.fs.reformer_mix.steam_inlet.temperature.fix(422)  # K
    m.fs.reformer_mix.steam_inlet.pressure.fix(206843)  # Pa, equal to 30 psia
    m.fs.reformer_mix.steam_inlet.mole_frac_comp.fix(0)
    m.fs.reformer_mix.steam_inlet.mole_frac_comp[0, 'H2O'].fix(1)

    # reformer outlet pressure
    m.fs.reformer.outlet.pressure.fix(137895)  # Pa, equal to 20 Psi
    m.fs.reformer.outlet.temperature.fix(1060.93)  # K, equal to 1450 F


def scale_flowsheet(m):
    # set NG_props default scaling
    m.fs.NG_props.set_default_scaling("flow_mol", 1e-3)
    m.fs.NG_props.set_default_scaling("flow_mol_phase", 1e-3)
    m.fs.NG_props.set_default_scaling("temperature", 1e-2)
    m.fs.NG_props.set_default_scaling("pressure", 1e-5)

    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e1)
    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e2, index="C2H6")
    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e2, index="C3H8")
    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e2, index="C4H10")

    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e1)
    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e2,
                                      index=('Vap', 'C2H6'))
    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e2,
                                      index=('Vap', 'C3H8'))
    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e2,
                                      index=('Vap', 'C4H10'))

    m.fs.NG_props.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.NG_props.set_default_scaling("entr_mol_phase", 1e-1)

    # set syn_props default scaling
    m.fs.syn_props.set_default_scaling("flow_mol", 1e-3)
    m.fs.syn_props.set_default_scaling("flow_mol_phase", 1e-3)
    m.fs.syn_props.set_default_scaling("temperature", 1e-2)
    m.fs.syn_props.set_default_scaling("pressure", 1e-5)
    m.fs.syn_props.set_default_scaling("mole_frac_comp", 1e1)
    m.fs.syn_props.set_default_scaling("mole_frac_phase_comp", 1e1)
    m.fs.syn_props.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.syn_props.set_default_scaling("entr_mol_phase", 1e-1)

    # set air_props default scaling
    m.fs.air_props.set_default_scaling("flow_mol", 1e-3)
    m.fs.air_props.set_default_scaling("flow_mol_phase", 1e-3)
    m.fs.air_props.set_default_scaling("temperature", 1e-2)
    m.fs.air_props.set_default_scaling("pressure", 1e-5)
    m.fs.air_props.set_default_scaling("mole_frac_comp", 1e1)
    m.fs.air_props.set_default_scaling("mole_frac_phase_comp", 1e1)
    m.fs.air_props.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.air_props.set_default_scaling("entr_mol_phase", 1e-3)

    iscale.set_scaling_factor(m.fs.prereformer.lagrange_mult, 1e-4)
    iscale.set_scaling_factor(m.fs.anode.lagrange_mult, 1e-4)

    iscale.calculate_scaling_factors(m)
    # apply scaling factors
    m.fs.NG_props.set_default_scaling("flow_mol", 1e-3)
    m.fs.NG_props.set_default_scaling("flow_mol_phase", 1e-3)
    m.fs.NG_props.set_default_scaling("temperature", 1e-2)
    m.fs.NG_props.set_default_scaling("pressure", 1e-5)

    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e1)
    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e2, index="C2H6")
    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e2, index="C3H8")
    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e2, index="C4H10")

    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e1)
    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e2,
                                      index=('Vap', 'C2H6'))
    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e2,
                                      index=('Vap', 'C3H8'))
    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e2,
                                      index=('Vap', 'C4H10'))

    m.fs.NG_props.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.NG_props.set_default_scaling("entr_mol_phase", 1e-1)

    iscale.set_scaling_factor(m.fs.reformer.lagrange_mult, 1e-4)

    iscale.calculate_scaling_factors(m)


def initialize_power_island(m):
    # cathode side
    m.fs.air_blower.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.cathode_hx.tube_inlet, m.fs.air_blower.outlet)

    # fix cathode inlet to initial guess
    m.fs.cathode.inlet.flow_mol[0] = 34174
    m.fs.cathode.inlet.temperature[0] = 892
    m.fs.cathode.inlet.pressure[0] = 105490
    m.fs.cathode.inlet.mole_frac_comp[0, 'H2O'] = 0.0109
    m.fs.cathode.inlet.mole_frac_comp[0, 'CO2'] = 0.0003
    m.fs.cathode.inlet.mole_frac_comp[0, 'N2'] = 0.8099
    m.fs.cathode.inlet.mole_frac_comp[0, 'O2'] = 0.1690
    m.fs.cathode.inlet.mole_frac_comp[0, 'Ar'] = 0.0099

    m.fs.cathode.initialize(outlvl=logging.INFO)

    m.fs.cathode.inlet.unfix()

    # cathode translator block
    copy_port_values(
        m.fs.cathode_translator.inlet, m.fs.cathode.ion_outlet)

    m.fs.cathode_translator.initialize()

    # rest of cathode side
    copy_port_values(
        m.fs.cathode_heat.inlet, m.fs.cathode.air_outlet)

    m.fs.cathode_heat.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.cathode_recycle.inlet, m.fs.cathode_heat.outlet)

    m.fs.cathode_recycle.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.cathode_hx.shell_inlet, m.fs.cathode_recycle.exhaust)

    m.fs.cathode_hx.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.cathode_blower.inlet, m.fs.cathode_recycle.recycle)

    m.fs.cathode_blower.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.cathode_mix.recycle, m.fs.cathode_blower.outlet)

    copy_port_values(
        m.fs.cathode_mix.feed, m.fs.cathode_hx.tube_outlet)

    m.fs.cathode_mix.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.cathode.inlet, m.fs.cathode_mix.outlet)

    # anode side
    # anode inlet is used as tear stream
    m.fs.anode.inlet.flow_mol[0] = 9702
    m.fs.anode.inlet.temperature[0] = 844
    m.fs.anode.inlet.pressure[0] = 137612
    m.fs.anode.inlet.mole_frac_comp[0, 'CH4'] = 0.058
    m.fs.anode.inlet.mole_frac_comp[0, 'CO'] = 0.044
    m.fs.anode.inlet.mole_frac_comp[0, 'CO2'] = 0.123
    m.fs.anode.inlet.mole_frac_comp[0, 'H2'] = 0.223
    m.fs.anode.inlet.mole_frac_comp[0, 'H2O'] = 0.183
    m.fs.anode.inlet.mole_frac_comp[0, 'N2'] = 0.195
    m.fs.anode.inlet.mole_frac_comp[0, 'O2'] = 0.172
    m.fs.anode.inlet.mole_frac_comp[0, 'Ar'] = 0.002

    m.fs.anode.lagrange_mult[0, "C"] = 51092
    m.fs.anode.lagrange_mult[0, "H"] = 78296
    m.fs.anode.lagrange_mult[0, "O"] = 291784

    m.fs.anode.outlet.mole_frac_comp[0, "O2"] = 0

    m.fs.anode.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.anode_recycle.inlet, m.fs.anode.outlet)

    m.fs.anode_recycle.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.anode_blower.inlet, m.fs.anode_recycle.recycle)

    m.fs.anode_blower.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.recycle_translator.inlet, m.fs.anode_blower.outlet)

    m.fs.recycle_translator.initialize()

    copy_port_values(
        m.fs.anode_mix.recycle, m.fs.recycle_translator.outlet)

    m.fs.anode_mix.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.anode_hx.tube_inlet, m.fs.anode_mix.outlet)

    copy_port_values(
        m.fs.anode_hx.shell_inlet, m.fs.anode_recycle.exhaust)

    m.fs.anode_hx.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.prereformer.inlet, m.fs.anode_hx.tube_outlet)

    copy_port_values(
        m.fs.prereformer.outlet, m.fs.prereformer.inlet)

    m.fs.prereformer.gibbs_scaling = 1e-4

    m.fs.prereformer.lagrange_mult[0, "C"] = 9707
    m.fs.prereformer.lagrange_mult[0, "H"] = 62744
    m.fs.prereformer.lagrange_mult[0, "O"] = 293569

    m.fs.prereformer.outlet.mole_frac_comp[0, "O2"] = 0
    m.fs.prereformer.outlet.mole_frac_comp[0, "Ar"] = 0.003
    m.fs.prereformer.outlet.mole_frac_comp[0, "C2H6"] = 0
    m.fs.prereformer.outlet.mole_frac_comp[0, "C3H8"] = 0
    m.fs.prereformer.outlet.mole_frac_comp[0, "C4H10"] = 0

    m.fs.prereformer.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.anode_translator.inlet, m.fs.prereformer.outlet)

    m.fs.anode_translator.initialize()

    copy_port_values(
        m.fs.fuel_cell_mix.fuel_inlet, m.fs.anode_translator.outlet)

    copy_port_values(
        m.fs.fuel_cell_mix.ion_inlet, m.fs.cathode_translator.outlet)

    m.fs.fuel_cell_mix.initialize(outlvl=logging.INFO)

    ##############################
    # Combustor and HRSG section #
    ##############################

    # cathode side
    copy_port_values(m.fs.cathode_exhaust_split.inlet,
                     m.fs.cathode_hx.shell_outlet)

    m.fs.cathode_exhaust_split.initialize(outlvl=logging.INFO)

    copy_port_values(m.fs.cathode_expander.inlet,
                     m.fs.cathode_exhaust_split.exhaust_outlet)

    m.fs.cathode_expander.initialize(outlvl=logging.INFO)

    copy_port_values(m.fs.cathode_HRSG.inlet,
                     m.fs.cathode_expander.outlet)

    m.fs.cathode_HRSG.initialize(outlvl=logging.INFO)

    copy_port_values(m.fs.cathode_exhaust_translator.inlet,
                     m.fs.cathode_exhaust_split.combustor_outlet)

    m.fs.cathode_exhaust_translator.initialize()

    # anode side
    copy_port_values(m.fs.combustor_mix.cathode_inlet,
                     m.fs.cathode_exhaust_translator.outlet)

    copy_port_values(m.fs.combustor_mix.anode_inlet,
                     m.fs.anode_hx.shell_outlet)

    m.fs.combustor_mix.initialize(outlvl=logging.INFO)

    copy_port_values(m.fs.combustor.inlet,
                     m.fs.combustor_mix.outlet)

    m.fs.combustor.initialize()

    copy_port_values(m.fs.combustor_expander.inlet,
                     m.fs.combustor.outlet)

    m.fs.combustor_expander.initialize()

    copy_port_values(m.fs.anode_HRSG.inlet,
                     m.fs.combustor_expander.outlet)

    m.fs.anode_HRSG.initialize()


def initialize_reformer(m):
    m.fs.reformer.inlet.flow_mol[0] = 2262  # mol/s
    m.fs.reformer.inlet.temperature[0] = 470  # K
    m.fs.reformer.inlet.pressure[0] = 203395  # Pa
    m.fs.reformer.inlet.mole_frac_comp[0, 'CH4'] = 0.191
    m.fs.reformer.inlet.mole_frac_comp[0, 'C2H6'] = 0.006
    m.fs.reformer.inlet.mole_frac_comp[0, 'C3H8'] = 0.002
    m.fs.reformer.inlet.mole_frac_comp[0, 'C4H10'] = 0.001
    m.fs.reformer.inlet.mole_frac_comp[0, 'H2'] = 0
    m.fs.reformer.inlet.mole_frac_comp[0, 'CO'] = 0
    m.fs.reformer.inlet.mole_frac_comp[0, 'CO2'] = 0.002
    m.fs.reformer.inlet.mole_frac_comp[0, 'H2O'] = 0.212
    m.fs.reformer.inlet.mole_frac_comp[0, 'N2'] = 0.458
    m.fs.reformer.inlet.mole_frac_comp[0, 'O2'] = 0.122
    m.fs.reformer.inlet.mole_frac_comp[0, 'Ar'] = 0.006

    m.fs.reformer.lagrange_mult[0, 'C'] = 39230
    m.fs.reformer.lagrange_mult[0, 'H'] = 81252
    m.fs.reformer.lagrange_mult[0, 'O'] = 315049

    m.fs.reformer.outlet.mole_frac_comp[0, 'O2'] = 0
    m.fs.reformer.outlet.mole_frac_comp[0, 'Ar'] = 0.004
    m.fs.reformer.outlet.mole_frac_comp[0, 'CH4'] = 0.0005
    m.fs.reformer.outlet.mole_frac_comp[0, 'C2H6'] = 0
    m.fs.reformer.outlet.mole_frac_comp[0, 'C3H8'] = 0
    m.fs.reformer.outlet.mole_frac_comp[0, 'C4H10'] = 0

    m.fs.reformer.initialize()

    # reformer recuperator
    copy_port_values(
        m.fs.reformer_recuperator.shell_inlet, m.fs.reformer.outlet)

    m.fs.reformer_recuperator.initialize()

    # NG expander
    copy_port_values(
        m.fs.NG_expander.inlet, m.fs.reformer_recuperator.tube_outlet)

    m.fs.NG_expander.initialize()

    # reformer bypass
    copy_port_values(
        m.fs.reformer_bypass.inlet, m.fs.NG_expander.outlet)

    m.fs.reformer_bypass.initialize()

    # air compressor train
    m.fs.air_compressor_s1.initialize()

    copy_port_values(
        m.fs.intercooler_s1.inlet, m.fs.air_compressor_s1.outlet)

    m.fs.intercooler_s1.initialize()

    copy_port_values(
        m.fs.air_compressor_s2.inlet, m.fs.intercooler_s1.outlet)

    m.fs.air_compressor_s2.initialize()

    copy_port_values(
        m.fs.intercooler_s2.inlet, m.fs.air_compressor_s2.outlet)

    m.fs.intercooler_s2.initialize()

    # reformer mixer
    copy_port_values(
        m.fs.reformer_mix.oxygen_inlet, m.fs.intercooler_s2.outlet)

    copy_port_values(
        m.fs.reformer_mix.gas_inlet, m.fs.reformer_bypass.reformer_outlet)

    m.fs.reformer_mix.initialize()

    # bypass rejoin
    copy_port_values(
        m.fs.bypass_rejoin.syngas_inlet,
        m.fs.reformer_recuperator.shell_outlet)

    copy_port_values(
        m.fs.bypass_rejoin.bypass_inlet, m.fs.reformer_bypass.bypass_outlet)

    m.fs.bypass_rejoin.initialize()


def connect_reformer_to_power_island(m):
    m.fs.anode_mix.feed.unfix()
    m.fs.SYN_CONNECT = Arc(
        source=m.fs.bypass_rejoin.outlet, destination=m.fs.anode_mix.feed)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)


def SOFC_ROM_setup(m):
    # create the ROM
    build_SOFC_ROM(m.fs)

    # build constraints connecting flowsheet to ROM input vars

    @m.fs.Constraint()
    def ROM_fuel_inlet_temperature(fs):
        return(fs.SOFC.fuel_temperature ==
               fs.anode_mix.feed.temperature[0]/pyunits.K - 273)

    @m.fs.Constraint()
    def ROM_air_inlet_temperature(fs):
        return(fs.SOFC.air_temperature ==
               fs.cathode.inlet.temperature[0]/pyunits.K - 273)

    @m.fs.Constraint()
    def ROM_air_recirculation(fs):
        return(fs.SOFC.air_recirculation ==
               fs.cathode_recycle.split_fraction[0, 'recycle'])

    @m.fs.Constraint()
    def ROM_OTC(fs):
        O_frac = (1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO']
                  + 2 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO2']
                  + 1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'H2O'])

        C_frac = (1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO']
                  + 1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO2']
                  + 1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CH4'])

        return fs.SOFC.OTC == O_frac/C_frac

    @m.fs.Constraint()
    def ROM_fuel_utilization(fs):
        full_O2_flow = (fs.anode_mix.feed.flow_mol[0] * (
            0.5 * fs.anode_mix.feed.mole_frac_comp[0, 'H2'] +
            0.5 * fs.anode_mix.feed.mole_frac_comp[0, 'CO'] +
            2.0 * fs.anode_mix.feed.mole_frac_comp[0, 'CH4'] +
            3.5 * fs.anode_mix.feed.mole_frac_comp[0, 'C2H6'] +
            5.0 * fs.anode_mix.feed.mole_frac_comp[0, 'C3H8'] +
            6.5 * fs.anode_mix.feed.mole_frac_comp[0, 'C4H10']))

        return (fs.SOFC.fuel_util*full_O2_flow ==
                fs.cathode.ion_outlet.flow_mol[0])

    @m.fs.Constraint()
    def ROM_air_utilization(fs):
        air_in = (fs.cathode_hx.tube_inlet.flow_mol[0] *
                  fs.cathode_hx.tube_inlet.mole_frac_comp[0, 'O2'])

        air_out = (fs.cathode_hx.shell_inlet.flow_mol[0] *
                   fs.cathode_hx.shell_inlet.mole_frac_comp[0, 'O2'])

        return fs.SOFC.air_util == 1 - air_out/air_in

    @m.fs.Constraint()
    def ROM_internal_reformation(fs):
        return (fs.SOFC.internal_reforming ==
                fs.reformer_bypass.split_fraction[0, 'bypass_outlet'])

    # we want to analyze this flowsheet in terms of ROM inputs, so fix ROM
    # inputs and unfix flowsheet vars that control inputs

    # current density
    m.fs.SOFC.current_density.fix(4000)

    # fuel temperature
    m.fs.reformer_recuperator.tube_outlet.temperature.unfix()
    m.fs.SOFC.fuel_temperature.fix(348.33)

    # internal reformation fraction
    m.fs.reformer_bypass.split_fraction[0, 'bypass_outlet'].unfix()
    m.fs.SOFC.internal_reforming.fix(0.6)

    # air temperature - constrained by max cell temperature
    m.fs.cathode_hx.area.unfix()
    m.fs.SOFC.max_cell_temperature.fix(750)

    # air recirculation fraction
    m.fs.cathode_recycle.split_fraction[0, 'recycle'].unfix()
    m.fs.SOFC.air_recirculation.fix(0.5)

    # oxygen to carbon ratio
    m.fs.anode_recycle.split_fraction[0, 'recycle'].unfix()
    m.fs.SOFC.OTC.fix(2.1)

    # fuel utilization
    m.fs.cathode.ion_outlet.flow_mol.unfix()
    m.fs.SOFC.fuel_util.fix(0.8)

    # air utilization - constrained by deltaT across cell
    m.fs.air_blower.inlet.flow_mol.unfix()
    m.fs.SOFC.deltaT_cell.fix(100)

    # pressure
    m.fs.SOFC.pressure.fix(1)

    # initialize ROM
    calculate_variable_from_constraint(m.fs.SOFC.air_temperature,
                                       m.fs.ROM_air_inlet_temperature)
    calculate_variable_from_constraint(m.fs.SOFC.air_util,
                                       m.fs.ROM_air_utilization)
    initialize_SOFC_ROM(m.fs.SOFC)

    # add constraints for power calculations
    m.fs.F = pyo.Param(initialize=96487, units=pyunits.C/pyunits.mol)
    m.fs.stack_current = pyo.Var(initialize=650, units=pyunits.MA)
    m.fs.stack_power = pyo.Var(initialize=600, units=pyunits.MW)

    @m.fs.Constraint()
    def stack_current_constraint(fs):
        return (fs.stack_current ==
                pyunits.convert(4*fs.F*m.fs.cathode.ion_outlet.flow_mol[0],
                                pyunits.MA))

    @m.fs.Constraint()
    def stack_power_constraint(fs):
        return fs.stack_power == fs.stack_current*fs.SOFC.stack_voltage

    # initialize power calculations
    calculate_variable_from_constraint(m.fs.stack_current,
                                       m.fs.stack_current_constraint)
    calculate_variable_from_constraint(m.fs.stack_power,
                                       m.fs.stack_power_constraint)


def add_SOFC_energy_balance(m):
    @m.fs.Constraint()
    def ROM_anode_outlet_temperature(fs):
        return (fs.SOFC.anode_outlet_temperature ==
                fs.anode.outlet.temperature[0]/pyunits.K - 273)

    m.fs.anode.outlet.temperature.unfix()

    @m.fs.Constraint()
    def SOFC_energy_balance(fs):
        return (-1*pyunits.convert(fs.anode.heat_duty[0], pyunits.MW) ==
                fs.stack_power +
                pyunits.convert(fs.cathode_heat.heat_duty[0], pyunits.MW))

    m.fs.cathode_heat.outlet.temperature.unfix()


def add_result_constraints(m):
    # total heat supplied to HRSG
    m.fs.HRSG_heat_duty = pyo.Var(initialize=300, units=pyunits.MW)

    @m.fs.Constraint()
    def HRSG_heat_duty_constraint(fs):
        return (-1*fs.HRSG_heat_duty ==
                pyunits.convert(fs.anode_HRSG.heat_duty[0], pyunits.MW) +
                pyunits.convert(fs.cathode_HRSG.heat_duty[0], pyunits.MW))

    # heat required for generating reformer steam
    m.fs.reformer_steam_heat = pyo.Var(initialize=25, units=pyunits.MW)

    # only valid when steam is heated to 422 K (149 C)
    @m.fs.Constraint()
    def reformer_steam_heat_constraint(fs):
        return (fs.reformer_steam_heat == 0.050735*pyunits.MJ/pyunits.mol *
                fs.reformer_mix.steam_inlet.flow_mol[0])

    # HRSG heat applied toward steam cycle
    m.fs.steam_cycle_heat = pyo.Var(initialize=300, units=pyunits.MW)

    @m.fs.Constraint()
    def steam_cycle_heat_constraint(fs):
        return (fs.steam_cycle_heat == fs.HRSG_heat_duty -
                fs.reformer_steam_heat -
                pyunits.convert(fs.reformer.heat_duty[0], pyunits.MW))

    # power generated by steam cycle
    m.fs.steam_cycle_power = pyo.Var(initialize=100, units=pyunits.MW)
    m.fs.steam_cycle_efficiency = pyo.Param(initialize=0.38, mutable=True)

    @m.fs.Constraint()
    def steam_cycle_power_constraint(fs):
        return (fs.steam_cycle_power ==
                fs.steam_cycle_heat * fs.steam_cycle_efficiency)

    # stack AC power
    m.fs.stack_power_AC = pyo.Var(initialize=550, units=pyunits.MW)
    m.fs.inverter_efficiency = pyo.Param(initialize=0.97, mutable=True)

    @m.fs.Constraint()
    def stack_AC_power_constraint(fs):
        return fs.stack_power_AC == fs.stack_power * fs.inverter_efficiency

    # gross plant power
    m.fs.gross_power = pyo.Var(initialize=670, units=pyunits.MW)

    @m.fs.Constraint()
    def gross_power_constraint(fs):
        return (fs.gross_power == fs.stack_power_AC +
                fs.steam_cycle_power +
                -1*pyunits.convert(fs.NG_expander.work_mechanical[0],
                                   pyunits.MW))

    # auxiliary load of the plant
    m.fs.auxiliary_load = pyo.Var(initialize=10, units=pyunits.MW)

    @m.fs.Constraint()
    def auxiliary_load_constraint(fs):
        return (fs.auxiliary_load == pyunits.convert(
                   (fs.air_blower.work_mechanical[0] +
                    fs.cathode_blower.work_mechanical[0] +
                    fs.anode_blower.work_mechanical[0] +
                    fs.air_compressor_s1.work_mechanical[0] +
                    fs.air_compressor_s2.work_mechanical[0]), pyunits.MW))

    # net plant power
    m.fs.net_power = pyo.Var(initialize=660, units=pyunits.MW)

    @m.fs.Constraint()
    def net_power_constraint(fs):
        return fs.net_power == fs.gross_power - fs.auxiliary_load

    # HHV efficiency
    m.fs.HHV_efficiency = pyo.Var(initialize=0.6)

    @m.fs.Constraint()
    def efficiency_rule(fs):
        NG_HHV = 908839.23*pyunits.J/pyunits.mol
        return (fs.HHV_efficiency == fs.net_power /
                pyunits.convert(
                    (NG_HHV * fs.reformer_recuperator.tube_inlet.flow_mol[0]),
                    pyunits.MW))

    # CO2 emissions in g/kWh
    m.fs.CO2_emissions = pyo.Var(initialize=300,
                                 units=pyunits.g/pyunits.kWh)

    @m.fs.Constraint()
    def CO2_emission_constraint(fs):
        mass_flow = (fs.anode_HRSG.outlet.flow_mol[0] *
                     fs.anode_HRSG.outlet.mole_frac_comp[0, 'CO2'] *
                     44.01*pyunits.g/pyunits.mol)
        return fs.CO2_emissions == pyunits.convert((mass_flow/fs.net_power),
                                                   (pyunits.g/pyunits.kWh))


def make_stream_dict(m):
    """Adds _streams to m, which contains a dictionary of streams for display

    Args:
        m (ConcreteModel): A Pyomo model from create_model()

    Returns:
        None
    """

    m._streams = OrderedDict(
        [
            ("FUEL_IN", m.fs.reformer_recuperator.tube_inlet),
            ("HOT_FUEL", m.fs.reformer_recuperator.tube_outlet),
            ("REF_IN", m.fs.reformer_mix.gas_inlet),
            ("AIR_IN", m.fs.reformer_mix.oxygen_inlet),
            ("STEAM_IN", m.fs.reformer_mix.steam_inlet),
            ("REF_OUT", m.fs.reformer.outlet),
            ("SYN_IN", m.fs.anode_mix.feed),
            ("ANO_HX_CI", m.fs.ANODE_MIXER),
            ("ANODE_IN", m.fs.fuel_cell_mix.fuel_inlet),
            ("ANODE_OUT", m.fs.anode.outlet),
            ("ANODE_REC", m.fs.ANODE_BLOWER),
            ("ANO_HX_HI", m.fs.ANODE_RECYCLE_HX),
            ("ANO_HX_HO", m.fs.anode_hx.shell_outlet),
            ("AIR", m.fs.air_blower.inlet),
            ("CATH_HX_CI", m.fs.AIR_BLOWER),
            ("CATH_HX_CO", m.fs.CATHODE_HX_COLD),
            ("CATH_IN", m.fs.CATHODE_MIXER),
            ("CATH_OUT", m.fs.CATHODE_HEAT),
            ("CATH_REC", m.fs.CATHODE_BLOWER),
            ("CATH_HX_HI", m.fs.CATHODE_RECYCLE_HX),
            ("CATH_HX_HO", m.fs.cathode_hx.shell_outlet),
            ("CATH_HRSG_IN", m.fs.cathode_HRSG.inlet),
            ("CATH_EXH", m.fs.cathode_HRSG.outlet),
            ("COMB_AIR", m.fs.combustor_mix.cathode_inlet),
            ("COMB_OUT", m.fs.combustor.outlet),
            ("ANOD_EXH", m.fs.anode_HRSG.outlet)])


def pfd_result(outfile, m, df):
    # a function for string formatting
    def fstr(value, decimals, unit=''):
        if decimals == 0:
            rounded_value = int(value)
        else:
            rounded_value = round(value, decimals)
        return "{:,}".format(rounded_value) + unit
    # set up tags
    tags = {}
    for i in df.index:
        tags[i + "_F"] = fstr(df.loc[i, "Total Molar Flowrate"], 0, ' mol/s')
        tags[i + "_T"] = fstr(df.loc[i, "Temperature"], 0, ' K')
        tags[i + "_P"] = fstr(df.loc[i, "Pressure"]/1000, 0, ' kPa')
    tags["current_density"] = fstr(pyo.value(m.fs.SOFC.current_density), 0)
    tags["fuel_temp"] = fstr(pyo.value(m.fs.SOFC.fuel_temperature), 1)
    tags["int_ref"] = fstr(pyo.value(m.fs.SOFC.internal_reforming), 2)
    tags["air_temp"] = fstr(pyo.value(m.fs.SOFC.air_temperature), 1)
    tags["air_recirc"] = fstr(pyo.value(m.fs.SOFC.air_recirculation), 2)
    tags["OTC"] = fstr(pyo.value(m.fs.SOFC.OTC), 3)
    tags["fuel_util"] = fstr(pyo.value(m.fs.SOFC.fuel_util), 3)
    tags["air_util"] = fstr(pyo.value(m.fs.SOFC.air_util), 3)
    tags["voltage"] = fstr(pyo.value(m.fs.SOFC.stack_voltage), 4)
    tags["power"] = fstr(pyo.value(m.fs.stack_power), 1)
    tags["recuperator_duty"] = fstr(
        (pyo.value(m.fs.reformer_recuperator.heat_duty[0])/1e6), 1)
    tags["anode_hx_duty"] = fstr(
        (pyo.value(m.fs.anode_hx.heat_duty[0])/1e6), 2)
    tags["cathode_hx_duty"] = fstr(
        (pyo.value(m.fs.cathode_hx.heat_duty[0])/1e6), 2)
    tags["anode_duty"] = fstr((pyo.value(m.fs.anode.heat_duty[0])/1e6), 2)
    tags["cathode_duty"] = fstr(
        (pyo.value(m.fs.cathode_heat.heat_duty[0])/1e6), 2)
    tags["stack_power_AC"] = fstr(pyo.value(m.fs.stack_power_AC), 1)
    tags["HRSG_power"] = fstr(pyo.value(m.fs.steam_cycle_power), 1)
    tags["expander_power"] = fstr(-1*pyo.value(
        m.fs.NG_expander.work_mechanical[0])/1e6, 1)
    tags["aux_load"] = fstr(pyo.value(m.fs.auxiliary_load), 1)
    tags["net_power"] = fstr(pyo.value(m.fs.net_power), 1)
    tags["thermal_input"] = fstr(pyo.value(
        m.fs.net_power/m.fs.HHV_efficiency), 1)
    tags["efficiency"] = fstr(pyo.value(m.fs.HHV_efficiency)*100, 2)
    tags["emissions"] = fstr(pyo.value(m.fs.CO2_emissions), 1)

    original_svg_file = os.path.join(this_file_dir(),
                                     "NGFC_results_template.svg")
    with open(original_svg_file, "r") as f:
        svg_tag(tags, f, outfile=outfile)


def main():
    # create model and flowsheet
    m = pyo.ConcreteModel(name='NGFC without carbon capture')
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # solver and options
    solver = pyo.SolverFactory("ipopt")
    solver.options = {'bound_push': 1e-16}

    if os.path.exists('NGFC_flowsheet_init.json.gz'):
        build_power_island(m)
        build_reformer(m)
        connect_reformer_to_power_island(m)
        SOFC_ROM_setup(m)
        add_SOFC_energy_balance(m)
        add_result_constraints(m)
        scale_flowsheet(m)
        ms.from_json(m, fname='NGFC_flowsheet_init.json.gz')

    else:
        build_power_island(m)
        build_reformer(m)
        set_power_island_inputs(m)
        set_reformer_inputs(m)
        scale_flowsheet(m)
        initialize_power_island(m)
        initialize_reformer(m)
        connect_reformer_to_power_island(m)
        SOFC_ROM_setup(m)
        add_SOFC_energy_balance(m)
        add_result_constraints(m)
        solver.solve(m, tee=True)
        ms.to_json(m, fname='NGFC_flowsheet_init.json.gz')

    # uncomment to report results
    # make_stream_dict(m)
    # df = create_stream_table_dataframe(streams=m._streams, orient="index")
    # pfd_result("NGFC_results.svg", m, df)

    return m
