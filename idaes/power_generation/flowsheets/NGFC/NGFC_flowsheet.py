##############################################################################
# The development of this flowsheet/code is funded by the ARPA-E DIFFERENTIATE
# project: “Machine Learning for Natural Gas to Electric Power System Design”
# Project number: DE-FOA-0002107-1625.
# This project is a collaborative effort between the Pacific Northwest National
# Laboratory, the National Energy Technology Laboratory, and the University of
# Washington to design NGFC systems with high efficiencies and low CO2
# emissions.
##############################################################################
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
from pyomo.network import Arc
from pyomo.common.fileutils import this_file_dir
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# IDAES Imports
from idaes.core.util import copy_port_values
import idaes.core.util.scaling as iscale
from idaes.core.util.misc import svg_tag

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

import logging

# Custom imports
from properties.natural_gas_prop import get_SOFC_properties, rxn_configuration
from ROM.SOFC_ROM import make_SOFC_ROM

# %% SOFC Power Island Section


def build_power_island(m):
    NG_config = get_SOFC_properties(
        components=['H2', 'CO', 'CH4', "C2H6", "C3H8", "C4H10",
                    "H2O", 'CO2',
                    'N2', 'O2', 'Ar'])
    m.fs.NG_props = GenericParameterBlock(default=NG_config)
    # apply scaling factors
    m.fs.NG_props.set_default_scaling("flow_mol", 1e-3)
    m.fs.NG_props.set_default_scaling("flow_mol_phase", 1e-3, index='Vap')
    m.fs.NG_props.set_default_scaling("temperature", 1e-2)
    m.fs.NG_props.set_default_scaling("pressure", 1e-5)
    m.fs.NG_props.set_default_scaling("mole_frac_comp", 10)
    m.fs.NG_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'H2'))
    m.fs.NG_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'CO'))
    m.fs.NG_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'CH4'))
    m.fs.NG_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'H2O'))
    m.fs.NG_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'CO2'))
    m.fs.NG_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'N2'))
    m.fs.NG_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'O2'))
    m.fs.NG_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'Ar'))
    m.fs.NG_props.set_default_scaling("enth_mol_phase", 1e-3)

    syn_config = get_SOFC_properties(
        components=["H2", "CO", "CH4",
                    "H2O", "CO2",
                    "N2", "O2", "Ar"])
    m.fs.syn_props = GenericParameterBlock(default=syn_config)
    # apply scaling factors
    m.fs.syn_props.set_default_scaling("flow_mol", 1e-3)
    m.fs.syn_props.set_default_scaling("flow_mol_phase", 1e-3, index='Vap')
    m.fs.syn_props.set_default_scaling("temperature", 1e-2)
    m.fs.syn_props.set_default_scaling("pressure", 1e-5)
    m.fs.syn_props.set_default_scaling("mole_frac_comp", 10)
    m.fs.syn_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'H2'))
    m.fs.syn_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'CO'))
    m.fs.syn_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'CH4'))
    m.fs.syn_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'H2O'))
    m.fs.syn_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'CO2'))
    m.fs.syn_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'N2'))
    m.fs.syn_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'O2'))
    m.fs.syn_props.set_default_scaling(
        "mole_frac_phase_comp", 10, index=('Vap', 'Ar'))
    m.fs.syn_props.set_default_scaling("enth_mol_phase", 1e-3)

    m.fs.rxn_props = GenericReactionParameterBlock(
        default={"property_package": m.fs.syn_props, **rxn_configuration})

    air_config = get_SOFC_properties(
        components=['H2O', 'CO2', 'N2', 'O2', 'Ar'])
    m.fs.air_props = GenericParameterBlock(default=air_config)
    # apply scaling factors
    m.fs.air_props.set_default_scaling("flow_mol", 1e-3)
    m.fs.air_props.set_default_scaling("flow_mol_phase", 1e-3)
    m.fs.air_props.set_default_scaling("temperature", 1e-2)
    m.fs.air_props.set_default_scaling("pressure", 1e-5)
    m.fs.air_props.set_default_scaling("mole_frac_comp", 1e1)
    m.fs.air_props.set_default_scaling("mole_frac_phase_comp", 1e1)
    m.fs.air_props.set_default_scaling("enth_mol_phase", 1e-3)

    # anode side units
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

    m.fs.fuel_cell_mix = Mixer(
        default={"inlet_list": ["fuel_inlet", "ion_inlet"],
                 "momentum_mixing_type": MomentumMixingType.none,
                 "property_package": m.fs.syn_props})

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

    # cathode side units
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

    m.fs.cathode = Separator(
        default={"outlet_list": ["air_outlet", "ion_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.air_props})

    m.fs.cathode_translator = Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": m.fs.air_props,
                 "outlet_property_package": m.fs.syn_props})

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

    # combustor and HRSG units
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
    # anode side arcs
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

    # cathode side arcs
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

    # combustor and HRSG arcs
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

    # build translator block constraints
    # anode translator
    @m.fs.anode_translator.Constraint(m.fs.time)
    def temp_constraint(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.anode_translator.Constraint(m.fs.time)
    def pressure_constraint(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.anode_translator.Constraint(m.fs.time)
    def flow_constraint(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.anode_translator.Constraint(m.fs.time,
                                      m.fs.syn_props.component_list)
    def mole_frac_constraint(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    # anode recycle translator
    @m.fs.recycle_translator.Constraint(m.fs.time)
    def temp_constraint(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.recycle_translator.Constraint(m.fs.time)
    def pressure_constraint(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.recycle_translator.Constraint(m.fs.time)
    def flow_constraint(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.recycle_translator.Constraint(m.fs.time,
                                        m.fs.syn_props.component_list)
    def mole_frac_constraint(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    for j in m.fs.NG_props.component_list:
        if j not in m.fs.syn_props.component_list:
            m.fs.recycle_translator.outlet.mole_frac_comp[0, j].fix(0)

    # cathode translator
    @m.fs.cathode_translator.Constraint(m.fs.time)
    def temp_constraint(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.cathode_translator.Constraint(m.fs.time)
    def pressure_constraint(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.cathode_translator.Constraint(m.fs.time)
    def flow_constraint(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    m.fs.cathode_translator.outlet.mole_frac_comp.fix(0)
    m.fs.cathode_translator.outlet.mole_frac_comp[0, 'O2'].fix(1)

    # cathode exhaust translator
    @m.fs.cathode_exhaust_translator.Constraint(m.fs.time)
    def temp_constraint(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.cathode_exhaust_translator.Constraint(m.fs.time)
    def pressure_constraint(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.cathode_exhaust_translator.Constraint(m.fs.time)
    def flow_constraint(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    air_species = ["N2", "O2", "Ar", "H2O", "CO2"]

    @m.fs.cathode_exhaust_translator.Constraint(m.fs.time, air_species)
    def mole_frac_constraint(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    m.fs.cathode_exhaust_translator.outlet.mole_frac_comp.fix(0)
    m.fs.cathode_exhaust_translator.outlet.mole_frac_comp[0, 'H2O'].unfix()
    m.fs.cathode_exhaust_translator.outlet.mole_frac_comp[0, 'CO2'].unfix()
    m.fs.cathode_exhaust_translator.outlet.mole_frac_comp[0, 'O2'].unfix()
    m.fs.cathode_exhaust_translator.outlet.mole_frac_comp[0, 'N2'].unfix()
    m.fs.cathode_exhaust_translator.outlet.mole_frac_comp[0, 'Ar'].unfix()


def set_power_island_inputs(m):
    ###################
    # anode side inputs
    ###################

    # syngas feed conditions
    m.fs.anode_mix.feed.flow_mol.fix(3642)  # mol/s
    m.fs.anode_mix.feed.temperature.fix(626.5)  # K
    m.fs.anode_mix.feed.pressure.fix(137888)  # Pa, equal to 20 psia
    m.fs.anode_mix.feed.mole_frac_comp[0, 'CH4'].fix(0.1890)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'C2H6'].fix(0.0061)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'C3H8'].fix(0.0013)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'C4H10'].fix(0.0008)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'CO'].fix(0.0913)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'CO2'].fix(0.0437)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'H2'].fix(0.2761)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'H2O'].fix(0.1113)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'N2'].fix(0.2877)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'O2'].fix(0.0000)
    m.fs.anode_mix.feed.mole_frac_comp[0, 'Ar'].fix(0.0035)

    # anode heat exchanger
    m.fs.anode_hx.tube.deltaP.fix(-137.9)  # equal to -0.02 psi
    m.fs.anode_hx.shell.deltaP.fix(-137.9)  # equal to -0.02 psi
    m.fs.anode_hx.area.fix(12664)  # m2
    m.fs.anode_hx.overall_heat_transfer_coefficient.fix(80)  # W/m2K

    # prereformer and anode
    m.fs.prereformer.heat_duty.fix(0)
    m.fs.prereformer.deltaP.fix(-137.9)  # equal to -0.02 psi

    m.fs.anode.outlet.temperature.fix(979.44)  # K, does not stay fixed
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
    m.fs.cathode_hx.area.fix(17024)  # m2
    m.fs.cathode_hx.overall_heat_transfer_coefficient.fix(80)  # W/m2K

    # cathode
    m.fs.cathode.ion_outlet.flow_mol.fix(1670)
    m.fs.cathode.split_fraction[0, 'ion_outlet', 'H2O'].fix(0)
    m.fs.cathode.split_fraction[0, 'ion_outlet', 'CO2'].fix(0)
    m.fs.cathode.split_fraction[0, 'ion_outlet', 'N2'].fix(0)
    m.fs.cathode.split_fraction[0, 'ion_outlet', 'Ar'].fix(0)

    m.fs.cathode_heat.outlet.temperature.fix(1000)  # K
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


def initialize_power_island(m):
    # cathode side
    m.fs.air_blower.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.cathode_hx.tube_inlet, m.fs.air_blower.outlet)

    # fix cathode inlet to initial guess
    m.fs.cathode.inlet.flow_mol.fix(34161)
    m.fs.cathode.inlet.temperature.fix(891)
    m.fs.cathode.inlet.pressure.fix(105490)
    m.fs.cathode.inlet.mole_frac_comp[0, 'H2O'].fix(0.0109)
    m.fs.cathode.inlet.mole_frac_comp[0, 'CO2'].fix(0.0003)
    m.fs.cathode.inlet.mole_frac_comp[0, 'N2'].fix(0.8099)
    m.fs.cathode.inlet.mole_frac_comp[0, 'O2'].fix(0.1690)
    m.fs.cathode.inlet.mole_frac_comp[0, 'Ar'].fix(0.0099)

    m.fs.cathode.initialize(outlvl=logging.INFO)

    m.fs.cathode.inlet.unfix()

    # cathode translator block
    copy_port_values(
        m.fs.cathode_translator.inlet, m.fs.cathode.ion_outlet)

    calculate_variable_from_constraint(
        m.fs.cathode_translator.outlet.flow_mol[0],
        m.fs.cathode_translator.flow_constraint[0])

    calculate_variable_from_constraint(
        m.fs.cathode_translator.outlet.temperature[0],
        m.fs.cathode_translator.temp_constraint[0])

    calculate_variable_from_constraint(
        m.fs.cathode_translator.outlet.pressure[0],
        m.fs.cathode_translator.pressure_constraint[0])

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
    # fix anode inlet to initial guess
    m.fs.anode.inlet.flow_mol.fix(9801)
    m.fs.anode.inlet.temperature.fix(845)
    m.fs.anode.inlet.pressure.fix(137612)
    m.fs.anode.inlet.mole_frac_comp[0, 'CH4'].fix(0.057)
    m.fs.anode.inlet.mole_frac_comp[0, 'CO'].fix(0.043)
    m.fs.anode.inlet.mole_frac_comp[0, 'CO2'].fix(0.122)
    m.fs.anode.inlet.mole_frac_comp[0, 'H2'].fix(0.222)
    m.fs.anode.inlet.mole_frac_comp[0, 'H2O'].fix(0.187)
    m.fs.anode.inlet.mole_frac_comp[0, 'N2'].fix(0.195)
    m.fs.anode.inlet.mole_frac_comp[0, 'O2'].fix(0.172)
    m.fs.anode.inlet.mole_frac_comp[0, 'Ar'].fix(0.002)

    m.fs.anode.initialize(outlvl=logging.INFO)

    m.fs.anode.inlet.unfix()

    copy_port_values(
        m.fs.anode_recycle.inlet, m.fs.anode.outlet)

    m.fs.anode_recycle.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.anode_blower.inlet, m.fs.anode_recycle.recycle)

    m.fs.anode_blower.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.recycle_translator.inlet, m.fs.anode_blower.outlet)

    # anode recycle translator
    calculate_variable_from_constraint(
        m.fs.recycle_translator.outlet.flow_mol[0],
        m.fs.recycle_translator.flow_constraint[0])

    calculate_variable_from_constraint(
        m.fs.recycle_translator.outlet.temperature[0],
        m.fs.recycle_translator.temp_constraint[0])

    calculate_variable_from_constraint(
        m.fs.recycle_translator.outlet.pressure[0],
        m.fs.recycle_translator.pressure_constraint[0])

    for key in m.fs.recycle_translator.mole_frac_constraint.keys():
        calculate_variable_from_constraint(
            m.fs.recycle_translator.outlet.mole_frac_comp[key],
            m.fs.recycle_translator.mole_frac_constraint[key])

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

    m.fs.prereformer.gibbs_scaling = 1e-5

    m.fs.prereformer.initialize(outlvl=logging.INFO)

    copy_port_values(
        m.fs.anode_translator.inlet, m.fs.prereformer.outlet)

    # anode translator
    calculate_variable_from_constraint(
        m.fs.anode_translator.outlet.flow_mol[0],
        m.fs.anode_translator.flow_constraint[0])

    calculate_variable_from_constraint(
        m.fs.anode_translator.outlet.temperature[0],
        m.fs.anode_translator.temp_constraint[0])

    calculate_variable_from_constraint(
        m.fs.anode_translator.outlet.pressure[0],
        m.fs.anode_translator.pressure_constraint[0])

    for key in m.fs.anode_translator.mole_frac_constraint.keys():
        calculate_variable_from_constraint(
            m.fs.anode_translator.outlet.mole_frac_comp[key],
            m.fs.anode_translator.mole_frac_constraint[key])

    m.fs.anode_translator.initialize()

    copy_port_values(
        m.fs.fuel_cell_mix.fuel_inlet, m.fs.anode_translator.outlet)

    copy_port_values(
        m.fs.fuel_cell_mix.ion_inlet, m.fs.cathode_translator.outlet)

    m.fs.fuel_cell_mix.initialize(outlvl=logging.INFO)

    ############################
    # Combustor and HRSG section
    ############################

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

    # cathode exhaust translator
    calculate_variable_from_constraint(
        m.fs.cathode_exhaust_translator.outlet.flow_mol[0],
        m.fs.cathode_exhaust_translator.flow_constraint[0])

    calculate_variable_from_constraint(
        m.fs.cathode_exhaust_translator.outlet.temperature[0],
        m.fs.cathode_exhaust_translator.temp_constraint[0])

    calculate_variable_from_constraint(
        m.fs.cathode_exhaust_translator.outlet.pressure[0],
        m.fs.cathode_exhaust_translator.pressure_constraint[0])

    for key in m.fs.cathode_exhaust_translator.mole_frac_constraint.keys():
        calculate_variable_from_constraint(
            m.fs.cathode_exhaust_translator.outlet.mole_frac_comp[key],
            m.fs.cathode_exhaust_translator.mole_frac_constraint[key])

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

# %% Reformer Section


def build_reformer(m):
    # build property package
    natural_gas_config = get_SOFC_properties(
        components=['CH4', 'CO', 'CO2', 'H2', 'H2O', 'C2H6', 'C3H8', 'C4H10',
                    'N2', 'O2', 'Ar'])
    m.fs.natural_gas_props = GenericParameterBlock(default=natural_gas_config)

    # build unit models
    m.fs.reformer_recuperator = HeatExchanger(
        default={"delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "shell": {"property_package": m.fs.natural_gas_props},
                 "tube": {"property_package": m.fs.natural_gas_props}})

    m.fs.NG_expander = PressureChanger(
        default={'compressor': False,
                 'property_package': m.fs.natural_gas_props,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic})

    m.fs.reformer_bypass = Separator(
        default={"outlet_list": ["reformer_outlet", "bypass_outlet"],
                 "property_package": m.fs.natural_gas_props})

    m.fs.air_compressor_s1 = PressureChanger(
        default={"compressor": True,
                 "property_package": m.fs.natural_gas_props,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    m.fs.intercooler_s1 = Heater(
        default={"property_package": m.fs.natural_gas_props,
                 "has_pressure_change": True})

    m.fs.air_compressor_s2 = PressureChanger(
        default={"compressor": True,
                 "property_package": m.fs.natural_gas_props,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    m.fs.intercooler_s2 = Heater(
        default={"property_package": m.fs.natural_gas_props,
                 "has_pressure_change": True})

    m.fs.reformer_mix = Mixer(
        default={"inlet_list": ["gas_inlet", "oxygen_inlet", "steam_inlet"],
                 "property_package": m.fs.natural_gas_props})

    m.fs.reformer = GibbsReactor(
        default={"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "inert_species": ["N2", "Ar"],
                 "property_package": m.fs.natural_gas_props})

    m.fs.bypass_rejoin = Mixer(
        default={"inlet_list": ["syngas_inlet", "bypass_inlet"],
                 "property_package": m.fs.natural_gas_props})

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
    m.fs.air_compressor_s1.inlet.flow_mol.fix(1332.9)  # mol/s
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
    m.fs.reformer_mix.steam_inlet.flow_mol.fix(464.77)  # mol/s
    m.fs.reformer_mix.steam_inlet.temperature.fix(422)  # K
    m.fs.reformer_mix.steam_inlet.pressure.fix(206843)  # Pa, equal to 30 psia
    m.fs.reformer_mix.steam_inlet.mole_frac_comp.fix(0)
    m.fs.reformer_mix.steam_inlet.mole_frac_comp[0, 'H2O'].fix(1)

    # reformer outlet pressure
    m.fs.reformer.outlet.pressure.fix(137895)  # Pa, equal to 20 Psi
    m.fs.reformer.outlet.temperature.fix(1060.93)  # K, equal to 1450 F


def initialize_reformer(m):
    m.fs.reformer.inlet.flow_mol.fix(2262)  # mol/s
    m.fs.reformer.inlet.temperature.fix(470)  # K
    m.fs.reformer.inlet.pressure.fix(203395)  # Pa
    m.fs.reformer.inlet.mole_frac_comp[0, 'CH4'].fix(0.191)
    m.fs.reformer.inlet.mole_frac_comp[0, 'C2H6'].fix(0.006)
    m.fs.reformer.inlet.mole_frac_comp[0, 'C3H8'].fix(0.002)
    m.fs.reformer.inlet.mole_frac_comp[0, 'C4H10'].fix(0.001)
    m.fs.reformer.inlet.mole_frac_comp[0, 'H2'].fix(0)
    m.fs.reformer.inlet.mole_frac_comp[0, 'CO'].fix(0)
    m.fs.reformer.inlet.mole_frac_comp[0, 'CO2'].fix(0.002)
    m.fs.reformer.inlet.mole_frac_comp[0, 'H2O'].fix(0.212)
    m.fs.reformer.inlet.mole_frac_comp[0, 'N2'].fix(0.458)
    m.fs.reformer.inlet.mole_frac_comp[0, 'O2'].fix(0.122)
    m.fs.reformer.inlet.mole_frac_comp[0, 'Ar'].fix(0.006)

    m.fs.reformer.gibbs_scaling = 1e-5

    m.fs.reformer.lagrange_mult[0, 'Ar'] = 217511
    m.fs.reformer.lagrange_mult[0, 'C'] = 35595
    m.fs.reformer.lagrange_mult[0, 'H'] = 79455
    m.fs.reformer.lagrange_mult[0, 'N'] = 111756
    m.fs.reformer.lagrange_mult[0, 'O'] = 315023

    m.fs.reformer.heat_duty[0] = 16748188

    ref_solver = pyo.SolverFactory('ipopt')
    ref_solver.options = {'tol': 1e-6, 'bound_push': 1e-16}

    m.fs.reformer.control_volume.properties_in.initialize()
    m.fs.reformer.control_volume.properties_out.initialize()
    ref_solver.solve(m.fs.reformer.control_volume, tee=True)
    ref_solver.solve(m.fs.reformer, tee=True)

    m.fs.reformer.inlet.unfix()

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


# %% SOFC ROM


def build_SOFC_ROM(m):
    # create the ROM
    m.fs.SOFC = pyo.Block()
    make_SOFC_ROM(m.fs.SOFC)

    # flowsheet - ROM interconnection constraints
    def fuel_temp_match(fs):
        return(fs.SOFC.fuel_temperature ==
               fs.anode_mix.feed.temperature[0] - 273)
    m.fs.fuel_temp_constraint = pyo.Constraint(rule=fuel_temp_match)

    def air_temp_match(fs):
        return(fs.SOFC.air_temperature ==
               fs.cathode_mix.outlet.temperature[0] - 273)
    m.fs.air_temp_constraint = pyo.Constraint(rule=air_temp_match)

    def air_recirc_match(fs):
        return(fs.SOFC.air_recirculation ==
               fs.cathode_recycle.split_fraction[0, 'recycle'])
    m.fs.air_recirc_constraint = pyo.Constraint(rule=air_recirc_match)

    def OTC_rule(fs):
        O_frac = (1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO']
                  + 2 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO2']
                  + 1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'H2O'])

        C_frac = (1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO']
                  + 1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO2']
                  + 1 * fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CH4'])

        return fs.SOFC.OTC == O_frac/C_frac
    m.fs.OTC_constraint = pyo.Constraint(rule=OTC_rule)

    def fuel_utilization_rule(fs):
        full_O2_flow = (fs.anode_mix.feed.flow_mol[0] * (
            0.5 * fs.anode_mix.feed.mole_frac_comp[0, 'H2'] +
            0.5 * fs.anode_mix.feed.mole_frac_comp[0, 'CO'] +
            2.0 * fs.anode_mix.feed.mole_frac_comp[0, 'CH4'] +
            3.5 * fs.anode_mix.feed.mole_frac_comp[0, 'C2H6'] +
            5.0 * fs.anode_mix.feed.mole_frac_comp[0, 'C3H8'] +
            6.5 * fs.anode_mix.feed.mole_frac_comp[0, 'C4H10']))

        return (fs.SOFC.fuel_util*full_O2_flow ==
                fs.cathode.ion_outlet.flow_mol[0])

    m.fs.fuel_util_constraint = pyo.Constraint(rule=fuel_utilization_rule)

    def air_utilization_rule(fs):
        air_in = (fs.cathode_hx.tube_inlet.flow_mol[0] *
                  fs.cathode_hx.tube_inlet.mole_frac_comp[0, 'O2'])

        air_out = (fs.cathode_hx.shell_inlet.flow_mol[0] *
                   fs.cathode_hx.shell_inlet.mole_frac_comp[0, 'O2'])

        return fs.SOFC.air_util == 1 - air_out/air_in

    m.fs.air_util_constraint = pyo.Constraint(rule=air_utilization_rule)

    if hasattr(m.fs, "reformer_bypass"):
        def internal_ref_rule(fs):
            return (fs.SOFC.internal_reforming ==
                    fs.reformer_bypass.split_fraction[0, 'bypass_outlet'])
        m.fs.internal_ref_constraint = pyo.Constraint(rule=internal_ref_rule)

    else:
        m.fs.SOFC.internal_reforming.fix(0.6)

    m.fs.SOFC.current_density.fix(4000)
    m.fs.SOFC.pressure.fix(1)

    # constraints for reporting results
    m.fs.stack_power = pyo.Var(initialize=600e6)

    def stack_power_rule(fs):  # power in MW
        F = 96487  # C/mol e-
        current = m.fs.cathode.ion_outlet.flow_mol[0]*4*F
        return fs.stack_power == current*fs.SOFC.stack_voltage
    m.fs.stack_power_constraint = pyo.Constraint(rule=stack_power_rule)

    def cathode_heat_target(fs):
        return -1*fs.anode.heat_duty[0] - fs.stack_power
    m.fs.cathode_heat_eq = pyo.Expression(rule=cathode_heat_target)


def add_anode_temp_constraint(m):
    m.fs.anode.outlet.temperature.unfix()

    def outlet_fuel_temp_rule(fs):
        return (fs.anode.outlet.temperature[0] - 273 ==
                fs.SOFC.anode_outlet_temperature)
    m.fs.outet_fuel_temp_constraint = pyo.Constraint(
        rule=outlet_fuel_temp_rule)


def add_cathode_heat_constraint(m):
    m.fs.cathode_heat.outlet.temperature.unfix()

    def cathode_heat_rule(fs):
        return (fs.cathode_heat.heat_duty[0]*1e-6 ==
                (-1*fs.anode.heat_duty[0] - fs.stack_power)*1e-6)
    m.fs.cathode_heat_constraint = pyo.Constraint(rule=cathode_heat_rule)


# %% Results


def add_result_constraints(m, ref=True):
    # heat supplied to HRSG
    m.fs.HRSG_heat_duty = pyo.Var(initialize=50e6)

    def HRSG_heat_duty_rule(fs):
        return (-1*fs.HRSG_heat_duty ==
                fs.anode_HRSG.heat_duty[0] + fs.cathode_HRSG.heat_duty[0])
    m.fs.HRSG_heat_duty_constraint = pyo.Constraint(rule=HRSG_heat_duty_rule)

    # heat required for generating reformer steam
    m.fs.reformer_steam_heat = pyo.Var(initialize=1e6)

    if ref:
        # only valid when steam is heated to 422 K (149 C)
        def reformer_steam_heat_rule(fs):
            return (fs.reformer_steam_heat ==
                    50735 * fs.reformer_mix.steam_inlet.flow_mol[0])
        m.fs.reformer_steam_heat_eq = pyo.Constraint(
            rule=reformer_steam_heat_rule)

    else:
        m.fs.reformer_steam_heat.fix(0)

    # HRSG heat applied toward steam cycle
    m.fs.steam_cycle_heat = pyo.Var(initialize=50e6)

    def steam_cycle_heat_rule(fs):
        return (fs.steam_cycle_heat == fs.HRSG_heat_duty
                - fs.reformer_steam_heat - fs.reformer.heat_duty[0])
    m.fs.steam_cycle_heat_eq = pyo.Constraint(rule=steam_cycle_heat_rule)

    # power generated by steam cycle
    m.fs.steam_cycle_power = pyo.Var(initialize=10e6)
    m.fs.steam_cycle_efficiency = pyo.Param(initialize=0.38,
                                            mutable=True)

    def steam_cycle_power_rule(fs):
        return (fs.steam_cycle_power ==
                fs.steam_cycle_heat * fs.steam_cycle_efficiency)
    m.fs.steam_cycle_power_eq = pyo.Constraint(rule=steam_cycle_power_rule)

    # stack AC power
    m.fs.stack_power_AC = pyo.Var(initialize=600e6)
    m.fs.inverter_efficiency = pyo.Param(initialize=0.97,
                                         mutable=True)

    def stack_AC_power_rule(fs):
        return fs.stack_power_AC == fs.stack_power * fs.inverter_efficiency
    m.fs.stack_AC_power_constraint = pyo.Constraint(rule=stack_AC_power_rule)

    # gross plant power
    m.fs.gross_power = pyo.Var(initialize=700e6)

    def gross_power_rule(fs):
        return (fs.gross_power == (fs.stack_power_AC
                                   + fs.steam_cycle_power
                                   + -1*fs.NG_expander.work_mechanical[0]))
    m.fs.gross_power_constraint = pyo.Constraint(rule=gross_power_rule)

    # auxiliary load of the plant
    m.fs.auxiliary_load = pyo.Var(initialize=10e6)

    def auxiliary_load_rule(fs):
        return (fs.auxiliary_load == fs.air_blower.work_mechanical[0]
                + fs.cathode_blower.work_mechanical[0]
                + fs.anode_blower.work_mechanical[0]
                + fs.air_compressor_s1.work_mechanical[0]
                + fs.air_compressor_s2.work_mechanical[0])
    m.fs.auxiliary_load_constraint = pyo.Constraint(rule=auxiliary_load_rule)

    # net plant power
    m.fs.net_power = pyo.Var(initialize=600e6)

    def net_power_rule(fs):
        return fs.net_power == fs.gross_power - fs.auxiliary_load
    m.fs.net_power_constraint = pyo.Constraint(rule=net_power_rule)

    # CO2 emissions in g/kWh
    m.fs.CO2_emissions = pyo.Var(initialize=1000)

    def CO2_emission_rule(fs):
        mass_flow = (fs.anode_HRSG.outlet.flow_mol[0] *
                     fs.anode_HRSG.outlet.mole_frac_comp[0, 'CO2'] * 44.01)
        return fs.CO2_emissions == mass_flow/(fs.net_power/3.6e6)
    m.fs.CO2_eq = pyo.Constraint(rule=CO2_emission_rule)

    # HHV efficiency
    m.fs.HHV_efficiency = pyo.Var(initialize=1)

    def efficiency_rule(fs):
        NG_HHV = 908839.23  # J/mol
        return (fs.HHV_efficiency ==
                fs.net_power /
                (NG_HHV * fs.reformer_recuperator.tube_inlet.flow_mol[0]))
    m.fs.efficiency_eq = pyo.Constraint(rule=efficiency_rule)

    # scale result variables
    iscale.set_scaling_factor(m.fs.HRSG_heat_duty, 1e-6)
    iscale.set_scaling_factor(m.fs.reformer_steam_heat, 1e-6)
    iscale.set_scaling_factor(m.fs.steam_cycle_heat, 1e-6)
    iscale.set_scaling_factor(m.fs.steam_cycle_power, 1e-6)
    iscale.set_scaling_factor(m.fs.stack_power_AC, 1e-6)
    iscale.set_scaling_factor(m.fs.gross_power, 1e-6)
    iscale.set_scaling_factor(m.fs.auxiliary_load, 1e-6)
    iscale.set_scaling_factor(m.fs.net_power, 1e-6)

    # scale result constraints
    iscale.constraint_scaling_transform(m.fs.HRSG_heat_duty_constraint, 1e-6)
    iscale.constraint_scaling_transform(m.fs.reformer_steam_heat_eq, 1e-6)
    iscale.constraint_scaling_transform(m.fs.steam_cycle_heat_eq, 1e-6)
    iscale.constraint_scaling_transform(m.fs.steam_cycle_power_eq, 1e-6)
    iscale.constraint_scaling_transform(m.fs.stack_AC_power_constraint, 1e-6)
    iscale.constraint_scaling_transform(m.fs.gross_power_constraint, 1e-6)
    iscale.constraint_scaling_transform(m.fs.auxiliary_load_constraint, 1e-6)
    iscale.constraint_scaling_transform(m.fs.net_power_constraint, 1e-6)


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
    tags["power"] = fstr((pyo.value(m.fs.stack_power)/1e6), 1)
    tags["recuperator_duty"] = fstr(
        (pyo.value(m.fs.reformer_recuperator.heat_duty[0])/1e6), 1)
    tags["anode_hx_duty"] = fstr(
        (pyo.value(m.fs.anode_hx.heat_duty[0])/1e6), 2)
    tags["cathode_hx_duty"] = fstr(
        (pyo.value(m.fs.cathode_hx.heat_duty[0])/1e6), 2)
    tags["anode_duty"] = fstr((pyo.value(m.fs.anode.heat_duty[0])/1e6), 2)
    tags["cathode_duty"] = fstr(
        (pyo.value(m.fs.cathode_heat.heat_duty[0])/1e6), 2)
    tags["stack_power_AC"] = fstr(pyo.value(m.fs.stack_power_AC)/1e6, 1)
    tags["HRSG_power"] = fstr(pyo.value(m.fs.steam_cycle_power)/1e6, 1)
    tags["expander_power"] = fstr(-1*pyo.value(
        m.fs.NG_expander.work_mechanical[0])/1e6, 1)
    tags["aux_load"] = fstr(pyo.value(m.fs.auxiliary_load)/1e6, 1)
    tags["net_power"] = fstr(pyo.value(m.fs.net_power)/1e6, 1)
    tags["thermal_input"] = fstr(pyo.value(
        m.fs.net_power/m.fs.HHV_efficiency)/1e6, 1)
    tags["efficiency"] = fstr(pyo.value(m.fs.HHV_efficiency)*100, 2)
    tags["emissions"] = fstr(pyo.value(m.fs.CO2_emissions), 1)

    original_svg_file = os.path.join(this_file_dir(),
                                     "NGFC_results_template.svg")
    with open(original_svg_file, "r") as f:
        svg_tag(tags, f, outfile=outfile)
