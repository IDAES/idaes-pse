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

import os
import sys
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
from idaes.core.util.tags import svg_tag
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
    Flash,
    Separator,
    Translator)
from idaes.generic_models.unit_models.heat_exchanger import \
    delta_temperature_underwood_callback
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
from idaes.generic_models.unit_models.separator import SplittingType
from idaes.generic_models.unit_models.mixer import MomentumMixingType

from idaes.power_generation.flowsheets.sofc.properties.CO2_H2O_Ideal_VLE_scaled import \
    configuration as CO2_H2O_VLE_config
from idaes.power_generation.flowsheets.sofc.properties.natural_gas_PR_scaled_units import \
    get_NG_properties, rxn_configuration

from idaes.power_generation.flowsheets.sofc.surrogates.cpu import CPU
from idaes.power_generation.flowsheets.sofc.surrogates.sofc_rom_builder \
    import build_SOFC_ROM, initialize_SOFC_ROM
from idaes.power_generation.flowsheets.sofc.surrogates.sofc_surrogate \
    import build_SOFC_SM, initialize_SOFC_SM

from idaes.power_generation.costing.power_plant_costing import \
    get_fixed_OM_costs, get_variable_OM_costs, initialize_fixed_OM_costs, \
    initialize_variable_OM_costs


def build_NGFC(m):
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # create property packages - 4 property packages and 1 reaction
    NG_config = get_NG_properties(
        components=['H2', 'CO', "H2O", 'CO2', 'CH4', "C2H6", "C3H8", "C4H10",
                    'N2', 'O2', 'Ar'])
    m.fs.NG_props = GenericParameterBlock(default=NG_config)

    syn_config = get_NG_properties(
        components=["H2", "CO", "H2O", "CO2", "CH4", "N2", "O2", "Ar"])
    m.fs.syn_props = GenericParameterBlock(default=syn_config)

    air_config = get_NG_properties(
        components=['H2O', 'CO2', 'N2', 'O2', 'Ar'])

    m.fs.air_props = GenericParameterBlock(default=air_config)

    m.fs.CO2_H2O_VLE = GenericParameterBlock(default=CO2_H2O_VLE_config)

    m.fs.rxn_props = GenericReactionParameterBlock(
        default={"property_package": m.fs.syn_props, **rxn_configuration})

    # build anode side units
    m.fs.anode_mix = Mixer(
        default={"inlet_list": ["feed_inlet", "recycle_inlet"],
                 "property_package": m.fs.NG_props})

    m.fs.anode_hx = HeatExchanger(
        default={"delta_temperature_callback":
                 delta_temperature_underwood_callback,
                 "tube": {"property_package": m.fs.NG_props,
                          "has_pressure_change": True},
                 "shell": {"property_package": m.fs.syn_props,
                           "has_pressure_change": True}})

    m.fs.prereformer = GibbsReactor(
        default={"has_heat_transfer": False,
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
        default={"outlet_list": ["exhaust_outlet", "recycle_outlet"],
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
        default={"inlet_list": ["air_inlet", "recycle_inlet"],
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
        default={"outlet_list": ["exhaust_outlet", "recycle_outlet"],
                 "property_package": m.fs.air_props})

    m.fs.cathode_blower = PressureChanger(
        default={'compressor': True,
                 'property_package': m.fs.air_props,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic})

    m.fs.cathode_HRSG = Heater(
        default={"has_pressure_change": True,
                 "property_package": m.fs.air_props})

    # build ASU, oxycombustor and CPU
    m.fs.air_compressor_s1 = PressureChanger(
        default={"compressor": True,
                 "property_package": m.fs.air_props,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    m.fs.intercooler_s1 = Heater(
        default={"property_package": m.fs.air_props,
                 "has_pressure_change": True})

    m.fs.air_compressor_s2 = PressureChanger(
        default={"compressor": True,
                 "property_package": m.fs.air_props,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    m.fs.intercooler_s2 = Heater(
        default={"property_package": m.fs.air_props,
                 "has_pressure_change": True})

    m.fs.ASU = Separator(
        default={"outlet_list": ["N2_outlet", "O2_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.air_props})

    m.fs.ASU_O2_outlet = Heater(
        default={"has_pressure_change": True,
                 "property_package": m.fs.air_props})

    m.fs.oxycombustor_translator = Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": m.fs.air_props,
                 "outlet_property_package": m.fs.syn_props})

    @m.fs.oxycombustor_translator.Constraint(m.fs.time)
    def oxycombustor_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.oxycombustor_translator.Constraint(m.fs.time)
    def oxycombustor_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.oxycombustor_translator.Constraint(m.fs.time)
    def oxycombustor_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.oxycombustor_translator.Constraint(m.fs.time,
                                             m.fs.air_props.component_list)
    def oxycombustor_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    for j in m.fs.syn_props.component_list:
        if j not in m.fs.air_props.component_list:
            m.fs.oxycombustor_translator.outlet.mole_frac_comp[0, j].fix(0)

    m.fs.combustor_mix = Mixer(
        default={"inlet_list": ["anode_inlet", "oxygen_inlet"],
                 "property_package": m.fs.syn_props})

    m.fs.oxycombustor = StoichiometricReactor(
        default={"has_heat_of_reaction": False,
                 "has_heat_transfer": False,
                 "has_pressure_change": True,
                 "property_package": m.fs.syn_props,
                 "reaction_package": m.fs.rxn_props})

    @m.fs.Constraint()
    def air_to_combustor(fs):
        XS = 1.09  # excess oxygen

        F_CH4 = (fs.combustor_mix.anode_inlet.flow_mol[0] *
                 fs.combustor_mix.anode_inlet.mole_frac_comp[(0, "CH4")])
        F_CO = (fs.combustor_mix.anode_inlet.flow_mol[0] *
                fs.combustor_mix.anode_inlet.mole_frac_comp[(0, "CO")])
        F_H2 = (fs.combustor_mix.anode_inlet.flow_mol[0] *
                fs.combustor_mix.anode_inlet.mole_frac_comp[(0, "H2")])
        O2_required = 2*F_CH4 + 0.5*F_CO + 0.5*F_H2

        O2_fed = (fs.combustor_mix.oxygen_inlet.flow_mol[0] *
                  fs.combustor_mix.oxygen_inlet.mole_frac_comp[(0, "O2")])
        return O2_fed == XS*O2_required

    m.fs.anode_HRSG = Heater(
        default={"has_pressure_change": True,
                 "property_package": m.fs.syn_props})

    m.fs.CPU_translator = Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": m.fs.syn_props,
                 "outlet_property_package": m.fs.CO2_H2O_VLE})

    @m.fs.CPU_translator.Constraint(m.fs.time)
    def CPU_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.CPU_translator.Constraint(m.fs.time)
    def CPU_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.CPU_translator.Constraint(m.fs.time)
    def CPU_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.CPU_translator.Constraint(m.fs.time, m.fs.CO2_H2O_VLE.component_list)
    def CPU_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    m.fs.condenser = Heater(
        default={"has_pressure_change": True,
                 "property_package": m.fs.CO2_H2O_VLE})

    m.fs.flash = Flash(
        default={"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "property_package": m.fs.CO2_H2O_VLE})

    m.fs.CPU = CPU()

    # build arcs to connect unit operations
    # arcs for anode side
    m.fs.ANODE_MIXER = Arc(
        source=m.fs.anode_mix.outlet,
        destination=m.fs.anode_hx.tube_inlet)

    m.fs.PREREF_IN = Arc(
        source=m.fs.anode_hx.tube_outlet,
        destination=m.fs.prereformer.inlet)

    m.fs.PREREF_OUT = Arc(
        source=m.fs.prereformer.outlet,
        destination=m.fs.anode_translator.inlet)

    m.fs.ANODE_TRANS_OUT = Arc(
        source=m.fs.anode_translator.outlet,
        destination=m.fs.fuel_cell_mix.fuel_inlet)

    m.fs.ANODE_IN = Arc(
        source=m.fs.fuel_cell_mix.outlet,
        destination=m.fs.anode.inlet)

    m.fs.ANODE = Arc(
        source=m.fs.anode.outlet,
        destination=m.fs.anode_recycle.inlet)

    m.fs.ANODE_RECYCLE_HX = Arc(
        source=m.fs.anode_recycle.exhaust_outlet,
        destination=m.fs.anode_hx.shell_inlet)

    m.fs.ANODE_RECYCLE = Arc(
        source=m.fs.anode_recycle.recycle_outlet,
        destination=m.fs.anode_blower.inlet)

    m.fs.ANODE_BLOWER = Arc(
        source=m.fs.anode_blower.outlet,
        destination=m.fs.recycle_translator.inlet)

    m.fs.REC_TRANS_OUT = Arc(
        source=m.fs.recycle_translator.outlet,
        destination=m.fs.anode_mix.recycle_inlet)

    # arcs for cathode side
    m.fs.AIR_BLOWER = Arc(
        source=m.fs.air_blower.outlet,
        destination=m.fs.cathode_hx.tube_inlet)

    m.fs.CATHODE_HX_COLD = Arc(
        source=m.fs.cathode_hx.tube_outlet,
        destination=m.fs.cathode_mix.air_inlet)

    m.fs.CATHODE_MIXER = Arc(
        source=m.fs.cathode_mix.outlet,
        destination=m.fs.cathode.inlet)

    m.fs.O2_IONS = Arc(
        source=m.fs.cathode.ion_outlet,
        destination=m.fs.cathode_translator.inlet)

    m.fs.O2_IONS_TRANSLATE = Arc(
        source=m.fs.cathode_translator.outlet,
        destination=m.fs.fuel_cell_mix.ion_inlet)

    m.fs.CATHODE = Arc(
        source=m.fs.cathode.air_outlet,
        destination=m.fs.cathode_heat.inlet)

    m.fs.CATHODE_HEAT = Arc(
        source=m.fs.cathode_heat.outlet,
        destination=m.fs.cathode_recycle.inlet)

    m.fs.CATHODE_RECYCLE_HX = Arc(
        source=m.fs.cathode_recycle.exhaust_outlet,
        destination=m.fs.cathode_hx.shell_inlet)

    m.fs.CATHODE_RECYCLE = Arc(
        source=m.fs.cathode_recycle.recycle_outlet,
        destination=m.fs.cathode_blower.inlet)

    m.fs.CATHODE_BLOWER = Arc(
        source=m.fs.cathode_blower.outlet,
        destination=m.fs.cathode_mix.recycle_inlet)

    m.fs.CATHODE_HRSG_IN = Arc(
        source=m.fs.cathode_hx.shell_outlet,
        destination=m.fs.cathode_HRSG.inlet)

    # arcs for ASU, oxycombustor and CPU
    m.fs.STAGE_1_OUT = Arc(
        source=m.fs.air_compressor_s1.outlet,
        destination=m.fs.intercooler_s1.inlet)

    m.fs.IC_1_OUT = Arc(
        source=m.fs.intercooler_s1.outlet,
        destination=m.fs.air_compressor_s2.inlet)

    m.fs.STAGE_2_OUT = Arc(
        source=m.fs.air_compressor_s2.outlet,
        destination=m.fs.intercooler_s2.inlet)

    m.fs.TO_ASU = Arc(
        source=m.fs.intercooler_s2.outlet,
        destination=m.fs.ASU.inlet)

    m.fs.ASU_OUT = Arc(
        source=m.fs.ASU.O2_outlet,
        destination=m.fs.ASU_O2_outlet.inlet)

    m.fs.TO_O2_TRANSLATOR = Arc(
        source=m.fs.ASU_O2_outlet.outlet,
        destination=m.fs.oxycombustor_translator.inlet)

    m.fs.ANODE_HX_HOT_OUT = Arc(
        source=m.fs.anode_hx.shell_outlet,
        destination=m.fs.combustor_mix.anode_inlet)

    m.fs.COMBUST_O2_IN = Arc(
        source=m.fs.oxycombustor_translator.outlet,
        destination=m.fs.combustor_mix.oxygen_inlet)

    m.fs.COMBUST_IN = Arc(
        source=m.fs.combustor_mix.outlet,
        destination=m.fs.oxycombustor.inlet)

    m.fs.COMBUST_OUT = Arc(
        source=m.fs.oxycombustor.outlet,
        destination=m.fs.anode_HRSG.inlet)

    m.fs.ANODE_HRSG_OUT = Arc(
        source=m.fs.anode_HRSG.outlet,
        destination=m.fs.CPU_translator.inlet)

    m.fs.CONDENSE_IN = Arc(
        source=m.fs.CPU_translator.outlet,
        destination=m.fs.condenser.inlet)

    m.fs.FLASH_IN = Arc(
        source=m.fs.condenser.outlet,
        destination=m.fs.flash.inlet)

    # converting from kmol/s to mol/s
    @m.fs.Constraint(m.fs.time)
    def CPU_inlet_F(fs, t):
        return fs.CPU.inlet.flow_mol[t] == fs.flash.vap_outlet.flow_mol[t]*1000

    @m.fs.Constraint(m.fs.time)
    def CPU_inlet_T(fs, t):
        return (fs.CPU.inlet.temperature[t] ==
                fs.flash.vap_outlet.temperature[t])

    @m.fs.Constraint(m.fs.time)
    def CPU_inlet_P(fs, t):
        return fs.CPU.inlet.pressure[t] == fs.flash.vap_outlet.pressure[t]*1000

    @m.fs.Constraint(m.fs.time, m.fs.CO2_H2O_VLE.component_list)
    def CPU_inlet_x(fs, t, j):
        return (fs.CPU.inlet.mole_frac_comp[t, j] ==
                fs.flash.vap_outlet.mole_frac_comp[t, j])

    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)


def set_NGFC_inputs(m):
    ###################
    # anode side inputs
    ###################
    # natural gas feed conditions
    m.fs.anode_mix.feed_inlet.flow_mol.fix(1.0946)  # kmol/s
    m.fs.anode_mix.feed_inlet.temperature.fix(288.15)  # K
    m.fs.anode_mix.feed_inlet.pressure.fix(137.895)  # kPa (20 psia)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'CH4'].fix(0.931)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'C2H6'].fix(0.032)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'C3H8'].fix(0.007)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'C4H10'].fix(0.004)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'CO'].fix(0)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'CO2'].fix(0.01)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'H2'].fix(0)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'H2O'].fix(0)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'N2'].fix(0.016)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'O2'].fix(0)
    m.fs.anode_mix.feed_inlet.mole_frac_comp[0, 'Ar'].fix(0)

    # anode heat exchanger
    m.fs.anode_hx.tube.deltaP.fix(-.1379)  # kPa (-0.02 psi)
    m.fs.anode_hx.shell.deltaP.fix(-.1379)  # kPa (-0.02 psi)
    m.fs.anode_hx.area.fix(12664)  # m2
    m.fs.anode_hx.overall_heat_transfer_coefficient.fix(80e-3)  # mW/m^2K

    # prereformer and anode
    m.fs.prereformer.deltaP.fix(-.1379)  # kPa (-0.02 psi)

    m.fs.anode.outlet.temperature.fix(1001.5)  # K, unfixed by ROM
    m.fs.anode.outlet.pressure.fix(137.137)  # kPa (19.89 psia)

    # anode recycle and blower
    m.fs.anode_recycle.split_fraction[0, 'recycle_outlet'].fix(0.627)  # unfixed by ROM

    m.fs.anode_blower.outlet.pressure.fix(137.888)  # kPa (20 psia)
    m.fs.anode_blower.efficiency_isentropic.fix(0.8)

    #####################
    # cathode side inputs
    #####################

    # air feed conditions
    m.fs.air_blower.inlet.flow_mol.fix(11.444)  # kmol/s, unfixed by ROM
    m.fs.air_blower.inlet.temperature.fix(288.15)  # K (59 F)
    m.fs.air_blower.inlet.pressure.fix(101.353)  # kPa (14.7 psia)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'H2O'].fix(0.0104)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'CO2'].fix(0.0003)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'N2'].fix(0.7722)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'O2'].fix(0.2077)
    m.fs.air_blower.inlet.mole_frac_comp[0, 'Ar'].fix(0.0094)

    # air blower
    m.fs.air_blower.outlet.pressure.fix(111.006)  # kPa (16.1 psia)
    m.fs.air_blower.efficiency_isentropic.fix(0.82)

    # cathode heat exchanger
    m.fs.cathode_hx.tube_outlet.pressure.fix(105.490)  # kPa (15.3 psia)
    m.fs.cathode_hx.shell.deltaP.fix(-1.379)  # kPa (-0.2 psi)
    m.fs.cathode_hx.area.fix(14098)  # m2, unfixed by ROM
    m.fs.cathode_hx.overall_heat_transfer_coefficient.fix(80e-3)  # mW/m2K

    # cathode
    m.fs.cathode.ion_outlet.flow_mol.fix(1.893)  # kmol/s, unfixed by ROM

    m.fs.cathode_heat.outlet.temperature.fix(972.5)  # K, unfixed by ROM
    m.fs.cathode_heat.outlet.pressure.fix(104.111)  # kPa (15.1 psia)

    # cathode recycle and blower
    m.fs.cathode_recycle.split_fraction[0, 'recycle_outlet'].fix(0.5)

    m.fs.cathode_blower.outlet.pressure.fix(105.490)  # kPa (15.3 psi)
    m.fs.cathode_blower.efficiency_isentropic.fix(0.8)

    m.fs.cathode_HRSG.outlet.temperature.fix(405.7)  # K (270 F)
    m.fs.cathode_HRSG.deltaP.fix(-1.379)  # kPa (-0.2 psi)

    #####################################
    # ASU, oxycombustor, and CPU inputs #
    #####################################
    # air to ASU
    m.fs.air_compressor_s1.inlet.flow_mol[0] = 1.809  # kmol/s
    m.fs.air_compressor_s1.inlet.temperature.fix(288.15)  # K
    m.fs.air_compressor_s1.inlet.pressure.fix(101.353)  # kPa (14.7 psia)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'CO2'].fix(0.0003)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'H2O'].fix(0.0104)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'N2'].fix(0.7722)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'O2'].fix(0.2077)
    m.fs.air_compressor_s1.inlet.mole_frac_comp[0, 'Ar'].fix(0.0094)

    # air compressors and intercoolers
    m.fs.air_compressor_s1.outlet.pressure.fix(234.422)  # kPa (34 psia)
    m.fs.air_compressor_s1.efficiency_isentropic.fix(0.84)

    m.fs.intercooler_s1.outlet.temperature.fix(310.93)  # K (100 F)
    m.fs.intercooler_s1.deltaP.fix(-3.447)  # kPa (-0.5 psi)

    m.fs.air_compressor_s2.outlet.pressure.fix(544.686)  # kPa (79 psia)
    m.fs.air_compressor_s2.efficiency_isentropic.fix(0.84)

    m.fs.intercooler_s2.outlet.temperature.fix(310.93)  # K (100 F)
    m.fs.intercooler_s2.deltaP.fix(-3.447)  # kPa (-0.5 psi)

    # air seperation unit
    m.fs.ASU.split_fraction[0, "O2_outlet", "CO2"].fix(0)
    m.fs.ASU.split_fraction[0, "O2_outlet", "H2O"].fix(0)
    m.fs.ASU.split_fraction[0, "O2_outlet", "N2"].fix(0.0005)
    m.fs.ASU.split_fraction[0, "O2_outlet", "O2"].fix(0.9691)
    m.fs.ASU.split_fraction[0, "O2_outlet", "Ar"].fix(0.0673)

    m.fs.ASU_O2_outlet.outlet.temperature.fix(299.82)  # K (80 F)
    m.fs.ASU_O2_outlet.outlet.pressure.fix(158.579)  # kPa (23 psia)

    m.fs.oxycombustor.outlet.temperature.setub(2000)
    m.fs.oxycombustor.deltaP.fix(-6.895)  # kPa (-1 psi)

    m.fs.oxycombustor.outlet.mole_frac_comp[0, "H2"].fix(0)
    m.fs.oxycombustor.outlet.mole_frac_comp[0, "CO"].fix(0)
    m.fs.oxycombustor.outlet.mole_frac_comp[0, "CH4"].fix(0)

    m.fs.anode_HRSG.inlet.temperature.setub(2000)
    m.fs.anode_HRSG.outlet.temperature.fix(405)  # K
    m.fs.anode_HRSG.deltaP.fix(-6.895)  # kPa (-1 psi)

    m.fs.condenser.outlet.temperature.fix(310.9)  # K (100 F)
    m.fs.condenser.deltaP.fix(-6.895)  # kPa (-1 psi)

    m.fs.flash.control_volume.properties_out[0].temperature.fix(310.9)  # K
    m.fs.flash.control_volume.properties_out[0].pressure.fix(101.3529)  # kPa (14.7 psia)


def scale_NGFC(m):
    # NG_props
    m.fs.NG_props.set_default_scaling("flow_mol", 1)
    m.fs.NG_props.set_default_scaling("flow_mol_phase", 1)
    m.fs.NG_props.set_default_scaling("temperature", 1e-2)
    m.fs.NG_props.set_default_scaling("pressure", 1e-2)
    m.fs.NG_props.set_default_scaling("mole_frac_comp", 1e2)
    m.fs.NG_props.set_default_scaling("mole_frac_phase_comp", 1e2)

    m.fs.NG_props.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.NG_props.set_default_scaling("entr_mol_phase", 1e-1)

    # syn_props
    m.fs.syn_props.set_default_scaling("flow_mol", 1)
    m.fs.syn_props.set_default_scaling("flow_mol_phase", 1)
    m.fs.syn_props.set_default_scaling("temperature", 1e-2)
    m.fs.syn_props.set_default_scaling("pressure", 1e-2)
    m.fs.syn_props.set_default_scaling("mole_frac_comp", 1e2)
    m.fs.syn_props.set_default_scaling("mole_frac_phase_comp", 1e2)

    m.fs.syn_props.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.syn_props.set_default_scaling("entr_mol_phase", 1e-1)

    # air_props
    m.fs.air_props.set_default_scaling("flow_mol", 1)
    m.fs.air_props.set_default_scaling("flow_mol_phase", 1)
    m.fs.air_props.set_default_scaling("temperature", 1e-2)
    m.fs.air_props.set_default_scaling("pressure", 1e-2)
    m.fs.air_props.set_default_scaling("mole_frac_comp", 1e2)
    m.fs.air_props.set_default_scaling("mole_frac_phase_comp", 1e2)

    m.fs.air_props.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.air_props.set_default_scaling("entr_mol_phase", 1e-1)

    # CO2_H2O_VLE
    m.fs.CO2_H2O_VLE.set_default_scaling("flow_mol", 1)
    m.fs.CO2_H2O_VLE.set_default_scaling("flow_mol_phase", 1)
    m.fs.CO2_H2O_VLE.set_default_scaling("temperature", 1e-2)
    m.fs.CO2_H2O_VLE.set_default_scaling("pressure", 1e-2)
    m.fs.CO2_H2O_VLE.set_default_scaling("mole_frac_comp", 1e2)
    m.fs.CO2_H2O_VLE.set_default_scaling("mole_frac_phase_comp", 1e2)

    m.fs.CO2_H2O_VLE.set_default_scaling("enth_mol_phase", 1e-3)
    m.fs.CO2_H2O_VLE.set_default_scaling("entr_mol_phase", 1e-1)

    # heat exchanger areas and overall heat transfer coefficiencts
    iscale.set_scaling_factor(m.fs.anode_hx.area, 1e-4)
    iscale.set_scaling_factor(m.fs.anode_hx.overall_heat_transfer_coefficient, 1e3)
    iscale.set_scaling_factor(m.fs.cathode_hx.area, 1e-4)
    iscale.set_scaling_factor(m.fs.cathode_hx.overall_heat_transfer_coefficient, 1e3)

    # control volume heats
    iscale.set_scaling_factor(m.fs.anode_hx.tube.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.anode_hx.shell.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.anode.control_volume.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.cathode_hx.tube.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.cathode_hx.shell.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.cathode_heat.control_volume.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.cathode_HRSG.control_volume.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.intercooler_s1.control_volume.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.intercooler_s2.control_volume.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.ASU_O2_outlet.control_volume.heat, 1e-1)
    iscale.set_scaling_factor(m.fs.anode_HRSG.control_volume.heat, 1e-5)
    iscale.set_scaling_factor(m.fs.condenser.control_volume.heat, 1e-4)
    iscale.set_scaling_factor(m.fs.flash.control_volume.heat, 1e-2)

    # work
    iscale.set_scaling_factor(m.fs.anode_blower.control_volume.work, 1e-2)
    iscale.set_scaling_factor(m.fs.air_blower.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.cathode_blower.control_volume.work, 1e-2)
    iscale.set_scaling_factor(m.fs.air_compressor_s1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.air_compressor_s2.control_volume.work, 1e-3)

    # reaction extents
    iscale.set_scaling_factor(m.fs.oxycombustor.control_volume.rate_reaction_extent[0, "R1"], 1e2)
    iscale.set_scaling_factor(m.fs.oxycombustor.control_volume.rate_reaction_extent[0, "R2"], 1e2)
    iscale.set_scaling_factor(m.fs.oxycombustor.control_volume.rate_reaction_extent[0, "R3"], 1e5)

    iscale.calculate_scaling_factors(m)


def initialize_NGFC(m):
    # cathode side
    m.fs.air_blower.initialize()

    copy_port_values(source=m.fs.air_blower.outlet,
                     destination=m.fs.cathode_hx.tube_inlet)

    # fix cathode inlet to initial guess
    m.fs.cathode.inlet.flow_mol[0] = 20.995  # kmol/s
    m.fs.cathode.inlet.temperature[0] = 869.55  # K
    m.fs.cathode.inlet.pressure[0] = 105.4895  # kPa
    m.fs.cathode.inlet.mole_frac_comp[0, 'H2O'] = 0.0113
    m.fs.cathode.inlet.mole_frac_comp[0, 'CO2'] = 0.0003
    m.fs.cathode.inlet.mole_frac_comp[0, 'N2'] = 0.8418
    m.fs.cathode.inlet.mole_frac_comp[0, 'O2'] = 0.1362
    m.fs.cathode.inlet.mole_frac_comp[0, 'Ar'] = 0.0102

    m.fs.cathode.initialize()

    # cathode translator block
    copy_port_values(source=m.fs.cathode.ion_outlet,
                     destination=m.fs.cathode_translator.inlet)

    m.fs.cathode_translator.initialize()

    # rest of cathode side
    copy_port_values(source=m.fs.cathode.air_outlet,
                     destination=m.fs.cathode_heat.inlet)

    m.fs.cathode_heat.initialize()

    copy_port_values(source=m.fs.cathode_heat.outlet,
                     destination=m.fs.cathode_recycle.inlet)

    m.fs.cathode_recycle.initialize()

    copy_port_values(source=m.fs.cathode_recycle.exhaust_outlet,
                     destination=m.fs.cathode_hx.shell_inlet)

    m.fs.cathode_hx.initialize()

    copy_port_values(source=m.fs.cathode_recycle.recycle_outlet,
                     destination=m.fs.cathode_blower.inlet)

    m.fs.cathode_blower.initialize()

    copy_port_values(source=m.fs.cathode_blower.outlet,
                     destination=m.fs.cathode_mix.recycle_inlet)

    copy_port_values(source=m.fs.cathode_hx.tube_outlet,
                     destination=m.fs.cathode_mix.air_inlet)

    m.fs.cathode_mix.initialize()

    copy_port_values(source=m.fs.cathode_hx.shell_outlet,
                     destination=m.fs.cathode_HRSG.inlet)

    m.fs.cathode_HRSG.initialize()

    # anode side
    # prereformer outlet tear stream
    m.fs.fuel_cell_mix.fuel_inlet.flow_mol[0] = 7.207  # kmol/s
    m.fs.fuel_cell_mix.fuel_inlet.temperature[0] = 795.6  # K
    m.fs.fuel_cell_mix.fuel_inlet.pressure[0] = 137.6  # kPa
    m.fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CH4'] = 0.1234
    m.fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO'] = 0.0425
    m.fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'CO2'] = 0.2581
    m.fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'H2'] = 0.2379
    m.fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'H2O'] = 0.3316
    m.fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'N2'] = 0.0065
    m.fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'O2'] = 0
    m.fs.fuel_cell_mix.fuel_inlet.mole_frac_comp[0, 'Ar'] = 0

    copy_port_values(source=m.fs.cathode_translator.outlet,
                     destination=m.fs.fuel_cell_mix.ion_inlet)

    m.fs.fuel_cell_mix.initialize()

    copy_port_values(source=m.fs.fuel_cell_mix.outlet,
                     destination=m.fs.anode.inlet)

    m.fs.anode.lagrange_mult[0, "C"] = 48722
    m.fs.anode.lagrange_mult[0, "H"] = 77156
    m.fs.anode.lagrange_mult[0, "O"] = 291729
    m.fs.anode.outlet.mole_frac_comp[0, "O2"] = 0
    m.fs.anode.gibbs_scaling = 1e-4

    m.fs.anode.initialize()

    copy_port_values(source=m.fs.anode.outlet,
                     destination=m.fs.anode_recycle.inlet)

    m.fs.anode_recycle.initialize()

    copy_port_values(source=m.fs.anode_recycle.recycle_outlet,
                     destination=m.fs.anode_blower.inlet)

    m.fs.anode_blower.initialize()

    copy_port_values(source=m.fs.anode_blower.outlet,
                     destination=m.fs.recycle_translator.inlet)

    m.fs.recycle_translator.initialize()

    copy_port_values(source=m.fs.recycle_translator.outlet,
                     destination=m.fs.anode_mix.recycle_inlet)

    m.fs.anode_mix.initialize()

    copy_port_values(source=m.fs.anode_mix.outlet,
                     destination=m.fs.anode_hx.tube_inlet)

    copy_port_values(source=m.fs.anode_recycle.exhaust_outlet,
                     destination=m.fs.anode_hx.shell_inlet)

    m.fs.anode_hx.initialize()

    copy_port_values(source=m.fs.anode_hx.tube_outlet,
                     destination=m.fs.prereformer.inlet)

    copy_port_values(source=m.fs.prereformer.inlet,
                     destination=m.fs.prereformer.outlet)

    m.fs.prereformer.lagrange_mult[0, "C"] = 9910
    m.fs.prereformer.lagrange_mult[0, "H"] = 62556
    m.fs.prereformer.lagrange_mult[0, "O"] = 293466
    m.fs.prereformer.outlet.mole_frac_comp[0, "O2"] = 0
    m.fs.prereformer.outlet.mole_frac_comp[0, "Ar"] = 0
    m.fs.prereformer.outlet.mole_frac_comp[0, "C2H6"] = 0
    m.fs.prereformer.outlet.mole_frac_comp[0, "C3H8"] = 0
    m.fs.prereformer.outlet.mole_frac_comp[0, "C4H10"] = 0
    m.fs.prereformer.gibbs_scaling = 1e-4

    m.fs.prereformer.initialize()

    copy_port_values(source=m.fs.prereformer.outlet,
                     destination=m.fs.anode_translator.inlet)

    m.fs.anode_translator.initialize()

    #############################
    # ASU, oxycombustor and CPU #
    #############################
    # air compressor train
    m.fs.air_compressor_s1.initialize()

    copy_port_values(source=m.fs.air_compressor_s1.outlet,
                     destination=m.fs.intercooler_s1.inlet)

    m.fs.intercooler_s1.initialize()

    copy_port_values(source=m.fs.intercooler_s1.outlet,
                     destination=m.fs.air_compressor_s2.inlet)

    m.fs.air_compressor_s2.initialize()

    copy_port_values(source=m.fs.air_compressor_s2.outlet,
                     destination=m.fs.intercooler_s2.inlet)

    m.fs.intercooler_s2.initialize()

    copy_port_values(source=m.fs.intercooler_s2.outlet,
                     destination=m.fs.ASU.inlet)

    m.fs.ASU.initialize()

    copy_port_values(source=m.fs.ASU.O2_outlet,
                     destination=m.fs.ASU_O2_outlet.inlet)

    m.fs.ASU_O2_outlet.initialize()

    copy_port_values(source=m.fs.ASU_O2_outlet.outlet,
                     destination=m.fs.oxycombustor_translator.inlet)

    m.fs.oxycombustor_translator.initialize()

    copy_port_values(source=m.fs.anode_hx.shell_outlet,
                     destination=m.fs.combustor_mix.anode_inlet)

    copy_port_values(source=m.fs.oxycombustor_translator.outlet,
                     destination=m.fs.combustor_mix.oxygen_inlet)

    m.fs.combustor_mix.initialize()

    copy_port_values(source=m.fs.combustor_mix.outlet,
                     destination=m.fs.oxycombustor.inlet)

    m.fs.oxycombustor.initialize()

    copy_port_values(source=m.fs.oxycombustor.outlet,
                     destination=m.fs.anode_HRSG.inlet)

    m.fs.anode_HRSG.initialize()

    copy_port_values(source=m.fs.anode_HRSG.outlet,
                     destination=m.fs.CPU_translator.inlet)

    m.fs.CPU_translator.initialize()

    copy_port_values(source=m.fs.CPU_translator.outlet,
                     destination=m.fs.condenser.inlet)

    m.fs.condenser.initialize()

    copy_port_values(source=m.fs.condenser.outlet,
                     destination=m.fs.flash.inlet)

    m.fs.flash.initialize()

    copy_port_values(source=m.fs.flash.vap_outlet,
                     destination=m.fs.CPU.inlet)

    calculate_variable_from_constraint(m.fs.CPU.inlet.flow_mol[0],
                                       m.fs.CPU_inlet_F[0])

    calculate_variable_from_constraint(m.fs.CPU.inlet.pressure[0],
                                       m.fs.CPU_inlet_P[0])

    m.fs.CPU.initialize()


def SOFC_ROM_setup(m, use_DNN=True):
    # create the ROM
    if use_DNN:  # option to use deep neural net sofc surrogate model
        build_SOFC_ROM(m.fs)
    else:  # option to use algebraic sofc surrogate model
        build_SOFC_SM(m.fs)

    # build constraints connecting flowsheet to ROM input vars

    @m.fs.Constraint()
    def ROM_fuel_inlet_temperature(fs):
        return(fs.SOFC.fuel_temperature ==
               fs.anode_mix.feed_inlet.temperature[0]/pyunits.K - 273.15)

    @m.fs.Constraint()
    def ROM_air_inlet_temperature(fs):
        return(fs.SOFC.air_temperature ==
               fs.cathode.inlet.temperature[0]/pyunits.K - 273.15)

    @m.fs.Constraint()
    def ROM_air_recirculation(fs):
        return(fs.SOFC.air_recirculation ==
               fs.cathode_recycle.split_fraction[0, 'recycle_outlet'])

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
        full_O2_flow = (fs.anode_mix.feed_inlet.flow_mol[0] * (
            0.5 * fs.anode_mix.feed_inlet.mole_frac_comp[0, 'H2'] +
            0.5 * fs.anode_mix.feed_inlet.mole_frac_comp[0, 'CO'] +
            2.0 * fs.anode_mix.feed_inlet.mole_frac_comp[0, 'CH4'] +
            3.5 * fs.anode_mix.feed_inlet.mole_frac_comp[0, 'C2H6'] +
            5.0 * fs.anode_mix.feed_inlet.mole_frac_comp[0, 'C3H8'] +
            6.5 * fs.anode_mix.feed_inlet.mole_frac_comp[0, 'C4H10']))

        return (fs.SOFC.fuel_util*full_O2_flow ==
                fs.cathode.ion_outlet.flow_mol[0])

    @m.fs.Constraint()
    def ROM_air_utilization(fs):
        air_in = (fs.cathode_hx.tube_inlet.flow_mol[0] *
                  fs.cathode_hx.tube_inlet.mole_frac_comp[0, 'O2'])

        air_out = (fs.cathode_hx.shell_inlet.flow_mol[0] *
                   fs.cathode_hx.shell_inlet.mole_frac_comp[0, 'O2'])

        return fs.SOFC.air_util == 1 - air_out/air_in

    # unfix flowsheet vars and fix ROM inputs

    # current density
    m.fs.SOFC.current_density.fix(4000)

    # internal reformation fraction
    m.fs.SOFC.internal_reforming.fix(1)

    # air temperature - constrained by max cell temperature
    m.fs.cathode_hx.area.unfix()
    m.fs.SOFC.max_cell_temperature.fix(750)

    # air recirculation fraction
    m.fs.cathode_recycle.split_fraction[0, 'recycle_outlet'].unfix()
    m.fs.SOFC.air_recirculation.fix(0.5)

    # oxygen to carbon ratio
    m.fs.anode_recycle.split_fraction[0, 'recycle_outlet'].unfix()
    m.fs.SOFC.OTC.fix(2.1)

    # fuel utilization
    m.fs.cathode.ion_outlet.flow_mol.unfix()
    m.fs.SOFC.fuel_util.fix(0.85)

    # air utilization - constrained by deltaT across cell
    m.fs.air_blower.inlet.flow_mol.unfix()
    m.fs.SOFC.deltaT_cell.fix(93)

    # pressure
    m.fs.SOFC.pressure.fix(1)

    # initialize ROM
    calculate_variable_from_constraint(m.fs.SOFC.fuel_temperature,
                                       m.fs.ROM_fuel_inlet_temperature)
    calculate_variable_from_constraint(m.fs.SOFC.air_temperature,
                                       m.fs.ROM_air_inlet_temperature)
    calculate_variable_from_constraint(m.fs.SOFC.air_util,
                                       m.fs.ROM_air_utilization)

    if use_DNN:
        initialize_SOFC_ROM(m.fs.SOFC)
    else:
        initialize_SOFC_SM(m.fs.SOFC)

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
                fs.anode.outlet.temperature[0]/pyunits.K - 273.15)

    m.fs.anode.outlet.temperature.unfix()

    @m.fs.Constraint()
    def SOFC_energy_balance(fs):
        return (-1*pyunits.convert(fs.anode.heat_duty[0], pyunits.MW) ==
                fs.stack_power +
                pyunits.convert(fs.cathode_heat.heat_duty[0], pyunits.MW))

    m.fs.cathode_heat.outlet.temperature.unfix()


def add_result_constraints(m):
    # total heat supplied to HRSG
    m.fs.HRSG_heat_duty = pyo.Var(m.fs.time, initialize=300, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def HRSG_heat_duty_constraint(fs, t):
        return (-1*fs.HRSG_heat_duty[t] ==
                pyunits.convert(fs.anode_HRSG.heat_duty[t], pyunits.MW) +
                pyunits.convert(fs.cathode_HRSG.heat_duty[t], pyunits.MW))

    # steam required for ASU
    m.fs.ASU_HP_steam_heat = pyo.Var(m.fs.time, initialize=4, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def ASU_HP_steam_constraint(fs, t):
        return (fs.ASU_HP_steam_heat[t] ==
                0.005381*1e3*pyunits.MJ/pyunits.mol *
                fs.ASU.O2_outlet.flow_mol[t])

    # HRSG heat applied toward steam cycle
    m.fs.steam_cycle_heat = pyo.Var(m.fs.time,
                                    initialize=300, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def steam_cycle_heat_constraint(fs, t):
        return (fs.steam_cycle_heat[t] == fs.HRSG_heat_duty[t]
                - fs.ASU_HP_steam_heat[t])

    # power generated by steam cycle
    m.fs.steam_cycle_efficiency = pyo.Param(initialize=0.381, mutable=True)
    m.fs.steam_cycle_power = pyo.Var(m.fs.time,
                                     initialize=100, units=pyunits.MW)
    m.fs.steam_cycle_loss = pyo.Var(m.fs.time,
                                    initialize=200, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def steam_cycle_power_constraint(fs, t):
        return (fs.steam_cycle_power[t] ==
                fs.steam_cycle_heat[t]*fs.steam_cycle_efficiency)

    @m.fs.Constraint(m.fs.time)
    def steam_cycle_loss_constraint(fs, t):
        return (fs.steam_cycle_loss[t] ==
                fs.steam_cycle_heat[t]*(1 - fs.steam_cycle_efficiency))

    # stack AC power
    m.fs.stack_power_AC = pyo.Var(m.fs.time, initialize=550, units=pyunits.MW)
    m.fs.inverter_efficiency = pyo.Param(initialize=0.97, mutable=True)

    @m.fs.Constraint(m.fs.time)
    def stack_AC_power_constraint(fs, t):
        return fs.stack_power_AC[t] == fs.stack_power * fs.inverter_efficiency

    # gross plant power
    m.fs.gross_power = pyo.Var(m.fs.time, initialize=670, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def gross_power_constraint(fs, t):
        return (fs.gross_power[t] == fs.stack_power_AC[t] +
                fs.steam_cycle_power[t])

    # boiler feedwater pump load - steam cycle not modeled, scaled from BB
    m.fs.feedwater_pump_work = pyo.Var(m.fs.time,
                                       initialize=2, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def feedwater_pump_work_constraint(fs, t):
        ref_steam_cycle_power = 262.8*pyunits.MW
        ref_pump_work = 4.83*pyunits.MW
        return (fs.feedwater_pump_work[t]*ref_steam_cycle_power ==
                ref_pump_work*fs.steam_cycle_power[t])

    # condensate pump load - steam cycle not modeled, scaled from BB
    m.fs.condensate_pump_work = pyo.Var(m.fs.time,
                                        initialize=0.1, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def condensate_pump_work_constraint(fs, t):
        ref_steam_cycle_power = 262.8*pyunits.MW
        ref_pump_work = 0.15*pyunits.MW
        return (fs.condensate_pump_work[t]*ref_steam_cycle_power ==
                ref_pump_work*fs.steam_cycle_power[t])

    # steam turbine auxiliary load - steam cycle not modeled, scaled from BB
    m.fs.steam_turbine_auxiliary = pyo.Var(m.fs.time,
                                           initialize=0.1, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def steam_turbine_auxiliary_constraint(fs, t):
        ref_steam_cycle_power = 262.8*pyunits.MW
        ref_turbine_aux = 0.2*pyunits.MW
        return (fs.steam_turbine_auxiliary[t]*ref_steam_cycle_power ==
                ref_turbine_aux*fs.steam_cycle_power[t])

    # miscellaneous BOP - scaled from BB based on NG mass flow
    m.fs.misc_BOP_load = pyo.Var(m.fs.time, initialize=0.5, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def misc_BOP_load_constraint(fs, t):
        NG_flow = pyunits.convert(
            m.fs.anode_mix.feed_inlet_state[t].flow_mass,
            pyunits.lb/pyunits.hr)
        ref_NG_flow = 148095*pyunits.lb/pyunits.hr
        ref_load = 0.396*pyunits.MW
        return m.fs.misc_BOP_load[t] == ref_load*NG_flow/ref_NG_flow

    # calculate cooling water flowrate
    m.fs.cooling_water_duty = pyo.Var(m.fs.time,
                                      initialize=450, units=pyunits.MW)
    m.fs.cooling_water_flowrate = pyo.Var(m.fs.time, initialize=3000,
                                          units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def cooling_water_duty_constraint(fs, t):
        return (fs.cooling_water_duty[t] ==
                fs.steam_cycle_loss[t] +
                fs.CPU.heat_duty[t]/1e6*pyunits.MW +
                7.327*pyunits.MW +  # 25 MMbtu/hr of additional heat
                pyunits.convert(-1*fs.intercooler_s1.heat_duty[t] +
                                -1*fs.intercooler_s2.heat_duty[t] +
                                -1*fs.condenser.heat_duty[t],
                                pyunits.MW))

    @m.fs.Constraint(m.fs.time)
    def cooling_water_flowrate_constraint(fs, t):
        heat_capacity = 75.38*pyunits.J/pyunits.mol/pyunits.K
        delta_T = 36*pyunits.K
        molar_mass = 18*pyunits.g/pyunits.mol
        return (fs.cooling_water_flowrate[t] ==
                pyunits.convert(
                    fs.cooling_water_duty[t]*molar_mass/heat_capacity/delta_T,
                    pyunits.lb/pyunits.s))

    # circulating pump work - not modeled, scaled from NGFC pathways study
    m.fs.circulating_pump_work = pyo.Var(m.fs.time,
                                         initialize=1, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def circulating_pump_work_constraint(fs, t):
        ref_flowrate = 18602.6*pyunits.lb/pyunits.s
        ref_load = 2.778*pyunits.MW
        return (fs.circulating_pump_work[t]*ref_flowrate ==
                ref_load*fs.cooling_water_flowrate[t])

    # cooling tower fan load - not modeled, scaled from NGFC pathways study
    m.fs.cooling_tower_load = pyo.Var(m.fs.time,
                                      initialize=1, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def cooling_tower_load_constraint(fs, t):
        ref_flowrate = 18602.6*pyunits.lb/pyunits.s
        ref_load = 1.4372*pyunits.MW
        return (fs.cooling_tower_load[t]*ref_flowrate ==
                ref_load*fs.cooling_water_flowrate[t])

    # auxiliary load of the plant
    m.fs.auxiliary_load = pyo.Var(m.fs.time, initialize=60, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def auxiliary_load_constraint(fs, t):
        return (fs.auxiliary_load[t] == fs.CPU.work[t]/1e6 +
                fs.feedwater_pump_work[t] + fs.condensate_pump_work[t] +
                fs.steam_turbine_auxiliary[t] + fs.misc_BOP_load[t] +
                fs.circulating_pump_work[t] + fs.cooling_tower_load[t] +
                pyunits.convert(
                   (fs.air_blower.work_mechanical[t]/0.95 +
                    fs.cathode_blower.work_mechanical[t]/0.95 +
                    fs.anode_blower.work_mechanical[t]/0.95 +
                    fs.air_compressor_s1.work_mechanical[t]/0.96 +
                    fs.air_compressor_s2.work_mechanical[t]/0.96), pyunits.MW))

    # transformer losses
    m.fs.transformer_losses = pyo.Var(m.fs.time,
                                      initialize=2, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def transformer_losses_constraint(fs, t):
        HV_aux = pyunits.convert(fs.air_compressor_s1.work_mechanical[t]/0.96 +
                                 fs.air_compressor_s2.work_mechanical[t]/0.96,
                                 pyunits.MW) + fs.CPU.work[t]/1e6
        MV_aux = fs.auxiliary_load[t] - HV_aux
        LV_aux = MV_aux*0.15

        HV_gen_loss = ((fs.gross_power[t] - fs.auxiliary_load[t]) + HV_aux)*.003
        HV_loss = HV_aux*.003
        MV_loss = MV_aux*.005
        LV_loss = LV_aux*.005

        return (fs.transformer_losses[t] ==
                HV_gen_loss + HV_loss + MV_loss + LV_loss)

    # net plant power
    m.fs.net_power = pyo.Var(m.fs.time, initialize=660, units=pyunits.MW)

    @m.fs.Constraint(m.fs.time)
    def net_power_constraint(fs, t):
        return (fs.net_power[t] ==
                fs.gross_power[t] - fs.auxiliary_load[t] -
                fs.transformer_losses[t])

    # HHV efficiency
    m.fs.HHV_efficiency = pyo.Var(m.fs.time, initialize=0.6)

    @m.fs.Constraint(m.fs.time)
    def efficiency_rule(fs, t):
        NG_HHV = 908839.23*pyunits.J/pyunits.mol
        return (fs.HHV_efficiency[t] == fs.net_power[t] /
                pyunits.convert(
                    (NG_HHV * fs.anode_mix.feed_inlet.flow_mol[t]),
                    pyunits.MW))

    # CO2 emissions in g/kWh
    m.fs.CO2_emissions = pyo.Var(m.fs.time, initialize=300,
                                 units=pyunits.g/pyunits.kWh)

    @m.fs.Constraint(m.fs.time)
    def CO2_emission_constraint(fs, t):
        mass_flow = (fs.CPU.vent.flow_mol[t] *
                     fs.CPU.vent.mole_frac_comp[t, 'CO2'] *
                     44.01*pyunits.g/pyunits.mol)
        return fs.CO2_emissions[t] == (
            pyunits.convert((mass_flow/fs.net_power[t]),
                            (pyunits.g/pyunits.kWh)))


def initialize_results(m):

    variables = [
        m.fs.HRSG_heat_duty,
        m.fs.ASU_HP_steam_heat,
        m.fs.steam_cycle_heat,
        m.fs.steam_cycle_power,
        m.fs.steam_cycle_loss,
        m.fs.stack_power_AC,
        m.fs.gross_power,
        m.fs.feedwater_pump_work,
        m.fs.condensate_pump_work,
        m.fs.steam_turbine_auxiliary,
        m.fs.misc_BOP_load,
        m.fs.cooling_water_duty,
        m.fs.cooling_water_flowrate,
        m.fs.circulating_pump_work,
        m.fs.cooling_tower_load,
        m.fs.auxiliary_load,
        m.fs.transformer_losses,
        m.fs.net_power,
        m.fs.HHV_efficiency,
        m.fs.CO2_emissions]

    constraints = [
        m.fs.HRSG_heat_duty_constraint,
        m.fs.ASU_HP_steam_constraint,
        m.fs.steam_cycle_heat_constraint,
        m.fs.steam_cycle_power_constraint,
        m.fs.steam_cycle_loss_constraint,
        m.fs.stack_AC_power_constraint,
        m.fs.gross_power_constraint,
        m.fs.feedwater_pump_work_constraint,
        m.fs.condensate_pump_work_constraint,
        m.fs.steam_turbine_auxiliary_constraint,
        m.fs.misc_BOP_load_constraint,
        m.fs.cooling_water_duty_constraint,
        m.fs.cooling_water_flowrate_constraint,
        m.fs.circulating_pump_work_constraint,
        m.fs.cooling_tower_load_constraint,
        m.fs.auxiliary_load_constraint,
        m.fs.transformer_losses_constraint,
        m.fs.net_power_constraint,
        m.fs.efficiency_rule,
        m.fs.CO2_emission_constraint]

    for v, c in zip(variables, constraints):
        for t in m.fs.time:
            calculate_variable_from_constraint(v[t], c[t])


def add_costing(m):
    # fixed O&M costs
    NGFC_TPC = 693.019  # MM$

    get_fixed_OM_costs(m, 650, tech=6, fixed_TPC=NGFC_TPC)

    m.fs.costing.stack_replacement_cost = pyo.Var(
        initialize=13.26,
        bounds=(0, 100),
        doc="stack replacement cost in $MM/yr")
    m.fs.costing.stack_replacement_cost.fix()

    # redefine total fixed OM cost to include stack replacement cost
    del m.fs.costing.total_fixed_OM_cost_rule

    @m.fs.costing.Constraint()
    def total_fixed_OM_cost_rule(c):
        return c.total_fixed_OM_cost == (c.annual_operating_labor_cost +
                                         c.maintenance_labor_cost +
                                         c.admin_and_support_labor_cost +
                                         c.property_taxes_and_insurance +
                                         c.stack_replacement_cost)
    # variable O&M costs

    # Natural Gas Fuel
    NG_HHV = 908839.23*pyunits.J/pyunits.mol

    @m.fs.Expression(m.fs.time)
    def NG_energy_rate(fs, t):
        return m.fs.anode_mix.feed_inlet.flow_mol[t]*NG_HHV

    # Water
    m.fs.condensate_flowrate = pyo.Var(m.fs.time, units=pyunits.lb/pyunits.s)

    # scaled based on BB Case 31A steam turbine power
    @m.fs.Constraint(m.fs.time)
    def condensate_flowrate_constraint(fs, t):
        ref_flowrate = 388.586*pyunits.lb/pyunits.s
        ref_power = 262.8*pyunits.MW
        return (fs.condensate_flowrate[t] == ref_flowrate *
                fs.steam_cycle_power[t]/ref_power)

    m.fs.BFW_circulation_rate = pyo.Var(m.fs.time,
                                        units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def BFW_circulation_rate_constraint(fs, t):
        return fs.BFW_circulation_rate[t] == fs.condensate_flowrate[t]

    m.fs.BFW_makeup_flow = pyo.Var(m.fs.time, units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def BFW_makeup_flow_constraint(fs, t):
        return fs.BFW_makeup_flow[t] == 0.01*fs.BFW_circulation_rate[t]

    m.fs.evaporation_losses = pyo.Var(m.fs.time, units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def evaporation_losses_constraint(fs, t):
        return fs.evaporation_losses[t] == 0.016*fs.cooling_water_flowrate[t]

    m.fs.blowdown_losses = pyo.Var(m.fs.time, units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def blowdown_losses_constraint(fs, t):
        cycles_of_concentration = 4
        return (fs.blowdown_losses[t] ==
                fs.evaporation_losses[t]/(cycles_of_concentration - 1))

    m.fs.total_wet_losses = pyo.Var(m.fs.time, units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def wet_losses_constraint(fs, t):
        return fs.total_wet_losses[t] == (fs.evaporation_losses[t] +
                                          fs.blowdown_losses[t])

    m.fs.water_demand = pyo.Var(m.fs.time, units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def water_demand_constraint(fs, t):
        return fs.water_demand[t] == (fs.BFW_makeup_flow[t] +
                                      fs.total_wet_losses[t])

    m.fs.water_recycle = pyo.Var(m.fs.time, units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def water_recycle_constraint(fs, t):
        molar_mass = 18*pyunits.g/pyunits.mol
        CPU_water = pyunits.convert(
            fs.CPU.water.flow_mol[t]*molar_mass,
            pyunits.lb/pyunits.s)
        flash_water = pyunits.convert(
            fs.flash.liq_outlet.flow_mol[t]*molar_mass,
            pyunits.lb/pyunits.s)
        return fs.water_recycle[t] == CPU_water + flash_water

    m.fs.raw_water_withdrawal = pyo.Var(m.fs.time, units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def raw_water_withdrawal_constraint(fs, t):
        return fs.raw_water_withdrawal[t] == (fs.water_demand[t] -
                                              fs.water_recycle[t])

    m.fs.process_water_discharge = pyo.Var(m.fs.time, initialize=9,
                                           units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def water_discharge_constraint(fs, t):
        return fs.process_water_discharge[t] == 0.9*fs.blowdown_losses[t]

    m.fs.raw_water_consumption = pyo.Var(m.fs.time,
                                         units=pyunits.lb/pyunits.s)

    @m.fs.Constraint(m.fs.time)
    def raw_water_consumption_constraint(fs, t):
        return fs.raw_water_consumption[t] == (fs.raw_water_withdrawal[t] -
                                               fs.process_water_discharge[t])

    # MU & WT Chem
    @m.fs.Expression(m.fs.time)
    def waste_treatment_use(fs, t):
        use_rate = 0.00297886*pyunits.lb/pyunits.gal
        density = 8.34*pyunits.lb/pyunits.gal
        return fs.water_demand[t]/density*use_rate

    # NG desulfur TDA Adsorbent
    @m.fs.Expression(m.fs.time)
    def desulfur_adsorbent_use(fs, t):
        NG_sulfur_ppm = 5/1e6
        sulfur_MW = 32*pyunits.g/pyunits.mol
        sulfur_capacity = 0.03
        return (m.fs.anode_mix.feed_inlet.flow_mol[t] *
                NG_sulfur_ppm * sulfur_MW / sulfur_capacity)

    # Methanation Catalyst
    @m.fs.Expression(m.fs.time)
    def methanation_catalyst_use(fs, t):
        syngas_flow = pyunits.convert(
            m.fs.anode_mix.feed_inlet_state[0].flow_mass, pyunits.lb/pyunits.hr)
        return (3.2*syngas_flow/611167 *
                pyunits.m**3/pyunits.day*pyunits.hr/pyunits.lb)

    resources = ["natural gas",
                 "water treatment chemicals",
                 "desulfur adsorbent",
                 "methanation catalyst"]

    rates = [m.fs.NG_energy_rate,
             m.fs.waste_treatment_use,
             m.fs.desulfur_adsorbent_use,
             m.fs.methanation_catalyst_use]

    prices = {"desulfur adsorbent": 6.0297*pyunits.USD/pyunits.lb,
              "methanation catalyst": 601.765*pyunits.USD/pyunits.m**3}

    get_variable_OM_costs(m, m.fs.net_power, resources, rates, prices)


def initialize_costing(m):
    variables = [
        m.fs.condensate_flowrate,
        m.fs.BFW_circulation_rate,
        m.fs.BFW_makeup_flow,
        m.fs.evaporation_losses,
        m.fs.blowdown_losses,
        m.fs.total_wet_losses,
        m.fs.water_demand,
        m.fs.water_recycle,
        m.fs.raw_water_withdrawal,
        m.fs.process_water_discharge,
        m.fs.raw_water_consumption,
        m.fs.costing.other_variable_costs]

    constraints = [
        m.fs.condensate_flowrate_constraint,
        m.fs.BFW_circulation_rate_constraint,
        m.fs.BFW_makeup_flow_constraint,
        m.fs.evaporation_losses_constraint,
        m.fs.blowdown_losses_constraint,
        m.fs.wet_losses_constraint,
        m.fs.water_demand_constraint,
        m.fs.water_recycle_constraint,
        m.fs.raw_water_withdrawal_constraint,
        m.fs.water_discharge_constraint,
        m.fs.raw_water_consumption_constraint,
        # m.fs.stack_replacement
        ]

    for v, c in zip(variables, constraints):
        for t in m.fs.time:
            calculate_variable_from_constraint(v[t], c[t])

    initialize_fixed_OM_costs(m)
    initialize_variable_OM_costs(m)


def make_stream_dict(m):
    m._streams = OrderedDict(
        [
            ("NG_IN", m.fs.anode_mix.feed_inlet),
            ("ANO_HX_CI", m.fs.anode_hx.tube_inlet),
            ("ANO_HX_CO", m.fs.anode_hx.tube_outlet),
            ("ANODE_IN", m.fs.fuel_cell_mix.fuel_inlet),
            ("ANODE_OUT", m.fs.anode.outlet),
            ("ANODE_REC", m.fs.anode_blower.outlet),
            ("ANO_HX_HI", m.fs.anode_hx.shell_inlet),
            ("ANO_HX_HO", m.fs.anode_hx.shell_outlet),
            ("ASU_IN", m.fs.air_compressor_s1.inlet),
            ("COMB_O2", m.fs.ASU_O2_outlet.outlet),
            ("COMB_OUT", m.fs.oxycombustor.outlet),
            ("COND_IN", m.fs.anode_HRSG.outlet),
            ("COND_OUT", m.fs.condenser.outlet),
            ("AIR", m.fs.air_blower.inlet),
            ("CATH_HX_CI", m.fs.cathode_hx.tube_inlet),
            ("CATH_HX_CO", m.fs.cathode_hx.tube_outlet),
            ("CATH_IN", m.fs.cathode.inlet),
            ("CATH_OUT", m.fs.cathode_heat.outlet),
            ("CATH_REC", m.fs.cathode_blower.outlet),
            ("CATH_HX_HI", m.fs.cathode_hx.shell_inlet),
            ("CATH_HX_HO", m.fs.cathode_hx.shell_outlet),
            ("CATH_EXH", m.fs.cathode_HRSG.outlet)
            ])


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
        tags[i + "_F"] = fstr(df.loc[i, "Total Molar Flowrate"], 2, ' kmol/s')
        tags[i + "_T"] = fstr(df.loc[i, "Temperature"], 0, ' K')
        tags[i + "_P"] = fstr(df.loc[i, "Pressure"], 0, ' kPa')

    streams = {"WATER_1": m.fs.flash.liq_outlet,
               "VAP_OUT": m.fs.flash.vap_outlet}

    for i in streams.keys():
        tags[i + "_F"] = fstr(pyo.value(streams[i].flow_mol[0]), 2, ' kmol/s')
        tags[i + "_T"] = fstr(pyo.value(streams[i].temperature[0]), 0, ' K')
        tags[i + "_P"] = fstr(pyo.value(streams[i].pressure[0]), 0, ' kPa')

    CPU_streams = {"WATER_2": m.fs.CPU.water,
                   "CO2_PURE": m.fs.CPU.pureco2}

    for i in CPU_streams.keys():
        tags[i + "_F"] = fstr(pyo.value(CPU_streams[i].flow_mol[0])/1000, 2, ' kmol/s')
        tags[i + "_T"] = fstr(pyo.value(CPU_streams[i].temperature[0]), 0, ' K')
        tags[i + "_P"] = fstr(pyo.value(CPU_streams[i].pressure[0])/1000, 0, ' kPa')

    tags["stack_power_AC"] = fstr(pyo.value(m.fs.stack_power_AC[0]), 1)
    tags["HRSG_power"] = fstr(pyo.value(m.fs.steam_cycle_power[0]), 1)
    tags["aux_load"] = fstr(pyo.value(m.fs.auxiliary_load[0]), 1)
    tags["net_power"] = fstr(pyo.value(m.fs.net_power[0]), 1)
    tags["thermal_input"] = fstr(pyo.value(
        m.fs.net_power[0]/m.fs.HHV_efficiency[0]), 1)
    tags["efficiency"] = fstr(pyo.value(m.fs.HHV_efficiency[0])*100, 2)
    tags["emissions"] = fstr(pyo.value(m.fs.CO2_emissions[0]), 1)

    original_svg_file = os.path.join(this_file_dir(),
                                     "sofc_results_template.svg")
    with open(original_svg_file, "r") as f:
        svg_tag(tags, f, outfile=outfile)


def report_results(m, outfile="sofc_results_report.svg"):
    make_stream_dict(m)
    df = create_stream_table_dataframe(streams=m._streams, orient="index")
    pfd_result(outfile, m, df)


def partial_load_setup(m):
    # fix heat exchanger area
    m.fs.cathode_hx.area.fix()
    m.fs.SOFC.deltaT_cell.unfix()

    # turndown parameters
    max_ng_flow = 1.0946
    max_current_density = 4000

    # add constraint relating ng_flow and current density
    @m.fs.Constraint(m.fs.time)
    def NG_current_relation(fs, t):
        return (m.fs.anode_mix.feed_inlet.flow_mol[t]/max_ng_flow ==
                m.fs.SOFC.current_density/max_current_density)

    # set new power output
    m.fs.anode_mix.feed_inlet.flow_mol.unfix()
    m.fs.SOFC.current_density.unfix()

    m.fs.net_power.fix(650)


def check_scaling(m):
    import idaes.core.util.scaling as iscale
    jac, nlp = iscale.get_jacobian(m, scaled=True)
    # print("Extreme Jacobian entries:")
    sourceFile = open('extreme_jacobian.txt', 'w')
    for i in iscale.extreme_jacobian_entries(
            jac=jac, nlp=nlp, small=1e-6, large=1e3):
        print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]", file=sourceFile)
    sourceFile.close()

    # print("Unscaled constraints:")
    sourceFile2 = open('unscaled_constraints.txt', 'w')
    for c in iscale.unscaled_constraints_generator(m):
        print(f"    {c}", file=sourceFile2)
    sourceFile2.close()

    sourceFile3 = open('constraints_with_scale_factor.txt', 'w')
    # print("Scaled constraints by factor:")
    for c, s in iscale.constraints_with_scale_factor_generator(m):
        print(f"    {c}, {s}", file=sourceFile3)
    sourceFile3.close()

    # print("Badly scaled variables:")
    sourceFile4 = open('badly_scaled_var.txt', 'w')
    for v, sv in iscale.badly_scaled_var_generator(
            m, large=1e3, small=1e-4, zero=1e-12):
        print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}",
              file=sourceFile4)
    sourceFile.close()
    print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")


def get_model(use_DNN=True):
    # create model and flowsheet
    m = pyo.ConcreteModel()

    # if use_DNN:
    #     init_fname = "sofc_dnn_init.json.gz"
    # else:
    #     init_fname = "sofc_surrogate_init.json.gz"

    # if os.path.exists(init_fname):
    #     build_NGFC(m)
    #     SOFC_ROM_setup(m, use_DNN=use_DNN)
    #     add_SOFC_energy_balance(m)
    #     add_result_constraints(m)
    #     add_costing(m)
    #     ms.from_json(m, fname=init_fname)
    #     # scale_NGFC(m)

    # else:
    # build model
    build_NGFC(m)
    set_NGFC_inputs(m)
    # scale_NGFC(m)
    initialize_NGFC(m)
    SOFC_ROM_setup(m, use_DNN=use_DNN)
    add_SOFC_energy_balance(m)
    # solve model
    solver = pyo.SolverFactory("ipopt")
    solver.options = {'bound_push': 1e-23}
    solver.solve(m, tee=True)
    # and results and costing
    check_scaling(m)
    add_result_constraints(m)
    initialize_results(m)
    add_costing(m)
    initialize_costing(m)
    # # save model
    # ms.to_json(m, fname=init_fname)

    return m

# %%       # ------------------------------------------------------------------
if __name__ == "__main__":
    m = get_model()