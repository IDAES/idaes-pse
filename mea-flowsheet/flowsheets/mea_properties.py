import copy
from pyomo.environ import value, units as pyunits
from pyomo.common.config import ConfigValue, ConfigDict, In, Bool, ListOf
from idaes.core import LiquidPhase, Component, ControlVolume0DBlock, ControlVolume1DBlock, StateBlock, FlowDirection
import idaes.core.util.scaling as iscale
from idaes.models.properties.modular_properties import GenericParameterBlock
from idaes.models_extra.column_models.properties.MEA_solvent import (
    configuration as liquid_config,
)
from idaes.models_extra.column_models.properties.MEA_vapor import (
    wet_co2 as vapor_config, flue_gas as flue_gas
)

parmest_parameters = {
    'VLE': {'bic_k_eq_coeff_1': 366.061867998774,
            'bic_k_eq_coeff_2': -13326.25411,
            'bic_k_eq_coeff_3': -55.68643292,
            'car_k_eq_coeff_1': 164.039636,
            'car_k_eq_coeff_2': -707.0056712,
            'car_k_eq_coeff_3': -26.40136817,
            'lwm_coeff_1': -2.076073001,
            'lwm_coeff_2': 0.037322205,
            'lwm_coeff_3': -0.00032721,
            'lwm_coeff_4': -0.111102655,
            },
    'surface_tension': {'surf_tens_CO2_coeff_1': -0.00589934906112609,
                       'surf_tens_CO2_coeff_2': 0.00175020536428591,
                       'surf_tens_CO2_coeff_3': 0.129650182728177,
                       'surf_tens_CO2_coeff_4': 0.0000126444768126308,
                       'surf_tens_CO2_coeff_5': -5.73954817199691E-06,
                       'surf_tens_CO2_coeff_6': -0.00018969005534195,
                       'surf_tens_F_coeff_a': 1070.65668317975,
                       'surf_tens_F_coeff_b': -2578.78134208703,
                       'surf_tens_F_coeff_c': 3399.24113311222,
                       'surf_tens_F_coeff_d': -2352.47410135319,
                       'surf_tens_F_coeff_e': 2960.24753687833,
                       'surf_tens_F_coeff_f': 3.06684894924048,
                       'surf_tens_F_coeff_g': -1.79435372759593,
                       'surf_tens_F_coeff_h': -7.2124219075848,
                       'surf_tens_F_coeff_i': 2.97502322396621,
                       'surf_tens_F_coeff_j': -10.5738529301824,
                       },
    'molar_volume': {'vol_mol_liq_comp_coeff_a': -10.5792012186177,
                     'vol_mol_liq_comp_coeff_b': -2.02049415703576,
                     'vol_mol_liq_comp_coeff_c': 3.1506793296904,
                     'vol_mol_liq_comp_coeff_d': 192.012600751473,
                     'vol_mol_liq_comp_coeff_e': -695.384861676286,
                     },
    'viscosity': {'visc_d_coeff_a': -0.0854041877181552,
                  'visc_d_coeff_b': 2.72913373574306,
                  'visc_d_coeff_c': 35.1158892542595,
                  'visc_d_coeff_d': 1805.52759876533,
                  'visc_d_coeff_e': 0.00716025669867574,
                  'visc_d_coeff_f': 0.0106488402285381,
                  'visc_d_coeff_g': -0.0854041877181552,
                  },
    }

state_bounds_default = ConfigDict(
    description="Dictionary containing state variable name as key, tuple containing lower "
                "bound, initial value, upper bound, and units as value.",
)

state_bounds_default.declare(
    "flow_mol",
    ConfigValue(
        domain=tuple,
        default=(0, 1e3, 1e6, pyunits.mol / pyunits.s)
    )
)

state_bounds_default.declare(
    "temperature",
    ConfigValue(
        domain=tuple,
        default=(273.15, 298.15, 500, pyunits.K)
    )
)

state_bounds_default.declare(
    "pressure",
    ConfigValue(
        domain=tuple,
        default=(1, 1e5, 1e8, pyunits.Pa)
    )
)



def remove_ions(mea_dict: dict):
    mea_dict["components"]["H2O"]["type"] = Component

    mea_dict["components"]["CO2"]["type"] = Component
    del mea_dict["components"]["CO2"]["henry_component"]["Liq"]["basis"]

    del mea_dict["components"]["MEA_+"]
    del mea_dict["components"]["MEACOO_-"]
    del mea_dict["components"]["HCO3_-"]

    mea_dict["phases"]["Liq"]["type"] = LiquidPhase
    del mea_dict["phases"]["Liq"]["equation_of_state_options"]

    del mea_dict["state_components"]
    del mea_dict["inherent_reactions"]


def MEALiquidParameterBlock(ions=True, state_bounds=None):
    config = copy.deepcopy(liquid_config)
    if not ions:
        remove_ions(config)
    if state_bounds is None:
        state_bounds = state_bounds_default
    for key, tup in state_bounds.items():
        config["state_bounds"][key] = tup

    params = GenericParameterBlock(**config)

    return params


def MEAVaporParameterBlock(state_bounds=None):
    config = copy.deepcopy(vapor_config)
    if state_bounds is None:
        state_bounds = state_bounds_default
    for key, tup in state_bounds.items():
        config["state_bounds"][key] = tup

    params = GenericParameterBlock(**config)

    return params

def FlueGasParameterBlock(state_bounds=None):
    config = copy.deepcopy(flue_gas)
    if state_bounds is None:
        state_bounds = state_bounds_default
    for key, tup in state_bounds.items():
        config["state_bounds"][key] = tup

    params = GenericParameterBlock(**config)

    return params


def scale_mea_liquid_params(params, ions, scaling_factor_flow_mol=3e-4):
    params.set_default_scaling("enth_mol_phase", 3e-4)
    params.set_default_scaling("pressure", 1e-5)
    params.set_default_scaling("temperature", 1)
    params.set_default_scaling("flow_mol", scaling_factor_flow_mol)
    params.set_default_scaling("flow_mol_phase", scaling_factor_flow_mol)

    params.set_default_scaling("flow_mass_phase", scaling_factor_flow_mol / 24e-3)  # MW mixture ~= 24 g/Mol
    params.set_default_scaling("dens_mol_phase", 1 / 43000)
    params.set_default_scaling("visc_d_phase", 700)
    params.set_default_scaling("log_k_eq", 1)

    mole_frac_scaling_factors = {
        "H2O": 1,
        "MEA": 10,
        "CO2": 20,
    }
    if ions:
        mole_frac_true_scaling_factors = {
            "CO2": 1e4,  # Could go to 1e4 or 3e4
            "H2O": 1,
            "HCO3_-": 1000,
            "MEA": 30,
            "MEACOO_-": 30,
            "MEA_+": 30,
        }
    for comp, sf_x in mole_frac_scaling_factors.items():
        params.set_default_scaling("mole_frac_comp", sf_x, index=comp)
        params.set_default_scaling("mole_frac_phase_comp", sf_x, index=("Liq", comp))
        params.set_default_scaling(
            "flow_mol_phase_comp",
            sf_x * scaling_factor_flow_mol,
            index=("Liq", comp)
        )

    if ions:
        for comp, sf_x in mole_frac_true_scaling_factors.items():
            params.set_default_scaling("mole_frac_phase_comp_true", sf_x, index=("Liq", comp))
            params.set_default_scaling(
                "flow_mol_phase_comp_true",
                sf_x * scaling_factor_flow_mol,
                index=("Liq", comp)
            )

        params.set_default_scaling("apparent_inherent_reaction_extent", 2e-2, index="bicarbonate")
        params.set_default_scaling("apparent_inherent_reaction_extent", 1e-3, index="carbamate")


def scale_mea_vapor_params(params, scaling_factor_flow_mol=3e-4):
    params.set_default_scaling("enth_mol_phase", 3e-4)
    params.set_default_scaling("pressure", 1e-5)
    params.set_default_scaling("temperature", 1)
    params.set_default_scaling("flow_mol", scaling_factor_flow_mol)
    params.set_default_scaling("flow_mol_phase", scaling_factor_flow_mol)

    params.set_default_scaling("flow_mol_phase_comp", 2 * scaling_factor_flow_mol)
    params.set_default_scaling("flow_mass_phase", scaling_factor_flow_mol / 24e-3)  # Say MW ~=24 g/mol
    params.set_default_scaling("visc_d_phase", 6e4)

def switch_liquid_to_parmest_params(params, ions):
    if ions:
        params.reaction_bicarbonate.k_eq_coeff_1.set_value(parmest_parameters['VLE']['bic_k_eq_coeff_1'])
        params.reaction_bicarbonate.k_eq_coeff_2.set_value(parmest_parameters['VLE']['bic_k_eq_coeff_2'])
        params.reaction_bicarbonate.k_eq_coeff_3.set_value(parmest_parameters['VLE']['bic_k_eq_coeff_3'])
        params.reaction_carbamate.k_eq_coeff_1.set_value(parmest_parameters['VLE']['car_k_eq_coeff_1'])
        params.reaction_carbamate.k_eq_coeff_2.set_value(parmest_parameters['VLE']['car_k_eq_coeff_2'])
        params.reaction_carbamate.k_eq_coeff_3.set_value(parmest_parameters['VLE']['car_k_eq_coeff_3'])
    params.CO2.lwm_coeff_1.set_value(parmest_parameters['VLE']['lwm_coeff_1'])
    params.CO2.lwm_coeff_2.set_value(parmest_parameters['VLE']['lwm_coeff_2'])
    params.CO2.lwm_coeff_3.set_value(parmest_parameters['VLE']['lwm_coeff_3'])
    params.CO2.lwm_coeff_4.set_value(parmest_parameters['VLE']['lwm_coeff_4'])

    params.Liq.surf_tens_CO2_coeff_1.set_value(
        parmest_parameters['surface_tension']['surf_tens_CO2_coeff_1'])
    params.Liq.surf_tens_CO2_coeff_2.set_value(
        parmest_parameters['surface_tension']['surf_tens_CO2_coeff_2'])
    params.Liq.surf_tens_CO2_coeff_3.set_value(
        parmest_parameters['surface_tension']['surf_tens_CO2_coeff_3'])
    params.Liq.surf_tens_CO2_coeff_4.set_value(
        parmest_parameters['surface_tension']['surf_tens_CO2_coeff_4'])
    params.Liq.surf_tens_CO2_coeff_5.set_value(
        parmest_parameters['surface_tension']['surf_tens_CO2_coeff_5'])
    params.Liq.surf_tens_CO2_coeff_6.set_value(
        parmest_parameters['surface_tension']['surf_tens_CO2_coeff_6'])
    params.Liq.surf_tens_F_coeff_a.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_a'])
    params.Liq.surf_tens_F_coeff_b.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_b'])
    params.Liq.surf_tens_F_coeff_c.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_c'])
    params.Liq.surf_tens_F_coeff_d.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_d'])
    params.Liq.surf_tens_F_coeff_e.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_e'])
    params.Liq.surf_tens_F_coeff_f.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_f'])
    params.Liq.surf_tens_F_coeff_g.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_g'])
    params.Liq.surf_tens_F_coeff_h.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_h'])
    params.Liq.surf_tens_F_coeff_i.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_i'])
    params.Liq.surf_tens_F_coeff_j.set_value(parmest_parameters['surface_tension']['surf_tens_F_coeff_j'])

    params.CO2.vol_mol_liq_comp_coeff_a.set_value(
        parmest_parameters['molar_volume']['vol_mol_liq_comp_coeff_a'])
    params.MEA.vol_mol_liq_comp_coeff_b.set_value(
        parmest_parameters['molar_volume']['vol_mol_liq_comp_coeff_b'])
    params.MEA.vol_mol_liq_comp_coeff_c.set_value(
        parmest_parameters['molar_volume']['vol_mol_liq_comp_coeff_c'])
    params.CO2.vol_mol_liq_comp_coeff_d.set_value(
        parmest_parameters['molar_volume']['vol_mol_liq_comp_coeff_d'])
    params.CO2.vol_mol_liq_comp_coeff_e.set_value(
        parmest_parameters['molar_volume']['vol_mol_liq_comp_coeff_e'])

    params.Liq.visc_d_coeff_a.set_value(parmest_parameters['viscosity']['visc_d_coeff_a'])
    params.Liq.visc_d_coeff_b.set_value(parmest_parameters['viscosity']['visc_d_coeff_b'])
    params.Liq.visc_d_coeff_c.set_value(parmest_parameters['viscosity']['visc_d_coeff_c'])
    params.Liq.visc_d_coeff_d.set_value(parmest_parameters['viscosity']['visc_d_coeff_d'])
    params.Liq.visc_d_coeff_e.set_value(parmest_parameters['viscosity']['visc_d_coeff_e'])
    params.Liq.visc_d_coeff_f.set_value(parmest_parameters['viscosity']['visc_d_coeff_f'])
    params.Liq.visc_d_coeff_g.set_value(parmest_parameters['viscosity']['visc_d_coeff_g'])

def switch_liquid_to_literature_params(params, ions):
    if ions:
        inherent_rxns = liquid_config["inherent_reactions"]
        bic_k_eq_coeff = inherent_rxns["bicarbonate"]["parameter_data"]["k_eq_coeff"]
        params.reaction_bicarbonate.k_eq_coeff_1.set_value(bic_k_eq_coeff["1"])
        params.reaction_bicarbonate.k_eq_coeff_2.set_value(bic_k_eq_coeff["2"])
        params.reaction_bicarbonate.k_eq_coeff_3.set_value(bic_k_eq_coeff["3"])
        car_k_eq_coeff = inherent_rxns["carbamate"]["parameter_data"]["k_eq_coeff"]
        params.reaction_carbamate.k_eq_coeff_1.set_value(car_k_eq_coeff["1"])
        params.reaction_carbamate.k_eq_coeff_2.set_value(car_k_eq_coeff["2"])
        params.reaction_carbamate.k_eq_coeff_3.set_value(car_k_eq_coeff["3"])

    CO2_param_data = liquid_config["components"]["CO2"]["parameter_data"]
    lwm_coeff = CO2_param_data["lwn_coeff"]
    params.CO2.lwm_coeff_1.set_value(lwm_coeff["1"])
    params.CO2.lwm_coeff_2.set_value(lwm_coeff["2"])
    params.CO2.lwm_coeff_3.set_value(lwm_coeff["3"])
    params.CO2.lwm_coeff_4.set_value(lwm_coeff["4"])

    vol_mol_liq_comp_coeff = CO2_param_data["vol_mol_liq_comp_coeff"]
    params.CO2.vol_mol_liq_comp_coeff_a.set_value(vol_mol_liq_comp_coeff["a"])
    params.CO2.vol_mol_liq_comp_coeff_b.set_value(vol_mol_liq_comp_coeff["b"])
    params.CO2.vol_mol_liq_comp_coeff_c.set_value(vol_mol_liq_comp_coeff["c"])
    params.CO2.vol_mol_liq_comp_coeff_d.set_value(vol_mol_liq_comp_coeff["d"])
    params.CO2.vol_mol_liq_comp_coeff_e.set_value(vol_mol_liq_comp_coeff["e"])

    liq_param_data = liquid_config["phases"]["Liq"]["parameter_data"]
    surf_tens_CO2_coeff = liq_param_data["surf_tens_CO2_coeff"]
    params.Liq.surf_tens_CO2_coeff_1.set_value(surf_tens_CO2_coeff["1"])
    params.Liq.surf_tens_CO2_coeff_2.set_value(surf_tens_CO2_coeff["2"])
    params.Liq.surf_tens_CO2_coeff_3.set_value(surf_tens_CO2_coeff["3"])
    params.Liq.surf_tens_CO2_coeff_4.set_value(surf_tens_CO2_coeff["4"])
    params.Liq.surf_tens_CO2_coeff_5.set_value(surf_tens_CO2_coeff["5"])
    params.Liq.surf_tens_CO2_coeff_6.set_value(surf_tens_CO2_coeff["6"])

    surf_tens_F_coeff = liq_param_data["surf_tens_F_coeff"]
    params.Liq.surf_tens_F_coeff_a.set_value(surf_tens_F_coeff["a"])
    params.Liq.surf_tens_F_coeff_b.set_value(surf_tens_F_coeff["b"])
    params.Liq.surf_tens_F_coeff_c.set_value(surf_tens_F_coeff["c"])
    params.Liq.surf_tens_F_coeff_d.set_value(surf_tens_F_coeff["d"])
    params.Liq.surf_tens_F_coeff_e.set_value(surf_tens_F_coeff["e"])
    params.Liq.surf_tens_F_coeff_f.set_value(surf_tens_F_coeff["f"])
    params.Liq.surf_tens_F_coeff_g.set_value(surf_tens_F_coeff["g"])
    params.Liq.surf_tens_F_coeff_h.set_value(surf_tens_F_coeff["h"])
    params.Liq.surf_tens_F_coeff_i.set_value(surf_tens_F_coeff["i"])
    params.Liq.surf_tens_F_coeff_j.set_value(surf_tens_F_coeff["j"])

    visc_d_coeff = liq_param_data["visc_d_coeff"]
    params.Liq.visc_d_coeff_a.set_value(visc_d_coeff["a"])
    params.Liq.visc_d_coeff_b.set_value(visc_d_coeff["b"])
    params.Liq.visc_d_coeff_c.set_value(visc_d_coeff["c"])
    params.Liq.visc_d_coeff_d.set_value(visc_d_coeff["d"])
    params.Liq.visc_d_coeff_e.set_value(visc_d_coeff["e"])
    params.Liq.visc_d_coeff_f.set_value(visc_d_coeff["f"])
    params.Liq.visc_d_coeff_g.set_value(visc_d_coeff["g"])

def initialize_inherent_reactions(blk):
    if isinstance(blk, ControlVolume1DBlock):
        if blk._flow_direction == FlowDirection.forward:
            source_idx = blk.length_domain.first()
        else:
            source_idx = blk.length_domain.last()
        source = blk.properties[blk.flowsheet().time.first(), source_idx]

        x_lim_reag = min(
            value(source.mole_frac_comp["MEA"]),
            value(source.mole_frac_comp["CO2"]),
        )

        F_lim_reag = value(source.flow_mol*x_lim_reag)

        for k in blk.properties.keys():
            blk.properties[k].apparent_inherent_reaction_extent["bicarbonate"].value = 0.2 * F_lim_reag
            blk.properties[k].apparent_inherent_reaction_extent["carbamate"].value = 0.1 * F_lim_reag
    elif isinstance(blk, ControlVolume0DBlock):
        source = blk.properties_in[blk.flowsheet().time.first()]
        x_lim_reag = min(
            value(source.mole_frac_comp["MEA"]),
            value(source.mole_frac_comp["CO2"]),
        )

        F_lim_reag = value(source.flow_mol*x_lim_reag)

        for t in blk.flowsheet().time:
            blk.properties_in[t].apparent_inherent_reaction_extent["bicarbonate"].value = 0.2 * F_lim_reag
            blk.properties_in[t].apparent_inherent_reaction_extent["carbamate"].value = 0.1 * F_lim_reag

            blk.properties_out[t].apparent_inherent_reaction_extent["bicarbonate"].value = 0.2 * F_lim_reag
            blk.properties_out[t].apparent_inherent_reaction_extent["carbamate"].value = 0.1 * F_lim_reag

    elif isinstance(blk, StateBlock):
        for sub_blk in blk.values():
            x_lim_reag = min(
                value(sub_blk.mole_frac_comp["MEA"]),
                value(sub_blk.mole_frac_comp["CO2"]),
            )

            F_lim_reag = value(sub_blk.flow_mol*x_lim_reag)


            sub_blk.apparent_inherent_reaction_extent["bicarbonate"].value = 0.2 * F_lim_reag
            sub_blk.apparent_inherent_reaction_extent["carbamate"].value = 0.1 * F_lim_reag
    else:
        raise ValueError("Must be passed control volume or state block")