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
Natural gas property package for the vapor phase using Peng-Robinson equation
of state.
"""
# Import Python libraries
import logging
import copy

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import VaporPhase, Component

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType

from idaes.generic_models.properties.core.pure import NIST
from idaes.generic_models.properties.core.pure import RPP4

from idaes.generic_models.properties.core.reactions.dh_rxn import \
    constant_dh_rxn
from idaes.generic_models.properties.core.reactions.rate_constant import \
    arrhenius
from idaes.generic_models.properties.core.reactions.rate_forms import \
    power_law_rate
from idaes.generic_models.properties.core.generic.generic_reaction import (
    ConcentrationForm)

# Set up logger
_log = logging.getLogger(__name__)

# Property Sources

# Source: NIST webbook
# Properties: Heat capacity coefficients for all species except ethane,
# propane, and butane. Reference enthalpies and entropies for all species.

# Source: The Properties of Gases and Liquids (1987)
# 4th edition, Chemical Engineering Series - Robert C. Reid
# Properties: Critical temperatures and pressures. Omega.
# Heat capacity coefficients for ethane, propane, and butane.


# returns a configuration dictionary for the list of specified components
def get_NG_properties(components='all'):
    configuration = {
        'components': {
            'H2': {'type': Component,
                   "elemental_composition": {"H": 2},
                   'enth_mol_ig_comp': NIST,
                   'entr_mol_ig_comp': NIST,
                   'parameter_data': {
                       'mw': (0.0020159, pyunits.kg/pyunits.mol),
                       'pressure_crit': (13e5, pyunits.Pa),
                       'temperature_crit': (33.2, pyunits.K),
                       'omega': -0.218,
                       'cp_mol_ig_comp_coeff': {'A': 33.066178,
                                                'B': -11.363417,
                                                'C': 11.432816,
                                                'D': -2.772874,
                                                'E': -0.158558,
                                                'F': -9.980797,
                                                'G': 172.707974,
                                                'H': 0.0},
                       'enth_mol_form_vap_comp_ref': (
                           0, pyunits.J/pyunits.mol),
                       'entr_mol_form_vap_comp_ref': (
                           130.7, pyunits.J/pyunits.mol/pyunits.K)}},

            'CO': {'type': Component,
                   "elemental_composition": {"C": 1, "O": 1},
                   'enth_mol_ig_comp': NIST,
                   'entr_mol_ig_comp': NIST,
                   'parameter_data': {
                       'mw': (0.0280101, pyunits.kg/pyunits.mol),
                       'pressure_crit': (35e5, pyunits.Pa),
                       'temperature_crit': (132.9, pyunits.K),
                       'omega': 0.066,
                       'cp_mol_ig_comp_coeff': {'A': 25.56759,
                                                'B': 6.09613,
                                                'C': 4.054656,
                                                'D': -2.671301,
                                                'E': 0.131021,
                                                'F': -118.0089,
                                                'G': 227.3665,
                                                'H': -110.5271},
                       'enth_mol_form_vap_comp_ref': (
                           -110530, pyunits.J/pyunits.mol),
                       'entr_mol_form_vap_comp_ref': (
                           197.7, pyunits.J/pyunits.mol/pyunits.K)}},

            'H2O': {'type': Component,
                    "elemental_composition": {"H": 2, "O": 1},
                    'enth_mol_ig_comp': NIST,
                    'entr_mol_ig_comp': NIST,
                    'parameter_data': {
                        'mw': (0.01801528, pyunits.kg/pyunits.mol),
                        'pressure_crit': (221.2e5, pyunits.Pa),
                        'temperature_crit': (647.3, pyunits.K),
                        'omega': 0.344,
                        'cp_mol_ig_comp_coeff': {'A': 30.092,
                                                 'B': 6.832514,
                                                 'C': 6.793435,
                                                 'D': -2.53448,
                                                 'E': 0.082139,
                                                 'F': -250.881,
                                                 'G': 223.3967,
                                                 'H': -241.8264},
                        'enth_mol_form_vap_comp_ref': (
                            -241826, pyunits.J/pyunits.mol),
                        'entr_mol_form_vap_comp_ref': (
                            188.8, pyunits.J/pyunits.mol/pyunits.K)}},

            'CO2': {'type': Component,
                    "elemental_composition": {"C": 1, "O": 2},
                    'enth_mol_ig_comp': NIST,
                    'entr_mol_ig_comp': NIST,
                    'parameter_data': {
                        'mw': (0.04401, pyunits.kg/pyunits.mol),
                        'pressure_crit': (73.8e5, pyunits.Pa),
                        'temperature_crit': (304.1, pyunits.K),
                        'omega': 0.239,
                        'cp_mol_ig_comp_coeff': {'A': 24.99735,
                                                 'B': 55.18696,
                                                 'C': -33.69137,
                                                 'D': 7.948387,
                                                 'E': -0.136638,
                                                 'F': -403.6075,
                                                 'G': 228.2431,
                                                 'H': -393.5224},
                        'enth_mol_form_vap_comp_ref': (
                            -393510, pyunits.J/pyunits.mol),
                        'entr_mol_form_vap_comp_ref': (
                            213.8, pyunits.J/pyunits.mol/pyunits.K)}},

            'O2': {'type': Component,
                   "elemental_composition": {"O": 2},
                   'enth_mol_ig_comp': NIST,
                   'entr_mol_ig_comp': NIST,
                   'parameter_data': {
                       'mw': (0.031998, pyunits.kg/pyunits.mol),
                       'pressure_crit': (50.4e5, pyunits.Pa),
                       'temperature_crit': (154.6, pyunits.K),
                       'omega': 0.025,
                       'cp_mol_ig_comp_coeff': {'A': 30.03235,
                                                'B': 8.772972,
                                                'C': -3.988133,
                                                'D': 0.788313,
                                                'E': -0.741599,
                                                'F': -11.32468,
                                                'G': 236.1663,
                                                'H': 0.0},
                       'enth_mol_form_vap_comp_ref': (
                           0, pyunits.J/pyunits.mol),
                       'entr_mol_form_vap_comp_ref': (
                           -205.2, pyunits.J/pyunits.mol/pyunits.K)}},

            'N2': {'type': Component,
                   "elemental_composition": {"N": 2},
                   'enth_mol_ig_comp': NIST,
                   'entr_mol_ig_comp': NIST,
                   'parameter_data': {
                       'mw': (0.0280134, pyunits.kg/pyunits.mol),
                       'pressure_crit': (33.9e5, pyunits.Pa),
                       'temperature_crit': (126.2, pyunits.K),
                       'omega': 0.039,
                       'cp_mol_ig_comp_coeff': {'A': 19.50583,
                                                'B': 19.88705,
                                                'C': -8.598535,
                                                'D': 1.369784,
                                                'E': 0.527601,
                                                'F': -4.935202,
                                                'G': 212.39,
                                                'H': 0.0},
                       'enth_mol_form_vap_comp_ref': (
                           0, pyunits.J/pyunits.mol),
                       'entr_mol_form_vap_comp_ref': (
                           -191.6, pyunits.J/pyunits.mol/pyunits.K)}},

            'Ar': {'type': Component,
                   "elemental_composition": {"Ar": 1},
                   'enth_mol_ig_comp': NIST,
                   'entr_mol_ig_comp': NIST,
                   'parameter_data': {
                       'mw': (0.039948, pyunits.kg/pyunits.mol),
                       'pressure_crit': (48.7e5, pyunits.Pa),
                       'temperature_crit': (150.8, pyunits.K),
                       'omega': 0.001,
                       'cp_mol_ig_comp_coeff': {'A': 20.786,
                                                'B': 2.825911e-07,
                                                'C': -1.464191e-07,
                                                'D': 1.092131e-08,
                                                'E': -3.661371e-08,
                                                'F': -6.19735,
                                                'G': 179.999,
                                                'H': 0.0},
                       'enth_mol_form_vap_comp_ref': (
                           0, pyunits.J/pyunits.mol),
                       'entr_mol_form_vap_comp_ref': (
                           -154.9, pyunits.J/pyunits.mol/pyunits.K)}},

            'CH4': {'type': Component,
                    "elemental_composition": {"C": 1, "H": 4},
                    'enth_mol_ig_comp': NIST,
                    'entr_mol_ig_comp': NIST,
                    'parameter_data': {
                        'mw': (0.0160425, pyunits.kg/pyunits.mol),
                        'pressure_crit': (46e5, pyunits.Pa),
                        'temperature_crit': (190.4, pyunits.K),
                        'omega': 0.011,
                        'cp_mol_ig_comp_coeff': {'A': -0.703029,
                                                 'B': 108.4773,
                                                 'C': -42.52157,
                                                 'D': 5.862788,
                                                 'E': 0.678565,
                                                 'F': -76.84376,
                                                 'G': 158.7163,
                                                 'H': -74.8731},
                        'enth_mol_form_vap_comp_ref': (
                            -74870, pyunits.J/pyunits.mol),
                        'entr_mol_form_vap_comp_ref': (
                            -186.3, pyunits.J/pyunits.mol/pyunits.K)}},

            'C2H6': {'type': Component,
                     "elemental_composition": {"C": 2, "H": 6},
                     'enth_mol_ig_comp': RPP4,
                     'entr_mol_ig_comp': RPP4,
                     'parameter_data': {
                         'mw': (0.030069, pyunits.kg/pyunits.mol),
                         'pressure_crit': (48.8e5, pyunits.Pa),
                         'temperature_crit': (305.4, pyunits.K),
                         'omega': 0.099,
                         'cp_mol_ig_comp_coeff': {'A': 5.409,
                                                  'B': 1.781e-1,
                                                  'C': -6.938e-5,
                                                  'D': 8.713e-9},
                         'enth_mol_form_vap_comp_ref': (
                             -84000, pyunits.J/pyunits.mol),
                         'entr_mol_form_vap_comp_ref': (
                             229.2, pyunits.J/pyunits.mol/pyunits.K)}},

            'C3H8': {'type': Component,
                     "elemental_composition": {"C": 3, "H": 8},
                     'enth_mol_ig_comp': RPP4,
                     'entr_mol_ig_comp': RPP4,
                     'parameter_data': {
                         'mw': (0.0320849, pyunits.kg/pyunits.mol),
                         'pressure_crit': (42.5e5, pyunits.Pa),
                         'temperature_crit': (369.8, pyunits.K),
                         'omega': 0.153,
                         'cp_mol_ig_comp_coeff': {'A': -4.224,
                                                  'B': 3.063e-1,
                                                  'C': -1.586e-4,
                                                  'D': 3.215e-8},
                         'enth_mol_form_vap_comp_ref': (
                             -104700, pyunits.J/pyunits.mol),
                         'entr_mol_form_vap_comp_ref': (
                             270.3, pyunits.J/pyunits.mol/pyunits.K)}},

            'C4H10': {'type': Component,
                      "elemental_composition": {"C": 4, "H": 10},
                      'enth_mol_ig_comp': RPP4,
                      'entr_mol_ig_comp': RPP4,
                      'parameter_data': {
                          'mw': (0.058122, pyunits.kg/pyunits.mol),
                          'pressure_crit': (38e5, pyunits.Pa),
                          'temperature_crit': (425.2, pyunits.K),
                          'omega': 0.199,
                          'cp_mol_ig_comp_coeff': {'A': 9.487,
                                                   'B': 3.313e-1,
                                                   'C': -1.108e-4,
                                                   'D': -2.822e-9},
                          'enth_mol_form_vap_comp_ref': (
                              -125600, pyunits.J/pyunits.mol),
                          'entr_mol_form_vap_comp_ref': (
                              310.0, pyunits.J/pyunits.mol/pyunits.K)}}},

        # Specifying phases
        "phases":  {'Vap': {"type": VaporPhase,
                            "equation_of_state": Cubic,
                            "equation_of_state_options":
                                {"type": CubicType.PR}}},

        # Set base units of measurement
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},

        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {"flow_mol": (0, 8000, 50000, pyunits.mol/pyunits.s),
                         "temperature": (273.15, 500, 2500, pyunits.K),
                         "pressure": (5e4, 1.3e5, 1e7, pyunits.Pa)},
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),

        "parameter_data": {"PR_kappa":
                           {('H2', 'H2'): 0, ('H2', 'CO'): 0,
                            ('H2', 'H2O'): 0, ('H2', 'CO2'): 0,
                            ('H2', 'O2'): 0, ('H2', 'N2'): 0,
                            ('H2', 'Ar'): 0, ('H2', 'CH4'): 0,
                            ('H2', 'C2H6'): 0, ('H2', 'C3H8'): 0,
                            ('H2', 'C4H10'): 0, ('CO', 'H2'): 0,
                            ('CO', 'CO'): 0, ('CO', 'H2O'): 0,
                            ('CO', 'CO2'): 0, ('CO', 'O2'): 0,
                            ('CO', 'N2'): 0, ('CO', 'Ar'): 0,
                            ('CO', 'CH4'): 0, ('CO', 'C2H6'): 0,
                            ('CO', 'C3H8'): 0, ('CO', 'C4H10'): 0,
                            ('H2O', 'H2'): 0, ('H2O', 'CO'): 0,
                            ('H2O', 'H2O'): 0, ('H2O', 'CO2'): 0,
                            ('H2O', 'O2'): 0, ('H2O', 'N2'): 0,
                            ('H2O', 'Ar'): 0, ('H2O', 'CH4'): 0,
                            ('H2O', 'C2H6'): 0, ('H2O', 'C3H8'): 0,
                            ('H2O', 'C4H10'): 0, ('CO2', 'H2'): 0,
                            ('CO2', 'CO'): 0, ('CO2', 'H2O'): 0,
                            ('CO2', 'CO2'): 0, ('CO2', 'O2'): 0,
                            ('CO2', 'N2'): 0, ('CO2', 'Ar'): 0,
                            ('CO2', 'CH4'): 0, ('CO2', 'C2H6'): 0,
                            ('CO2', 'C3H8'): 0, ('CO2', 'C4H10'): 0,
                            ('O2', 'H2'): 0, ('O2', 'CO'): 0,
                            ('O2', 'H2O'): 0, ('O2', 'CO2'): 0,
                            ('O2', 'O2'): 0, ('O2', 'N2'): 0,
                            ('O2', 'Ar'): 0, ('O2', 'CH4'): 0,
                            ('O2', 'C2H6'): 0, ('O2', 'C3H8'): 0,
                            ('O2', 'C4H10'): 0, ('N2', 'H2'): 0,
                            ('N2', 'CO'): 0, ('N2', 'H2O'): 0,
                            ('N2', 'CO2'): 0, ('N2', 'O2'): 0,
                            ('N2', 'N2'): 0, ('N2', 'Ar'): 0,
                            ('N2', 'CH4'): 0, ('N2', 'C2H6'): 0,
                            ('N2', 'C3H8'): 0, ('N2', 'C4H10'): 0,
                            ('Ar', 'H2'): 0, ('Ar', 'CO'): 0,
                            ('Ar', 'H2O'): 0, ('Ar', 'CO2'): 0,
                            ('Ar', 'O2'): 0, ('Ar', 'N2'): 0,
                            ('Ar', 'Ar'): 0, ('Ar', 'CH4'): 0,
                            ('Ar', 'C2H6'): 0, ('Ar', 'C3H8'): 0,
                            ('Ar', 'C4H10'): 0, ('CH4', 'H2'): 0,
                            ('CH4', 'CO'): 0, ('CH4', 'H2O'): 0,
                            ('CH4', 'CO2'): 0, ('CH4', 'O2'): 0,
                            ('CH4', 'N2'): 0, ('CH4', 'Ar'): 0,
                            ('CH4', 'CH4'): 0, ('CH4', 'C2H6'): 0,
                            ('CH4', 'C3H8'): 0, ('CH4', 'C4H10'): 0,
                            ('C2H6', 'H2'): 0, ('C2H6', 'CO'): 0,
                            ('C2H6', 'H2O'): 0, ('C2H6', 'CO2'): 0,
                            ('C2H6', 'O2'): 0, ('C2H6', 'N2'): 0,
                            ('C2H6', 'Ar'): 0, ('C2H6', 'CH4'): 0,
                            ('C2H6', 'C2H6'): 0, ('C2H6', 'C3H8'): 0,
                            ('C2H6', 'C4H10'): 0, ('C3H8', 'H2'): 0,
                            ('C3H8', 'CO'): 0, ('C3H8', 'H2O'): 0,
                            ('C3H8', 'CO2'): 0, ('C3H8', 'O2'): 0,
                            ('C3H8', 'N2'): 0, ('C3H8', 'Ar'): 0,
                            ('C3H8', 'CH4'): 0, ('C3H8', 'C2H6'): 0,
                            ('C3H8', 'C3H8'): 0, ('C3H8', 'C4H10'): 0,
                            ('C4H10', 'H2'): 0, ('C4H10', 'CO'): 0,
                            ('C4H10', 'H2O'): 0, ('C4H10', 'CO2'): 0,
                            ('C4H10', 'O2'): 0, ('C4H10', 'N2'): 0,
                            ('C4H10', 'Ar'): 0, ('C4H10', 'CH4'): 0,
                            ('C4H10', 'C2H6'): 0,
                            ('C4H10', 'C3H8'): 0,
                            ('C4H10', 'C4H10'): 0}}}

    # logic to select a subset of components
    if components == 'all':
        return configuration
    else:
        # collect required components
        comp_dict = {}
        for key in components:
            comp_dict[key] = configuration['components'][key]
        # create a copy of configuration and overwrite components
        config_copy = copy.copy(configuration)
        config_copy['components'] = comp_dict
        return config_copy


# stoichiometry for natural gas combustion
rxn_configuration = {
    "rate_reactions": {
        "R1": {"stoichiometry": {("Vap", "H2"): -1,
                                 ("Vap", "O2"): -0.5,
                                 ("Vap", "H2O"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": 0,
                   "arrhenius_const": 0,
                   "energy_activation": 0}},
        "R2": {"stoichiometry": {("Vap", "CO"): -1,
                                 ("Vap", "O2"): -0.5,
                                 ("Vap", "CO2"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": 0,
                   "arrhenius_const": 0,
                   "energy_activation": 0}},
        "R3": {"stoichiometry": {("Vap", "CH4"): -1,
                                 ("Vap", "O2"): -2,
                                 ("Vap", "H2O"): 2,
                                 ("Vap", "CO2"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": 0,
                   "arrhenius_const": 0,
                   "energy_activation": 0}}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K}}
