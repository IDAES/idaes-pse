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
Phase equilibrium package for methanol synthesis using ideal VLE.
Author: Brandon Paul
"""
from pyomo.environ import units as pyunits

from idaes.models.properties.modular_properties.base.generic_reaction import (
        ConcentrationForm)
from idaes.models.properties.modular_properties.reactions.dh_rxn import \
    constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_forms import \
    power_law_rate
from idaes.models.properties.modular_properties.reactions.rate_constant import \
    arrhenius

# [1] Reaction properties and stoichiometric coefficients obtained from
# Nieminen, H.; Laari, A.; Koiranen, T. CO2 Hydrogenation to Methanol
# by a Liquid-Phase Process with Alcoholic Solvents: A Techno-Economic
# Analysis. Processes 2019, 7, 405. https://doi.org/10.3390/pr7070405

# Methane reformation to syngas: CH4 + O2 => CO2 + 2H2
# Syngas conversion to methanol: CO + 2H2 => CH3OH (assumed sole reaction here)


config_dict = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "rate_reactions": {
        "R1": {"stoichiometry": {("Vap", "CO"): -1,
                                 ("Vap", "H2"): -2,
                                 ("Vap", "CH3OH"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "concentration_form": ConcentrationForm.moleFraction,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "parameter_data": {
                    "reaction_order": {("Vap", "CO"): 1,
                                       ("Vap", "H2"): 2},
                    "arrhenius_const": (3.77287e19, pyunits.mol / pyunits.m**3
                                        / pyunits.s),
                    "dh_rxn_ref": (-90640, pyunits.J/pyunits.mol),
                    "energy_activation": (109.2e3, pyunits.J/pyunits.mol)}}}}
