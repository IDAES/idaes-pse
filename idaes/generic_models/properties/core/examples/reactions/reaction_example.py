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
Example of defining a reaction package for the follwoing reaction system:

R1)  A + B --> 2C  Power law rate reaction
R2)  B + C <-> D  Power law equilibrium reaction

Author: Andrew Lee
"""
from idaes.core import LiquidPhase, Component

from idaes.generic_models.properties.core.state_definitions import FcTP
from idaes.generic_models.properties.core.eos.ideal import Ideal

from idaes.generic_models.properties.core.reactions.dh_rxn import \
    constant_dh_rxn
from idaes.generic_models.properties.core.reactions.rate_constant import \
    arrhenius
from idaes.generic_models.properties.core.reactions.rate_forms import \
    mole_frac_power_law_rate
from idaes.generic_models.properties.core.reactions.equilibrium_constant import \
    van_t_hoff
from idaes.generic_models.properties.core.reactions.equilibrium_forms import \
    mole_frac_power_law_equil


# First, create a thermophsyical property definition that will be used
# with the reactions

# For this example, the thermophsycial properties will only define components
# and phases, but in practice users would also need to define some properties
thermo_configuration = {
    "components": {
        'A': {"type": Component},
        'B': {"type": Component},
        'C': {"type": Component},
        'D': {"type": Component}},
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Ideal}},
    "state_definition": FcTP,
    "state_bounds": {"flow_mol": (0, 1000),
                     "temperature": (273.15, 450),
                     "pressure": (5e4, 1e6)},
    "pressure_ref": 1e5,
    "temperature_ref": 300}


# Next, create the reaction property definition which describes the system on
# reactions to be modeled.
rxn_configuration = {
    "rate_reactions": {
        "R1": {"stoichiometry": {("Liq", "A"): -1,
                                 ("Liq", "B"): -1,
                                 ("Liq", "C"): 2},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": mole_frac_power_law_rate,
               "parameter_data": {
                   "dh_rxn_ref": -10000,
                   "arrhenius_const": 1,
                   "energy_activation": 1000}}},
    "equilibrium_reactions": {
        "R2": {"stoichiometry": {("Liq", "B"): -1,
                                 ("Liq", "C"): -1,
                                 ("Liq", "D"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "equilibrium_constant": van_t_hoff,
               "equilibrium_form": mole_frac_power_law_equil,
               "parameter_data": {
                   "dh_rxn_ref": -20000,
                   "k_eq_ref": 100,
                   "T_eq_ref": 350}}}}
