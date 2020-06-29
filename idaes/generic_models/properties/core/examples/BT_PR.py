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
Benzene-Toluene phase equilibrium package using ideal liquid and vapor.

Example property package using the Generic Property Package Framework.
This exmample shows how to set up a property package to do benzene-toluene
phase equilibrium in the generic framework using ideal liquid and vapor
assumptions along with methods drawn from the pre-built IDAES property
libraries.
"""
# Import Python libraries
import logging

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import \
        LogBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import log_fugacity

import idaes.generic_models.properties.core.pure.RPP as RPP

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019

configuration = {
    # Specifying components
    "components": {
        'benzene': {"type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "pressure_sat_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": 78.1136E-3,  # [1]
                        "pressure_crit": 48.9e5,  # [1]
                        "temperature_crit": 562.2,  # [1]
                        "omega": 0.212,  # [1]
                        "cp_mol_ig_comp_coeff": {'A': -3.392E1,  # [1]
                                                 'B': 4.739E-1,
                                                 'C': -3.017E-4,
                                                 'D': 7.130E-8},
                        "enth_mol_form_vap_comp_ref": 82.9e3,  # [3]
                        "entr_mol_form_vap_comp_ref": -269,  # [3]
                        "pressure_sat_comp_coeff": {'A': -6.98273,  # [1]
                                                    'B': 1.33213,
                                                    'C': -2.62863,
                                                    'D': -3.33399}}},
        'toluene': {"type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "pressure_sat_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": 92.1405E-3,  # [1]
                        "pressure_crit": 41e5,  # [1]
                        "temperature_crit": 591.8,  # [1]
                        "omega": 0.263,  # [1]
                        "cp_mol_ig_comp_coeff": {'A': -2.435E1,
                                                 'B': 5.125E-1,
                                                 'C': -2.765E-4,
                                                 'D': 4.911E-8},
                        "enth_mol_form_vap_comp_ref": 50.1e3,  # [3]
                        "entr_mol_form_vap_comp_ref": -321,  # [3]
                        "pressure_sat_comp_coeff": {'A': -7.28607,  # [1]
                                                    'B': 1.38091,
                                                    'C': -2.83433,
                                                    'D': -2.79168}}}},

    # Specifying phases
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.PR}},
                'Vap': {"type": VaporPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.PR}}},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 1000),
                     "temperature": (273.15, 500),
                     "pressure": (5e4, 1e6)},
    "pressure_ref": 101325,
    "temperature_ref": 298.15,

    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {"PR_kappa": {("benzene", "benzene"): 0.000,
                                    ("benzene", "toluene"): 0.000,
                                    ("toluene", "benzene"): 0.000,
                                    ("toluene", "toluene"): 0.000}}}
