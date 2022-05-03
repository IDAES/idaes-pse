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
# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import VaporPhase, Component
from idaes.core.base.phases import PhaseType as PT
from idaes.models.properties.modular_properties.state_definitions import FPhx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil.forms import \
    fugacity

from idaes.models.properties.modular_properties.pure.RPP4 import RPP4
from idaes.models.properties.modular_properties.pure.NIST import NIST

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal syngas methanol system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [2] Perry's Chemical Engineers' Handbook 7th Ed.
# [3] NIST Chemistry WebBook, https://webbook.nist.gov/chemistry/
#     Retrieved 21st October, 2021

config_dict = {
    # Specifying components
    "components": {
        'CH4':
            {"type": Component,
             "elemental_composition": {"C": 1, "H": 4, "O": 0},
             "enth_mol_ig_comp": RPP4,
             "entr_mol_ig_comp": RPP4,
             "valid_phase_types": [PT.vaporPhase],
             "parameter_data": {
                 "mw": (16.043E-3, pyunits.kg/pyunits.mol),  # [1]
                 "pressure_crit": (46.0e5, pyunits.Pa),  # [1]
                 "temperature_crit": (190.4, pyunits.K),  # [1]
                 "cp_mol_ig_comp_coeff": {
                     'A': (1.925E1, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                     'B': (5.213E-2, pyunits.J/pyunits.mol/pyunits.K**2),
                     'C': (1.197E-5, pyunits.J/pyunits.mol/pyunits.K**3),
                     'D': (-1.132E-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                 "entr_mol_form_vap_comp_ref": (
                     186.25, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                 "enth_mol_form_vap_comp_ref": (
                     -74.6e3, pyunits.J/pyunits.mol)}},  # [3]
        'CO':
            {"type": Component,
             "elemental_composition": {"C": 1, "H": 0, "O": 1},
             "enth_mol_ig_comp": RPP4,
             "entr_mol_ig_comp": RPP4,
             "valid_phase_types": [PT.vaporPhase],
             "parameter_data": {
                 "mw": (28.010E-3, pyunits.kg/pyunits.mol),  # [1]
                 "pressure_crit": (35.0e5, pyunits.Pa),  # [1]
                 "temperature_crit": (132.9, pyunits.K),  # [1]
                 "cp_mol_ig_comp_coeff": {
                     'A': (3.087E1, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                     'B': (-1.285E-2, pyunits.J/pyunits.mol/pyunits.K**2),
                     'C': (2.789E-5, pyunits.J/pyunits.mol/pyunits.K**3),
                     'D': (-1.272E-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                 "entr_mol_form_vap_comp_ref": (
                     197.66, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                 "enth_mol_form_vap_comp_ref": (
                     -110.53e3, pyunits.J/pyunits.mol)}},  # [3]
        'H2':
            {"type": Component,
             "elemental_composition": {"C": 0, "H": 2, "O": 0},
             "enth_mol_ig_comp": RPP4,
             "entr_mol_ig_comp": RPP4,
             "valid_phase_types": [PT.vaporPhase],
             "parameter_data": {
                 "mw": (2.016E-3, pyunits.kg/pyunits.mol),  # [1]
                 "pressure_crit": (12.9e5, pyunits.Pa),  # [1]
                 "temperature_crit": (33.0, pyunits.K),  # [1]
                 "cp_mol_ig_comp_coeff": {
                     'A': (2.714E1, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                     'B': (9.274E-3, pyunits.J/pyunits.mol/pyunits.K**2),
                     'C': (-1.381E-5, pyunits.J/pyunits.mol/pyunits.K**3),
                     'D': (7.645E-9, pyunits.J/pyunits.mol/pyunits.K**4)},
                 "entr_mol_form_vap_comp_ref": (
                     130.68, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                 "enth_mol_form_vap_comp_ref": (
                     0.0, pyunits.J/pyunits.mol)}},  # [3]
        'CH3OH':
            {"type": Component,
             "elemental_composition": {"C": 1, "H": 4, "O": 1},
             "enth_mol_ig_comp": RPP4,
             "entr_mol_ig_comp": RPP4,
             "valid_phase_types": [PT.vaporPhase],
             "pressure_sat_comp": NIST,
             "phase_equilibrium_form": {("Vap"): fugacity},
             "parameter_data": {
                 "mw": (32.042E-3, pyunits.kg/pyunits.mol),  # [1]
                 "pressure_crit": (80.9e5, pyunits.Pa),  # [1]
                 "temperature_crit": (512.6, pyunits.K),  # [1]
                 "cp_mol_ig_comp_coeff": {
                     'A': (2.115E1, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                     'B': (7.092E-2, pyunits.J/pyunits.mol/pyunits.K**2),
                     'C': (2.587E-5, pyunits.J/pyunits.mol/pyunits.K**3),
                     'D': (-2.852E-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                 "entr_mol_form_vap_comp_ref": (
                     239.81, pyunits.J/pyunits.mol/pyunits.K),
                 "enth_mol_form_vap_comp_ref": (
                     -205e3, pyunits.J/pyunits.mol),  # [3]
                 "pressure_sat_comp_coeff": {'A': (5.15853, None),  # [3]
                                             'B': (1569.613, pyunits.K),
                                             'C': (-34.846, pyunits.K)}}}},

    # Specifying phases
    "phases":  {'Vap': {"type": VaporPhase,
                        "equation_of_state": Ideal}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FPhx,
    "state_bounds": {"flow_mol": (1e-10, 100, 1e10, pyunits.mol/pyunits.s),
                     "enth_mol": (-1e10, 100, 1e10, pyunits.J/pyunits.mol),
                     "temperature": (198.15, 298.15, 512.75, pyunits.K),
                     "pressure": (1e-10, 1e5, 1e10, pyunits.Pa)},
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K)}
