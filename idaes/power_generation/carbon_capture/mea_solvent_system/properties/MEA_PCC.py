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
Property Package for Monoethanolamine (MEA)-based Post Combustion CO2 Capture(PCC)
using the Generic Property Package Framework.This includes three configuration
dictionaries containing parameters required by property  methods  pre-built in
IDAES property libraries.
The configuration dictionaries are:
1. configuration_vapor_abs which  contains parameters for the vapor phase in the
absorber.
2. configuration_vapor_reg which  contains parameters for the vapor phase in the
stripper/regenerator,reboiler, and condenser.
3. configuration_liquid which  contains parameters for the liquid phase
applicable to the absorber,stripper,reboiler,plate heat exchanger and condenser

"""

# Import Python libraries
import logging
import copy

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.pure import RPP4
from idaes.generic_models.properties.core.pure import NIST

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Data Sources:
# [1] Gmehling, J., Kleiber, M., Kolbe, B., & Rarey, J. (2012).
#     Chemical thermodynamics for process simulation. 1st edition. Wiley-VCH
# [2]
# [3]
# [4]

absorber_vapor_phase_components = ['CO2','H2O','O2','N2']
stripper_vapor_phase_components = ['CO2','H2O']

configuration_vapor_abs = {
    # Specifying components
    "components": {
        "N2": {"type": Component,
                     "enth_mol_ig_comp": RPP4,
                     "parameter_data": {
                         "mw": (28.014E-3, pyunits.kg/pyunits.mol),  # [1]
                         "pressure_crit": (33.958E5, pyunits.Pa),  # [1]
                         "temperature_crit": (126.192, pyunits.K),  # [1]
                         "omega": 0.0372,  # [1]
                         "cp_mol_ig_comp_coeff": { # [1]
                             "A": (31.128,
                                   pyunits.J/pyunits.mol/pyunits.K),
                             "B": (-13.556E-3,
                                   pyunits.J/pyunits.mol/pyunits.K**2),
                             "C": (26.77E-6,
                                   pyunits.J/pyunits.mol/pyunits.K**3),
                             "D": (-11.673E-9,
                                   pyunits.J/pyunits.mol/pyunits.K**4)},
                         "enth_mol_form_vap_comp_ref": (
                             0.0, pyunits.J/pyunits.mol)}},

        "O2": {"type": Component,
                   "enth_mol_ig_comp": RPP4,
                   "parameter_data": {
                       "mw": (31.999E-3, pyunits.kg/pyunits.mol), # [1]
                       "pressure_crit": (50.464e5, pyunits.Pa),    # [1]
                       "temperature_crit": (154.599, pyunits.K),   # [1]
                       "omega": 0.0221,  # [1]
                       "cp_mol_ig_comp_coeff": { #[1]
                           "A": (28.087, pyunits.J/pyunits.mol/pyunits.K),
                           "B": (-0.004E-3,
                                 pyunits.J/pyunits.mol/pyunits.K**2),
                           "C": (17.447E-6, pyunits.J/pyunits.mol/pyunits.K**3),
                           "D": (-10.644E-9,
                                 pyunits.J/pyunits.mol/pyunits.K**4)},
                       "enth_mol_form_vap_comp_ref": (
                           0.0, pyunits.J/pyunits.mol)}},

        "H2O": {"type": Component,
                   "enth_mol_ig_comp": RPP4,
                   "parameter_data": {
                       "mw": (18.015E-3, pyunits.kg/pyunits.mol),  # [1]
                       "pressure_crit": (220.64e5, pyunits.Pa),    # [1]
                       "temperature_crit": (647.096, pyunits.K),   # [1]
                       "omega": 0.3443,  # [1]
                       "cp_mol_ig_comp_coeff": { # [1]
                           "A": (32.22, pyunits.J/pyunits.mol/pyunits.K),
                           "B": (1.923E3,
                                 pyunits.J/pyunits.mol/pyunits.K**2),
                           "C": (10.548E-6, pyunits.J/pyunits.mol/pyunits.K**3),
                           "D": (-3.594E-9,
                                 pyunits.J/pyunits.mol/pyunits.K**4)},
                       "enth_mol_form_vap_comp_ref": (
                           -241820.0, pyunits.J/pyunits.mol)}},

        "CO2": {"type": Component,
                   "enth_mol_ig_comp": RPP4,
                   "parameter_data": {
                       "mw": (44.009E-3, pyunits.kg/pyunits.mol),  # [1]
                       "pressure_crit": (73.773e5, pyunits.Pa),    # [1]
                       "temperature_crit": (304.128, pyunits.K),   # [1]
                       "omega": 0.2236,  # [1]
                       "cp_mol_ig_comp_coeff": { #[1]
                           "A": (32.22, pyunits.J/pyunits.mol/pyunits.K),
                           "B": (1.923E3,
                                 pyunits.J/pyunits.mol/pyunits.K**2),
                           "C": (10.548E-6, pyunits.J/pyunits.mol/pyunits.K**3),
                           "D": (-3.594E-9,
                                 pyunits.J/pyunits.mol/pyunits.K**4)},
                       "enth_mol_form_vap_comp_ref": (
                           -393500, pyunits.J/pyunits.mol)}}},

    # Specifying phases
    "phases": {"Vap": {"type": VaporPhase,
                       "equation_of_state": Ideal}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K)
    }


# configuration for stripper/regenerator
configuration_vapor_reg = copy.deepcopy(configuration_vapor_abs)

# remove components not in stripper
for comp in absorber_vapor_phase_components:
    if not comp in stripper_vapor_phase_components:
        configuration_vapor_reg['components'].pop(comp)


