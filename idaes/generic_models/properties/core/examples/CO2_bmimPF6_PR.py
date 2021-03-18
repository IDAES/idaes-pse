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
CO2-bmimPF6 phase equilibrium package using ideal liquid and vapor.

Example property package using the Generic Property Package Framework.
This exmample shows how to set up a property package to do CO2-Ionic Liquid
phase equilibrium in the generic framework using Peng-robinsons EOS
assumptions along with methods drawn from the pre-built IDAES property
libraries.
"""
# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.core.phases import PhaseType as PT
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import \
        LogBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import log_fugacity


import idaes.generic_models.properties.core.pure.RPP4 as RPP4
import idaes.generic_models.properties.core.pure.NIST as NIST


# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [2] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019
# [3] Separation of CO2 and H2S using room-temperature ionic liquid [bmim][PF6]
#     Mark B.Shiflett, A.Yokozeki, 2010
# [4] Thermodynamic Properties of Imidazolium-Based Ionic Liquids:  Densities, Heat Capacities, and Enthalpies of Fusion of [bmim][PF6] and [bmim][NTf2]
#     Jacobo Troncoso, Claudio A. Cerdeiriña, Yolanda A. Sanmamed, Luís Romaní, and Luís Paulo N. Rebelo, 2006

configuration = {
    # Specifying components
    "components": {
        "bmimPF6": {"type": Component,
                    "enth_mol_ig_comp": RPP4,
                    "entr_mol_ig_comp": RPP4,
                    "valid_phase_types": PT.liquidPhase,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (284.18E-3, pyunits.kg/pyunits.mol),  # [3]
                        "pressure_crit": (24e5, pyunits.Pa),  # [3]
                        "temperature_crit": (860, pyunits.K),  # [3]
                        "omega": 0.7917,  # [4]
                        "cp_mol_ig_comp_coeff": {
                            'A': (3259.5745, pyunits.J/pyunits.mol/pyunits.K),  # [4]
                            'B': (-28.5610, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (0.09354, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (-0.0001000673, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "enth_mol_form_vap_comp_ref": (
                            0.145, pyunits.J/pyunits.mol),  # [4]
                        "entr_mol_form_vap_comp_ref": (
                            137.5, pyunits.J/pyunits.mol/pyunits.K)}},

        "carbon_dioxide": {"type": Component,
                  "enth_mol_ig_comp": RPP4,
                  "entr_mol_ig_comp": RPP4,
                  "pressure_sat_comp": NIST,
                  "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                  "parameter_data": {
                      "mw": (44.010E-3, pyunits.kg/pyunits.mol),  # [1]
                      "pressure_crit": (71.8e5, pyunits.Pa),  # [1]
                      "temperature_crit": (304.1, pyunits.K),  # [1]
                      "omega": 0.239,  # [1]
                      "cp_mol_ig_comp_coeff": {
                          "A": (1.980E1, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                          "B": (7.344E-2, pyunits.J/pyunits.mol/pyunits.K**2),
                          "C": (-5.602E-5, pyunits.J/pyunits.mol/pyunits.K**3),
                          "D": (1.715E-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                      "enth_mol_form_vap_comp_ref": (
                          -393.51, pyunits.J/pyunits.mol),  # [2]
                      "entr_mol_form_vap_comp_ref": (
                          -269, pyunits.J/pyunits.mol/pyunits.K),  # [2]
                      "pressure_sat_comp_coeff": {"A": (6.81228, None),  # [2]
                                                  "B": (1301.679, pyunits.K),
                                                  "C": (-3.494, pyunits.K)}}}},

    # Specifying phases
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.PR}},
                'Vap': {"type": VaporPhase,
                        "equation_of_state": Cubic,
                        "equation_of_state_options": {
                            "type": CubicType.PR}}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                     "temperature": (10, 300, 500, pyunits.K),
                     "pressure": (5e-4, 1e5, 1e10, pyunits.Pa)},
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),

    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {"PR_kappa": {("bmimPF6", "bmimPF6"): 0.000,
                                    ("bmimPF6", "carbon_dioxide"): -0.4071,
                                    ("carbon_dioxide", "carbon_dioxide"): 0.000,
                                    ("carbon_dioxide", "bmimPF6"): 0.0206}}}
