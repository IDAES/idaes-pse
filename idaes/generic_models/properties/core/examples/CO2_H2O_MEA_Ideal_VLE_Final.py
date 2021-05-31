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
Property package for MEA Carbon Capture System - Stripper Reboiler

CO2-H2O-MEA System
Thermodynamic Model: Ideal liquid and vapor phase VLE with non-volatiles
"""
import sys
import os
import logging

# Import units package from Pyomo
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.core.phases import PhaseType as PT
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import \
        IdealBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import fugacity
import idaes.generic_models.properties.core.pure.Perrys as Perrys
import idaes.generic_models.properties.core.pure.NIST as NIST

# Import custom methods required to build the property package
from idaes.generic_models.properties.core.examples import CustomMethods_Final,\
    CustomMethods2_Final

# Set up logger
_log = logging.getLogger(__name__)

# Configuration dictionary for an ideal CO2-H2O-MEA system
# Assumptions
# 1. Ideal gas-liquid system, based on apparent components
# 2. MEA is non-volatile
# 3. Properties are applicable at standard pressure (1 atm)

# Reference Temperature : 298.15 K
# Reference Pressure : 101325 Pa   

# Data Sources:

# [1] NIST Webbook, https://webbook.nist.gov/
      # Converted from bar to Pa
      
# [2] Perry's Chemical Engineers' Handbook 7th Ed.

# [3] Model Development, Validation, and Optimization of an 
      # MEA-Based Post-Combustion CO2 Capture Process under Part-Load and 
      # Variable Capture Operations: Paul Akula, John Eslick, Debangsu Bhattacharyya, 
      # and David C. Miller, I&EC research article - Supporting Information
      
# [4] A Predictive Thermodynamic Model for an Aqueous Blend of Potassium 
      # Carbonate, Piperazine, and Monoethanolamine for Carbon Dioxide 
      # Capture from Flue Gas: Marcus Douglas Hilliard
      
# [5] Uncertainty Quantification of Property Models: 
      # Methodology and Its Application to CO2-Loaded Aqueous MEA Solutions 
      # Joshua C. Morgan, Debangsu Bhattacharya, Charles Tong, David Miller
      
# [6] A Semi-Empirical Model for Estimating the Heat Capacity of Aqueous
      # Solutions of Alkanolamines for CO2 Capture:
      # Elvis O. Agbonghae,* Kevin J. Hughes, Derek B. Ingham, Lin Ma, 
      # and Mohamed Pourkashanian
      # Energy Technology and Innovation Initiative (ETII), University of Leeds
      # Leeds, LS2 9JT, United Kingdom - I&EC research article 


configuration = {
    # Specifying components
    "components": {
        'H2O': {"type": Component,
                "dens_mol_liq_comp": CustomMethods_Final,
                "enth_mol_liq_comp": Perrys,
                "enth_mol_ig_comp": NIST,
                "pressure_sat_comp": NIST,
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                "parameter_data": {
                    "mw": (18.0153E-3, pyunits.kg/pyunits.mol),  # [1]
                    "pressure_crit": (220.64E5, pyunits.Pa),  # [1]
                    "temperature_crit": (647, pyunits.K),  # [1]
                    "mol_vol_liq_comp_coeff": { # [5]
                        '1': (-3.2484e-6, pyunits.g*pyunits.cm**-3*pyunits.K**-2),  
                        '2': (0.00165, pyunits.g*pyunits.cm**-3*pyunits.K**-1),
                        '3': (0.793,pyunits.g*pyunits.cm**-3)},
                    "cp_mol_ig_comp_coeff": {
                        'A': (30.09200, pyunits.J/pyunits.mol/pyunits.K), # [1] temperature range 500 K- 1700 K
                        'B': (6.832514, pyunits.J*pyunits.mol**-1*\
                              pyunits.K**-1*pyunits.kiloK**-1),
                        'C': (6.793435, pyunits.J*pyunits.mol**-1*\
                              pyunits.K**-1*pyunits.kiloK**-2),
                        'D': (-2.534480, pyunits.J*pyunits.mol**-1*\
                              pyunits.K**-1*pyunits.kiloK**-3),
                        'E': (0.082139, pyunits.J*pyunits.mol**-1*\
                              pyunits.K**-1*pyunits.kiloK**2),
                        'F': (-250.8810, pyunits.kJ/pyunits.mol),
                        'G': (223.3967, pyunits.J/pyunits.mol/pyunits.K),
                        'H': (0, pyunits.kJ/pyunits.mol)},
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K), # [2] pg 2-174, temperature range 273.16 K - 533.15 K
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "enth_mol_form_liq_comp_ref": (
                        -285.83E3, pyunits.J/pyunits.mol),  # [1]
                    "pressure_sat_comp_coeff": {
                        'A': (3.55959, None),  # [1], temperature range 379 K - 573 K
                        'B': (643.748, pyunits.K),
                        'C': (-198.043, pyunits.K)}}},
        'MEA': {"type": Component,
                "valid_phase_types": PT.liquidPhase,
                "dens_mol_liq_comp": CustomMethods_Final,
                "enth_mol_liq_comp": CustomMethods_Final,
                "parameter_data": {
                    "mw": (61.0831E-3, pyunits.kg/pyunits.mol), # [1] 
                    "pressure_crit": (80.3E5, pyunits.Pa), # [1]
                    "temperature_crit": (671.4, pyunits.K), # [1]
                    "mol_vol_liq_comp_coeff": { # [5]
                        '1': (-5.35162E-7, pyunits.g*pyunits.cm**-3*\
                              pyunits.K**-2),  
                        '2': (-4.51417E-4, pyunits.g*pyunits.cm**-3*\
                              pyunits.K**-1),
                        '3': (1.19451,pyunits.g*pyunits.cm**-3)},
                    "cp_mol_liq_comp_coeff": { # [6]
                        'a0': (78.2498, pyunits.J*pyunits.mol**-1*\
                              pyunits.K**-1),  
                        'a1': (0.293, pyunits.J*pyunits.mol**-1*\
                              pyunits.K**-2)},
                    "enth_mol_form_liq_comp_ref": ( # [1]
                        -507.5E3, pyunits.J/pyunits.mol)}},
        'CO2': {"type": Component,
                "enth_mol_liq_comp": CustomMethods2_Final, # [4]
                "enth_mol_ig_comp": NIST, 
                "pressure_sat_comp": CustomMethods2_Final,
                "phase_equilibrium_form": {("Vap", "Liq"): \
                                           CustomMethods2_Final.fugacity},
                "parameter_data": {
                   "mw": (44.0095E-3, pyunits.kg/pyunits.mol),  # [1]
                   "pressure_crit": (73.825E5, pyunits.Pa),  # [1]
                   "temperature_crit": (304.23, pyunits.K),  # [1]
                   "cp_mol_ig_comp_coeff": {                 # [1], temperature range 298 K - 1200 K
                       'A': (24.99735, pyunits.J/pyunits.mol/pyunits.K),
                       'B': (55.18696, pyunits.J*pyunits.mol**-1*pyunits.K**-1*\
                             pyunits.kiloK**-1),
                       'C': (-33.69137, pyunits.J*pyunits.mol**-1*pyunits.K**-1*\
                             pyunits.kiloK**-2),
                       'D': (7.948387, pyunits.J*pyunits.mol**-1*pyunits.K**-1*\
                             pyunits.kiloK**-3),
                       'E': (-0.136638, pyunits.J*pyunits.mol**-1*pyunits.K**-1*\
                             pyunits.kiloK**2),
                       'F': (-403.6075, pyunits.kJ/pyunits.mol),
                       'G': (228.2431, pyunits.J/pyunits.mol/pyunits.K),
                       'H': (0, pyunits.kJ/pyunits.mol)}
                         }}},

    # Specifying phases
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Ideal},
                'Vap': {"type": VaporPhase,
                        "equation_of_state": Ideal}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 92.0817, 100, pyunits.mol/pyunits.s),
                     "temperature": (350, 380, 500, pyunits.K),
                     "pressure": (5E4, 183430, 1e7, pyunits.Pa)
                     },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),

    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
    "bubble_dew_method": IdealBubbleDew}
