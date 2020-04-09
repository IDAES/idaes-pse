##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
from idaes.core import declare_process_block_class, LiquidPhase, VaporPhase

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterData)

from idaes.generic_models.properties.core.state_definitions import FTPx
import idaes.generic_models.properties.core.eos.ideal as ideal
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import (
        bubble_temp_ideal,
        dew_temp_ideal,
        bubble_press_ideal,
        dew_press_ideal)
from idaes.generic_models.properties.core.phase_equil.forms import fugacity

import idaes.generic_models.properties.core.pure.Perrys as Perrys
import idaes.generic_models.properties.core.pure.RPP as RPP

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("BTIdealParameterBlock")
class BTIdealParameterData(GenericParameterData):
    def configure(self):
        '''
        Method to set construction arguments.
        '''
        # ---------------------------------------------------------------------
        # Set config arguments

        # Data Sources:
        # [1] The Properties of Gases and Liquids (1987)
        #     4th edition, Chemical Engineering Series - Robert C. Reid
        # [2] Perry's Chemical Engineers' Handbook 7th Ed.
        #     Converted to J/mol.K, mol/m^3
        # [3] Engineering Toolbox, https://www.engineeringtoolbox.com
        #     Retrieved 1st December, 2019
        self.config.components = {
            'benzene': {"dens_mol_liq_comp": Perrys,
                        "enth_mol_liq_comp": Perrys,
                        "enth_mol_ig_comp": RPP,
                        "pressure_sat_comp": RPP,
                        "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                        "parameter_data": {
                            "mw": 78.1136E-3,  # [1]
                            "pressure_crit": 48.9e5,  # [1]
                            "temperature_crit": 562.2,  # [1]
                            "dens_mol_liq_comp_coeff": {'1': 1.0162*1e3,  # [2] pg. 2-98
                                                        '2': 0.2655,
                                                        '3': 562.16,
                                                        '4': 0.28212},
                            "cp_mol_ig_comp_coeff": {'A': -3.392E1,  # [1]
                                                     'B': 4.739E-1,
                                                     'C': -3.017E-4,
                                                     'D': 7.130E-8},
                            "cp_mol_liq_comp_coeff": {'1': 1.29E2,  # [2]
                                                      '2': -1.7E-1,
                                                      '3': 6.48E-4,
                                                      '4': 0,
                                                      '5': 0},
                            "enth_mol_form_liq_comp_ref": 49.0e3,  # [3]
                            "enth_mol_form_vap_comp_ref": 82.9e3,  # [3]
                            "pressure_sat_comp_coeff": {'A': -6.98273,  # [1]
                                                        'B': 1.33213,
                                                        'C': -2.62863,
                                                        'D': -3.33399}}},
            'toluene': {"dens_mol_liq_comp": Perrys,
                        "enth_mol_liq_comp": Perrys,
                        "enth_mol_ig_comp": RPP,
                        "pressure_sat_comp": RPP,
                        "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                        "parameter_data": {
                            "mw": 92.1405E-3,  # [1]
                            "pressure_crit": 41e5,  # [1]
                            "temperature_crit": 591.8,  # [1]
                            "dens_mol_liq_comp_coeff": {'1': 0.8488*1e3,  # [2] pg. 2-98
                                                        '2': 0.26655,
                                                        '3': 591.8,
                                                        '4': 0.2878},
                            "cp_mol_ig_comp_coeff": {'A': -2.435E1,
                                                     'B': 5.125E-1,
                                                     'C': -2.765E-4,
                                                     'D': 4.911E-8},
                            "cp_mol_liq_comp_coeff": {'1': 1.40E2,  # [2]
                                                      '2': -1.52E-1,
                                                      '3': 6.95E-4,
                                                      '4': 0,
                                                      '5': 0},
                            "enth_mol_form_liq_comp_ref": 12.0e3,  # [3]
                            "enth_mol_form_vap_comp_ref": 50.1e3,  # [3]
                            "pressure_sat_comp_coeff": {'A': -7.28607,  # [1]
                                                        'B': 1.38091,
                                                        'C': -2.83433,
                                                        'D': -2.79168}}}}
        self.config.phases = {
            'Liq': {"type": LiquidPhase,
                    "equation_of_state": ideal},
            'Vap': {"type": VaporPhase,
                    "equation_of_state": ideal}}

        self.config.state_definition = FTPx
        self.config.state_bounds = {"flow_mol": (0, 1000),
                                    "temperature": (273.15, 450),
                                    "pressure": (5e4, 1e6)}
        self.config.pressure_ref = 1e5
        self.config.temperature_ref = 300

        self.config.phases_in_equilibrium = [("Vap", "Liq")]
        self.config.phase_equilibrium_state = {
            ("Vap", "Liq"): smooth_VLE}

        self.config.temperature_bubble = bubble_temp_ideal
        self.config.temperature_dew = dew_temp_ideal
        self.config.pressure_bubble = bubble_press_ideal
        self.config.pressure_dew = dew_press_ideal
