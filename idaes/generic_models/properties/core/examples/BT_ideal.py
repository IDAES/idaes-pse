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

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals

# Import IDAES cores
from idaes.core import declare_process_block_class, LiquidPhase, VaporPhase

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterData)

from idaes.generic_models.properties.core.state_definitions import FTPx
import idaes.generic_models.properties.core.eos.ideal as ideal
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
# from idaes.generic_models.properties.core.phase_equil.bubble_dew import (
#         bubble_temp_ideal,
#         dew_temp_ideal,
#         bubble_press_ideal,
#         dew_press_ideal)

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
        self.config.components = {
            'benzene': {"enth_mol_liq_comp": Perrys,
                        "enth_mol_ig_comp": RPP},
            'toluene': {"enth_mol_liq_comp": Perrys,
                        "enth_mol_ig_comp": RPP}}
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

        self.config.phase_equilibrium_formulation = smooth_VLE
        self.config.phases_in_equilibrium = [("Vap", "Liq")]

        # self.config.temperature_bubble = bubble_temp_ideal
        # self.config.temperature_dew = dew_temp_ideal
        # self.config.pressure_bubble = bubble_press_ideal
        # self.config.pressure_dew = dew_press_ideal

        # self.config.dens_mol_liq_comp = Perrys
        # self.config.entr_mol_liq_comp = Perrys
        # self.config.entr_mol_ig_comp = RPP
        # self.config.pressure_sat_comp = RPP

        # # Constants for ideal gas specific enthalpy
        # # Source: The Properties of Gases and Liquids (1987)
        # # 4th edition, Chemical Engineering Series - Robert C. Reid
        # cp_mol_ig_comp_coeff_data = {('benzene', 'A'): -3.392E1,
        #                              ('benzene', 'B'): 4.739E-1,
        #                              ('benzene', 'C'): -3.017E-4,
        #                              ('benzene', 'D'): 7.130E-8,
        #                              ('toluene', 'A'): -2.435E1,
        #                              ('toluene', 'B'): 5.125E-1,
        #                              ('toluene', 'C'): -2.765E-4,
        #                              ('toluene', 'D'): 4.911E-8}

        # self.cp_mol_ig_comp_coeff = Param(
        #         self.component_list,
        #         ['A', 'B', 'C', 'D'],
        #         mutable=False,
        #         initialize=cp_mol_ig_comp_coeff_data,
        #         doc="Parameters for ideal gas heat capacity [J/mol.K]")

        # # Constants for liquid phase specific enthalpy
        # # Source: Perry's Chemical Engineers' Handbook 7th Ed.
        # # Units converted to J/mol.K
        # cp_mol_liq_comp_coeff_data = {('benzene', '1'): 1.29E2,
        #                               ('benzene', '2'): -1.7E-1,
        #                               ('benzene', '3'): 6.48E-4,
        #                               ('benzene', '4'): 0,
        #                               ('benzene', '5'): 0,
        #                               ('toluene', '1'): 1.40E2,
        #                               ('toluene', '2'): -1.52E-1,
        #                               ('toluene', '3'): 6.95E-4,
        #                               ('toluene', '4'): 0,
        #                               ('toluene', '5'): 0}

        # self.cp_mol_liq_comp_coeff = Param(
        #         self.component_list,
        #         ['1', '2', '3', '4', '5'],
        #         mutable=False,
        #         initialize=cp_mol_liq_comp_coeff_data,
        #         doc="Parameters for liquid cp [J/mol.K]")
