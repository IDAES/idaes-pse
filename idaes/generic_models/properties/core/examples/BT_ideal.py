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
from pyomo.environ import Var

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

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        self.config.components = {
            'benzene': {"dens_mol_liq_comp": Perrys,
                        "enth_mol_liq_comp": Perrys,
                        "enth_mol_ig_comp": RPP,
                        "entr_mol_liq_comp": Perrys,
                        "entr_mol_ig_comp": RPP,
                        "pressure_sat_comp": RPP,
                        "phase_equilibrium_form": fugacity,
                        "mw": 78.1136E-3,
                        "pressure_crit": 48.9e5,
                        "temperature_crit": 562.2},
            'toluene': {"dens_mol_liq_comp": Perrys,
                        "enth_mol_liq_comp": Perrys,
                        "enth_mol_ig_comp": RPP,
                        "entr_mol_liq_comp": Perrys,
                        "entr_mol_ig_comp": RPP,
                        "pressure_sat_comp": RPP,
                        "phase_equilibrium_form": fugacity,
                        "mw": 92.1405E-3,
                        "pressure_crit": 41e5,
                        "temperature_crit": 591.8}}
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

        self.config.temperature_bubble = bubble_temp_ideal
        self.config.temperature_dew = dew_temp_ideal
        self.config.pressure_bubble = bubble_press_ideal
        self.config.pressure_dew = dew_press_ideal

    def parameters(self):
        # Constants for ideal gas specific enthalpy
        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        self.benzene.cp_mol_ig_comp_coeff = Var(
                ['A', 'B', 'C', 'D'],
                initialize={'A': -3.392E1,
                            'B': 4.739E-1,
                            'C': -3.017E-4,
                            'D': 7.130E-8},
                doc="Parameters for ideal gas heat capacity [J/mol.K]")

        self.toluene.cp_mol_ig_comp_coeff = Var(
                ['A', 'B', 'C', 'D'],
                initialize={'A': -2.435E1,
                            'B': 5.125E-1,
                            'C': -2.765E-4,
                            'D': 4.911E-8},
                doc="Parameters for ideal gas heat capacity [J/mol.K]")

        # Constants for liquid phase specific enthalpy
        # Source: Perry's Chemical Engineers' Handbook 7th Ed.
        # Units converted to J/mol.K
        self.benzene.cp_mol_liq_comp_coeff = Var(
                ['1', '2', '3', '4', '5'],
                initialize={'1': 1.29E2,
                            '2': -1.7E-1,
                            '3': 6.48E-4,
                            '4': 0,
                            '5': 0},
                doc="Parameters for liquid cp [J/mol.K]")

        self.toluene.cp_mol_liq_comp_coeff = Var(
                ['1', '2', '3', '4', '5'],
                initialize={'1': 1.40E2,
                            '2': -1.52E-1,
                            '3': 6.95E-4,
                            '4': 0,
                            '5': 0},
                doc="Parameters for liquid cp [J/mol.K]")

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        self.benzene.pressure_sat_comp_coeff = Var(
            ['A', 'B', 'C', 'D'],
            initialize={'A': -6.98273,
                        'B': 1.33213,
                        'C': -2.62863,
                        'D': -3.33399},
            doc="Parameters for saturation pressure [Pa]")

        self.toluene.pressure_sat_comp_coeff = Var(
            ['A', 'B', 'C', 'D'],
            initialize={'A': -7.28607,
                        'B': 1.38091,
                        'C': -2.83433,
                        'D': -2.79168},
            doc="Parameters for saturation pressure [Pa]")

        # Source: "Perry's Chemical Engineers' Handbook by Robert H. Perry"
        # 7th Edition, pg. 2-98
        # Units converted to mol/m^3
        self.benzene.dens_mol_liq_comp_coeff = Var(
            ['1', '2', '3', '4'],
            initialize={'1': 1.0162*1e3,
                        '2': 0.2655,
                        '3': 562.16,
                        '4': 0.28212},
            doc="Coefficients for liquid molar densities [mol/m^3]")

        self.toluene.dens_mol_liq_comp_coeff = Var(
            ['1', '2', '3', '4'],
            initialize={'1': 0.8488*1e3,
                        '2': 0.26655,
                        '3': 591.8,
                        '4': 0.2878},
            doc="Coefficients for liquid molar densities [mol/m^3]")

        # Standard heats of formation
        # Source: NIST Webbook, https://webbook.nist.gov
        # Retrieved 1st December 2019
        self.benzene.enth_mol_form_phase_comp_ref = Var(
                self.phase_list,
                initialize={'Liq': 49.0e3, 'Vap': 82.9e3},
                doc="Molar heat of formation @ Tref [J/mol]")

        self.toluene.enth_mol_form_phase_comp_ref = Var(
                self.phase_list,
                initialize={'Liq': 12.0e3, 'Vap': 50.1e3},
                doc="Molar heat of formation @ Tref [J/mol]")

        # Standard entropy of formation
        # Source: Engineering Toolbox, https://www.engineeringtoolbox.com
        # Retrieved 1st December, 2019
        self.benzene.entr_mol_phase_comp_ref = Var(
                self.phase_list,
                initialize={'Liq': -173, 'Vap': -269},
                doc="Standard molar entropy @ Tref [J/mol.K]")

        self.toluene.entr_mol_phase_comp_ref = Var(
                self.phase_list,
                initialize={'Liq': -220, 'Vap': -321},
                doc="Standard molar entropy @ Tref [J/mol.K]")
