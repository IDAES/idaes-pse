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

"""
# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals, Set

# Import IDAES cores
from idaes.core import declare_process_block_class

from idaes.property_models.core.generic.generic_property import (
        GenericParameterData)

from idaes.property_models.core.state_definitions import FTPx
import idaes.property_models.core.eos.ideal as ideal
from idaes.property_models.core.phase_equil import smooth_VLE
from idaes.property_models.core.generic.bubble_dew import (bubble_temp_ideal,
                                                           dew_temp_ideal,
                                                           bubble_press_ideal,
                                                           dew_press_ideal)

import idaes.property_models.core.pure.Perrys as Perrys
import idaes.property_models.core.pure.RPP as RPP

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("BTIdealParameterBlock")
class BTIdealParameterData(GenericParameterData):

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(BTIdealParameterData, self).build()

        # ---------------------------------------------------------------------
        # Set config arguments
        self.config.state_definition = FTPx

        self.config.state_bounds = {"flow_mol": (0, 1000),
                                    "temperature": (273.15, 450),
                                    "pressure": (5e4, 1e6)}

        self.config.equation_of_state = {"Vap": ideal,
                                         "Liq": ideal}

        self.config.phase_equilibrium_formulation = smooth_VLE

        self.config.bubble_temperature = bubble_temp_ideal
        self.config.dew_temperature = dew_temp_ideal
        self.config.bubble_pressure = bubble_press_ideal
        self.config.dew_pressure = dew_press_ideal

        self.config.dens_mol_comp_liq = Perrys
        self.config.enth_mol_comp_liq = Perrys
        self.config.enth_mol_comp_ig = RPP
        self.config.entr_mol_comp_liq = Perrys
        self.config.entr_mol_comp_ig = RPP
        self.config.pressure_sat_comp = RPP
        # ---------------------------------------------------------------------
        self.component_list = Set(initialize=['benzene', 'toluene'])

        self.phase_list = Set(initialize=["Vap", "Liq"], ordered=True)

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        # List of phase equilibrium index
        self.phase_equilibrium_idx = Set(initialize=[1, 2])

        self.phase_equilibrium_list = \
            {1: ["benzene", ("Vap", "Liq")],
             2: ["toluene", ("Vap", "Liq")]}

        # Thermodynamic reference state
        self.pressure_ref = Param(mutable=True,
                                  default=101325,
                                  doc='Reference pressure [Pa]')
        self.temperature_ref = Param(mutable=True,
                                     default=298.15,
                                     doc='Reference temperature [K]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_crit_data = {'benzene': 48.9e5,
                              'toluene': 41e5}

        self.pressure_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=pressure_crit_data,
            doc='Critical pressure [Pa]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        temperature_crit_data = {'benzene': 562.2,
                                 'toluene': 591.8}

        self.temperature_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=temperature_crit_data,
            doc='Critical temperature [K]')

        # Gas Constant
        self.gas_const = Param(within=NonNegativeReals,
                               mutable=False,
                               default=8.314,
                               doc='Gas constant [J/mol.K]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        mw_comp_data = {'benzene': 78.1136E-3,
                        'toluene': 92.1405E-3}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=mw_comp_data,
                             doc="Molecular weight [kg/mol]")

        # Constants for ideal gas specific enthalpy
        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        cp_ig_coeff_data = {('benzene', 'A'): -3.392E1,
                            ('benzene', 'B'): 4.739E-1,
                            ('benzene', 'C'): -3.017E-4,
                            ('benzene', 'D'): 7.130E-8,
                            ('toluene', 'A'): -2.435E1,
                            ('toluene', 'B'): 5.125E-1,
                            ('toluene', 'C'): -2.765E-4,
                            ('toluene', 'D'): 4.911E-8}

        self.cp_ig_coeff = Param(self.component_list,
                                 ['A', 'B', 'C', 'D'],
                                 mutable=False,
                                 initialize=cp_ig_coeff_data,
                                 doc="Parameters for ideal gas heat capacity")

        # Constants for liquid phase specific enthalpy
        # Source: Perry's Chemical Engineers Handbook 7th Ed.
        # Units converted to J/mol.K
        cp_liq_coeff_data = {('benzene', '1'): 1.29E2,
                             ('benzene', '2'): -1.7E-1,
                             ('benzene', '3'): 6.48E-4,
                             ('benzene', '4'): 0,
                             ('benzene', '5'): 0,
                             ('toluene', '1'): 1.40E2,
                             ('toluene', '2'): -1.52E-1,
                             ('toluene', '3'): 6.95E-4,
                             ('toluene', '4'): 0,
                             ('toluene', '5'): 0}

        self.cp_liq_coeff = Param(self.component_list,
                                  ['1', '2', '3', '4', '5'],
                                  mutable=False,
                                  initialize=cp_liq_coeff_data,
                                  doc="Parameters for liquid cp")

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_sat_coeff_data = {('benzene', 'A'): -6.98273,
                                   ('benzene', 'B'): 1.33213,
                                   ('benzene', 'C'): -2.62863,
                                   ('benzene', 'D'): -3.33399,
                                   ('toluene', 'A'): -7.28607,
                                   ('toluene', 'B'): 1.38091,
                                   ('toluene', 'C'): -2.83433,
                                   ('toluene', 'D'): -2.79168}

        self.pressure_sat_coeff = Param(
            self.component_list,
            ['A', 'B', 'C', 'D'],
            mutable=False,
            initialize=pressure_sat_coeff_data,
            doc="parameters to compute Cp_comp")

        # Source: "Perry's Chemical Engineers Handbook by Robert H. Perry"
        # 7th Edition, pg. 2-98
        # Units converted to mol/m^3
        dens_mol_liq_coeff_data = {('benzene', '1'): 1.0162*1e3,
                                   ('benzene', '2'): 0.2655,
                                   ('benzene', '3'): 562.16,
                                   ('benzene', '4'): 0.28212,
                                   ('toluene', '1'): 0.8488*1e3,
                                   ('toluene', '2'): 0.26655,
                                   ('toluene', '3'): 591.8,
                                   ('toluene', '4'): 0.2878}
        self.dens_mol_liq_coeff = Param(
            self.component_list,
            ['1', '2', '3', '4'],
            mutable=False,
            initialize=dens_mol_liq_coeff_data,
            doc="Coefficients for calculating liquid molar densities")

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        dh_vap_data = {'benzene': 3.377e4,
                       'toluene': 3.8262e4}

        self.dh_vap_ref = Param(
                self.component_list,
                mutable=False,
                initialize=dh_vap_data,
                doc="Molar heat of vaporization @ Tref [J/mol]")

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        ds_vap_data = {'benzene': 3.377e4/298.15,
                       'toluene': 3.8262e4/298.15}

        self.ds_vap_ref = Param(
                self.component_list,
                mutable=False,
                initialize=ds_vap_data,
                doc="Molar entropy of vaporization @ Tref [J/mol.K]")
