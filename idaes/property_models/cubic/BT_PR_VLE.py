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
Example Peng-Robinson parameter block for the VLE calucations for an N2-O2-Ar
system.
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Reals, Param, NonNegativeReals, Set

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.core.util.misc import extract_data

from idaes.property_models.cubic.cubic_prop_pack_VLE import (
        CubicParameterData, CubicEoS)


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("BTParameterBlock")
class BTParameterData(CubicParameterData):

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(BTParameterData, self).build()

        self.cubic_type = CubicEoS.PR

        self.component_list = Set(initialize=['benzene', 'toluene'])

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        # List of phase equilibrium index
        self.phase_equilibrium_idx = Set(initialize=[1, 2, 3])

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
            initialize=extract_data(pressure_crit_data),
            doc='Critical pressure [Pa]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        temperature_crit_data = {'benzene': 562.2,
                                 'toluene': 591.8}

        self.temperature_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_crit_data),
            doc='Critical temperature [K]')

        # Pitzer acentricity factor (from Prop. Gases & Liquids)
        omega_data = {'benzene': 0.212,
                      'toluene': 0.263}

        self.omega = Param(
            self.component_list,
            within=Reals,
            mutable=False,
            initialize=extract_data(omega_data),
            doc='Acentricity Factor')

        # Peng-Robinson binary interaction parameters
        kappa_data = {
            ('benzene', 'benzene'): 0.0000, ('benzene', 'toluene'): 0.0000,
            ('toluene', 'benzene'): 0.0000, ('toluene', 'toluene'): 0.0000}

        self.kappa = Param(
            self.component_list,
            self.component_list,
            within=Reals,
            mutable=False,
            initialize=extract_data(kappa_data),
            doc='Peng-Robinson binary interaction parameters')

        # Gas Constant
        self.gas_const = Param(within=NonNegativeReals,
                               mutable=False,
                               default=8.314,
                               doc='Gas Constant [J/mol.K]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        mw_comp_data = {'benzene': 78.1136E-3,
                        'toluene': 92.1405E-3}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             doc="molecular weight Kg/mol")

        # Constants for specific heat capacity, enthalpy
        # Sources: The Properties of Gases and Liquids (1987)
        #         4th edition, Chemical Engineering Series - Robert C. Reid
        cp_ig_data = {('benzene', '1'): -3.392E1,
                      ('benzene', '2'): 4.739E-1,
                      ('benzene', '3'): -3.017E-4,
                      ('benzene', '4'): 7.130E-8,
                      ('benzene', '5'): 0,
                      ('toluene', '1'): -2.435E1,
                      ('toluene', '2'): 5.125E-1,
                      ('toluene', '3'): -2.765E-4,
                      ('toluene', '4'): 4.911E-8,
                      ('toluene', '5'): 0}

        self.cp_ig = Param(self.component_list,
                           ['1', '2', '3', '4', '5'],
                           mutable=False,
                           initialize=extract_data(cp_ig_data),
                           doc="Parameters to compute cp_comp")

        # Antoine coefficients for ideal vapour (units: bar, K)
        # This is needed for initial guesses of bubble and dew points
        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        antoine_data = {('benzene', '1'): 4.202,
                        ('benzene', '2'): 1322,
                        ('benzene', '3'): -38.56,
                        ('toluene', '1'): 4.216,
                        ('toluene', '2'): 1435,
                        ('toluene', '3'): -43.33}

        self.antoine = Param(self.component_list,
                             ['1', '2', '3'],
                             mutable=False,
                             initialize=extract_data(antoine_data),
                             doc="Antoine Parameters to pressure_sat")
