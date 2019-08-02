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
Example Peng-Robinson parameter block for the VLE calucations for H2O.
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


@declare_process_block_class("H2OParameterBlock")
class H2OParameterData(CubicParameterData):

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(H2OParameterData, self).build()

        self.cubic_type = CubicEoS.PR

        self.component_list = Set(initialize=['H2O'])

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        # List of phase equilibrium index
        self.phase_equilibrium_idx = Set(initialize=[1])

        self.phase_equilibrium_list = \
            {1: ["H2O", ("Vap", "Liq")]}

        # List of all chemical elements that constitute the chemical species
        self.elem = ['H', 'O']

        # Elemental composition of all species
        self.elem_comp = {'H2O': {'H': 2, 'O': 1}}

        # Thermodynamic reference state
        self.pressure_ref = Param(mutable=True,
                                  default=101325,
                                  doc='Reference pressure [Pa]')
        self.temperature_ref = Param(mutable=True,
                                     default=298.15,
                                     doc='Reference temperature [K]')

        pressure_crit_data = {'H2O': 22120000}

        self.pressure_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(pressure_crit_data),
            doc='Critical pressure [Pa]')

        temperature_crit_data = {'H2O': 647.3}

        self.temperature_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_crit_data),
            doc='Critical temperature [K]')

        # Pitzer acentricity factor (from Prop. Gases & Liquids)
        omega_data = {'H2O': 0.344}

        self.omega = Param(
            self.component_list,
            within=Reals,
            mutable=False,
            initialize=extract_data(omega_data),
            doc='Acentricity Factor')

        # Peng-Robinson binary interaction parameters
        kappa_data = {
            ('H2O', 'H2O'): 0.0000}

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
        mw_comp_data = {'H2O': 0.018015}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             doc="molecular weight Kg/mol")

        # Constants for specific heat capacity, enthalpy
        cp_ig_data = {('H2O', '1'): 32.24,
                      ('H2O', '2'): 1.924e-3,
                      ('H2O', '3'): 1.055e-5,
                      ('H2O', '4'): -3.596e-9,
                      ('H2O', '5'): 0}

        self.cp_ig = Param(self.component_list,
                           ['1', '2', '3', '4', '5'],
                           mutable=False,
                           initialize=extract_data(cp_ig_data),
                           doc="Parameters to compute cp_comp")

        # Antoine coefficients for ideal vapour (units: bar, K)
        # This is needed for initial guesses of bubble and dew points
        antoine_data = {('H2O', '1'): 3.55959,
                        ('H2O', '2'): 643.748,
                        ('H2O', '3'): -198.043}

        self.antoine = Param(self.component_list,
                             ['1', '2', '3'],
                             mutable=False,
                             initialize=extract_data(antoine_data),
                             doc="Antoine Parameters to pressure_sat")
