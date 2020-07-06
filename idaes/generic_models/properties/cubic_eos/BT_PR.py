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
Example Peng-Robinson parameter block for the VLE calucations for a
benzene-toluene system.

Unless otherwise noted, parameters are from:
"The Properties of Gases and Liquids, 4th Edition", Reid, Prausnitz and Poling,
McGraw-Hill, 1987
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Pyomo libraries
from pyomo.environ import Reals, Param, NonNegativeReals, Set

# Import IDAES cores
from idaes.core import declare_process_block_class, Component
from idaes.core.util.misc import extract_data

from idaes.generic_models.properties.cubic_eos.cubic_prop_pack import (
        CubicParameterData, CubicEoS)
from idaes.logger import getIdaesLogger


# Set up logger
_log = getIdaesLogger(__name__)


@declare_process_block_class("BTParameterBlock")
class BTParameterData(CubicParameterData):

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(BTParameterData, self).build()

        self.cubic_type = CubicEoS.PR

        # Add Component objects
        self.benzene = Component()
        self.toluene = Component()

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

        # Critical Properties
        pressure_crit_data = {'benzene': 48.9e5,
                              'toluene': 41.0e5}

        self.pressure_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(pressure_crit_data),
            doc='Critical pressure [Pa]')

        temperature_crit_data = {'benzene': 562.2,
                                 'toluene': 591.8}

        self.temperature_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_crit_data),
            doc='Critical temperature [K]')

        # Pitzer acentricity factor
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

        # Molecular Weights
        mw_comp_data = {'benzene': 78.1136E-3,
                        'toluene': 92.1405E-3}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             doc="molecular weight kg/mol")

        # Constants for specific heat capacity, enthalpy and entropy
        cp_ig_data = {('benzene', '1'): -3.392E1,
                      ('benzene', '2'): 4.739E-1,
                      ('benzene', '3'): -3.017E-4,
                      ('benzene', '4'): 7.130E-8,
                      ('toluene', '1'): -2.435E1,
                      ('toluene', '2'): 5.125E-1,
                      ('toluene', '3'): -2.765E-4,
                      ('toluene', '4'): 4.911E-8}

        self.cp_ig = Param(self.component_list,
                           ['1', '2', '3', '4'],
                           mutable=False,
                           initialize=extract_data(cp_ig_data),
                           doc="Parameters to compute cp_comp")

        # Standard heats of formation
        # Source: NIST Webbook, https://webbook.nist.gov
        # Retrieved 25th September 2019
        dh_form_data = {'benzene': 82.9e3,
                        'toluene': 50.1e3}

        self.enth_mol_form_ref = Param(self.component_list,
                                       mutable=False,
                                       initialize=extract_data(dh_form_data),
                                       doc="Standard heats of formation")

        # Standard entropy of formation
        # Source: Engineering Toolbox, https://www.engineeringtoolbox.com
        # Retrieved 25th September, 2019
        ds_form_data = {'benzene': -269,
                        'toluene': -321}

        self.entr_mol_form_ref = Param(self.component_list,
                                       mutable=False,
                                       initialize=extract_data(ds_form_data),
                                       doc="Standard entropy of formation")

        # Antoine coefficients for ideal vapour (units: bar, K)
        # This is needed for initial guesses of bubble and dew points
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
