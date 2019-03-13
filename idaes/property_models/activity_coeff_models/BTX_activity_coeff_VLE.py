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
Example property package for the VLE calucations for a Benzene-Toluene-o-Xylene
system. If using the activity coefficient models (NRTL or Wilson), the user is
expected to provide the paramters necessary for these models. Please note that
these parameters are declared as variables here to allow for use in a parameter
estimation problem if the VLE data is available.
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals, Set

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.core.util.misc import extract_data

from idaes.property_models.activity_coeff_models.activity_coeff_prop_pack \
    import ActivityCoeffParameterData

# Some more inforation about this module
__author__ = "Jaffer Ghouse"
__version__ = "0.0.1"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("BTXParameterBlock")
class BTXParameterData(ActivityCoeffParameterData):

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(BTXParameterData, self).build()

        self.component_list_master = Set(initialize=['benzene',
                                                     'toluene',
                                                     'o-xylene'])

        # Component list - a list of component identifiers
        # NOTE: User needs to update this list; can be a subset or
        # equal to the master component list
        self.component_list = Set(initialize=['benzene', 'toluene'])

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        # List of phase equilibrium index
        self.phase_equilibrium_idx_master = Set(initialize=[1, 2, 3])

        self.phase_equilibrium_idx = Set(initialize=[1, 2])

        self.phase_equilibrium_list_master = \
            {1: ["benzene", ("Vap", "Liq")],
             2: ["toluene", ("Vap", "Liq")],
             3: ["o-xylene", ("Vap", "Liq")]}

        self.phase_equilibrium_list = \
            {1: ["benzene", ("Vap", "Liq")],
             2: ["toluene", ("Vap", "Liq")]}

        # Thermodynamic reference state
        self.pressure_reference = Param(mutable=True,
                                        default=101325,
                                        doc='Reference pressure [Pa]')
        self.temperature_reference = Param(mutable=True,
                                           default=298.15,
                                           doc='Reference temperature [K]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_critical_data = {'benzene': 48.9e5,
                                  'toluene': 41e5,
                                  'o-xylene': 37.3e5
                                  }

        self.pressure_critical = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(pressure_critical_data),
            doc='Critical pressure [Pa]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        temperature_critical_data = {'benzene': 562.2,
                                     'toluene': 591.8,
                                     'o-xylene': 630.3
                                     }

        self.temperature_critical = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_critical_data),
            doc='Critical temperature [K]')

        # Gas Constant
        self.gas_const = Param(within=NonNegativeReals,
                               mutable=False,
                               default=8.314,
                               doc='Gas Constant [J/mol.K]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        mw_comp_data = {'benzene': 78.1136E-3,
                        'toluene': 92.1405E-3,
                        'o-xylene': 106.167e-3}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             doc="molecular weight Kg/mol")

        # Constants for specific heat capacity, enthalpy
        # Sources: The Properties of Gases and Liquids (1987)
        #         4th edition, Chemical Engineering Series - Robert C. Reid
        #         Perry's Chemical Engineers Handbook
        #         - Robert H. Perry (Cp_liq)
        CpIG_data = {('Liq', 'benzene', '1'): 1.29E5,
                     ('Liq', 'benzene', '2'): -1.7E2,
                     ('Liq', 'benzene', '3'): 6.48E-1,
                     ('Liq', 'benzene', '4'): 0,
                     ('Liq', 'benzene', '5'): 0,
                     ('Vap', 'benzene', '1'): -3.392E1,
                     ('Vap', 'benzene', '2'): 4.739E-1,
                     ('Vap', 'benzene', '3'): -3.017E-4,
                     ('Vap', 'benzene', '4'): 7.130E-8,
                     ('Vap', 'benzene', '5'): 0,
                     ('Liq', 'toluene', '1'): 1.40E5,
                     ('Liq', 'toluene', '2'): -1.52E2,
                     ('Liq', 'toluene', '3'): 6.95E-1,
                     ('Liq', 'toluene', '4'): 0,
                     ('Liq', 'toluene', '5'): 0,
                     ('Vap', 'toluene', '1'): -2.435E1,
                     ('Vap', 'toluene', '2'): 5.125E-1,
                     ('Vap', 'toluene', '3'): -2.765E-4,
                     ('Vap', 'toluene', '4'): 4.911E-8,
                     ('Vap', 'toluene', '5'): 0,
                     ('Liq', 'o-xylene', '1'): 3.65e4,
                     ('Liq', 'o-xylene', '2'): 1.0175e3,
                     ('Liq', 'o-xylene', '3'): -2.63,
                     ('Liq', 'o-xylene', '4'): 3.02e-3,
                     ('Liq', 'o-xylene', '5'): 0,
                     ('Vap', 'o-xylene', '1'): -1.585e-1,
                     ('Vap', 'o-xylene', '2'): 5.962e-1,
                     ('Vap', 'o-xylene', '3'): -3.443e-4,
                     ('Vap', 'o-xylene', '4'): 7.528E-8,
                     ('Vap', 'o-xylene', '5'): 0}

        self.CpIG = Param(self.phase_list, self.component_list,
                          ['1', '2', '3', '4', '5'],
                          mutable=False,
                          initialize=extract_data(CpIG_data),
                          doc="parameters to compute Cp_comp")

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_sat_coeff_data = {('benzene', 'A'): -6.98273,
                                   ('benzene', 'B'): 1.33213,
                                   ('benzene', 'C'): -2.62863,
                                   ('benzene', 'D'): -3.33399,
                                   ('toluene', 'A'): -7.28607,
                                   ('toluene', 'B'): 1.38091,
                                   ('toluene', 'C'): -2.83433,
                                   ('toluene', 'D'): -2.79168,
                                   ('o-xylene', 'A'): -7.53357,
                                   ('o-xylene', 'B'): 1.40968,
                                   ('o-xylene', 'C'): -3.10985,
                                   ('o-xylene', 'D'): -2.85992}

        self.pressure_sat_coeff = Param(
            self.component_list,
            ['A', 'B', 'C', 'D'],
            mutable=False,
            initialize=extract_data(pressure_sat_coeff_data),
            doc="parameters to compute Cp_comp")

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        dh_vap = {'benzene': 3.377e4, 'toluene': 3.8262e4,
                  'o-xylene': 4.34584e4}

        self.dh_vap = Param(self.component_list,
                            mutable=False,
                            initialize=extract_data(dh_vap),
                            doc="heat of vaporization")
