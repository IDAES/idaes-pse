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
Example property package for the VLE calucations for a Benzene-Toluene-o-Xylene
system. If using the activity coefficient models (NRTL or Wilson), the user is
expected to provide the paramters necessary for these models. Please note that
these parameters are declared as variables here to allow for use in a parameter
estimation problem if the VLE data is available.
"""

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals, Set

# Import IDAES cores
from idaes.core import declare_process_block_class, Component
from idaes.core.util.misc import extract_data

from idaes.generic_models.properties.activity_coeff_models.activity_coeff_prop_pack \
    import ActivityCoeffParameterData
from idaes.logger import getIdaesLogger

# Some more inforation about this module
__author__ = "Jaffer Ghouse"
__version__ = "0.0.1"


# Set up logger
_log = getIdaesLogger(__name__)


@declare_process_block_class("BTXParameterBlock")
class BTXParameterData(ActivityCoeffParameterData):

    def build(self):
        '''
        Callable method for Block construction.
        '''
        self.component_list_master = Set(initialize=['benzene',
                                                     'toluene',
                                                     'o-xylene'])

        # Create component objects
        # NOTE: User needs to update this list; can be a subset or
        # equal to the master component list
        self.benzene = Component()
        self.toluene = Component()

        super(BTXParameterData, self).build()

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
        CpIG_data = {('Liq', 'benzene', 'A'): 1.29E5,
                     ('Liq', 'benzene', 'B'): -1.7E2,
                     ('Liq', 'benzene', 'C'): 6.48E-1,
                     ('Liq', 'benzene', 'D'): 0,
                     ('Liq', 'benzene', 'E'): 0,
                     ('Vap', 'benzene', 'A'): -3.392E1,
                     ('Vap', 'benzene', 'B'): 4.739E-1,
                     ('Vap', 'benzene', 'C'): -3.017E-4,
                     ('Vap', 'benzene', 'D'): 7.130E-8,
                     ('Vap', 'benzene', 'E'): 0,
                     ('Liq', 'toluene', 'A'): 1.40E5,
                     ('Liq', 'toluene', 'B'): -1.52E2,
                     ('Liq', 'toluene', 'C'): 6.95E-1,
                     ('Liq', 'toluene', 'D'): 0,
                     ('Liq', 'toluene', 'E'): 0,
                     ('Vap', 'toluene', 'A'): -2.435E1,
                     ('Vap', 'toluene', 'B'): 5.125E-1,
                     ('Vap', 'toluene', 'C'): -2.765E-4,
                     ('Vap', 'toluene', 'D'): 4.911E-8,
                     ('Vap', 'toluene', 'E'): 0,
                     ('Liq', 'o-xylene', 'A'): 3.65e4,
                     ('Liq', 'o-xylene', 'B'): 1.0175e3,
                     ('Liq', 'o-xylene', 'C'): -2.63,
                     ('Liq', 'o-xylene', 'D'): 3.02e-3,
                     ('Liq', 'o-xylene', 'E'): 0,
                     ('Vap', 'o-xylene', 'A'): -1.585e-1,
                     ('Vap', 'o-xylene', 'B'): 5.962e-1,
                     ('Vap', 'o-xylene', 'C'): -3.443e-4,
                     ('Vap', 'o-xylene', 'D'): 7.528E-8,
                     ('Vap', 'o-xylene', 'E'): 0}

        self.CpIG = Param(self.phase_list, self.component_list,
                          ['A', 'B', 'C', 'D', 'E'],
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

        # Standard heats of formation
        # Source: NIST Webbook, https://webbook.nist.gov
        # Retrieved 25th September 2019
        dh_form_data = {('Vap', 'benzene'): 82.9e3,
                        ('Vap', 'toluene'): 50.1e3,
                        ('Vap', 'o-xylene'): 19.0e3,
                        ('Liq', 'benzene'): 49.0e3,
                        ('Liq', 'toluene'): 12.0e3,
                        ('Liq', 'o-xylene'): -24.4e3}

        self.dh_form = Param(self.phase_list,
                             self.component_list,
                             mutable=False,
                             initialize=extract_data(dh_form_data),
                             doc="Standard heats of formation [J/mol]")

        # Standard entropy of formation
        # Source: Engineering Toolbox, https://www.engineeringtoolbox.com
        # o-xylene from NIST Webbook, https://webbook.nist.gov
        # Retrieved 9th October, 2019
        ds_form_data = {('Vap', 'benzene'): -269,
                        ('Vap', 'toluene'): -321,
                        ('Vap', 'o-xylene'): -353.6,
                        ('Liq', 'benzene'): -173,
                        ('Liq', 'toluene'): -220,
                        ('Liq', 'o-xylene'): -246}

        self.ds_form = Param(self.phase_list,
                             self.component_list,
                             mutable=False,
                             initialize=extract_data(ds_form_data),
                             doc="Standard entropy of formation [J/mol.K]")
