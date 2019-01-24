##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Example property package for the VLE calucations for a Benzene-Toluene
system.
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals, Set

# Import IDAES cores
from idaes.core import declare_process_block_class, PhysicalParameterBase

from idaes.property_models.ideal_prop_pack_VLE import IdealStateBlock

# Some more inforation about this module
__author__ = "Jaffer Ghouse"
__version__ = "0.0.1"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("PhysicalParameterBlock")
class PhysicalParameterData(PhysicalParameterBase):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PhysicalParameterData, self).build()

        self.state_block_class = IdealStateBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Liq', 'Vap'])

        self.component_list_master = Set(initialize=['benzene',
                                                     'toluene',
                                                     'o-xylene'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['benzene', 'toluene',
                                              'o-xylene'])

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        self.phase_equilibrium_idx_master = Set(initialize=[1, 2, 3])

        self.phase_equilibrium_idx = Set(initialize=[1, 2, 3])

        # Reaction Stoichiometry
        self.phase_equilibrium_list_master = \
            {1: ["benzene", ("Vap", "Liq")],
             2: ["toluene", ("Vap", "Liq")],
             3: ["o-xylene", ("Vap", "Liq")]}

        self.phase_equilibrium_list = \
            {1: ["benzene", ("Vap", "Liq")],
             2: ["toluene", ("Vap", "Liq")],
             3: ["o-xylene", ("Vap", "Liq")]}

        # Thermodynamic reference state
        self.pressure_reference = Param(mutable=True,
                                        default=101325,
                                        doc='Reference pressure [Pa]')
        self.temperature_reference = Param(mutable=True,
                                           default=298.15,
                                           doc='Reference temperature [K]')

        pressure_critical_data = {'benzene': 48.9e5,
                                  'toluene': 41e5,
                                  'o-xylene': 37.3e5
                                  }

        def rule_pressure_critical(self, i):
            return pressure_critical_data[i]
        self.pressure_critical = Param(self.component_list,
                                       within=NonNegativeReals,
                                       mutable=False,
                                       initialize=rule_pressure_critical,
                                       doc='Critical pressure [Pa]')

        temperature_critical_data = {'benzene': 562.2,
                                     'toluene': 591.8,
                                     'o-xylene': 630.3
                                     }

        def rule_temperature_critical(self, i):
            return temperature_critical_data[i]
        self.temperature_critical = Param(self.component_list,
                                          within=NonNegativeReals,
                                          mutable=False,
                                          initialize=rule_temperature_critical,
                                          doc='Critical temperature [K]')

        # Gas Constant
        self.gas_constant = Param(within=NonNegativeReals,
                                  mutable=False,
                                  default=8.314,
                                  doc='Gas Constant [J/mol.K]')

        # Molecular weights
        mw_comp_data = {'benzene': 78.1136E-3,
                        'toluene': 92.1405E-3,
                        'o-xylene': 106.167e-3}

        def rule_mw(self, i):
            return mw_comp_data[i]
        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=rule_mw,
                             doc="molecular weight Kg/mol")

        # Constants for specific heat capacity, enthalpy, entropy
        # calculations for ideal gas (from NIST)
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

        def rule_cp_ig(self, i, j, k):
            return CpIG_data[i, j, k]
        self.CpIG = Param(self.phase_list, self.component_list,
                          ['1', '2', '3', '4', '5'],
                          mutable=False,
                          initialize=rule_cp_ig,
                          doc="parameters to compute Cp_comp")

        vapor_pressure_coeff_data = {('benzene', 'A'): -6.98273,
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

        def rule_vap_pressure(self, i, j):
            return vapor_pressure_coeff_data[i, j]
        self.vapor_pressure_coeff = Param(self.component_list,
                                          ['A', 'B', 'C', 'D'],
                                          mutable=False,
                                          initialize=rule_vap_pressure,
                                          doc="parameters to compute Cp_comp")

        delH_vap = {'benzene': 3.377e4, 'toluene': 3.8262e4,
                    'o-xylene': 4.34584e4}

        def rule_hvap(self, i):
            return delH_vap[i]
        self.delH_vap = Param(self.component_list,
                              mutable=False,
                              initialize=rule_hvap,
                              doc="heat of vaporization")

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mol': {'method': None, 'units': 'mol/s'},
             'mole_frac': {'method': None, 'units': 'no unit'},
             'temperature': {'method': None, 'units': 'K'},
             'pressure': {'method': None, 'units': 'Pa'},
             'flow_mol_phase': {'method': None, 'units': 'mol/s'},
             'density_mol': {'method': '_density_mol',
                             'units': 'mol/m^3'},
             'vapor_pressure': {'method': '_vapor_pressure', 'units': 'Pa'},
             'mole_frac_phase': {'method': '_mole_frac_phase',
                                 'units': 'no unit'},
             'enthalpy_comp_liq': {'method': '_enthalpy_comp_liq',
                                   'units': 'J/mol'},
             'enthalpy_comp_vap': {'method': '_enthalpy_comp_vap',
                                   'units': 'J/mol'},
             'enthalpy_liq': {'method': '_enthalpy_liq',
                              'units': 'J/mol'},
             'enthalpy_vap': {'method': '_enthalpy_vap',
                              'units': 'J/mol'}})

        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})
