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


@declare_process_block_class("ASUParameterBlock")
class ASUParameterData(CubicParameterData):

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(ASUParameterData, self).build()
        
        self.cubic_type = CubicEoS.PR

        self.component_list = Set(initialize=['N2', 'O2', 'Ar'])

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        # List of phase equilibrium index
        self.phase_equilibrium_idx = Set(initialize=[1, 2, 3])

        self.phase_equilibrium_list = \
            {1: ["N2", ("Vap", "Liq")],
             2: ["O2", ("Vap", "Liq")],
             3: ["Ar", ("Vap", "Liq")]}

        # List of all chemical elements that constitute the chemical species 
        self.elem = ['N', 'O', 'Ar']

        # Elemental composition of all species
        self.elem_comp = {'N2'  : {'N':2, 'O':0, 'Ar':0},
                          'O2'  : {'N':0, 'O':2, 'Ar':0},
                          'Ar'  : {'N':0, 'O':0, 'Ar':1},
        }

        # Thermodynamic reference state
        self.pressure_ref = Param(mutable=True,
                                  default=101325,
                                  doc='Reference pressure [Pa]')
        self.temperature_ref = Param(mutable=True,
                                     default=298.15,
                                     doc='Reference temperature [K]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_crit_data = {'N2' : 3394387.5,
                              'O2' : 5045985.0,
                              'Ar' : 4873732.5}

        self.pressure_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(pressure_crit_data),
            doc='Critical pressure [Pa]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        temperature_crit_data = {'N2' : 126.20,
                                 'O2' : 154.58,
                                 'Ar' : 150.86}

        self.temperature_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_crit_data),
            doc='Critical temperature [K]')



        # Pitzer acentricity factor (from Prop. Gases & Liquids)
        omega_data = {'N2' : 0.040,
                      'O2' : 0.021,
                      'Ar' :-0.004}

        self.omega = Param(
            self.component_list,
            within=Reals,
            mutable=False,
            initialize=extract_data(omega_data),
            doc='Acentricity Factor')

        # Peng-Robinson binary interaction parameters
        kappa_data = {
            ('N2', 'N2'): 0.0000, ('N2', 'O2'):-0.0119, ('N2', 'Ar'):-0.0026,
            ('O2', 'N2'):-0.0119, ('O2', 'O2'): 0.0000, ('O2', 'Ar'): 0.0104,
            ('Ar', 'N2'):-0.0026, ('Ar', 'O2'): 0.0104, ('Ar', 'Ar'): 0.0000}

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
        mw_comp_data = {'Ar' : 0.04401,
                        'N2' : 0.028014,
                        'O2' : 0.03199806}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             doc="molecular weight Kg/mol")

        # Constants for specific heat capacity, enthalpy
        # Sources: The Properties of Gases and Liquids (1987)
        #         4th edition, Chemical Engineering Series - Robert C. Reid
        cp_ig_data = {('N2', '1'): 31.128960,
                      ('N2', '2'): -1.356e-2,
                      ('N2', '3'): 2.678e-5,
                      ('N2', '4'): -1.167e-8,
                      ('N2', '5'): 0,
                      ('O2', '1'): 28.087192,
                      ('O2', '2'): -3.678e-6,
                      ('O2', '3'): 1.745e-5,
                      ('O2', '4'): -1.064e-8,
                      ('O2', '5'): 0,
                      ('Ar', '1'): 20.790296,
                      ('Ar', '2'): -3.209e-5,
                      ('Ar', '3'): 5.163e-8,
                      ('Ar', '4'): 0,
                      ('Ar', '5'): 0}

        self.cp_ig = Param(self.component_list,
                           ['1', '2', '3', '4', '5'],
                           mutable=False,
                           initialize=extract_data(cp_ig_data),
                           doc="Parameters to compute cp_comp")
