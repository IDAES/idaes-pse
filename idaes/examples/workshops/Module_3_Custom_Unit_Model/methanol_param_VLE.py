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
Example property package for the VLE calucations for the methanol synthesis
problem from Turkay & Grossmann. The parameters and correlations are from the
paper.
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals, Set
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import declare_process_block_class, PhysicalParameterBlock
from idaes.core.util.constants import Constants as const

from .methanol_state_block_VLE import IdealStateBlock

# Some more inforation about this module
__author__ = "Jaffer Ghouse"
__version__ = "0.0.1"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("PhysicalParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class.
    """

    # Config block for the _IdealStateBlock
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare("valid_phase", ConfigValue(
        default=('Vap', 'Liq'),
        domain=In(['Liq', 'Vap', ('Vap', 'Liq'), ('Liq', 'Vap')]),
        description="Flag indicating the valid phase",
        doc="""Flag indicating the valid phase for a given set of
    conditions, and thus corresponding constraints  should be included,
    **default** - ('Vap', 'Liq').
    **Valid values:** {
    **'Liq'** - Liquid only,
    **'Vap'** - Vapor only,
    **('Vap', 'Liq')** - Vapor-liquid equilibrium,
    **('Liq', 'Vap')** - Vapor-liquid equilibrium,}"""))

    CONFIG.declare("Cp", ConfigValue(
        default=0.035,
        domain=float,
        description="Constant pressure heat capacity in MJ/(kgmol K)",
        doc="""Value for the constant pressure heat capacity,
        **default** = 0.035 MJ/(kgmol K)"""))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PhysicalParameterData, self).build()

        self.state_block_class = IdealStateBlock

        # List of valid phases in property package
        # List of valid phases in property package
        if self.config.valid_phase == ('Liq', 'Vap') or \
                self.config.valid_phase == ('Vap', 'Liq'):
            self.phase_list = Set(initialize=['Liq', 'Vap'],
                                  ordered=True)
        elif self.config.valid_phase == 'Liq':
            self.phase_list = Set(initialize=['Liq'])
        else:
            self.phase_list = Set(initialize=['Vap'])

        self.component_list = Set(initialize=['CH4',
                                              'CO',
                                              'H2',
                                              'CH3OH'])

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        self.phase_equilibrium_idx = Set(initialize=[1, 2, 3, 4])

        self.phase_equilibrium_list = \
            {1: ["CH4", ("Vap", "Liq")],
             2: ["CO", ("Vap", "Liq")],
             3: ["H2", ("Vap", "Liq")],
             4: ["CH3OH", ("Vap", "Liq")]}

        self.vapor_pressure_coeff = {('CH4', 'A'): 15.2243,
                                     ('CH4', 'B'): 897.84,
                                     ('CH4', 'C'): -7.16,
                                     ('CO', 'A'): 14.3686,
                                     ('CO', 'B'): 530.22,
                                     ('CO', 'C'): -13.15,
                                     ('H2', 'A'): 13.6333,
                                     ('H2', 'B'): 164.9,
                                     ('H2', 'C'): 3.19,
                                     ('CH3OH', 'A'): 18.5875,
                                     ('CH3OH', 'B'): 3626.55,
                                     ('CH3OH', 'C'): -34.29}

        Cp = self.config.Cp
        Cv = Cp - const.gas_constant*1e-3
        gamma = Cp / Cv

        self.gamma = Param(within=NonNegativeReals, mutable=True, default=gamma, doc='Ratio of Cp to Cv')

        self.Cp = Param(within=NonNegativeReals, mutable=True, default=Cp, doc='Constant pressure heat capacity [MJ/(kgmol K)]')

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mol': {'method': None, 'units': 'kgmol/s'},
             'mole_frac': {'method': None, 'units': 'no unit'},
             'temperature': {'method': None, 'units': '100K'},
             'pressure': {'method': None, 'units': 'MPa'},
             'flow_mol_phase': {'method': None, 'units': 'kgmol/s'},
             'density_mol': {'method': '_density_mol',
                             'units': 'kgmol/m^3'},
             'vapor_pressure': {'method': '_vapor_pressure', 'units': 'MPa'},
             'mole_frac_phase': {'method': '_mole_frac_phase',
                                 'units': 'no unit'},
             'enthalpy_comp_liq': {'method': '_enthalpy_comp_liq',
                                   'units': 'MJ/kgmol'},
             'enthalpy_comp_vap': {'method': '_enthalpy_comp_vap',
                                   'units': 'MJ/kgmol'},
             'enthalpy_liq': {'method': '_enthalpy_liq',
                              'units': 'MJ/kgmol'},
             'enthalpy_vap': {'method': '_enthalpy_vap',
                              'units': 'MJ/kgmol'}})

        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'kg',
                               'amount': 'kgmol',
                               'temperature': '100K',
                               'energy': 'MJ',
                               'holdup': 'kgmol'})
