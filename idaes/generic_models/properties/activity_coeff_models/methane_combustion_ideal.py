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
from pyomo.environ import Param, Set
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBlock,
                        Component)
from idaes.core.util.misc import extract_data

from idaes.generic_models.properties.activity_coeff_models.activity_coeff_prop_pack \
    import ActivityCoeffParameterData
from idaes.logger import getIdaesLogger


# Some more inforation about this module
__author__ = "Andrew Lee, Jaffer Ghouse"
__version__ = "0.0.1"


# Set up logger
_log = getIdaesLogger(__name__)


@declare_process_block_class("MethaneParameterBlock")
class MethaneParameterData(ActivityCoeffParameterData):
    # Methane combstion only considers and ideal vapor phase, so need to
    # overload the user-selection of activity coefficient model and valid
    # phases. Do this by creating our own Config block with limited choices.
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare("activity_coeff_model", ConfigValue(
        default="Ideal",
        domain=In(["Ideal"]),
        description="Methane combustion supports ideal gas only"))

    CONFIG.declare("state_vars", ConfigValue(
        default="FTPz",
        domain=In(["FTPz", "FcTP"]),
        description="Flag indicating the choice for state variables",
        doc="""Flag indicating the choice for state variables to be used
    for the state block, and thus corresponding constraints  should be
    included,
    **default** - FTPz
    **Valid values:** {
    **"FTPx"** - Total flow, Temperature, Pressure and Mole fraction,
    **"FcTP"** - Component flow, Temperature and Pressure}"""))

    CONFIG.declare("valid_phase", ConfigValue(
        default="Vap",
        domain=In(["Vap"]),
        description="Methane combustion supports ideal gas only"))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(MethaneParameterData, self).build()

        # Component list - a list of component identifiers
        self.H2 = Component()
        self.N2 = Component()
        self.O2 = Component()
        self.CH4 = Component()
        self.CO = Component()
        self.CO2 = Component()
        self.H2O = Component()
        self.NH3 = Component()

        # List of all chemical elements that constitute the chemical species
        self.element_list = Set(initialize=['H', 'N', 'O', 'C'])

        # Elemental composition of all species
        self.element_comp = {'H2': {'H': 2, 'N': 0, 'O': 0, 'C': 0},
                             'N2': {'H': 0, 'N': 2, 'O': 0, 'C': 0},
                             'O2': {'H': 0, 'N': 0, 'O': 2, 'C': 0},
                             'CH4': {'H': 4, 'N': 0, 'O': 0, 'C': 1},
                             'CO': {'H': 0, 'N': 0, 'O': 1, 'C': 1},
                             'CO2': {'H': 0, 'N': 0, 'O': 2, 'C': 1},
                             'H2O': {'H': 2, 'N': 0, 'O': 1, 'C': 0},
                             'NH3': {'H': 3, 'N': 1, 'O': 0, 'C': 0}}

        # Thermodynamic reference state
        self.pressure_reference = Param(mutable=True,
                                        default=101325,
                                        doc='Reference pressure [Pa]')
        self.temperature_reference = Param(mutable=True,
                                           default=1500,
                                           doc='Reference temperature [K]')

        # Constants for specific heat capacity, enthalpy
        # Sources: The Properties of Gases and Liquids (1987)
        #         4th edition, Chemical Engineering Series - Robert C. Reid
        CpIG_data = {('Vap', 'CH4', 'A'): 1.925e1,
                     ('Vap', 'CH4', 'B'): 5.213e-2,
                     ('Vap', 'CH4', 'C'): 1.197e-5,
                     ('Vap', 'CH4', 'D'): -1.132e-8,
                     ('Vap', 'CH4', 'E'): 0,
                     ('Vap', 'CO', 'A'): 3.087e1,
                     ('Vap', 'CO', 'B'): -1.285e-2,
                     ('Vap', 'CO', 'C'): 2.789e-5,
                     ('Vap', 'CO', 'D'): -1.272e-8,
                     ('Vap', 'CO', 'E'): 0,
                     ('Vap', 'CO2', 'A'): 1.980e1,
                     ('Vap', 'CO2', 'B'): 7.344e-2,
                     ('Vap', 'CO2', 'C'): -5.602e-5,
                     ('Vap', 'CO2', 'D'): 1.715e-8,
                     ('Vap', 'CO2', 'E'): 0,
                     ('Vap', 'H2', 'A'): 2.714e1,
                     ('Vap', 'H2', 'B'): 9.274e-3,
                     ('Vap', 'H2', 'C'): -1.381e-5,
                     ('Vap', 'H2', 'D'): 7.645e-9,
                     ('Vap', 'H2', 'E'): 0,
                     ('Vap', 'H2O', 'A'): 3.224e1,
                     ('Vap', 'H2O', 'B'): 1.924e-3,
                     ('Vap', 'H2O', 'C'): 1.055e-5,
                     ('Vap', 'H2O', 'D'): -3.596e-9,
                     ('Vap', 'H2O', 'E'): 0,
                     ('Vap', 'N2', 'A'): 3.115e1,
                     ('Vap', 'N2', 'B'): -1.357e-2,
                     ('Vap', 'N2', 'C'): 2.680e-5,
                     ('Vap', 'N2', 'D'): -1.168e-8,
                     ('Vap', 'N2', 'E'): 0,
                     ('Vap', 'NH3', 'A'): 2.731e1,
                     ('Vap', 'NH3', 'B'): 2.383e-2,
                     ('Vap', 'NH3', 'C'): 1.707e-5,
                     ('Vap', 'NH3', 'D'): -1.185e-8,
                     ('Vap', 'NH3', 'E'): 0,
                     ('Vap', 'O2', 'A'): 2.811e1,
                     ('Vap', 'O2', 'B'): -3.680e-6,
                     ('Vap', 'O2', 'C'): 1.746e-5,
                     ('Vap', 'O2', 'D'): -1.065e-8,
                     ('Vap', 'O2', 'E'): 0}

        self.CpIG = Param(self.phase_list, self.component_list,
                          ['A', 'B', 'C', 'D', 'E'],
                          mutable=False,
                          initialize=extract_data(CpIG_data),
                          doc="parameters to compute Cp_comp")

        # Source: NIST Webbook, 9th October 2019
        dh_form = {("Vap", "CH4"): -74600,
                   ("Vap", "CO"): -110530,
                   ("Vap", "CO2"): -393520,
                   ("Vap", "H2"): 0,
                   ("Vap", "H2O"): -241830,
                   ("Vap", "N2"): 0,
                   ("Vap", "NH3"): -45900,
                   ("Vap", "O2"): 0}

        self.dh_form = Param(self.phase_list,
                             self.component_list,
                             mutable=False,
                             initialize=extract_data(dh_form),
                             doc="Heats of formation (J/mol)")

        # Source: NIST Webbook, 9th October 2019
        ds_form = {("Vap", "CH4"): 186.25,
                   ("Vap", "CO"): 197.66,
                   ("Vap", "CO2"): 213.79,
                   ("Vap", "H2"): 130.68,
                   ("Vap", "H2O"): 188.84,
                   ("Vap", "N2"): 191.61,
                   ("Vap", "NH3"): 192.77,
                   ("Vap", "O2"): 205.15}

        self.ds_form = Param(self.phase_list,
                             self.component_list,
                             mutable=False,
                             initialize=extract_data(ds_form),
                             doc="Entropies of formation (J/mol.K)")
