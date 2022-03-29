#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Example property package for the VLE calucations for a Benzene-Toluene-o-Xylene
system. If using the activity coefficient models (NRTL or Wilson), the user is
expected to provide the paramters necessary for these models. Please note that
these parameters are declared as variables here to allow for use in a parameter
estimation problem if the VLE data is available.
"""

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals, Set, units as pyunits

# Import IDAES cores
from idaes.core import declare_process_block_class, Component
from idaes.core.util.misc import extract_data

from idaes.models.properties.activity_coeff_models.activity_coeff_prop_pack import (
    ActivityCoeffParameterData,
)
from idaes.logger import getIdaesLogger

# Some more inforation about this module
__author__ = "Jaffer Ghouse"
__version__ = "0.0.1"


# Set up logger
_log = getIdaesLogger(__name__)


@declare_process_block_class("BTXParameterBlock")
class BTXParameterData(ActivityCoeffParameterData):
    def build(self):
        """
        Callable method for Block construction.
        """
        self.component_list_master = Set(initialize=["benzene", "toluene", "o-xylene"])

        # Create component objects
        # NOTE: User needs to update this list; can be a subset or
        # equal to the master component list
        self.benzene = Component()
        self.toluene = Component()

        super(BTXParameterData, self).build()

        # List of phase equilibrium index
        self.phase_equilibrium_idx_master = Set(initialize=[1, 2, 3])

        self.phase_equilibrium_idx = Set(initialize=[1, 2])

        self.phase_equilibrium_list_master = {
            1: ["benzene", ("Vap", "Liq")],
            2: ["toluene", ("Vap", "Liq")],
            3: ["o-xylene", ("Vap", "Liq")],
        }

        self.phase_equilibrium_list = {
            1: ["benzene", ("Vap", "Liq")],
            2: ["toluene", ("Vap", "Liq")],
        }

        # Thermodynamic reference state
        self.pressure_reference = Param(
            mutable=True,
            default=101325,
            doc="Reference pressure [Pa]",
            units=pyunits.Pa,
        )
        self.temperature_reference = Param(
            mutable=True,
            default=298.15,
            doc="Reference temperature [K]",
            units=pyunits.K,
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_critical_data = {
            "benzene": 48.9e5,
            "toluene": 41e5,
            "o-xylene": 37.3e5,
        }

        self.pressure_critical = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(pressure_critical_data),
            doc="Critical pressure [Pa]",
            units=pyunits.Pa,
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        temperature_critical_data = {
            "benzene": 562.2,
            "toluene": 591.8,
            "o-xylene": 630.3,
        }

        self.temperature_critical = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_critical_data),
            doc="Critical temperature [K]",
            units=pyunits.K,
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        mw_comp_data = {
            "benzene": 78.1136e-3,
            "toluene": 92.1405e-3,
            "o-xylene": 106.167e-3,
        }

        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize=extract_data(mw_comp_data),
            doc="molecular weight kg/mol",
            units=pyunits.kg / pyunits.mol,
        )

        # Constants for specific heat capacity, enthalpy
        # Sources: The Properties of Gases and Liquids (1987)
        #         4th edition, Chemical Engineering Series - Robert C. Reid
        #         Perry's Chemical Engineers Handbook
        #         - Robert H. Perry (Cp_liq)
        Cp_Liq_A_data = {("benzene"): 1.29e5, ("toluene"): 1.40e5, ("o-xylene"): 3.65e4}
        Cp_Liq_B_data = {
            ("benzene"): -1.7e2,
            ("toluene"): -1.52e2,
            ("o-xylene"): 1.0175e3,
        }
        Cp_Liq_C_data = {
            ("benzene"): 6.48e-1,
            ("toluene"): 6.95e-1,
            ("o-xylene"): -2.63,
        }
        Cp_Liq_D_data = {("benzene"): 0, ("toluene"): 0, ("o-xylene"): 3.02e-3}
        Cp_Liq_E_data = {("benzene"): 0, ("toluene"): 0, ("o-xylene"): 0}

        self.cp_mol_liq_comp_coeff_A = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K,
            doc="Liquid phase Cp parameter A",
            initialize=extract_data(Cp_Liq_A_data),
        )

        self.cp_mol_liq_comp_coeff_B = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K**2,
            doc="Liquid phase Cp parameter B",
            initialize=extract_data(Cp_Liq_B_data),
        )

        self.cp_mol_liq_comp_coeff_C = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K**3,
            doc="Liquid phase Cp parameter C",
            initialize=extract_data(Cp_Liq_C_data),
        )

        self.cp_mol_liq_comp_coeff_D = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K**4,
            doc="Liquid phase Cp parameter D",
            initialize=extract_data(Cp_Liq_D_data),
        )

        self.cp_mol_liq_comp_coeff_E = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K**5,
            doc="Liquid phase Cp parameter E",
            initialize=extract_data(Cp_Liq_E_data),
        )

        Cp_Vap_A_data = {
            ("benzene"): -3.392e1,
            ("toluene"): -2.435e1,
            ("o-xylene"): -1.585e-1,
        }
        Cp_Vap_B_data = {
            ("benzene"): 4.739e-1,
            ("toluene"): 5.125e-1,
            ("o-xylene"): 5.962e-1,
        }
        Cp_Vap_C_data = {
            ("benzene"): -3.017e-4,
            ("toluene"): -2.765e-4,
            ("o-xylene"): -3.443e-4,
        }
        Cp_Vap_D_data = {
            ("benzene"): 7.130e-8,
            ("toluene"): 4.911e-8,
            ("o-xylene"): 7.528e-8,
        }
        Cp_Vap_E_data = {("benzene"): 0, ("toluene"): 0, ("o-xylene"): 0}

        self.cp_mol_vap_comp_coeff_A = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Vapor phase Cp parameter A",
            initialize=extract_data(Cp_Vap_A_data),
        )

        self.cp_mol_vap_comp_coeff_B = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
            doc="Vapor phase Cp parameter B",
            initialize=extract_data(Cp_Vap_B_data),
        )

        self.cp_mol_vap_comp_coeff_C = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K**3,
            doc="Vapor phase Cp parameter C",
            initialize=extract_data(Cp_Vap_C_data),
        )

        self.cp_mol_vap_comp_coeff_D = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K**4,
            doc="Vapor phase Cp parameter D",
            initialize=extract_data(Cp_Vap_D_data),
        )

        self.cp_mol_vap_comp_coeff_E = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K**5,
            doc="Vapor phase Cp parameter E",
            initialize=extract_data(Cp_Vap_E_data),
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_sat_coeff_data = {
            ("benzene", "A"): -6.98273,
            ("benzene", "B"): 1.33213,
            ("benzene", "C"): -2.62863,
            ("benzene", "D"): -3.33399,
            ("toluene", "A"): -7.28607,
            ("toluene", "B"): 1.38091,
            ("toluene", "C"): -2.83433,
            ("toluene", "D"): -2.79168,
            ("o-xylene", "A"): -7.53357,
            ("o-xylene", "B"): 1.40968,
            ("o-xylene", "C"): -3.10985,
            ("o-xylene", "D"): -2.85992,
        }

        self.pressure_sat_coeff = Param(
            self.component_list,
            ["A", "B", "C", "D"],
            mutable=False,
            initialize=extract_data(pressure_sat_coeff_data),
            doc="parameters to compute P_sat",
        )

        # Standard heats of formation
        # Source: NIST Webbook, https://webbook.nist.gov
        # Retrieved 25th September 2019
        dh_form_data = {
            ("Vap", "benzene"): 82.9e3,
            ("Vap", "toluene"): 50.1e3,
            ("Vap", "o-xylene"): 19.0e3,
            ("Liq", "benzene"): 49.0e3,
            ("Liq", "toluene"): 12.0e3,
            ("Liq", "o-xylene"): -24.4e3,
        }

        self.dh_form = Param(
            self.phase_list,
            self.component_list,
            mutable=False,
            initialize=extract_data(dh_form_data),
            doc="Standard heats of formation [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        # Standard entropy of formation
        # Source: Engineering Toolbox, https://www.engineeringtoolbox.com
        # o-xylene from NIST Webbook, https://webbook.nist.gov
        # Retrieved 9th October, 2019
        ds_form_data = {
            ("Vap", "benzene"): -269,
            ("Vap", "toluene"): -321,
            ("Vap", "o-xylene"): -353.6,
            ("Liq", "benzene"): -173,
            ("Liq", "toluene"): -220,
            ("Liq", "o-xylene"): -246,
        }

        self.ds_form = Param(
            self.phase_list,
            self.component_list,
            mutable=False,
            initialize=extract_data(ds_form_data),
            doc="Standard entropy of formation [J/mol.K]",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
