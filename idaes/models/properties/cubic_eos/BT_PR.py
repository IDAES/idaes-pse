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
Example Peng-Robinson parameter block for the VLE calucations for a
benzene-toluene system.

Unless otherwise noted, parameters are from:
"The Properties of Gases and Liquids, 4th Edition", Reid, Prausnitz and Poling,
McGraw-Hill, 1987
"""
# Import Pyomo libraries
from pyomo.environ import Reals, Param, NonNegativeReals, Set, units as pyunits

# Import IDAES cores
from idaes.core import declare_process_block_class, Component
from idaes.core.util.misc import extract_data

from idaes.models.properties.cubic_eos.cubic_prop_pack import (
    CubicParameterData,
    CubicEoS,
)
from idaes.logger import getIdaesLogger


# Set up logger
_log = getIdaesLogger(__name__)


@declare_process_block_class("BTParameterBlock")
class BTParameterData(CubicParameterData):
    def build(self):
        """
        Callable method for Block construction.
        """
        super(BTParameterData, self).build()

        self.cubic_type = CubicEoS.PR

        # Add Component objects
        self.benzene = Component(elemental_composition={"C": 6, "H": 6})
        self.toluene = Component(elemental_composition={"C": 7, "H": 8})

        # List of phase equilibrium index
        self.phase_equilibrium_idx = Set(initialize=[1, 2])

        self.phase_equilibrium_list = {
            1: ["benzene", ("Vap", "Liq")],
            2: ["toluene", ("Vap", "Liq")],
        }

        # Thermodynamic reference state
        self.pressure_ref = Param(
            mutable=True,
            default=101325,
            doc="Reference pressure [Pa]",
            units=pyunits.Pa,
        )
        self.temperature_ref = Param(
            mutable=True,
            default=298.15,
            doc="Reference temperature [K]",
            units=pyunits.K,
        )

        # Critical Properties
        pressure_crit_data = {"benzene": 48.9e5, "toluene": 41.0e5}

        self.pressure_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(pressure_crit_data),
            doc="Critical pressure [Pa]",
            units=pyunits.Pa,
        )

        temperature_crit_data = {"benzene": 562.2, "toluene": 591.8}

        self.temperature_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_crit_data),
            doc="Critical temperature [K]",
            units=pyunits.K,
        )

        # Pitzer acentricity factor
        omega_data = {"benzene": 0.212, "toluene": 0.263}

        self.omega = Param(
            self.component_list,
            within=Reals,
            mutable=False,
            initialize=extract_data(omega_data),
            doc="Acentricity Factor",
        )

        # Peng-Robinson binary interaction parameters
        kappa_data = {
            ("benzene", "benzene"): 0.0000,
            ("benzene", "toluene"): 0.0000,
            ("toluene", "benzene"): 0.0000,
            ("toluene", "toluene"): 0.0000,
        }

        self.kappa = Param(
            self.component_list,
            self.component_list,
            within=Reals,
            mutable=False,
            initialize=extract_data(kappa_data),
            doc="Peng-Robinson binary interaction parameters",
        )

        # Molecular Weights
        mw_comp_data = {"benzene": 78.1136e-3, "toluene": 92.1405e-3}

        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize=extract_data(mw_comp_data),
            doc="molecular weight kg/mol",
            units=pyunits.kg / pyunits.mol,
        )

        # Constants for specific heat capacity, enthalpy and entropy
        self.cp_mol_ig_comp_coeff_1 = Param(
            self.component_list,
            mutable=False,
            initialize={"benzene": -3.392e1, "toluene": -2.435e1},
            doc="Parameter 1 to compute cp_mol_comp",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )

        self.cp_mol_ig_comp_coeff_2 = Param(
            self.component_list,
            mutable=False,
            initialize={"benzene": 4.739e-1, "toluene": 5.125e-1},
            doc="Parameter 2 to compute cp_mol_comp",
            units=pyunits.J / pyunits.mol / pyunits.K**2,
        )

        self.cp_mol_ig_comp_coeff_3 = Param(
            self.component_list,
            mutable=False,
            initialize={"benzene": -3.017e-4, "toluene": -2.765e-4},
            doc="Parameter 3 to compute cp_mol_comp",
            units=pyunits.J / pyunits.mol / pyunits.K**3,
        )

        self.cp_mol_ig_comp_coeff_4 = Param(
            self.component_list,
            mutable=False,
            initialize={"benzene": 7.130e-8, "toluene": 4.911e-8},
            doc="Parameter 4 to compute cp_mol_comp",
            units=pyunits.J / pyunits.mol / pyunits.K**4,
        )

        # Standard heats of formation
        # Source: NIST Webbook, https://webbook.nist.gov
        # Retrieved 25th September 2019
        dh_form_data = {"benzene": 82.9e3, "toluene": 50.1e3}

        self.enth_mol_form_ref = Param(
            self.component_list,
            mutable=False,
            initialize=extract_data(dh_form_data),
            doc="Standard heats of formation",
            units=pyunits.J / pyunits.mol,
        )

        # Standard entropy of formation
        # Source: Engineering Toolbox, https://www.engineeringtoolbox.com
        # Retrieved 25th September, 2019
        ds_form_data = {"benzene": -269, "toluene": -321}

        self.entr_mol_form_ref = Param(
            self.component_list,
            mutable=False,
            initialize=extract_data(ds_form_data),
            doc="Standard entropy of formation",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )

        # Antoine coefficients for ideal vapour (units: bar, K)
        # This is needed for initial guesses of bubble and dew points
        self.antoine_coeff_A = Param(
            self.component_list,
            mutable=False,
            initialize={"benzene": 4.202, "toluene": 4.216},
            doc="Antoine A Parameter to calculate pressure_sat",
            units=pyunits.dimensionless,
        )

        self.antoine_coeff_B = Param(
            self.component_list,
            mutable=False,
            initialize={"benzene": 1322, "toluene": 1435},
            doc="Antoine B Parameter to calculate pressure_sat",
            units=pyunits.K,
        )

        self.antoine_coeff_C = Param(
            self.component_list,
            mutable=False,
            initialize={"benzene": -38.56, "toluene": -43.33},
            doc="Antoine C Parameter to calculate pressure_sat",
            units=pyunits.K,
        )
