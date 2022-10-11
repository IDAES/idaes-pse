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
Method to set constant pure component properties:

"""
from pyomo.environ import log, Var, units as pyunits
import pyomo.environ as pyo
from idaes.models.properties.modular_properties.base.utility import get_method
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError


# ------------------------------------------------------------------------------------
# Gas viscosity at low pressures from Lennard Jones paramters. Note that LJ parameters
# are underdetermined when estimated from viscosity data (see The Indeterminacy of
# the Values of Potential Parameters as Derived from Transport and Virial Coefficient
# by Reichenberg D., 1973 for more information) so it's important to use LJ parameters
# from the same source.
class Eucken(object):
    class therm_cond_phase_comp(object):
        @staticmethod
        def build_parameters(cobj, p):
            # pobj = cobj.parent_block().get_phase(p)
            # if not pobj.is_vapor_phase:
            #     raise ConfigurationError(f"The Eucken method works only for vapor phases, not phase {p}")

            # Properties of Gases and Liquids 10-3-1: Three common values for f_int: f_int=1, which is Eucken
            # classic, f_int=1.32, which is "Modified Eucken", and f_int=1.15, which is a compromise suggested by
            # Stiel, L. I., and G. Thodos:AIChE J., 10: 26 (1964)
            if not hasattr(cobj, "f_int_eucken"):
                cobj.f_int_eucken = Var(
                    doc="Dimensionless factor in Eucken formula associated with internal degrees of freedom",
                    units=pyunits.dimensionless,
                )
                set_param_from_config(cobj, param="f_int_eucken")

        @staticmethod
        def return_expression(b, cobj, p, T):
            # Properties of Gases and Liquids, Eq. 10-3.2
            units = b.params.get_metadata().derived_units

            M = pyunits.convert(cobj.mw, pyunits.kg / pyunits.mol)
            R = pyunits.convert(Constants.gas_constant, units["heat_capacity_mole"])
            f_int = cobj.f_int_eucken
            try:
                cp_mol_ig_comp = cobj.config.cp_mol_ig_comp
            except AttributeError:
                raise ConfigurationError(
                    f"Cannot find method to calculate cp_mol_ig_comp for component "
                    f"{cobj.local_name}."
                )

            if hasattr(cp_mol_ig_comp, "return_expression"):
                cp_func = cp_mol_ig_comp.return_expression
            else:
                cp_func = cp_mol_ig_comp.cp_mol_ig_comp.return_expression

            therm_cond = (
                b.visc_d_phase_comp[p, cobj.local_name]
                / M
                * (f_int * cp_func(b, cobj, T) + (15 / 4 - 5 * f_int / 2) * R)
            )

            return pyunits.convert(therm_cond, units["thermal_conductivity"])
