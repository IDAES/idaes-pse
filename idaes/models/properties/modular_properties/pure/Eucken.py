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
Eucken correlation to calculate low pressure pure component gas thermal conductivity,
(therm_cond_phase_comp) as outlined in The Properties of Gases and Liquids, 5th Ed.,
Section 10-3-1. Not particularly accurate, with errors up to 30% reported, especially
for polar compounds. The component-specific parameter f_int_eucken allows the user to
adjust the correlation on a component-by-component basis.
"""
from pyomo.environ import Var, units as pyunits
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError


class Eucken(object):
    """Class to implement Eucken correlation for gas-phase thermal conductivity"""

    class therm_cond_phase_comp(object):
        """Eucken correlation for thermal conductivity"""

        @staticmethod
        def build_parameters(cobj, p):
            """Builds dimensionless parameter f_int"""
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
            """Returns expression for therm_cond_phase_comp"""
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

            if not hasattr(b, "_visc_d_phase_comp"):
                b._make_visc_d_phase_comp()  # pylint: disable=protected-access

            if hasattr(cp_mol_ig_comp, "return_expression"):
                cp_func = cp_mol_ig_comp.return_expression
            else:
                cp_func = cp_mol_ig_comp.cp_mol_ig_comp.return_expression

            therm_cond = (
                b._visc_d_phase_comp[  # pylint: disable=protected-access
                    p, cobj.local_name
                ]
                / M
                * (f_int * cp_func(b, cobj, T) + (15 / 4 - 5 * f_int / 2) * R)
            )

            return pyunits.convert(therm_cond, units["thermal_conductivity"])
