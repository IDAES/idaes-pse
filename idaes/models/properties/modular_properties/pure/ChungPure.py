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
Pure component gas viscosities at low pressures from critical properties,
acentric factor, and dipole moment. Also requires an association factor, which
can be found for some highly polar substances in Table 9-1 in the Properties
of Gases and Liquids, 5th Ed. Chung et al. might also have additional factors
in some of their papers. If unknown, set the association factor to zero.
"""
from pyomo.environ import Var, units as pyunits
import pyomo.environ as pyo

from idaes.core.util.misc import set_param_from_config
from idaes.models.properties.modular_properties.pure.ChapmanEnskog import (
    collision_integral_neufeld_callback,
)


class ChungViscosityPure(object):
    """Container class for physical property calculations via methods of Chung et al."""

    @staticmethod
    def build_common_parameters(cobj):
        """Build dipole moment and association factor"""
        units = cobj.parent_block().get_metadata().derived_units
        if not hasattr(cobj, "dipole_moment"):
            cobj.dipole_moment = Var(
                doc="Molecular dipole moment",
                units=units["current"] * units["time"] * units["length"],
            )
            set_param_from_config(cobj, param="dipole_moment")

        if not hasattr(cobj, "association_factor_chung"):
            cobj.association_factor_chung = Var(
                doc="Correction term used for highly polar components",
                units=pyunits.dimensionless,
            )
            set_param_from_config(cobj, param="association_factor_chung")

    class visc_d_phase_comp(object):
        """
        Implementation of pure component gas dynamic viscosity via the method of Chung et al. as described in The Properties
         of Gases and Liquids, Section 9-4-2.
        """

        @staticmethod
        def build_parameters(cobj, p):
            """Build common parameters and viscosity collision integral callback"""
            ChungViscosityPure.build_common_parameters(cobj)
            if not hasattr(cobj, "viscosity_collision_integral_callback"):
                cobj.viscosity_collision_integral_callback = (
                    collision_integral_neufeld_callback
                )

        @staticmethod
        def return_expression(b, cobj, p, T):
            """Return expression for visc_d_phase_comp"""
            # Properties of Gases and Liquids 5th Ed., Section 9-4-2
            units = b.params.get_metadata().derived_units

            T = pyunits.convert(T, to_units=pyunits.K)
            T_crit = pyunits.convert(cobj.temperature_crit, to_units=pyunits.K)
            V_crit = 1 / pyunits.convert(
                cobj.dens_mol_crit, to_units=pyunits.mol / pyunits.mL
            )
            M = pyunits.convert(cobj.mw, pyunits.g / pyunits.mol)
            mu = pyunits.convert(cobj.dipole_moment, pyunits.debye)
            omega = cobj.omega
            kappa = cobj.association_factor_chung
            Omega = cobj.viscosity_collision_integral_callback(
                1.2593 * T / T_crit
            )  # Eq. 9-4.10
            # Eq. 9-4.12
            mu_r = (
                131.3
                * pyo.sqrt(pyunits.mL / pyunits.mol * pyunits.K)
                / pyunits.debye
                * (mu / pyo.sqrt(V_crit * T_crit))
            )
            # Eq. 9-4.11
            Fc = 1 - 0.2756 * omega + 0.059035 * mu_r**4 + kappa
            C = (
                40.785
                * pyunits.micropoise
                * (pyunits.mL / pyunits.mol) ** (2 / 3)
                / pyo.sqrt(pyunits.g / pyunits.mol * pyunits.K)
            )
            # Eq. 9-4.10
            visc = C * Fc * pyo.sqrt(M * T) / (V_crit ** (2 / 3) * Omega)

            return pyunits.convert(visc, units["dynamic_viscosity"])
