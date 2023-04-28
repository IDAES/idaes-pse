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
Gas viscosity at low pressures from Chapman Enskog theory using Lennard Jones parameters.
Note that LJ parameters are underdetermined when estimated from viscosity data (see
The Indeterminacy of the Values of Potential Parameters as Derived from Transport and
Virial Coefficient by Reichenberg D., 1973 for more information) so it's important to
use LJ parameters from the same source.
"""
from pyomo.environ import Var, units as pyunits
import pyomo.environ as pyo

from idaes.core.util.misc import set_param_from_config


class ChapmanEnskogLennardJones(object):
    """Class implementing physical properties from Chapman Enskog theory"""

    @staticmethod
    def build_lennard_jones_parameters(cobj):
        """Method to build Lennard Jones parameters"""
        units = cobj.parent_block().get_metadata().derived_units
        if not hasattr(cobj, "lennard_jones_sigma"):
            cobj.lennard_jones_sigma = Var(
                doc="Parameter sigma from Lennard Jones potential",
                units=units["length"],
            )
            set_param_from_config(cobj, param="lennard_jones_sigma")

        if not hasattr(cobj, "lennard_jones_epsilon_reduced"):
            cobj.lennard_jones_epsilon_reduced = Var(
                doc="Parameter epsilon from Lennard Jones potential reduced by Boltzmann constant",
                units=units["temperature"],
            )
            set_param_from_config(cobj, param="lennard_jones_epsilon_reduced")

    class visc_d_phase_comp(object):
        """Implementation of pure component dynamic viscosity from Chapman Enskog Theory"""

        @staticmethod
        def build_parameters(cobj, p):
            """Build Lennard Jones parameters and add callback for viscosity collision integral"""
            ChapmanEnskogLennardJones.build_lennard_jones_parameters(cobj)
            if not hasattr(cobj, "viscosity_collision_integral_callback"):
                cobj.viscosity_collision_integral_callback = (
                    collision_integral_neufeld_callback
                )

        @staticmethod
        def return_expression(b, cobj, p, T):
            """Return expression for visc_d_phase_comp"""
            # Properties of Gases and Liquids, Eq. 9.3.9
            units = b.params.get_metadata().derived_units

            T = pyunits.convert(T, to_units=pyunits.K)
            sigma = pyunits.convert(cobj.lennard_jones_sigma, pyunits.angstrom)
            M = pyunits.convert(cobj.mw, pyunits.g / pyunits.mol)
            T_dim = T / pyunits.convert(
                cobj.lennard_jones_epsilon_reduced, to_units=pyunits.K
            )
            Omega = cobj.viscosity_collision_integral_callback(T_dim)

            C = (
                26.69
                * pyunits.micropoise
                * pyunits.angstrom**2
                / pyo.sqrt(pyunits.g / pyunits.mol * pyunits.K)
            )
            visc = C * pyo.sqrt(M * T) / (sigma**2 * Omega)

            return pyunits.convert(visc, units["dynamic_viscosity"])


def collision_integral_kim_ross_callback(T_dim):
    """Equation 9.4.5 from Properties of Gases and Liquids.

    .. math::

      \\Omega(T^*) =1.604/\\sqrt{T^*}

    """
    return 1.604 / pyo.sqrt(T_dim)


def collision_integral_neufeld_callback(T_dim):
    """Equation 9.4.3 from Properties of Gases and Liquids.

    .. math::

      \\Omega(T^*) = A(T^*)^{-B} +  C\\exp(-DT^*) + E\\exp(-FT^*)

    in which :math:`A=1.16145`, :math:`B=0.14874`, :math:`C=0.52487`, :math:`D=0.77320`, :math:`E=2.16178`, and :math:`F=2.43787`.
    """
    A = 1.16145
    B = 0.14874
    C = 0.52487
    D = 0.77320
    E = 2.16178
    F = 2.43787
    return A * T_dim**-B + C * pyo.exp(-D * T_dim) + E * pyo.exp(-F * T_dim)
