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
Wassiljew-Mason-Saxena rules to obtain mixture thermal conductivity of a low pressure
gas from component thermal conductivities. Has a synergy with the Wilke viscosity
mixing rules because they both rely on the same phi_ij terms.

Method taken from Properties of Gases and Liquids, 5th Ed., Sections 10-6-1 and 10-6-2
"""

from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import (
    ViscosityWilke,
)


class ThermalConductivityWMS(object):
    """Class to use the Wassiljew-Mason-Saxena rules to obtain mixture thermal conductivity"""

    class therm_cond_phase(object):
        """Implement gas mixture thermal conductivity using Wassiljew-Mason-Saxena rules"""

        @staticmethod
        def build_parameters(pobj):
            """The parameters needed by this method are the same as the ViscosityWilke method"""
            ViscosityWilke.build_parameters(pobj)

        @staticmethod
        def return_expression(b, p):
            """Return expression for therm_cond_phase"""
            # Properties of Gases and Liquids, Eq. 9-5.14
            # and Eq. 10-6.4
            ViscosityWilke.build_phi_ij(b, p)
            if not hasattr(b, "_therm_cond_phase_comp"):
                b._make_therm_cond_phase_comp()  # pylint: disable=protected-access

            # Properties of Gases and Liquids, Eq. 10-6.2
            return sum(
                [
                    b.mole_frac_phase_comp[p, i]
                    * b._therm_cond_phase_comp[p, i]  # pylint: disable=protected-access
                    / sum(
                        [
                            b.mole_frac_phase_comp[p, j] * b.visc_d_phi_ij[i, j]
                            for j in b.components_in_phase(p)
                        ]
                    )
                    for i in b.components_in_phase(p)
                ]
            )
