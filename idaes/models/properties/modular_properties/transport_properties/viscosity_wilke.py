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
Wilke rules to obtain mixture viscosity of a low pressure gas from pure component
viscosities given by visc_d_phase_comp. See Section 9-5-2 of The Properties of
Gases and Liquids, 5th Ed., for a detailed description of the method. Table 9-4
shows this method producing mixture viscosity within 5% of true viscosity, with
the notable exception of hydrogen-containing systems which have some errors up
to 12%. Those error figures are using experimentally-determined pure component
viscosity; additional error can accumulate when those viscosities must be
estimated.
"""
import pyomo.environ as pyo


class ViscosityWilke(object):
    """Class to implement Wilke method to calculate mixture dynamic viscosity"""

    @staticmethod
    def build_parameters(pobj):
        """
        Ensures there's a method to build phi_ij associated with the phase. The problem: the
        transport_properties_options dictionary may not exist, and, if it does exist, it may not have
        an explicit method assigned for the callback. Therefore create the dictionary if it doesn't exist
        and ensure the correct field is populated with the default callback.
        """
        # Have to use get method because transport_properties_options dictionary may not exist
        if (
            pobj.config.transport_property_options.get("viscosity_phi_ij_callback")
            is None
        ):
            pobj.config.transport_property_options[
                "viscosity_phi_ij_callback"
            ] = wilke_phi_ij_callback

    @staticmethod
    def build_phi_ij(b, p):
        """Method to build mixing function phi_ij"""
        pobj = b.params.get_phase(p)

        if not hasattr(b, "_visc_d_phase_comp"):
            b._make_visc_d_phase_comp()  # pylint: disable=protected-access

        if not hasattr(b, "visc_d_phi_ij"):
            mw_dict = {
                k: b.params.get_component(k).mw for k in b.components_in_phase(p)
            }

            def phi_rule(blk, i, j):
                return pobj.config.transport_property_options[
                    "viscosity_phi_ij_callback"
                ](blk, i, j, p, mw_dict)

            b.visc_d_phi_ij = pyo.Expression(
                b.components_in_phase(p),
                b.components_in_phase(p),
                rule=phi_rule,
                doc="Intermediate quantity for calculating gas mixture viscosity and thermal conductivity",
            )

    class visc_d_phase(object):
        """Method to construct gas mixture dynamic viscosity"""

        @staticmethod
        def build_parameters(pobj):
            """Build parameters"""
            ViscosityWilke.build_parameters(pobj)

        @staticmethod
        def return_expression(b, p):
            """Return expression for visc_d_phase"""
            # Properties of Gases and Liquids, Eq. 9-5.14
            ViscosityWilke.build_phi_ij(b, p)

            return sum(
                [
                    b.mole_frac_phase_comp[p, i]
                    * b._visc_d_phase_comp[p, i]  # pylint: disable=protected-access
                    / sum(
                        [
                            b.mole_frac_phase_comp[p, j] * b.visc_d_phi_ij[i, j]
                            for j in b.components_in_phase(p)
                        ]
                    )
                    for i in b.components_in_phase(p)
                ]
            )


def wilke_phi_ij_callback(b, i, j, pname, mw_dict):
    """Equation 9.5.14 from Properties of Gases and Liquids.

    .. math::

      \\phi_{ij} = \\frac{[1 + (\\mu_i/\\mu_j)^{1/2}(M_j/M_i)^{1/4}]^2}{[8(1+M_i/M_j)]^{1/2}}

    in which :math:`M_i` is the molecular weight of component :math:`i`.
    """
    visc_i = b._visc_d_phase_comp[pname, i]  # pylint: disable=protected-access
    visc_j = b._visc_d_phase_comp[pname, j]  # pylint: disable=protected-access
    return (
        1 + pyo.sqrt(visc_i / visc_j) * (mw_dict[j] / mw_dict[i]) ** 0.25
    ) ** 2 / pyo.sqrt(8 * (1 + mw_dict[i] / mw_dict[j]))


def herring_zimmer_phi_ij_callback(b, i, j, pname, mw_dict):
    """Equation 9.5.17 from Properties of Gases and Liquids.

    .. math::

      \\phi_{ij} = \\left(\\frac{M_j}{M_i}\\right)^{1/2}

    in which :math:`M_i` is the molecular weight of component :math:`i`.
    """
    return pyo.sqrt(mw_dict[j] / mw_dict[i])
