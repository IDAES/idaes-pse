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
from pyomo.environ import log, Var, units as pyunits
import pyomo.environ as pyo

from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants
from idaes.models.properties.modular_properties.base.utility import get_component_object
from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import (
    ViscosityWilkePhase,
    wilke_phi_ij_callback,
    herring_zimmer_phi_ij_callback,
)


class ThermalConductivityWMSPhase(object):
    class therm_cond_phase(object):
        @staticmethod
        def build_parameters(pobj):
            ViscosityWilkePhase.build_parameters(pobj)

        @staticmethod
        def return_expression(b, pobj):
            # Properties of Gases and Liquids, Eq. 9-5.14
            # and Eq. 10-6.4
            ViscosityWilkePhase.build_phi_ij(b, pobj)

            pname = pobj.local_name

            # Properties of Gases and Liquids, Eq. 10-6.2
            return sum(
                [
                    b.mole_frac_phase_comp[pname, i]
                    * b.therm_cond_phase_comp[pname, i]
                    / sum(
                        [
                            b.mole_frac_phase_comp[pname, j] * b.visc_d_phi_ij[i, j]
                            for j in b.components_in_phase(pname)
                        ]
                    )
                    for i in b.components_in_phase(pname)
                ]
            )
