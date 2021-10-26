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
Methods for calculating fugacity of Henry's Law components

For now, only mole fraction basis (Kpx) form is fully supported. The remainder
is prototype code
"""
from enum import Enum

from pyomo.environ import log, Var

from idaes.core.util.exceptions import BurntToast


class HenryType(Enum):
    # TODO: ADd more forms as needed
    Hcp = 1
    Kpc = 2
    Hxp = 3
    Kpx = 4


# Define a method for returning vapor pressure of Henry components
def henry_pressure(b, p, j, T=None):
    henry_def = b.params.get_component(j).config.henry_component[p]

    if T is None:
        henry = b.henry[p, j]
    else:
        henry = henry_def["method"](b, p, j, T)
    # Need to get the appropriate concentration term
    # TODO: Add support for true and apparent bases

    if henry_def["type"] == HenryType.Hcp:
        h_press = b.conc_mole_phase_comp[p, j]/henry
    elif henry_def["type"] == HenryType.Kpc:
        h_press = b.conc_mole_phase_comp[p, j]*henry
    elif henry_def["type"] == HenryType.Hxp:
        h_press = b.mole_frac_phase_comp[p, j]/henry
    elif henry_def["type"] == HenryType.Kpx:
        h_press = b.mole_frac_phase_comp[p, j]*henry
    else:
        raise BurntToast(
            f"Unrecognized value for HenryType {henry_def['type']}")

    return h_press


def log_henry_pressure(b, p, j, T=None):
    henry_def = b.params.get_component(j).config.henry_component[p]

    if T is None:
        henry = b.henry[p, j]
    else:
        henry = henry_def["method"](b, p, j, T)

    # Need to get the appropriate concentration term
    # TODO: Add support for true and apparent bases

    if henry_def["type"] == HenryType.Hcp:
        h_press = b.log_conc_mole_phase_comp[p, j] - log(henry)
    elif henry_def["type"] == HenryType.Kpc:
        h_press = b.log_conc_mole_phase_comp[p, j] + log(henry)
    elif henry_def["type"] == HenryType.Hxp:
        h_press = b.log_mole_frac_phase_comp[p, j] - log(henry)
    elif henry_def["type"] == HenryType.Kpx:
        h_press = b.log_mole_frac_phase_comp[p, j] + log(henry)
    else:
        raise BurntToast(
            f"Unrecognized value for HenryType {henry_def['type']}")

    return h_press


# Define units for Henry's constant
def henry_units(henry_type, units):
    if henry_type == HenryType.Hcp:
        h_units = units["density_mole"]/units["pressure"]
    elif henry_type == HenryType.Kpc:
        h_units = units["pressure"]/units["density_mole"]
    elif henry_type == HenryType.Hxp:
        h_units = units["pressure"]**-1
    elif henry_type == HenryType.Kpx:
        h_units = units["pressure"]
    else:
        raise BurntToast(f"Unrecognized value for HenryType {henry_type}")

    return h_units


class ConstantH():

    @staticmethod
    def build_parameters(cobj, p, h_type):
        b = cobj.parent_block()
        units = b.get_metadata().derived_units
        h_units = henry_units(h_type, units)

        cobj.add_component(
            "henry_ref_"+p,
            Var(initialize=cobj.config.parameter_data["henry_ref"][p],
                doc="Henry coefficient (mole fraction basis) at reference "
                "state for phase "+p,
                units=h_units))

    @staticmethod
    def return_expression(b, p, j, T=None):
        cobj = b.params.get_component(j)
        H = getattr(cobj, "henry_ref_"+p)

        return H

    @staticmethod
    def dT_expression(b, p, j, T=None):
        return 0
