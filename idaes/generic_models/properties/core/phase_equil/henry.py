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

from idaes.generic_models.properties.core.generic.utility import (
    StateIndex)
from idaes.core.util.exceptions import BurntToast


class HenryType(Enum):
    """
    'Henry Constant' types are numbered 1-50 (i.e. H = conc/pressure)
    'Henry Volatiltiy' types are numbered 51-100 (i.e. K = pressure/conc)
    We use this fact to simplify determining wheterh to multiply or divide
    by Henry's constant
    Any differnt forms can use values 101+, but will need ot add the custom
    code in if branches where necessary
    """
    # TODO: Add more forms as needed
    Hcp = 1
    Kpc = 51
    Hxp = 2
    Kpx = 52


# Define a method for getting the correct concentration term
# TODO : Look to see if this can be unified with get_concentration_term
# in generic/utilties
def get_henry_concentration_term(blk, henry_dict, log=False):
    if log:
        pre = "log_"
    else:
        pre = ""

    if hasattr(blk.params, "_electrolyte") and blk.params._electrolyte:
        if henry_dict["basis"] == StateIndex.true:
            sub = "_true"
        else:
            sub = "_apparent"
    else:
        sub = ""

    henry_type = henry_dict["type"]
    if henry_type == HenryType.Hcp or henry_type == HenryType.Kpc:
        conc_type = "conc_mol_phase_comp"
    elif henry_type == HenryType.Hxp or henry_type == HenryType.Kpx:
        conc_type = "mole_frac_phase_comp"
    else:
        raise BurntToast(
            f"Unrecognized value for HenryType {henry_type}")

    return getattr(blk, pre+conc_type+sub)


# Define a method for returning vapor pressure of Henry components
def henry_pressure(b, p, j, T=None):
    henry_def = b.params.get_component(j).config.henry_component[p]

    if T is None:
        henry = b.henry[p, j]
    else:
        henry = henry_def["method"].return_expression(b, p, j, T)
    # Need to get the appropriate concentration term
    # TODO: Add support for true and apparent bases

    h_conc = get_henry_concentration_term(b, henry_def, log=False)[p, j]

    if henry_def["type"].value <= 50:
        # H = c/P type
        h_press = h_conc/henry
    elif henry_def["type"].value <= 100:
        # K = P/c type
        h_press = h_conc*henry
    else:
        raise BurntToast(
            f"Unrecognized value for HenryType Enum {henry_def['type']}.")

    return h_press


def log_henry_pressure(b, p, j, T=None):
    henry_def = b.params.get_component(j).config.henry_component[p]

    # TODO: Should use a log henry var/expression
    if T is None:
        henry = b.henry[p, j]
    else:
        henry = henry_def["method"].return_expression(b, p, j, T)

    # Need to get the appropriate concentration term
    # TODO: Add support for true and apparent bases

    h_conc = get_henry_concentration_term(b, henry_def, log=True)[p, j]

    if henry_def["type"].value <= 50:
        # H = c/P type
        log_h_press = h_conc - log(henry)
    elif henry_def["type"].value <= 100:
        # K = P/c type
        log_h_press = h_conc + log(henry)
    else:
        raise BurntToast(
            f"Unrecognized value for HenryType Enum {henry_def['type']}.")

    return log_h_press


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

    # TODO: Need a return log expression method too

    @staticmethod
    def dT_expression(b, p, j, T=None):
        return 0
