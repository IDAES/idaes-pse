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

from idaes.models.properties.modular_properties.base.utility import StateIndex
from idaes.core.util.exceptions import ConfigurationError


class HenryType(Enum):
    """
    'Henry Constant' types are numbered 1-50 (i.e. H = conc/pressure)
    'Henry Volatiltiy' types are numbered 51-100 (i.e. K = pressure/conc)
    We use this fact to simplify determining wheterh to multiply or divide
    by Henry's constant
    Any different forms can use values 101+, but will need to add the custom
    code in if branches where necessary
    """

    # TODO: Add more forms as needed
    Hcp = 1
    Kpc = 51
    Hxp = 2
    Kpx = 52
    Dummy = 999  # To test error handling


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
        _raise_henry_type_error(henry_type)

    return getattr(blk, pre + conc_type + sub)


# TODO pressure -> fugacity
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
        h_press = h_conc / henry
    elif henry_def["type"].value <= 100:
        # K = P/c type
        h_press = h_conc * henry
    else:
        _raise_henry_type_error(henry_def["type"])

    return h_press


# TODO pressure -> fugacity
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
        _raise_henry_type_error(henry_def["type"])

    return log_h_press


def henry_equilibrium_ratio(b, p, j):
    """
    Returns vapor/liquid equilibrium mole ratio of Henry component j at the
    temperature and pressure of block b.

    Arguments:
        b: Property block for which this calculation is taking place
        p: Liquid phase for which this ratio is being calculated
        j: Henry component in phase p

    Returns:
        Molar ratio of component j in the vapor phase to that in liquid phase p

    Notes:
        If Henry's law is defined in terms of mole fraction, this method
        gives a meaningful answer whether or not the liquid phase composition
        has been set. If it is defined in terms of concentration, a reasonable
        value is given only if a reasonable liquid phase composition has been
        specified.

        This function returns a value regardless of whether or not block b
        is experiencing phase equilibrium at its current conditions.
    """

    henry_def = b.params.get_component(j).config.henry_component[p]
    henry_constant = b.henry[p, j]
    if henry_def["type"] == HenryType.Hcp:
        henry_constant /= b.dens_mol_phase[p]
    elif henry_def["type"] == HenryType.Kpc:
        henry_constant *= b.dens_mol_phase[p]

    if henry_def["type"].value <= 50:
        # H = c/P type
        return 1 / (henry_constant * b.pressure)
    elif henry_def["type"].value <= 100:
        # K = P/c type
        return henry_constant / b.pressure
    else:
        _raise_henry_type_error(henry_def["type"])


# Define units for Henry's constant
def henry_units(henry_type, units):
    if henry_type == HenryType.Hcp:
        h_units = units.DENSITY_MOLE / units.PRESSURE
    elif henry_type == HenryType.Kpc:
        h_units = units.PRESSURE / units.DENSITY_MOLE
    elif henry_type == HenryType.Hxp:
        h_units = units.PRESSURE**-1
    elif henry_type == HenryType.Kpx:
        h_units = units.PRESSURE
    else:
        _raise_henry_type_error(henry_type)

    return h_units


class ConstantH:
    @staticmethod
    def build_parameters(cobj, p, h_type):
        b = cobj.parent_block()
        units = b.get_metadata().derived_units
        h_units = henry_units(h_type, units)

        cobj.add_component(
            "henry_ref_" + p,
            Var(
                initialize=cobj.config.parameter_data["henry_ref"][p],
                doc="Henry coefficient (mole fraction basis) at reference "
                "state for phase " + p,
                units=h_units,
            ),
        )

    @staticmethod
    def return_expression(b, p, j, T=None):
        cobj = b.params.get_component(j)
        H = getattr(cobj, "henry_ref_" + p)

        return H

    # TODO: Need a return log expression method too

    @staticmethod
    def dT_expression(b, p, j, T=None):
        return 0


def _raise_henry_type_error(henry_type):
    raise ConfigurationError(f"Unrecognized value for HenryType {henry_type}")
