#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""Generic Helmholtz EOS Functions and Parameters
"""
__author__ = "John Eslick"

import enum
import ctypes
import os

from matplotlib import pyplot as plt
import numpy as np

import pyomo.environ as pyo
from pyomo.common.fileutils import find_library
from pyomo.common.config import ConfigValue, In
import idaes
from idaes.core.util.exceptions import ConfigurationError
from idaes.core import declare_process_block_class
from idaes.core import (
    PhysicalParameterBlock,
    LiquidPhase,
    VaporPhase,
    Phase,
    Component,
)
from idaes.models.properties.general_helmholtz.components import (
    get_transport_module,
    component_registered,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

_data_dir = os.path.join(idaes.bin_directory, "helm_data")
_data_dir = os.path.join(_data_dir, "")

try:
    # When compiling these, I don't bother changing the extension based on OS,
    # so the file name is always ends in .so. It's fine.
    _flib = find_library("general_helmholtz_external.so")
    ctypes.cdll.LoadLibrary(_flib)
except:
    _flib = None


def helmholtz_available():
    """Returns True if the shared library is installed and loads properly
    otherwise returns False
    """
    if _flib is None:
        return False
    if not os.path.exists(_data_dir):
        _log.error(f"The Helmholtz EoS data directory {_data_dir} does not exist.")
        return False
    return True


helmholtz_data_dir = _data_dir


class StateVars(enum.Enum):
    """
    State variable set options
    """

    PH = 1  # Pressure, Enthalpy
    PS = 2  # Pressure, Entropy
    PU = 3  # Pressure, Internal Energy
    TPX = 4  # Temperature, Pressure, Quality


class PhaseType(enum.Enum):
    """
    Ways to present phases to the framework
    """

    MIX = 1  # Looks like a single phase called mixed with a vapor fraction
    LG = 2  # Looks like two phases vapor and liquid
    L = 3  # Assume only liquid is present
    G = 4  # Assume only vapor is present


class AmountBasis(enum.Enum):
    """
    Mass or mole basis
    """

    MOLE = 1
    MASS = 2


dimensionless = pyo.units.dimensionless


_external_function_map = {
    # Thermo properties as a function of delta and tau
    "p_func": {  # pressure
        "fname": "p",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "p(comp, delta, tau)",
    },
    "u_func": {  # internal energy
        "fname": "u",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "u(comp, delta, tau)",
    },
    "s_func": {  # entropy
        "fname": "s",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "s(comp, delta, tau)",
    },
    "h_func": {  # enthalpy
        "fname": "h",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "h(comp, delta, tau)",
    },
    "g_func": {  # Gibbs free energy
        "fname": "g",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "g(comp, delta, tau)",
    },
    "f_func": {  # Helmholtz free energy
        "fname": "f",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "f(comp, delta, tau)",
    },
    "cv_func": {  # constant volume heat capacity
        "fname": "cv",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "cv(comp, delta, tau)",
    },
    "cp_func": {  # constant pressure heat capacity
        "fname": "cp",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "cp(comp, delta, tau)",
    },
    "w_func": {  # speed of sound
        "fname": "w",
        "units": pyo.units.m / pyo.units.s,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "w(comp, delta, tau)",
    },
    "v_func": {  # specific volume
        "fname": "v",
        "units": pyo.units.m**3 / pyo.units.kg,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "v(comp, delta, tau)",
    },
    # Functions of (h, p)
    "u_hp_func": {  # internal energy
        "fname": "u_hp",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "u(comp, enthalpy, pressure)",
    },
    "s_hp_func": {  # entropy
        "fname": "s_hp",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "s(comp, enthalpy, pressure)",
    },
    "g_hp_func": {  # Gibbs free energy
        "fname": "g_hp",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "g(comp, enthalpy, pressure)",
    },
    "f_hp_func": {  # Helmholtz free energy
        "fname": "f_hp",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "f(comp, enthalpy, pressure)",
    },
    "cv_hp_func": {  # constant volume heat capacity
        "fname": "cv_hp",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "cv(comp, enthalpy, pressure)",
    },
    "cp_hp_func": {  # constant pressure heat capacity
        "fname": "cp_hp",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "cp(comp, enthalpy, pressure)",
    },
    "w_hp_func": {  # speed of sound
        "fname": "w_hp",
        "units": pyo.units.m / pyo.units.s,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "w(comp, enthalpy, pressure)",
    },
    "v_hp_func": {  # specific volume
        "fname": "v_hp",
        "units": pyo.units.m**3 / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "v(comp, enthalpy, pressure)",
    },
    # Functions of (s, p)
    "u_sp_func": {  # internal energy
        "fname": "u_sp",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "u(comp, entropy, pressure)",
    },
    "h_sp_func": {  # enthalpy
        "fname": "h_sp",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "h(comp, entropy, pressure)",
    },
    "g_sp_func": {  # Gibbs free energy
        "fname": "g_sp",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "g(comp, entropy, pressure)",
    },
    "f_sp_func": {  # Helmholtz free energy
        "fname": "f_sp",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "f(comp, entropy, pressure)",
    },
    "cv_sp_func": {  # constant volume heat capacity
        "fname": "cv_sp",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "cv(comp, entropy, pressure)",
    },
    "cp_sp_func": {  # constant pressure heat capacity
        "fname": "cp_sp",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "cp(comp, entropy, pressure)",
    },
    "w_sp_func": {  # speed of sound
        "fname": "w_sp",
        "units": pyo.units.m / pyo.units.s,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "w(comp, entropy, pressure)",
    },
    "v_sp_func": {  # specific volume
        "fname": "v_sp",
        "units": pyo.units.m**3 / pyo.units.kg,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "v(comp, entropy, pressure)",
    },
    # Functions of (u, p)
    "h_up_func": {  # enthalpy
        "fname": "h_up",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "h(comp, internal energy, pressure)",
    },
    "s_up_func": {  # entropy
        "fname": "s_up",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "s(comp, internal energy, pressure)",
    },
    "g_up_func": {  # Gibbs free energy
        "fname": "g_up",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "g(comp, internal energy, pressure)",
    },
    "f_up_func": {  # Helmholtz free energy
        "fname": "f_up",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "f(comp, internal energy, pressure)",
    },
    "cv_up_func": {  # constant volume heat capacity
        "fname": "cv_up",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "cv(comp, internal energy, pressure)",
    },
    "cp_up_func": {  # constant pressure heat capacity
        "fname": "cp_up",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "cp(comp, internal energy, pressure)",
    },
    "w_up_func": {  # speed of sound
        "fname": "w_up",
        "units": pyo.units.m / pyo.units.s,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "w(comp, internal energy, pressure)",
    },
    "v_up_func": {  # specific volume
        "fname": "v_up",
        "units": pyo.units.m**3 / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "v(comp, internal energy, pressure)",
    },
    # Dimensionless Helmholtz energy to calculate other thermo properties
    "phi0_func": {  # ideal part
        "fname": "phi0",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
        "doc": "phi0(comp, delta, tau)",
    },
    "phi0_d_func": {  # ideal part derivative wrt delta
        "fname": "phi0_d",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_dd_func": {  # ideal part second derivative wrt delta
        "fname": "phi0_dd",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_t_func": {  # ideal part derivative wrt tau
        "fname": "phi0_t",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_dt_func": {  # ideal part second derivative wrt delta and tau
        "fname": "phi0_dt",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phi0_tt_func": {  # ideal part second derivative wrt tau
        "fname": "phi0_tt",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_func": {  # residual part
        "fname": "phir",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_d_func": {  # residual part derivative wrt delta
        "fname": "phir_d",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_dd_func": {  # residual part second derivative wrt delta
        "fname": "phir_dd",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_t_func": {  # residual part derivative wrt tau
        "fname": "phir_t",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_dt_func": {  # residual part second derivative wrt delta and tau
        "fname": "phir_dt",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    "phir_tt_func": {  # residual part second derivative wrt tau
        "fname": "phir_tt",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless, dimensionless],
    },
    # Phase specific functions of pressure and tau
    "hvpt_func": {  # vapor enthalpy
        "fname": "hvpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
        "doc": "h_v(p, tau)",
        "latex_symbol": "\h_v",
    },
    "hlpt_func": {  # liquid enthalpy
        "fname": "hlpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "svpt_func": {  # vapor entropy
        "fname": "svpt",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "slpt_func": {  # liquid entropy
        "fname": "slpt",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "uvpt_func": {  # vapor internal energy
        "fname": "uvpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "ulpt_func": {  # liquid internal energy
        "fname": "ulpt",
        "units": pyo.units.kJ / pyo.units.kg,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "delta_liq_func": {  # liquid density
        "fname": "delta_liq",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    "delta_vap_func": {  # vapor density
        "fname": "delta_vap",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kPa, dimensionless],
    },
    # state variable change functions
    "tau_func": {  # tau as a function of h, p
        "fname": "tau",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
    },
    "vf_func": {  # vapor fraction as a function of h, p
        "fname": "vf",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kJ / pyo.units.kg, pyo.units.kPa],
        "doc": "x(comp, h, p)",
    },
    "taus_func": {  # tau as a function of s, p
        "fname": "taus",
        "units": dimensionless,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
    },
    "vfs_func": {  # vapor fraction as a function of s, p
        "fname": "vfs",
        "units": dimensionless,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg / pyo.units.K,
            pyo.units.kPa,
        ],
        "doc": "x(comp, s, p)",
    },
    "tauu_func": {  # tau as a function of u, p
        "fname": "tauu",
        "units": dimensionless,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg,
            pyo.units.kPa,
        ],
    },
    "vfu_func": {  # vapor fraction as a function of u, p
        "fname": "vfu",
        "units": dimensionless,
        "arg_units": [
            dimensionless,
            pyo.units.kJ / pyo.units.kg,
            pyo.units.kPa,
        ],
    },
    # saturation curve as a function of tau
    "p_sat_func": {
        "fname": "p_sat",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless, dimensionless],
    },
    "delta_sat_v_func": {
        "fname": "delta_sat_v",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless],
    },
    "delta_sat_l_func": {
        "fname": "delta_sat_l",
        "units": dimensionless,
        "arg_units": [dimensionless, dimensionless],
    },
    # saturation curve as a function of p
    "tau_sat_func": {
        "fname": "tau_sat",
        "units": dimensionless,
        "arg_units": [dimensionless, pyo.units.kPa],
    },
    # Parameters (these functions take no arguments)
    "mw_func": {
        "fname": "mw",
        "units": pyo.units.g / pyo.units.mol,
        "arg_units": [dimensionless],
    },
    "sgc_func": {
        "fname": "sgc",
        "units": pyo.units.kJ / pyo.units.kg / pyo.units.K,
        "arg_units": [dimensionless],
    },
    "t_star_func": {
        "fname": "t_star",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
    },
    "rho_star_func": {
        "fname": "rho_star",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [dimensionless],
    },
    "pc_func": {
        "fname": "pc",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless],
    },
    "tc_func": {
        "fname": "tc",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
    },
    "rhoc_func": {
        "fname": "rhoc",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [dimensionless],
    },
    "pt_func": {
        "fname": "pt",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless],
    },
    "tt_func": {
        "fname": "tt",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
    },
    "rhot_l_func": {
        "fname": "rhot_l",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [dimensionless],
    },
    "rhot_v_func": {
        "fname": "rhot_v",
        "units": pyo.units.kg / pyo.units.m**3,
        "arg_units": [dimensionless],
    },
    "pmin_func": {
        "fname": "pmin",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless],
    },
    "tmin_func": {
        "fname": "tmin",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
    },
    "pmax_func": {
        "fname": "pmax",
        "units": pyo.units.kPa,
        "arg_units": [dimensionless],
    },
    "tmax_func": {
        "fname": "tmax",
        "units": pyo.units.K,
        "arg_units": [dimensionless],
    },
}


def add_helmholtz_external_functions(blk, names=None):
    """Add a Helmholtz EoS function to a Pyomo Block.

    Args:
        blk: block to add function to
        names: if None, add all functions, if a string add the single function,
            otherwise a list of function names.

    Returns:
        None
    """
    if names is None:
        names = _external_function_map.keys()
    if isinstance(names, str):
        names = [names]
    for name in names:
        if hasattr(blk, name):
            continue
        fdict = _external_function_map[name]
        setattr(
            blk,
            name,
            pyo.ExternalFunction(
                library=_flib,
                function=fdict["fname"],
                units=fdict["units"],
                arg_units=fdict["arg_units"],
                doc=fdict.get("doc", None),
            ),
        )


class HelmholtzThermoExpressions(object):
    """Class to write thermodynamic property expressions.  Take one of these
    possible sets of state variables: {h, p}, {u, p}, {s, p}, {s, T}, {T, x},
    {P, x}, or {T, P, x}, and return an expression for a thermo property.
    This works by converting the given state variables to temperature, density,
    and vapor fraction expressions then using those to write an expression for
    requested property. This writes expressions in a way that looks like a
    thermodynamic property function.
    """

    def __init__(self, blk, parameters, amount_basis=None):
        """Create a new thermodynamic property expression writer class.

        Args:
            blk: the block to attach the external functions to
            parameters: property parameter block

        Returns:
            HelmholtzThermoExpressions
        """
        if not helmholtz_available():
            raise RuntimeError("Helmholtz EoS external functions not available")
        if amount_basis is None:
            amount_basis = parameters.config.amount_basis
        self.param = parameters
        self.blk = blk
        self.amount_basis = amount_basis

    @staticmethod
    def _sv_str(**kwargs):
        a = [x for x in kwargs if kwargs[x] is not None]
        return ", ".join(a)

    def add_funcs(self, names=None):
        add_helmholtz_external_functions(self.blk, names=names)

    def _state_vars(self, **kwargs):
        """
        This function takes the state variable args, and if they are in
        one of the sets (h, p), (s, p) or (u, p), it returns
        (enum for state variable set, block with external function objects,
        the component string, the first state variable with units for
        external function, and the second state variable with units for
        the external functions).  If the state variable set is not one
        of the above, this function returns None for all but the block and
        component.
        """
        c = self.param.pure_component  # string for chemical component
        blk = self.blk  # block with external functions
        if "p" in kwargs and kwargs["p"] is not None:
            p = kwargs["p"] * self.param.uc["Pa to kPa"]
            if "h" in kwargs and kwargs["h"] is not None:
                if self.amount_basis == AmountBasis.MOLE:
                    h = kwargs["h"] * self.param.uc["J/mol to kJ/kg"]
                else:
                    h = kwargs["h"] * self.param.uc["J/kg to kJ/kg"]
                return StateVars.PH, blk, c, h, p
            if "s" in kwargs and kwargs["s"] is not None:
                if self.amount_basis == AmountBasis.MOLE:
                    s = kwargs["s"] * self.param.uc["J/mol/K to kJ/kg/K"]
                else:
                    s = kwargs["s"] * self.param.uc["J/kg/K to kJ/kg/K"]
                return StateVars.PS, blk, c, s, p
            if "u" in kwargs and kwargs["u"] is not None:
                if self.amount_basis == AmountBasis.MOLE:
                    u = kwargs["u"] * self.param.uc["J/mol to kJ/kg"]
                else:
                    u = kwargs["u"] * self.param.uc["J/kg to kJ/kg"]
                return StateVars.PU, blk, c, u, p
        return None, blk, c, None, None

    def basic_calculations(
        self, h=None, s=None, p=None, T=None, u=None, x=None, tau=None
    ):
        """This function is called as the basis for most thermo expression
        writer functions.  It takes the given state variables and returns
        expressions for liquid density, vapor density, vapor fraction and
        temperature, which can be used to write an expression for any thermo
        quantity.
        """
        c = self.param.pure_component

        # 1.) convert units to those expected by external functions
        if self.amount_basis == AmountBasis.MOLE:
            if h is not None:
                h *= self.param.uc["J/mol to kJ/kg"]
            if u is not None:
                u *= self.param.uc["J/mol to kJ/kg"]
            if s is not None:
                s *= self.param.uc["J/mol/K to kJ/kg/K"]
        else:
            if h is not None:
                h *= self.param.uc["J/kg to kJ/kg"]
            if u is not None:
                u *= self.param.uc["J/kg to kJ/kg"]
            if s is not None:
                s *= self.param.uc["J/kg/K to kJ/kg/K"]
        if p is not None:
            p *= self.param.uc["Pa to kPa"]
        if T is not None:
            tau = self.param.temperature_star / T
        if tau is not None:
            T = self.param.temperature_star / tau
        # 2.) find the block with the external functions
        blk = self.blk

        # 3.) Take given state variables and convert to density, T, and x
        if h is not None and p is not None:
            # h, p
            self.add_funcs(names=["tau_func", "vf_func"])
            tau = blk.tau_func(c, h, p, _data_dir)
            if x is None:
                x = blk.vf_func(c, h, p, _data_dir)
        elif s is not None and p is not None:
            # s, p
            self.add_funcs(names=["taus_func", "vfs_func"])
            tau = blk.taus_func(c, s, p, _data_dir)
            if x is None:
                x = blk.vfs_func(c, s, p, _data_dir)
        elif u is not None and p is not None:
            # u, p
            self.add_funcs(names=["tauu_func", "vfu_func"])
            tau = blk.tauu_func(c, u, p, _data_dir)
            if x is None:
                x = blk.vfu_func(c, u, p, _data_dir)
        elif x is not None and T is not None and p is not None:
            # T, P, x (okay, but I hope you know what you're doing)
            pass
        elif x is not None and p is not None:
            # x, p
            self.add_funcs(names=["tau_sat_func"])
            tau = blk.tau_sat_func(c, p, _data_dir)
        elif x is not None and T is not None:
            # x, T
            self.add_funcs(names=["p_sat_func"])
            p = blk.p_sat_func(c, tau, _data_dir)
        else:
            m = "This choice of state variables ({}) is not yet supported.".format(
                self._sv_str(h=h, s=s, p=p, T=T, u=u, x=x)
            )
            _log.error(m)
            raise NotImplementedError(m)

        # 4.) Calculate density
        self.add_funcs(names=["delta_liq_func", "delta_vap_func"])
        delta_liq = blk.delta_liq_func(c, p, tau, _data_dir)
        delta_vap = blk.delta_vap_func(c, p, tau, _data_dir)
        # 5.) From here its straight forward to calculate any property
        return blk, delta_liq, delta_vap, tau, x, c

    def s(self, **kwargs):
        """Mixed phase entropy"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["s_hp_func"])
            s = blk.s_hp_func(c, u, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["s_up_func"])
            s = blk.s_up_func(c, u, p, _data_dir)
        elif sv == StateVars.PS:
            s = u
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            self.add_funcs(names=["s_func"])
            s = (
                blk.s_func(c, delta_liq, tau, _data_dir) * (1 - x)
                + blk.s_func(c, delta_vap, tau, _data_dir) * x
            )
        if self.amount_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def s_liq(self, **kwargs):
        """Liquid phase entropy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["s_func"])
        s = blk.s_func(c, delta_liq, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def s_vap(self, **kwargs):
        """Vapor phase entropy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["s_func"])
        s = blk.s_func(c, delta_vap, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def h(self, **kwargs):
        """Mixed phase enthalpy"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            h = u
        elif sv == StateVars.PU:
            self.add_funcs(names=["h_up_func"])
            h = blk.h_up_func(c, u, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["h_sp_func"])
            h = blk.h_sp_func(c, u, p, _data_dir)
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            self.add_funcs(names=["h_func"])
            h = (
                blk.h_func(c, delta_liq, tau, _data_dir) * (1 - x)
                + blk.h_func(c, delta_vap, tau, _data_dir) * x
            )
        if self.amount_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg to J/mol"]
        return h * self.param.uc["kJ/kg to J/kg"]

    def h_liq(self, **kwargs):
        """Liquid phase enthalpy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["h_func"])
        h = blk.h_func(c, delta_liq, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg to J/mol"]
        return h * self.param.uc["kJ/kg to J/kg"]

    def h_vap(self, **kwargs):
        """Vapor phase enthalpy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["h_func"])
        h = blk.h_func(c, delta_vap, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg to J/mol"]
        return h * self.param.uc["kJ/kg to J/kg"]

    def u(self, **kwargs):
        """Mixed phase internal energy"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["u_hp_func"])
            u = blk.u_hp_func(c, u, p, _data_dir)
        elif sv == StateVars.PU:
            pass
        elif sv == StateVars.PS:
            self.add_funcs(names=["u_sp_func"])
            u = blk.u_sp_func(c, u, p, _data_dir)
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            self.add_funcs(names=["u_func"])
            u = (
                blk.u_func(c, delta_liq, tau, _data_dir) * (1 - x)
                + blk.u_func(c, delta_vap, tau, _data_dir) * x
            )
        if self.amount_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg to J/mol"]
        return u * self.param.uc["kJ/kg to J/kg"]

    def u_liq(self, **kwargs):
        """Liquid phase internal energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["u_func"])
        u = blk.u_func(c, delta_liq, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg to J/mol"]
        return u * self.param.uc["kJ/kg to J/kg"]

    def u_vap(self, **kwargs):
        """Vapor phase internal energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["u_func"])
        u = blk.u_func(c, delta_vap, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg to J/mol"]
        return u * self.param.uc["kJ/kg to J/kg"]

    def g(self, **kwargs):
        """Mixed phase Gibbs free energy"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["g_hp_func"])
            g = blk.g_hp_func(c, u, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["g_up_func"])
            g = blk.g_up_func(c, u, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["g_sp_func"])
            g = blk.g_sp_func(c, u, p, _data_dir)
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            self.add_funcs(names=["g_func"])
            g = (
                blk.g_func(c, delta_liq, tau, _data_dir) * (1 - x)
                + blk.g_func(c, delta_vap, tau, _data_dir) * x
            )
        if self.amount_basis == AmountBasis.MOLE:
            return g * self.param.uc["kJ/kg to J/mol"]
        return g * self.param.uc["kJ/kg to J/kg"]

    def g_liq(self, **kwargs):
        """Liquid phase Gibbs free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["g_func"])
        g = blk.g_func(c, delta_liq, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return g * self.param.uc["kJ/kg to J/mol"]
        return g * self.param.uc["kJ/kg to J/kg"]

    def g_vap(self, **kwargs):
        """Vapor phase Gibbs free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["g_func"])
        g = blk.g_func(c, delta_vap, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return g * self.param.uc["kJ/kg to J/mol"]
        return g * self.param.uc["kJ/kg to J/kg"]

    def f(self, **kwargs):
        """Mixed phase Helmholtz free energy"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["f_hp_func"])
            f = blk.f_hp_func(c, u, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["f_up_func"])
            f = blk.f_up_func(c, u, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["f_sp_func"])
            f = blk.f_sp_func(c, u, p, _data_dir)
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            self.add_funcs(names=["f_func"])
            f = (
                blk.f_func(c, delta_liq, tau, _data_dir) * (1 - x)
                + blk.f_func(c, delta_vap, tau, _data_dir) * x
            )
        if self.amount_basis == AmountBasis.MOLE:
            return f * self.param.uc["kJ/kg to J/mol"]
        return f * self.param.uc["kJ/kg to J/kg"]

    def f_liq(self, **kwargs):
        """Liquid phase Helmholtz free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["f_func"])
        f = blk.f_func(c, delta_liq, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return f * self.param.uc["kJ/kg to J/mol"]
        return f * self.param.uc["kJ/kg to J/kg"]

    def f_vap(self, **kwargs):
        """Vapor phase Helmholtz free energy"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["f_func"])
        f = blk.f_func(c, delta_vap, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return f * self.param.uc["kJ/kg to J/mol"]
        return f * self.param.uc["kJ/kg to J/kg"]

    def p(self, **kwargs):
        """Pressure"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv is not None:
            pass
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            self.add_funcs(names=["p_func"])
            # The following line looks a bit weird, but it is okay.  When in the
            # two-phase region the pressure for both phases is the same
            p = (
                blk.p_func(c, delta_liq, tau, _data_dir) * (1 - x)
                + blk.p_func(c, delta_vap, tau, _data_dir) * x
            )
        return p * self.param.uc["kPa to Pa"]

    def v_mol(self, **kwargs):
        """Mixed phase molar volume"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["v_hp_func"])
            v = blk.v_hp_func(c, u, p, _data_dir) * self.param.mw
        elif sv == StateVars.PU:
            self.add_funcs(names=["v_up_func"])
            v = blk.v_up_func(c, u, p, _data_dir) * self.param.mw
        elif sv == StateVars.PS:
            self.add_funcs(names=["v_sp_func"])
            v = blk.v_sp_func(c, u, p, _data_dir) * self.param.mw
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            v = (
                ((1 - x) / delta_liq + x / delta_vap)
                / self.param.dens_mass_star
                * self.param.mw
            )
        return v

    def v_mol_liq(self, **kwargs):
        """Liquid phase molar volume"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = self.param.mw / delta_liq / self.param.dens_mass_star
        return v

    def v_mol_vap(self, **kwargs):
        """Vapor phase molar volume"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = self.param.mw / delta_vap / self.param.dens_mass_star
        return v

    def v_mass(self, **kwargs):
        """Mixed phase molar volume"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["v_hp_func"])
            v = blk.v_hp_func(c, u, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["v_up_func"])
            v = blk.v_up_func(c, u, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["v_sp_func"])
            v = blk.v_sp_func(c, u, p, _data_dir)
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            v = ((1 - x) / delta_liq + x / delta_vap) / self.param.dens_mass_star
        return v

    def v_mass_liq(self, **kwargs):
        """Liquid phase molar volume"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = 1 / delta_liq / self.param.dens_mass_star
        return v

    def v_mass_vap(self, **kwargs):
        """Vapor phase molar volume"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        v = 1 / delta_vap / self.param.dens_mass_star
        return v

    def x(self, **kwargs):
        """Vapor faction"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return x

    def T(self, **kwargs):
        """Temperature"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["tau_func"])
            tau = blk.tau_func(c, u, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["tauu_func"])
            tau = blk.tauu_func(c, u, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["taus_func"])
            tau = blk.taus_func(c, u, p, _data_dir)
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return self.param.temperature_star / tau

    def tau(self, **kwargs):
        """Critical Temperature (K)/Temperature (K)"""
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["tau_func"])
            tau = blk.tau_func(c, u, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["tauu_func"])
            tau = blk.tauu_func(c, u, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["taus_func"])
            tau = blk.taus_func(c, u, p, _data_dir)
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return tau

    def delta_liq(self, **kwargs):
        """Return liquid phase reduced density (dens/critical dens) expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq

    def rho_liq(self, **kwargs):
        """Return liquid phase mass density expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq * self.param.dens_mass_star

    def rho_mol_liq(self, **kwargs):
        """Return liquid phase molar density expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_liq * self.param.dens_mass_star / self.param.mw

    def delta_vap(self, **kwargs):
        """Return vapor phase reduced density (dens/critical dens) expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap

    def rho_vap(self, **kwargs):
        """Return vapor phase mass density expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap * self.param.dens_mass_star

    def rho_mol_vap(self, **kwargs):
        """Return vapor phase molar density expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        return delta_vap * self.param.dens_mass_star / self.param.mw

    def cv_liq(self, **kwargs):
        """Return liquid phase constant volume heat capacity expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        cv = blk.cv_func(c, delta_liq, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return cv * self.param.uc["kJ/kg/K to J/mol/K"]
        return cv * self.param.uc["kJ/kg/K to J/kg/K"]

    def cv_mol_liq(self, **kwargs):
        """Backward Compatibility; Return liquid phase molar cv expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        cv = blk.cv_func(c, delta_liq, tau, _data_dir)
        return cv * self.param.uc["kJ/kg/K to J/mol/K"]

    def cv_vap(self, **kwargs):
        """Return vapor phase constant volume heat capacity expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        cv = blk.cv_func(c, delta_vap, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return cv * self.param.uc["kJ/kg/K to J/mol/K"]
        return cv * self.param.uc["kJ/kg/K to J/kg/K"]

    def cv_mol_vap(self, **kwargs):
        """Backward Compatibility; Return vapor phase molar cv expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cv_func"])
        cv = blk.cv_func(c, delta_vap, tau, _data_dir)
        return cv * self.param.uc["kJ/kg/K to J/mol/K"]

    def cp_liq(self, **kwargs):
        """Return liquid phase constant volume heat capacity expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cp_func"])
        cp = blk.cp_func(c, delta_liq, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return cp * self.param.uc["kJ/kg/K to J/mol/K"]
        return cp * self.param.uc["kJ/kg/K to J/kg/K"]

    def cp_mol_liq(self, **kwargs):
        """Backward Compatibility; Return liquid phase molar cp expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cp_func"])
        cp = blk.cp_func(c, delta_liq, tau, _data_dir)
        return cp * self.param.uc["kJ/kg/K to J/mol/K"]

    def cp_vap(self, **kwargs):
        """Return vapor phase constant volume heat capacity expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cp_func"])
        cp = blk.cp_func(c, delta_vap, tau, _data_dir)
        if self.amount_basis == AmountBasis.MOLE:
            return cp * self.param.uc["kJ/kg/K to J/mol/K"]
        return cp * self.param.uc["kJ/kg/K to J/kg/K"]

    def cp_mol_vap(self, **kwargs):
        """Backward Compatibility; Return liquid phase molar cp expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["cp_func"])
        cp = blk.cp_func(c, delta_vap, tau, _data_dir)
        return cp * self.param.uc["kJ/kg/K to J/mol/K"]

    def w(self, **kwargs):
        """Return speed of sound expression, this may not make sense
        in the two phase region
        """
        sv, blk, c, u, p = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["w_hp_func"])
            return blk.w_hp_func(c, u, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["w_up_func"])
            return blk.w_up_func(c, u, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["w_sp_func"])
            return blk.w_sp_func(c, u, p, _data_dir)
        else:
            blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
            self.add_funcs(names=["w_func"])
            return (
                blk.w_func(c, delta_liq, tau, _data_dir) * (1 - x)
                + blk.w_func(c, delta_vap, tau, _data_dir) * x
            )

    def w_liq(self, **kwargs):
        """Return liquid phase speed of sound expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["w_func"])
        return blk.w_func(c, delta_liq, tau, _data_dir)

    def w_vap(self, **kwargs):
        """Return vapor phase speed of sound expression"""
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        self.add_funcs(names=["w_func"])
        return blk.w_func(c, delta_vap, tau, _data_dir)

    def viscosity_liq(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        tmod = get_transport_module(c)
        if tmod is None:
            raise RuntimeError(f"Transport properties not available for {c}")
        return tmod._viscosity(self.param, delta_liq, tau, blk)

    def viscosity_vap(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        tmod = get_transport_module(c)
        if tmod is None:
            raise RuntimeError(f"Transport properties not available for {c}")
        return tmod._viscosity(self.param, delta_vap, tau, blk)

    def thermal_conductivity_liq(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        tmod = get_transport_module(c)
        if tmod is None:
            raise RuntimeError(f"Transport properties not available for {c}")
        return tmod._thermal_conductivity(self.param, delta_liq, tau, blk)

    def thermal_conductivity_vap(self, **kwargs):
        blk, delta_liq, delta_vap, tau, x, c = self.basic_calculations(**kwargs)
        tmod = get_transport_module(c)
        if tmod is None:
            raise RuntimeError(f"Transport properties not available for {c}")
        return tmod._thermal_conductivity(self.param, delta_vap, tau, blk)

    def p_sat(self, T=None, tau=None):
        """Return saturation pressure as a function of T or tau"""
        if T is not None:
            tau = self.param.temperature_star / T
        elif tau is not None:
            pass
        else:
            raise RuntimeError("p_sat expression requires either T or tau arg")
        self.add_funcs(names=["p_sat_func"])
        return (
            self.blk.p_sat_func(self.param.pure_component, tau, _data_dir)
            * self.param.uc["kPa to Pa"]
        )

    def delta_liq_sat(self, T=None, tau=None):
        """Return saturation pressure as a function of T or tau"""
        if T is not None:
            tau = self.param.temperature_star / T
        elif tau is not None:
            pass
        else:
            raise RuntimeError("delta_liq_sat expression requires either T or tau arg")
        self.add_funcs(names=["delta_sat_l_func"])
        return self.blk.delta_sat_l_func(self.param.pure_component, tau, _data_dir)

    def delta_vap_sat(self, T=None, tau=None):
        """Return saturation pressure as a function of T or tau"""
        if T is not None:
            tau = self.param.temperature_star / T
        elif tau is not None:
            pass
        else:
            raise RuntimeError("delta_vap_sat expression requires either T or tau arg")
        self.add_funcs(names=["delta_sat_v_func"])
        return self.blk.delta_sat_v_func(self.param.pure_component, tau, _data_dir)

    def T_sat(self, p):
        """Return saturation temperature as a function of p"""
        p *= self.param.uc["Pa to kPa"]
        self.add_funcs(names=["tau_sat_func"])
        return self.param.temperature_star / self.blk.tau_sat_func(
            self.param.pure_component, p, _data_dir
        )

    def tau_sat(self, p):
        """Return saturation tau as a function of p"""
        p *= self.param.uc["Pa to kPa"]
        self.add_funcs(names=["tau_sat_func"])
        return self.blk.tau_sat_func(self.param.pure_component, p, _data_dir)


@declare_process_block_class("HelmholtzParameterBlock")
class HelmholtzParameterBlockData(PhysicalParameterBlock):
    """
    This is a base class for Helmholtz equations of state using IDAES standard
    Helmholtz EOS external functions written in C++.
    """

    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "pure_component",
        ConfigValue(
            default=None,
            domain=str,
            description="Pure chemical component",
            doc="Pure component to calculate properties for",
        ),
    )
    CONFIG.declare(
        "phase_presentation",
        ConfigValue(
            default=PhaseType.MIX,
            domain=In(PhaseType),
            description="Set the way phases are presented to models",
            doc="""Set the way phases are presented to models. The MIX option
appears to the framework to be a mixed phase containing liquid and/or vapor.
The mixed option can simplify calculations at the unit model level since it can
be treated as a single phase, but unit models such as flash vessels will not be
able to treat the phases independently. The LG option presents as two separate
phases to the framework. The L or G options can be used if it is known for sure
that only one phase is present.
**default** - PhaseType.MIX
**Valid values:** {
**PhaseType.MIX** - Present a mixed phase with liquid and/or vapor,
**PhaseType.LG** - Present a liquid and vapor phase,
**PhaseType.L** - Assume only liquid can be present,
**PhaseType.G** - Assume only vapor can be present}""",
        ),
    )

    CONFIG.declare(
        "state_vars",
        ConfigValue(
            default=StateVars.PH,
            domain=In(StateVars),
            description="State variable set",
            doc="""The set of state variables to use. Depending on the use, one
state variable set or another may be better computationally. Usually pressure
and enthalpy are the best choice because they are well behaved during a phase
change.
**default** - StateVars.PH
**Valid values:** {
**StateVars.PH** - Pressure-Enthalpy,
**StateVars.TPX** - Temperature-Pressure-Quality}""",
        ),
    )

    CONFIG.declare(
        "amount_basis",
        ConfigValue(
            default=AmountBasis.MOLE,
            domain=In(AmountBasis),
            description="State variable set",
            doc="""The set of state variables to use. Depending on the use, one
state variable set or another may be better computationally. Usually pressure
and enthalpy are the best choice because they are well behaved during a phase
change.
**default** - StateVars.PH
**Valid values:** {
**StateVars.PH** - Pressure-Enthalpy,
**StateVars.TPX** - Temperature-Pressure-Quality}""",
        ),
    )

    def available(self):
        """Returns True if the shared library is installed and loads properly
        otherwise returns False
        """
        return helmholtz_available()

    def _suh_tpx(
        self,
        T=None,
        p=None,
        x=None,
        units=None,
        amount_basis=None,
        with_units=False,
        prop="h",
    ):
        """
        Convenience function to calculate enthalpy from temperature and either
        pressure or vapor fraction. This function can be used for inlet streams and
        initialization where temperature is known instead of enthalpy.
        User must provide values for one of these sets of values: {T, P}, {T, x},
        or {P, x}.
        Args:
            T: Temperature
            P: Pressure, None if saturated
            x: Vapor fraction [mol vapor/mol total] (between 0 and 1), None if
                superheated or sub-cooled
            units: The units to report the result in, if None use the default
                units appropriate for the amount basis.
            amount_basis (AmountBasis): Whether to use a mass or mole basis
            with_units (bool): if True return an expression with units
        Returns:
            Total molar enthalpy.
        """
        if amount_basis is None:
            amount_basis = self.config.amount_basis
        if units is None:
            if prop in ["h", "u"]:
                if amount_basis == AmountBasis.MOLE:
                    units = pyo.units.J / pyo.units.mol
                else:
                    units = pyo.units.J / pyo.units.kg
            else:
                if amount_basis == AmountBasis.MOLE:
                    units = pyo.units.J / pyo.units.mol / pyo.units.K
                else:
                    units = pyo.units.J / pyo.units.kg / pyo.units.K

        te = HelmholtzThermoExpressions(self, self, amount_basis=amount_basis)
        tmin = pyo.value(self.temperature_min)
        tmax = pyo.value(self.temperature_max)
        pmin = pyo.value(self.pressure_min)
        pmax = pyo.value(self.pressure_max)

        if not sum((p is None, T is None, x is None)) == 1:
            raise RuntimeError(
                "htpx function must be provided exactly two of the arguments T, p, x"
            )
        if T is not None:
            T = pyo.units.convert(T, to_units=pyo.units.K)
            if not tmin <= pyo.value(T) <= tmax:
                raise RuntimeError(f"T = {pyo.value(T)}, ({tmin} K <= T <= {tmax} K)")
        if x is not None:
            if not 0 <= pyo.value(x) <= 1:
                raise RuntimeError(f"x = {pyo.value(x)}, (0 K <= x <= 1)")
        if p is not None:
            p = pyo.units.convert(p, to_units=pyo.units.Pa)
            if not pmin <= pyo.value(p) <= pmax:
                raise RuntimeError(
                    f"p = {pyo.value(p)}, ({pmin} kPa <= p <= {pmax} kPa)"
                )
            if T is not None:
                # P, T may be under-specified, but assume you know it's clearly a
                # vapor or liquid so figure out which and set x.
                psat = te.p_sat(T)
                if pyo.value(p) < pyo.value(psat):
                    x = 1
                else:
                    x = 0
        if prop == "h":
            if with_units:
                return pyo.value(pyo.units.convert(te.h(T=T, p=p, x=x), units)) * units
            return pyo.value(pyo.units.convert(te.h(T=T, p=p, x=x), units))
        elif prop == "s":
            if with_units:
                return pyo.value(pyo.units.convert(te.s(T=T, p=p, x=x), units)) * units
            return pyo.value(pyo.units.convert(te.s(T=T, p=p, x=x), units))
        elif prop == "u":
            if with_units:
                return pyo.value(pyo.units.convert(te.u(T=T, p=p, x=x), units)) * units
            return pyo.value(pyo.units.convert(te.u(T=T, p=p, x=x), units))

    def htpx(
        self,
        T=None,
        p=None,
        x=None,
        units=None,
        amount_basis=None,
        with_units=False,
    ):
        return self._suh_tpx(
            T=T,
            p=p,
            x=x,
            units=units,
            amount_basis=amount_basis,
            with_units=with_units,
            prop="h",
        )

    def stpx(
        self,
        T=None,
        p=None,
        x=None,
        units=None,
        amount_basis=None,
        with_units=False,
    ):
        return self._suh_tpx(
            T=T,
            p=p,
            x=x,
            units=units,
            amount_basis=amount_basis,
            with_units=with_units,
            prop="s",
        )

    def utpx(
        self,
        T=None,
        p=None,
        x=None,
        units=None,
        amount_basis=None,
        with_units=False,
    ):
        return self._suh_tpx(
            T=T,
            p=p,
            x=x,
            units=units,
            amount_basis=amount_basis,
            with_units=with_units,
            prop="u",
        )

    def _set_default_scaling(self):
        """Set default scaling parameters to be used if not otherwise set"""

        # TODO<jce> leaving flow to not break things, but plan to remove it and
        #     make the user specify with no default.
        self.set_default_scaling("flow_mol", 1e-4)
        self.set_default_scaling("flow_mass", 1)
        self.set_default_scaling("flow_mol_comp", 1e-4)
        self.set_default_scaling("flow_mass_comp", 1)
        self.set_default_scaling("flow_vol", 100)
        self.set_default_scaling("flow_mass", 1)

        # Set some scaling with reasonable a priori values
        self.set_default_scaling("temperature_crit", 1e-2)
        self.set_default_scaling("temperature_star", 1e-2)
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("enth_mass", 1e-3)
        self.set_default_scaling("temperature", 1e-1)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("vapor_frac", 1e1)
        self.set_default_scaling("dens_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e1)  # Not liq
        self.set_default_scaling("pressure_crit", 1e-6)
        self.set_default_scaling("dens_mass_crit", 1e-2)
        self.set_default_scaling("mw", 1e3)
        self.set_default_scaling("temperature_sat", 1e-2)
        self.set_default_scaling("temperature_red", 1)
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("dens_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e1)
        self.set_default_scaling("dens_phase_red", 1)
        self.set_default_scaling("enth_mol_sat_phase", 1e-2, index="Liq")
        self.set_default_scaling("enth_mol_sat_phase", 1e-4, index="Vap")
        self.set_default_scaling("enth_mass_sat_phase", 1e-2, index="Liq")
        self.set_default_scaling("enth_mass_sat_phase", 1e-4, index="Vap")
        self.set_default_scaling("dh_vap_mol", 1e-4)
        self.set_default_scaling("dh_vap_mass", 1e-4)
        self.set_default_scaling("energy_internal_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("energy_internal_mol_phase", 1e-4)
        self.set_default_scaling("energy_internal_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("energy_internal_mass_phase", 1e-4)
        self.set_default_scaling("enth_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("enth_mol_phase", 1e-4, index="Vap")
        self.set_default_scaling("entr_mol_phase", 1e-1, index="Liq")
        self.set_default_scaling("entr_mol_phase", 1e-1, index="Vap")
        self.set_default_scaling("enth_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("enth_mass_phase", 1e-4, index="Vap")
        self.set_default_scaling("entr_mass_phase", 1e-1, index="Liq")
        self.set_default_scaling("entr_mass_phase", 1e-1, index="Vap")
        self.set_default_scaling("cp_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("cp_mol_phase", 1e-2, index="Vap")
        self.set_default_scaling("cv_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("cv_mol_phase", 1e-2, index="Vap")
        self.set_default_scaling("cp_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("cp_mass_phase", 1e-2, index="Vap")
        self.set_default_scaling("cv_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("cv_mass_phase", 1e-2, index="Vap")
        self.set_default_scaling("dens_mol_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mol_phase", 1e-4, index="Vap")
        self.set_default_scaling("phase_frac", 10)
        self.set_default_scaling("energy_internal_mol", 1e-3)
        self.set_default_scaling("energy_internal_mass", 1e-3)
        self.set_default_scaling("entr_mol", 1e-1)
        self.set_default_scaling("entr_mass", 1e-1)
        self.set_default_scaling("cp_mol", 1e-2)
        self.set_default_scaling("cv_mol", 1e-2)
        self.set_default_scaling("cp_mass", 1e-2)
        self.set_default_scaling("cv_mass", 1e-2)
        self.set_default_scaling("dens_mass", 1)
        self.set_default_scaling("dens_mol", 1e-3)
        self.set_default_scaling("heat_capacity_ratio", 1e1)

    def _create_component_and_phase_objects(self):
        # Create chemical component objects
        pp = self.config.phase_presentation
        for c in self.component_list:
            setattr(self, str(c), Component(_component_list_exists=True))
        # Create phase objects
        if pp == PhaseType.MIX:
            self.private_phase_list = pyo.Set(initialize=["Vap", "Liq"])
            self.Mix = Phase()
        elif pp == PhaseType.LG:
            self.private_phase_list = pyo.Set(initialize=["Vap", "Liq"])
            self.Liq = LiquidPhase()
            self.Vap = VaporPhase()
        elif pp == PhaseType.L:
            self.private_phase_list = pyo.Set(initialize=["Liq"])
            self.Liq = LiquidPhase()
        elif pp == PhaseType.G:
            self.private_phase_list = pyo.Set(initialize=["Vap"])
            self.Vap = VaporPhase()

    def build(self):
        if not self.available():
            raise RuntimeError("Helmholtz EoS external functions not available")
        super().build()
        # Check if the specified component is supported
        if not component_registered(self.config.pure_component):
            raise ConfigurationError(
                f"Component {self.config.pure_component} not supported."
            )
        # This is imported here to avoid a circular import
        from idaes.models.properties.general_helmholtz.helmholtz_state import (
            HelmholtzStateBlock,
        )

        self._state_block_class = HelmholtzStateBlock
        # set the component_list as required for the generic IDAES properties
        self.component_list = pyo.Set(initialize=[self.config.pure_component])
        # since this a only a pure component package, have a specific
        # pure_component attribute
        self.pure_component = self.config.pure_component
        # State var set
        self.state_vars = self.config.state_vars
        # Phase equilibrium description
        self.phase_equilibrium_idx = (pyo.Set(initialize=[1]),)
        self.phase_equilibrium_list = ({1: ["H2O", ("Vap", "Liq")]},)
        # Add default scaling factors for property expressions or variables
        self._set_default_scaling()
        # Add idaes component and phase objects
        self._create_component_and_phase_objects()
        # To ensure consistency pull parameters from external functions
        add_helmholtz_external_functions(  # Add parameter external functions
            self,
            [
                # Constants
                "sgc_func",  # specific gas constant
                "mw_func",  # molecular weight
                # Critical properties
                "pc_func",  # Critical pressure
                "tc_func",  # Critical temperature
                "t_star_func",  # Critical temperature
                "rhoc_func",  # Critical density
                "rho_star_func",  # Critical temperature
                # triple point properties
                "pt_func",  # triple point pressure
                "tt_func",  # triple point temperature
                "rhot_l_func",  # triple point liquid density
                "rhot_v_func",  # triple point vapor density
                # triple point properties
                "pt_func",  # triple point pressure
                "tt_func",  # triple point temperature
                "rhot_l_func",  # triple point liquid density
                "rhot_v_func",  # triple point vapor density
                # Bounds
                "pmin_func",  # pmin
                "tmin_func",  # tmin
                "pmax_func",  # pmax
                "tmax_func",  # tmax
                #
                "hlpt_func",
                "slpt_func",
                "ulpt_func",
                "hvpt_func",
                "svpt_func",
                "uvpt_func",
            ],
        )
        # The parameters are constants and we don't want to call the external
        # functions more than once, so define Pyomo parameters with the values,
        # use the external function here to avoid defining the parameters twice
        pu = pyo.units
        cmp = self.pure_component
        self.add_param(
            "mw",
            pu.convert(self.mw_func(cmp, _data_dir), pu.kg / pu.mol),
        )
        self.add_param(
            "sgc",
            pyo.units.convert(self.sgc_func(cmp, _data_dir), pu.J / pu.kg / pu.K),
        )
        self.add_param(
            "sgc_mol",
            pu.convert(self.sgc * self.mw, pu.J / pu.mol / pu.K),
        )
        self.add_param(
            "pressure_crit",
            pu.convert(self.pc_func(cmp, _data_dir), pu.Pa),
        )
        self.add_param(
            "pressure_trip",
            pu.convert(self.pt_func(cmp, _data_dir), pu.Pa),
        )
        self.add_param(
            "pressure_min",
            pu.convert(self.pmin_func(cmp, _data_dir), pu.Pa),
        )
        self.add_param(
            "pressure_max",
            pu.convert(self.pmax_func(cmp, _data_dir), pu.Pa),
        )
        self.add_param(
            "default_pressure_value",
            pu.convert((self.pressure_crit + self.pressure_trip) / 2.0, pu.Pa),
        )
        self.default_pressure_bounds = (self.pressure_min, self.pressure_max)
        self.add_param(
            "temperature_star",
            pu.convert(self.t_star_func(cmp, _data_dir), pu.K),
        )
        self.add_param(
            "temperature_crit",
            pu.convert(self.tc_func(cmp, _data_dir), pu.K),
        )
        self.add_param(
            "temperature_trip",
            pu.convert(self.tt_func(cmp, _data_dir), pu.K),
        )
        self.add_param(
            "temperature_min",
            pu.convert(self.tmin_func(cmp, _data_dir), pu.K),
        )
        self.add_param(
            "temperature_max",
            pu.convert(self.tmax_func(cmp, _data_dir), pu.K),
        )
        self.add_param(
            "default_temperature_value",
            pu.convert((self.temperature_crit + self.temperature_trip) / 2.0, pu.K),
        )
        self.default_temperature_bounds = (self.temperature_min, self.temperature_max)
        self.add_param(
            "dens_mass_star",
            pu.convert(self.rho_star_func(cmp, _data_dir), pu.kg / pu.m**3),
        )
        self.add_param(
            "dens_mol_star",
            pu.convert(self.dens_mass_star / self.mw, pu.mol / pu.m**3),
        )
        self.add_param(
            "dens_mass_crit",
            pu.convert(self.rhoc_func(cmp, _data_dir), pu.kg / pu.m**3),
        )
        self.add_param(
            "dens_mol_crit",
            pu.convert(self.dens_mass_crit / self.mw, pu.mol / pu.m**3),
        )

        self.uc = {
            "J/mol to kJ/kg": (pyo.units.kJ / 1000 / pyo.units.J) / self.mw,
            "J/mol to J/kg": 1.0 / self.mw,
            "kJ/kg to J/mol": (pyo.units.J * 1000 / pyo.units.kJ) * self.mw,
            "J/mol/K to kJ/kg/K": (pyo.units.kJ / 1000 / pyo.units.J) / self.mw,
            "J/mol/K to J/kg/K": 1 / self.mw,
            "kJ/kg/K to J/mol/K": (pyo.units.J * 1000 / pyo.units.kJ) * self.mw,
            "J/kg to kJ/kg": (pyo.units.kJ / 1000 / pyo.units.J),
            "J/kg to J/mol": self.mw,
            "J/kg/K to J/mol/K": self.mw,
            "kJ/kg to J/kg": (pyo.units.J * 1000 / pyo.units.kJ),
            "J/kg/K to kJ/kg/K": (pyo.units.kJ / 1000 / pyo.units.J),
            "kJ/kg/K to J/kg/K": (pyo.units.J * 1000 / pyo.units.kJ),
            "kPa to Pa": (pyo.units.Pa * 1000 / pyo.units.kPa),
            "Pa to kPa": (pyo.units.kPa / 1000 / pyo.units.Pa),
            "kg/m3 to mol/m3": 1 / self.mw,
        }
        self.add_param(
            "enthalpy_mol_min",
            self.hlpt_func(
                cmp,
                self.pressure_trip * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_trip,
                _data_dir,
            )
            * self.uc["kJ/kg to J/mol"],
        )
        self.add_param(
            "enthalpy_mass_min",
            self.hlpt_func(
                cmp,
                self.pressure_trip * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_trip,
                _data_dir,
            )
            * self.uc["kJ/kg to J/kg"],
        )
        self.add_param(
            "enthalpy_mol_max",
            self.hlpt_func(
                cmp,
                self.pressure_max * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_max,
                _data_dir,
            )
            * self.uc["kJ/kg to J/mol"],
        )
        self.add_param(
            "enthalpy_mass_max",
            self.hlpt_func(
                cmp,
                self.pressure_max * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_max,
                _data_dir,
            )
            * self.uc["kJ/kg to J/kg"],
        )
        self.add_param("default_enthalpy_mass_value", self.enthalpy_mass_min)
        self.add_param("default_enthalpy_mol_value", self.enthalpy_mol_min)
        self.default_enthalpy_mass_bounds = (
            self.enthalpy_mass_min,
            self.enthalpy_mass_max,
        )
        self.default_enthalpy_mol_bounds = (
            self.enthalpy_mol_min,
            self.enthalpy_mol_max,
        )
        self.add_param(
            "entropy_mol_min",
            self.slpt_func(
                cmp,
                self.pressure_trip * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_trip,
                _data_dir,
            )
            * self.uc["kJ/kg/K to J/mol/K"],
        )
        self.add_param(
            "entropy_mass_min",
            self.slpt_func(
                cmp,
                self.pressure_trip * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_trip,
                _data_dir,
            )
            * self.uc["kJ/kg/K to J/kg/K"],
        )
        self.add_param(
            "entropy_mol_max",
            self.slpt_func(
                cmp,
                self.pressure_max * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_max,
                _data_dir,
            )
            * self.uc["kJ/kg/K to J/mol/K"],
        )
        self.add_param(
            "entropy_mass_max",
            self.slpt_func(
                cmp,
                self.pressure_max * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_max,
                _data_dir,
            )
            * self.uc["kJ/kg/K to J/kg/K"],
        )
        self.add_param("default_entropy_mass_value", self.entropy_mass_min)
        self.add_param("default_entropy_mol_value", self.entropy_mol_min)
        self.default_entropy_mass_bounds = (
            self.entropy_mass_min,
            self.entropy_mass_max,
        )
        self.default_entropy_mol_bounds = (
            self.entropy_mol_min,
            self.entropy_mol_max,
        )
        self.add_param(
            "energy_internal_mol_min",
            self.ulpt_func(
                cmp,
                self.pressure_trip * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_trip,
                _data_dir,
            )
            * self.uc["kJ/kg to J/mol"],
        )
        self.add_param(
            "energy_internal_mass_min",
            self.ulpt_func(
                cmp,
                self.pressure_trip * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_trip,
                _data_dir,
            )
            * self.uc["kJ/kg to J/kg"],
        )
        self.add_param(
            "energy_internal_mol_max",
            self.ulpt_func(
                cmp,
                self.pressure_max * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_max,
                _data_dir,
            )
            * self.uc["kJ/kg to J/mol"],
        )
        self.add_param(
            "energy_internal_mass_max",
            self.ulpt_func(
                cmp,
                self.pressure_max * self.uc["Pa to kPa"],
                self.temperature_star / self.temperature_max,
                _data_dir,
            )
            * self.uc["kJ/kg to J/kg"],
        )
        self.add_param(
            "default_energy_internal_mass_value", self.energy_internal_mass_min
        )
        self.add_param(
            "default_energy_internal_mol_value", self.energy_internal_mol_min
        )
        self.default_energy_internal_mass_bounds = (
            self.energy_internal_mass_min,
            self.energy_internal_mass_max,
        )
        self.default_energy_internal_mol_bounds = (
            self.energy_internal_mol_min,
            self.energy_internal_mol_max,
        )

        # Smoothing arameters for TPX complimentarity form
        if self.config.state_vars == StateVars.TPX:
            self.smoothing_pressure_over = pyo.Param(
                mutable=True,
                initialize=1e-4,
                doc="Smooth max parameter (pressure over)",
                units=pyo.units.kPa,
            )

            self.smoothing_pressure_under = pyo.Param(
                mutable=True,
                initialize=1e-4,
                doc="Smooth max parameter (pressure under)",
                units=pyo.units.kPa,
            )

    def add_param(self, name, expr):
        self.add_component(
            name,
            pyo.Param(
                initialize=pyo.value(expr),
                units=pyo.units.get_units(expr),
                domain=pyo.Reals,
            ),
        )

    def initialize(self, *args, **kwargs):
        pass

    def dome_data(
        self,
        amount_basis=None,
        pressure_unit=pyo.units.kPa,
        energy_unit=pyo.units.kJ,
        mass_unit=pyo.units.kg,
        mol_unit=pyo.units.kmol,
        n=60,
    ):
        """Get data to plot the two-phase dome or saturation curve.  This data
        can be used to plot the 2 phase dome for p-h and t-s diagrams and the
        saturation curve on the p-t diagram.

        Args:
            amount_bases (AmountBasis): Mass or mole basis. Get from parameter
                block if None.
            pressure_unit (PyomoUnit): Pressure units of measure
            energy_unit (PyomoUnit): Energy units of measure
            mass_unit (PyomoUnit): Mass units of measure
            mol_unit (PyomoUnit): Mole unit of measure

        Returns:
            dict: dictonary with the keys {'T', 'tau', 'p', 'delta_liq',
                'delta_vap', 'h_liq', 'h_vap', 's_liq', 's_vap'} each a list of
                numbers corresponding to states along the two-phase dome.
        """

        # Attach the external functions needed to this parameter block (self) as
        # Pyomo ExternalFunction.
        add_helmholtz_external_functions(
            self,
            [
                "p_sat_func",
                "delta_sat_l_func",
                "delta_sat_v_func",
                "h_func",
                "s_func",
            ],
        )
        # Set the amount basis, if not provided by the user get it from this
        # paraemeter block (self)
        if amount_basis is None:
            amount_basis = self.config.amount_basis

        # We'll do the calcuation based on temeprautre.  We want more points
        # near the critical point, so use log space to get that. The lower bound
        # for temperature is the triple point and the upper bound is the
        # critical point.  Here temperature is as tau = T*/T.
        tau_c = pyo.value(self.temperature_star / self.temperature_crit)
        tau_t = pyo.value(self.temperature_star / self.temperature_trip)
        tau_dist_vec = np.logspace(-5, 0, n)
        tau_vec = [tau_c + td * (tau_t - tau_c) for td in tau_dist_vec]
        tau_vec = [tau_c] + tau_vec

        # Get pressures by calling the p_sat function from the temperature vector
        p_vec = [
            pyo.value(
                pyo.units.convert(
                    self.p_sat_func(self.pure_component, tau, _data_dir), pressure_unit
                )
            )
            for tau in tau_vec
        ]

        # Get the reduced density vector for the liquid, from this and
        # tau we can clculate all the rest of the properties for sat liquid
        delta_sat_l_vec = [
            pyo.value(self.delta_sat_l_func(self.pure_component, tau, _data_dir))
            for tau in tau_vec
        ]
        # Get the reduced density vector for the liquid, from this and
        # tau we can clculate all the rest of the properties for sat liquid
        delta_sat_v_vec = [
            pyo.value(self.delta_sat_v_func(self.pure_component, tau, _data_dir))
            for tau in tau_vec
        ]
        if amount_basis == AmountBasis.MOLE:
            s_liq_vec = [
                pyo.value(
                    pyo.units.convert(
                        self.s_func(self.pure_component, delta, tau, _data_dir)
                        * self.uc["kJ/kg/K to J/mol/K"],
                        energy_unit / mol_unit / pyo.units.K,
                    )
                )
                for delta, tau in zip(delta_sat_l_vec, tau_vec)
            ]
            s_vap_vec = [
                pyo.value(
                    pyo.units.convert(
                        self.s_func(self.pure_component, delta, tau, _data_dir)
                        * self.uc["kJ/kg/K to J/mol/K"],
                        energy_unit / mol_unit / pyo.units.K,
                    )
                )
                for delta, tau in zip(delta_sat_v_vec, tau_vec)
            ]
            h_liq_vec = [
                pyo.value(
                    pyo.units.convert(
                        self.h_func(self.pure_component, delta, tau, _data_dir)
                        * self.uc["kJ/kg to J/mol"],
                        energy_unit / mol_unit,
                    )
                )
                for delta, tau in zip(delta_sat_l_vec, tau_vec)
            ]
            h_vap_vec = [
                pyo.value(
                    pyo.units.convert(
                        self.h_func(self.pure_component, delta, tau, _data_dir)
                        * self.uc["kJ/kg to J/mol"],
                        energy_unit / mol_unit,
                    )
                )
                for delta, tau in zip(delta_sat_v_vec, tau_vec)
            ]
        else:
            s_liq_vec = [
                pyo.value(
                    pyo.units.convert(
                        self.s_func(self.pure_component, delta, tau, _data_dir),
                        energy_unit / mass_unit / pyo.units.K,
                    )
                )
                for delta, tau in zip(delta_sat_l_vec, tau_vec)
            ]
            s_vap_vec = [
                pyo.value(
                    pyo.units.convert(
                        self.s_func(self.pure_component, delta, tau, _data_dir),
                        energy_unit / mass_unit / pyo.units.K,
                    )
                )
                for delta, tau in zip(delta_sat_v_vec, tau_vec)
            ]
            h_liq_vec = [
                pyo.value(
                    pyo.units.convert(
                        self.h_func(self.pure_component, delta, tau, _data_dir),
                        energy_unit / mass_unit,
                    )
                )
                for delta, tau in zip(delta_sat_l_vec, tau_vec)
            ]
            h_vap_vec = [
                pyo.value(
                    pyo.units.convert(
                        self.h_func(self.pure_component, delta, tau, _data_dir),
                        energy_unit / mass_unit,
                    )
                )
                for delta, tau in zip(delta_sat_v_vec, tau_vec)
            ]

        return {
            "T": [pyo.value(self.temperature_star / tau) for tau in tau_vec],
            "tau": tau_vec,
            "p": p_vec,
            "delta_liq": delta_sat_l_vec,
            "delta_vap": delta_sat_v_vec,
            "h_liq": h_liq_vec,
            "h_vap": h_vap_vec,
            "s_liq": s_liq_vec,
            "s_vap": s_vap_vec,
        }

    def isotherms(self, temperatures):
        """Get isotherm data for a P-H diagram.

        Args:
            temperatures: A list of temperatures

        Returns:
            dict: The keys are temperatures the values are dicts with "p", "h",
                "s", and "delta" data for the isotherm.
        """
        add_helmholtz_external_functions(
            self,
            [
                "p_sat_func",
                "delta_liq_func",
                "delta_vap_func",
                "h_func",
                "s_func",
            ],
        )
        pt = pyo.value(pyo.units.convert(self.pressure_trip, pyo.units.kPa))
        pc = pyo.value(pyo.units.convert(self.pressure_crit, pyo.units.kPa))
        pmax = pyo.value(pyo.units.convert(self.pressure_max, pyo.units.kPa))
        d = {}

        def _pvec(d, tau, p1, p2=None, phase="sat"):
            if phase == "sat":
                p_vec = [p1, p1]
            elif phase == "liq":
                dist_vec = np.logspace(-4, -0.25, 20)
                vec = [p1 + pd * (p2 - p1) for pd in dist_vec]
                p_vec = [p1] + vec
            else:
                dist_vec = np.logspace(-4, -0.3, 20)
                vec = [p1 + pd * (p2 - p1) for pd in dist_vec]
                p_vec = [p1] + vec + np.linspace(vec[-1], p2, 10).tolist()
            if phase == "liq" or phase == "sc":
                delta = [
                    pyo.value(
                        self.delta_liq_func(self.pure_component, p, tau, _data_dir)
                    )
                    for p in p_vec
                ]
            elif phase == "vap":
                delta = [
                    pyo.value(
                        self.delta_vap_func(self.pure_component, p, tau, _data_dir)
                    )
                    for p in p_vec
                ]
            else:  # sat
                delta = [
                    pyo.value(
                        self.delta_liq_func(
                            self.pure_component, p_vec[0], tau, _data_dir
                        )
                    ),
                    pyo.value(
                        self.delta_vap_func(
                            self.pure_component, p_vec[1], tau, _data_dir
                        )
                    ),
                ]
            h_vec = [
                pyo.value(self.h_func(self.pure_component, dv, tau, _data_dir))
                for dv in delta
            ]
            s_vec = [
                pyo.value(self.s_func(self.pure_component, dv, tau, _data_dir))
                for dv in delta
            ]
            d["p"] = p_vec
            d["delta"] = delta
            d["h"] = h_vec
            d["s"] = s_vec

        for t in temperatures:
            d[t] = {}
            d2 = d[t]
            tau = self.temperature_star / t

            for key in ["liq", "vap", "sat", "sc"]:
                d2[key] = {}

            if t >= pyo.value(self.temperature_crit):
                _pvec(d2["sc"], tau, pc, pmax, "sc")
                _pvec(d2["vap"], tau, pc, pt, "vap")
            else:
                p_sat = pyo.value(self.p_sat_func(self.pure_component, tau, _data_dir))
                _pvec(d2["liq"], tau, p_sat, pmax, "liq")
                _pvec(d2["sat"], tau, p_sat, p_sat, "sat")
                _pvec(d2["vap"], tau, p_sat, pt, "vap")

        return d

    def ph_diagram(
        self,
        ylim=None,
        xlim=None,
        points=None,
        figsize=None,
        dpi=None,
        isotherms=None,
        isotherms_line_format=None,
        isotherms_label=True,
    ):
        if points is None:
            points = {}

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        if ylim is not None:
            ax.set_ylim(ylim)
        if xlim is not None:
            ax.set_xlim(xlim)

        # Get some parameters for plot limits
        pc = pyo.value(pyo.units.convert(self.pressure_crit, pyo.units.kPa))
        pt = pyo.value(pyo.units.convert(self.pressure_trip, pyo.units.kPa))
        tc = pyo.value(self.temperature_crit)
        ts = pyo.value(self.temperature_star)
        tt = pyo.value(self.temperature_trip)
        pmax = pyo.value(pyo.units.convert(self.pressure_max, pyo.units.kPa))

        dome = self.dome_data()

        # plot saturaion curves use log scale for pressure
        ax.set_yscale("log")
        ax.plot(dome["h_liq"], dome["p"], c="b", label="sat liquid")
        ax.plot(dome["h_vap"], dome["p"], c="r", label="sat vapor")

        # Temperatures for isotherms
        t_vec = np.linspace(
            pyo.value(self.temperature_trip), pyo.value(self.temperature_crit), 6
        )
        p = {}
        h_l = {}
        h_v = {}

        if isotherms is not None:  # plot isotherms
            if isinstance(isotherms, np.ndarray):
                t_vec = isotherms.tolist()
            elif isinstance(isotherms, (list, tuple)):
                t_vec = isotherms
            else:
                t_vec = np.linspace(
                    pyo.value(self.temperature_trip),
                    pyo.value(self.temperature_crit),
                    6,
                ).tolist()
            isotherm_data = self.isotherms(temperatures=t_vec)
            if isotherms_line_format is None:
                isotherms_line_format = {"c": "g"}
            for t, dat in isotherm_data.items():
                for phase in ["sc", "liq", "vap", "sat"]:
                    if dat[phase]:
                        ax.plot(
                            dat[phase]["h"], dat[phase]["p"], **isotherms_line_format
                        )
            if isotherms_label:
                for t, dat in isotherm_data.items():
                    if dat["sat"]:
                        ppos = dat["sat"]["p"][0]
                        hpos = (dat["sat"]["h"][0] + dat["sat"]["h"][1]) / 2.0
                        angle = -30
                    else:
                        ppos = dat["vap"]["p"][0]
                        hpos = dat["vap"]["h"][0]
                        angle = -30
                    ax.text(
                        hpos,
                        ppos,
                        f"{t:.2f}  K",
                        fontsize="small",
                        ha="center",
                        va="center",
                        rotation=angle,
                    )
        x = []
        y = []
        for p, v in points.items():
            ax.scatter([v[0]], [v[1]])
            ax.text(v[0], v[1], p, ha="center", fontsize="xx-large")
            x.append(v[0])
            y.append(v[1])
        if len(x) > 1:  # This closes the loop
            x.append(x[0])
            y.append(y[0])
        ax.plot(x, y, c="black")

        # Titles
        ax.set_title(f"P-H Diagram for {self.pure_component}")
        ax.set_xlabel("Enthalpy (kJ/kg)")
        ax.set_ylabel("Pressure (kPa)")
        return fig, ax

    def ts_diagram(self, ylim=None, xlim=None, points=None, figsize=None, dpi=None):
        # Add external functions needed to plot PH-diagram
        if points is None:
            points = {}

        add_helmholtz_external_functions(
            self,
            [
                "p_sat_func",
                "tau_sat_func",
                "delta_sat_l_func",
                "delta_sat_v_func",
                "s_func",
                "slpt_func",
                "svpt_func",
                "delta_liq_func",
                "delta_vap_func",
            ],
        )
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        if ylim is not None:
            ax.set_ylim(ylim)
        if xlim is not None:
            ax.set_xlim(xlim)

        # Get some parameters for plot limits
        pc = pyo.value(pyo.units.convert(self.pressure_crit, pyo.units.kPa))
        pt = pyo.value(pyo.units.convert(self.pressure_trip, pyo.units.kPa))
        tc = pyo.value(self.temperature_crit)
        ts = pyo.value(self.temperature_star)
        tt = pyo.value(self.temperature_trip)
        pmax = pyo.value(pyo.units.convert(self.pressure_max, pyo.units.kPa))

        dome = self.dome_data()

        # plot saturaion curves use log scale for pressure
        # plt.yscale("log")
        ax.plot(dome["s_liq"], dome["T"], c="b", label="sat liquid")
        ax.plot(dome["s_vap"], dome["T"], c="r", label="sat vapor")

        x = []
        y = []
        for p, v in points.items():
            ax.scatter([v[0]], [v[1]])
            ax.text(v[0], v[1], p, ha="center", fontsize="xx-large")
            x.append(v[0])
            y.append(v[1])
        if len(x) > 1:
            x.append(x[0])
            y.append(y[0])
        plt.plot(x, y, c="black")

        ax.set_title(f"T-S Diagram for {self.pure_component}")
        ax.set_xlabel("Entropy (kJ/kg/K)")
        ax.set_ylabel("Temperature (K)")
        return fig, ax

    # TODO: points argument is unused
    def pt_diagram(self, ylim=None, xlim=None, points=None, figsize=None, dpi=None):
        # Add external functions needed to plot PH-diagram
        add_helmholtz_external_functions(
            self,
            [
                "p_sat_func",
                "tau_sat_func",
                "delta_sat_l_func",
                "delta_sat_v_func",
                "s_func",
                "slpt_func",
                "svpt_func",
                "delta_liq_func",
                "delta_vap_func",
            ],
        )
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        if ylim is not None:
            ax.set_ylim(ylim)
        if xlim is not None:
            ax.set_xlim(xlim)

        # Get some parameters for plot limits
        pc = pyo.value(pyo.units.convert(self.pressure_crit, pyo.units.kPa))
        pt = pyo.value(pyo.units.convert(self.pressure_trip, pyo.units.kPa))
        tc = pyo.value(self.temperature_crit)
        ts = pyo.value(self.temperature_star)
        tt = pyo.value(self.temperature_trip)
        pmax = pyo.value(pyo.units.convert(self.pressure_max, pyo.units.kPa))

        # Temperatures to plot the saturation from triple point to critical
        # in the form of tau (Tc/T)
        tau_sat_vec = np.linspace(
            1,
            pyo.value(self.temperature_star) / (pyo.value(self.temperature_trip)),
            200,
        )
        p_sat_vec = [None] * len(tau_sat_vec)
        delta_sat_v_vec = [None] * len(tau_sat_vec)
        delta_sat_l_vec = [None] * len(tau_sat_vec)
        h_sat_v_vec = [None] * len(tau_sat_vec)
        h_sat_l_vec = [None] * len(tau_sat_vec)

        for i, tau in enumerate(tau_sat_vec):
            p_sat_vec[i] = pyo.value(
                self.p_sat_func(self.pure_component, tau, _data_dir)
            )
            delta_sat_l_vec[i] = pyo.value(
                self.delta_sat_l_func(self.pure_component, tau, _data_dir)
            )
            delta_sat_v_vec[i] = pyo.value(
                self.delta_sat_v_func(self.pure_component, tau, _data_dir)
            )

        # plot saturaion curves use log scale for pressure
        ax.set_yscale("log")
        ax.plot(ts / tau_sat_vec, p_sat_vec, c="m", label="sat")
        ax.plot([tc, tc], [pc, pc * 10], c="m", label="sat")
        ax.plot([tc, tc * 1.1], [pc, pc], c="m", label="sat")

        # Points for critical point and triple point
        # plt.scatter([tc], [pc])
        # plt.scatter([tt], [pt])

        ax.set_title(f"P-T Diagram for {self.pure_component}")
        ax.set_ylabel("Pressure (kPa)")
        ax.set_xlabel("Temperature (K)")
        return fig, ax

    # In case you can't rember the axis order in the diagrams
    hp_diagram = ph_diagram
    tp_diagram = pt_diagram
    st_diagram = ts_diagram

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "temperature_crit": {"method": None, "units": "K"},
                "pressure_crit": {"method": None, "units": "Pa"},
                "dens_mass_crit": {"method": None, "units": "kg/m^3"},
                "dens_mol_crit": {"method": None, "units": "mol/m^3"},
                "mw": {"method": None, "units": "kg/mol"},
                "temperature_sat": {"method": "None", "units": "K"},
                "flow_mol": {"method": None, "units": "mol/s"},
                "flow_mass": {"method": None, "units": "kg/s"},
                "flow_vol": {"method": None, "units": "m^3/s"},
                "temperature": {"method": None, "units": "K"},
                "pressure": {"method": None, "units": "Pa"},
                "dens_mass_phase": {"method": None, "units": "kg/m^3"},
                "temperature_red": {"method": None, "units": None},
                "pressure_sat": {"method": None, "units": "kPa"},
                "energy_internal_mol_phase": {"method": None, "units": "J/mol"},
                "enth_mol_phase": {"method": None, "units": "J/mol"},
                "entr_mol_phase": {"method": None, "units": "J/mol.K"},
                "cp_mol_phase": {"method": None, "units": "J/mol.K"},
                "cv_mol_phase": {"method": None, "units": "J/mol.K"},
                "energy_internal_mass_phase": {"method": None, "units": "J/kg"},
                "enth_mass_phase": {"method": None, "units": "J/kg"},
                "entr_mass_phase": {"method": None, "units": "J/kg.K"},
                "cp_mass_phase": {"method": None, "units": "J/kg.K"},
                "cv_mass_phase": {"method": None, "units": "J/kg.K"},
                "dens_mol_phase": {"method": None, "units": "mol/m^3"},
                "therm_cond_phase": {"method": None, "units": "W/m.K"},
                "visc_d_phase": {"method": None, "units": "Pa.s"},
                "visc_k_phase": {"method": None, "units": "m^2/s"},
                "phase_frac": {"method": None, "units": None},
                "flow_mol_comp": {"method": None, "units": "mol/s"},
                "flow_mass_comp": {"method": None, "units": "kg/s"},
                "energy_internal_mol": {"method": None, "units": "J/mol"},
                "enth_mol": {"method": None, "units": "J/mol"},
                "entr_mol": {"method": None, "units": "J/mol.K"},
                "cp_mol": {"method": None, "units": "J/mol.K"},
                "cv_mol": {"method": None, "units": "J/mol.K"},
                "energy_internal_mass": {"method": None, "units": "J/kg"},
                "enth_mass": {"method": None, "units": "J/kg"},
                "entr_mass": {"method": None, "units": "J/kg.K"},
                "cp_mass": {"method": None, "units": "J/kg.K"},
                "cv_mass": {"method": None, "units": "J/kg.K"},
                "heat_capacity_ratio": {"method": None, "units": None},
                "dens_mass": {"method": None, "units": "kg/m^3"},
                "dens_mol": {"method": None, "units": "mol/m^3"},
            }
        )

        obj.define_custom_properties(
            {
                "temperature_star": {
                    "method": None,
                    "units": obj.derived_units.TEMPERATURE,
                },
                "dens_mass_star": {
                    "method": None,
                    "units": obj.derived_units.DENSITY_MASS,
                },
                "dens_mol_star": {
                    "method": None,
                    "units": obj.derived_units.DENSITY_MOLE,
                },
                "specific_gas_constant": {
                    "method": None,
                    "units": obj.derived_units.ENTROPY_MASS,
                },
                "speed_sound_phase": {
                    "method": None,
                    "units": obj.derived_units.VELOCITY,
                },
                "vapor_frac": {"method": None, "units": pyo.units.dimensionless},
                "dh_vap_mol": {"method": None, "units": obj.derived_units.ENERGY_MOLE},
                "dh_vap_mass": {"method": None, "units": obj.derived_units.ENERGY_MASS},
            }
        )

        obj.add_default_units(
            {
                "time": pyo.units.s,
                "length": pyo.units.m,
                "mass": pyo.units.kg,
                "amount": pyo.units.mol,
                "temperature": pyo.units.K,
            }
        )
