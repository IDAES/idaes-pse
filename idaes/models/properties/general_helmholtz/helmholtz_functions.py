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
    viscosity_available,
    thermal_conductivity_available,
    surface_tension_available,
    component_registered,
)
from idaes.models.properties.general_helmholtz.components.parameters import (
    get_parameter_path,
    auto_register,
)
import idaes.logger as idaeslog
from idaes.models.properties.general_helmholtz.helmholtz_functions_map import (
    external_function_map as _external_function_map,
)


_log = idaeslog.getLogger(__name__)

_data_dir = get_parameter_path()
_data_dir = os.path.join(_data_dir, "")
auto_register()

# General Helmholtz functions return variables for all phases,
# but single phase properties do not need all of these.
# pylint: disable=W0612

try:
    _flib = find_library("general_helmholtz_external")
    ctypes.cdll.LoadLibrary(_flib)
except Exception:  # pylint: disable=W0703
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
    Enum, state variable set options.

    * PH: Pressure and enthalpy
    * PS: Pressure and entropy
    * PU: Pressure and internal energy
    * TPX: Temperature, pressure, and quality
    """

    PH = 1  # Pressure, Enthalpy
    PS = 2  # Pressure, Entropy
    PU = 3  # Pressure, Internal Energy
    TPX = 4  # Temperature, Pressure, Quality


class PhaseType(enum.Enum):
    """
    Enum, possible phases and presentation.

    * MIX: Two phase presented and a single phase to framework
    * LG: Two phases
    * L: Liquid only
    * G: Vapor only
    """

    MIX = 1  # Looks like a single phase called mixed with a vapor fraction
    LG = 2  # Looks like two phases vapor and liquid
    L = 3  # Assume only liquid is present
    G = 4  # Assume only vapor is present


class AmountBasis(enum.Enum):
    """
    Enum, mass or mole basis

    * MOLE: Amount is measured in moles
    * MASS: Amount is measured in mass
    """

    MOLE = 1
    MASS = 2


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

    Args:
        blk (Block): block to attach the external functions to
        parameters (HelmholtzParameterBlock): property parameter block
        amount_basis (AmountBasis|None): If none get the amount basis
            from the parameter block, otherwise use this amount basis for
            inputs.

    Returns:
        HelmholtzThermoExpressions
    """

    def __init__(self, blk, parameters, amount_basis=None):
        """Create a new thermodynamic property expression writer class.

        Args:
            blk (Block): block to attach the external functions to
            parameters (HelmholtzParameterBlock): property parameter block
            amount_basis (AmountBasis|None): If none get the amount basis
                from the parameter block, otherwise use this amount basis for
                inputs.

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

    def add_funcs(self, names=None):
        """Add external functions for the block expressions will be written
        from.

        Args:
            name (str): function name
        """
        add_helmholtz_external_functions(self.blk, names=names)

    @staticmethod
    def _validate_args(kwargs):
        """Make sure the provided arguments in kwargs is valid. For now
        the supported options for state variables are:
        (h, p), (s, p), (u, p), (T, p, x), and for single phase only (T, p)
        For sat properties, (T, x) and (p, x) are also valid, Psat or Tsat
        will be added automatically as the missing sate variable.
        """
        new_kwargs = {a: kwargs[a] for a in kwargs if kwargs[a] is not None}
        kwset = set(new_kwargs.keys())
        if len(new_kwargs) == 3:
            if {"T", "p", "x"} == set(new_kwargs):
                return new_kwargs, kwset
        if len(new_kwargs) != 2:
            raise RuntimeError(
                f"Helmholtz EoS expression writer: state variable set {kwset} is not supported."
            )
        if kwset in [
            {"T", "p"},
            {"h", "p"},
            {"s", "p"},
            {"u", "p"},
            {"T", "x"},
            {"p", "x"},
        ]:
            return new_kwargs, kwset
        raise RuntimeError(
            f"Helmholtz EoS expression writer: state variable set {kwset} is not supported."
        )

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
        result_basis = kwargs.pop("result_basis", self.amount_basis)
        convert_args = kwargs.pop("convert_args", True)
        kwargs, kwset = self._validate_args(kwargs)
        c = self.param.pure_component  # string for chemical component
        blk = self.blk  # block with external functions
        if "p" in kwargs:
            if convert_args:
                p = kwargs["p"] * self.param.uc["Pa to kPa"]
            else:
                p = kwargs["p"]
        if "x" in kwargs:
            x = kwargs["x"]
        else:
            x = None
        if "T" in kwargs:
            T = kwargs["T"]
        if {"h", "p"} == kwset:
            if convert_args:
                if self.amount_basis == AmountBasis.MOLE:
                    h = kwargs["h"] * self.param.uc["J/mol to kJ/kg"]
                else:
                    h = kwargs["h"] * self.param.uc["J/kg to kJ/kg"]
            else:
                h = kwargs["h"]
            return StateVars.PH, result_basis, blk, c, h, p, None
        if {"s", "p"} == kwset:
            if convert_args:
                if self.amount_basis == AmountBasis.MOLE:
                    s = kwargs["s"] * self.param.uc["J/mol/K to kJ/kg/K"]
                else:
                    s = kwargs["s"] * self.param.uc["J/kg/K to kJ/kg/K"]
            else:
                s = kwargs["s"]
            return StateVars.PS, result_basis, blk, c, s, p, None
        if {"u", "p"} == kwset:
            if convert_args:
                if self.amount_basis == AmountBasis.MOLE:
                    u = kwargs["u"] * self.param.uc["J/mol to kJ/kg"]
                else:
                    u = kwargs["u"] * self.param.uc["J/kg to kJ/kg"]
            else:
                u = kwargs["u"]
            return StateVars.PU, result_basis, blk, c, u, p, None
        if {"T", "x"} == kwset:
            self.add_funcs(names=["p_sat_t_func"])
            p = blk.p_sat_t_func(c, T, _data_dir)
            return StateVars.TPX, result_basis, blk, c, T, p, x
        if {"p", "x"} == kwset:
            self.add_funcs(names=["t_sat_func"])
            T = blk.t_sat_func(c, p, _data_dir)
            return StateVars.TPX, result_basis, blk, c, T, p, x
        if kwset in [{"T", "p"}, {"T", "p", "x"}]:
            return StateVars.TPX, result_basis, blk, c, T, p, x

    def _generic_prop(
        self,
        hp_func,
        up_func,
        sp_func,
        tp_liq_func,
        tp_vap_func,
        mass_uc,
        mole_uc,
        **kwargs,
    ):
        """mixed phase property"""
        sv, result_basis, blk, c, state1, p, x = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=[hp_func])
            prop = getattr(blk, hp_func)(c, state1, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=[up_func])
            prop = getattr(blk, up_func)(c, state1, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=[sp_func])
            prop = getattr(blk, sp_func)(c, state1, p, _data_dir)
        else:
            if x is None:
                raise RuntimeError("Cannot calculate vapor fraction from T and p")
            self.add_funcs(names=[tp_liq_func])
            self.add_funcs(names=[tp_vap_func])
            prop = (
                getattr(blk, tp_liq_func)(c, state1, p, _data_dir) * (1 - x)
                + getattr(blk, tp_vap_func)(c, state1, p, _data_dir) * x
            )
        if result_basis == AmountBasis.MOLE:
            return prop * mole_uc
        return prop * mass_uc

    def _generic_prop_phase(
        self,
        hp_func,
        up_func,
        sp_func,
        tp_func,
        mass_uc,
        mole_uc,
        **kwargs,
    ):
        """single phase property"""
        sv, result_basis, blk, c, state1, p, x = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=[hp_func])
            prop = getattr(blk, hp_func)(c, state1, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=[up_func])
            prop = getattr(blk, up_func)(c, state1, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=[sp_func])
            prop = getattr(blk, sp_func)(c, state1, p, _data_dir)
        else:
            self.add_funcs(names=[tp_func])
            prop = getattr(blk, tp_func)(c, state1, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return prop * mole_uc
        return prop * mass_uc

    def p(self, **kwargs):
        """Pressure"""
        # The only way we'd be calculating pressure is when we have T, x. In
        # that case, it has been calculated already, so just return pressure.
        sv, result_basis, blk, c, u, p, x = self._state_vars(**kwargs)
        return p * self.param.uc["kPa to Pa"]

    def T(self, **kwargs):
        """Temperature"""
        sv, results_basis, blk, c, state1, p, x = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["t_hp_func"])
            t = blk.t_hp_func(c, state1, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["t_up_func"])
            t = blk.t_up_func(c, state1, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["t_sp_func"])
            t = blk.t_sp_func(c, state1, p, _data_dir)
        else:
            t = state1
        return t

    def tau(self, **kwargs):
        """Critical Temperature (K)/Temperature (K)"""
        sv, result_basis, blk, c, state1, p, x = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["tau_hp_func"])
            tau = blk.tau_hp_func(c, state1, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["tau_up_func"])
            tau = blk.tau_up_func(c, state1, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["tau_sp_func"])
            tau = blk.tau_sp_func(c, state1, p, _data_dir)
        else:
            tau = self.param.temperature_star / state1
        return tau

    def x(self, **kwargs):
        """Vapor fraction"""
        sv, result_basis, blk, c, state1, p, x = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["vf_hp_func"])
            x = blk.vf_hp_func(c, state1, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["vf_up_func"])
            x = blk.vf_up_func(c, state1, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["vf_up_func"])
            x = blk.vf_up_func(c, state1, p, _data_dir)
        elif sv == StateVars.TPX:
            if x is None:
                raise RuntimeError("There is no nice way to calculate x from T and P")
        return x

    def s(self, **kwargs):
        """Mixed phase entropy"""
        sv, result_basis, blk, c, state1, p, x = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["s_hp_func"])
            s = blk.s_hp_func(c, state1, p, _data_dir)
        elif sv == StateVars.PU:
            self.add_funcs(names=["s_up_func"])
            s = blk.s_up_func(c, state1, p, _data_dir)
        elif sv == StateVars.PS:
            s = state1
        elif sv == StateVars.TPX:
            if x is None:
                raise RuntimeError("Cannot calculate vapor fraction from T and p")
            self.add_funcs(names=["s_vap_tp_func"])
            self.add_funcs(names=["s_liq_tp_func"])
            s = (
                blk.s_liq_tp_func(c, state1, p, _data_dir) * (1 - x)
                + blk.s_vap_tp_func(c, state1, p, _data_dir) * x
            )
        if result_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def s_liq(self, **kwargs):
        """Liquid phase entropy"""
        return self._generic_prop_phase(
            hp_func="s_liq_hp_func",
            up_func="s_liq_up_func",
            sp_func="s_liq_sp_func",
            tp_func="s_liq_tp_func",
            mass_uc=self.param.uc["kJ/kg/K to J/kg/K"],
            mole_uc=self.param.uc["kJ/kg/K to J/mol/K"],
            **kwargs,
        )

    def s_vap(self, **kwargs):
        """Vapor phase entropy"""
        return self._generic_prop_phase(
            hp_func="s_vap_hp_func",
            up_func="s_vap_up_func",
            sp_func="s_vap_sp_func",
            tp_func="s_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg/K to J/kg/K"],
            mole_uc=self.param.uc["kJ/kg/K to J/mol/K"],
            **kwargs,
        )

    def h(self, **kwargs):
        """Mixed phase enthalpy"""
        sv, result_basis, blk, c, state1, p, x = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            h = state1
        elif sv == StateVars.PU:
            self.add_funcs(names=["h_up_func"])
            h = blk.h_up_func(c, state1, p, _data_dir)
        elif sv == StateVars.PS:
            self.add_funcs(names=["h_sp_func"])
            h = blk.h_sp_func(c, state1, p, _data_dir)
        elif sv == StateVars.TPX:
            if x is None:
                raise RuntimeError("Cannot calculate vapor fraction from T and p")
            self.add_funcs(names=["h_vap_tp_func"])
            self.add_funcs(names=["h_liq_tp_func"])
            h = (
                blk.h_liq_tp_func(c, state1, p, _data_dir) * (1 - x)
                + blk.h_vap_tp_func(c, state1, p, _data_dir) * x
            )
        if result_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg to J/mol"]
        return h * self.param.uc["kJ/kg to J/kg"]

    def h_liq(self, **kwargs):
        """Liquid phase enthalpy"""
        return self._generic_prop_phase(
            hp_func="h_liq_hp_func",
            up_func="h_liq_up_func",
            sp_func="h_liq_sp_func",
            tp_func="h_liq_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def h_vap(self, **kwargs):
        """Vapor phase enthalpy"""
        return self._generic_prop_phase(
            hp_func="h_vap_hp_func",
            up_func="h_vap_up_func",
            sp_func="h_vap_sp_func",
            tp_func="h_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def u(self, **kwargs):
        """Mixed phase internal energy"""
        sv, result_basis, blk, c, state1, p, x = self._state_vars(**kwargs)
        if sv == StateVars.PH:
            self.add_funcs(names=["u_hp_func"])
            u = blk.u_hp_func(c, state1, p, _data_dir)
        elif sv == StateVars.PU:
            u = state1
        elif sv == StateVars.PS:
            self.add_funcs(names=["u_sp_func"])
            u = blk.u_sp_func(c, state1, p, _data_dir)
        elif sv == StateVars.TPX:
            if x is None:
                raise RuntimeError("Cannot calculate vapor fraction from T and p")
            self.add_funcs(names=["u_vap_tp_func"])
            self.add_funcs(names=["u_liq_tp_func"])
            u = (
                blk.u_liq_tp_func(c, state1, p, _data_dir) * (1 - x)
                + blk.u_vap_tp_func(c, state1, p, _data_dir) * x
            )
        if result_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg to J/mol"]
        return u * self.param.uc["kJ/kg to J/kg"]

    def u_liq(self, **kwargs):
        """Liquid phase internal energy"""
        return self._generic_prop_phase(
            hp_func="u_liq_hp_func",
            up_func="u_liq_up_func",
            sp_func="u_liq_sp_func",
            tp_func="u_liq_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def u_vap(self, **kwargs):
        """Vapor phase internal energy"""
        return self._generic_prop_phase(
            hp_func="u_vap_hp_func",
            up_func="u_vap_up_func",
            sp_func="u_vap_sp_func",
            tp_func="u_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def g(self, **kwargs):
        """Mixed phase Gibbs free energy"""
        return self._generic_prop(
            hp_func="g_hp_func",
            up_func="g_up_func",
            sp_func="g_sp_func",
            tp_liq_func="g_liq_tp_func",
            tp_vap_func="g_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def g_liq(self, **kwargs):
        """Liquid phase Gibbs free energy"""
        return self._generic_prop_phase(
            hp_func="g_liq_hp_func",
            up_func="g_liq_up_func",
            sp_func="g_liq_sp_func",
            tp_func="g_liq_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def g_vap(self, **kwargs):
        """Vapor phase Gibbs free energy"""
        return self._generic_prop_phase(
            hp_func="g_vap_hp_func",
            up_func="g_vap_up_func",
            sp_func="g_vap_sp_func",
            tp_func="g_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def f(self, **kwargs):
        """Mixed phase Helmholtz free energy"""
        return self._generic_prop(
            hp_func="f_hp_func",
            up_func="f_up_func",
            sp_func="f_sp_func",
            tp_liq_func="f_liq_tp_func",
            tp_vap_func="f_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def f_liq(self, **kwargs):
        """Liquid phase Helmholtz free energy"""
        return self._generic_prop_phase(
            hp_func="f_liq_hp_func",
            up_func="f_liq_up_func",
            sp_func="f_liq_sp_func",
            tp_func="f_liq_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def f_vap(self, **kwargs):
        """Vapor phase Helmholtz free energy"""
        return self._generic_prop_phase(
            hp_func="f_vap_hp_func",
            up_func="f_vap_up_func",
            sp_func="f_vap_sp_func",
            tp_func="f_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg to J/kg"],
            mole_uc=self.param.uc["kJ/kg to J/mol"],
            **kwargs,
        )

    def v(self, **kwargs):
        """Mixed phase volume"""
        return self._generic_prop(
            hp_func="v_hp_func",
            up_func="v_up_func",
            sp_func="v_sp_func",
            tp_liq_func="v_liq_tp_func",
            tp_vap_func="v_vap_tp_func",
            mass_uc=1,
            mole_uc=self.param.uc["m3/kg to m3/mol"],
            **kwargs,
        )

    def v_liq(self, **kwargs):
        """Liquid phase volume"""
        return self._generic_prop_phase(
            hp_func="v_liq_hp_func",
            up_func="v_liq_up_func",
            sp_func="v_liq_sp_func",
            tp_func="v_liq_tp_func",
            mass_uc=1,
            mole_uc=self.param.uc["m3/kg to m3/mol"],
            **kwargs,
        )

    def v_vap(self, **kwargs):
        """Vapor phase volume"""
        return self._generic_prop_phase(
            hp_func="v_vap_hp_func",
            up_func="v_vap_up_func",
            sp_func="v_vap_sp_func",
            tp_func="v_vap_tp_func",
            mass_uc=1,
            mole_uc=self.param.uc["m3/kg to m3/mol"],
            **kwargs,
        )

    def rho(self, **kwargs):
        """Mixed phase density"""
        return 1 / self.v(**kwargs)

    def rho_liq(self, **kwargs):
        """Liquid density"""
        return 1 / self.v_liq(**kwargs)

    def rho_vap(self, **kwargs):
        """Vapor Density"""
        return 1 / self.v_vap(**kwargs)

    def cp(self, **kwargs):
        """Mixed phase isobaric heat capacity"""
        return self._generic_prop(
            hp_func="cp_hp_func",
            up_func="cp_up_func",
            sp_func="cp_sp_func",
            tp_liq_func="cp_liq_tp_func",
            tp_vap_func="cp_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg/K to J/kg/K"],
            mole_uc=self.param.uc["kJ/kg/K to J/mol/K"],
            **kwargs,
        )

    def cp_liq(self, **kwargs):
        """Liquid phase isobaric heat capacity"""
        return self._generic_prop_phase(
            hp_func="cp_liq_hp_func",
            up_func="cp_liq_up_func",
            sp_func="cp_liq_sp_func",
            tp_func="cp_liq_tp_func",
            mass_uc=self.param.uc["kJ/kg/K to J/kg/K"],
            mole_uc=self.param.uc["kJ/kg/K to J/mol/K"],
            **kwargs,
        )

    def cp_vap(self, **kwargs):
        """Vapor phase isobaric heat capacity"""
        return self._generic_prop_phase(
            hp_func="cp_vap_hp_func",
            up_func="cp_vap_up_func",
            sp_func="cp_vap_sp_func",
            tp_func="cp_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg/K to J/kg/K"],
            mole_uc=self.param.uc["kJ/kg/K to J/mol/K"],
            **kwargs,
        )

    def cv(self, **kwargs):
        """Mixed phase isochoric heat capacity"""
        return self._generic_prop(
            hp_func="cv_hp_func",
            up_func="cv_up_func",
            sp_func="cv_sp_func",
            tp_liq_func="cv_liq_tp_func",
            tp_vap_func="cv_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg/K to J/kg/K"],
            mole_uc=self.param.uc["kJ/kg/K to J/mol/K"],
            **kwargs,
        )

    def cv_liq(self, **kwargs):
        """Liquid phase isochoric heat capacity"""
        return self._generic_prop_phase(
            hp_func="cv_liq_hp_func",
            up_func="cv_liq_up_func",
            sp_func="cv_liq_sp_func",
            tp_func="cv_liq_tp_func",
            mass_uc=self.param.uc["kJ/kg/K to J/kg/K"],
            mole_uc=self.param.uc["kJ/kg/K to J/mol/K"],
            **kwargs,
        )

    def cv_vap(self, **kwargs):
        """Vapor phase isochoric heat capacity"""
        return self._generic_prop_phase(
            hp_func="cv_vap_hp_func",
            up_func="cv_vap_up_func",
            sp_func="cv_vap_sp_func",
            tp_func="cv_vap_tp_func",
            mass_uc=self.param.uc["kJ/kg/K to J/kg/K"],
            mole_uc=self.param.uc["kJ/kg/K to J/mol/K"],
            **kwargs,
        )

    def w(self, **kwargs):
        """Mixed phase speed of sound (value for two-phase is not meaningful)
        use the liquid or vapor version if two phases are expected.
        """
        return self._generic_prop(
            hp_func="w_hp_func",
            up_func="w_up_func",
            sp_func="w_sp_func",
            tp_liq_func="w_liq_tp_func",
            tp_vap_func="w_vap_tp_func",
            mass_uc=1,
            mole_uc=1,
            **kwargs,
        )

    def w_liq(self, **kwargs):
        """Liquid phase speed of sound"""
        return self._generic_prop_phase(
            hp_func="w_liq_hp_func",
            up_func="w_liq_up_func",
            sp_func="w_liq_sp_func",
            tp_func="w_liq_tp_func",
            mass_uc=1,
            mole_uc=1,
            **kwargs,
        )

    def w_vap(self, **kwargs):
        """Vapor phase speed of sound"""
        return self._generic_prop_phase(
            hp_func="w_vap_hp_func",
            up_func="w_vap_up_func",
            sp_func="w_vap_sp_func",
            tp_func="w_vap_tp_func",
            mass_uc=1,
            mole_uc=1,
            **kwargs,
        )

    def viscosity(self, **kwargs):
        """Mixed phase viscosity (value for two-phase is not meaningful)
        use the liquid or vapor version if two phases are expected.
        """
        if not viscosity_available(self.param.pure_component):
            raise RuntimeError(
                f"Viscosity not available for {self.param.pure_component}"
            )
        return self._generic_prop(
            hp_func="mu_hp_func",
            up_func="mu_up_func",
            sp_func="mu_sp_func",
            tp_liq_func="mu_liq_tp_func",
            tp_vap_func="mu_vap_tp_func",
            mass_uc=self.param.uc["uPa to Pa"],
            mole_uc=self.param.uc["uPa to Pa"],
            **kwargs,
        )

    def viscosity_liq(self, **kwargs):
        """Liquid phase viscosity"""
        if not viscosity_available(self.param.pure_component):
            raise RuntimeError(
                f"Viscosity not available for {self.param.pure_component}"
            )
        return self._generic_prop_phase(
            hp_func="mu_liq_hp_func",
            up_func="mu_liq_up_func",
            sp_func="mu_liq_sp_func",
            tp_func="mu_liq_tp_func",
            mass_uc=self.param.uc["uPa to Pa"],
            mole_uc=self.param.uc["uPa to Pa"],
            **kwargs,
        )

    def viscosity_vap(self, **kwargs):
        """Vapor phase viscosity"""
        if not viscosity_available(self.param.pure_component):
            raise RuntimeError(
                f"Viscosity not available for {self.param.pure_component}"
            )
        return self._generic_prop_phase(
            hp_func="mu_vap_hp_func",
            up_func="mu_vap_up_func",
            sp_func="mu_vap_sp_func",
            tp_func="mu_vap_tp_func",
            mass_uc=self.param.uc["uPa to Pa"],
            mole_uc=self.param.uc["uPa to Pa"],
            **kwargs,
        )

    def thermal_conductivity(self, **kwargs):
        """Mixed phase thermal conductivity (value for two-phase is not meaningful)
        use the liquid or vapor version if two phases are expected.
        """
        if not thermal_conductivity_available(self.param.pure_component):
            raise RuntimeError(
                f"Thermal conductivity not available for {self.param.pure_component}"
            )
        return self._generic_prop(
            hp_func="lambda_hp_func",
            up_func="lambda_up_func",
            sp_func="lambda_sp_func",
            tp_liq_func="lambda_liq_tp_func",
            tp_vap_func="lambda_vap_tp_func",
            mass_uc=self.param.uc["mW to W"],
            mole_uc=self.param.uc["mW to W"],
            **kwargs,
        )

    def thermal_conductivity_liq(self, **kwargs):
        """Liquid phase thermal conductivity"""
        if not thermal_conductivity_available(self.param.pure_component):
            raise RuntimeError(
                f"Thermal conductivity not available for {self.param.pure_component}"
            )
        return self._generic_prop_phase(
            hp_func="lambda_liq_hp_func",
            up_func="lambda_liq_up_func",
            sp_func="lambda_liq_sp_func",
            tp_func="lambda_liq_tp_func",
            mass_uc=self.param.uc["mW to W"],
            mole_uc=self.param.uc["mW to W"],
            **kwargs,
        )

    def thermal_conductivity_vap(self, **kwargs):
        """Vapor phase thermal conductivity"""
        if not thermal_conductivity_available(self.param.pure_component):
            raise RuntimeError(
                f"Thermal conductivity not available for {self.param.pure_component}"
            )
        return self._generic_prop_phase(
            hp_func="lambda_vap_hp_func",
            up_func="lambda_vap_up_func",
            sp_func="lambda_vap_sp_func",
            tp_func="lambda_vap_tp_func",
            mass_uc=self.param.uc["mW to W"],
            mole_uc=self.param.uc["mW to W"],
            **kwargs,
        )

    def surface_tension(self, **kwargs):
        """Surface tension, this is only meaningful for two-phase region"""
        if not surface_tension_available(self.param.pure_component):
            raise RuntimeError(
                f"Surface tension not available for {self.param.pure_component}"
            )
        return self._generic_prop(
            hp_func="sigma_hp_func",
            up_func="sigma_up_func",
            sp_func="sigma_sp_func",
            tp_liq_func="sigma_liq_tp_func",
            tp_vap_func="sigma_vap_tp_func",
            mass_uc=self.param.uc["mN to N"],
            mole_uc=self.param.uc["mN to N"],
            **kwargs,
        )

    def p_sat(self, T):
        """Return saturation pressure as a function of T or tau"""
        self.add_funcs(names=["p_sat_t_func"])
        return (
            self.blk.p_sat_t_func(self.param.pure_component, T, _data_dir)
            * self.param.uc["kPa to Pa"]
        )

    def T_sat(self, p, convert_args=True):
        """Return saturation temperature as a function of p"""
        if convert_args:
            p *= self.param.uc["Pa to kPa"]
        self.add_funcs(names=["t_sat_func"])
        return self.blk.t_sat_func(self.param.pure_component, p, _data_dir)

    def h_vap_sat(self, T=None, p=None, result_basis=None, convert_args=True):
        """Return saturation vapor enthalpy as a function of T or p"""
        if result_basis is None:
            result_basis = self.amount_basis
        if T is not None:
            self.add_funcs(names=["h_vap_sat_t_func"])
            h = self.blk.h_vap_sat_t_func(self.param.pure_component, T, _data_dir)
        elif p is not None:
            self.add_funcs(names=["h_vap_sat_p_func"])
            if convert_args:
                p *= self.param.uc["Pa to kPa"]
            h = self.blk.h_vap_sat_p_func(self.param.pure_component, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg to J/mol"]
        return h * self.param.uc["kJ/kg to J/kg"]

    def h_liq_sat(self, T=None, p=None, result_basis=None, convert_args=True):
        """Return saturation liquid enthalpy as a function of T or p"""
        if result_basis is None:
            result_basis = self.amount_basis
        if T is not None:
            self.add_funcs(names=["h_liq_sat_t_func"])
            h = self.blk.h_liq_sat_t_func(self.param.pure_component, T, _data_dir)
        elif p is not None:
            self.add_funcs(names=["h_liq_sat_p_func"])
            if convert_args:
                p *= self.param.uc["Pa to kPa"]
            h = self.blk.h_liq_sat_p_func(self.param.pure_component, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return h * self.param.uc["kJ/kg to J/mol"]
        return h * self.param.uc["kJ/kg to J/kg"]

    def s_vap_sat(self, T=None, p=None, result_basis=None, convert_args=True):
        """Return saturation vapor entropy as a function of T or p"""
        if result_basis is None:
            result_basis = self.amount_basis
        if T is not None:
            self.add_funcs(names=["s_vap_sat_t_func"])
            s = self.blk.s_vap_sat_t_func(self.param.pure_component, T, _data_dir)
        elif p is not None:
            self.add_funcs(names=["s_vap_sat_p_func"])
            if convert_args:
                p *= self.param.uc["Pa to kPa"]
            s = self.blk.s_vap_sat_p_func(self.param.pure_component, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def s_liq_sat(self, T=None, p=None, result_basis=None, convert_args=True):
        """Return saturation liquid entropy as a function of T or p"""
        if result_basis is None:
            result_basis = self.amount_basis
        if T is not None:
            self.add_funcs(names=["s_liq_sat_t_func"])
            s = self.blk.s_liq_sat_t_func(self.param.pure_component, T, _data_dir)
        elif p is not None:
            self.add_funcs(names=["s_liq_sat_p_func"])
            if convert_args:
                p *= self.param.uc["Pa to kPa"]
            s = self.blk.s_liq_sat_p_func(self.param.pure_component, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return s * self.param.uc["kJ/kg/K to J/mol/K"]
        return s * self.param.uc["kJ/kg/K to J/kg/K"]

    def u_vap_sat(self, T=None, p=None, result_basis=None, convert_args=True):
        """Return saturation vapor internal energy as a function of T or p"""
        if result_basis is None:
            result_basis = self.amount_basis
        if T is not None:
            self.add_funcs(names=["u_vap_sat_t_func"])
            u = self.blk.u_vap_sat_t_func(self.param.pure_component, T, _data_dir)
        elif p is not None:
            self.add_funcs(names=["u_vap_sat_p_func"])
            if convert_args:
                p *= self.param.uc["Pa to kPa"]
            u = self.blk.u_vap_sat_p_func(self.param.pure_component, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg to J/mol"]
        return u * self.param.uc["kJ/kg to J/kg"]

    def u_liq_sat(self, T=None, p=None, result_basis=None, convert_args=True):
        """Return saturation liquid internal energy as a function of T or p"""
        if result_basis is None:
            result_basis = self.amount_basis
        if T is not None:
            self.add_funcs(names=["u_liq_sat_t_func"])
            u = self.blk.u_liq_sat_t_func(self.param.pure_component, T, _data_dir)
        elif p is not None:
            self.add_funcs(names=["u_liq_sat_p_func"])
            if convert_args:
                p *= self.param.uc["Pa to kPa"]
            u = self.blk.u_liq_sat_p_func(self.param.pure_component, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return u * self.param.uc["kJ/kg to J/mol"]
        return u * self.param.uc["kJ/kg to J/kg"]

    def v_vap_sat(self, T=None, p=None, result_basis=None, convert_args=True):
        """Return saturation vapor volume as a function of T or p"""
        if result_basis is None:
            result_basis = self.amount_basis
        if T is not None:
            self.add_funcs(names=["v_vap_sat_t_func"])
            v = self.blk.v_vap_sat_t_func(self.param.pure_component, T, _data_dir)
        elif p is not None:
            self.add_funcs(names=["v_vap_sat_p_func"])
            if convert_args:
                p *= self.param.uc["Pa to kPa"]
            v = self.blk.v_vap_sat_p_func(self.param.pure_component, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return v * self.param.uc["m3/kg to m3/mol"]
        return v

    def v_liq_sat(self, T=None, p=None, result_basis=None, convert_args=True):
        """Return saturation liquid volume as a function of T or p"""
        if result_basis is None:
            result_basis = self.amount_basis
        if T is not None:
            self.add_funcs(names=["v_liq_sat_t_func"])
            v = self.blk.v_liq_sat_t_func(self.param.pure_component, T, _data_dir)
        elif p is not None:
            self.add_funcs(names=["v_liq_sat_p_func"])
            if convert_args:
                p *= self.param.uc["Pa to kPa"]
            v = self.blk.v_liq_sat_p_func(self.param.pure_component, p, _data_dir)
        if result_basis == AmountBasis.MOLE:
            return v * self.param.uc["m3/kg to m3/mol"]
        return v


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
            doc="(str) Pure component for which to calculate properties",
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
**StateVars.PS** - Pressure-Entropy,
**StateVars.PU** - Pressure-Internal Energy,
**StateVars.TPX** - Temperature-Pressure-Quality}""",
        ),
    )

    CONFIG.declare(
        "amount_basis",
        ConfigValue(
            default=AmountBasis.MOLE,
            domain=In(AmountBasis),
            description="Quantities on a mass or mole basis",
            doc="""The amount basis (mass or mole) for quantities
**default** - AmountBasis.mole
**Valid values:** {
**AmountBasis.mole** - use mole units (mol),
**AmountBasis.mass** - use mass units (kg)}""",
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
        Convenience function to calculate a state variable from temperature and either
        pressure or vapor fraction. This function can be used for inlet streams and
        initialization where temperature is known instead of a state variable.
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
            prop (str): h, s, or u
        Returns:
            Total selected state variable.
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
        """
        Convenience method to calculate enthalpy from temperature and either
        pressure or vapor fraction. This function can be used for inlet streams and
        initialization where temperature is known instead of enthalpy.
        User must provide values for one of these sets of values: {T, P}, {T, x},
        or {P, x}.

        Args:
            T (float): Temperature
            P (float): Pressure, None if saturated
            x (float): Vapor fraction [mol vapor/mol total] (between 0 and 1), None if superheated or sub-cooled
            units (Units): The units to report the result in, if None use the default units appropriate for the amount basis.
            amount_basis (AmountBasis): Whether to use a mass or mole basis
            with_units (bool): if True return an expression with units

        Returns:
            float: Specific or molar enthalpy
        """
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
        """
        Convenience method to calculate entropy from temperature and either
        pressure or vapor fraction. This function can be used for inlet streams and
        initialization where temperature is known instead of entropy.
        User must provide values for one of these sets of values: {T, P}, {T, x},
        or {P, x}.

        Args:
            T (float): Temperature
            P (float): Pressure, None if saturated
            x (float): Vapor fraction [mol vapor/mol total] (between 0 and 1), None if superheated or sub-cooled
            units (Units): The units to report the result in, if None use the default units appropriate for the amount basis.
            amount_basis (AmountBasis): Whether to use a mass or mole basis
            with_units (bool): if True return an expression with units
        Returns:
            float: Specific or molar entropy
        """
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
        """
        Convenience method to calculate internal energy from temperature and either
        pressure or vapor fraction. This function can be used for inlet streams and
        initialization where temperature is known instead of internal energy.
        User must provide values for one of these sets of values: {T, P}, {T, x},
        or {P, x}.

        Args:
            T (float): Temperature
            P (float): Pressure, None if saturated
            x (float): Vapor fraction [mol vapor/mol total] (between 0 and 1), None if superheated or sub-cooled
            units (Units): The units to report the result in, if None use the default units appropriate for the amount basis.
            amount_basis (AmountBasis): Whether to use a mass or mole basis
            with_units (bool): if True return an expression with units

        Returns:
            float: Specific or molar internal energy
        """
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
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("dens_mass_phase", 1e-2, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e1)
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
        """Populate the parameter block"""
        if not self.available():
            raise RuntimeError("Helmholtz EoS external functions not available")
        super().build()
        # Check if the specified component is supported
        if not component_registered(self.config.pure_component):
            raise ConfigurationError(
                f"Component {self.config.pure_component} not supported."
            )
        # This is imported here to avoid a circular import
        # pylint: disable-next=import-outside-toplevel
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
            "m3/kg to m3/mol": self.mw,
            "uPa to Pa": 1e-6 * pyo.units.Pa / pyo.units.uPa,
            "mW to W": 1e-3 * pyo.units.W / pyo.units.mW,
            "mN to N": 1e-3 * pyo.units.N / pyo.units.mN,
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

        # Smoothing parameters for TPX complimentarity form
        if self.config.state_vars == StateVars.TPX:
            self.smoothing_pressure_over = pyo.Param(
                mutable=True,
                initialize=1e-4,
                doc="Smooth max parameter (pressure over)",
                units=pyo.units.Pa,
            )

            self.smoothing_pressure_under = pyo.Param(
                mutable=True,
                initialize=1e-4,
                doc="Smooth max parameter (pressure under)",
                units=pyo.units.Pa,
            )

    def add_param(self, name, expr):
        """Add a parameter to the block.

        Args:
            name (str): parameter name
            expr (expression): Pyomo expression for parameter value
        """
        self.add_component(
            name,
            pyo.Param(
                initialize=pyo.value(expr),
                units=pyo.units.get_units(expr),
                domain=pyo.Reals,
            ),
        )

    def initialize(self, *args, **kwargs):
        """No initialization required here. This method is included for
        compatibility.
        """

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
            dict: dictionary with the keys {'T', 'tau', 'p', 'delta_liq',
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
        # parameter block (self)
        if amount_basis is None:
            amount_basis = self.config.amount_basis

        # We'll do the calculation based on temperature.  We want more points
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
        # tau we can calculate all the rest of the properties for sat liquid
        delta_sat_l_vec = [
            pyo.value(self.delta_sat_l_func(self.pure_component, tau, _data_dir))
            for tau in tau_vec
        ]
        # Get the reduced density vector for the liquid, from this and
        # tau we can calculate all the rest of the properties for sat liquid
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
        """Create a enthalpy-pressure diagram using Matplotlib

        Args:
            ylim (tuple): lower and upper limits for pressure axis
            xlim (tuple): lower and upper limits for enthalpy axis
            points (dict): dict of tuples points to label on the plot
            figsize (tuple): figure size
            dpi (int): figure dots per inch
            isotherms (list|None): list of temperatures for plotting isotherms
            isotherms_line_format (str|None): line format for isotherms
            isotherms_label (bool): if true label isotherms

        Returns:
            (figure, axis)
        """
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
        """Create a entropy-temperautre diagram using Matplotlib

        Args:
            ylim (tuple): lower and upper limits for temperature axis
            xlim (tuple): lower and upper limits for entropy axis
            points (dict): dict of tuples points to label on the plot
            figsize (tuple): figure size
            dpi (int): figure dots per inch

        Returns:
            (figure, axis)
        """
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

    def pt_diagram(self, ylim=None, xlim=None, figsize=None, dpi=None):
        """Create a pressure-teperature diagram using Matplotlib

        Args:
            ylim (tuple): lower and upper limits for pressure axis
            xlim (tuple): lower and upper limits for temperature axis
            figsize (tuple): figure size
            dpi (int): figure dots per inch

        Returns:
            (figure, axis)
        """
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

    # In case you can't remember the axis order in the diagrams
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
                "surf_tens": {"method": None, "units": "N/m"},
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
                "speed_sound": {
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
