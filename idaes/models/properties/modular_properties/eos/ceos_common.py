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

import os
import enum
from pyomo.environ import ExternalFunction, units as pyunits
from idaes import bin_directory

cubic_so_path = os.path.join(bin_directory, "cubic_roots.so")


def cubic_roots_available():
    """Make sure the compiled cubic root functions are available. Yes, in
    Windows the .so extention is still used.
    """
    return os.path.isfile(cubic_so_path)


class CubicType(enum.Enum):
    PR = 0
    SRK = 1


EoS_param = {
    CubicType.PR: {"u": 2, "w": -1, "omegaA": 0.45724, "coeff_b": 0.07780},
    CubicType.SRK: {"u": 1, "w": 0, "omegaA": 0.42748, "coeff_b": 0.08664},
}


class _ExternalFunctionSpecs(object):
    """Stores the information for creating an ExternalFunction object including.
    the library, function, and units of measure for the arguments and the return
    values.
    """

    def __init__(self, func, lib, units=None, arg_units=None):
        self.func = func
        self.lib = lib
        self.units = units
        self.arg_units = arg_units

    def kwargs(self):
        return {
            "library": self.lib,
            "function": self.func,
            "units": self.units,
            "arg_units": self.arg_units,
        }


class CubicThermoExpressions(object):
    """This class provides some standard expressions for cubic equations of
    state.
    """

    _external = {
        "compress_fact_liq_func": _ExternalFunctionSpecs(
            "cubic_root_l",
            cubic_so_path,
            units=pyunits.dimensionless,
            arg_units=[
                pyunits.dimensionless,
                pyunits.dimensionless,
                pyunits.dimensionless,
            ],
        ),
        "compress_fact_vap_func": _ExternalFunctionSpecs(
            "cubic_root_h",
            cubic_so_path,
            units=pyunits.dimensionless,
            arg_units=[
                pyunits.dimensionless,
                pyunits.dimensionless,
                pyunits.dimensionless,
            ],
        ),
    }

    def __init__(self, blk):
        self.blk = blk

    def add_funcs(self, names=None):
        if names is None:
            names = []
        elif isinstance(names, str):
            names = [names]
        for name in names:
            if hasattr(self.blk, name):
                continue
            setattr(self.blk, name, ExternalFunction(**self._external[name].kwargs()))

    def z_liq(self, eos=None, u=None, w=None, A=None, B=None, b=None, c=None, d=None):
        self.add_funcs(names="compress_fact_liq_func")
        if b is not None and c is not None and d is not None:
            pass
        else:
            if eos is not None:
                u = EoS_param[eos]["u"]
                w = EoS_param[eos]["w"]
            if u is None or w is None or A is None or B is None:
                raise RuntimeError(
                    "Please supply args: eos, {u, w, A, B}, or {b, c, d}"
                )
            b = -(1.0 + B - u * B)
            c = A + w * B**2 - u * B - u * B**2
            d = -A * B - w * B**2 - w * B**3
        return self.blk.compress_fact_liq_func(b, c, d)

    def z_vap(self, eos=None, u=None, w=None, A=None, B=None, b=None, c=None, d=None):
        self.add_funcs(names="compress_fact_vap_func")
        if b is not None and c is not None and d is not None:
            pass
        else:
            if eos is not None:
                u = EoS_param[eos]["u"]
                w = EoS_param[eos]["w"]
            if u is None or w is None or A is None or B is None:
                raise RuntimeError(
                    "Please supply args: eos, {u, w, A, B}, or {b, c, d}"
                )
            b = -(1.0 + B - u * B)
            c = A + w * B**2 - u * B - u * B**2
            d = -A * B - w * B**2 - w * B**3
        return self.blk.compress_fact_vap_func(b, c, d)
