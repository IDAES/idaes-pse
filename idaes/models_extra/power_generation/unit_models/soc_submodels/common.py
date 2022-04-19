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

__author__ = "John Eslick, Douglas Allan"

import copy
import enum

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.dae import DerivativeVar
import pyomo.environ as pyo
from pyomo.network import Port


from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.math import safe_log
from idaes.core.util import get_solver
import idaes.core.util.model_statistics as mstat

from idaes.core.util.misc import VarLikeExpression
import idaes.logger as idaeslog

_constR = 8.3145 * pyo.units.J / pyo.units.mol / pyo.units.K  # or Pa*m3/K/mol
_constF = 96485 * pyo.units.coulomb / pyo.units.mol
_safe_log_eps = 1e-9
_safe_sqrt_eps = 1e-9

class CV_Bound(enum.Enum):
    EXTRAPOLATE = 1
    NODE_VALUE = 2


class CV_Interpolation(enum.Enum):
    UDS = 1  # Upwind difference scheme, exit face same as node center
    CDS = 2  # Linear interpolation from upstream and downstream node centers
    QUICK = 3  # Quadratic upwind interpolation, quadratic from two upwind
    # centers and one downwind


def _set_default_factor(c, s):
    for i in c:
        if iscale.get_scaling_factor(c[i]) is None:
            iscale.set_scaling_factor(c[i], s)


def _set_scaling_factor_if_none(c, s):
    if iscale.get_scaling_factor(c) is None:
        iscale.set_scaling_factor(c, s)


def _set_and_get_scaling_factor(c, s):
    _set_scaling_factor_if_none(c, s)
    return iscale.get_scaling_factor(c)


def _set_if_unfixed(v, val):
    if not v.fixed:
        v.set_value(pyo.value(val))


def _create_if_none(blk, var, idx_set, units):
    attr = getattr(blk.config, var)
    if attr is None:
        if idx_set is None:
            setattr(
                blk,
                var,
                pyo.Var(
                    initialize=0,
                    units=units,
                ),
            )
        else:
            setattr(
                blk,
                var,
                pyo.Var(
                    *idx_set,
                    initialize=0,
                    units=units,
                ),
            )
    else:
        setattr(blk, var, pyo.Reference(attr))

def _init_solve_block(blk, solver, log):
    if not mstat.degrees_of_freedom(blk) == 0:
        raise InitializationError(
            f"{blk.name} encountered a nonsquare problem during "
            "initialization. Check whether all parameter values are fixed and "
            "keep in mind that the SOC submodel initialization methods assume "
            "that certain degrees of freedom have been fixed by the cell "
            "initialization method."
        )
    with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
        results = solver.solve(blk, tee=slc.tee)
    if not pyo.check_optimal_termination(results):
        raise InitializationError(
            f"{blk.name} failed to initialize successfully. Please check "
            f"the output logs for more information."
        )


def _interpolate_channel(
    iz, ifaces, nodes, faces, phi_func, phi_inlet, method, opposite_flow
):
    """PRIVATE Function: Interpolate faces of control volumes in 1D

    Args:
        iz: index of the face
        ifaces: set of face indexes
        nodes: set of node z locations
        faces: set of face z locations
        phi_func: function that returns an expression for the quantity to be
            interpolated at the node as a function of node index
        phi_inlet: expression for the value of the quantity to be interpolated
            at the channel inlet face (not the inlet to the CV)
        method: interpolation method
        opposite_flow: if True assume velocity is negative

    Returns:
        expression for phi at face[iz]
    """
    # I don't always need these, but it doesn't take long to calculate them
    if not opposite_flow:
        izu = iz - 1  # adjacent node upstream of the face
        izuu = iz - 2  # node upstream adacjent to node upstream adjacent to face
        izd = iz  # downstream node adjacent to face
    else:
        izu = iz  # adjacent node upstream of the face
        izuu = iz + 1  # node upstream adacjent to node upstream adjacent to face
        izd = iz - 1  # downstream node adjacent to face
    if iz == ifaces.first() and not opposite_flow:
        return phi_inlet
    if iz == ifaces.last() and opposite_flow:
        return phi_inlet
    if method == CV_Interpolation.UDS:
        return phi_func(izu)
    if method == CV_Interpolation.CDS:
        zu = nodes.at(izu)
        if (opposite_flow and iz == ifaces.first()) or (
            not opposite_flow and iz == ifaces.last()
        ):
            izd = izuu
        zd = nodes.at(izd)
        zf = faces.at(iz)
        lambf = (zd - zf) / (zd - zu)
        return (1 - lambf) * phi_func(izu) + lambf * phi_func(izd)


def _interpolate_2D(
    ic,
    ifaces,
    nodes,
    faces,
    phi_func,
    phi_bound_0,
    phi_bound_1,
    derivative=False,
    method=CV_Interpolation.CDS,
):
    """PRIVATE Function: Interpolate faces of control volumes in 1D

    Args:
        ic: face index in x or z direction matching specified direction
        ifaces: set of face indexes
        nodes: set of node locations in specified direction
        faces: set of face locations in specified direction
        phi_func: function that returns an expression for the quantity to be
            interpolated at a node cneter as a function of node index
        phi_bound_0: expression for the value of the quantity to be interpolated
            at the 0 bound
        phi_bound_1: expression for the value of the quantity to be interpolated
            at the 1 bound
        derivative: If True estimate derivative
        method: interpolation method currently only CDS is supported
        direction: direction to interpolate

    Returns:
        expression for phi at face
    """
    if method != CV_Interpolation.CDS:
        raise RuntimeError(
            "SOFC/SOEC sections modeled in 2D do not have bulk "
            "flow, so currently only CDS interpolation is supported"
        )
    icu = ic - 1
    icd = ic
    if ic == ifaces.first():
        if not isinstance(phi_bound_0, CV_Bound):
            return phi_bound_0
        icu = ic + 1
    if ic == ifaces.last():
        if not isinstance(phi_bound_1, CV_Bound):
            return phi_bound_1
        icd = ic - 2
    # For now CDS is the only option
    cu = nodes.at(icu)
    cd = nodes.at(icd)
    if not derivative:
        cf = faces.at(ic)
        lambf = (cd - cf) / (cd - cu)
        return (1 - lambf) * phi_func(icu) + lambf * phi_func(icd)
    else:
        # Since we are doing linear interpolation derivative is the slope
        # between node centers even if they are not evenly spaced
        return (phi_func(icd) - phi_func(icu)) / (cd - cu)

_element_list = ["Ar", "H", "O", "C", "N", "S"]
_species_list = [
    "Ar",
    "CH4",
    "CO",
    "CO2",
    "C2H4",
    "C2H6",
    "C3H8",
    "H2",
    "H2O",
    "H2S",
    "NO",
    "N2",
    "N2O",
    "O2",
    "SO2",
]
_element_dict = {
    "Ar": {"Ar": 1},
    "H": {"CH4": 4, "C2H4": 4, "C2H6": 6, "C3H8": 8, "H2": 2, "H2O": 2, "H2S": 2},
    "O": {"CO": 1, "CO2": 2, "H2O": 1, "NO": 1, "N2O": 1, "O2": 2, "SO2": 2},
    "C": {"CH4": 1, "CO": 1, "CO2": 1, "C2H4": 4, "C2H6": 2, "C3H8": 3},
    "N": {"NO": 1, "N2": 2, "N2O": 2},
    "S": {"H2S": 1, "SO2": 1},
}

for value in _element_dict.values():
    for species in _species_list:
        if species not in value.keys():
            value[species] = 0

# Parameters for binary diffusion coefficients from:
#  "Properties of Gases and Liquids" 5th Ed.
bin_diff_sigma = {
    "Ar": 3.542,
    "Air": 3.711,
    "CH4": 3.758,
    "CO": 3.690,
    "CO2": 3.941,
    "C2H4": 4.163,
    "C2H6": 4.443,
    "C3H8": 5.118,
    "H2": 2.827,
    "H2O": 2.641,
    "H2S": 3.623,
    "NO": 3.492,
    "N2": 3.798,
    "N2O": 3.826,
    "O2": 3.467,
    "SO2": 4.112,
}
bin_diff_epsok = {
    "Ar": 93.3,
    "Air": 78.6,
    "CH4": 148.6,
    "CO": 91.7,
    "CO2": 195.2,
    "C2H4": 224.7,
    "C2H6": 215.7,
    "C3H8": 237.1,
    "H2": 59.7,
    "H2O": 809.1,
    "H2S": 301.1,
    "NO": 116.7,
    "N2": 71.4,
    "N2O": 232.4,
    "O2": 106.7,
    "SO2": 335.4,
}
bin_diff_M = {
    "Ar": 39.948,
    "Air": 28.96,
    "CH4": 16.0425,
    "CO": 28.0101,
    "CO2": 44.0095,
    "C2H4": 28.0532,
    "C2H6": 30.0690,
    "C3H8": 44.0956,
    "H2": 2.01588,
    "H2O": 18.0153,
    "H2S": 34.081,
    "NO": 30.0061,
    "N2": 28.0134,
    "N2O": 44.0128,
    "O2": 31.9988,
    "SO2": 64.064,
}

# Shomate equation parameters
# NIST Webbook
h_params = {
    "Ar": {
        "A": 20.78600,
        "B": 2.825911e-7,
        "C": -1.464191e-7,
        "D": 1.092131e-8,
        "E": -3.661371e-8,
        "F": -6.197350,
        "G": 179.9990,
        "H": 0.000000,
    },
    "CH4": {
        "A": -0.703029,
        "B": 108.4773,
        "C": -42.52157,
        "D": 5.862788,
        "E": 0.678565,
        "F": -76.84376,
        "G": 158.7163,
        "H": -74.87310,
    },
    "CO": {
        "A": 25.56759,
        "B": 6.096130,
        "C": 4.054656,
        "D": -2.671301,
        "E": 0.131021,
        "F": -118.0089,
        "G": 227.3665,
        "H": -110.5271,
    },
    "CO2": {
        "A": 24.99735,
        "B": 55.18696,
        "C": -33.69137,
        "D": 7.948387,
        "E": -0.136638,
        "F": -403.6075,
        "G": 228.2431,
        "H": -393.5224,
    },
    "C2H4": {
        "A": -6.387880,
        "B": 184.4019,
        "C": -112.9718,
        "D": 28.49593,
        "E": 0.315540,
        "F": 48.17332,
        "G": 163.1568,
        "H": 52.46694,
    },
    "H2": {  # 1000K-2500K
        "A": 33.066178,
        "B": -11.363417,
        "C": 11.432816,
        "D": -2.772874,
        "E": -0.158558,
        "F": -9.980797,
        "G": 172.707974,
        "H": 0.0,
    },
    "H2O": {  # 500K-1700K
        "A": 30.09200,
        "B": 6.832514,
        "C": 6.793435,
        "D": -2.534480,
        "E": 0.082139,
        "F": -250.8810,
        "G": 223.3967,
        "H": -241.8264,
    },
    "H2S": {
        "A": 26.88412,
        "B": 18.67809,
        "C": 3.434203,
        "D": -3.378702,
        "E": 0.135882,
        "F": -28.91211,
        "G": 233.3747,
        "H": -20.50202,
    },
    "N2": {
        "A": 19.50583,
        "B": 19.88705,
        "C": -8.598535,
        "D": 1.369784,
        "E": 0.527601,
        "F": -4.935202,
        "G": 212.3900,
        "H": 0.0,
    },
    "O2": {  # 700 to 2000
        "A": 30.03235,
        "B": 8.772972,
        "C": -3.988133,
        "D": 0.788313,
        "E": -0.741599,
        "F": -11.32468,
        "G": 236.1663,
        "H": 0.0,
    },
    "SO2": {
        "A": 21.43049,
        "B": 74.35094,
        "C": -57.75217,
        "D": 16.35534,
        "E": 0.086731,
        "F": -305.7688,
        "G": 254.8872,
        "H": -296.8422,
    },
    "Vac": {
        "A": 0.0,
        "B": 0.0,
        "C": 0.0,
        "D": 0.0,
        "E": 0.0,
        "F": 0.0,
        "G": 0.0,
        "H": 0.0,
    },
    "O^2-": {  # Chosen to match O2 coeffs besides enthalpy and entropy of formation
        "A": 0.5 * 30.03235,
        "B": 0.5 * 8.772972,
        "C": 0.5 * -3.988133,
        "D": 0.5 * 0.788313,
        "E": 0.5 * -0.741599,
        "F": 0.5 * -11.32468,
        "G": 0.5 * 236.1663,  # 0,
        "H": 0,  # -236.0,
    },
}


def _binary_diffusion_coefficient_expr(temperature, p, c1, c2):
    mab = 2 * (1.0 / bin_diff_M[c1] + 1.0 / bin_diff_M[c2]) ** (-1)
    sab = (bin_diff_sigma[c1] + bin_diff_sigma[c2]) / 2.0
    epsok = (bin_diff_epsok[c1] * bin_diff_epsok[c2]) ** 0.5
    tr = temperature / epsok
    a = 1.06036
    b = 0.15610
    c = 0.19300
    d = 0.47635
    e = 1.03587
    f = 1.52996
    g = 1.76474
    h = 3.89411
    omega = (
        a / tr**b + c / pyo.exp(d * tr) + e / pyo.exp(f * tr) + g / pyo.exp(h * tr)
    )
    cm2_to_m2 = 0.01 * 0.01
    Pa_to_bar = 1e-5
    return (
        0.002666
        * cm2_to_m2
        * temperature ** (3 / 2)
        / p
        / Pa_to_bar
        / mab**0.5
        / sab**2
        / omega
    )


def _comp_enthalpy_expr(temperature, comp):
    # ideal gas enthalpy
    d = h_params[comp]
    t = temperature / 1000.0 / pyo.units.K
    return (
        1000
        * (
            d["A"] * t
            + d["B"] * t**2 / 2.0
            + d["C"] * t**3 / 3.0
            + d["D"] * t**4 / 4.0
            - d["E"] / t
            + d["F"]
        )
        * pyo.units.J
        / pyo.units.mol
    )


def _comp_int_energy_expr(temperature, comp):
    # ideal gas enthalpy
    d = h_params[comp]
    t = temperature / 1000.0
    return _comp_enthalpy_expr(temperature, comp) - _constR * temperature


def _comp_entropy_expr(temperature, comp):
    # ideal gas enthalpy
    d = h_params[comp]
    t = temperature / 1000.0 / pyo.units.K
    return (
        (
            d["A"] * safe_log(t, eps=_safe_log_eps)
            + d["B"] * t
            + d["C"] * t**2 / 2.0
            + d["D"] * t**3 / 3.0
            - d["E"] / 2.0 / t**2
            + d["G"]
        )
        * pyo.units.J
        / pyo.units.mol
        / pyo.units.K
    )


def _submodel_boilerplate_config(CONFIG):
    CONFIG.declare(
        "cv_zfaces",
        ConfigValue(description="CV boundary set, should start with 0 and end with 1."),
    )
    CONFIG.declare(
        "length_z",
        ConfigValue(
            default=None, description="Length in the direction of flow (z-direction)"
        ),
    )
    CONFIG.declare(
        "length_y", ConfigValue(default=None, description="Width of cell (y-direction)")
    )
    CONFIG.declare(
        "current_density",
        ConfigValue(default=None, description="Optional current_density variable"),
    )
    CONFIG.declare(
        "include_temperature_x_thermo",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=True,
            description="Whether to consider temperature variations in "
            "x direction in thermodynamic equations",
        ),
    )


def _submodel_boilerplate_create_if_none(unit):
    tset = unit.flowsheet().config.time
    iznodes = unit.iznodes
    _create_if_none(unit, "length_z", idx_set=None, units=pyo.units.m)
    _create_if_none(unit, "length_y", idx_set=None, units=pyo.units.m)
    _create_if_none(
        unit,
        "current_density",
        idx_set=(tset, iznodes),
        units=pyo.units.A / pyo.units.m**2,
    )


def _thermal_boundary_conditions_config(CONFIG, thin):
    CONFIG.declare(
        "temperature_z",
        ConfigValue(
            default=None, description="Temperature as indexed by time " "and z"
        ),
    )
    if thin:
        CONFIG.declare(
            "Dtemp",
            ConfigValue(
                default=None,
                description="Deviation of temperature " "from temperature_z",
            ),
        )
    else:
        CONFIG.declare(
            "Dtemp_x0",
            ConfigValue(
                default=None,
                description="Deviation of temperature at x=0 " "from temperature_z",
            ),
        )
        CONFIG.declare(
            "Dtemp_x1",
            ConfigValue(
                default=None,
                description="Deviation of temperature at x=1 " "from temperature_z",
            ),
        )
    CONFIG.declare(
        "qflux_x0",
        ConfigValue(
            default=None, description="Heat flux through x=0 " "(positive is in)"
        ),
    )
    CONFIG.declare(
        "qflux_x1",
        ConfigValue(
            default=None, description="Heat flux through x=1 " "(positive is out)"
        ),
    )


def _create_thermal_boundary_conditions_if_none(unit, thin):
    tset = unit.flowsheet().config.time
    include_temp_x_thermo = unit.config.include_temperature_x_thermo
    iznodes = unit.iznodes

    _create_if_none(unit, "temperature_z", idx_set=(tset, iznodes), units=pyo.units.K)

    if thin:
        _create_if_none(unit, "Dtemp", idx_set=(tset, iznodes), units=pyo.units.K)

        @unit.Expression(tset, iznodes)
        def temperature(b, t, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.Dtemp[t, iz]
            else:
                return b.temperature_z[t, iz]

    else:
        _create_if_none(unit, "Dtemp_x0", idx_set=(tset, iznodes), units=pyo.units.K)

        @unit.Expression(tset, iznodes)
        def temperature_x0(b, t, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.Dtemp_x0[t, iz]
            else:
                return b.temperature_z[t, iz]

        _create_if_none(unit, "Dtemp_x1", idx_set=(tset, iznodes), units=pyo.units.K)

        @unit.Expression(tset, iznodes)
        def temperature_x1(b, t, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.Dtemp_x1[t, iz]
            else:
                return b.temperature_z[t, iz]

    _create_if_none(
        unit, "qflux_x0", idx_set=(tset, iznodes), units=pyo.units.W / pyo.units.m**2
    )

    _create_if_none(
        unit, "qflux_x1", idx_set=(tset, iznodes), units=pyo.units.W / pyo.units.m**2
    )


def _material_boundary_conditions_config(CONFIG, thin):
    if thin:
        CONFIG.declare(
            "xflux",
            ConfigValue(default=None, description="Variable for material flux"),
        )
        CONFIG.declare(
            "Dconc",
            ConfigValue(
                default=None,
                description="Deviation of concentration from channel bulk "
                "concentration",
            ),
        )
    else:
        CONFIG.declare(
            "Dconc_x0",
            ConfigValue(
                default=None,
                description="Deviation of concentration at x= 0 from "
                "channel bulk concentration",
            ),
        )
        CONFIG.declare(
            "Dconc_x1",
            ConfigValue(
                default=None,
                description="Deviation of concentration at x= 0 from "
                "channel bulk concentration",
            ),
        )
        CONFIG.declare(
            "xflux_x0",
            ConfigValue(default=None, description="Variable for material flux at x=0"),
        )
        CONFIG.declare(
            "xflux_x1",
            ConfigValue(default=None, description="Variable for material flux at x=1"),
        )


def _create_material_boundary_conditions_if_none(unit, thin):
    tset = unit.flowsheet().config.time
    iznodes = unit.iznodes
    comps = unit.comps
    if thin:
        _create_if_none(
            unit,
            "Dconc",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / pyo.units.m**3,
        )
        _create_if_none(
            unit,
            "xflux",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / (pyo.units.s * pyo.units.m**2),
        )
    else:
        _create_if_none(
            unit,
            "Dconc_x0",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / pyo.units.m**3,
        )
        _create_if_none(
            unit,
            "Dconc_x1",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / pyo.units.m**3,
        )
        _create_if_none(
            unit,
            "xflux_x0",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / (pyo.units.s * pyo.units.m**2),
        )
        _create_if_none(
            unit,
            "xflux_x1",
            idx_set=(tset, iznodes, comps),
            units=pyo.units.mol / (pyo.units.s * pyo.units.m**2),
        )
