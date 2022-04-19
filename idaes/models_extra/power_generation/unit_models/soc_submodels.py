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
_scaling_factor_pressure = 1e-4
_scaling_factor_temperature_thermo = 1e-2
_scaling_factor_temperature_channel = 1e2
_scaling_factor_temperature_conduction = 1e2
_scaling_factor_concentration = 1
_scaling_factor_Dconcentration = 10
_scaling_factor_energy = 1e-3
_safe_log_eps = 1e-9
_safe_sqrt_eps = 1e-9
_Tmax = None


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


def _init_dof_error(blk):
    if not mstat.degrees_of_freedom(blk) == 0:
        raise InitializationError(
            f"{blk.name} encountered a nonsquare optimization problem during "
            "initialization. Check whether all parameter values are fixed and "
            "keep in mind that the SOC submodel initialization methods assume "
            "that certain degrees of freedom have been fixed by the cell "
            "initialization method."
        )


def _init_solve_block(blk, solver, log):
    _init_dof_error(blk)
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


def contour_grid_data(
    var, time, xnodes, znodes, left=None, right=None, top=None, bottom=None
):
    if left is None and right is None:
        nz = len(znodes)
    elif left is not None and right is not None:
        nz = len(znodes) + 2
    else:
        nz = len(znodes) + 1

    if top is None and bottom is None:
        nx = len(xnodes)
    elif top is not None and bottom is not None:
        nx = len(xnodes) + 2
    else:
        nx = len(xnodes) + 1

    data_z = [None] * nx
    data_x = [None] * nx
    data_w = [None] * nx
    hg = {}
    xg = {}
    zg = {}
    for it, t in enumerate(time):
        i = 0
        zg[it] = [None] * nz
        xg[it] = [None] * nz
        hg[it] = [None] * nz
        for iz, z in enumerate(znodes):
            for ix, x in enumerate(xnodes):
                data_z[ix] = z
                data_x[ix] = x
                data_w[ix] = pyo.value(var[t, ix + 1, iz + 1])
                i += 1
            zg[it][iz] = copy.copy(data_z)
            xg[it][iz] = copy.copy(data_x)
            hg[it][iz] = copy.copy(data_w)
        # Add top and/or bottom
        for iz, z in enumerate(znodes):
            if bottom is not None:
                zg[it][iz].pop()
                xg[it][iz].pop()
                hg[it][iz].pop()
                zg[it][iz].insert(0, z)
                xg[it][iz].insert(0, 0.0)
                hg[it][iz].insert(0, pyo.value(bottom[t, iz + 1]))
            if top is not None:
                zg[it][iz][-1] = z
                xg[it][iz][-1] = 1.0
                hg[it][iz][-1] = pyo.value(top[t, iz + 1])
        if left is not None:
            zg[it].pop()
            xg[it].pop()
            hg[it].pop()
            zg[it].insert(0, copy.copy(data_z))
            xg[it].insert(0, copy.copy(data_z))
            hg[it].insert(0, copy.copy(data_z))
            if bottom is not None:
                offset = 1
                zg[it][0][0] = 0.0
                xg[it][0][0] = 0.0
                hg[it][0][0] = None
            else:
                offset = 0
            if top is not None:
                offset = 1
                zg[it][0][-1] = 0.0
                xg[it][0][-1] = 1.0
                hg[it][0][-1] = None
            for ix, x in enumerate(xnodes):
                zg[it][0][ix + offset] = 0.0
                xg[it][0][ix + offset] = x
                hg[it][0][ix + offset] = pyo.value(left[t, ix + 1])
        if right is not None:
            zg[it][-1] = copy.copy(data_z)
            xg[it][-1] = copy.copy(data_z)
            hg[it][-1] = copy.copy(data_z)
            if bottom is not None:
                offset = 1
                zg[it][-1][0] = 1.0
                xg[it][-1][0] = 0.0
                hg[it][-1][0] = None
            else:
                offset = 0
            if top is not None:
                offset = 1
                zg[it][-1][-1] = 1.0
                xg[it][-1][-1] = 1.0
                hg[it][-1][-1] = None
            for ix, x in enumerate(xnodes):
                zg[it][-1][ix + offset] = 1.0
                xg[it][-1][ix + offset] = x
                hg[it][-1][ix + offset] = pyo.value(right[t, ix + 1])

    return zg, xg, hg


element_list = ["Ar", "H", "O", "C", "N", "S"]
species_list = [
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
element_dict = {
    "Ar": {"Ar": 1},
    "H": {"CH4": 4, "C2H4": 4, "C2H6": 6, "C3H8": 8, "H2": 2, "H2O": 2, "H2S": 2},
    "O": {"CO": 1, "CO2": 2, "H2O": 1, "NO": 1, "N2O": 1, "O2": 2, "SO2": 2},
    "C": {"CH4": 1, "CO": 1, "CO2": 1, "C2H4": 4, "C2H6": 2, "C3H8": 3},
    "N": {"NO": 1, "N2": 2, "N2O": 2},
    "S": {"H2S": 1, "SO2": 1},
}

for value in element_dict.values():
    for species in species_list:
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


def binary_diffusion_coefficient_expr(temperature, p, c1, c2):
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


def comp_enthalpy_expr(temperature, comp):
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


def comp_int_energy_expr(temperature, comp):
    # ideal gas enthalpy
    d = h_params[comp]
    t = temperature / 1000.0
    return comp_enthalpy_expr(temperature, comp) - _constR * temperature


def comp_entropy_expr(temperature, comp):
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


@declare_process_block_class("SocChannel")
class SocChannelData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
                **default** = useDefault.
                **Valid values:** {
                **useDefault** - get flag from parent (default = False),
                **True** - set as a dynamic model,
                **False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(domain=In([useDefault, True, False]), default=useDefault),
    )
    CONFIG.declare(
        "comp_list",
        ConfigValue(default=["H2", "H2O"], description="List of components"),
    )

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
        "include_temperature_x_thermo",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=True,
            description="Whether to consider temperature variations in "
            "x direction in thermodynamic equations",
        ),
    )
    _thermal_boundary_conditions_config(CONFIG, thin=False)
    CONFIG.declare(
        "opposite_flow",
        ConfigValue(default=False, description="If True assume velocity is negative"),
    )
    CONFIG.declare(
        "below_electrode",
        ConfigValue(
            domain=In([True, False]),
            description="Decides whether or not to create material "
            "flux terms above or below the channel.",
        ),
    )
    CONFIG.declare(
        "interpolation_scheme",
        ConfigValue(
            default=CV_Interpolation.UDS,
            description="Method used to interpolate face values",
        ),
    )

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        self.comps = pyo.Set(initialize=self.config.comp_list)
        # z coordinates for nodes and faces
        self.zfaces = pyo.Set(initialize=self.config.cv_zfaces)
        self.znodes = pyo.Set(
            initialize=[
                (self.zfaces.at(i) + self.zfaces.at(i + 1)) / 2.0
                for i in range(1, len(self.zfaces))
            ]
        )
        # This sets provide an integer index for nodes and faces
        self.izfaces = pyo.Set(initialize=range(1, len(self.zfaces) + 1))
        self.iznodes = pyo.Set(initialize=range(1, len(self.znodes) + 1))

        # Space saving aliases
        comps = self.comps
        izfaces = self.izfaces
        iznodes = self.iznodes
        zfaces = self.zfaces
        znodes = self.znodes
        include_temp_x_thermo = self.config.include_temperature_x_thermo

        _create_if_none(self, "length_z", idx_set=None, units=pyo.units.m)
        _create_if_none(self, "length_y", idx_set=None, units=pyo.units.m)

        _create_thermal_boundary_conditions_if_none(self, thin=False)

        if not self.config.below_electrode:
            self.Dconc_x0 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Deviation of concentration at electrode surface from that "
                "in the bulk channel",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3,
            )
            self.xflux_x0 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Material flux from electrode surface to channel "
                "(positive is in)",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**2 / pyo.units.s,
            )

            @self.Expression(tset, iznodes, comps)
            def Dconc_x1(b, t, iz, j):
                return 0

            @self.Expression(tset, iznodes, comps)
            def xflux_x1(b, t, iz, j):
                return 0

        else:
            self.Dconc_x1 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Deviation of concentration at electrode surface from that "
                "in the bulk channel",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3,
            )
            self.xflux_x1 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Material flux from channel to electrode surface "
                "(positive is out)",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**2 / pyo.units.s,
            )

            @self.Expression(tset, iznodes, comps)
            def Dconc_x0(b, t, iz, j):
                return 0

            @self.Expression(tset, iznodes, comps)
            def xflux_x0(b, t, iz, j):
                return 0

        # Channel thickness AKA length in the x direction is specific to the
        # channel so local variable here is the only option
        self.length_x = pyo.Var(
            doc="Thickness from interconnect to electrode (x-direction)",
            units=pyo.units.m,
        )
        self.heat_transfer_coefficient = pyo.Var(
            tset,
            iznodes,
            doc="Local channel heat transfer coefficient",
            initialize=500,
            units=pyo.units.J / pyo.units.m**2 / pyo.units.s / pyo.units.K,
        )
        self.flow_mol = pyo.Var(
            tset,
            iznodes,
            doc="Molar flow in the z-direction through faces",
            units=pyo.units.mol / pyo.units.s,
            bounds=(-1e-9, None),
        )
        self.conc = pyo.Var(
            tset,
            iznodes,
            comps,
            doc="Component concentration at node centers",
            units=pyo.units.mol / pyo.units.m**3,
            bounds=(0, None),
        )
        self.Dtemp = pyo.Var(
            tset,
            iznodes,
            doc="Deviation of temperature at node centers " "from temperature_z",
            units=pyo.units.K,
            bounds=(-1000, 1000),
        )
        self.enth_mol = pyo.Var(
            tset,
            iznodes,
            doc="Molar enthalpy at node centers",
            units=pyo.units.J / pyo.units.mol,
        )
        if self.config.has_holdup:
            self.int_energy_mol = pyo.Var(
                tset,
                iznodes,
                doc="Molar internal energy at node centers",
                units=pyo.units.J / pyo.units.mol,
            )
            self.int_energy_density = pyo.Var(
                tset,
                iznodes,
                doc="Molar internal energy density at node centers",
                units=pyo.units.J / pyo.units.m**3,
            )
        self.velocity = pyo.Var(
            tset,
            iznodes,
            doc="Fluid velocity at node centers",
            units=pyo.units.m / pyo.units.s,
        )
        self.pressure = pyo.Var(
            tset,
            iznodes,
            doc="Pressure at node centers",
            units=pyo.units.Pa,
            bounds=(0, None),
        )
        self.mole_frac_comp = pyo.Var(
            tset,
            iznodes,
            comps,
            doc="Component mole fraction at node centers",
            bounds=(0, None),
            units=pyo.units.dimensionless,
        )
        self.flow_mol_inlet = pyo.Var(
            tset,
            doc="Inlet face molar flow rate",
            bounds=(0, None),
            units=pyo.units.mol / pyo.units.s,
        )
        self.pressure_inlet = pyo.Var(
            tset, doc="Inlet pressure", bounds=(0, None), units=pyo.units.Pa
        )
        self.temperature_inlet = pyo.Var(
            tset, doc="Inlet temperature", bounds=(600, _Tmax), units=pyo.units.K
        )
        self.temperature_outlet = pyo.Var(
            tset, doc="Outlet temperature", bounds=(600, _Tmax), units=pyo.units.K
        )
        self.mole_frac_comp_inlet = pyo.Var(
            tset,
            comps,
            doc="Inlet compoent mole fractions",
            bounds=(0, 1),
            units=pyo.units.dimensionless,
        )

        # Add time derivative varaible if steady state use const 0.
        if dynamic:
            self.dcdt = DerivativeVar(
                self.conc,
                wrt=tset,
                initialize=0,
                doc="Component concentration time derivative",
            )
        else:
            self.dcdt = pyo.Param(
                tset,
                iznodes,
                comps,
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3 / pyo.units.s,
            )

        # Add time derivative varaible if steady state use const 0.
        if dynamic:
            self.dcedt = DerivativeVar(
                self.int_energy_density,
                wrt=tset,
                initialize=0,
                doc="Internal energy density time derivative",
            )
        else:
            self.dcedt = pyo.Param(
                tset,
                iznodes,
                initialize=0,
                units=pyo.units.J / pyo.units.m**3 / pyo.units.s,
            )

        @self.Expression()
        def flow_area(b):
            return b.length_x[None] * b.length_y[None]

        @self.Expression(iznodes)
        def dz(b, iz):
            return b.zfaces.at(iz + 1) - b.zfaces.at(iz)

        @self.Expression(iznodes)
        def node_volume(b, iz):
            return b.length_x[None] * b.length_y[None] * b.length_z[None] * b.dz[iz]

        @self.Expression(iznodes)
        def xface_area(b, iz):
            return b.length_z[None] * b.length_y[None] * b.dz[iz]

        @self.Expression(tset, iznodes)
        def temperature(b, t, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.Dtemp[t, iz]
            else:
                return b.temperature_z[t, iz]

        @self.Expression(tset, iznodes)
        def volume_molar(b, t, iz):
            return _constR * b.temperature[t, iz] / b.pressure[t, iz]

        @self.Expression(tset)
        def volume_molar_inlet(b, t):
            return _constR * b.temperature_inlet[t] / b.pressure_inlet[t]

        # TODO maybe replace with variable-constraint pair?
        @self.Expression(tset)
        def enth_mol_inlet(b, t):
            return sum(
                comp_enthalpy_expr(b.temperature_inlet[t], i)
                * b.mole_frac_comp_inlet[t, i]
                for i in comps
            )

        # TODO maybe replace with variable-constraint pair?
        @self.Expression(tset, iznodes, comps)
        def diff_eff_coeff(b, t, iz, i):
            T = b.temperature[t, iz]
            P = b.pressure[t, iz]
            x = b.mole_frac_comp
            bfun = binary_diffusion_coefficient_expr
            return (1.0 - x[t, iz, i]) / sum(
                x[t, iz, j] / bfun(T, P, i, j) for j in comps if i != j
            )

        @self.Expression(tset, iznodes, comps)
        def mass_transfer_coeff(b, t, iz, i):
            # Quick and dirty approximation based on Ficks law through a thin
            # film of length L_x/2. For small concentration gradients
            # this will (hopefully) be enough
            return 2 * b.diff_eff_coeff[t, iz, i] / b.length_x

        if not self.config.below_electrode:

            @self.Constraint(tset, iznodes, comps)
            def xflux_x0_eqn(b, t, iz, i):
                return (
                    b.xflux_x0[t, iz, i]
                    == b.mass_transfer_coeff[t, iz, i] * b.Dconc_x0[t, iz, i]
                )

        else:

            @self.Constraint(tset, iznodes, comps)
            def xflux_x1_eqn(b, t, iz, i):
                return (
                    b.xflux_x1[t, iz, i]
                    == -b.mass_transfer_coeff[t, iz, i] * b.Dconc_x1[t, iz, i]
                )

        @self.Constraint(tset, iznodes)
        def flow_mol_eqn(b, t, iz):
            # either way the flow goes, want the flow rate to be positive, but
            # in the opposite flow cases want flux and velocity to be negative
            if self.config.opposite_flow:
                return (
                    b.flow_mol[t, iz]
                    == -b.flow_area * b.velocity[t, iz] / b.volume_molar[t, iz]
                )
            return (
                b.flow_mol[t, iz]
                == b.flow_area * b.velocity[t, iz] / b.volume_molar[t, iz]
            )

        @self.Constraint(tset, iznodes)
        def constant_pressure_eqn(b, t, iz):
            return b.pressure[t, iz] == b.pressure_inlet[t]

        @self.Constraint(tset, iznodes, comps)
        def conc_eqn(b, t, iz, i):
            return (
                b.conc[t, iz, i] * b.temperature[t, iz] * _constR
                == b.pressure[t, iz] * b.mole_frac_comp[t, iz, i]
            )

        @self.Constraint(tset, iznodes)
        def enth_mol_eqn(b, t, iz):
            return b.enth_mol[t, iz] == sum(
                comp_enthalpy_expr(b.temperature[t, iz], i) * b.mole_frac_comp[t, iz, i]
                for i in comps
            )

        if self.config.has_holdup:

            @self.Constraint(tset, iznodes)
            def int_energy_mol_eqn(b, t, iz):
                return b.int_energy_mol[t, iz] == sum(
                    comp_int_energy_expr(b.temperature[t, iz], i)
                    * b.mole_frac_comp[t, iz, i]
                    for i in comps
                )

            @self.Constraint(tset, iznodes)
            def int_energy_density_eqn(b, t, iz):
                return (
                    b.int_energy_density[t, iz]
                    == b.int_energy_mol[t, iz] / b.volume_molar[t, iz]
                )

        @self.Constraint(tset, iznodes)
        def mole_frac_eqn(b, t, iz):
            return 1 == sum(b.mole_frac_comp[t, iz, i] for i in comps)

        @self.Expression(tset, comps)
        def flow_mol_comp_inlet(b, t, i):
            return b.flow_mol_inlet[t] * b.mole_frac_comp_inlet[t, i]

        @self.Expression(tset, comps)
        def zflux_inlet(b, t, i):
            # either way the flow goes, want the flow rate to be positive, but
            # in the opposite flow cases want flux and velocity to be negative
            if self.config.opposite_flow:
                return -b.flow_mol_inlet[t] / b.flow_area * b.mole_frac_comp_inlet[t, i]
            return b.flow_mol_inlet[t] / b.flow_area * b.mole_frac_comp_inlet[t, i]

        @self.Expression(tset)
        def zflux_enth_inlet(b, t):
            # either way the flow goes, want the flow rate to be positive, but
            # in the opposite flow cases want flux and velocity to be negative
            if self.config.opposite_flow:
                return -b.flow_mol_inlet[t] / b.flow_area * b.enth_mol_inlet[t]
            return b.flow_mol_inlet[t] / b.flow_area * b.enth_mol_inlet[t]

        @self.Expression(tset, izfaces, comps)
        def zflux(b, t, iz, i):
            return _interpolate_channel(
                iz=iz,
                ifaces=izfaces,
                nodes=znodes,
                faces=zfaces,
                phi_func=lambda iface: b.velocity[t, iface] * b.conc[t, iface, i],
                phi_inlet=b.zflux_inlet[t, i],
                method=self.config.interpolation_scheme,
                opposite_flow=self.config.opposite_flow,
            )

        @self.Expression(tset, izfaces)
        def zflux_enth(b, t, iz):
            return _interpolate_channel(
                iz=iz,
                ifaces=izfaces,
                nodes=znodes,
                faces=zfaces,
                phi_func=lambda iface: b.velocity[t, iface]
                / b.volume_molar[t, iface]
                * b.enth_mol[t, iface],
                phi_inlet=b.zflux_enth_inlet[t],
                method=self.config.interpolation_scheme,
                opposite_flow=self.config.opposite_flow,
            )

        @self.Expression(tset, izfaces)
        def pressure_face(b, t, iz):
            # Although I'm currently assuming no pressure drop in the channel
            # and don't have a momentum balance, this will let me estimate the
            # outlet pressure in a way that will let me add in a momentum balance
            # later
            return _interpolate_channel(
                iz=iz,
                ifaces=izfaces,
                nodes=znodes,
                faces=zfaces,
                phi_func=lambda iface: b.pressure[t, iface],
                phi_inlet=b.pressure_inlet[t],
                method=self.config.interpolation_scheme,
                opposite_flow=self.config.opposite_flow,
            )

        @self.Constraint(tset, iznodes, comps)
        def material_balance_eqn(b, t, iz, i):
            # if t == tset.first() and dynamic:
            #     return pyo.Constraint.Skip
            return b.dcdt[t, iz, i] * b.node_volume[iz] == b.flow_area * (
                b.zflux[t, iz, i] - b.zflux[t, iz + 1, i]
            ) + b.xface_area[iz] * (b.xflux_x0[t, iz, i] - b.xflux_x1[t, iz, i])

        if dynamic:
            self.material_balance_eqn[tset.first(), :, :].deactivate()

        @self.Constraint(tset, iznodes)
        def energy_balance_eqn(b, t, iz):
            return (
                b.dcedt[t, iz] * b.node_volume[iz]
                == b.flow_area * (b.zflux_enth[t, iz] - b.zflux_enth[t, iz + 1])
                + b.xface_area[iz]
                * sum(
                    (b.xflux_x0[t, iz, i] - b.xflux_x1[t, iz, i])
                    * comp_enthalpy_expr(b.temperature_x1[t, iz], i)
                    for i in comps
                )
                + b.qflux_x0[t, iz] * b.xface_area[iz]
                - b.qflux_x1[t, iz] * b.xface_area[iz]
            )

        if dynamic:
            self.energy_balance_eqn[tset.first(), :].deactivate()

        @self.Constraint(tset, iznodes)
        def temperature_x0_eqn(b, t, iz):
            return self.qflux_x0[t, iz] == self.heat_transfer_coefficient[t, iz] * (
                self.Dtemp_x0[t, iz] - self.Dtemp[t, iz]
            )

        @self.Constraint(tset, iznodes)
        def temperature_x1_eqn(b, t, iz):
            return self.qflux_x1[t, iz] == self.heat_transfer_coefficient[t, iz] * (
                self.Dtemp[t, iz] - self.Dtemp_x1[t, iz]
            )

        # For convenience define outlet expressions
        if self.config.opposite_flow:
            izfout = self.izfout = izfaces.first()
            iznout = self.iznout = iznodes.first()
        else:
            izfout = self.izfout = izfaces.last()
            iznout = self.iznout = iznodes.last()

        @self.Expression(tset, comps)
        def flow_mol_comp_outlet(b, t, i):
            if self.config.opposite_flow:
                return -b.zflux[t, izfout, i] * b.flow_area
            else:
                return b.zflux[t, izfout, i] * b.flow_area

        # @self.Expression(tset)
        def rule_flow_mol_outlet(b, t):
            return sum(b.flow_mol_comp_outlet[t, i] for i in comps)

        # @self.Expression(tset)
        def rule_pressure_outlet(b, t):
            return b.pressure_face[t, izfout]

        # @self.Expression(tset, comps)
        def rule_mole_frac_comp_outlet(b, t, i):
            return b.flow_mol_comp_outlet[t, i] / b.flow_mol_outlet[t]

        self.flow_mol_outlet = VarLikeExpression(tset, rule=rule_flow_mol_outlet)
        self.pressure_outlet = VarLikeExpression(tset, rule=rule_pressure_outlet)
        self.mole_frac_comp_outlet = VarLikeExpression(
            tset, comps, rule=rule_mole_frac_comp_outlet
        )

        @self.Expression(tset)
        def enth_mol_outlet(b, t):
            if self.config.opposite_flow:
                return -b.zflux_enth[t, izfout] * b.flow_area / b.flow_mol_outlet[t]
            else:
                return b.zflux_enth[t, izfout] * b.flow_area / b.flow_mol_outlet[t]

        # know enthalpy need a constraint to back calculate temperature
        @self.Constraint(tset)
        def temperature_outlet_eqn(b, t):
            return b.enth_mol_outlet[t] == sum(
                comp_enthalpy_expr(b.temperature_outlet[t], i)
                * b.mole_frac_comp_outlet[t, i]
                for i in comps
            )

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        xflux_guess=None,
        qflux_x1_guess=None,
        qflux_x0_guess=None,
        velocity_guess=None,
    ):
        # At present, this method does not fix inlet variables because they are
        # fixed at the cell level instead.
        # TODO Add ports to submodel instead?

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        tset = self.flowsheet().config.time
        t0 = tset.first()

        for t in tset:
            _set_if_unfixed(self.temperature_outlet[t], self.temperature_inlet[t])
            for iz in self.iznodes:
                _set_if_unfixed(self.temperature_z[t, iz], self.temperature_inlet[t])
                _set_if_unfixed(self.Dtemp_x0[t, iz], 0)
                _set_if_unfixed(self.Dtemp[t, iz], 0)
                _set_if_unfixed(self.Dtemp_x1[t, iz], 0)
                _set_if_unfixed(self.pressure[t, iz], self.pressure_inlet[t])
                _set_if_unfixed(self.flow_mol[t, iz], self.flow_mol_inlet[t])
                for i in self.config.comp_list:
                    _set_if_unfixed(
                        self.mole_frac_comp[t, iz, i], self.mole_frac_comp_inlet[t, i]
                    )
                    _set_if_unfixed(
                        self.conc[t, iz, i],
                        self.pressure[t, iz]
                        * self.mole_frac_comp[t, iz, i]
                        / self.temperature[t, iz]
                        / _constR,
                    )
                _set_if_unfixed(
                    self.velocity[t, iz],
                    self.flow_mol[t, iz] / self.flow_area * self.volume_molar[t, iz],
                )
                _set_if_unfixed(
                    self.enth_mol[t, iz],
                    sum(
                        comp_enthalpy_expr(self.temperature[t, iz], i)
                        * self.mole_frac_comp[t, iz, i]
                        for i in self.config.comp_list
                    ),
                )
                if self.config.has_holdup:
                    _set_if_unfixed(
                        self.int_energy_mol[t, iz],
                        sum(
                            comp_int_energy_expr(self.temperature[t, iz], i)
                            * self.mole_frac_comp[t, iz, i]
                            for i in self.config.comp_list
                        ),
                    )
                    _set_if_unfixed(
                        self.int_energy_density[t, iz],
                        self.int_energy_mol[t, iz] / self.volume_molar[t, iz],
                    )
        slvr = get_solver(solver, optarg)
        _init_solve_block(self, slvr, solve_log)

    def calculate_scaling_factors(self):
        pass

    def model_check(self, steady_state=True):
        if not steady_state:
            # Mass and energy conservation equations steady state only at present
            return

        comp_set = set(self.comps)
        elements_present = set()

        for element in element_list:
            include_element = False
            for species in species_list:
                # Floating point equality take warning!
                if species in comp_set and element_dict[element][species] != 0:
                    include_element = True
            if include_element:
                elements_present.add(element)

        for t in self.flowsheet().config.time:
            for element in element_list:
                if element not in elements_present:
                    continue
                sum_in = sum(
                    element_dict[element][j] * self.flow_mol_comp_inlet[t, j]
                    for j in self.comps
                )
                sum_out = sum(
                    element_dict[element][j] * self.flow_mol_comp_outlet[t, j]
                    for j in self.comps
                )
                for iz in self.iznodes:
                    sum_in += sum(
                        element_dict[element][j]
                        * self.xflux_x0[t, iz, j]
                        * self.xface_area[iz]
                        for j in self.comps
                    )
                    sum_out += sum(
                        element_dict[element][j]
                        * self.xflux_x1[t, iz, j]
                        * self.xface_area[iz]
                        for j in self.comps
                    )
                normal = max(
                    pyo.value(sum_in), pyo.value(sum_out), 1e-9
                )  # FIXME justify this number
                fraction_change = pyo.value((sum_out - sum_in) / normal)
                if abs(fraction_change) > 1e-5:
                    raise RuntimeError(
                        f"{element} is not being conserved in {self.name}; "
                        f"fractional change {fraction_change}."
                    )
            enth_in = self.enth_mol_inlet[t] * self.flow_mol_inlet[t]
            enth_out = self.enth_mol_outlet[t] * self.flow_mol_outlet[t]

            for iz in self.iznodes:
                enth_in += self.xface_area[iz] * (
                    self.qflux_x0[t, iz]
                    + sum(
                        comp_enthalpy_expr(self.temperature_x0[t, iz], j)
                        * self.xflux_x0[t, iz, j]
                        for j in self.comps
                    )
                )
                enth_out += self.xface_area[iz] * (
                    self.qflux_x1[t, iz]
                    + sum(
                        comp_enthalpy_expr(self.temperature_x1[t, iz], j)
                        * self.xflux_x1[t, iz, j]
                        for j in self.comps
                    )
                )

            normal = max(
                pyo.value(abs(enth_in)), pyo.value(abs(enth_out)), 1e-4
            )  # FIXME justify this number
            fraction_change = pyo.value((enth_out - enth_in) / normal)
            if abs(fraction_change) > 3e-3:
                raise RuntimeError(
                    f"Energy is not being conserved in {self.name}; "
                    f"fractional change {fraction_change}"
                )

    def recursive_scaling(self):
        gsf = iscale.get_scaling_factor
        ssf = _set_scaling_factor_if_none
        sgsf = _set_and_get_scaling_factor
        sdf = _set_default_factor
        cst = lambda c, s: iscale.constraint_scaling_transform(c, s, overwrite=False)
        sR = 1e-1  # Scaling factor for R
        # sD = 5e3 # Heuristic scaling factor for diffusion coefficient
        sD = 1e4
        sy_def = 10  # Mole frac comp scaling
        sh = 1e-2  # Heat xfer coeff
        # sh = 1
        sH = 1e-4  # Enthalpy/int energy
        sLx = sgsf(self.length_x, 1 / self.length_x.value)
        sLy = 1 / self.length_y[None].value
        sLz = len(self.iznodes) / self.length_z[None].value

        for t in self.flowsheet().time:
            sT = sgsf(self.temperature_inlet[t], 1e-2)
            ssf(self.temperature_outlet[t], sT)
            sP = sgsf(self.pressure_inlet[t], 1e-4)

            s_flow_mol = sgsf(self.flow_mol_inlet[t], 1e3)
            sy_in = {}
            for j in self.comps:
                sy_in[j] = sgsf(self.mole_frac_comp_inlet[t, j], sy_def)

            for iz in self.iznodes:
                # These should have been scaled by the cell-level method, so
                # notify the user if they're using a standalone channel
                # and forgot to scale these
                if not self.temperature_z[t, iz].is_reference():
                    gsf(self.temperature_z[t, iz], warning=True)
                gsf(self.qflux_x0[t, iz], warning=True)
                gsf(self.qflux_x1[t, iz], warning=True)

                s_flow_mol = sgsf(self.flow_mol[t, iz], s_flow_mol)
                sT = sgsf(self.temperature_z[t, iz], sT)
                sP = sgsf(self.pressure[t, iz], sP)
                cst(self.constant_pressure_eqn[t, iz], sP)

                cst(self.flow_mol_eqn[t, iz], s_flow_mol)
                sV = sR * sT / sP
                ssf(self.velocity[t, iz], sV * s_flow_mol / (sLx * sLy))

                sH = sgsf(self.enth_mol[t, iz], sH)
                cst(self.enth_mol_eqn[t, iz], sH)
                cst(self.energy_balance_eqn[t, iz], sH * s_flow_mol)

                if self.config.has_holdup:
                    sU = sgsf(self.int_energy_mol[t, iz], sH)
                    cst(self.int_energy_mol_eqn[t, iz], sU)

                    s_rho_U = sgsf(self.int_energy_density[t, iz], sU / sV)
                    cst(self.int_energy_density_eqn[t, iz], s_rho_U)

                sq0 = sgsf(self.qflux_x0[t, iz], 1e-2)
                cst(self.temperature_x0_eqn[t, iz], sq0)
                sq1 = sgsf(self.qflux_x1[t, iz], 1e-2)
                cst(self.temperature_x1_eqn[t, iz], sq1)
                sq = min(sq0, sq1)

                sDT = sgsf(self.Dtemp[t, iz], sq / sh)
                ssf(self.Dtemp_x0[t, iz], sDT)
                ssf(self.Dtemp_x1[t, iz], sDT)

                # Pointless other than making a record that this equation
                # is well-scaled by default
                cst(self.mole_frac_eqn[t, iz], 1)

                for j in self.comps:
                    # These should have been scaled by the cell-level method, so
                    # notify the user if they're using a standalone channel
                    # and forgot to scale these
                    if self.config.below_electrode:
                        gsf(self.xflux_x1[t, iz, j], warning=True)
                    else:
                        gsf(self.xflux_x0[t, iz, j], warning=True)

                    sy = sgsf(self.mole_frac_comp[t, iz, j], sy_in[j])
                    cst(self.material_balance_eqn[t, iz, j], s_flow_mol * sy)

                    ssf(self.conc[t, iz, j], sy * sP / (sR * sT))
                    cst(self.conc_eqn[t, iz, j], sy * sP)

                    if hasattr(self, "xflux_x0_eqn"):
                        sXflux = gsf(
                            self.xflux_x0[t, iz, j], default=1e-1, warning=True
                        )
                        cst(self.xflux_x0_eqn[t, iz, j], sXflux)
                        ssf(self.Dconc_x0[t, iz, j], sLx * sXflux / sD)
                    if hasattr(self, "xflux_x1_eqn"):
                        sXflux = gsf(
                            self.xflux_x1[t, iz, j], default=1e-1, warning=True
                        )
                        cst(self.xflux_x1_eqn[t, iz, j], sXflux)
                        ssf(self.Dconc_x1[t, iz, j], sLx * sXflux / sD)

            cst(self.temperature_outlet_eqn[t], sH)


@declare_process_block_class("SocElectrode")
class SocElectrodeData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
                    **default** = useDefault.
                    **Valid values:** {
                    **useDefault** - get flag from parent (default = False),
                    **True** - set as a dynamic model,
                    **False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(domain=In([useDefault, True, False]), default=useDefault),
    )
    CONFIG.declare(
        "cv_xfaces",
        ConfigValue(
            description="CV x-boundary set, should start with 0 and end with 1."
        ),
    )
    CONFIG.declare(
        "comp_list",
        ConfigValue(default=["H2", "H2O"], description="List of components"),
    )
    CONFIG.declare(
        "conc_ref",
        ConfigValue(
            default=None,
            description="Variable for the component concentration in bulk channel ",
        ),
    )
    CONFIG.declare(
        "dconc_refdt",
        ConfigValue(
            default=None,
            description="Variable for time derivative of the "
            "component concentration in the bulk channel",
        ),
    )
    _submodel_boilerplate_config(CONFIG)
    _thermal_boundary_conditions_config(CONFIG, thin=False)
    _material_boundary_conditions_config(CONFIG, thin=False)

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        self.comps = pyo.Set(initialize=self.config.comp_list)
        # z coordinates for nodes and faces
        self.zfaces = pyo.Set(initialize=self.config.cv_zfaces)
        self.znodes = pyo.Set(
            initialize=[
                (self.zfaces.at(i) + self.zfaces.at(i + 1)) / 2.0
                for i in range(1, len(self.zfaces))
            ]
        )
        self.xfaces = pyo.Set(initialize=self.config.cv_xfaces)
        self.xnodes = pyo.Set(
            initialize=[
                (self.xfaces.at(i) + self.xfaces.at(i + 1)) / 2.0
                for i in range(1, len(self.xfaces))
            ]
        )
        # This sets provide an integer index for nodes and faces
        self.izfaces = pyo.Set(initialize=range(1, len(self.zfaces) + 1))
        self.iznodes = pyo.Set(initialize=range(1, len(self.znodes) + 1))
        self.ixfaces = pyo.Set(initialize=range(1, len(self.xfaces) + 1))
        self.ixnodes = pyo.Set(initialize=range(1, len(self.xnodes) + 1))

        # Space saving aliases
        comps = self.comps
        izfaces = self.izfaces
        iznodes = self.iznodes
        ixfaces = self.ixfaces
        ixnodes = self.ixnodes
        zfaces = self.zfaces
        znodes = self.znodes
        xfaces = self.xfaces
        xnodes = self.xnodes
        include_temp_x_thermo = self.config.include_temperature_x_thermo

        # Electrode thickness AKA length in the x direction is specific to the
        # electrode so local variable here is the only option
        self.length_x = pyo.Var(
            doc="Thickness of the electrode (x-direction)",
            units=pyo.units.m,
        )

        if self.config.dynamic:
            # If we're dynamic, the user needs to either provide both the
            # reference concentration and its derivative, or provide neither.
            assert (
                self.config.conc_ref is None and self.config.dconc_refdt is None
            ) or (
                self.config.conc_ref is not None and self.config.dconc_refdt is not None
            )

        if self.config.conc_ref is None:
            self.conc_ref = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Concentration of components in the channel bulk",
                initialize=0.0,
                units=pyo.units.mol / pyo.units.m**3,
            )
            if self.config.dynamic:
                self.dconc_refdt = DerivativeVar(
                    self.conc_ref,
                    wrt=tset,
                    doc="Derivative of concentration of components in the channel bulk",
                    initialize=0.0,
                    units=pyo.units.mol / (pyo.units.s * pyo.units.m**3),
                )
            else:
                self.dconc_refdt = pyo.Param(
                    tset,
                    iznodes,
                    comps,
                    doc="Derivative of concentration of components in the channel bulk",
                    initialize=0.0,
                    units=pyo.units.mol / (pyo.units.s * pyo.units.m**3),
                )
        else:
            self.conc_ref = pyo.Reference(self.config.conc_ref)
            if self.config.dconc_refdt.ctype == DerivativeVar:
                self.dconc_refdt = pyo.Reference(self.config.dconc_refdt, ctype=pyo.Var)
            else:
                self.dconc_refdt = pyo.Reference(self.config.dconc_refdt)
        _submodel_boilerplate_create_if_none(self)
        _create_thermal_boundary_conditions_if_none(self, thin=False)
        _create_material_boundary_conditions_if_none(self, thin=False)

        self.porosity = pyo.Var(
            initialize=0.50, doc="Electrode porosity", units=pyo.units.dimensionless
        )
        self.tortuosity = pyo.Var(
            initialize=2.0, doc="Electrode tortuosity", units=pyo.units.dimensionless
        )

        self.Dtemp = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            doc="Deviation of temperature at node centers " "from temperature_z",
            units=pyo.units.K,
            bounds=(-1000, 1000),
        )
        self.Dconc = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            comps,
            doc="Deviation of component concentration at node centers " "from conc_ref",
            units=pyo.units.mol / pyo.units.m**3,
            bounds=(-100, 100),
        )
        self.enth_mol = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            doc="Molar enthalpy at node centers",
            units=pyo.units.J / pyo.units.mol,
        )
        if self.config.has_holdup:
            self.int_energy_mol = pyo.Var(
                tset,
                ixnodes,
                iznodes,
                doc="Fluid molar internal energy at node centers",
                units=pyo.units.J / pyo.units.mol,
            )
            self.int_energy_density = pyo.Var(
                tset,
                ixnodes,
                iznodes,
                doc="Fluid molar internal energy density at node centers",
                units=pyo.units.J / pyo.units.m**3,
            )
            self.int_energy_density_solid = pyo.Var(
                tset,
                ixnodes,
                iznodes,
                doc="Internal energy density of solid electrode",
                units=pyo.units.J / pyo.units.m**3,
            )
        self.pressure = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            doc="Pressure at node centers",
            units=pyo.units.Pa,
            bounds=(0, None),
        )

        # Assume the the electrode gas phase and solid are same temp
        self.mole_frac_comp = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            comps,
            doc="Component mole fraction at node centers",
            initialize=1 / len(comps),
            bounds=(0, 1),
        )

        self.resistivity_log_preexponential_factor = pyo.Var(
            doc="Logarithm of resistivity preexponential factor " "in units of ohm*m",
            units=pyo.units.dimensionless,
        )
        self.resistivity_thermal_exponent_dividend = pyo.Var(
            doc="Parameter divided by temperature in exponential", units=pyo.units.K
        )

        # Parameters
        self.solid_heat_capacity = pyo.Var()
        self.solid_density = pyo.Var()
        self.solid_thermal_conductivity = pyo.Var()

        # Add time derivative varaible if steady state use const 0.
        if dynamic:
            self.dDconcdt = DerivativeVar(
                self.Dconc,
                wrt=tset,
                initialize=0,
                doc="Component concentration time derivative in deviation " "variable",
            )
        else:
            self.dDconcdt = pyo.Param(
                tset,
                ixnodes,
                iznodes,
                comps,
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3 / pyo.units.s,
            )
        # Add time derivative varaible if steady state use const 0.
        if dynamic:
            self.dcedt = DerivativeVar(
                self.int_energy_density,
                wrt=tset,
                initialize=0,
                doc="Internal energy density time derivative",
            )
        else:
            self.dcedt = pyo.Param(
                tset,
                ixnodes,
                iznodes,
                initialize=0,
                units=pyo.units.W / pyo.units.m**3,
            )
        # Add time derivative varaible if steady state use const 0.
        if dynamic:
            self.dcedt_solid = DerivativeVar(
                self.int_energy_density_solid,
                wrt=tset,
                initialize=0,
                doc="Internal energy density time derivative",
            )
        else:
            self.dcedt_solid = pyo.Param(
                tset,
                ixnodes,
                iznodes,
                initialize=0,
                units=pyo.units.W / pyo.units.m**3,
            )

        @self.Expression(iznodes)
        def dz(b, iz):
            return b.zfaces.at(iz + 1) - b.zfaces.at(iz)

        @self.Expression(ixnodes)
        def dx(b, ix):
            return b.zfaces.at(ix + 1) - b.zfaces.at(ix)

        @self.Expression(ixnodes, iznodes)
        def node_volume(b, ix, iz):
            return (
                b.length_x[None]
                * b.length_y[None]
                * b.length_z[None]
                * b.dz[iz]
                * b.dx[ix]
            )

        @self.Expression(ixnodes)
        def zface_area(b, ix):
            return b.length_y[None] * b.length_x[None] * b.dx[ix]

        @self.Expression(iznodes)
        def xface_area(b, iz):
            return b.length_y[None] * b.length_z[None] * b.dz[iz]

        @self.Expression(tset, ixnodes, iznodes)
        def temperature(b, t, ix, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.Dtemp[t, ix, iz]
            else:
                return b.temperature_z[t, iz]

        @self.Expression(tset, iznodes, comps)
        def conc_x0(b, t, iz, j):
            return b.conc_ref[t, iz, j] + b.Dconc_x0[t, iz, j]

        @self.Expression(tset, ixnodes, iznodes, comps)
        def conc(b, t, ix, iz, j):
            return b.conc_ref[t, iz, j] + b.Dconc[t, ix, iz, j]

        @self.Expression(tset, iznodes, comps)
        def conc_x1(b, t, iz, j):
            return b.conc_ref[t, iz, j] + b.Dconc_x1[t, iz, j]

        @self.Expression(tset, ixnodes, iznodes, comps)
        def dcdt(b, t, ix, iz, j):
            return b.dconc_refdt[t, iz, j] + b.dDconcdt[t, ix, iz, j]

        @self.Expression(tset, ixnodes, iznodes)
        def volume_molar(b, t, ix, iz):
            return _constR * b.temperature[t, ix, iz] / b.pressure[t, ix, iz]

        @self.Constraint(tset, ixnodes, iznodes, comps)
        def conc_eqn(b, t, ix, iz, i):
            return (
                b.conc[t, ix, iz, i] * b.temperature[t, ix, iz] * _constR
                == b.pressure[t, ix, iz] * b.mole_frac_comp[t, ix, iz, i]
            )

        @self.Constraint(tset, ixnodes, iznodes)
        def enth_mol_eqn(b, t, ix, iz):
            return b.enth_mol[t, ix, iz] == sum(
                comp_enthalpy_expr(b.temperature[t, ix, iz], i)
                * b.mole_frac_comp[t, ix, iz, i]
                for i in comps
            )

        if self.config.has_holdup:
            # For the vapor phase
            @self.Constraint(tset, ixnodes, iznodes)
            def int_energy_mol_eqn(b, t, ix, iz):
                return b.int_energy_mol[t, ix, iz] == sum(
                    comp_int_energy_expr(b.temperature[t, ix, iz], i)
                    * b.mole_frac_comp[t, ix, iz, i]
                    for i in comps
                )

            @self.Constraint(tset, ixnodes, iznodes)
            def int_energy_density_eqn(b, t, ix, iz):
                return (
                    b.int_energy_density[t, ix, iz]
                    == b.int_energy_mol[t, ix, iz] / b.volume_molar[t, ix, iz]
                )

            @self.Constraint(tset, ixnodes, iznodes)
            def int_energy_density_solid_eqn(b, t, ix, iz):
                return b.int_energy_density_solid[
                    t, ix, iz
                ] == b.solid_heat_capacity * b.solid_density * (
                    b.temperature[t, ix, iz] - 1000 * pyo.units.K
                )

        # @self.Constraint(tset, ixnodes, iznodes)
        # def int_energy_density_eqn(b, t, ix, iz):
        #     return (
        #         b.int_energy_density[t, ix, iz]
        #         == b.porosity * b.int_energy_mol[t, ix, iz] / b.volume_molar[t, ix, iz]
        #         + (1 - b.porosity) * b.int_energy_density_solid[t, ix, iz]
        #     )

        @self.Constraint(tset, ixnodes, iznodes)
        def mole_frac_eqn(b, t, ix, iz):
            return 1 == sum(b.mole_frac_comp[t, ix, iz, i] for i in comps)

        @self.Expression(tset, ixnodes, iznodes, comps)
        def diff_eff_coeff(b, t, ix, iz, i):
            T = b.temperature[t, ix, iz]
            P = b.pressure[t, ix, iz]
            x = b.mole_frac_comp
            bfun = binary_diffusion_coefficient_expr
            return (
                b.porosity
                / b.tortuosity
                * (1.0 - x[t, ix, iz, i])
                / sum(x[t, ix, iz, j] / bfun(T, P, i, j) for j in comps if i != j)
            )

        @self.Expression(tset, ixfaces, iznodes, comps)
        def dcdx(b, t, ix, iz, i):
            return _interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=xnodes,
                faces=xfaces,
                phi_func=lambda ixf: b.Dconc[t, ixf, iz, i] / b.length_x,
                phi_bound_0=(b.Dconc[t, ixnodes.first(), iz, i] - b.Dconc_x0[t, iz, i])
                / (xnodes.first() - xfaces.first())
                / b.length_x,
                phi_bound_1=(b.Dconc_x1[t, iz, i] - b.Dconc[t, ixnodes.last(), iz, i])
                / (xfaces.last() - xnodes.last())
                / b.length_x,
                derivative=True,
            )

        @self.Expression(tset, ixnodes, izfaces, comps)
        def dcdz(b, t, ix, iz, i):
            return _interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=znodes,
                faces=zfaces,
                phi_func=lambda izf: b.conc[t, ix, izf, i] / b.length_z[None],
                phi_bound_0=0,  # solid wall no flux
                phi_bound_1=0,  # solid wall no flux
                derivative=True,
            )

        @self.Expression(tset, ixfaces, iznodes)
        def dTdx(b, t, ix, iz):
            return _interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=xnodes,
                faces=xfaces,
                phi_func=lambda ixf: b.Dtemp[t, ixf, iz] / b.length_x,
                phi_bound_0=(b.Dtemp[t, ixnodes.first(), iz] - b.Dtemp_x0[t, iz])
                / (xnodes.first() - xfaces.first())
                / b.length_x,
                phi_bound_1=(b.Dtemp_x1[t, iz] - b.Dtemp[t, ixnodes.last(), iz])
                / (xfaces.last() - xnodes.last())
                / b.length_x,
                derivative=True,
            )

        @self.Expression(tset, ixnodes, izfaces)
        def dTdz(b, t, ix, iz):
            return _interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=znodes,
                faces=zfaces,
                phi_func=lambda izf: b.temperature[t, ix, izf] / b.length_z[None],
                phi_bound_0=0,
                phi_bound_1=0,
                derivative=True,
            )

        @self.Expression(tset, ixfaces, iznodes, comps)
        def diff_eff_coeff_xfaces(b, t, ix, iz, i):
            return _interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=xnodes,
                faces=xfaces,
                phi_func=lambda ixf: b.diff_eff_coeff[t, ixf, iz, i],
                # TODO we can probably use conc_x0 and conc_x1 now
                phi_bound_0=b.diff_eff_coeff[
                    t, ixnodes.first(), iz, i
                ],  # use node value
                phi_bound_1=b.diff_eff_coeff[
                    t, ixnodes.last(), iz, i
                ],  # use node value
                derivative=False,
            )

        @self.Expression(tset, ixnodes, izfaces, comps)
        def diff_eff_coeff_zfaces(b, t, ix, iz, i):
            return _interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=znodes,
                faces=zfaces,
                phi_func=lambda izf: b.diff_eff_coeff[t, ix, izf, i],
                phi_bound_0=0,  # solid wall no flux
                phi_bound_1=0,  # solid wall no flux
                derivative=False,
            )

        @self.Expression(tset, ixfaces, iznodes)
        def temperature_xfaces(b, t, ix, iz):
            return _interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=xnodes,
                faces=xfaces,
                phi_func=lambda ixf: b.temperature[t, ixf, iz],
                phi_bound_0=b.temperature_x0[t, iz],
                phi_bound_1=b.temperature_x1[t, iz],
                derivative=False,
            )

        @self.Expression(tset, ixnodes, izfaces)
        def temperature_zfaces(b, t, ix, iz):
            return _interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=znodes,
                faces=zfaces,
                phi_func=lambda izf: b.temperature[t, ix, izf],
                phi_bound_0=b.temperature[t, ix, iznodes.first()],
                phi_bound_1=b.temperature[t, ix, iznodes.last()],
                derivative=False,
            )

        @self.Expression(tset, ixfaces, iznodes, comps)
        def xflux(b, t, ix, iz, i):
            return -b.dcdx[t, ix, iz, i] * b.diff_eff_coeff_xfaces[t, ix, iz, i]

        @self.Constraint(tset, iznodes, comps)
        def xflux_x0_eqn(b, t, iz, i):
            return b.xflux[t, ixfaces.first(), iz, i] == b.xflux_x0[t, iz, i]

        @self.Constraint(tset, iznodes, comps)
        def xflux_x1_eqn(b, t, iz, i):
            return b.xflux[t, ixfaces.last(), iz, i] == b.xflux_x1[t, iz, i]

        @self.Expression(tset, ixnodes, izfaces, comps)
        def zflux(b, t, ix, iz, i):
            return -b.dcdz[t, ix, iz, i] * b.diff_eff_coeff_zfaces[t, ix, iz, i]

        @self.Expression(tset, ixfaces, iznodes)
        def qxflux(b, t, ix, iz):
            return -(1 - b.porosity) * b.solid_thermal_conductivity * b.dTdx[t, ix, iz]

        @self.Expression(tset, ixnodes, izfaces)
        def qzflux(b, t, ix, iz):
            return -(1 - b.porosity) * b.solid_thermal_conductivity * b.dTdz[t, ix, iz]

        @self.Constraint(tset, iznodes)
        def qflux_x0_eqn(b, t, iz):
            return b.qflux_x0[t, iz] == b.qxflux[t, ixfaces.first(), iz]

        @self.Constraint(tset, iznodes)
        def qflux_x1_eqn(b, t, iz):
            return b.qflux_x1[t, iz] == b.qxflux[t, ixfaces.last(), iz]

        @self.Expression(tset, ixnodes, iznodes)
        def resistivity(b, t, ix, iz):
            return (
                pyo.units.ohm
                * pyo.units.m
                * pyo.exp(
                    b.resistivity_log_preexponential_factor
                    + b.resistivity_thermal_exponent_dividend / b.temperature[t, ix, iz]
                )
            )

        @self.Expression(tset, ixnodes, iznodes)
        def resistance(b, t, ix, iz):
            return (
                b.resistivity[t, ix, iz]
                * b.length_x
                * b.dx[ix]
                / b.xface_area[iz]
                / (1 - b.porosity)
            )

        @self.Expression(tset, iznodes)
        def current(b, t, iz):
            return b.current_density[t, iz] * b.xface_area[iz]

        @self.Expression(tset, ixnodes, iznodes)
        def voltage_drop(b, t, ix, iz):
            return b.current[t, iz] * b.resistance[t, ix, iz]

        @self.Expression(tset, iznodes)
        def resistance_total(b, t, iz):
            return sum(b.resistance[t, ix, iz] for ix in ixnodes)

        @self.Expression(tset, iznodes)
        def voltage_drop_total(b, t, iz):
            return sum(b.voltage_drop[t, ix, iz] for ix in ixnodes)

        @self.Constraint(tset, ixnodes, iznodes, comps)
        def material_balance_eqn(b, t, ix, iz, i):
            return b.node_volume[ix, iz] * b.dcdt[t, ix, iz, i] == b.xface_area[iz] * (
                b.xflux[t, ix, iz, i] - b.xflux[t, ix + 1, iz, i]
            ) + b.zface_area[ix] * (b.zflux[t, ix, iz, i] - b.zflux[t, ix, iz + 1, i])

        if dynamic:
            self.material_balance_eqn[tset.first(), :, :, :].deactivate()

        @self.Expression(tset, ixnodes, iznodes)
        def joule_heating(b, t, ix, iz):
            return b.current[t, iz] * b.voltage_drop[t, ix, iz]

        @self.Constraint(tset, ixnodes, iznodes)
        def energy_balance_solid_eqn(b, t, ix, iz):
            return (
                b.node_volume[ix, iz]
                * (
                    b.porosity * b.dcedt[t, ix, iz]
                    + (1 - b.porosity) * b.dcedt_solid[t, ix, iz]
                )
                == b.xface_area[iz] * (b.qxflux[t, ix, iz] - b.qxflux[t, ix + 1, iz])
                + b.zface_area[ix] * (b.qzflux[t, ix, iz] - b.qzflux[t, ix, iz + 1])
                + b.joule_heating[t, ix, iz]
                # For mass flux heat transfer include exchange with channel
                # probably make little differece, but want to ensure the energy
                # balance closes
                + b.xface_area[iz]
                * sum(
                    b.xflux[t, ix, iz, i]
                    * comp_enthalpy_expr(b.temperature_xfaces[t, ix, iz], i)
                    for i in comps
                )
                - b.xface_area[iz]
                * sum(
                    b.xflux[t, ix + 1, iz, i]
                    * comp_enthalpy_expr(b.temperature_xfaces[t, ix + 1, iz], i)
                    for i in comps
                )
                + b.zface_area[ix]
                * sum(
                    b.zflux[t, ix, iz, i]
                    * comp_enthalpy_expr(b.temperature_zfaces[t, ix, iz], i)
                    for i in comps
                )
                - b.zface_area[ix]
                * sum(
                    b.zflux[t, ix, iz + 1, i]
                    * comp_enthalpy_expr(b.temperature_zfaces[t, ix, iz + 1], i)
                    for i in comps
                )
            )

        if dynamic:
            self.energy_balance_solid_eqn[tset.first(), :, :].deactivate()

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        temperature_guess=None,
        pressure_guess=None,
        mole_frac_guess=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        comps = self.config.comp_list
        # Use conc_x0 and the temperature guess to start filling in initial
        # guess values.
        for t in self.flowsheet().time:
            for iz in self.iznodes:
                for i in comps:
                    _set_if_unfixed(
                        self.Dconc[t, self.ixnodes.first(), iz, i],
                        self.Dconc_x0[t, iz, i],
                    )
                for ix in self.ixnodes:
                    if temperature_guess is not None:
                        _set_if_unfixed(self.temperature_z[t, iz], temperature_guess)
                        _set_if_unfixed(self.Dtemp[t, ix, iz], 0)
                        _set_if_unfixed(self.Dtemp_x0[t, iz], 0)
                        _set_if_unfixed(self.Dtemp_x1[t, iz], 0)
                    if pressure_guess is not None:
                        _set_if_unfixed(self.pressure[t, ix, iz], pressure_guess)
                    for i in comps:
                        _set_if_unfixed(self.Dconc[t, ix, iz, i], 0)
                    mol_dens = pyo.value(sum(self.conc[t, ix, iz, i] for i in comps))
                    _set_if_unfixed(
                        self.pressure[t, ix, iz],
                        _constR * self.temperature[t, ix, iz] * mol_dens,
                    )
                    for i in comps:
                        _set_if_unfixed(
                            self.mole_frac_comp[t, ix, iz, i],
                            self.conc[t, ix, iz, i] / mol_dens,
                        )
                    _set_if_unfixed(
                        self.enth_mol[t, ix, iz],
                        sum(
                            comp_enthalpy_expr(self.temperature[t, ix, iz], i)
                            * self.mole_frac_comp[t, ix, iz, i]
                            for i in comps
                        ),
                    )
                    if self.config.has_holdup:
                        _set_if_unfixed(
                            self.int_energy_mol[t, ix, iz],
                            sum(
                                comp_int_energy_expr(self.temperature[t, ix, iz], i)
                                * self.mole_frac_comp[t, ix, iz, i]
                                for i in comps
                            ),
                        )
                        _set_if_unfixed(
                            self.int_energy_density[t, ix, iz],
                            self.int_energy_mol[t, ix, iz]
                            / self.volume_molar[t, ix, iz],
                        )
                        # _set_if_unfixed(
                        #    self.int_energy_density_solid[t, ix, iz],
                        #    self.solid_heat_capacity * self.solid_density * (self.temperature[t, ix, iz] - 1000 * pyo.units.K)
                        # )
                for i in comps:
                    _set_if_unfixed(
                        self.Dconc_x1[t, iz, i],
                        self.Dconc[t, self.ixnodes.last(), iz, i],
                    )

        slvr = get_solver(solver, optarg)
        _init_solve_block(self, slvr, solve_log)

    def calculate_scaling_factors(self):
        pass

    def model_check(self, steady_state=True):
        comp_set = set(self.comps)
        elements_present = set()

        for element in element_list:
            include_element = False
            for species in species_list:
                # Floating point equality take warning!
                if species in comp_set and element_dict[element][species] != 0:
                    include_element = True
            if include_element:
                elements_present.add(element)

        if not steady_state:
            # Mass and energy conservation equations steady state only at present
            return
        for t in self.flowsheet().config.time:
            for element in element_list:
                if element not in elements_present:
                    continue
                sum_in = 0
                sum_out = 0
                for iz in self.iznodes:
                    sum_in += sum(
                        element_dict[element][j]
                        * self.xflux_x0[t, iz, j]
                        * self.xface_area[iz]
                        for j in self.comps
                    )
                    sum_out += sum(
                        element_dict[element][j]
                        * self.xflux_x1[t, iz, j]
                        * self.xface_area[iz]
                        for j in self.comps
                    )
                normal = max(
                    pyo.value(sum_in), pyo.value(sum_out), 1e-8
                )  # FIXME justify this number
                fraction_change = pyo.value((sum_out - sum_in) / normal)
                if abs(fraction_change) > 3e-3:
                    raise RuntimeError(
                        f"{element} is not being conserved in {self.name}; "
                        f"fractional change {fraction_change}"
                    )
            enth_in = 0
            enth_out = 0

            for iz in self.iznodes:
                enth_in += self.xface_area[iz] * (
                    self.qflux_x0[t, iz]
                    + sum(
                        comp_enthalpy_expr(self.temperature_x0[t, iz], j)
                        * self.xflux_x0[t, iz, j]
                        for j in self.comps
                    )
                )
                enth_out += self.xface_area[iz] * (
                    self.qflux_x1[t, iz]
                    + sum(
                        comp_enthalpy_expr(self.temperature_x1[t, iz], j)
                        * self.xflux_x1[t, iz, j]
                        for j in self.comps
                    )
                )
            total_joule_heating = sum(
                sum(self.joule_heating[t, ix, iz] for ix in self.ixnodes)
                for iz in self.iznodes
            )

            normal = max(
                pyo.value(abs(enth_in)),
                pyo.value(abs(enth_out)),
                pyo.value(abs(total_joule_heating)),
                1e-4,
            )  # FIXME justify this number
            fraction_change = pyo.value(
                (enth_out - enth_in - total_joule_heating) / normal
            )
            if abs(fraction_change) > 3e-3:
                raise RuntimeError(
                    f"Energy is not being conserved in {self.name}; "
                    f"fractional change {fraction_change}"
                )

    def recursive_scaling(self):
        gsf = iscale.get_scaling_factor
        ssf = _set_scaling_factor_if_none
        sgsf = _set_and_get_scaling_factor
        sdf = _set_default_factor
        cst = lambda c, s: iscale.constraint_scaling_transform(c, s, overwrite=False)
        sR = 1e-1  # Scaling factor for R
        sD = 1e5  # Heuristic scaling factor for diffusion coefficient
        sy_def = 10  # Mole frac comp scaling
        sh = 1e-2  # Heat xfer coeff
        sH = 1e-4  # Enthalpy/int energy
        sk = 10  # Fudge factor to scale Dtemp
        sLx = sgsf(self.length_x, len(self.ixnodes) / self.length_x.value)
        sLy = 1 / self.length_y[None].value
        sLz = len(self.iznodes) / self.length_z[None].value

        for t in self.flowsheet().time:
            for iz in self.iznodes:
                if not self.temperature_z[t, iz].is_reference():
                    sT = sgsf(self.temperature_z[t, iz], 1e-2)

                if self.qflux_x0[t, iz].is_reference():
                    sq0 = gsf(self.qflux_x0[t, iz].referent, default=1e-2)
                else:
                    sq0 = sgsf(self.qflux_x0[t, iz], 1e-2)
                cst(self.qflux_x0_eqn[t, iz], sq0)
                if not self.Dtemp_x0.is_reference():
                    ssf(self.Dtemp_x0, sq0 * sLx / sk)

                if self.qflux_x1[t, iz].is_reference():
                    sq1 = gsf(self.qflux_x1[t, iz].referent, default=1e-2)
                else:
                    sq1 = sgsf(self.qflux_x1[t, iz], 1e-2)
                cst(self.qflux_x1_eqn[t, iz], sq1)
                if not self.Dtemp_x1.is_reference():
                    ssf(self.Dtemp_x1, sq1 * sLx / sk)

                sqx = min(sq0, sq1)
                sqz = 10 * sqx  # Heuristic

                sxflux = {}
                for j in self.comps:
                    # ssf(self.conc_ref[t,iz,j],sy_def*1e-4/(sR*sT))

                    if self.xflux_x0[t, iz, j].is_reference():
                        sxflux0 = gsf(self.xflux_x0[t, iz, j].referent, default=1e-2)
                    else:
                        sxflux0 = sgsf(self.xflux_x0[t, iz, j], 1e-2)
                    cst(self.xflux_x0_eqn[t, iz, j], sxflux0)
                    if not self.Dconc_x0[t, iz, j].is_reference():
                        ssf(self.Dconc_x0[t, iz, j], sxflux0 * sLx / sD)

                    if self.xflux_x1[t, iz, j].is_reference():
                        sxflux1 = gsf(self.xflux_x1[t, iz, j].referent, default=1e-2)
                    else:
                        sxflux1 = sgsf(self.xflux_x1[t, iz, j], 1e-2)
                    cst(self.xflux_x1_eqn[t, iz, j], sxflux1)

                    if not self.Dconc_x1[t, iz, j].is_reference():
                        ssf(self.Dconc_x1[t, iz, j], sxflux1 * sLx / sD)

                    sxflux[j] = min(sxflux0, sxflux1)

                for ix in self.ixnodes:
                    sP = sgsf(self.pressure[t, ix, iz], 1e-4)
                    sV = sR * sT / sP

                    sDT = sgsf(self.Dtemp[t, ix, iz], sqx * sLx / sk)
                    sH = sgsf(self.enth_mol[t, ix, iz], sH)
                    cst(self.enth_mol_eqn[t, ix, iz], sH)

                    if self.config.has_holdup:
                        sU = sgsf(self.int_energy_mol[t, ix, iz], sH)
                        cst(self.int_energy_mol_eqn[t, ix, iz], sU)

                        s_rho_U = sgsf(self.int_energy_density[t, ix, iz], sU / sV)
                        cst(self.int_energy_density_eqn[t, ix, iz], s_rho_U)

                        s_rho_U_solid = sgsf(
                            self.int_energy_density_solid[t, ix, iz],
                            1
                            / (
                                self.solid_heat_capacity.value
                                * self.solid_density.value
                                * sDT
                            ),
                        )
                        cst(self.int_energy_density_solid_eqn[t, ix, iz], s_rho_U_solid)

                    cst(self.mole_frac_eqn[t, ix, iz], 1)
                    cst(self.energy_balance_solid_eqn[t, ix, iz], sqx * sLy * sLz)

                    for j in self.comps:
                        sy = sgsf(self.mole_frac_comp[t, ix, iz, j], sy_def)
                        cst(self.conc_eqn[t, ix, iz, j], sy * sP)
                        ssf(self.Dconc[t, ix, iz, j], sxflux[j] * sLx / sD)

                        cst(
                            self.material_balance_eqn[t, ix, iz, j],
                            sxflux[j] * sLy * sLz,
                        )


@declare_process_block_class("SocTriplePhaseBoundary")
class SocTriplePhaseBoundaryData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
            **default** = useDefault.
            **Valid values:** {
            **useDefault** - get flag from parent (default = False),
            **True** - set as a dynamic model,
            **False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(domain=In([False]), default=False),
    )
    CONFIG.declare(
        "comp_list",
        ConfigValue(default=["H2", "H2O"], description="List of components"),
    )
    CONFIG.declare(
        "tpb_stoich_dict",
        ConfigValue(
            description="Stochiometry coefficients for component reactions on "
            "the triple phase boundary.",
        ),
    )
    CONFIG.declare(
        "conc_ref",
        ConfigValue(
            default=None,
            description="Variable for the component concentration in bulk channel ",
        ),
    )
    CONFIG.declare(
        "below_electrolyte",
        ConfigValue(
            domain=In([True, False]),
            description="Variable for the component concentration in bulk channel ",
        ),
    )

    _submodel_boilerplate_config(CONFIG)
    _thermal_boundary_conditions_config(CONFIG, thin=True)
    _material_boundary_conditions_config(CONFIG, thin=True)

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        # z coordinates for nodes and faces
        zfaces = self.config.cv_zfaces
        self.znodes = pyo.Set(
            initialize=[
                (zfaces[i] + zfaces[i + 1]) / 2.0 for i in range(len(zfaces) - 1)
            ]
        )
        comps = self.comps = pyo.Set(initialize=self.config.comp_list)
        self.tpb_stoich = copy.copy(self.config.tpb_stoich_dict)
        # TODO maybe let user specify inert species directly? Floating point
        # equalities make me nervous---Doug
        self.inert_comps = {j for j, coeff in self.tpb_stoich.items() if coeff == 0}
        self.reacting_comps = {
            j for j, coeff in self.tpb_stoich.items() if j not in self.inert_comps
        }
        self.reacting_gases = {j for j in comps if j not in self.inert_comps}

        iznodes = self.iznodes = pyo.Set(initialize=range(1, len(self.znodes) + 1))

        _submodel_boilerplate_create_if_none(self)
        _create_thermal_boundary_conditions_if_none(self, thin=True)
        _create_material_boundary_conditions_if_none(self, thin=True)

        _create_if_none(
            self,
            "conc_ref",
            idx_set=(tset, iznodes),
            units=pyo.units.mol / pyo.units.m**3,
        )

        self.mole_frac_comp = pyo.Var(
            tset,
            iznodes,
            comps,
            initialize=1 / len(comps),
            units=pyo.units.dimensionless,
            bounds=(0, 1),
        )

        self.log_mole_frac_comp = pyo.Var(
            tset,
            iznodes,
            comps,
            initialize=-1,
            units=pyo.units.dimensionless,
            bounds=(None, 0),
        )

        self.activation_potential = pyo.Var(
            tset,
            iznodes,
            initialize=1,
            units=pyo.units.V,
        )

        self.activation_potential_alpha1 = pyo.Var(
            initialize=0.5,
            units=pyo.units.dimensionless,
        )

        self.activation_potential_alpha2 = pyo.Var(
            initialize=0.5,
            units=pyo.units.dimensionless,
        )

        self.exchange_current_exponent_comp = pyo.Var(
            self.reacting_gases,
            initialize=1,
            units=pyo.units.dimensionless,
            bounds=(0, None),
        )

        self.exchange_current_log_preexponential_factor = pyo.Var(
            initialize=1, units=(pyo.units.amp / pyo.units.m**2), bounds=(0, None)
        )

        self.exchange_current_activation_energy = pyo.Var(
            initialize=0, units=pyo.units.J / pyo.units.mol, bounds=(0, None)
        )

        @self.Expression(tset, iznodes, comps)
        def conc(b, t, iz, j):
            return b.conc_ref[t, iz, j] + b.Dconc[t, iz, j]

        @self.Expression(tset, iznodes)
        def pressure(b, t, iz):
            return sum(b.conc[t, iz, i] for i in comps) * _constR * b.temperature[t, iz]

        # mole_frac_comp must be a variable because we want IPOPT to enforce
        # a lower bound of 0 in order to avoid AMPL errors, etc.
        @self.Constraint(tset, iznodes, comps)
        def mole_frac_comp_eqn(b, t, iz, j):
            return b.mole_frac_comp[t, iz, j] == b.conc[t, iz, j] / sum(
                b.conc[t, iz, i] for i in comps
            )

        @self.Constraint(tset, iznodes, comps)
        def log_mole_frac_comp_eqn(b, t, iz, j):
            return b.mole_frac_comp[t, iz, j] == pyo.exp(b.log_mole_frac_comp[t, iz, j])

        @self.Expression(tset, iznodes)
        def ds_rxn(b, t, iz):
            T = b.temperature[t, iz]
            P = b.pressure[t, iz]
            P_ref = 1e5 * pyo.units.Pa
            log_y_j = b.log_mole_frac_comp
            nu_j = b.tpb_stoich
            # Any j not in comps is assumed to not be vapor phase
            pressure_exponent = sum(nu_j[j] for j in b.reacting_gases)
            if abs(pressure_exponent) < 1e-6:
                out_expr = 0
            else:
                out_expr = -_constR * pressure_exponent * pyo.log(P / P_ref)
            return out_expr + (
                sum(nu_j[j] * comp_entropy_expr(T, j) for j in b.reacting_comps)
                - _constR
                * sum(
                    nu_j[j] * log_y_j[t, iz, j]
                    for j in b.comps
                    if j not in b.inert_comps
                )
            )

        @self.Expression(tset, iznodes)
        def dh_rxn(b, t, iz):
            T = b.temperature[t, iz]
            nu_j = b.tpb_stoich
            return sum(nu_j[j] * comp_enthalpy_expr(T, j) for j in b.reacting_comps)

        @self.Expression(tset, iznodes)
        def dg_rxn(b, t, iz):
            return b.dh_rxn[t, iz] - b.temperature[t, iz] * b.ds_rxn[t, iz]

        @self.Expression(tset, iznodes)
        def nernst_potential(b, t, iz):
            return -b.dg_rxn[t, iz] / _constF

        @self.Expression(tset, iznodes)
        def log_exchange_current_density(b, t, iz):
            T = b.temperature[t, iz]
            log_k = b.exchange_current_log_preexponential_factor[None]
            expo = b.exchange_current_exponent_comp
            E_A = b.exchange_current_activation_energy[None]
            out = log_k - E_A / (_constR * T)
            for j in b.reacting_gases:
                out += expo[j] * b.log_mole_frac_comp[t, iz, j]
            return out

        # Butler Volmer equation
        @self.Constraint(tset, iznodes)
        def activation_potential_eqn(b, t, iz):
            i = b.current_density[t, iz]
            log_i0 = b.log_exchange_current_density[t, iz]
            eta = b.activation_potential[t, iz]
            T = b.temperature[t, iz]
            alpha_1 = b.activation_potential_alpha1[None]
            alpha_2 = b.activation_potential_alpha2[None]
            exp_expr = _constF * eta / (_constR * T)
            return i == pyo.exp(log_i0 + alpha_1 * exp_expr) - pyo.exp(
                log_i0 - alpha_2 * exp_expr
            )

        @self.Expression(tset, iznodes)
        def reaction_rate_per_unit_area(b, t, iz):
            # Assuming there are no current leaks, the reaction rate can be
            # calculated directly from the current density
            return b.current_density[t, iz] / _constF

        # Put this expression in to prepare for a contact resistance term
        @self.Expression(tset, iznodes)
        def voltage_drop_total(b, t, iz):
            return b.activation_potential[t, iz]

        @self.Constraint(tset, iznodes)
        def qflux_eqn(b, t, iz):
            return (
                b.qflux_x1[t, iz]
                == b.qflux_x0[t, iz]
                + b.current_density[t, iz]
                * b.voltage_drop_total[t, iz]  # Resistive heating
                - b.reaction_rate_per_unit_area[t, iz]  # Reversible heat of reaction
                * b.temperature[t, iz]
                * b.ds_rxn[t, iz]
            )

        @self.Constraint(tset, iznodes, comps)
        def xflux_eqn(b, t, iz, j):
            if b.config.below_electrolyte:
                return (
                    b.xflux[t, iz, j]
                    == -b.reaction_rate_per_unit_area[t, iz] * b.tpb_stoich[j]
                )
            else:
                return (
                    b.xflux[t, iz, j]
                    == b.reaction_rate_per_unit_area[t, iz] * b.tpb_stoich[j]
                )

    def initialize_build(
        self, outlvl=idaeslog.NOTSET, solver=None, optarg=None, fix_x0=False
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        self.Dtemp.fix()
        self.conc_ref.fix()
        self.Dconc.fix()
        if fix_x0:
            self.qflux_x0.fix()
        else:
            self.qflux_x1.fix()

        for t in self.flowsheet().time:
            for iz in self.iznodes:
                denom = pyo.value(sum(self.conc[t, iz, j] for j in self.comps))
                for j in self.comps:
                    self.mole_frac_comp[t, iz, j].value = pyo.value(
                        self.conc[t, iz, j] / denom
                    )
                    self.log_mole_frac_comp[t, iz, j].value = pyo.value(
                        pyo.log(self.mole_frac_comp[t, iz, j])
                    )

        slvr = get_solver(solver, optarg)
        _init_solve_block(self, slvr, solve_log)

        self.Dtemp.unfix()
        self.conc_ref.unfix()
        self.Dconc.unfix()
        if fix_x0:
            self.qflux_x0.unfix()
        else:
            self.qflux_x1.unfix()

    def calculate_scaling_factors(self):
        pass

    def recursive_scaling(self):
        gsf = iscale.get_scaling_factor
        ssf = _set_scaling_factor_if_none
        sgsf = _set_and_get_scaling_factor
        sdf = _set_default_factor
        cst = lambda c, s: iscale.constraint_scaling_transform(c, s, overwrite=False)
        sR = 1e-1  # Scaling factor for R
        sy_def = 10  # Mole frac comp scaling
        # sLy = sgsf(self.length_y,1/self.length_y[None].value)
        # sLz = sgsf(self.length_z,len(self.iznodes)/self.length_z[None].value)
        sLy = 1 / self.length_y[None].value
        sLz = len(self.iznodes) / self.length_z[None].value

        for t in self.flowsheet().time:
            for iz in self.iznodes:
                ssf(self.activation_potential[t, iz], 10)
                if self.current_density[t, iz].is_reference():
                    si = gsf(self.current_density[t, iz].referent, default=1e-2)
                else:
                    si = gsf(self.current_density[t, iz], default=1e-2, warning=True)
                # TODO come back when I come up with a good formulation
                cst(self.activation_potential_eqn[t, iz], si)
                if self.qflux_x0[t, iz].is_reference():
                    gsf(self.qflux_x0[t, iz].referent, default=1e-2)
                else:
                    sqx0 = sgsf(self.qflux_x0[t, iz], 1e-2)
                if self.qflux_x1[t, iz].is_reference():
                    sqx1 = gsf(self.qflux_x1[t, iz].referent, 1e-2)
                else:
                    sqx1 = sgsf(self.qflux_x1[t, iz], 1e-2)
                sqx = min(sqx0, sqx1)
                cst(self.qflux_eqn[t, iz], sqx)
                for j in self.comps:
                    # TODO Come back with good formulation for trace components
                    # and scale DConc and Cref
                    sy = sgsf(self.mole_frac_comp[t, iz, j], sy_def)
                    ssf(self.log_mole_frac_comp[t, iz, j], 1)
                    cst(self.mole_frac_comp_eqn[t, iz, j], sy)
                    cst(self.log_mole_frac_comp_eqn[t, iz, j], sy)
                    sxflux = sgsf(self.xflux[t, iz, j], 1e-2)
                    cst(self.xflux_eqn[t, iz, j], sxflux)


@declare_process_block_class("SocContactResistor")
class SocContactResistorData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(domain=In([False]), default=False),
    )

    _submodel_boilerplate_config(CONFIG)
    _thermal_boundary_conditions_config(CONFIG, thin=True)

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        # z coordinates for nodes and faces
        zfaces = self.config.cv_zfaces
        self.znodes = pyo.Set(
            initialize=[
                (zfaces[i] + zfaces[i + 1]) / 2.0 for i in range(len(zfaces) - 1)
            ]
        )
        # This sets provide an integer index for nodes and faces
        iznodes = self.iznodes = pyo.Set(initialize=range(1, len(self.znodes) + 1))

        _submodel_boilerplate_create_if_none(self)
        _create_thermal_boundary_conditions_if_none(self, thin=True)

        # Preexponential factor needs to be given in units of  omh*m**2
        self.log_preexponential_factor = pyo.Var(
            initialize=-50, units=pyo.units.dimensionless
        )
        # units=pyo.units.ohm*pyo.units.m**2,)
        self.thermal_exponent_dividend = pyo.Var(initialize=0, units=pyo.units.K)
        self.contact_fraction = pyo.Var(
            initialize=1, units=pyo.units.dimensionless, bounds=(0, 1)
        )

        @self.Expression(tset, iznodes)
        def contact_resistance(b, t, iz):
            return (
                1
                * pyo.units.ohm
                * pyo.units.m**2
                / b.contact_fraction
                * pyo.exp(
                    b.log_preexponential_factor
                    + b.thermal_exponent_dividend / (b.temperature[t, iz])
                )
            )

        @self.Expression(tset, iznodes)
        def voltage_drop(b, t, iz):
            return b.contact_resistance[t, iz] * b.current_density[t, iz]

        @self.Expression(tset, iznodes)
        def joule_heating_flux(b, t, iz):
            return b.voltage_drop[t, iz] * b.current_density[t, iz]

        @self.Constraint(tset, iznodes)
        def qflux_eqn(b, t, iz):
            return b.qflux_x1[t, iz] == b.qflux_x0[t, iz] + b.joule_heating_flux[t, iz]

    def initialize_build(
        self, outlvl=idaeslog.NOTSET, solver=None, optarg=None, fix_qflux_x0=True
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        self.Dtemp.fix()
        if fix_qflux_x0:
            self.qflux_x0.fix()
        else:
            self.qflux_x1.fix()

        slvr = get_solver(solver, optarg)
        _init_solve_block(self, slvr, solve_log)

        self.Dtemp.unfix()
        if fix_qflux_x0:
            self.qflux_x0.unfix()
        else:
            self.qflux_x1.unfix()

    def calculate_scaling_factors(self):
        pass

    def recursive_scaling(self):
        gsf = iscale.get_scaling_factor
        for i, c in self.qflux_eqn.items():
            if self.qflux_x0[i].is_reference():
                sq0 = gsf(self.qflux_x0[i].referent)
            else:
                sq0 = gsf(self.qflux_x0[i], warning=True)
            if self.qflux_x1[i].is_reference():
                sq1 = gsf(self.qflux_x1[i].referent)
            else:
                sq1 = gsf(self.qflux_x1[i], warning=True)
            sq = min(sq0, sq1)
            iscale.constraint_scaling_transform(c, sq, overwrite=False)


@declare_process_block_class("SocConductiveSlab")
class SocConductiveSlabData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
                **default** = useDefault.
                **Valid values:** {
                **useDefault** - get flag from parent (default = False),
                **True** - set as a dynamic model,
                **False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(domain=In([useDefault, True, False]), default=useDefault),
    )
    CONFIG.declare(
        "cv_xfaces",
        ConfigValue(
            description="CV x-boundary set, should start with 0 and end with 1."
        ),
    )

    _submodel_boilerplate_config(CONFIG)
    _thermal_boundary_conditions_config(CONFIG, thin=False)

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        # z coordinates for nodes and faces
        self.zfaces = pyo.Set(initialize=self.config.cv_zfaces)
        self.znodes = pyo.Set(
            initialize=[
                (self.zfaces.at(i) + self.zfaces.at(i + 1)) / 2.0
                for i in range(1, len(self.zfaces))
            ]
        )
        self.xfaces = pyo.Set(initialize=self.config.cv_xfaces)
        self.xnodes = pyo.Set(
            initialize=[
                (self.xfaces.at(i) + self.xfaces.at(i + 1)) / 2.0
                for i in range(1, len(self.xfaces))
            ]
        )
        # This sets provide an integer index for nodes and faces
        self.izfaces = pyo.Set(initialize=range(1, len(self.zfaces) + 1))
        self.iznodes = pyo.Set(initialize=range(1, len(self.znodes) + 1))
        self.ixfaces = pyo.Set(initialize=range(1, len(self.xfaces) + 1))
        self.ixnodes = pyo.Set(initialize=range(1, len(self.xnodes) + 1))

        # Space saving aliases
        izfaces = self.izfaces
        iznodes = self.iznodes
        ixfaces = self.ixfaces
        ixnodes = self.ixnodes
        zfaces = self.zfaces
        znodes = self.znodes
        xfaces = self.xfaces
        xnodes = self.xnodes

        # Channel thickness AKA length in the x direction is specific to the
        # channel so local variable here is the only option
        self.length_x = pyo.Var(
            doc="Thickness of the electrode (x-direction)",
            units=pyo.units.m,
        )

        _submodel_boilerplate_create_if_none(self)
        _create_thermal_boundary_conditions_if_none(self, thin=False)

        self.Dtemp = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            initialize=0,
            doc="Temperature at node centers",
            units=pyo.units.K,
        )
        if self.config.has_holdup:
            self.int_energy_density_solid = pyo.Var(
                tset,
                ixnodes,
                iznodes,
                doc="Internal energy density of solid electrode",
                units=pyo.units.J / pyo.units.m**3,
            )
        self.resistivity_log_preexponential_factor = pyo.Var(
            doc="Logarithm of resistivity preexponential factor " "in units of ohm*m",
            units=pyo.units.dimensionless,
        )
        self.resistivity_thermal_exponent_dividend = pyo.Var(
            doc="Parameter divided by temperature in exponential", units=pyo.units.K
        )

        # Parameters
        self.heat_capacity = pyo.Var()
        self.density = pyo.Var()
        self.thermal_conductivity = pyo.Var()

        @self.Expression(tset, ixnodes, iznodes)
        def temperature(b, t, ix, iz):
            if b.config.include_temperature_x_thermo:
                return b.temperature_z[t, iz] + b.Dtemp[t, ix, iz]
            else:
                return b.temperature_z[t, iz]

        if self.config.has_holdup:

            @self.Constraint(tset, ixnodes, iznodes)
            def int_energy_density_solid_eqn(b, t, ix, iz):
                return b.int_energy_density_solid[
                    t, ix, iz
                ] == b.heat_capacity * b.density * (
                    b.temperature[t, ix, iz] - 1000 * pyo.units.K
                )

        if dynamic:
            self.dcedt_solid = DerivativeVar(
                self.int_energy_density_solid,
                wrt=tset,
                initialize=0,
                doc="Internal energy density time derivative",
            )
        else:
            self.dcedt_solid = pyo.Param(
                tset,
                ixnodes,
                iznodes,
                initialize=0,
                units=pyo.units.W / pyo.units.m**3,
            )

        @self.Expression(iznodes)
        def dz(b, iz):
            return b.zfaces.at(iz + 1) - b.zfaces.at(iz)

        @self.Expression(ixnodes)
        def dx(b, ix):
            return b.zfaces.at(ix + 1) - b.zfaces.at(ix)

        @self.Expression(ixnodes, iznodes)
        def node_volume(b, ix, iz):
            return (
                b.length_x[None]
                * b.length_y[None]
                * b.length_z[None]
                * b.dz[iz]
                * b.dx[ix]
            )

        @self.Expression(ixnodes)
        def zface_area(b, ix):
            return b.length_y[None] * b.length_x[None] * b.dx[ix]

        @self.Expression(iznodes)
        def xface_area(b, iz):
            return b.length_y[None] * b.length_z[None] * b.dz[iz]

        @self.Expression(tset, iznodes)
        def current(b, t, iz):
            return b.current_density[t, iz] * b.xface_area[iz]

        @self.Expression(tset, ixfaces, iznodes)
        def dTdx(b, t, ix, iz):
            return _interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=xnodes,
                faces=xfaces,
                phi_func=lambda ixf: b.Dtemp[t, ixf, iz] / b.length_x,
                phi_bound_0=(b.Dtemp_x0[t, iz] - b.Dtemp[t, ixnodes.first(), iz])
                / (xfaces.first() - xnodes.first())
                / b.length_x,
                phi_bound_1=(b.Dtemp[t, ixnodes.last(), iz] - b.Dtemp_x1[t, iz])
                / (xnodes.last() - xfaces.last())
                / b.length_x,
                derivative=True,
            )

        @self.Expression(tset, ixnodes, izfaces)
        def dTdz(b, t, ix, iz):
            return _interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=znodes,
                faces=zfaces,
                phi_func=lambda izf: b.temperature[t, ix, izf] / b.length_z[None],
                phi_bound_0=0,
                phi_bound_1=0,
                derivative=True,
            )

        @self.Expression(tset, ixfaces, iznodes)
        def qxflux(b, t, ix, iz):
            return -b.thermal_conductivity * b.dTdx[t, ix, iz]

        @self.Constraint(tset, iznodes)
        def qflux_x0_eqn(b, t, iz):
            return b.qflux_x0[t, iz] == b.qxflux[t, ixfaces.first(), iz]

        @self.Constraint(tset, iznodes)
        def qflux_x1_eqn(b, t, iz):
            return b.qflux_x1[t, iz] == b.qxflux[t, ixfaces.last(), iz]

        @self.Expression(tset, ixnodes, izfaces)
        def qzflux(b, t, ix, iz):
            return -b.thermal_conductivity * b.dTdz[t, ix, iz]

        @self.Expression(tset, ixnodes, iznodes)
        def resistivity(b, t, ix, iz):
            return (
                pyo.units.ohm
                * pyo.units.m
                * pyo.exp(
                    b.resistivity_log_preexponential_factor
                    + b.resistivity_thermal_exponent_dividend / b.temperature[t, ix, iz]
                )
            )

        @self.Expression(tset, ixnodes, iznodes)
        def resistance(b, t, ix, iz):
            return b.resistivity[t, ix, iz] * b.length_x * b.dx[ix] / b.xface_area[iz]

        @self.Expression(tset, ixnodes, iznodes)
        def voltage_drop(b, t, ix, iz):
            return b.current[t, iz] * b.resistance[t, ix, iz]

        @self.Expression(tset, iznodes)
        def resistance_total(b, t, iz):
            return sum(b.resistance[t, ix, iz] for ix in ixnodes)

        @self.Expression(tset, iznodes)
        def voltage_drop_total(b, t, iz):
            return sum(b.voltage_drop[t, ix, iz] for ix in ixnodes)

        @self.Expression(tset, ixnodes, iznodes)
        def joule_heating(b, t, ix, iz):
            # current_density is the current density so have to multiply it be Area I**2 = i**2*A**2
            # R = rho * dx / Area / (1-porosity) heating = I**2*R
            return b.current[t, iz] * b.voltage_drop[t, ix, iz]

        @self.Constraint(tset, ixnodes, iznodes)
        def energy_balance_solid_eqn(b, t, ix, iz):
            return (
                b.node_volume[ix, iz] * b.dcedt_solid[t, ix, iz]
                == b.xface_area[iz] * (b.qxflux[t, ix, iz] - b.qxflux[t, ix + 1, iz])
                + b.zface_area[ix] * (b.qzflux[t, ix, iz] - b.qzflux[t, ix, iz + 1])
                + b.joule_heating[t, ix, iz]
            )

        if dynamic:
            self.energy_balance_solid_eqn[tset.first(), :, :].deactivate()

    def calculate_scaling_factors(self):
        pass

    def model_check(self):
        pass

    def recursive_scaling(self):
        gsf = iscale.get_scaling_factor
        ssf = _set_scaling_factor_if_none
        sgsf = _set_and_get_scaling_factor
        sdf = _set_default_factor
        cst = lambda c, s: iscale.constraint_scaling_transform(c, s, overwrite=False)
        sR = 1e-1  # Scaling factor for R
        sD = 1e4  # Heuristic scaling factor for diffusion coefficient
        sy_def = 10  # Mole frac comp scaling
        sh = 1e-2  # Heat xfer coeff
        sH = 1e-4  # Enthalpy/int energy
        sk = 1  # Thermal conductivity is ~1
        sLx = sgsf(self.length_x, len(self.ixnodes) / self.length_x.value)
        # sLy = sgsf(self.length_y,1/self.length_y[None].value)
        # sLz = sgsf(self.length_z,len(self.iznodes)/self.length_z[None].value)
        sLy = 1 / self.length_y[None].value
        sLz = len(self.iznodes) / self.length_z[None].value

        for t in self.flowsheet().time:
            for iz in self.iznodes:
                if not self.temperature_z[t, iz].is_reference():
                    sT = sgsf(self.temperature_z[t, iz], 1e-2)

                if self.qflux_x0[t, iz].is_reference():
                    sq0 = gsf(self.qflux_x0[t, iz].referent, default=1e-2)
                else:
                    sq0 = sgsf(self.qflux_x0[t, iz], 1e-2)
                cst(self.qflux_x0_eqn[t, iz], sq0)
                if not self.Dtemp_x0.is_reference():
                    ssf(self.Dtemp_x0, sq0 * sLx / sk)

                if self.qflux_x1[t, iz].is_reference():
                    sq1 = gsf(self.qflux_x1[t, iz].referent, default=1e-2)
                else:
                    sq1 = sgsf(self.qflux_x1[t, iz], 1e-2)
                cst(self.qflux_x1_eqn[t, iz], sq1)
                if not self.Dtemp_x1.is_reference():
                    ssf(self.Dtemp_x1, sq1 * sLx / sk)

                sqx = min(sq0, sq1)
                sqz = 10 * sqx  # Heuristic

                for ix in self.ixnodes:
                    sDT = sgsf(self.Dtemp[t, ix, iz], sqx * sLx / sk)

                    if self.config.has_holdup:
                        s_rho_U_solid = sgsf(
                            self.int_energy_density_solid[t, ix, iz],
                            1 / (self.heat_capacity.value * self.density.value * sDT),
                        )
                        cst(self.int_energy_density_solid_eqn[t, ix, iz], s_rho_U_solid)

                    cst(self.energy_balance_solid_eqn[t, ix, iz], sqx * sLy * sLz)


@declare_process_block_class("SolidOxideCell")
class SolidOxideCellData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
                **default** = useDefault.
                **Valid values:** {
                **useDefault** - get flag from parent (default = False),
                **True** - set as a dynamic model,
                **False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(domain=In([useDefault, True, False]), default=useDefault),
    )
    CONFIG.declare(
        "cv_zfaces",
        ConfigValue(
            description="CV z-boundary set, should start with 0 and end with 1."
        ),
    )
    CONFIG.declare(
        "cv_xfaces_fuel_electrode",
        ConfigValue(
            description="CV x-boundary set for electrode producing or consuming"
            "fuel, should start with 0 and end with 1."
        ),
    )
    CONFIG.declare(
        "cv_xfaces_oxygen_electrode",
        ConfigValue(
            description="CV x-boundary set for electrode producing or consuming"
            "oxygen, should start with 0 and end with 1."
        ),
    )
    CONFIG.declare(
        "cv_xfaces_electrolyte",
        ConfigValue(
            description="CV x-boundary set for electrolyte, should start with "
            "0 and end with 1."
        ),
    )
    CONFIG.declare(
        "fuel_comps",
        ConfigValue(
            default=["H2", "H2O"], description="List of components in fuel stream"
        ),
    )
    CONFIG.declare(
        "oxygen_comps",
        ConfigValue(
            default=["O2", "H2O"], description="List of components in the oxygen stream"
        ),
    )
    # TODO do we want them to provide an electron term in stoich dict?
    # improve documentation
    CONFIG.declare(
        "fuel_tpb_stoich_dict",
        ConfigValue(
            default={"H2": -0.5, "H2O": 0.5},
            description="Dictionary describing the stoichiometry of a "
            "reaction to produce one electron in the fuel "
            "electrode.",
        ),
    )
    CONFIG.declare(
        "oxygen_tpb_stoich_dict",
        ConfigValue(
            default={"O2": -0.25, "H2O": 0.0},
            description="Dictionary describing the stoichiometry of a "
            "reaction to consume one electron in the oxygen "
            "electrode.",
        ),
    )
    CONFIG.declare(
        "flow_pattern",
        ConfigValue(
            default=HeatExchangerFlowPattern.countercurrent,
            domain=In(HeatExchangerFlowPattern),
            description="Co-current or counter-current flow pattern",
        ),
    )
    CONFIG.declare(
        "flux_through_interconnect",
        ConfigValue(
            default=False,
            description="If True write periodic constraint "
            "to model flux through interconnect.",
        ),
    )
    CONFIG.declare(
        "interpolation_scheme",
        ConfigValue(
            default=CV_Interpolation.UDS,
            description="Method used to interpolate face values",
        ),
    )
    # Setting this to false caused initialization issues, so I'm forcing it to
    # be true until I figure out whether those issues can be fixed ---Doug
    CONFIG.declare(
        "include_temperature_x_thermo",
        ConfigValue(
            domain=In([True]),
            default=True,
            description="Whether to consider temperature variations in "
            "x direction in thermodynamic equations",
        ),
    )

    def build(self):
        super().build()
        has_holdup = self.config.has_holdup
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        t0 = tset.first()

        self.fuel_comps = pyo.Set(initialize=self.config.fuel_comps)
        self.oxygen_comps = pyo.Set(initialize=self.config.oxygen_comps)
        zfaces = self.zfaces = self.config.cv_zfaces

        iznodes = self.iznodes = pyo.Set(initialize=range(1, len(zfaces)))
        self.current_density = pyo.Var(
            tset, iznodes, initialize=0, units=pyo.units.A / pyo.units.m**2
        )
        self.potential = pyo.Var(tset, initialize=1.25, units=pyo.units.V)

        include_temp_x_thermo = self.config.include_temperature_x_thermo

        self.temperature_z = pyo.Var(
            tset,
            iznodes,
            doc="Temperature indexed by z",
            initialize=1000,
            units=pyo.units.K,
            bounds=(500, _Tmax),
        )
        self.length_z = pyo.Var(
            doc="Length of cell in direction parallel " "to channel flow.",
            initialize=0.05,
            units=pyo.units.m,
            bounds=(0, None),
        )
        self.length_y = pyo.Var(
            doc="Length of cell in direction perpendicular "
            "to both channel flow and current flow.",
            initialize=0.05,
            units=pyo.units.m,
            bounds=(0, None),
        )

        if self.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
            opposite_flow = False
        elif self.config.flow_pattern == HeatExchangerFlowPattern.countercurrent:
            opposite_flow = True
        else:
            raise ConfigurationError(
                "{} SolidOxideCell supports only cocurrent and "
                "countercurrent flow patterns, but flow_type configuration"
                " argument was set to {}.".format(self.name, self.config.flow_pattern)
            )

        self.fuel_chan = SocChannel(
            default={
                "has_holdup": has_holdup,
                "dynamic": dynamic,
                "cv_zfaces": zfaces,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "temperature_z": self.temperature_z,
                "interpolation_scheme": CV_Interpolation.UDS,
                "opposite_flow": False,
                "below_electrode": True,
                "comp_list": self.fuel_comps,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.contact_interconnect_fuel_flow_mesh = SocContactResistor(
            default={
                "has_holdup": False,
                "dynamic": False,
                "cv_zfaces": zfaces,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "temperature_z": self.temperature_z,
                "Dtemp": self.fuel_chan.Dtemp_x0,
                "qflux_x1": self.fuel_chan.qflux_x0,
                "current_density": self.current_density,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.contact_flow_mesh_fuel_electrode = SocContactResistor(
            default={
                "has_holdup": False,
                "dynamic": False,
                "cv_zfaces": zfaces,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "temperature_z": self.temperature_z,
                "Dtemp": self.fuel_chan.Dtemp_x1,
                "qflux_x0": self.fuel_chan.qflux_x1,
                "current_density": self.current_density,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.oxygen_chan = SocChannel(
            default={
                "has_holdup": has_holdup,
                "dynamic": dynamic,
                "cv_zfaces": zfaces,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "temperature_z": self.temperature_z,
                "interpolation_scheme": CV_Interpolation.UDS,
                "below_electrode": False,
                "opposite_flow": opposite_flow,
                "comp_list": self.oxygen_comps,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.contact_interconnect_oxygen_flow_mesh = SocContactResistor(
            default={
                "has_holdup": False,
                "dynamic": False,
                "cv_zfaces": zfaces,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "temperature_z": self.temperature_z,
                "Dtemp": self.oxygen_chan.Dtemp_x1,
                "qflux_x0": self.oxygen_chan.qflux_x1,
                "current_density": self.current_density,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.contact_flow_mesh_oxygen_electrode = SocContactResistor(
            default={
                "has_holdup": False,
                "dynamic": False,
                "cv_zfaces": zfaces,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "temperature_z": self.temperature_z,
                "Dtemp": self.oxygen_chan.Dtemp_x0,
                "qflux_x1": self.oxygen_chan.qflux_x0,
                "current_density": self.current_density,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.fuel_electrode = SocElectrode(
            default={
                "has_holdup": has_holdup,
                "dynamic": dynamic,
                "cv_zfaces": zfaces,
                "cv_xfaces": self.config.cv_xfaces_fuel_electrode,
                "comp_list": self.fuel_comps,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "conc_ref": self.fuel_chan.conc,
                "dconc_refdt": self.fuel_chan.dcdt,
                "Dconc_x0": self.fuel_chan.Dconc_x1,
                "xflux_x0": self.fuel_chan.xflux_x1,
                "qflux_x0": self.contact_flow_mesh_fuel_electrode.qflux_x1,
                "temperature_z": self.temperature_z,
                "Dtemp_x0": self.fuel_chan.Dtemp_x1,
                "current_density": self.current_density,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.oxygen_electrode = SocElectrode(
            default={
                "has_holdup": has_holdup,
                "dynamic": dynamic,
                "cv_zfaces": zfaces,
                "cv_xfaces": self.config.cv_xfaces_oxygen_electrode,
                "comp_list": self.oxygen_chan.comps,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "conc_ref": self.oxygen_chan.conc,
                "dconc_refdt": self.oxygen_chan.dcdt,
                "Dconc_x1": self.oxygen_chan.Dconc_x0,
                "xflux_x1": self.oxygen_chan.xflux_x0,
                "qflux_x1": self.contact_flow_mesh_oxygen_electrode.qflux_x0,
                "temperature_z": self.temperature_z,
                "Dtemp_x1": self.oxygen_chan.Dtemp_x0,
                "current_density": self.current_density,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.fuel_tpb = SocTriplePhaseBoundary(
            default={
                "has_holdup": False,
                "dynamic": False,
                "cv_zfaces": zfaces,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "comp_list": self.fuel_comps,
                "tpb_stoich_dict": self.config.fuel_tpb_stoich_dict,
                "current_density": self.current_density,
                "temperature_z": self.temperature_z,
                "Dtemp": self.fuel_electrode.Dtemp_x1,
                "qflux_x0": self.fuel_electrode.qflux_x1,
                "conc_ref": self.fuel_chan.conc,
                "Dconc": self.fuel_electrode.Dconc_x1,
                "xflux": self.fuel_electrode.xflux_x1,
                "include_temperature_x_thermo": include_temp_x_thermo,
                "below_electrolyte": True,
            }
        )
        self.oxygen_tpb = SocTriplePhaseBoundary(
            default={
                "has_holdup": False,
                "dynamic": False,
                "cv_zfaces": zfaces,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "comp_list": self.oxygen_comps,
                "tpb_stoich_dict": self.config.oxygen_tpb_stoich_dict,
                "current_density": self.current_density,
                "temperature_z": self.temperature_z,
                "Dtemp": self.oxygen_electrode.Dtemp_x0,
                "qflux_x1": self.oxygen_electrode.qflux_x0,
                "conc_ref": self.oxygen_chan.conc,
                "Dconc": self.oxygen_electrode.Dconc_x0,
                "xflux": self.oxygen_electrode.xflux_x0,
                "include_temperature_x_thermo": include_temp_x_thermo,
                "below_electrolyte": False,
            }
        )
        self.electrolyte = SocConductiveSlab(
            default={
                "has_holdup": has_holdup,
                "dynamic": dynamic,
                "cv_zfaces": zfaces,
                "cv_xfaces": self.config.cv_xfaces_electrolyte,
                "length_z": self.length_z,
                "length_y": self.length_y,
                "temperature_z": self.temperature_z,
                "current_density": self.current_density,
                "Dtemp_x0": self.fuel_electrode.Dtemp_x1,
                "Dtemp_x1": self.oxygen_electrode.Dtemp_x0,
                "qflux_x0": self.fuel_tpb.qflux_x1,
                "qflux_x1": self.oxygen_tpb.qflux_x0,
                "include_temperature_x_thermo": include_temp_x_thermo,
            }
        )
        self.state_vars = {"flow_mol", "mole_frac_comp", "temperature", "pressure"}
        for chan, alias in zip([self.fuel_chan, self.oxygen_chan], ["fuel", "oxygen"]):
            setattr(
                self,
                alias + "_inlet",
                Port(
                    initialize={
                        var: getattr(chan, var + "_inlet") for var in self.state_vars
                    }
                ),
            )
            setattr(
                self,
                alias + "_outlet",
                Port(
                    initialize={
                        var: getattr(chan, var + "_outlet") for var in self.state_vars
                    }
                ),
            )

        @self.Expression(tset, iznodes)
        def eta_contact(b, t, iz):
            return (
                b.contact_interconnect_fuel_flow_mesh.voltage_drop[t, iz]
                + b.contact_flow_mesh_fuel_electrode.voltage_drop[t, iz]
                + b.contact_interconnect_oxygen_flow_mesh.voltage_drop[t, iz]
                + b.contact_flow_mesh_oxygen_electrode.voltage_drop[t, iz]
            )

        @self.Expression(tset, iznodes)
        def eta_ohm(b, t, iz):
            return (
                b.electrolyte.voltage_drop_total[t, iz]
                + b.fuel_electrode.voltage_drop_total[t, iz]
                + b.oxygen_electrode.voltage_drop_total[t, iz]
                + b.eta_contact[t, iz]
            )

        if self.config.flux_through_interconnect:
            raise NotImplementedError(
                "Flux through interconnect has not been " "implemented yet"
            )
        else:
            self.contact_interconnect_fuel_flow_mesh.qflux_x0.value = 0
            self.contact_interconnect_oxygen_flow_mesh.qflux_x1.value = 0

            @self.Constraint(tset, iznodes)
            def no_qflux_fuel_interconnect_eqn(b, t, iz):
                return 0 == b.contact_interconnect_fuel_flow_mesh.qflux_x0[t, iz]

            @self.Constraint(tset, iznodes)
            def no_qflux_oxygen_interconnect_eqn(b, t, iz):
                return 0 == b.contact_interconnect_oxygen_flow_mesh.qflux_x1[t, iz]

        @self.Constraint(tset, iznodes)
        def mean_temperature_eqn(b, t, iz):
            return 0 == b.fuel_chan.Dtemp[t, iz] + b.oxygen_chan.Dtemp[t, iz]

        if dynamic:
            self.mean_temperature_eqn[tset.first(), :].deactivate()

        @self.Constraint(tset, iznodes)
        def potential_eqn(b, t, iz):
            return b.potential[t] == b.fuel_tpb.nernst_potential[
                t, iz
            ] + b.oxygen_tpb.nernst_potential[t, iz] - (
                b.eta_ohm[t, iz]
                + b.fuel_tpb.voltage_drop_total[t, iz]
                + b.oxygen_tpb.voltage_drop_total[t, iz]
            )

        # This is net flow of power *into* the cell. In fuel cell mode, this will
        # be negative
        @self.Expression(tset)
        def electrical_work(b, t):
            return -sum(
                b.potential[t] * b.electrolyte.current[t, iz]
                for iz in b.electrolyte.iznodes
            )

    def initialize_build(
        self,
        current_density_guess,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        xflux_guess=None,
        qflux_x1_guess=None,
        qflux_x0_guess=None,
        temperature_guess=None,
        velocity_guess=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        tset = self.flowsheet().config.time
        # t0 = tset.first()

        unfix_inlet_var = {}

        # Save fixedness status of inlet variables and fix them for initialization
        for chan, comps in zip(
            ["fuel", "oxygen"], [self.fuel_comps, self.oxygen_comps]
        ):
            pname = chan + "_inlet"
            p = getattr(self, pname)
            unfix_inlet_var[pname] = {}
            for varname in self.state_vars:
                unfix_inlet_var[pname][varname] = {}
                var = getattr(p, varname)
                for t in tset:
                    if varname == "mole_frac_comp":
                        unfix_inlet_var[pname][varname][t] = {}
                        for j in comps:
                            unfix_inlet_var[pname][varname][t][j] = not var[t, j].fixed
                            var[t, j].fix()
                    else:
                        unfix_inlet_var[pname][varname][t] = not var[t].fixed
                        var[t].fix()

        unfix_potential = {}
        slvr = get_solver(solver, optarg)
        for t in tset:
            unfix_potential[t] = not self.potential[t].fixed
        self.potential_eqn.deactivate()
        self.current_density.fix(current_density_guess)

        self.temperature_z.fix(temperature_guess)
        self.mean_temperature_eqn.deactivate()

        self.contact_interconnect_fuel_flow_mesh.initialize_build(
            outlvl=outlvl, solver=solver, optarg=optarg, fix_qflux_x0=True
        )
        self.contact_interconnect_oxygen_flow_mesh.initialize_build(
            outlvl=outlvl, solver=solver, optarg=optarg, fix_qflux_x0=False
        )

        # Reset the fluxes to zero in case there are stale values
        init_log.info_high("Initializing Fuel Channel")
        self.fuel_chan.xflux_x1.fix(0)
        self.fuel_chan.qflux_x0.fix(0)
        self.fuel_chan.qflux_x1.fix(0)
        self.fuel_chan.initialize_build(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )
        self.fuel_chan.xflux_x1.unfix()
        self.fuel_chan.qflux_x0.unfix()
        self.fuel_chan.qflux_x1.unfix()

        init_log.info_high("Initializing Oxygen Channel")
        self.oxygen_chan.xflux_x0.fix(0)
        self.oxygen_chan.qflux_x0.fix(0)
        self.oxygen_chan.qflux_x1.fix(0)
        self.oxygen_chan.initialize_build(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )
        self.oxygen_chan.xflux_x0.unfix()
        self.oxygen_chan.qflux_x0.unfix()
        self.oxygen_chan.qflux_x1.unfix()

        init_log.info_high("Initializing Contact Resistors")

        self.contact_flow_mesh_fuel_electrode.initialize_build(fix_qflux_x0=True)
        self.contact_flow_mesh_oxygen_electrode.initialize_build(fix_qflux_x0=False)

        init_log.info_high("Calculating reaction rate at current guess")
        self.fuel_tpb.Dtemp.fix(0)
        self.fuel_tpb.Dconc.fix(0)
        self.fuel_tpb.conc_ref.fix()
        self.fuel_tpb.qflux_x1.fix(0)
        for t in self.flowsheet().time:
            for iz in self.fuel_tpb.iznodes:
                denom = pyo.value(
                    sum(self.fuel_tpb.conc[t, iz, j] for j in self.fuel_tpb.comps)
                )
                for j in self.fuel_tpb.comps:
                    self.fuel_tpb.mole_frac_comp[t, iz, j].value = pyo.value(
                        self.fuel_tpb.conc[t, iz, j] / denom
                    )
                    self.fuel_tpb.log_mole_frac_comp[t, iz, j].value = pyo.value(
                        pyo.log(self.fuel_tpb.mole_frac_comp[t, iz, j])
                    )

        _init_solve_block(self.fuel_tpb, slvr, solve_log)

        self.fuel_tpb.Dtemp.unfix()
        self.fuel_tpb.Dconc.unfix()
        self.fuel_tpb.conc_ref.unfix()
        self.fuel_tpb.qflux_x1.unfix()

        self.oxygen_tpb.Dtemp.fix(0)
        self.oxygen_tpb.Dconc.fix(0)
        self.oxygen_tpb.conc_ref.fix()
        self.oxygen_tpb.qflux_x0.fix(0)
        for t in self.flowsheet().time:
            for iz in self.oxygen_tpb.iznodes:
                denom = pyo.value(
                    sum(self.oxygen_tpb.conc[t, iz, j] for j in self.oxygen_tpb.comps)
                )
                for j in self.oxygen_tpb.comps:
                    self.oxygen_tpb.mole_frac_comp[t, iz, j].value = pyo.value(
                        self.oxygen_tpb.conc[t, iz, j] / denom
                    )
                    self.oxygen_tpb.log_mole_frac_comp[t, iz, j].value = pyo.value(
                        pyo.log(self.oxygen_tpb.mole_frac_comp[t, iz, j])
                    )

        _init_solve_block(self.oxygen_tpb, slvr, solve_log)

        self.oxygen_tpb.Dtemp.unfix()
        self.oxygen_tpb.Dconc.unfix()
        self.oxygen_tpb.conc_ref.unfix()
        self.oxygen_tpb.qflux_x0.unfix()

        init_log.info_high("Initializing Fuel Electrode")
        self.fuel_electrode.conc_ref.fix()
        self.fuel_electrode.Dconc_x0.fix()
        self.fuel_electrode.Dtemp_x0.fix()
        self.fuel_electrode.qflux_x0.fix()
        # self.fuel_electrode.Dconc_x0.fix()
        self.fuel_electrode.xflux_x1.fix()

        self.fuel_electrode.initialize_build(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            temperature_guess=temperature_guess,
        )

        self.fuel_electrode.conc_ref.unfix()
        self.fuel_electrode.Dconc_x0.unfix()
        self.fuel_electrode.Dtemp_x0.unfix()
        self.fuel_electrode.qflux_x0.unfix()
        # self.fuel_electrode.Dconc_x0.unfix()
        self.fuel_electrode.xflux_x1.unfix()

        init_log.info_high("Initializing Oxygen Electrode")
        self.oxygen_electrode.conc_ref.fix()
        self.oxygen_electrode.Dconc_x1.fix()
        self.oxygen_electrode.Dtemp_x1.fix()
        self.oxygen_electrode.qflux_x1.fix()
        self.oxygen_electrode.xflux_x0.fix()

        self.oxygen_electrode.initialize_build(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
            temperature_guess=temperature_guess,
        )

        self.oxygen_electrode.conc_ref.unfix()
        self.oxygen_electrode.Dconc_x1.unfix()
        self.oxygen_electrode.Dtemp_x1.unfix()
        self.oxygen_electrode.qflux_x1.unfix()
        self.oxygen_electrode.xflux_x0.unfix()

        init_log.info_high("Initializing Triple Phase Boundaries")
        self.fuel_tpb.initialize_build(
            outlvl=outlvl, solver=solver, optarg=optarg, fix_x0=True
        )
        self.oxygen_tpb.initialize_build(
            outlvl=outlvl, solver=solver, optarg=optarg, fix_x0=False
        )

        self.temperature_z.unfix()
        self.mean_temperature_eqn.activate()

        init_log.info_high("Solving cell with fixed current density")
        _init_solve_block(self, slvr, solve_log)

        self.potential_eqn.activate()
        # TODO---does the user ever have a reason to fix current density
        # besides initialization and initial conditions?
        self.current_density.unfix()
        self.potential.fix()

        init_log.info_high("Solving cell with potential equations active")
        _init_solve_block(self, slvr, solve_log)

        # Unfix any inlet variables that were fixed by initialization
        # To be honest, using a state block would probably have been less work
        for chan, comps in zip(
            ["fuel", "oxygen"], [self.fuel_comps, self.oxygen_comps]
        ):
            pname = chan + "_inlet"
            p = getattr(self, pname)
            for varname in self.state_vars:
                var = getattr(p, varname)
                for t in tset:
                    if varname == "mole_frac_comp":
                        for j in comps:
                            if unfix_inlet_var[pname][varname][t][j]:
                                var[t, j].unfix()
                    else:
                        if unfix_inlet_var[pname][varname][t]:
                            var[t].unfix()

        for t in tset:
            if unfix_potential[t]:
                self.potential[t].unfix()
        init_log.info(f"{self.name} initialization completed successfully.")

    def calculate_scaling_factors(self):
        self.recursive_scaling()

    def model_check(self, steady_state=True):
        self.fuel_chan.model_check()
        self.oxygen_chan.model_check()
        self.fuel_electrode.model_check()
        self.oxygen_electrode.model_check()
        self.electrolyte.model_check()

        # Make sure arguments to safe_log and safe_sqrt
        # are sufficiently large at solution
        for expr in [
            self.temperature_z,
            self.fuel_electrode.temperature_x1,
            self.oxygen_electrode.temperature_x1,
        ]:
            for T in expr.values():
                assert pyo.value(T) / 1000 > _safe_log_eps * 100
        print("No problems with safe_math functions in electrochemistry.")

        comp_set = set(self.fuel_comps)
        comp_set = comp_set.union(self.oxygen_comps)
        elements_present = set()

        for element in element_list:
            include_element = False
            for species in species_list:
                # Floating point equality take warning!
                if species in comp_set and element_dict[element][species] != 0:
                    include_element = True
            if include_element:
                elements_present.add(element)

        if not steady_state:
            # Mass and energy conservation equations steady state only at present
            return
        for t in self.flowsheet().config.time:
            for element in element_list:
                if element not in elements_present:
                    continue
                sum_in = 0
                sum_out = 0
                for chan, comps in zip(
                    [self.fuel_chan, self.oxygen_chan],
                    [self.fuel_comps, self.oxygen_comps],
                ):
                    sum_in += sum(
                        element_dict[element][j] * chan.flow_mol_comp_inlet[t, j]
                        for j in comps
                    )
                    sum_out += sum(
                        element_dict[element][j] * chan.flow_mol_comp_inlet[t, j]
                        for j in comps
                    )
                fraction_change = pyo.value((sum_out - sum_in) / sum_in)
                if abs(fraction_change) > 1e-5:
                    raise RuntimeError(
                        f"{element} is not being conserved {self.name}; "
                        f"fractional change {fraction_change}"
                    )

            enth_in = (
                self.fuel_chan.enth_mol_inlet[t] * self.fuel_chan.flow_mol_inlet[t]
                + self.oxygen_chan.enth_mol_inlet[t]
                * self.oxygen_chan.flow_mol_inlet[t]
            )
            enth_out = (
                self.fuel_chan.enth_mol_outlet[t] * self.fuel_chan.flow_mol_outlet[t]
                + self.oxygen_chan.enth_mol_outlet[t]
                * self.oxygen_chan.flow_mol_outlet[t]
            )
            normal = max(
                pyo.value(abs(enth_in)),
                pyo.value(abs(enth_out)),
                pyo.value(abs(self.electrical_work[t])),
                1e-3,
            )  # FIXME justify this value
            fraction_change = pyo.value(
                (enth_out - enth_in - self.electrical_work[t]) / normal
            )
            if abs(fraction_change) > 3e-3:
                raise RuntimeError(
                    f"Energy is not being conserved in {self.name}; "
                    f"fractional change {fraction_change}"
                )

    def recursive_scaling(self):
        gsf = iscale.get_scaling_factor
        ssf = iscale.set_scaling_factor
        cst = iscale.constraint_scaling_transform
        sdf = _set_default_factor

        sy_def = 10
        s_inert_flux = 1e4

        sdf(self.current_density, 1e-2)
        sdf(self.temperature_z, 1e-2)
        sdf(self.potential, 1)
        sdf(self.fuel_inlet.temperature, 1e-2)
        sdf(self.fuel_inlet.pressure, 1e-4)
        sdf(self.fuel_inlet.flow_mol, 1e5)
        sdf(self.fuel_inlet.mole_frac_comp, sy_def)
        sdf(self.oxygen_inlet.temperature, 1e-2)
        sdf(self.oxygen_inlet.pressure, 1e-4)
        sdf(self.oxygen_inlet.flow_mol, 1e5)
        sdf(self.oxygen_inlet.mole_frac_comp, sy_def)
        sdf(self.length_z, 1 / self.length_z.value)
        sdf(self.length_y, 1 / self.length_y.value)

        iscale.propagate_indexed_component_scaling_factors(self)

        # Need to scale xfluxes by component because inerts have much smaller
        # fluxes than actively reacting species.
        # TODO Revisit when reforming equations are added

        for t in self.flowsheet().time:
            sy_in_fuel = {}
            for j in self.fuel_comps:
                sy_in_fuel[j] = gsf(self.fuel_inlet.mole_frac_comp[t, j], sy_def)
            sy_in_oxygen = {}
            for j in self.oxygen_comps:
                sy_in_oxygen[j] = gsf(self.oxygen_inlet.mole_frac_comp[t, j], sy_def)
            for iz in self.iznodes:
                s_react_flux = gsf(self.current_density[t, iz]) * pyo.value(_constF)
                for j in self.fuel_comps:
                    if abs(self.fuel_tpb.tpb_stoich[j]) > 1e-3:
                        s_flux_j = sy_in_fuel[j] * s_react_flux
                    else:
                        s_flux_j = sy_in_fuel[j] * s_inert_flux
                    for var in [self.fuel_chan.xflux_x1, self.fuel_electrode.xflux_x1]:
                        if gsf(var[t, iz, j]) is None:
                            ssf(var[t, iz, j], s_flux_j)
                    ssf(self.fuel_tpb.mole_frac_comp[t, iz, j], sy_in_fuel[j])
                for j in self.oxygen_comps:
                    if abs(self.oxygen_tpb.tpb_stoich[j]) > 1e-3:
                        s_flux_j = sy_in_oxygen[j] * s_react_flux
                    else:
                        s_flux_j = sy_in_oxygen[j] * s_inert_flux
                    # s_flux_j = sy_in_oxygen[j]*s_mat_flux / max(abs(self.oxygen_tpb.tpb_stoich[j]),0.25)
                    for var in [
                        self.oxygen_chan.xflux_x0,
                        self.oxygen_electrode.xflux_x0,
                    ]:
                        if gsf(var[t, iz, j]) is None:
                            ssf(var[t, iz, j], s_flux_j)
                    ssf(self.oxygen_tpb.mole_frac_comp[t, iz, j], sy_in_oxygen[j])

                s_q_flux = s_react_flux * 1e-4  # Chosen heuristically based on TdS_rxn
                # s_q_flux = s_react_flux*1e-6
                for submodel in [
                    self.fuel_chan,
                    self.contact_interconnect_fuel_flow_mesh,
                    self.contact_flow_mesh_fuel_electrode,
                    self.fuel_tpb,
                    self.oxygen_chan,
                    self.contact_interconnect_oxygen_flow_mesh,
                    self.contact_flow_mesh_oxygen_electrode,
                    self.oxygen_tpb,
                    self.electrolyte,
                ]:
                    if gsf(submodel.qflux_x0[t, iz]) is None:
                        ssf(submodel.qflux_x0[t, iz], s_q_flux)
                    if gsf(submodel.qflux_x1[t, iz]) is None:
                        ssf(submodel.qflux_x1[t, iz], s_q_flux)
                if not self.config.flux_through_interconnect:
                    sq = gsf(
                        self.contact_interconnect_fuel_flow_mesh.qflux_x0[t, iz],
                        default=s_q_flux,
                    )
                    cst(self.no_qflux_fuel_interconnect_eqn[t, iz], sq, overwrite=False)
                    sq = gsf(
                        self.contact_interconnect_oxygen_flow_mesh.qflux_x1[t, iz],
                        default=s_q_flux,
                    )
                    cst(
                        self.no_qflux_oxygen_interconnect_eqn[t, iz],
                        sq,
                        overwrite=False,
                    )
        for idx, con in self.mean_temperature_eqn.items():
            cst(con, 1, overwrite=False)

        for submodel in [
            self.fuel_chan,
            self.contact_interconnect_fuel_flow_mesh,
            self.contact_flow_mesh_fuel_electrode,
            self.fuel_electrode,
            self.fuel_tpb,
            self.oxygen_chan,
            self.contact_interconnect_oxygen_flow_mesh,
            self.contact_flow_mesh_oxygen_electrode,
            self.oxygen_electrode,
            self.oxygen_tpb,
            self.electrolyte,
        ]:
            submodel.recursive_scaling()
