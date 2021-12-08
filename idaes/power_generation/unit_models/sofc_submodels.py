##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

import copy
import numpy as np
from scipy.interpolate import griddata
from itertools import combinations
import enum

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.dae import ContinuousSet, DerivativeVar
import pyomo.environ as pyo

from idaes.core.util.config import is_physical_parameter_block
from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.solvers import use_idaes_solver_configuration_defaults
import idaes
from idaes.core.util.math import safe_log
from idaes.core.util import get_solver
from pyomo.network import Port


_constR = 8.3145 * pyo.units.J / pyo.units.mol / pyo.units.K  # or Pa*m3/K/mol
_constF = 96485 * pyo.units.coulomb / pyo.units.mol


class CV_Direction(enum.Enum):
    Z = 1
    X = 2
    Y = 3


class CV_Bound(enum.Enum):
    EXTRAPOLATE = 1


class CV_Interpolation(enum.Enum):
    UDS = 1  # Upwind difference scheme, exit face same as node center
    CDS = 2  # Linear interpolation from upstream and downstream node centers
    QUICK = 3  # Quadratic upwind interpolation, quadratic from two upwind
    # centers and one downwind


def _set_default_factor(c, s):
    for i in c:
        if iscale.get_scaling_factor(c[i]) is None:
            iscale.set_scaling_factor(c[i], s)


def _set_if_unfixed(v, val):
    if not v.fixed:
        v.set_value(pyo.value(val))


def _interpolate_channel(
    iz, ifaces, nodes, faces, phi_func, phi_inlet, method, oposite_flow
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
        oposite_flow: if True assume velocity is negative

    Returns:
        expression for phi at face[iz]
    """
    # I don't always need these, but it doesn't take long to calculate them
    if not oposite_flow:
        izu = iz - 1  # adjacent node upstream of the face
        izuu = iz - 2  # node upstream adacjent to node upstream adjacent to face
        izd = iz  # downstream node adjacent to face
    else:
        izu = iz  # adjacent node upstream of the face
        izuu = iz + 1  # node upstream adacjent to node upstream adjacent to face
        izd = iz - 1  # downstream node adjacent to face
    if iz == ifaces.first() and not oposite_flow:
        return phi_inlet
    if iz == ifaces.last() and oposite_flow:
        return phi_inlet
    if method == CV_Interpolation.UDS:
        return phi_func(izu)
    if method == CV_Interpolation.CDS:
        zu = nodes.at(izu)
        if (oposite_flow and iz == ifaces.first()) or (
            not oposite_flow and iz == ifaces.last()
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
    direction=CV_Direction.Z,
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


def contour_grid_data(var, time, xnodes, znodes):
    nx = len(xnodes)
    nz = len(znodes)
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
    return zg, xg, hg


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
        a / tr ** b + b / pyo.exp(d * tr) + e / pyo.exp(f * tr) + g / pyo.exp(h * tr)
    )
    cm2_to_m2 = 0.01 * 0.01
    Pa_to_bar = 1e-5
    return (
        0.002666
        * cm2_to_m2
        * temperature ** (3 / 2)
        / p
        / Pa_to_bar
        / mab ** 0.5
        / sab ** 2
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
            + d["B"] * t ** 2 / 2.0
            + d["C"] * t ** 3 / 3.0
            + d["D"] * t ** 4 / 4.0
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
            d["A"] * pyo.log(t)
            + d["B"] * t
            + d["C"] * t ** 2 / 2.0
            + d["D"] * t ** 3 / 3.0
            - d["E"] / 2.0 / t ** 2
            + d["G"]
        )
        * pyo.units.J
        / pyo.units.mol
        / pyo.units.K
    )


@declare_process_block_class("SofcChannel")
class SofcChannelData(UnitModelBlockData):
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
        ConfigValue(description="CV boundary set, should start with 0 and end with 1."),
    )
    CONFIG.declare(
        "comp_list",
        ConfigValue(default=["H2", "H2O"], description="List of components"),
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
        "oposite_flow",
        ConfigValue(default=False, description="If True assume velocity is negative"),
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
        is_dynamic = self.config.dynamic
        time_units = self.flowsheet().time_units
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

        # Since the length and width of the cell are general for all the parts
        # of the cell, provide the option of just referencing cell level
        # variables
        if self.config.length_z is None:
            self.length_z = pyo.Var(
                initialize=0.25,
                doc="Length in the direction of flow (z-direction)",
                units=pyo.units.m,
            )
        else:
            self.length_z = pyo.Reference(self.config.length_z)

        if self.config.length_y is None:
            self.length_y = pyo.Var(
                initialize=0.25, doc="Width of cell (y-direction)", units=pyo.units.m
            )
        else:
            self.length_y = pyo.Reference(self.config.length_y)

        # Channel thickness AKA length in the x direction is specific to the
        # channel so local variable here is the only option
        self.length_x = pyo.Var(
            doc="Thickness from interconnect to electrode (x-direction)",
            units=pyo.units.m,
        )

        # X direction flux variables
        self.xflux = pyo.Var(
            tset,
            iznodes,
            comps,
            doc="Component flux to electrode (positive is out)",
            initialize=0,
            units=pyo.units.mol / pyo.units.m ** 2 / time_units,
        )
        self.qflux_x1 = pyo.Var(
            tset,
            iznodes,
            doc="Heat flux to electrode (positive is out)",
            initialize=0,
            units=pyo.units.J / pyo.units.m ** 2 / time_units,
        )
        self.qflux_x0 = pyo.Var(
            tset,
            iznodes,
            doc="Heat flux to interconnect (positive is out)",
            initialize=0,
            units=pyo.units.J / pyo.units.m ** 2 / time_units,
        )
        self.htc = pyo.Var(
            tset,
            iznodes,
            doc="Local channel heat fransfer coefficient",
            initialize=500,
            units=pyo.units.J / pyo.units.m ** 2 / time_units / pyo.units.K,
        )
        #
        self.flow_mol = pyo.Var(
            tset,
            iznodes,
            doc="Molar flow in the z-direction through faces",
            units=pyo.units.mol / time_units,
        )
        self.conc = pyo.Var(
            tset,
            iznodes,
            comps,
            doc="Component concentration at node centers",
            units=pyo.units.mol / pyo.units.m ** 3,
        )
        self.enth_mol = pyo.Var(
            tset,
            iznodes,
            doc="Molar enthalpy at node centers",
            units=pyo.units.J / pyo.units.mol,
        )
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
            units=pyo.units.J / pyo.units.m ** 3,
        )
        self.velocity = pyo.Var(
            tset,
            iznodes,
            doc="Fluid velocity at node centers",
            units=pyo.units.m / time_units,
        )
        self.pressure = pyo.Var(
            tset, iznodes, doc="Pressure at node centers", units=pyo.units.Pa
        )
        self.temperature = pyo.Var(
            tset, iznodes, doc="Bulk temperature at node centers", units=pyo.units.K
        )
        self.temperature_el = pyo.Var(
            tset, iznodes, doc="Temperature at electrode interface", units=pyo.units.K
        )
        self.mole_frac_comp = pyo.Var(
            tset, iznodes, comps, doc="Component mole fraction at node centers"
        )
        self.flow_mol_inlet = pyo.Var(tset, doc="Inlet face molar flow rate")
        self.pressure_inlet = pyo.Var(tset, doc="Inlet pressure")
        self.temperature_inlet = pyo.Var(tset, doc="Inlet temperature")
        self.temperature_outlet = pyo.Var(tset, doc="Outlet temperature")
        self.mole_frac_comp_inlet = pyo.Var(
            tset, comps, doc="Inlet compoent mole fractions"
        )

        # Add time derivative varaible if steady state use const 0.
        if is_dynamic:
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
                units=pyo.units.mol / pyo.units.m ** 3 / time_units,
            )

        # Add time derivative varaible if steady state use const 0.
        if is_dynamic:
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
                units=pyo.units.J / pyo.units.m ** 3 / time_units,
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
        def electrode_area(b, iz):
            return b.length_z[None] * b.length_y[None] * b.dz[iz]

        @self.Expression(tset, iznodes)
        def volume_molar(b, t, iz):
            return _constR * b.temperature[t, iz] / b.pressure[t, iz]

        @self.Expression(tset)
        def volume_molar_inlet(b, t):
            return _constR * b.temperature_inlet[t] / b.pressure_inlet[t]

        @self.Expression(tset)
        def enth_mol_inlet(b, t):
            return sum(
                comp_enthalpy_expr(b.temperature_inlet[t], i)
                * b.mole_frac_comp_inlet[t, i]
                for i in comps
            )

        @self.Constraint(tset, iznodes)
        def flow_mol_eqn(b, t, iz):
            # either way the flow goes, want the flow rate to be positive, but
            # in the oposite flow cases want flux and velocity to be negative
            if self.config.oposite_flow:
                return (
                    b.flow_mol[t, iz]
                    == -b.flow_area * b.velocity[t, iz] / b.volume_molar[t, iz]
                )
            return (
                b.flow_mol[t, iz]
                == b.flow_area * b.velocity[t, iz] / b.volume_molar[t, iz]
            )

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
        def zflux_inlet(b, t, i):
            # either way the flow goes, want the flow rate to be positive, but
            # in the oposite flow cases want flux and velocity to be negative
            if self.config.oposite_flow:
                return -b.flow_mol_inlet[t] / b.flow_area * b.mole_frac_comp_inlet[t, i]
            return b.flow_mol_inlet[t] / b.flow_area * b.mole_frac_comp_inlet[t, i]

        @self.Expression(tset)
        def zflux_enth_inlet(b, t):
            # either way the flow goes, want the flow rate to be positive, but
            # in the oposite flow cases want flux and velocity to be negative
            if self.config.oposite_flow:
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
                oposite_flow=self.config.oposite_flow,
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
                oposite_flow=self.config.oposite_flow,
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
                oposite_flow=self.config.oposite_flow,
            )

        @self.Constraint(tset, iznodes, comps)
        def mass_balance_eqn(b, t, iz, i):
            if t == tset.first() and is_dynamic:
                return pyo.Constraint.Skip
            else:
                return (
                    b.dcdt[t, iz, i] * b.node_volume[iz]
                    == b.flow_area * (b.zflux[t, iz, i] - b.zflux[t, iz + 1, i])
                    - b.electrode_area[iz] * b.xflux[t, iz, i]
                )

        @self.Constraint(tset, iznodes)
        def energy_balance_eqn(b, t, iz):
            if t == tset.first() and is_dynamic:
                return pyo.Constraint.Skip
            else:
                return (
                    b.dcedt[t, iz] * b.node_volume[iz]
                    == b.flow_area * (b.zflux_enth[t, iz] - b.zflux_enth[t, iz + 1])
                    - b.electrode_area[iz]
                    * sum(
                        b.xflux[t, iz, i] * comp_enthalpy_expr(b.temperature_el[t, iz], i)
                        for i in comps
                    )
                    - b.qflux_x1[t, iz] * b.electrode_area[iz]
                    - b.qflux_x0[t, iz] * b.electrode_area[iz]
                )

        @self.Constraint(tset, iznodes)
        def temperature_el_eqn(b, t, iz):
            return self.qflux_x1[t, iz] == self.htc[t, iz]*(self.temperature[t, iz] - self.temperature_el[t, iz])

        # For convenience define outlet expressions
        if self.config.oposite_flow:
            izfout = izfaces.first()
            iznout = iznodes.first()
        else:
            izfout = izfaces.last()
            iznout = iznodes.last()

        @self.Expression(tset, comps)
        def flow_mol_comp_outlet(b, t, i):
            return b.zflux[t, izfout, i] * b.flow_area

        @self.Expression(tset)
        def flow_mol_outlet(b, t):
            return sum(b.flow_mol_comp_outlet[t, i] for i in comps)

        @self.Expression(tset)
        def pressure_outlet(b, t):
            return b.pressure_face[t, izfout]

        @self.Expression(tset, comps)
        def mole_frac_comp_outlet(b, t, i):
            return b.flow_mol_comp_outlet[t, i] / b.flow_mol_outlet[t]

        @self.Expression(tset)
        def enth_mol_outlet(b, t):
            return b.zflux_enth[t, izfout] * b.flow_area / b.flow_mol_outlet[t]

        # know enthalpy need a constraint to back calculate temperature
        @self.Constraint(tset)
        def temperature_outlet_eqn(b, t):
            return b.enth_mol_outlet[t] == sum(
                comp_enthalpy_expr(b.temperature_outlet[t], i)
                * b.mole_frac_comp_outlet[t, i]
                for i in comps
            )

    def initialize(
        self,
        outlvl=idaeslog.DEBUG,
        solver=None,
        optarg=None,
        xflux_guess=None,
        qflux_x1_guess=None,
        qflux_x0_guess=None,
        velocity_guess=None,
    ):
        tset = self.flowsheet().config.time
        t0 = tset.first()

        for t in tset:
            _set_if_unfixed(self.temperature_outlet[t], self.temperature_inlet[t])
            for iz in self.iznodes:
                _set_if_unfixed(self.temperature[t, iz], self.temperature_inlet[t])
                _set_if_unfixed(self.temperature_el[t, iz], self.temperature_inlet[t])
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

        solver = get_solver(solver, optarg)
        solver.solve(self, tee=True)

    def calculate_scaling_factors(self):
        # Base default scaling on typical conditions and dimensions
        _set_default_factor(self.pressure, 1e-5)
        _set_default_factor(self.temperature, 1e-2)
        _set_default_factor(self.conc, 1e-1)
        _set_default_factor(self.mole_frac_comp, 10)
        _set_default_factor(self.flow_mol, 1e4)
        _set_default_factor(self.velocity, 10)
        _set_default_factor(self.xflux, 100)
        iscale.propagate_indexed_component_scaling_factors(self)
        s_length_y = 1 / pyo.value(self.length_y[None])
        s_length_z = 1 / pyo.value(self.length_z[None])
        s_length_x = 1 / pyo.value(self.length_x[None])
        _set_default_factor(self.length_x, s_length_x)
        _set_default_factor(self.length_y, s_length_y)
        _set_default_factor(self.length_z, s_length_z)
        s_deltaz = len(self.iznodes) / s_length_z

        # estimate conc scaling
        for i in self.conc:
            if iscale.get_scaling_factor(self.conc[i]) is None:
                sp = iscale.get_scaling_factor(self.pressure[i[0], i[1]])
                sx = iscale.get_scaling_factor(self.mole_frac_comp[i])
                st = iscale.get_scaling_factor(self.temperature[i[0], i[1]])
                sr = 1  # scale for gas constant
                iscale.set_scaling_factor(self.conc[i], sp * sx / sr / st)

        for i in self.conc_eqn:
            sp = iscale.get_scaling_factor(self.pressure[i[0], i[1]])
            sx = iscale.get_scaling_factor(self.mole_frac_comp[i])
            iscale.constraint_scaling_transform(self.conc_eqn[i], sp * sx)

        for i, c in self.flow_mol_eqn.items():
            s = iscale.get_scaling_factor(self.flow_mol[i])
            iscale.constraint_scaling_transform(c, s)

        for i in self.mole_frac_eqn:
            iscale.constraint_scaling_transform(self.mole_frac_eqn[i], 10)

        # for i in self.p_eqn:
        #    sP = iscale.get_scaling_factor(self.pressure[i])
        #    iscale.constraint_scaling_transform(self.p_eqn[i], sP)

        for t, iz, i in self.mass_balance_eqn:
            sc = iscale.get_scaling_factor(self.conc[t, iz, i])
            sv = iscale.get_scaling_factor(self.velocity[t, iz])
            iscale.constraint_scaling_transform(
                self.mass_balance_eqn[t, iz, i], sc * sv / s_deltaz
            )


@declare_process_block_class("SofcElectrode")
class SofcElectrodeData(UnitModelBlockData):
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
        "length_z",
        ConfigValue(default=None, description="Length in the z-direction"),
    )
    CONFIG.declare(
        "length_y", ConfigValue(default=None, description="Width of cell (y-direction)")
    )
    CONFIG.declare(
        "current_density",
        ConfigValue(default=None, description="Optional current_density variable"),
    )
    CONFIG.declare(
        "conc_x0",
        ConfigValue(
            default=None, description="Variable for the component concentration at x=0"
        ),
    )
    CONFIG.declare(
        "xflux_x0",
        ConfigValue(default=None, description="Variable for fluid flux at x=0"),
    )
    CONFIG.declare(
        "qflux_x0",
        ConfigValue(default=None, description="Variable for heat flux at x=0"),
    )
    CONFIG.declare(
        "channel_tempeature",
        ConfigValue(default=None, description="Variable for channel temperature"),
    )
    CONFIG.declare(
        "channel_htc",
        ConfigValue(
            default=None, description="Variable for channel heat transfer coefficent"
        ),
    )
    CONFIG.declare(
        "xflux_fluid_enth_x0",
        ConfigValue(
            default=None, description="Variable for fluid enthalpy flux at x=0"
        ),
    )
    CONFIG.declare(
        "tpb_stoich_dict",
        ConfigValue(
            default=None,
            description="Stochiometry coefficients for component reactions on "
            "the tripple phase boundary.",
        ),
    )

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        is_dynamic = self.config.dynamic
        time_units = self.flowsheet().time_units
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

        self.tpb_stoich = copy.copy(self.config.tpb_stoich_dict)

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

        # Since the length and width of the cell are general for all the parts
        # of the cell, provide the option of just referencing cell level
        # variables
        if self.config.length_z is None:
            self.length_z = pyo.Var(
                initialize=0.25,
                doc="Length in the direction z-direction",
                units=pyo.units.m,
            )
        else:
            self.length_z = pyo.Reference(self.config.length_z)

        if self.config.length_y is None:
            self.length_y = pyo.Var(
                initialize=0.25, doc="Width of cell (y-direction)", units=pyo.units.m
            )
        else:
            self.length_y = pyo.Reference(self.config.length_y)

        # Channel thickness AKA length in the x direction is specific to the
        # channel so local variable here is the only option
        self.length_x = pyo.Var(
            doc="Thickness of the electrode (x-direction)",
            units=pyo.units.m,
        )

        # The concentration along the x=0 boundary comes from the channel, so
        # optionally make a refernce to that external variable.  Otherwise
        # just create a new varaible.  Generally the concentration is calculated
        # at the nodes and interpolated at the faces, so this doesn't overlap
        # with the conc variable, since this is at the x = 0 face
        if self.config.conc_x0 is None:
            self.conc_x0 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Concentration of components at the x=0 boundary",
                initialize=0.0,
                units=pyo.units.mol / pyo.units.m ** 3,
            )
        else:
            self.conc_x0 = pyo.Reference(self.config.conc_x0)

        if self.config.current_density is None:
            self.current_density = pyo.Var(
                tset,
                iznodes,
                doc="Concentration of components at the x=0 boundary",
                initialize=0.0,
                units=pyo.units.amps / pyo.units.m ** 2,
            )
        else:
            self.current_density = pyo.Reference(self.config.current_density)

        #
        self.porosity = pyo.Var(
            initialize=0.50, doc="Electrode porosity", units=pyo.units.dimensionless
        )
        self.tortuosity = pyo.Var(
            initialize=2.0, doc="Electrode tortuosity", units=pyo.units.dimensionless
        )
        self.kec = pyo.Var(units=pyo.units.A)
        self.eec = pyo.Var(units=pyo.units.J/pyo.units.mol)
        self.alpha = pyo.Var()

        # flux variables
        self.zflux = pyo.Var(
            tset,
            ixnodes,
            izfaces,
            comps,
            doc="Component mole flux in z-direction",
            initialize=0,
            units=pyo.units.mol / pyo.units.m ** 2 / time_units,
        )
        self.xflux = pyo.Var(
            tset,
            ixfaces,
            iznodes,
            comps,
            doc="Component mole flux in x-direction",
            initialize=0,
            units=pyo.units.mol / pyo.units.m ** 2 / time_units,
        )

        #
        self.conc = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            comps,
            doc="Component concentration at node centers",
            units=pyo.units.mol / pyo.units.m ** 3,
        )
        self.conc_x1 = pyo.Var(
            tset,
            iznodes,
            comps,
            doc="Concentration of components at the x=1 boundary",
            initialize=0.0,
            units=pyo.units.mol / pyo.units.m ** 3,
        )
        self.enth_mol = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            doc="Molar enthalpy at node centers",
            units=pyo.units.J / pyo.units.mol,
        )
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
            units=pyo.units.J / pyo.units.m ** 3,
        )
        self.int_energy_density_solid = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            doc="Internal energy density of solid electrode",
            units=pyo.units.J / pyo.units.m ** 3,
        )
        self.pressure = pyo.Var(
            tset, ixnodes, iznodes, doc="Pressure at node centers", units=pyo.units.Pa
        )
        # Assume the the electrode gas phase and solid are same temp
        self.temperature = pyo.Var(
            tset, ixnodes, iznodes, doc="Temperature at node centers", units=pyo.units.K
        )
        self.temperature_x0 = pyo.Var(
            tset, iznodes, doc="Temperature at x=0 bound", units=pyo.units.K
        )
        self.temperature_x1 = pyo.Var(
            tset, iznodes, doc="Temperature at x=1 bound", units=pyo.units.K
        )
        self.mole_frac_comp = pyo.Var(
            tset, ixnodes, iznodes, comps, doc="Component mole fraction at node centers"
        )
        self.qflux_x1 = pyo.Var(
            tset,
            iznodes,
            doc="Heat flux from electrode to electrolyte",
            units=pyo.units.J / pyo.units.m ** 2 / time_units,
        )
        self.qrate_per_area_gen_tbp = pyo.Var(
            tset,
            iznodes,
            doc="Heat generated at tripple phase bound by other than reaction",
            units=pyo.units.J / pyo.units.m ** 2 / time_units,
        )
        self.a_res = pyo.Var(
            doc="Resistance preexponential parameter", units=pyo.units.ohm * pyo.units.m
        )
        self.b_res = pyo.Var(doc="Resistance parameter", units=pyo.units.K)

        # Parameters
        self.electrode_heat_capacity = pyo.Var()
        self.electrode_density = pyo.Var()
        self.electrode_thermal_conductivity = pyo.Var()

        # Add time derivative varaible if steady state use const 0.
        if is_dynamic:
            self.dcdt = DerivativeVar(
                self.conc,
                wrt=tset,
                initialize=0,
                doc="Component concentration time derivative",
            )
        else:
            self.dcdt = pyo.Param(
                tset,
                ixnodes,
                iznodes,
                comps,
                initialize=0,
                units=pyo.units.mol / pyo.units.m ** 3 / time_units,
            )
        # Add time derivative varaible if steady state use const 0.
        if is_dynamic:
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
                units=pyo.units.J / pyo.units.m ** 3 / time_units,
            )
        # Add time derivative varaible if steady state use const 0.
        if is_dynamic:
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
                units=pyo.units.J / pyo.units.m ** 3 / time_units,
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

        @self.Expression(tset, iznodes)
        def enth_mol_x0(b, t, iz):
            return sum(
                comp_enthalpy_expr(b.temperature_x0[t, iz], i) * b.conc_x0[t, iz, i]
                for i in comps
            ) / sum(b.conc_x0[t, iz, i] for i in comps)

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
            ] == b.electrode_heat_capacity * b.electrode_density * (
                b.temperature[t, ix, iz] - 1000 * pyo.units.K
            )

        @self.Constraint(tset, ixnodes, iznodes)
        def mole_frac_eqn(b, t, ix, iz):
            return 1 == sum(b.mole_frac_comp[t, ix, iz, i] for i in comps)

        @self.Expression(tset, iznodes, comps)
        def mole_frac_comp_x1(b, t, iz, i):
            return b.conc_x1[t, iz, i] / sum(b.conc_x1[t, iz, j] for j in comps)

        @self.Expression(tset, iznodes)
        def pressure_x1(b, t, iz):
            return (
                _constR
                * sum(b.conc_x1[t, iz, j] for j in comps)
                * b.temperature_x1[t, iz]
            )

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
                phi_func=lambda ixf: b.conc[t, ixf, iz, i] / b.length_x,
                phi_bound_0=(b.conc_x0[t, iz, i] - b.conc[t, ixnodes.first(), iz, i])
                / (xfaces.first() - xnodes.first())
                / b.length_x,
                phi_bound_1=(b.conc_x1[t, iz, i] - b.conc[t, ixnodes.last(), iz, i])
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
                phi_func=lambda izf: b.conc[t, ix, izf, i] / b.length_z,
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
                phi_func=lambda ixf: b.temperature[t, ixf, iz] / b.length_x,
                phi_bound_0=(
                    b.temperature_x0[t, iz] - b.temperature[t, ixnodes.first(), iz]
                )
                / (xfaces.first() - xnodes.first())
                / b.length_x,
                phi_bound_1=(
                    b.temperature[t, ixnodes.last(), iz] - b.temperature_x1[t, iz]
                )
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
                phi_func=lambda izf: b.temperature[t, ix, izf] / b.length_z,
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
                phi_bound_1=0,  # reactions at this bound this is not used
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

        @self.Constraint(tset, iznodes, comps)
        def conc_x1_eqn(b, t, iz, i):
            return (
                b.xflux[t, ixfaces.last(), iz, i]
                == -b.dcdx[t, ixfaces.last(), iz, i]
                * b.diff_eff_coeff_xfaces[t, ixfaces.last(), iz, i]
            )

        @self.Constraint(tset, ixfaces, iznodes, comps)
        def xflux_eqn(b, t, ix, iz, i):
            if ix == ixfaces.last():
                if b.tpb_stoich[i] == 0:
                    return b.xflux[t, ix, iz, i] == 0
                else:
                    return (
                        b.xflux[t, ix, iz, i]
                        == -b.tpb_stoich[i] * b.current_density[t, iz] / _constF
                    )
            return (
                b.xflux[t, ix, iz, i]
                == -b.dcdx[t, ix, iz, i] * b.diff_eff_coeff_xfaces[t, ix, iz, i]
            )

        if self.config.xflux_x0 is not None:
            # if this is provided the xflux at the channel bound should be the
            # same as the flux out of the channel
            self.xflux_x0 = pyo.Reference(self.config.xflux_x0)

            @self.Constraint(tset, iznodes, comps)
            def xflux_from_channel_eqn(b, t, iz, i):
                return b.xflux_x0[t, iz, i] == b.xflux[t, 1, iz, i]

        if (
            self.config.channel_tempeature is not None
            and self.config.channel_htc is not None
            and self.config.qflux_x0 is not None
        ):
            # if this is provided the xflux at the channel bound should be the
            # same as the flux out of the channel
            self.channel_htc = pyo.Reference(self.config.channel_htc)
            self.channel_tempeature = pyo.Reference(self.config.channel_tempeature)
            self.qflux_x0 = pyo.Reference(self.config.qflux_x0)

            @self.Constraint(tset, iznodes)
            def qflux_from_channel_htc_eqn(b, t, iz):
                return b.qflux_x0[t, iz] == b.channel_htc[t, iz] * (
                    b.channel_tempeature[t, iz] - b.temperature_x0[t, iz]
                )

        @self.Expression(tset, ixfaces, iznodes)
        def qxflux(b, t, ix, iz):
            return (
                -(1 - b.porosity) * b.electrode_thermal_conductivity * b.dTdx[t, ix, iz]
            )

        @self.Expression(tset, ixnodes, izfaces)
        def qzflux(b, t, ix, iz):
            return (
                -(1 - b.porosity) * b.electrode_thermal_conductivity * b.dTdz[t, ix, iz]
            )

        @self.Constraint(tset, iznodes)
        def qflux_x0_eqn(b, t, iz):
            return self.qflux_x0[t, iz] == b.qxflux[t, ixfaces.first(), iz]

        @self.Expression(tset, iznodes)
        def h2_tpb_delta_s(b, t, iz):
            return (
                comp_entropy_expr(b.temperature_x1[t, iz], "H2O")
                - comp_entropy_expr(b.temperature_x1[t, iz], "H2")
                - 0.5 * comp_entropy_expr(b.temperature_x1[t, iz], "O2")
            )

        @self.Expression(tset, iznodes)
        def h2_tpb_delta_h(b, t, iz):
            return (
                comp_enthalpy_expr(b.temperature_x1[t, iz], "H2O")
                - comp_enthalpy_expr(b.temperature_x1[t, iz], "H2")
                - 0.5 * comp_enthalpy_expr(b.temperature_x1[t, iz], "O2")
            )

        @self.Expression(tset, iznodes)
        def h2_tpb_delta_g(b, t, iz):
            return b.h2_tpb_delta_h[t, iz] - b.h2_tpb_delta_s[t, iz] * b.temperature_x1[t, iz]

        @self.Expression(tset, iznodes)
        def qflux_rxn_tpb(b, t, iz):
            if "H2" in comps:
                return (
                    - b.xflux[t, ixfaces.last(), iz, "H2"]
                    * b.temperature_x1[t, iz]
                    * b.h2_tpb_delta_s[t, iz]
                )
            else:
                return 0 * pyo.units.J / pyo.units.m ** 2 / time_units

        @self.Expression(tset, iznodes)
        def qflux_tpb(b, t, iz):
            return b.qrate_per_area_gen_tbp[t, iz] + b.qflux_rxn_tpb[t, iz]

        @self.Constraint(tset, iznodes)
        def qflux_x1_eqn(b, t, iz):
            return (
                -b.qxflux[t, ixfaces.last(), iz] + b.qflux_x1[t, iz]
                == b.qflux_tpb[t, iz]
            )

        @self.Constraint(tset, ixnodes, izfaces, comps)
        def zflux_eqn(b, t, ix, iz, i):
            return (
                b.zflux[t, ix, iz, i]
                == -b.dcdz[t, ix, iz, i] * b.diff_eff_coeff_zfaces[t, ix, iz, i]
            )

        @self.Expression(tset, ixnodes, iznodes)
        def resistivity(b, t, ix, iz):
            return b.a_res * pyo.exp(b.b_res / b.temperature[t, ix, iz])

        @self.Expression(tset, ixnodes, iznodes)
        def resistance(b, t, ix, iz):
            return b.resistivity[t, ix, iz] * b.length_x * b.dx[ix] / b.xface_area[iz] / (1 - b.porosity)

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
        def mass_balance_eqn(b, t, ix, iz, i):
            if is_dynamic and t == tset.first():
                return pyo.Constraint.Skip
            return b.node_volume[ix, iz] * b.dcdt[t, ix, iz, i] == b.xface_area[iz] * (
                b.xflux[t, ix, iz, i] - b.xflux[t, ix + 1, iz, i]
            ) + b.zface_area[ix] * (b.zflux[t, ix, iz, i] - b.zflux[t, ix, iz + 1, i])

        @self.Expression(tset, ixnodes, iznodes)
        def joule_heating(b, t, ix, iz):
            return b.current[t, iz] * b.voltage_drop[t, ix, iz]

        @self.Constraint(tset, ixnodes, iznodes)
        def energy_balance_solid_eqn(b, t, ix, iz):
            if is_dynamic and t == tset.first():
                return pyo.Constraint.Skip
            return (
                b.node_volume[ix, iz] * (1 - b.porosity) * b.dcedt_solid[t, ix, iz]
                == b.xface_area[iz] * (b.qxflux[t, ix, iz] - b.qxflux[t, ix + 1, iz])
                + b.zface_area[ix] * (b.qzflux[t, ix, iz] - b.qzflux[t, ix, iz + 1])
                + b.joule_heating[t, ix, iz]
                # For mass flux heat transfer include exchange with channel
                # probably make little differece, but want to ensure the energy
                # balance closes
                + b.xface_area[iz]
                    * sum(
                        b.xflux[t, ix, iz, i] * comp_enthalpy_expr(b.temperature_xfaces[t, ix, iz], i)
                        for i in comps
                    )
                - b.xface_area[iz]
                    * sum(
                        b.xflux[t, ix + 1, iz, i] * comp_enthalpy_expr(b.temperature_xfaces[t, ix + 1, iz], i)
                        for i in comps
                    )
                + b.zface_area[ix]
                    * sum(
                        b.zflux[t, ix, iz, i] * comp_enthalpy_expr(b.temperature_zfaces[t, ix, iz], i)
                        for i in comps
                    )
                - b.zface_area[ix]
                    * sum(
                        b.zflux[t, ix, iz+1, i] * comp_enthalpy_expr(b.temperature_zfaces[t, ix, iz+1], i)
                        for i in comps
                    )
            )

    def initialize(
        self,
        outlvl=idaeslog.DEBUG,
        solver=None,
        optarg=None,
        temperature_guess=None,
        pressure_guess=None,
        mole_frac_guess=None,
    ):
        comps = self.config.comp_list
        # Use conc_x0 and the temperature guess to start filling in initial
        # guess values.
        for t in self.flowsheet().time:
            for iz in self.iznodes:
                for ix in self.ixnodes:
                    if temperature_guess is not None:
                        _set_if_unfixed(self.temperature[t, ix, iz], temperature_guess)
                        _set_if_unfixed(self.temperature_x0[t, iz], temperature_guess)
                        _set_if_unfixed(self.temperature_x1[t, iz], temperature_guess)
                    if pressure_guess is not None:
                        _set_if_unfixed(self.pressure[t, ix, iz], pressure_guess)
                    for i in comps:
                        _set_if_unfixed(self.conc[t, ix, iz, i], self.conc_x0[t, iz, i])
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
                        self.int_energy_mol[t, ix, iz] / self.volume_molar[t, ix, iz],
                    )
                    # _set_if_unfixed(
                    #    self.int_energy_density_solid[t, ix, iz],
                    #    self.electrode_heat_capacity * self.electrode_density * (self.temperature[t, ix, iz] - 1000 * pyo.units.K)
                    # )
        # solver = get_solver(solver, optarg)
        # solver.solve(self, tee=True)

    def calculate_scaling_factors(self):
        # Base default scaling on typical conditions and dimensions
        _set_default_factor(self.pressure, 1e-5)
        _set_default_factor(self.int_energy_density_solid, 1e-10)


@declare_process_block_class("SofcElectrolyte")
class SofcElectrolyteData(UnitModelBlockData):
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
        "cv_xfaces",
        ConfigValue(
            description="CV x-boundary set, should start with 0 and end with 1."
        ),
    )
    CONFIG.declare(
        "length_z",
        ConfigValue(default=None, description="Length in the z-direction"),
    )
    CONFIG.declare(
        "length_y", ConfigValue(default=None, description="Width of cell (y-direction)")
    )
    CONFIG.declare(
        "current_density",
        ConfigValue(default=None, description="Optional current_density variable"),
    )
    CONFIG.declare(
        "temperature_x0", ConfigValue(default=None, description="Temperature at x=0")
    )
    CONFIG.declare(
        "temperature_x1", ConfigValue(default=None, description="Temperature at x=1")
    )
    CONFIG.declare(
        "qflux_x0", ConfigValue(default=None, description="Temperature at x=0")
    )
    CONFIG.declare(
        "qflux_x1", ConfigValue(default=None, description="Temperature at x=1")
    )

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        is_dynamic = self.config.dynamic
        time_units = self.flowsheet().time_units
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

        # Since the length and width of the cell are general for all the parts
        # of the cell, provide the option of just referencing cell level
        # variables
        if self.config.length_z is None:
            self.length_z = pyo.Var(
                initialize=0.25,
                doc="Length in the direction z-direction",
                units=pyo.units.m,
            )
        else:
            self.length_z = pyo.Reference(self.config.length_z)

        if self.config.length_y is None:
            self.length_y = pyo.Var(
                initialize=0.25, doc="Width of cell (y-direction)", units=pyo.units.m
            )
        else:
            self.length_y = pyo.Reference(self.config.length_y)

        # Channel thickness AKA length in the x direction is specific to the
        # channel so local variable here is the only option
        self.length_x = pyo.Var(
            doc="Thickness of the electrode (x-direction)",
            units=pyo.units.m,
        )

        if self.config.current_density is None:
            self.current_density = pyo.Var(
                tset,
                iznodes,
                doc="Concentration of components at the x=0 boundary",
                initialize=0.0,
                units=pyo.units.amps / pyo.units.m ** 2,
            )
        else:
            self.current_density = pyo.Reference(self.config.current_density)

        if self.config.temperature_x0 is None:
            self.temperature_x0 = pyo.Var(
                tset,
                iznodes,
                doc="Temperature at the x=0 boundary",
                initialize=1000,
                units=pyo.units.K,
            )
        else:
            self.temperature_x0 = pyo.Reference(self.config.temperature_x0)

        if self.config.temperature_x1 is None:
            self.temperature_x1 = pyo.Var(
                tset,
                iznodes,
                doc="Temperature at the x=1 boundary",
                initialize=1000,
                units=pyo.units.K,
            )
        else:
            self.temperature_x1 = pyo.Reference(self.config.temperature_x1)

        if self.config.qflux_x0 is None:
            self.qflux_x0 = pyo.Var(
                tset,
                iznodes,
                doc="Heat flux at the x=0 boundary",
                initialize=1000,
                units=pyo.units.J / pyo.units.M ** 2 / time_units,
            )
        else:
            self.qflux_x0 = pyo.Reference(self.config.qflux_x0)

        if self.config.qflux_x1 is None:
            self.qflux_x1 = pyo.Var(
                tset,
                iznodes,
                doc="Heat flux at the x=1 boundary",
                initialize=1000,
                units=pyo.units.J / pyo.units.M ** 2 / time_units,
            )
        else:
            self.qflux_x1 = pyo.Reference(self.config.qflux_x1)

        self.temperature = pyo.Var(
            tset, ixnodes, iznodes, doc="Temperature at node centers", units=pyo.units.K
        )
        self.int_energy_density_solid = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            doc="Internal energy density of solid electrode",
            units=pyo.units.J / pyo.units.m ** 3,
        )
        self.a_res = pyo.Var(
            doc="Resistance preexponential parameter", units=pyo.units.ohm * pyo.units.m
        )
        self.b_res = pyo.Var(doc="Resistance parameter", units=pyo.units.K)

        # Parameters
        self.heat_capacity = pyo.Var()
        self.density = pyo.Var()
        self.thermal_conductivity = pyo.Var()

        @self.Constraint(tset, ixnodes, iznodes)
        def int_energy_density_solid_eqn(b, t, ix, iz):
            return b.int_energy_density_solid[
                t, ix, iz
            ] == b.heat_capacity * b.density * (
                b.temperature[t, ix, iz] - 1000 * pyo.units.K
            )

        if is_dynamic:
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
                units=pyo.units.J / pyo.units.m ** 3 / time_units,
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

        @self.Expression(tset, ixfaces, iznodes)
        def dTdx(b, t, ix, iz):
            return _interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=xnodes,
                faces=xfaces,
                phi_func=lambda ixf: b.temperature[t, ixf, iz] / b.length_x,
                phi_bound_0=(
                    b.temperature_x0[t, iz] - b.temperature[t, ixnodes.first(), iz]
                )
                / (xfaces.first() - xnodes.first())
                / b.length_x,
                phi_bound_1=(
                    b.temperature[t, ixnodes.last(), iz] - b.temperature_x1[t, iz]
                )
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
                phi_func=lambda izf: b.temperature[t, ix, izf] / b.length_z,
                phi_bound_0=0,
                phi_bound_1=0,
                derivative=True,
            )

        @self.Expression(tset, ixfaces, iznodes)
        def qxflux(b, t, ix, iz):
            return -b.thermal_conductivity * b.dTdx[t, ix, iz]

        @self.Constraint(tset, iznodes)
        def qflux_x1_eqn(b, t, iz):
            return -b.qflux_x1[t, iz] == b.qxflux[t, ixfaces.last(), iz]

        @self.Constraint(tset, iznodes)
        def qflux_x0_eqn(b, t, iz):
            return b.qflux_x0[t, iz] == b.qxflux[t, ixfaces.first(), iz]

        @self.Expression(tset, ixnodes, izfaces)
        def qzflux(b, t, ix, iz):
            return -b.thermal_conductivity * b.dTdz[t, ix, iz]

        @self.Expression(tset, ixnodes, iznodes)
        def resistivity(b, t, ix, iz):
            return b.a_res * pyo.exp(b.b_res / b.temperature[t, ix, iz])

        @self.Expression(tset, ixnodes, iznodes)
        def resistance(b, t, ix, iz):
            return b.resistivity[t, ix, iz] * b.length_x * b.dx[ix] / b.xface_area[iz]

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

        @self.Expression(tset, ixnodes, iznodes)
        def joule_heating(b, t, ix, iz):
            # current_density is the current density so have to multiply it be Area I**2 = i**2*A**2
            # R = rho * dx / Area / (1-porosity) heating = I**2*R
            return b.current[t, iz] * b.voltage_drop[t, ix, iz]

        @self.Constraint(tset, ixnodes, iznodes)
        def energy_balance_solid_eqn(b, t, ix, iz):
            if is_dynamic and t == tset.first():
                return pyo.Constraint.Skip
            return (
                b.node_volume[ix, iz] * b.dcedt_solid[t, ix, iz]
                == b.xface_area[iz] * (b.qxflux[t, ix, iz] - b.qxflux[t, ix + 1, iz])
                + b.zface_area[ix] * (b.qzflux[t, ix, iz] - b.qzflux[t, ix, iz + 1])
                + b.joule_heating[t, ix, iz]
            )

    def calculate_scaling_factors(self):
        # Base default scaling on typical conditions and dimensions
        _set_default_factor(self.int_energy_density_solid, 1e-10)


def use_channel():
    from idaes.core import FlowsheetBlock

    dynamic = True
    time_nfe = 10
    time_set = [0, 10] if dynamic else [0]

    zfaces = [0.0, 0.1, 0.2, 0.3, 0.7, 0.8, 0.9, 1.0]
    # zfaces = np.linspace(0, 1, 40).tolist()
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        default={
            "dynamic": dynamic,
            "time_set": time_set,
            "time_units": pyo.units.s,
        }
    )
    m.fs.chan = SofcChannel(
        default={
            "cv_zfaces": zfaces,
            "interpolation_scheme": CV_Interpolation.UDS,
            "oposite_flow": True,
        }
    )

    if dynamic:
        pyo.TransformationFactory("dae.finite_difference").apply_to(
            m.fs, nfe=time_nfe, wrt=m.fs.time, scheme="BACKWARD"
        )
        m.fs.chan.temperature[0, :].fix(1023.15)
        m.fs.chan.flow_mol[0, :].fix(1)
        m.fs.chan.mole_frac_comp[0, :, "H2"].fix(0.7)

    m.fs.chan.temperature_inlet.fix(1023.15)

    m.fs.chan.pressure_inlet.fix(20e5)
    m.fs.chan.pressure.fix(20e5)

    m.fs.chan.flow_mol_inlet.fix(10)
    m.fs.chan.mole_frac_comp_inlet[:, "H2"].fix(0.7)
    m.fs.chan.mole_frac_comp_inlet[:, "H2O"].fix(0.3)

    m.fs.chan.xflux[:, :, "H2"].fix(100)
    m.fs.chan.xflux[:, :, "H2O"].fix(-100)
    m.fs.chan.qflux_x1[:, :].fix(0)
    m.fs.chan.qflux_x0[:, :].fix(0)

    m.fs.chan.length_x.fix(0.01)
    m.fs.chan.length_y.fix(0.01)
    m.fs.chan.length_z.fix(2)

    m.fs.chan.initialize()
    solver = pyo.SolverFactory("ipopt")
    solver.solve(m, tee=True, options={"tol": 1e-6})

    m.fs.chan.display()
    m.fs.chan.temperature.display()

    return m


def use_elecrode():
    import matplotlib.pyplot as plt
    from matplotlib import cm as color_map
    from idaes.core import FlowsheetBlock
    import idaes.core.plugins

    dynamic = False
    time_nfe = 15
    time_set = [0, 10] if dynamic else [0]

    zfaces = np.linspace(0, 1, 11).tolist()
    xfaces_electrode = [0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0]
    xfaces_electrolyte = [0.0, 0.25, 0.5, 0.75, 1.0]
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        default={
            "dynamic": dynamic,
            "time_set": time_set,
            "time_units": pyo.units.s,
        }
    )
    m.fs.fuel_chan = SofcChannel(
        default={
            "cv_zfaces": zfaces,
            "interpolation_scheme": CV_Interpolation.UDS,
            "oposite_flow": False,
            "comp_list": ["H2", "H2O"],
        }
    )
    m.fs.oxygen_chan = SofcChannel(
        default={
            "cv_zfaces": zfaces,
            "interpolation_scheme": CV_Interpolation.UDS,
            "oposite_flow": False,
            "comp_list": ["O2", "H2O"],
        }
    )
    m.fs.fuel_electrode = SofcElectrode(
        default={
            "cv_zfaces": zfaces,
            "cv_xfaces": xfaces_electrode,
            "tpb_stoich_dict": {"H2": -0.5, "H2O": 0.5},
            "comp_list": ["H2", "H2O"],
            "conc_x0": m.fs.fuel_chan.conc,
            "xflux_x0": m.fs.fuel_chan.xflux,
            "qflux_x0": m.fs.fuel_chan.qflux_x1,
            "channel_tempeature": m.fs.fuel_chan.temperature,
            "channel_htc": m.fs.fuel_chan.htc,
        }
    )
    m.fs.oxygen_electrode = SofcElectrode(
        default={
            "cv_zfaces": zfaces,
            "cv_xfaces": xfaces_electrode,
            "tpb_stoich_dict": {"O2": -0.25, "H2O": 0.0},
            "comp_list": ["O2", "H2O"],
            "current_density": m.fs.fuel_electrode.current_density,
            "conc_x0": m.fs.oxygen_chan.conc,
            "xflux_x0": m.fs.oxygen_chan.xflux,
            "qflux_x0": m.fs.oxygen_chan.qflux_x1,
            "channel_tempeature": m.fs.oxygen_chan.temperature,
            "channel_htc": m.fs.oxygen_chan.htc,
        }
    )
    m.fs.electrolyte = SofcElectrolyte(
        default={
            "cv_zfaces": zfaces,
            "cv_xfaces": xfaces_electrolyte,
            "current_density": m.fs.fuel_electrode.current_density,
            "temperature_x0": m.fs.fuel_electrode.temperature_x1,
            "qflux_x0": m.fs.fuel_electrode.qflux_x1,
            "temperature_x1": m.fs.oxygen_electrode.temperature_x1,
            "qflux_x1": m.fs.oxygen_electrode.qflux_x1,
        }
    )

    if dynamic:
        pyo.TransformationFactory("dae.finite_difference").apply_to(
            m.fs, nfe=time_nfe, wrt=m.fs.time, scheme="BACKWARD"
        )
        m.fs.fuel_chan.temperature[0, :].fix(1023.15)
        m.fs.fuel_chan.flow_mol[0, :].fix(1e-5)
        m.fs.fuel_chan.mole_frac_comp[0, :, "H2"].fix(0.1)

        m.fs.oxygen_chan.temperature[0, :].fix(1023.15)
        m.fs.oxygen_chan.flow_mol[0, :].fix(1e-5)
        m.fs.oxygen_chan.mole_frac_comp[0, :, "O2"].fix(0.10)

        m.fs.fuel_electrode.conc[0, :, :, :].fix(
            pyo.value(1.02e5 / _constR / 1023.15 / 2)
        )
        m.fs.fuel_electrode.temperature[0, :, :].fix(1023.15)

        m.fs.oxygen_electrode.conc[0, :, :, :].fix(
            pyo.value(1.02e5 / _constR / 1023.15 / 2)
        )
        m.fs.oxygen_electrode.temperature[0, :, :].fix(1023.15)

        m.fs.electrolyte.temperature[0, :, :].fix(1023.15)
        m.fs.fuel_electrode.temperature_x1[0, :].fix(1023.15)
        m.fs.oxygen_electrode.temperature_x1[0, :].fix(1023.15)

    m.fs.fuel_chan.temperature_inlet.fix(1023.15)
    m.fs.fuel_chan.pressure_inlet.fix(1.02e5)
    m.fs.fuel_chan.pressure.fix(1.02e5)
    m.fs.fuel_chan.flow_mol_inlet.fix(1e-4)
    m.fs.fuel_chan.mole_frac_comp_inlet[:, "H2"].fix(0.10)
    m.fs.fuel_chan.mole_frac_comp_inlet[:, "H2O"].fix(0.90)
    m.fs.fuel_chan.xflux[:, :, "H2"].set_value(0.0)
    m.fs.fuel_chan.xflux[:, :, "H2O"].set_value(0.0)
    m.fs.fuel_chan.qflux_x1[:, :].fix(0)
    m.fs.fuel_chan.qflux_x0[:, :].fix(0)
    m.fs.fuel_chan.length_x.fix(0.002)
    m.fs.fuel_chan.length_y.fix(0.05)
    m.fs.fuel_chan.length_z.fix(0.05)
    m.fs.fuel_chan.htc.fix(100)

    m.fs.oxygen_chan.temperature_inlet.fix(1023.15)
    m.fs.oxygen_chan.pressure_inlet.fix(1.02e5)
    m.fs.oxygen_chan.pressure.fix(1.02e5)
    m.fs.oxygen_chan.flow_mol_inlet.fix(1e-4)
    m.fs.oxygen_chan.mole_frac_comp_inlet[:, "O2"].fix(0.1)
    m.fs.oxygen_chan.mole_frac_comp_inlet[:, "H2O"].fix(0.9)
    m.fs.oxygen_chan.xflux[:, :, "O2"].set_value(0.0)
    m.fs.oxygen_chan.xflux[:, :, "H2O"].set_value(0.0)
    m.fs.oxygen_chan.qflux_x1[:, :].fix(0)
    m.fs.oxygen_chan.qflux_x0[:, :].fix(0)
    m.fs.oxygen_chan.length_x.fix(0.002)
    m.fs.oxygen_chan.length_y.fix(0.05)
    m.fs.oxygen_chan.length_z.fix(0.05)
    m.fs.oxygen_chan.htc.fix(100)

    m.fs.fuel_electrode.qrate_per_area_gen_tbp.fix(0)
    m.fs.fuel_electrode.pressure[:, :, :].set_value(1.02e5)
    m.fs.fuel_electrode.current_density.fix(-3000)
    m.fs.fuel_electrode.length_x.fix(750e-6)
    m.fs.fuel_electrode.length_y.fix(0.05)
    m.fs.fuel_electrode.length_z.fix(0.05)
    m.fs.fuel_electrode.porosity.fix(0.30)
    m.fs.fuel_electrode.tortuosity.fix(3.0)
    m.fs.fuel_electrode.electrode_heat_capacity.fix(430)
    m.fs.fuel_electrode.electrode_density.fix(3030)
    m.fs.fuel_electrode.electrode_thermal_conductivity.fix(1.6)
    m.fs.fuel_electrode.a_res.fix(2.98e-5)
    m.fs.fuel_electrode.b_res.fix(-1392.0)

    m.fs.oxygen_electrode.qrate_per_area_gen_tbp.fix(0)
    m.fs.oxygen_electrode.pressure[:, :, :].set_value(1.02e5)
    m.fs.oxygen_electrode.length_x.fix(40e-6)
    m.fs.oxygen_electrode.length_y.fix(0.05)
    m.fs.oxygen_electrode.length_z.fix(0.05)
    m.fs.oxygen_electrode.porosity.fix(0.35)
    m.fs.oxygen_electrode.tortuosity.fix(3.0)
    m.fs.oxygen_electrode.electrode_heat_capacity.fix(430)
    m.fs.oxygen_electrode.electrode_density.fix(3030)
    m.fs.oxygen_electrode.electrode_thermal_conductivity.fix(1.6)
    m.fs.oxygen_electrode.a_res.fix(8.115e-5)
    m.fs.oxygen_electrode.b_res.fix(600.0)

    m.fs.electrolyte.temperature[:, :, :].set_value(1023.15)
    m.fs.electrolyte.length_x.fix(30e-6)
    m.fs.electrolyte.length_y.fix(0.05)
    m.fs.electrolyte.length_z.fix(0.05)
    m.fs.electrolyte.heat_capacity.fix(430)
    m.fs.electrolyte.density.fix(3030)
    m.fs.electrolyte.thermal_conductivity.fix(1.6)
    m.fs.electrolyte.a_res.fix(2.94e-5)
    m.fs.electrolyte.b_res.fix(10350.0)

    iscale.calculate_scaling_factors(m)
    use_idaes_solver_configuration_defaults()

    m.fs.fuel_chan.xflux.fix()
    m.fs.oxygen_chan.xflux.fix()
    m.fs.fuel_chan.initialize()
    m.fs.oxygen_chan.initialize()

    # m.fs.fuel_electrode.qflux_x1.fix(0)
    # m.fs.oxygen_electrode.qflux_x1.fix(0)
    m.fs.fuel_electrode.initialize(temperature_guess=1023.15)
    m.fs.oxygen_electrode.initialize(temperature_guess=1023.15)

    m.fs.fuel_chan.xflux.unfix()
    m.fs.oxygen_chan.xflux.unfix()

    # m.fs.fuel_electrode.qflux_x1.unfix()
    # m.fs.oxygen_electrode.qflux_x1.unfix()
    m.fs.fuel_electrode.channel_tempeature.unfix()
    m.fs.oxygen_electrode.channel_tempeature.unfix()
    m.fs.fuel_electrode.temperature_x0.unfix()
    m.fs.fuel_chan.qflux_x1.unfix()
    m.fs.oxygen_electrode.temperature_x0.unfix()
    m.fs.oxygen_chan.qflux_x1.unfix()

    # see = pyo.TransformationFactory("simple_equality_eliminator")
    # see.apply_to(m)
    solver = pyo.SolverFactory("ipopt")
    # see.revert()
    solver.solve(
        m,
        tee=True,
        symbolic_solver_labels=True,
        options={"tol": 1e-6, "halt_on_ampl_error": "no"},
    )

    m.fs.potential = pyo.Var(m.fs.time, m.fs.fuel_electrode.iznodes, initialize=1.2, units=pyo.units.V)
    m.fs.potential_cell = pyo.Var(m.fs.time, initialize=1.25, units=pyo.units.V)

    m.fs.fuel_electrode.kec.fix(1.35e10)
    m.fs.fuel_electrode.eec.fix(110000)
    m.fs.fuel_electrode.alpha.fix(0.4)
    m.fs.oxygen_electrode.kec.fix(8.7e7*300)
    m.fs.oxygen_electrode.eec.fix(120000)
    m.fs.oxygen_electrode.alpha.fix(0.5)

    @m.fs.Expression(m.fs.time, m.fs.fuel_electrode.iznodes)
    def potential_nernst(b, t, iz):
        T = b.fuel_electrode.temperature_x1[t, iz]
        dg = b.fuel_electrode.h2_tpb_delta_g[t, iz]
        poe = b.oxygen_electrode.pressure_x1[t, iz]
        yh2 = b.fuel_electrode.mole_frac_comp_x1[t, iz, "H2"]
        yh2o = b.fuel_electrode.mole_frac_comp_x1[t, iz, "H2O"]
        yo2 = b.oxygen_electrode.mole_frac_comp_x1[t, iz, "O2"]
        return (
            -dg/2/_constF
            + _constR
            * T
            / 2.0
            / _constF
            * pyo.log(yh2 / yh2o * (poe * yo2 / 1e5) ** 0.5)
        )

    @m.fs.Expression(m.fs.time, m.fs.fuel_electrode.iznodes)
    def eta_ohm(b, t, iz):
        return (
            b.electrolyte.voltage_drop_total[t, iz]
            + b.fuel_electrode.voltage_drop_total[t, iz]
            + b.oxygen_electrode.voltage_drop_total[t, iz]
        )

    @m.fs.Expression(m.fs.time, m.fs.fuel_electrode.iznodes)
    def fe_ecd(b, t, iz):
        k = b.fuel_electrode.kec
        E = b.fuel_electrode.eec
        yh2 = b.fuel_electrode.mole_frac_comp_x1[t, iz, "H2"]
        yh2o = b.fuel_electrode.mole_frac_comp_x1[t, iz, "H2O"]
        T = b.fuel_electrode.temperature_x1[t, iz]
        A = b.fuel_electrode.xface_area[iz] * (1 - m.fs.fuel_electrode.porosity)
        return k * A * yh2 * yh2o * pyo.exp(-E / _constR / T)

    @m.fs.Expression(m.fs.time, m.fs.oxygen_electrode.iznodes)
    def oe_ecd(b, t, iz):
        k = b.oxygen_electrode.kec
        E = b.oxygen_electrode.eec
        yo2 = b.oxygen_electrode.mole_frac_comp_x1[t, iz, "O2"]
        T = b.oxygen_electrode.temperature_x1[t, iz]
        A = b.oxygen_electrode.xface_area[iz] * (1 - m.fs.oxygen_electrode.porosity)
        return k * A * yo2 ** 0.25 * pyo.exp(-E / _constR / T)

    @m.fs.Expression(m.fs.time, m.fs.oxygen_electrode.iznodes)
    def eta_fe(b, t, iz):
        T = b.fuel_electrode.temperature_x1[t, iz]
        ecd = b.fe_ecd[t, iz]
        I = b.fuel_electrode.current[t, iz]
        alpha = b.fuel_electrode.alpha
        # this has the same sign as current, so is negative in soec mode
        return _constR * T / alpha / _constF * pyo.asinh(I / ecd / 2.0)

    @m.fs.Expression(m.fs.time, m.fs.oxygen_electrode.iznodes)
    def eta_oe(b, t, iz):
        T = b.oxygen_electrode.temperature_x1[t, iz]
        ecd = b.oe_ecd[t, iz]
        I = b.oxygen_electrode.current[t, iz]
        alpha = b.oxygen_electrode.alpha
        # this has the same sign as current, so is negative in soec mode
        return _constR * T / alpha / _constF * pyo.asinh(I / ecd / 2.0)


    @m.fs.Constraint(m.fs.time, m.fs.oxygen_electrode.iznodes)
    def potential_eqn(b, t, iz):
        return b.potential_cell[t] == b.potential_nernst[t, iz] - (
            b.eta_ohm[t, iz] + b.eta_fe[t, iz] + b.eta_oe[t, iz]
        )
    m.fs.fuel_electrode.current_density.unfix()
    m.fs.potential_cell.fix(1.29)

    solver.solve(
        m,
        tee=True,
        symbolic_solver_labels=True,
        options={"tol": 1e-6, "halt_on_ampl_error": "no"},
    )

    @m.fs.Constraint(m.fs.time, m.fs.oxygen_electrode.iznodes)
    def activation_heat_oe_eqn(b, t, iz):
        oe = b.oxygen_electrode
        return oe.qrate_per_area_gen_tbp[t, iz] == oe.current[t, iz] * b.eta_oe[t, iz] / oe.xface_area[iz] #/ (1 - oe.porosity)

    @m.fs.Constraint(m.fs.time, m.fs.fuel_electrode.iznodes)
    def activation_heat_fe_eqn(b, t, iz):
        fe = b.fuel_electrode
        return fe.qrate_per_area_gen_tbp[t, iz] == fe.current[t, iz] * b.eta_fe[t, iz] / fe.xface_area[iz] #/ (1 - fe.porosity)

    m.fs.fuel_electrode.qrate_per_area_gen_tbp.unfix()
    m.fs.oxygen_electrode.qrate_per_area_gen_tbp.unfix()

    solver.solve(
        m,
        tee=True,
        symbolic_solver_labels=True,
        options={"tol": 1e-6, "halt_on_ampl_error": "no"},
    )


    z, x, h = contour_grid_data(
        #var=pyo.Reference(m.fs.fuel_electrode.mole_frac_comp[:, :, :, "H2"]),
        # var=m.fs.fuel_electrode.pressure,
        # var=m.fs.fuel_electrode.temperature,
        var=m.fs.electrolyte.temperature,
        time=m.fs.time,
        #xnodes=m.fs.fuel_electrode.xnodes,
        #znodes=m.fs.fuel_electrode.znodes,
        xnodes=m.fs.electrolyte.xnodes,
        znodes=m.fs.electrolyte.znodes,
    )

    fig, ax = plt.subplots()
    # levels = np.linspace(800, 1400, 40)
    levels = 40

    def animate(i):
        ax.clear()
        ax.set_title("Mole Fraction H$_2$")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        img = ax.contourf(z[i], x[i], h[i], levels=levels, cmap="RdYlBu_r")
        if animate.first:
            plt.colorbar(img, ax=ax)
            animate.first = False

    animate.first = True
    import matplotlib.animation as animation

    ani = animation.FuncAnimation(fig, animate, len(m.fs.time), interval=100)
    plt.show()

    ani.save("animation.gif", writer=animation.PillowWriter(fps=2))

    # m.fs.fuel_electrode.xflux_x0.display()
    # m.fs.fuel_electrode.conc_x0.display()
    # m.fs.fuel_electrode.current_density.display()
    # m.fs.fuel_electrode.conc.display()
    # m.fs.fuel_chan.temperature.display()
    # m.fs.oxygen_chan.temperature.display()
    # m.fs.fuel_electrode.mole_frac_comp.display()
    # m.fs.fuel_chan.mole_frac_comp.display()
    # m.fs.fuel_electrode.int_energy_density_solid.display()
    # m.fs.fuel_electrode.temperature.display()
    # m.fs.electrolyte.temperature.display()
    return m


if __name__ == "__main__":
    m = use_elecrode()

    import idaes.core.util.model_statistics as mstat

    # m.fs.chan.temperature.display()
    # m.fs.chan.mole_frac_comp.display()

    m.fs.fuel_chan.velocity.display()
    m.fs.oxygen_chan.velocity.display()
    m.fs.fuel_electrode.qflux_tpb.display()
    m.fs.fuel_electrode.h2_tpb_delta_s.display()
    m.fs.fuel_electrode.conc_x1.display()

    m.fs.potential_nernst.display()
    m.fs.fe_ecd.display()
    m.fs.oe_ecd.display()
    m.fs.potential.display()
    m.fs.fuel_electrode.current_density.display()
    m.fs.electrolyte.resistance_total.display()
    m.fs.fuel_electrode.resistance_total.display()
    m.fs.oxygen_electrode.resistance_total.display()
    m.fs.eta_oe.display()
    m.fs.eta_fe.display()
    m.fs.eta_ohm.display()
    m.fs.fe_ecd.display()
    m.fs.oe_ecd.display()

    m.fs.fuel_electrode.mole_frac_comp_x1.display()
    m.fs.oxygen_electrode.mole_frac_comp_x1.display()
    m.fs.fuel_electrode.temperature_x1.display()
    m.fs.oxygen_electrode.temperature_x1.display()
    m.fs.fuel_electrode.pressure_x1.display()
    m.fs.oxygen_electrode.pressure_x1.display()

    print(mstat.degrees_of_freedom(m))
    print(binary_diffusion_coefficient_expr(590, 1e5, "CO2", "N2"))
    print(pyo.value(comp_enthalpy_expr(1023.15, "H2O")) + 241.8264 * 1000)

    """
    check_scaling=False
    if check_scaling:
        jac, nlp = iscale.get_jacobian(m, scaled=True)
        print("Extreme Jacobian entries:")
        for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, large=100):
            print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]")
        #print("Unscaled constraints:")
        #for c in iscale.unscaled_constraints_generator(m):
        #    print(f"    {c}")
        #print("Scaled constraints by factor:")
        #for c, s in iscale.constraints_with_scale_factor_generator(m):
        #    print(f"    {c}, {s}")
        print("Badly scaled variables:")
        for v, sv in iscale.badly_scaled_var_generator(
            m, large=1e2, small=1e-2, zero=1e-12
        ):
            print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}")
        print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")
    """

    for temperature in [273, 298, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]:
        ds = pyo.value(
            comp_entropy_expr(temperature, "H2O")
            - comp_entropy_expr(temperature, "H2")
            - 0.5 * comp_entropy_expr(temperature, "O2")
        )
        dh = pyo.value(
            comp_enthalpy_expr(temperature, "H2O")
            - comp_enthalpy_expr(temperature, "H2")
            - 0.5 * comp_enthalpy_expr(temperature, "O2")
        )
        dg = dh - temperature*ds
        dg_approx = pyo.value(-2*_constF*(1.253 - 0.00024516*temperature))

        print(f"{temperature}, {ds}, {dh}, {dg}, {dg_approx}")


    m.fs.fuel_electrode.qrate_per_area_gen_tbp.display()
    m.fs.oxygen_electrode.qrate_per_area_gen_tbp.display()
    m.fs.fuel_electrode.current_density.display()

    m.fs.fuel_electrode.qflux_tpb.display()
    m.fs.oxygen_electrode.qflux_tpb.display()
    m.fs.fuel_electrode.qflux_rxn_tpb.display()
    m.fs.oxygen_electrode.qflux_rxn_tpb.display()
    m.fs.fuel_electrode.qrate_per_area_gen_tbp.display()
    m.fs.oxygen_electrode.qrate_per_area_gen_tbp.display()


    m.fs.fuel_chan.temperature.display()
    m.fs.oxygen_chan.temperature.display()
    m.fs.fuel_electrode.temperature_x1.display()
    m.fs.oxygen_electrode.temperature_x1.display()

    dhfc = pyo.value(
        m.fs.fuel_chan.flow_area*(
            m.fs.fuel_chan.zflux_enth[0, m.fs.fuel_chan.izfaces.last()]
            - m.fs.fuel_chan.zflux_enth[0, m.fs.fuel_chan.izfaces.first()])
    )

    dhoc = pyo.value(
        m.fs.oxygen_chan.flow_area*(
            m.fs.oxygen_chan.zflux_enth[0, m.fs.oxygen_chan.izfaces.last()]
            - m.fs.oxygen_chan.zflux_enth[0, m.fs.oxygen_chan.izfaces.first()])
    )
    dmfc = pyo.value(
        m.fs.fuel_chan.flow_area*sum(
            m.fs.fuel_chan.zflux[0, m.fs.fuel_chan.izfaces.last(), i]*bin_diff_M[i]
            - m.fs.fuel_chan.zflux[0, m.fs.fuel_chan.izfaces.first(), i]*bin_diff_M[i]
            for i in m.fs.fuel_chan.config.comp_list
        )
    )
    dmoc = pyo.value(
        m.fs.oxygen_chan.flow_area*sum(
            m.fs.oxygen_chan.zflux[0, m.fs.oxygen_chan.izfaces.last(), i]*bin_diff_M[i]
            - m.fs.oxygen_chan.zflux[0, m.fs.oxygen_chan.izfaces.first(), i]*bin_diff_M[i]
            for i in m.fs.oxygen_chan.config.comp_list
        )
    )

    print(f"Mass Change: {dmoc + dmfc}")
    print(f"Enthalpy Change {dhoc + dhfc}")
    p = 0
    for iz in m.fs.oxygen_chan.iznodes:
        p += pyo.value(m.fs.potential_cell[0]*m.fs.fuel_electrode.current[0, iz])
    print(f"Total Electric Power: {-p}")

    m.fs.fuel_chan.temperature.display()
    m.fs.fuel_chan.temperature_el.display()
    m.fs.fuel_electrode.temperature_x0.display()

    m.fs.oxygen_chan.temperature.display()
    m.fs.oxygen_chan.temperature_el.display()
    m.fs.oxygen_electrode.temperature_x0.display()
