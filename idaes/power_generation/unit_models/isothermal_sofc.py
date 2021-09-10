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

import numpy as np
from itertools import combinations

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.dae import ContinuousSet, DerivativeVar
import pyomo.environ as pyo

from idaes.core.util.config import is_physical_parameter_block
from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.solvers import use_idaes_solver_configuration_defaults
import idaes
from idaes.core.util.math import safe_log
from idaes.core.util import get_solver
from pyomo.network import Port

# def _entropy_expr(temperature, x, comp_set):
#    return sum(_comp_entropy_expr(temperature, i)*x[i] + \
#        _constR*x[i]*safe_log(x[i]) for i in comp_set)

_constR = 8.3145  # J/mol/K or Pa*m3/K/mol
_constF = 96485  # C/mol

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
        "A": 18.563083,
        "B": 12.257357,
        "C": -2.859786,
        "D": 0.268238,
        "E": 1.977990,
        "F": -1.147438,
        "G": 156.288133,
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
    cm2_to_m2 = 0.0001
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
    t = temperature / 1000.0
    return 1000 * (
        d["A"] * t
        + d["B"] * t ** 2 / 2.0
        + d["C"] * t ** 3 / 3.0
        + d["D"] * t ** 4 / 4.0
        - d["E"] / t
        + d["F"]
    )  # J


def comp_int_energy_expr(temperature, comp):
    # ideal gas enthalpy
    d = h_params[comp]
    t = temperature / 1000.0
    return comp_enthalpy_expr(temperature, comp) - _constR * temperature


def comp_entropy_expr(temperature, comp):
    # ideal gas enthalpy
    d = h_params[comp]
    t = temperature / 1000.0
    return (
        d["A"] * pyo.log(t)
        + d["B"] * t
        + d["C"] * t ** 2 / 2.0
        + d["D"] * t ** 3 / 3.0
        - d["E"] / 2.0 / t ** 2
        + d["G"]
    )


@declare_process_block_class("IsothermalSofcChannel")
class IsothermalSofcChannelData(UnitModelBlockData):
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
        "zset", ConfigValue(description="Evenly spaced symetric CV boundary set")
    )
    CONFIG.declare(
        "comp_list",
        ConfigValue(default=["H2", "H2O"], description="List of components"),
    )
    CONFIG.declare(
        "length", ConfigValue(default=None, description="Lenght variable from parent")
    )
    CONFIG.declare(
        "width", ConfigValue(default=None, description="Width variable from parent")
    )

    def build(self):
        super().build()

        is_dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        zset = self.config.zset
        comps = self.config.comp_list

        self.izset = pyo.Set(initialize=range(0, len(zset)))
        self.izset_cv = pyo.Set(initialize=range(1, len(zset)))
        izset = self.izset
        izset_cv = self.izset_cv

        if self.config.length is None:
            self.length = pyo.Var(initialize=0.25)
        else:
            self.length = pyo.Reference(self.config.length)

        if self.config.width is None:
            self.width = pyo.Var(initialize=0.25)
        else:
            self.width = pyo.Reference(self.config.width)

        self.conc = pyo.Var(tset, izset, comps, doc="Component concentration", units=pyo.units.mol/pyo.units.m**3)
        self.pressure = pyo.Var(tset, izset, doc="Pressure", units=pyo.units.Pa)
        self.temperature = pyo.Var(tset, izset, doc="Temperature", units=pyo.units.K)
        self.thickness = pyo.Var(doc="Thickness in x direction", units=pyo.units.m)
        self.mole_frac_comp = pyo.Var(tset, izset, comps, doc="Component mole fraction")
        self.velocity = pyo.Var(tset, izset, doc="Gas velocity", units=pyo.units.m/pyo.units.s)
        self.flow_mol = pyo.Var(tset, izset, units=pyo.units.mol/pyo.units.s)
        self.xflux = pyo.Var(
            tset, izset, comps, doc="Component flux to electrode", initialize=0
        )
        if is_dynamic:
            self.dcdt = DerivativeVar(
                self.conc,
                wrt=tset,
                initialize=0,
                doc="Component concentration time derivative",
            )
        else:
            self.dcdt = 0

        @self.Expression()
        def flow_area(b):
            return b.thickness * b.width[None]

        @self.Expression(tset, izset)
        def volume_molar(b, t, iz):
            return _constR * b.temperature[t, iz] / b.pressure[t, iz]

        @self.Constraint(tset, izset)
        def flow_mol_eqn(b, t, iz):
            return (
                b.flow_mol[t, iz]
                == b.flow_area * b.velocity[t, iz] / b.volume_molar[t, iz]
            )

        @self.Constraint(tset, izset, comps)
        def conc_eqn(b, t, iz, i):
            return (
                b.conc[t, iz, i] * b.temperature[t, iz] * _constR
                == b.pressure[t, iz] * b.mole_frac_comp[t, iz, i]
            )

        @self.Constraint(tset, izset_cv)
        def mole_frac_eqn(b, t, iz):
            return 1 == sum(b.mole_frac_comp[t, iz, i] for i in comps)

        @self.Constraint(tset, izset)
        def p_eqn(b, t, iz):
            if iz == 0:
                return pyo.Constraint.Skip
            return b.pressure[t, iz] == b.pressure[t, 0]

        @self.Constraint(tset, izset_cv, comps)
        def mass_balance_eqn(b, t, iz, i):
            dz = (zset[iz] - zset[iz - 1]) * self.length[None]
            if t == tset.first() and is_dynamic:
                return pyo.Constraint.Skip
            else:
                return (
                    self.dcdt
                    == (
                        b.conc[t, iz - 1, i] * b.velocity[t, iz - 1]
                        - b.conc[t, iz, i] * b.velocity[t, iz]
                    )
                    / dz
                    - b.xflux[t, iz, i] / b.thickness
                )

    def initialize(
        self,
        outlvl=idaeslog.DEBUG,
        solver=None,
        optarg=None,
        xflux_guess=0,
        velocity_guess=None,
    ):
        # Get a reaconable value in for concentration, ...
        if velocity_guess is None:
            t0 = self.flowsheet().config.time.first()
            velocity_guess = pyo.value(
                self.flow_mol[t0, 0] * self.volume_molar[t0, 0] / self.flow_area
            )
        self.xflux.fix(xflux_guess)
        fix_flow = {}
        for i, v in self.flow_mol.items():
            fix_flow[i] = v.fixed
        self.flow_mol.unfix()
        self.velocity[:, 0].fix(velocity_guess)
        for t in self.flowsheet().config.time:
            for z in self.izset_cv:
                for i in self.config.comp_list:
                    self.mole_frac_comp[t, z, i].value = self.mole_frac_comp[
                        t, 0, i
                    ].value
                self.velocity[t, z].value = self.velocity[t, 0].value
                self.pressure[t, z].value = self.pressure[t, 0].value
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        slvr = get_solver(solver, optarg)
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = slvr.solve(self, tee=slc.tee)
        self.xflux.unfix()
        self.velocity.unfix()
        for i, v in fix_flow.items():
            if v:
                self.flow_mol[i].fix()

    def calculate_scaling_factors(self):
        # Base default scaling on typical conditions and dimensions
        def _set_default_factor(c, s):
            for i in c:
                if iscale.get_scaling_factor(c[i]) is None:
                    iscale.set_scaling_factor(c[i], s)

        _set_default_factor(self.pressure, 1e-5)
        _set_default_factor(self.temperature, 1e-2)
        _set_default_factor(self.thickness, 1e5)
        _set_default_factor(self.conc, 1)
        _set_default_factor(self.mole_frac_comp, 10)
        _set_default_factor(self.flow_mol, 1e4)
        _set_default_factor(self.velocity, 10)
        _set_default_factor(self.xflux, 100)
        iscale.propagate_indexed_component_scaling_factors(self)
        s_width = 1 / pyo.value(self.width[None])
        s_length = 1 / pyo.value(self.length[None])
        s_deltaz = len(self.izset_cv) * s_length

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

        for i in self.p_eqn:
            sP = iscale.get_scaling_factor(self.pressure[i])
            iscale.constraint_scaling_transform(self.p_eqn[i], sP)

        for t, iz, i in self.mass_balance_eqn:
            sc = iscale.get_scaling_factor(self.conc[t, iz, i])
            sv = iscale.get_scaling_factor(self.velocity[t, iz])
            iscale.constraint_scaling_transform(
                self.mass_balance_eqn[t, iz, i], sc * sv / s_deltaz
            )




@declare_process_block_class("IsothermalSofcElectrode")
class IsothermalSofcElectrodeData(UnitModelBlockData):
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
        "zset", ConfigValue(description="Evenly spaced symetric CV boundary set")
    )
    CONFIG.declare(
        "nx",
        ConfigValue(
            default=10, domain=int, description="Number of cells in x-direction"
        ),
    )
    CONFIG.declare(
        "comp_list",
        ConfigValue(default=["H2", "H2O"], description="List of components"),
    )
    CONFIG.declare(
        "stoic_dict", ConfigValue(default={}, description="List of components")
    )
    CONFIG.declare(
        "length", ConfigValue(default=None, description="Lenght variable from parent")
    )
    CONFIG.declare(
        "width", ConfigValue(default=None, description="Width variable from parent")
    )
    CONFIG.declare(
        "current", ConfigValue(default=None, description="Current variable from parent")
    )
    CONFIG.declare(
        "channel_conc",
        ConfigValue(default=None, description="Channel concentration variable"),
    )
    CONFIG.declare(
        "reversed",
        ConfigValue(default=False, description="Runs oposite of main z domain"),
    )

    def build(self):
        super().build()

        is_dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        zset = self.config.zset
        nx = self.config.nx
        comps = self.comps = self.config.comp_list

        self.xset = np.linspace(0, 1, nx + 1).tolist()
        self.izset = pyo.Set(initialize=range(0, len(zset)))
        self.ixset = pyo.Set(initialize=range(0, nx + 1))
        self.ixset_cv = pyo.Set(initialize=range(1, nx + 1))
        self.izset_cv = pyo.Set(initialize=range(1, len(zset)))
        xset = self.xset
        izset = self.izset
        ixset = self.ixset
        izset_cv = self.izset_cv
        ixset_cv = self.ixset_cv
        self.current = pyo.Reference(self.config.current)

        if self.config.length is None:
            self.length = pyo.Var(initialize=0.25)
        else:
            self.length = pyo.Reference(self.config.length)

        if self.config.width is None:
            self.width = pyo.Var(initialize=0.25)
        else:
            self.width = pyo.Reference(self.config.width)

        self.conc = pyo.Var(tset, izset_cv, ixset, comps)
        self.pressure = pyo.Var(tset, izset_cv, ixset)
        self.temperature = pyo.Var(tset, izset_cv, ixset)
        self.mole_frac_comp = pyo.Var(tset, izset_cv, ixset, comps)
        self.thickness = pyo.Var()
        self.porosity = pyo.Var(doc="Electrode porosity")
        self.tortuosity = pyo.Var(doc="Electrode tortuosity")
        self.k_res = pyo.Var()
        self.E_res = pyo.Var()
        if is_dynamic and has_holdup:
            self.dcdt = DerivativeVar(
                self.conc,
                wrt=tset,
                initialize=0,
                doc="Component concentration time derivative",
            )
        else:
            self.dcdt = 0

        @self.Constraint(tset, izset_cv, ixset, comps)
        def conc_eqn(b, t, iz, ix, i):
            return (
                b.conc[t, iz, ix, i] * b.temperature[t, iz, ix] * _constR
                == b.pressure[t, iz, ix] * b.mole_frac_comp[t, iz, ix, i]
            )

        @self.Constraint(tset, izset_cv, ixset)
        def mole_frac_eqn(b, t, iz, ix):
            return 1 == sum(b.mole_frac_comp[t, iz, ix, i] for i in comps)

        @self.Expression(tset, izset_cv, ixset, comps)
        def diff_coeff_im(b, t, iz, ix, i):
            T = b.temperature[t, iz, ix]
            P = b.pressure[t, iz, ix]
            x = b.mole_frac_comp
            bfun = binary_diffusion_coefficient_expr
            return (1.0 - x[t, iz, ix, i]) / sum(
                x[t, iz, ix, j] / bfun(T, P, i, j) for j in comps if i != j
            )

        @self.Expression(tset, izset_cv, ixset, comps)
        def diff_eff_coeff(b, t, iz, ix, i):
            return b.porosity / b.tortuosity * b.diff_coeff_im[t, iz, ix, i]

        @self.Expression(tset, izset_cv)
        def resistance(b, t, iz):
            T = b.temperature[t, iz, nx]
            A = self.width[None] * (zset[iz] - zset[iz - 1]) * self.length[None]
            return b.thickness * b.k_res * pyo.exp(b.E_res / T) / A

        @self.Expression(tset, izset_cv, ixset, comps)
        def dcdx(b, t, iz, ix, i):
            if ix == nx:
                c = self.config.stoic_dict[i]
                dz = (zset[iz] - zset[iz - 1]) * self.length[None]
                if c == 0:
                    return 0
                elif not self.config.reversed:
                    return (
                        c
                        * b.current[t, iz]
                        / _constF
                        / b.diff_eff_coeff[t, iz, ix, i]
                        / dz
                        / b.width[None]
                    )
                else:
                    return (
                        c
                        * b.current[t, 1 + izset_cv.last() - iz]
                        / _constF
                        / b.diff_eff_coeff[t, iz, ix, i]
                        / dz
                        / b.width[None]
                    )
            else:
                dx = (xset[ix + 1] - xset[ix]) * self.thickness
                return (b.conc[t, iz, ix + 1, i] - b.conc[t, iz, ix, i]) / dx

        @self.Expression(tset, izset_cv, ixset, comps)
        def dcdz(b, t, iz, ix, i):
            if iz == izset_cv.last():
                return 0
            else:
                dz = (zset[iz + 1] - zset[iz]) * self.length[None]
                return (b.conc[t, iz, ix, i] - b.conc[t, iz + 1, ix, i]) / dz

        @self.Constraint(tset, izset_cv, ixset, comps)
        def mass_balance_eqn(b, t, iz, ix, i):
            dz = (zset[iz] - zset[iz - 1]) * b.length[None]
            dx = (xset[ix] - xset[ix - 1]) * b.thickness
            if t == tset.first() and is_dynamic:
                return pyo.Constraint.Skip
            elif ix == 0:  # bound between the channel and electrode, same conc
                # This isn't a mass balance, x = 0 is a bound control volumes
                # indexing starts at 1. This is a boundary condition and it
                # does make the flux into the electrode the same as the flux
                # out of the channel.
                return b.conc[t, iz, 0, i] == self.config.channel_conc[t, iz, i]
            elif iz == 1:
                return self.dcdt == b.diff_eff_coeff[t, iz, ix, i] * (
                    -b.dcdz[t, iz, ix, i] / dz
                    + (  # nothing in at z = 0 (is wall)
                        b.dcdx[t, iz, ix - 1, i] - b.dcdx[t, iz, ix, i]
                    )
                    / dx
                )
            else:
                return self.dcdt == b.diff_eff_coeff[t, iz, ix, i] * (
                    (b.dcdz[t, iz - 1, ix, i] - b.dcdz[t, iz, ix, i]) / dz
                    + (b.dcdx[t, iz, ix - 1, i] - b.dcdx[t, iz, ix, i]) / dx
                )

    def initialize(
        self,
        outlvl=idaeslog.DEBUG,
        solver=None,
        optarg=None,
        comp_guess={},
        pressure_guess=1e5,
    ):
        for t in self.flowsheet().config.time:
            for x in self.ixset:
                for z in self.izset_cv:
                    for i in self.config.comp_list:
                        if i in comp_guess:
                            self.mole_frac_comp[t, z, x, i] = comp_guess[i]
                    self.pressure[t, z, x] = pressure_guess

    def calculate_scaling_factors(self):
        def _set_default_factor(c, s):
            for i in c:
                if iscale.get_scaling_factor(c[i]) is None:
                    iscale.set_scaling_factor(c[i], s)

        _set_default_factor(self.pressure, 1e-5)
        _set_default_factor(self.temperature, 1e-2)
        _set_default_factor(self.thickness, 1e5)
        _set_default_factor(self.conc, 1)
        _set_default_factor(self.mole_frac_comp, 10)
        # _set_default_factor(self.dcdx, 10)
        # _set_default_factor(self.dcdz, 10)
        iscale.propagate_indexed_component_scaling_factors(self)
        s_width = 1 / pyo.value(self.width[None])
        s_length = 1 / pyo.value(self.length[None])
        s_thickness = 1 / pyo.value(self.thickness)
        s_deltaz = len(self.izset_cv) * s_length
        s_deltax = len(self.ixset_cv) * s_thickness
        for i in self.conc_eqn:
            sp = iscale.get_scaling_factor(self.pressure[i[0], i[1], i[2]])
            sx = iscale.get_scaling_factor(self.mole_frac_comp[i])
            iscale.constraint_scaling_transform(self.conc_eqn[i], sp * sx)

        for i in self.mole_frac_eqn:
            iscale.constraint_scaling_transform(self.mole_frac_eqn[i], 10)

        for i in self.mass_balance_eqn:
            sc = iscale.get_scaling_factor(self.conc[i])
            sd = iscale.get_scaling_factor(self.dcdx[i], default=100)
            if i[2] == 0:  # x = 0
                iscale.constraint_scaling_transform(self.mass_balance_eqn[i], sc)
            else:
                iscale.constraint_scaling_transform(
                    self.mass_balance_eqn[i], sd / s_deltax
                )


@declare_process_block_class("IsothermalSofcElectrolyte")
class IsothermalSofcElectrolyteData(UnitModelBlockData):
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
        "zset", ConfigValue(description="Evenly spaced symetric CV boundary set")
    )
    CONFIG.declare(
        "length", ConfigValue(default=None, description="Lenght variable from parent")
    )
    CONFIG.declare(
        "width", ConfigValue(default=None, description="Width variable from parent")
    )

    def build(self):
        #
        super().build()
        tset = self.flowsheet().config.time
        zset = self.config.zset
        self.izset = pyo.Set(initialize=range(0, len(zset)))
        self.izset_cv = pyo.Set(initialize=range(1, len(zset)))
        self.thickness = pyo.Var()
        self.temperature = pyo.Var(tset, self.izset)


        if self.config.length is None:
            self.length = pyo.Var(initialize=0.25)
        else:
            self.length = pyo.Reference(self.config.length)

        if self.config.width is None:
            self.width = pyo.Var(initialize=0.25)
        else:
            self.width = pyo.Reference(self.config.width)

        self.k_res = pyo.Var(initialize=2.94e-5)
        self.E_res = pyo.Var(initialize=10350)

        @self.Expression(tset, self.izset_cv)
        def resistance(b, t, iz):
            T = b.temperature[t, iz]
            A = b.width[None] * (zset[iz] - zset[iz - 1]) * b.length[None]
            return b.thickness * b.k_res * pyo.exp(b.E_res / T) / A


@declare_process_block_class("IsothermalSofc")
class IsothermalSofcData(UnitModelBlockData):
    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([useDefault, True, False]),
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates wether this model will be dynamic,
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
        "nz",
        ConfigValue(
            default=20, domain=int, description="Number of cells in z-direction"
        ),
    )

    CONFIG.declare(
        "nxfe",
        ConfigValue(
            default=20,
            domain=int,
            description="Number of cells in x-direction of fuel electrode",
        ),
    )

    CONFIG.declare(
        "nxae",
        ConfigValue(
            default=20,
            domain=int,
            description="Number of cells in x-direction of air electrode",
        ),
    )

    CONFIG.declare(
        "fuel_side_comp_list",
        ConfigValue(default=["H2", "H2O"], description="List of components"),
    )

    CONFIG.declare(
        "air_side_comp_list",
        ConfigValue(default=["N2", "O2"], description="List of components"),
    )

    CONFIG.declare(
        "fuel_side_stoich",
        ConfigValue(default={"H2O": 0.5, "H2": -0.5}, description="List of components"),
    )

    CONFIG.declare(
        "air_side_stoich",
        ConfigValue(default={"O2": -0.25, "N2": 0}, description="List of components"),
    )

    CONFIG.declare(
        "soec",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="If False, SOFC model, if True, SOEC model",
        ),
    )

    def build(self):
        #
        super().build()
        # Parts:
        # fe = fuel electrode (anode)
        # ae = air electrode (cathode)
        # fc = fuel channel
        # ac = air channel
        # el = electrolyte

        # Get some config stuff
        tset = self.flowsheet().config.time
        is_dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup
        nz = self.config.nz
        nxfe = self.config.nxfe
        nxae = self.config.nxae

        # Set up space sets will use dimensionless lenghts from 0 to 1
        self.zset = np.linspace(0, 1, nz + 1).tolist()
        self.izset = pyo.Set(initialize=range(0, len(self.zset)))
        self.izset_cv = pyo.Set(initialize=range(1, len(self.zset)))

        # Easier aliases
        zset = self.zset
        izset = self.izset
        izset_cv = self.izset_cv

        # Variables
        self.length = pyo.Var(initialize=0.25)
        self.width = pyo.Var(initialize=0.25)
        self.h2_utilization = pyo.Var(tset, initialize=0.75)
        self.o2_stoichs = pyo.Var(tset, initialize=8.0)

        # Eletrochemical vars
        self.E_cell = pyo.Var(tset, units=pyo.units.V)
        self.current = pyo.Var(tset, izset_cv, units=pyo.units.A)
        self.total_current = pyo.Var(tset, units=pyo.units.A)
        self.eta_fe = pyo.Var(tset, izset_cv, initialize=2000 * _constR / _constF)
        self.eta_ae = pyo.Var(tset, izset_cv, initialize=2000 * _constR / _constF)
        self.k_fe = pyo.Var()
        self.k_ae = pyo.Var()
        self.alpha_fe = pyo.Var()
        self.alpha_ae = pyo.Var()
        self.eact_fe = pyo.Var()
        self.eact_ae = pyo.Var()
        self.heat_duty = pyo.Var(tset)

        # multiple cells ports
        self.n_cells = pyo.Var(initialize=20e6)
        self.n_cells.fix()
        self.mult_flow_mol_fc_inlet = pyo.Var(tset)
        self.mult_flow_mol_fc_outlet = pyo.Var(tset)
        self.mult_flow_mol_ac_inlet = pyo.Var(tset)
        self.mult_flow_mol_ac_outlet = pyo.Var(tset)


        self.el = IsothermalSofcElectrolyte(
            default={
                "dynamic": is_dynamic,
                "has_holdup": has_holdup,
                "zset": self.zset,
                "length": self.length,
                "width": self.width,
            }
        )

        self.fc = IsothermalSofcChannel(
            default={
                "dynamic": is_dynamic,
                "has_holdup": has_holdup,
                "zset": self.zset,
                "comp_list": self.config.fuel_side_comp_list,
                "length": self.length,
                "width": self.width,
            }
        )
        self.ac = IsothermalSofcChannel(
            default={
                "dynamic": is_dynamic,
                "has_holdup": has_holdup,
                "zset": self.zset,
                "comp_list": self.config.air_side_comp_list,
                "length": self.length,
                "width": self.width,
            }
        )
        self.fe = IsothermalSofcElectrode(
            default={
                "dynamic": is_dynamic,
                "has_holdup": has_holdup,
                "zset": self.zset,
                "nx": self.config.nxfe,
                "comp_list": self.config.fuel_side_comp_list,
                "stoic_dict": self.config.fuel_side_stoich,
                "current": self.current,
                "reversed": False,
                "length": self.length,
                "width": self.width,
                "channel_conc": self.fc.conc,
            }
        )
        self.ae = IsothermalSofcElectrode(
            default={
                "dynamic": is_dynamic,
                "has_holdup": has_holdup,
                "zset": self.zset,
                "nx": self.config.nxae,
                "comp_list": self.config.air_side_comp_list,
                "stoic_dict": self.config.air_side_stoich,
                "current": self.current,
                "reversed": True,
                "length": self.length,
                "width": self.width,
                "channel_conc": self.ac.conc,
            }
        )

        @self.Expression(tset, izset_cv)
        def E_nerst(b, t, iz):
            T = b.fe.temperature[t, iz, b.fe.ixset.last()]
            Pae = b.ae.pressure[t, 1 + nz - iz, b.ae.ixset.last()]
            yH2 = b.fe.mole_frac_comp[t, iz, b.fe.ixset.last(), "H2"]
            yO2 = b.ae.mole_frac_comp[t, 1 + nz - iz, b.ae.ixset.last(), "O2"]
            yH2O = b.fe.mole_frac_comp[t, iz, b.fe.ixset.last(), "H2O"]
            return (
                1.253
                - 0.00024516 * T
                + _constR
                * T
                / 2.0
                / _constF
                * pyo.log(yH2 / yH2O * (Pae * yO2 / 1e5) ** 0.5)
            )

        @self.Expression(tset, izset_cv)
        def eta_ohm(b, t, iz):
            return b.current[t, iz] * (
                b.el.resistance[t, iz] + b.fe.resistance[t, iz] + b.ae.resistance[t, 1 + nz - iz]
            )

        @self.Expression(tset, izset_cv)
        def fe_ecd(b, t, iz):
            k = self.k_fe
            E = self.eact_fe
            yh2 = self.fe.mole_frac_comp[t, iz, nxfe, "H2"]
            yh2o = self.fe.mole_frac_comp[t, iz, nxfe, "H2O"]
            A = self.width * (zset[iz] - zset[iz - 1]) * self.length * self.fe.porosity
            T = b.fe.temperature[t, iz, nxfe]
            return k * A * yh2 * yh2o * pyo.exp(-E / _constR / T)

        @self.Expression(tset, izset_cv)
        def eta_fe_expr(b, t, iz):
            T = b.fe.temperature[t, iz, nxfe]
            ecd = b.fe_ecd[t, iz]
            I = b.current[t, iz]
            alpha = b.alpha_fe
            return _constR * T / alpha / _constF * pyo.asinh(I / ecd / 2.0)

        @self.Expression(tset, izset_cv)
        def ae_ecd(b, t, iz):
            k = self.k_ae
            E = self.eact_ae
            yo2 = self.ae.mole_frac_comp[t, iz, nxae, "O2"]
            A = self.width * (zset[iz] - zset[iz - 1]) * self.length * self.ae.porosity
            T = b.ae.temperature[t, iz, nxae]
            return k * A * yo2 ** 0.25 * pyo.exp(-E / _constR / T)

        @self.Expression(tset, izset_cv)
        def eta_ae_expr(b, t, iz):
            T = b.ae.temperature[t, iz, nxfe]
            ecd = b.ae_ecd[t, iz]
            I = b.current[t, 1 + nz - iz]
            alpha = b.alpha_ae
            return _constR * T / alpha / _constF * pyo.asinh(I / ecd / 2.0)

        @self.Constraint(tset, izset_cv)
        def potential_eqn(b, t, iz):
            return b.E_cell[t] == b.E_nerst[t, iz] - (
                b.eta_ohm[t, iz] + b.eta_fe_expr[t, iz] + b.eta_ae_expr[t, 1 + nz - iz]
            )

        @self.Constraint(tset, izset_cv, self.fc.config.comp_list)
        def fc_xflux_eqn(b, t, iz, i):
            return b.fc.xflux[t, iz, i] == -(
                b.fe.diff_eff_coeff[t, iz, 0, i] * b.fe.dcdx[t, iz, 0, i]
            )

        @self.Constraint(tset, izset_cv, self.ac.config.comp_list)
        def ac_xflux_eqn(b, t, iz, i):
            return b.ac.xflux[t, iz, i] == -(
                b.ae.diff_eff_coeff[t, iz, 0, i] * b.ae.dcdx[t, iz, 0, i]
            )

        @self.Expression(tset)
        def total_current_expr(b, t):
            return sum(b.current[t, iz] for iz in izset_cv)

        @self.Constraint(tset)
        def total_current_eqn(b, t):
            return b.total_current[t] == b.total_current_expr[t]

        @self.Expression(tset)
        def power(b, t):
            return -b.E_cell[t] * b.total_current[t]

        @self.Expression(tset)
        def total_power(b, t):
            return -b.E_cell[t] * b.total_current[t] * b.n_cells

        @self.Expression(tset)
        def h2_utilization_expr(b, t):
            return (
                (
                    b.fc.flow_mol[t, 0] * b.fc.mole_frac_comp[t, 0, "H2"]
                    - b.fc.flow_mol[t, izset.last()] * b.fc.mole_frac_comp[t, izset.last(), "H2"]
                )
                / b.fc.flow_mol[t, 0]
                / b.fc.mole_frac_comp[t, 0, "H2"]
            )

        @self.Expression(tset)
        def h2_generation_expr(b, t):
            return (
                b.fc.flow_mol[t, izset.last()] * b.fc.mole_frac_comp[t, izset.last(), "H2"]
                - b.fc.flow_mol[t, 0] * b.fc.mole_frac_comp[t, 0, "H2"]
            )

        @self.Expression(tset)
        def power_per_h2_generation_expr(b, t):
            return b.power[t]/b.h2_generation_expr[t]

        if not self.config.soec:

            @self.Expression(tset)
            def o2_stoichs_expr(b, t):
                return (
                    2
                    * b.ac.flow_mol[t, 0]
                    * b.ac.mole_frac_comp[t, 0, "O2"]
                    / b.fc.flow_mol[t, 0]
                    / b.fc.mole_frac_comp[t, 0, "H2"]
                )

            @self.Constraint(tset)
            def h2_utilization_eqn(b, t):
                return b.h2_utilization[t] == b.h2_utilization_expr[t]

            @self.Constraint(tset)
            def o2_stoichs_eqn(b, t):
                return b.o2_stoichs[t] == b.o2_stoichs_expr[t]

        else:
            pass

        @self.Expression(tset)
        def h_fuel_in(b, t):
            return b.fc.flow_mol[t, 0] * sum(
                comp_enthalpy_expr(b.fc.temperature[t, 0], i)
                * b.fc.mole_frac_comp[t, 0, i]
                for i in self.fc.config.comp_list
            )

        @self.Expression(tset)
        def h_air_in(b, t):
            return b.ac.flow_mol[t, 0] * sum(
                comp_enthalpy_expr(b.ac.temperature[t, 0], i)
                * b.ac.mole_frac_comp[t, 0, i]
                for i in self.ac.config.comp_list
            )

        @self.Expression(tset)
        def h_fuel_out(b, t):
            return b.fc.flow_mol[t, izset.last()] * sum(
                comp_enthalpy_expr(b.fc.temperature[t, nz], i)
                * b.fc.mole_frac_comp[t, nz, i]
                for i in self.fc.config.comp_list
            )

        @self.Expression(tset)
        def h_air_out(b, t):
            return b.ac.flow_mol[t, izset.last()] * sum(
                comp_enthalpy_expr(b.ac.temperature[t, nz], i)
                * b.ac.mole_frac_comp[t, nz, i]
                for i in self.ac.config.comp_list
            )

        @self.Expression(tset)
        def deltah_therm(b, t):
            return b.h_fuel_out[t] + b.h_air_out[t] - b.h_fuel_in[t] - b.h_air_in[t]

        @self.Constraint(tset)
        def enth_bal(b, t):
            return  b.deltah_therm[t] + b.heat_duty[t] == b.power[t]

        @self.Constraint(tset)
        def mult_flow_mol_fc_inlet_eqn(b, t):
            return b.mult_flow_mol_fc_inlet[t] == b.n_cells * self.fc.flow_mol[t, 0]

        @self.Constraint(tset)
        def mult_flow_mol_fc_outlet_eqn(b, t):
            return b.mult_flow_mol_fc_outlet[t] == b.n_cells * self.fc.flow_mol[t, nz]

        @self.Constraint(tset)
        def mult_flow_mol_ac_inlet_eqn(b, t):
            return b.mult_flow_mol_ac_inlet[t] == b.n_cells * self.ac.flow_mol[t, 0]

        @self.Constraint(tset)
        def mult_flow_mol_ac_outlet_eqn(b, t):
            return b.mult_flow_mol_ac_outlet[t] == b.n_cells * self.ac.flow_mol[t, nz]

        self.inlet_ac_flow_mol_ref = pyo.Reference(self.ac.flow_mol[:, 0])
        self.inlet_ac_temperature_ref = pyo.Reference(self.ac.temperature[:, 0])
        self.inlet_ac_pressure_ref = pyo.Reference(self.ac.pressure[:, 0])
        self.inlet_ac_mole_frac_comp_ref = pyo.Reference(self.ac.mole_frac_comp[:, 0, :])
        self.inlet_ac = Port()
        self.inlet_ac.add(self.inlet_ac_flow_mol_ref, "flow_mol")
        self.inlet_ac.add(self.inlet_ac_temperature_ref, "temperature")
        self.inlet_ac.add(self.inlet_ac_pressure_ref, "pressure")
        self.inlet_ac.add(self.inlet_ac_mole_frac_comp_ref, "mole_frac_comp")
        self.inlet_ac_mult = Port()
        self.inlet_ac_mult.add(self.mult_flow_mol_ac_inlet, "flow_mol")
        self.inlet_ac_mult.add(self.inlet_ac_temperature_ref, "temperature")
        self.inlet_ac_mult.add(self.inlet_ac_pressure_ref, "pressure")
        self.inlet_ac_mult.add(self.inlet_ac_mole_frac_comp_ref, "mole_frac_comp")

        self.outlet_ac_flow_mol_ref = pyo.Reference(self.ac.flow_mol[:, nz])
        self.outlet_ac_temperature_ref = pyo.Reference(self.ac.temperature[:, nz])
        self.outlet_ac_pressure_ref = pyo.Reference(self.ac.pressure[:, nz])
        self.outlet_ac_mole_frac_comp_ref = pyo.Reference(self.ac.mole_frac_comp[:, nz, :])
        self.outlet_ac = Port()
        self.outlet_ac.add(self.outlet_ac_flow_mol_ref, "flow_mol")
        self.outlet_ac.add(self.outlet_ac_temperature_ref, "temperature")
        self.outlet_ac.add(self.outlet_ac_pressure_ref, "pressure")
        self.outlet_ac.add(self.outlet_ac_mole_frac_comp_ref, "mole_frac_comp")
        self.outlet_ac_mult = Port()
        self.outlet_ac_mult.add(self.mult_flow_mol_ac_outlet, "flow_mol")
        self.outlet_ac_mult.add(self.outlet_ac_temperature_ref, "temperature")
        self.outlet_ac_mult.add(self.outlet_ac_pressure_ref, "pressure")
        self.outlet_ac_mult.add(self.outlet_ac_mole_frac_comp_ref, "mole_frac_comp")

        self.inlet_fc_flow_mol_ref = pyo.Reference(self.fc.flow_mol[:, 0])
        self.inlet_fc_temperature_ref = pyo.Reference(self.fc.temperature[:, 0])
        self.inlet_fc_pressure_ref = pyo.Reference(self.fc.pressure[:, 0])
        self.inlet_fc_mole_frac_comp_ref = pyo.Reference(self.fc.mole_frac_comp[:, 0, :])
        self.inlet_fc = Port()
        self.inlet_fc.add(self.inlet_fc_flow_mol_ref, "flow_mol")
        self.inlet_fc.add(self.inlet_fc_temperature_ref, "temperature")
        self.inlet_fc.add(self.inlet_fc_pressure_ref, "pressure")
        self.inlet_fc.add(self.inlet_fc_mole_frac_comp_ref, "mole_frac_comp")
        self.inlet_fc_mult = Port()
        self.inlet_fc_mult.add(self.mult_flow_mol_fc_inlet, "flow_mol")
        self.inlet_fc_mult.add(self.inlet_fc_temperature_ref, "temperature")
        self.inlet_fc_mult.add(self.inlet_fc_pressure_ref, "pressure")
        self.inlet_fc_mult.add(self.inlet_fc_mole_frac_comp_ref, "mole_frac_comp")

        self.outlet_fc_flow_mol_ref = pyo.Reference(self.fc.flow_mol[:, nz])
        self.outlet_fc_temperature_ref = pyo.Reference(self.fc.temperature[:, nz])
        self.outlet_fc_pressure_ref = pyo.Reference(self.fc.pressure[:, nz])
        self.outlet_fc_mole_frac_comp_ref = pyo.Reference(self.fc.mole_frac_comp[:, nz, :])
        self.outlet_fc = Port()
        self.outlet_fc.add(self.outlet_fc_flow_mol_ref, "flow_mol")
        self.outlet_fc.add(self.outlet_fc_temperature_ref, "temperature")
        self.outlet_fc.add(self.outlet_fc_pressure_ref, "pressure")
        self.outlet_fc.add(self.outlet_fc_mole_frac_comp_ref, "mole_frac_comp")
        self.outlet_fc_mult = Port()
        self.outlet_fc_mult.add(self.mult_flow_mol_fc_outlet, "flow_mol")
        self.outlet_fc_mult.add(self.outlet_fc_temperature_ref, "temperature")
        self.outlet_fc_mult.add(self.outlet_fc_pressure_ref, "pressure")
        self.outlet_fc_mult.add(self.outlet_fc_mole_frac_comp_ref, "mole_frac_comp")

    def initialize(self, outlvl=idaeslog.DEBUG, solver=None, optarg=None, soec=False):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        t0 = self.flowsheet().time.first()
        self.fc.initialize(outlvl=outlvl, optarg=optarg)
        self.ac.initialize(outlvl=outlvl, optarg=optarg)
        fe_cg = {}
        for i in self.fc.config.comp_list:
            fe_cg[i] = self.fc.mole_frac_comp[t0, 0, i].value
        ae_cg = {}
        for i in self.ac.config.comp_list:
            ae_cg[i] = self.ac.mole_frac_comp[t0, 0, i].value
        self.fe.initialize(
            outlvl=outlvl,
            optarg=optarg,
            pressure_guess=self.fc.pressure[t0, 0].value,
            comp_guess=fe_cg,
        )
        self.ae.initialize(
            outlvl=outlvl,
            optarg=optarg,
            pressure_guess=self.ac.pressure[t0, 0].value,
            comp_guess=ae_cg,
        )
        slvr = get_solver(solver, optarg)
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = slvr.solve(self, tee=slc.tee)

    def calculate_scaling_factors(self):
        # For the most part, will try to set scaling factors based on
        # dimensions and units of measure, and leave little to the user
        s_width = 1 / pyo.value(self.width)
        s_length = 1 / pyo.value(self.length)
        s_deltaz = len(self.izset_cv) * s_length
        s_deltax_fe = len(self.fe.xset) / pyo.value(self.fe.thickness)
        s_deltax_ae = len(self.ae.xset) / pyo.value(self.ae.thickness)
        iscale.set_scaling_factor(self.width, s_width)
        iscale.set_scaling_factor(self.length, s_length)
        iscale.set_scaling_factor(self.h2_utilization, 1.0)
        iscale.set_scaling_factor(self.o2_stoichs, 1.0)
        iscale.set_scaling_factor(self.el.temperature, 1e-2)
        iscale.set_scaling_factor(self.E_cell, 1)
        iscale.set_scaling_factor(self.current, s_deltaz * s_width / 300.0)
        iscale.set_scaling_factor(self.eta_fe, 100)
        iscale.set_scaling_factor(self.eta_ae, 100)
        iscale.set_scaling_factor(self.mult_flow_mol_fc_inlet, 1e-3)
        iscale.set_scaling_factor(self.mult_flow_mol_ac_inlet, 1e-3)
        iscale.set_scaling_factor(self.mult_flow_mol_fc_outlet, 1e-3)
        iscale.set_scaling_factor(self.mult_flow_mol_ac_outlet, 1e-3)

        for i, c in self.ac_xflux_eqn.items():
            s = iscale.get_scaling_factor(self.ac.xflux[i])
            iscale.constraint_scaling_transform(c, s)

        for i, c in self.fc_xflux_eqn.items():
            s = iscale.get_scaling_factor(self.fc.xflux[i])
            iscale.constraint_scaling_transform(c, s)

        for i, c in self.mult_flow_mol_fc_inlet_eqn.items():
            s = iscale.get_scaling_factor(self.mult_flow_mol_fc_inlet[i])
            iscale.constraint_scaling_transform(c, s)

        for i, c in self.mult_flow_mol_ac_inlet_eqn.items():
            s = iscale.get_scaling_factor(self.mult_flow_mol_ac_inlet[i])
            iscale.constraint_scaling_transform(c, s)

        for i, c in self.mult_flow_mol_fc_outlet_eqn.items():
            s = iscale.get_scaling_factor(self.mult_flow_mol_fc_outlet[i])
            iscale.constraint_scaling_transform(c, s)

        for i, c in self.mult_flow_mol_ac_outlet_eqn.items():
            s = iscale.get_scaling_factor(self.mult_flow_mol_ac_outlet[i])
            iscale.constraint_scaling_transform(c, s)


def sofc_example():
    from idaes.core import FlowsheetBlock
    from idaes.power_generation.properties.natural_gas_PR import get_prop
    import numpy as np
    import matplotlib.pyplot as plt

    nz = 20
    nxfe = 10
    nxae = 10
    dynamic = False
    check_scaling = False

    m = pyo.ConcreteModel("Test SOFC")
    if dynamic:
        fsc = {"dynamic": dynamic, "time_set": [0, 25], "time_units": pyo.units.s}
    else:
        fsc = {"dynamic": dynamic}
    m.fs = FlowsheetBlock(default=fsc)
    m.fs.sofc = IsothermalSofc(
        default={"dynamic": dynamic, "nz": nz, "nxfe": nxfe, "nxae": nxae}
    )
    if dynamic:
        disc = pyo.TransformationFactory("dae.finite_difference")
        disc.apply_to(m, nfe=5, wrt=m.fs.time, scheme="BACKWARD")

    # Scaling
    for i in m.fs.sofc.ac.xflux:
        if i[2] == "N2":
            iscale.set_scaling_factor(m.fs.sofc.ac.xflux[i], 1000)

    # Model inputs
    m.fs.sofc.E_cell[0].fix(0.8)
    m.fs.sofc.el.thickness.fix(8e-6)
    m.fs.sofc.fe.thickness.fix(1e-3)
    m.fs.sofc.ae.thickness.fix(20e-6)
    m.fs.sofc.length.fix(0.05)
    m.fs.sofc.width.fix(0.05)
    m.fs.sofc.k_ae.fix(1e10)
    m.fs.sofc.eact_ae.fix(120000)
    m.fs.sofc.alpha_ae.fix(0.5)
    m.fs.sofc.k_fe.fix(6.1e9)
    m.fs.sofc.eact_fe.fix(110000)
    m.fs.sofc.alpha_fe.fix(0.5)
    m.fs.sofc.fe.k_res.fix(2.98e-5)
    m.fs.sofc.fe.E_res.fix(-1392)
    m.fs.sofc.ae.k_res.fix(8.114e-5)
    m.fs.sofc.ae.E_res.fix(600)
    m.fs.sofc.el.k_res.fix(2.94e-5)
    m.fs.sofc.el.E_res.fix(10350)
    m.fs.sofc.fc.thickness.fix(0.002)
    m.fs.sofc.ac.thickness.fix(0.002)
    m.fs.sofc.fe.porosity.fix(0.48)
    m.fs.sofc.fe.tortuosity.fix(5.4)
    m.fs.sofc.ae.porosity.fix(0.48)
    m.fs.sofc.ae.tortuosity.fix(5.4)
    m.fs.sofc.o2_stoichs.fix(8.0 * 0.75)
    m.fs.sofc.h2_utilization.fix(0.75)
    temperature = 1073.15
    m.fs.sofc.el.temperature.fix(temperature)
    m.fs.sofc.fc.temperature.fix(temperature)
    m.fs.sofc.ac.temperature.fix(temperature)
    m.fs.sofc.fe.temperature.fix(temperature)
    m.fs.sofc.ae.temperature.fix(temperature)

    m.fs.sofc.fc.pressure[:, 0].fix(1e5)
    m.fs.sofc.fc.flow_mol[:, 0] = 9e-5 # don't fix calculated for utilization
    m.fs.sofc.fc.mole_frac_comp[:, 0, "H2O"].fix(0.03)
    m.fs.sofc.fc.mole_frac_comp[:, 0, "H2"].fix(0.97)

    m.fs.sofc.ac.pressure[:, 0].fix(1e5)
    m.fs.sofc.ac.flow_mol[:, 0] = 1.2e-3 # calculated from stoichs
    m.fs.sofc.ac.mole_frac_comp[:, 0, "O2"].fix(0.21)
    m.fs.sofc.ac.mole_frac_comp[:, 0, "N2"].fix(0.79)

    iscale.calculate_scaling_factors(m)
    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    idaes.cfg.ipopt["options"]["tol"] = 1e-8
    # due to a lot of component mole fractions being on their lower bound of 0
    # bound push result in much longer solve times, so set it low.
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-10
    idaes.cfg.ipopt["options"]["max_iter"] = 150
    iscale.calculate_scaling_factors(m)

    print("Init")
    m.fs.sofc.initialize()

    solver = pyo.SolverFactory("ipopt")
    # solver.options["halt_on_ampl_error"] = "yes"
    # solver.solve(m, tee=True, symbolic_solver_labels=True)
    solver.solve(m, tee=True)

    print(f"Potential: {pyo.value(m.fs.sofc.E_cell[0])} V")
    print(f"Current: {pyo.value(m.fs.sofc.total_current[0])} A")
    print(
        f"Current Density: {pyo.value(m.fs.sofc.total_current[0]/m.fs.sofc.length/m.fs.sofc.width)} A/m^2"
    )
    print(f"Power: {pyo.value(m.fs.sofc.power[0])} W")
    print(f"Heat duty: {pyo.value(m.fs.sofc.heat_duty[0])} W")
    print(f"Cell temperature: {pyo.value(m.fs.sofc.ac.temperature[0, nz]):.1f}K")
    print(f"Fuel inlet enthalpy: {pyo.value(m.fs.sofc.h_fuel_in[0])} W")
    print(f"Fuel outlet enthalpy: {pyo.value(m.fs.sofc.h_fuel_out[0])} W")
    print(f"Air inlet enthalpy: {pyo.value(m.fs.sofc.h_air_in[0])} W")
    print(f"Air outlet enthalpy: {pyo.value(m.fs.sofc.h_air_out[0])} W")
    print(f"DeltaH_therm: {pyo.value(m.fs.sofc.deltah_therm[0])} W")
    print(f"Flow_mol FC inlet: {pyo.value(m.fs.sofc.fc.flow_mol[0, 0])} mol/s")
    print(f"Flow_mol FC outlet: {pyo.value(m.fs.sofc.fc.flow_mol[0, nz])} mol/s")
    print(f"Flow_mol AC inlet: {pyo.value(m.fs.sofc.ac.flow_mol[0, 0])} mol/s")
    print(f"Flow_mol AC outlet: {pyo.value(m.fs.sofc.ac.flow_mol[0, nz])} mol/s")
    print(f"H2 Utilization: {pyo.value(m.fs.sofc.h2_utilization[0])}")
    print(f"O2 Stoichs: {pyo.value(m.fs.sofc.o2_stoichs[0])}")
    z = [
        (m.fs.sofc.zset[iz] + m.fs.sofc.zset[iz - 1]) / 2.0 for iz in m.fs.sofc.izset_cv
    ]
    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(m.fs.sofc.fc.mole_frac_comp[t, iz, "H2"])
            for iz in m.fs.sofc.izset_cv
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(
        f"Fuel Channel H$_2$ (E$_{{cell}}$ = {pyo.value(m.fs.sofc.E_cell[0]):.4f}V)"
    )
    plt.xlabel("z/L")
    plt.ylabel("$x_{H_2}$")
    plt.legend()
    plt.show()

    z = [
        (m.fs.sofc.zset[iz] + m.fs.sofc.zset[iz - 1]) / 2.0 for iz in m.fs.sofc.izset_cv
    ]
    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(m.fs.sofc.fc.mole_frac_comp[t, iz, "H2O"])
            for iz in m.fs.sofc.izset_cv
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(
        f"Fuel Channel H$_2$O (E$_{{cell}}$ = {pyo.value(m.fs.sofc.E_cell[0]):.4f}V)"
    )
    plt.xlabel("z/L")
    plt.ylabel("$x_{H_2O}$")
    plt.legend()
    plt.show()

    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(
                m.fs.sofc.current[t, iz]
                / m.fs.sofc.length
                / m.fs.sofc.width
                / (m.fs.sofc.zset[iz] - m.fs.sofc.zset[iz - 1])
            )
            for iz in m.fs.sofc.izset_cv
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(f"Current Density (E$_{{cell}}$ = {pyo.value(m.fs.sofc.E_cell[0])}V)")
    plt.xlabel("z/L")
    plt.ylabel("Current (A/m$^2$)")
    plt.legend()
    plt.show()

    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [pyo.value(m.fs.sofc.E_nerst[t, iz]) for iz in m.fs.sofc.izset_cv]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(f"E$_{{nerst}}$ (E$_{{cell}}$ = {pyo.value(m.fs.sofc.E_cell[0])}V)")
    plt.xlabel("z/L")
    plt.ylabel("E$_{nerst}$ (V)")
    plt.legend()
    plt.show()

    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(m.fs.sofc.ac.mole_frac_comp[t, iz, "O2"])
            for iz in reversed(m.fs.sofc.izset_cv)
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(f"Air Channel O$_2$ (E$_{{cell}}$ = {pyo.value(m.fs.sofc.E_cell[0])}V)")
    plt.xlabel("z/L")
    plt.ylabel("$x_{O_2}$")
    plt.legend()
    plt.show()

    """
    for ix in m.fs.sofc.fe_ixset:

        v = {t:None for t in m.fs.time}
        for t in v:
            v[t] = [pyo.value(m.fs.sofc.fe_conc[t, iz, ix, "H2"]) for iz in m.fs.sofc.izset_cv]

        for t in v:
            plt.plot(m.fs.sofc.zset[1:], v[t], label=f"t = {t}")

        plt.title("Fuel Electrode H2")
        plt.legend()
        plt.show()

    for ix in m.fs.sofc.ae_ixset:

        v = {t:None for t in m.fs.time}
        for t in v:
            v[t] = [pyo.value(m.fs.sofc.ae_conc[t, iz, ix, "O2"]) for iz in m.fs.sofc.izset_cv]

        for t in v:
            plt.plot(m.fs.sofc.zset[1:], v[t], label=f"t = {t}")

        plt.title("Air Electrode O2")
        plt.legend()
        plt.show()
    """
    m.fs.sofc.k_ae.unfix()
    m.obj = pyo.Objective(
        expr=(m.fs.sofc.total_current[0] / m.fs.sofc.length / m.fs.sofc.width - 6000)
        ** 2
    )
    solver.solve(m, tee=True)
    m.fs.sofc.k_ae.fix()

    if check_scaling:
        jac, nlp = iscale.get_jacobian(m, scaled=True)
        print("Extreme Jacobian entries:")
        for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, large=100):
            print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]")
        print("Unscaled constraints:")
        for c in iscale.unscaled_constraints_generator(m):
            print(f"    {c}")
        print("Scaled constraints by factor:")
        for c, s in iscale.constraints_with_scale_factor_generator(m):
            print(f"    {c}, {s}")
        print("Badly scaled variables:")
        for v, sv in iscale.badly_scaled_var_generator(
            m, large=1e2, small=1e-2, zero=1e-12
        ):
            print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}")
        print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")

    """

    dat = {
        800:{"e":[], "i":[]},
        750:{"e":[], "i":[]},
        700:{"e":[], "i":[]},
        650:{"e":[], "i":[]},
    }
    et = {
        800:[0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.91],
        750:[0.93, 0.92, 0.91, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50],
        700:[0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.91, 0.92, 0.93, 0.94],
        650:[0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20],
    }
    for T in reversed(sorted(dat)):
        for e in et[T]:
            m.fs.sofc.E_cell[0].fix(e)
            m.fs.sofc.fc.temperature.fix(T+273.15)
            m.fs.sofc.ac.temperature.fix(T+273.15)
            m.fs.sofc.fe.temperature.fix(T+273.15)
            m.fs.sofc.ae.temperature.fix(T+273.15)
            m.fs.sofc.el.temperature.fix(T+273.15)
            print("\n\nStart")
            print(T)
            print(e)
            solver.solve(m, tee=True)
            dat[T]["e"].append(e)
            dat[T]["i"].append(pyo.value(m.fs.sofc.total_current[0]/m.fs.sofc.length/m.fs.sofc.width))

    plt.plot(dat[800]["i"], dat[800]["e"], label="T = 800C")
    plt.plot(dat[750]["i"], dat[750]["e"], label="T = 750C")
    plt.plot(dat[700]["i"], dat[700]["e"], label="T = 700C")
    plt.plot(dat[650]["i"], dat[650]["e"], label="T = 650C")
    plt.legend()
    plt.xlabel("Current Density (A/m$^2$)")
    plt.ylabel("Voltage (V)")
    plt.ylim([0, 1.2])
    plt.xlim([0, 10000])
    plt.grid(axis='x', color='0.95')
    plt.grid(axis='y', color='0.95')
    plt.show()

    print(f"k_ae = {pyo.value(m.fs.sofc.k_ae):.4e}")
    print(f"k_fe = {pyo.value(m.fs.sofc.k_fe):.4e}")
    print(f"t, h, u, s")
    for t in [800, 900, 1000, 1100, 1200, 1300]:
        h = pyo.value(comp_enthalpy_expr(t, "O2"))
        s = pyo.value(comp_entropy_expr(t, "O2"))
        u = pyo.value(comp_int_energy_expr(t, "O2"))
        print(f"{t}, {h}, {u}, {s}")
    """
    return m


def soec_example():
    from idaes.core import FlowsheetBlock
    from idaes.power_generation.properties.natural_gas_PR import get_prop
    import numpy as np
    import matplotlib.pyplot as plt

    nz = 20
    nxfe = 10
    nxae = 10
    dynamic = False
    check_scaling = False

    m = pyo.ConcreteModel("Test SOEC")
    if dynamic:
        fsc = {"dynamic": dynamic, "time_set": [0, 25], "time_units": pyo.units.s}
    else:
        fsc = {"dynamic": dynamic}
    m.fs = FlowsheetBlock(default=fsc)
    m.fs.soec = IsothermalSofc(
        default={
            "dynamic": dynamic,
            "nz": nz,
            "nxfe": nxfe,
            "nxae": nxae,
            "soec": True,
            "air_side_comp_list": ["H2O", "O2"],
            "fuel_side_comp_list": ["H2O", "H2"],
            "air_side_stoich": {"H2O": 0, "O2": -0.25},
        }
    )
    if dynamic:
        disc = pyo.TransformationFactory("dae.finite_difference")
        disc.apply_to(m, nfe=5, wrt=m.fs.time, scheme="BACKWARD")

    # Model inputs
    m.fs.soec.E_cell[0].fix(1.287)
    m.fs.soec.el.thickness.fix(8e-6)
    m.fs.soec.fe.thickness.fix(1e-3)
    m.fs.soec.ae.thickness.fix(20e-6)
    m.fs.soec.length.fix(0.05)
    m.fs.soec.width.fix(0.05)
    m.fs.soec.k_ae.fix(8.7e7*3)
    m.fs.soec.eact_ae.fix(120000)
    m.fs.soec.alpha_ae.fix(0.5)
    m.fs.soec.k_fe.fix(1.35e10)
    m.fs.soec.eact_fe.fix(110000)
    m.fs.soec.alpha_fe.fix(0.4)
    m.fs.soec.fe.k_res.fix(2.98e-5)
    m.fs.soec.fe.E_res.fix(-1392)
    m.fs.soec.ae.k_res.fix(8.114e-5)
    m.fs.soec.ae.E_res.fix(600)
    m.fs.soec.el.k_res.fix(2.94e-5)
    m.fs.soec.el.E_res.fix(10350)
    m.fs.soec.fc.thickness.fix(0.002)
    m.fs.soec.ac.thickness.fix(0.002)
    m.fs.soec.fe.porosity.fix(0.48)
    m.fs.soec.fe.tortuosity.fix(5.4)
    m.fs.soec.ae.porosity.fix(0.48)
    m.fs.soec.ae.tortuosity.fix(5.4)
    temperature = 1073.15
    m.fs.soec.el.temperature.fix(temperature)
    m.fs.soec.fc.temperature.fix(temperature)
    m.fs.soec.ac.temperature.fix(temperature)
    m.fs.soec.fe.temperature.fix(temperature)
    m.fs.soec.ae.temperature.fix(temperature)

    m.fs.soec.fc.pressure[:, 0].fix(1e5)
    m.fs.soec.fc.flow_mol[:, 0].fix(1e-4)
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2O"].fix(0.7)
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2"].fix(0.30)

    m.fs.soec.ac.pressure[:, 0].fix(1e5)
    m.fs.soec.ac.flow_mol[:, 0].fix(1e-4)
    m.fs.soec.ac.mole_frac_comp[:, 0, "O2"].fix(0.1)
    m.fs.soec.ac.mole_frac_comp[:, 0, "H2O"].fix(0.9)

    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    idaes.cfg.ipopt["options"]["tol"] = 1e-8
    # due to a lot of component mole fractions being on their lower bound of 0
    # bound push result in much longer solve times, so set it low.
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-10
    idaes.cfg.ipopt["options"]["max_iter"] = 150
    iscale.calculate_scaling_factors(m)

    m.fs.soec.initialize()

    solver = pyo.SolverFactory("ipopt")
    # solver.options["halt_on_ampl_error"] = "yes"
    # solver.solve(m, tee=True, symbolic_solver_labels=True)
    solver.solve(m, tee=True)

    @m.fs.Constraint(m.fs.time)
    def heat_duty_zero_eqn(b, t):
        return b.soec.heat_duty[t] == 0
    m.fs.soec.E_cell.unfix()
    solver.solve(m, tee=True)

    z = [
        (m.fs.soec.zset[iz] + m.fs.soec.zset[iz - 1]) / 2.0 for iz in m.fs.soec.izset_cv
    ]
    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(m.fs.soec.fc.mole_frac_comp[t, iz, "H2"])
            for iz in m.fs.soec.izset_cv
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(
        f"Fuel Channel H$_2$ (E$_{{cell}}$ = {pyo.value(m.fs.soec.E_cell[0]):.4f}V)"
    )
    plt.xlabel("z/L")
    plt.ylabel("$x_{H_2}$")
    plt.legend()
    plt.show()

    z = [
        (m.fs.soec.zset[iz] + m.fs.soec.zset[iz - 1]) / 2.0 for iz in m.fs.soec.izset_cv
    ]
    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(m.fs.soec.fc.mole_frac_comp[t, iz, "H2O"])
            for iz in m.fs.soec.izset_cv
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(
        f"Fuel Channel H$_2$O (E$_{{cell}}$ = {pyo.value(m.fs.soec.E_cell[0]):.4f}V)"
    )
    plt.xlabel("z/L")
    plt.ylabel("$x_{H_2O}$")
    plt.legend()
    plt.show()

    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(
                m.fs.soec.current[t, iz]
                / m.fs.soec.length
                / m.fs.soec.width
                / (m.fs.soec.zset[iz] - m.fs.soec.zset[iz - 1])
            )
            for iz in m.fs.soec.izset_cv
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(f"Current Density (E$_{{cell}}$ = {pyo.value(m.fs.soec.E_cell[0])}V)")
    plt.xlabel("z/L")
    plt.ylabel("Current (A/m$^2$)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [pyo.value(m.fs.soec.E_nerst[t, iz]) for iz in m.fs.soec.izset_cv]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(f"E$_{{nerst}}$ (E$_{{cell}}$ = {pyo.value(m.fs.soec.E_cell[0])}V)")
    plt.xlabel("z/L")
    plt.ylabel("E$_{nerst}$ (V)")
    plt.legend()
    plt.show()

    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(m.fs.soec.ac.mole_frac_comp[t, iz, "O2"])
            for iz in reversed(m.fs.soec.izset_cv)
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(f"Air Channel O$_2$ (E$_{{cell}}$ = {pyo.value(m.fs.soec.E_cell[0])}V)")
    plt.xlabel("z/L")
    plt.ylabel("$x_{O_2}$")
    plt.legend()
    plt.show()

    print(f"Potential: {pyo.value(m.fs.soec.E_cell[0])} V")
    print(f"Current: {pyo.value(m.fs.soec.total_current[0])} A")
    print(
        f"Current Density: {pyo.value(m.fs.soec.total_current[0]/m.fs.soec.length/m.fs.soec.width)} A/m^2"
    )
    print(f"Power: {pyo.value(m.fs.soec.power[0])} W")
    print(f"H2 generated: {pyo.value(m.fs.soec.h2_generation_expr[0])} mol/s")
    print(f"Power/H2 generated: {pyo.value(m.fs.soec.power_per_h2_generation_expr[0])} J/mol H2")
    print(f"Power/H2 generated: {pyo.value(m.fs.soec.power_per_h2_generation_expr[0])/2e-3*1e-6} MJ/kg H2")
    print(f"Heat duty: {pyo.value(m.fs.soec.heat_duty[0])} W")
    print(f"Cell temperature: {pyo.value(m.fs.soec.ac.temperature[t,nz]):.1f} K")
    print(f"Fuel inlet enthalpy: {pyo.value(m.fs.soec.h_fuel_in[0])} W")
    print(f"Fuel outlet enthalpy: {pyo.value(m.fs.soec.h_fuel_out[0])} W")
    print(f"Air inlet enthalpy: {pyo.value(m.fs.soec.h_air_in[0])} W")
    print(f"Air outlet enthalpy: {pyo.value(m.fs.soec.h_air_out[0])} W")
    print(f"DeltaH_therm: {pyo.value(m.fs.soec.deltah_therm[0])} W")
    print(f"Flow_mol FC inlet: {pyo.value(m.fs.soec.fc.flow_mol[0, 0])} mol/s")
    print(f"Flow_mol FC outlet: {pyo.value(m.fs.soec.fc.flow_mol[0, nz])} mol/s")
    print(f"Flow_mol AC inlet: {pyo.value(m.fs.soec.ac.flow_mol[0, 0])} mol/s")
    print(f"Flow_mol AC outlet: {pyo.value(m.fs.soec.ac.flow_mol[0, nz])} mol/s")

    m.fs.soec.inlet_ac.display()
    m.fs.soec.outlet_ac.display()
    m.fs.soec.inlet_fc.display()
    m.fs.soec.outlet_fc.display()
    return m


if __name__ == "__main__":
    m = soec_example()
