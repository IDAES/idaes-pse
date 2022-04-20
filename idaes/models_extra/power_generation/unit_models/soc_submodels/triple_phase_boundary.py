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

from pyomo.common.config import ConfigBlock, ConfigValue, In
import pyomo.environ as pyo


from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
from idaes.models_extra.power_generation.unit_models.soc_submodels.common import (
    _constR, _constF, _set_if_unfixed, _species_list, _element_list, _element_dict
)
import idaes.core.util.scaling as iscale
from idaes.core.util import get_solver

import idaes.logger as idaeslog

@declare_process_block_class("SocTriplePhaseBoundary")
class SocTriplePhaseBoundaryData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag",
            doc="No capacities or holdups, so no internal dynamics",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(domain=In([False]), default=False),
    )
    CONFIG.declare(
        "component_list",
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

    common._submodel_boilerplate_config(CONFIG)
    common._thermal_boundary_conditions_config(CONFIG, thin=True)
    common._material_boundary_conditions_config(CONFIG, thin=True)

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        tset = self.flowsheet().config.time
        # z coordinates for nodes and faces
        zfaces = self.config.cv_zfaces
        self.znodes = pyo.Set(
            initialize=[
                (zfaces[i] + zfaces[i + 1]) / 2.0 for i in range(len(zfaces) - 1)
            ]
        )
        comps = self.component_list = pyo.Set(
            initialize=self.config.component_list,
            ordered=True,
            doc="Set of all gas-phase components present in submodel",
        )
        self.tpb_stoich = copy.copy(self.config.tpb_stoich_dict)
        # TODO maybe let user specify inert species directly? Floating point
        # equalities make me nervous---Doug
        self.inert_component_list = pyo.Set(
            initialize=[j for j, coeff in self.tpb_stoich.items() if coeff == 0],
            ordered=True,
            doc="Set of components that do not react at triple phase boundary"
        )
        self.reacting_component_list =pyo.Set(
            initialize=[
                j for j, coeff in self.tpb_stoich.items() 
                if j not in self.inert_component_list
            ],
            ordered=True,
            doc="Set of components (gas-phase and solid) that react at triple "
            "phase boundary"
        )
        self.reacting_gas_list = pyo.Set(
            initialize=[j for j in comps if j not in self.inert_component_list],
            ordered=True,
            doc="Set of gas-phase components that react at triple phase boundary"
        )

        iznodes = self.iznodes = pyo.Set(initialize=range(1, len(self.znodes) + 1))

        common._submodel_boilerplate_create_if_none(self)
        common._create_thermal_boundary_conditions_if_none(self, thin=True)
        common._create_material_boundary_conditions_if_none(self, thin=True)

        common._create_if_none(
            self,
            "conc_ref",
            idx_set=(tset, iznodes, comps),
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
            self.reacting_gas_list,
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
            pressure_exponent = sum(nu_j[j] for j in b.reacting_gas_list)
            if abs(pressure_exponent) < 1e-6:
                out_expr = 0
            else:
                out_expr = -_constR * pressure_exponent * pyo.log(P / P_ref)
            return out_expr + (
                sum(nu_j[j] * common._comp_entropy_expr(T, j)
                    for j in b.reacting_component_list)
                - _constR
                * sum(
                    nu_j[j] * log_y_j[t, iz, j]
                    for j in b.reacting_gas_list
                    # TODO verify that excluding solids is correct
                )
            )

        @self.Expression(tset, iznodes)
        def dh_rxn(b, t, iz):
            T = b.temperature[t, iz]
            nu_j = b.tpb_stoich
            return sum(nu_j[j] * common._comp_enthalpy_expr(T, j) 
                       for j in b.reacting_component_list)

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
            for j in b.reacting_gas_list:
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
                denom = pyo.value(sum(self.conc[t, iz, j] for j in self.component_list))
                for j in self.component_list:
                    self.mole_frac_comp[t, iz, j].value = pyo.value(
                        self.conc[t, iz, j] / denom
                    )
                    self.log_mole_frac_comp[t, iz, j].value = pyo.value(
                        pyo.log(self.mole_frac_comp[t, iz, j])
                    )

        slvr = get_solver(solver, optarg)
        common._init_solve_block(self, slvr, solve_log)

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
        ssf = common._set_scaling_factor_if_none
        sgsf = common._set_and_get_scaling_factor
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
                for j in self.component_list:
                    # TODO Come back with good formulation for trace components
                    # and scale DConc and Cref
                    sy = sgsf(self.mole_frac_comp[t, iz, j], sy_def)
                    ssf(self.log_mole_frac_comp[t, iz, j], 1)
                    cst(self.mole_frac_comp_eqn[t, iz, j], sy)
                    cst(self.log_mole_frac_comp_eqn[t, iz, j], sy)
                    sxflux = sgsf(self.xflux[t, iz, j], 1e-2)
                    cst(self.xflux_eqn[t, iz, j], sxflux)