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

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.dae import DerivativeVar
import pyomo.environ as pyo


from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
from idaes.models_extra.power_generation.unit_models.soc_submodels.common import (
    _constR, _set_if_unfixed, _species_list, _element_list, _element_dict
)
import idaes.core.util.scaling as iscale
from idaes.core.util import get_solver

import idaes.logger as idaeslog

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

    common._submodel_boilerplate_config(CONFIG)
    common._thermal_boundary_conditions_config(CONFIG, thin=False)

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

        common._submodel_boilerplate_create_if_none(self)
        common._create_thermal_boundary_conditions_if_none(self, thin=False)

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
            return b.xfaces.at(ix + 1) - b.xfaces.at(ix)

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
            return common._interpolate_2D(
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
            return common._interpolate_2D(
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
        ssf = common._set_scaling_factor_if_none
        sgsf = common._set_and_get_scaling_factor
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