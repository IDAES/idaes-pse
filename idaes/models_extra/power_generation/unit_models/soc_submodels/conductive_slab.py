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
Two-dimensional model of a slab of some conductive material. Presently it is
used to for the SOC electrolyte, but it can also be used for the interconnect.
It is discretized using a finite volume method with centered differences for
interpolation, and contains only thermal transport and resistive heating terms.

Boundary variables:
    - ``temperature_deviation_x0[t, iz]``
    - ``temperature_deviation_x1[t, iz]``
    - ``heat_flux_x0[t, iz]``
    - ``heat_flux_x1[t, iz]``

Instances of ``Var`` that must be fixed:
    - ``length_x``: Thickness of conductive slab
    - ``heat_capacity``: Specific heat capacity of slab on a mass-basis (i.e., in J/(kg * K))
    - ``density``: Mass density of slab
    - ``thermal_conductivity``: Thermal conductivity of slab
    - ``resistivity_log_preexponential_factor``: Natural logarithm of resistivity preexponential factor in ohm * m
    - ``resistivity_thermal_exponent_dividend``: Parameter divided by temperature in resistivity equation, in K.
      Would be something like (reduced) activation energy, but it can be both negative and positive.
"""
__author__ = "John Eslick, Douglas Allan"

from pyomo.common.config import ConfigValue
from pyomo.dae import DerivativeVar
import pyomo.environ as pyo


from idaes.core import declare_process_block_class, UnitModelBlockData
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog


@declare_process_block_class("SocConductiveSlab")
class SocConductiveSlabData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "control_volume_xfaces",
        ConfigValue(
            description="List containing coordinates of control volume faces "
            "in x direction. Coordinates must start with zero, be strictly "
            "increasing, and end with one"
        ),
    )

    common._submodel_boilerplate_config(CONFIG)
    common._thermal_boundary_conditions_config(CONFIG, thin=False)

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        # Set up node and face sets and get integer indices for them
        izfaces, iznodes = common._face_initializer(
            self, self.config.control_volume_zfaces, "z"
        )
        ixfaces, ixnodes = common._face_initializer(
            self, self.config.control_volume_xfaces, "x"
        )
        # Channel thickness AKA length in the x direction is specific to the
        # channel so local variable here is the only option
        self.length_x = pyo.Var(
            doc="Thickness of the electrode (x-direction)",
            units=pyo.units.m,
        )

        common._submodel_boilerplate_create_if_none(self)
        common._create_thermal_boundary_conditions_if_none(self, thin=False)

        self.temperature_deviation_x = pyo.Var(
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
        self.heat_capacity = pyo.Var(
            initialize=200,
            doc="Heat capacity of slab (mass basis)",
            units=pyo.units.J / pyo.units.kg / pyo.units.K,
        )
        self.density = pyo.Var(
            initialize=1000,
            doc="Density of slab",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.thermal_conductivity = pyo.Var(
            initialize=80,
            doc="Thermal conductivity of slab",
            units=pyo.units.W / pyo.units.m / pyo.units.K,
        )

        @self.Expression(tset, ixnodes, iznodes)
        def temperature(b, t, ix, iz):
            if b.config.include_temperature_x_thermo:
                return b.temperature_z[t, iz] + b.temperature_deviation_x[t, ix, iz]
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
                nodes=b.xnodes,
                faces=b.xfaces,
                phi_func=lambda ixf: b.temperature_deviation_x[t, ixf, iz] / b.length_x,
                phi_bound_0=(
                    b.temperature_deviation_x0[t, iz]
                    - b.temperature_deviation_x[t, ixnodes.first(), iz]
                )
                / (b.xfaces.first() - b.xnodes.first())
                / b.length_x,
                phi_bound_1=(
                    b.temperature_deviation_x[t, ixnodes.last(), iz]
                    - b.temperature_deviation_x1[t, iz]
                )
                / (b.xnodes.last() - b.xfaces.last())
                / b.length_x,
                derivative=True,
            )

        @self.Expression(tset, ixnodes, izfaces)
        def dTdz(b, t, ix, iz):
            return common._interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=b.znodes,
                faces=b.zfaces,
                phi_func=lambda izf: b.temperature[t, ix, izf] / b.length_z[None],
                phi_bound_0=0,
                phi_bound_1=0,
                derivative=True,
            )

        @self.Expression(tset, ixfaces, iznodes)
        def heat_flux_x(b, t, ix, iz):
            return -b.thermal_conductivity * b.dTdx[t, ix, iz]

        @self.Constraint(tset, iznodes)
        def heat_flux_x0_eqn(b, t, iz):
            return b.heat_flux_x0[t, iz] == b.heat_flux_x[t, ixfaces.first(), iz]

        @self.Constraint(tset, iznodes)
        def heat_flux_x1_eqn(b, t, iz):
            return b.heat_flux_x1[t, iz] == b.heat_flux_x[t, ixfaces.last(), iz]

        @self.Expression(tset, ixnodes, izfaces)
        def heat_flux_z(b, t, ix, iz):
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
                == b.xface_area[iz]
                * (b.heat_flux_x[t, ix, iz] - b.heat_flux_x[t, ix + 1, iz])
                + b.zface_area[ix]
                * (b.heat_flux_z[t, ix, iz] - b.heat_flux_z[t, ix, iz + 1])
                + b.joule_heating[t, ix, iz]
            )

        if dynamic:
            self.energy_balance_solid_eqn[tset.first(), :, :].deactivate()

    def calculate_scaling_factors(self):
        pass

    def initialize_build(self, outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        # We've got initial estimates of temperature and heat flux at both ends of the slab. The fluxes are actually the
        # important things to preserve, but we cannot fix both of them (BVPs with two sets of Neumann boundary
        # conditions are ill-defined).
        raise NotImplementedError(
            "Initialization for the conductive_slab unit model is not implemented because there is no obvious set "
            "of boundary conditions to fix during solid_oxide_cell initialization, and it is not meant to be "
            "initialized in isolation."
        )

    def model_check(self):
        pass

    def recursive_scaling(self):
        gsf = iscale.get_scaling_factor

        def ssf(c, s):
            iscale.set_scaling_factor(c, s, overwrite=False)

        sgsf = iscale.set_and_get_scaling_factor

        def cst(c, s):
            iscale.constraint_scaling_transform(c, s, overwrite=False)

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

                if self.heat_flux_x0[t, iz].is_reference():
                    sq0 = gsf(self.heat_flux_x0[t, iz].referent, default=1e-2)
                else:
                    sq0 = sgsf(self.heat_flux_x0[t, iz], 1e-2)
                cst(self.heat_flux_x0_eqn[t, iz], sq0)
                if not self.temperature_deviation_x0.is_reference():
                    ssf(self.temperature_deviation_x0, sq0 * sLx / sk)

                if self.heat_flux_x1[t, iz].is_reference():
                    sq1 = gsf(self.heat_flux_x1[t, iz].referent, default=1e-2)
                else:
                    sq1 = sgsf(self.heat_flux_x1[t, iz], 1e-2)
                cst(self.heat_flux_x1_eqn[t, iz], sq1)
                if not self.temperature_deviation_x1.is_reference():
                    ssf(self.temperature_deviation_x1, sq1 * sLx / sk)

                sqx = min(sq0, sq1)
                sqz = 10 * sqx  # Heuristic

                for ix in self.ixnodes:
                    sDT = sgsf(self.temperature_deviation_x[t, ix, iz], sqx * sLx / sk)

                    if self.config.has_holdup:
                        s_rho_U_solid = sgsf(
                            self.int_energy_density_solid[t, ix, iz],
                            1 / (self.heat_capacity.value * self.density.value * sDT),
                        )
                        cst(self.int_energy_density_solid_eqn[t, ix, iz], s_rho_U_solid)

                    cst(self.energy_balance_solid_eqn[t, ix, iz], sqx * sLy * sLz)
