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
Two-dimensional model of a porous, conductive medium, such as electrodes used in
a solid oxide fuel/electrolytic cell. Equations regarding the reactions that
occur at the triple phase boundary are located in the TriplePhaseBoundary class.
Material and heat transfer PDEs are  discretized in both x and z directions using
a finite volume method. Quantities are interpolated using a centered-difference
scheme. Presently only conventional diffusion is considered, because the Knudsen
number is much less than one for reactive species (H2, H2O, and O2). The energy
balance considers thermal condition from the solid phase, enthalpy flows from
gas species, and resistive heating from the current flowing through the medium.

If variables are not provided for boundary conditions at the ``x=0`` or ``x=1``
edges, new variables will be created for those terms. No-flux boundary conditions
are used along the ``z=0``and ``z=1`` edges.

Boundary variables:
    - ``temperature_deviation_x0[t, iz]``
    - ``temperature_deviation_x1[t, iz]``
    - ``heat_flux_x0[t, iz]``
    - ``heat_flux_x1[t, iz]``
    - ``conc_mol_comp_deviation_x0[t, iz]``
    - ``conc_mol_comp_deviation_x1[t, iz]``
    - ``material_flux_x0[t, iz]``
    - ``material_flux_x1[t, iz]``

Instances of ``Var`` that must be fixed:
    - ``length_x``: Thickness of porous conductive slab
    - ``porosity``: Porosity or void fraction of slab
    - ``tortuosity``: Tortuosity of slab
    - ``solid_heat_capacity``: Specific heat capacity of solid phase of slab on a mass-basis (i.e., in J/(kg * K))
    - ``solid_density``: Mass density of the solid component of the slab. Note that literature
      values that include void areas must be corrected by dividing by one minus porosity
    - ``solid_thermal_conductivity``: Thermal conductivity of solid phase of slab. Again, literature
      values may need to be divided by one minus porosity
    - ``resistivity_log_preexponential_factor``: Natural logarithm of resistivity preexponential factor in ohm * m.
      Again, literature values may need to be divided by one minus porosity
    - ``resistivity_thermal_exponent_dividend``: Parameter divided by temperature in resistivity equation, in K.
      Would be something like (reduced) activation energy, but it can be both negative and positive.
"""
__author__ = "John Eslick, Douglas Allan"

from pyomo.common.config import ConfigValue, Bool, ListOf
from pyomo.dae import DerivativeVar
import pyomo.environ as pyo


from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.util.constants import Constants
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
from idaes.models_extra.power_generation.unit_models.soc_submodels.common import (
    _set_if_unfixed,
    _gas_species_list,
    _element_list,
    _element_dict,
)
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


@declare_process_block_class("PorousConductiveSlab")
class PorousConductiveSlabData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "has_gas_holdup",
        ConfigValue(
            domain=Bool,
            default=False,
            description="If True, make holdup terms corresponding to gas phase",
        ),
    )
    CONFIG.declare(
        "control_volume_xfaces",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates of control volume faces "
            "in x direction. Coordinates must start with zero, be strictly "
            "increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "component_list",
        ConfigValue(
            domain=common._SubsetOf(_gas_species_list), description="List of components"
        ),
    )
    CONFIG.declare(
        "conc_mol_comp_ref",
        ConfigValue(
            default=None,
            description="Variable for the component concentration in bulk channel ",
        ),
    )
    CONFIG.declare(
        "dconc_mol_comp_refdt",
        ConfigValue(
            default=None,
            description="Variable for time derivative of the "
            "component concentration in the bulk channel",
        ),
    )
    common._submodel_boilerplate_config(CONFIG)
    common._thermal_boundary_conditions_config(CONFIG, thin=False)
    common._material_boundary_conditions_config(CONFIG, thin=False)

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        self.component_list = pyo.Set(
            initialize=self.config.component_list,
            ordered=True,
            doc="Set of all gas-phase components present in submodel",
        )
        # Set up node and face sets and get integer indices for them
        izfaces, iznodes = common._face_initializer(
            self, self.config.control_volume_zfaces, "z"
        )
        ixfaces, ixnodes = common._face_initializer(
            self, self.config.control_volume_xfaces, "x"
        )

        # Space saving aliases
        comps = self.component_list
        include_temp_x_thermo = self.config.include_temperature_x_thermo

        # Electrode thickness AKA length in the x direction is specific to the
        # electrode so local variable here is the only option
        self.length_x = pyo.Var(
            doc="Thickness of the slab (x-direction)",
            units=pyo.units.m,
        )

        if self.config.dynamic:
            # If we're dynamic, the user needs to either provide both the
            # reference concentration and its derivative, or provide neither.
            assert (
                self.config.conc_mol_comp_ref is None
                and self.config.dconc_mol_comp_refdt is None
            ) or (
                self.config.conc_mol_comp_ref is not None
                and self.config.dconc_mol_comp_refdt is not None
            )

        if self.config.has_gas_holdup and not self.config.has_holdup:
            raise ConfigurationError(
                "Creating gas holdup terms while not creating other holdup terms is not supported."
            )

        if self.config.conc_mol_comp_ref is None:
            self.conc_mol_comp_ref = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Concentration of components in the channel bulk",
                initialize=0.0,
                units=pyo.units.mol / pyo.units.m**3,
            )
            if self.config.dynamic and self.config.has_gas_holdup:
                self.dconc_mol_comp_refdt = DerivativeVar(
                    self.conc_mol_comp_ref,
                    wrt=tset,
                    doc="Derivative of concentration of components in the channel bulk",
                    initialize=0.0,
                    units=pyo.units.mol / (pyo.units.s * pyo.units.m**3),
                )
            else:
                self.dconc_mol_comp_refdt = pyo.Param(
                    tset,
                    iznodes,
                    comps,
                    doc="Derivative of concentration of components in the channel bulk",
                    initialize=0.0,
                    units=pyo.units.mol / (pyo.units.s * pyo.units.m**3),
                )
        else:
            self.conc_mol_comp_ref = pyo.Reference(self.config.conc_mol_comp_ref)
            if self.config.dconc_mol_comp_refdt.ctype == DerivativeVar:
                self.dconc_mol_comp_refdt = pyo.Reference(
                    self.config.dconc_mol_comp_refdt, ctype=pyo.Var
                )
            else:
                self.dconc_mol_comp_refdt = pyo.Reference(
                    self.config.dconc_mol_comp_refdt
                )
        common._submodel_boilerplate_create_if_none(self)
        common._create_thermal_boundary_conditions_if_none(self, thin=False)
        common._create_material_boundary_conditions_if_none(self, thin=False)

        self.porosity = pyo.Var(
            initialize=0.50, doc="Electrode porosity", units=pyo.units.dimensionless
        )
        self.tortuosity = pyo.Var(
            initialize=2.0, doc="Electrode tortuosity", units=pyo.units.dimensionless
        )

        self.temperature_deviation_x = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            doc="Deviation of temperature at node centers from temperature_z",
            units=pyo.units.K,
            bounds=(-1000, 1000),
        )
        self.conc_mol_comp_deviation_x = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            comps,
            doc="Deviation of component concentration at node centers "
            "from conc_mol_comp_ref",
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
            if self.config.has_gas_holdup:
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
        self.diff_eff_coeff = pyo.Var(
            tset,
            ixnodes,
            iznodes,
            comps,
            doc="Effective Fick's law diffusion coefficient at node centers",
            initialize=2e-5,
            bounds=(0, None),
            units=pyo.units.m**2 / pyo.units.s,
        )
        self.resistivity_log_preexponential_factor = pyo.Var(
            doc="Logarithm of resistivity preexponential factor " "in units of ohm*m",
            units=pyo.units.dimensionless,
        )
        self.resistivity_thermal_exponent_dividend = pyo.Var(
            doc="Parameter divided by temperature in exponential", units=pyo.units.K
        )

        # Parameters
        self.solid_heat_capacity = pyo.Var(
            initialize=200,
            doc="Heat capacity of solid part of porous slab (mass basis)",
            units=pyo.units.J / pyo.units.kg / pyo.units.K,
        )
        self.solid_density = pyo.Var(
            initialize=1000,
            doc="Density of solid part of porous slab",
            units=pyo.units.kg / pyo.units.m**3,
        )
        self.solid_thermal_conductivity = pyo.Var(
            initialize=80,
            doc="Thermal conductivity of solid part of porous slab",
            units=pyo.units.W / pyo.units.m / pyo.units.K,
        )

        # Add time derivative varaible if steady state use const 0.
        if dynamic and self.config.has_gas_holdup:
            self.dconc_mol_comp_deviation_xdt = DerivativeVar(
                self.conc_mol_comp_deviation_x,
                wrt=tset,
                initialize=0,
                doc="Component concentration time derivative in deviation " "variable",
            )
        else:
            self.dconc_mol_comp_deviation_xdt = pyo.Param(
                tset,
                ixnodes,
                iznodes,
                comps,
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3 / pyo.units.s,
            )
        # Add time derivative variable if steady state use const 0.
        if dynamic and self.config.has_gas_holdup:
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
        # Add time derivative variable if steady state use const 0.
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

        @self.Expression(tset, ixnodes, iznodes)
        def temperature(b, t, ix, iz):
            if include_temp_x_thermo:
                return b.temperature_z[t, iz] + b.temperature_deviation_x[t, ix, iz]
            else:
                return b.temperature_z[t, iz]

        @self.Expression(tset, iznodes, comps)
        def conc_mol_comp_x0(b, t, iz, j):
            return (
                b.conc_mol_comp_ref[t, iz, j] + b.conc_mol_comp_deviation_x0[t, iz, j]
            )

        @self.Expression(tset, ixnodes, iznodes, comps)
        def conc_mol_comp(b, t, ix, iz, j):
            return (
                b.conc_mol_comp_ref[t, iz, j]
                + b.conc_mol_comp_deviation_x[t, ix, iz, j]
            )

        @self.Expression(tset, iznodes, comps)
        def conc_mol_comp_x1(b, t, iz, j):
            return (
                b.conc_mol_comp_ref[t, iz, j] + b.conc_mol_comp_deviation_x1[t, iz, j]
            )

        @self.Expression(tset, ixnodes, iznodes, comps)
        def dconc_mol_compdt(b, t, ix, iz, j):
            return (
                b.dconc_mol_comp_refdt[t, iz, j]
                + b.dconc_mol_comp_deviation_xdt[t, ix, iz, j]
            )

        @self.Expression(tset, ixnodes, iznodes)
        def vol_mol(b, t, ix, iz):
            return (
                Constants.gas_constant
                * b.temperature[t, ix, iz]
                / b.pressure[t, ix, iz]
            )

        @self.Constraint(tset, ixnodes, iznodes, comps)
        def conc_mol_comp_eqn(b, t, ix, iz, i):
            return (
                b.conc_mol_comp[t, ix, iz, i]
                * b.temperature[t, ix, iz]
                * Constants.gas_constant
                == b.pressure[t, ix, iz] * b.mole_frac_comp[t, ix, iz, i]
            )

        @self.Constraint(tset, ixnodes, iznodes)
        def enth_mol_eqn(b, t, ix, iz):
            return b.enth_mol[t, ix, iz] == sum(
                common._comp_enthalpy_expr(b.temperature[t, ix, iz], i)
                * b.mole_frac_comp[t, ix, iz, i]
                for i in comps
            )

        if self.config.has_holdup:
            if self.config.has_gas_holdup:
                # For the vapor phase
                @self.Constraint(tset, ixnodes, iznodes)
                def int_energy_mol_eqn(b, t, ix, iz):
                    return b.int_energy_mol[t, ix, iz] == sum(
                        common._comp_int_energy_expr(b.temperature[t, ix, iz], i)
                        * b.mole_frac_comp[t, ix, iz, i]
                        for i in comps
                    )

                @self.Constraint(tset, ixnodes, iznodes)
                def int_energy_density_eqn(b, t, ix, iz):
                    return (
                        b.int_energy_density[t, ix, iz]
                        == b.int_energy_mol[t, ix, iz] / b.vol_mol[t, ix, iz]
                    )

            @self.Constraint(tset, ixnodes, iznodes)
            def int_energy_density_solid_eqn(b, t, ix, iz):
                return b.int_energy_density_solid[
                    t, ix, iz
                ] == b.solid_heat_capacity * b.solid_density * (
                    b.temperature[t, ix, iz] - 1000 * pyo.units.K
                )

        @self.Constraint(tset, ixnodes, iznodes)
        def mole_frac_comp_eqn(b, t, ix, iz):
            return 1 == sum(b.mole_frac_comp[t, ix, iz, i] for i in comps)

        @self.Constraint(tset, ixnodes, iznodes, comps)
        def diff_eff_coeff_eqn(b, t, ix, iz, i):
            T = b.temperature[t, ix, iz]
            P = b.pressure[t, ix, iz]
            x = b.mole_frac_comp
            bfun = common._binary_diffusion_coefficient_expr
            return b.diff_eff_coeff[t, ix, iz, i] == b.porosity / b.tortuosity * (
                1.0 - x[t, ix, iz, i]
            ) / sum(x[t, ix, iz, j] / bfun(T, P, i, j) for j in comps if i != j)

        @self.Expression(tset, ixfaces, iznodes, comps)
        def dcdx(b, t, ix, iz, i):
            return common._interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=b.xnodes,
                faces=b.xfaces,
                phi_func=lambda ixf: b.conc_mol_comp_deviation_x[t, ixf, iz, i]
                / b.length_x,
                phi_bound_0=(
                    b.conc_mol_comp_deviation_x[t, ixnodes.first(), iz, i]
                    - b.conc_mol_comp_deviation_x0[t, iz, i]
                )
                / (b.xnodes.first() - b.xfaces.first())
                / b.length_x,
                phi_bound_1=(
                    b.conc_mol_comp_deviation_x1[t, iz, i]
                    - b.conc_mol_comp_deviation_x[t, ixnodes.last(), iz, i]
                )
                / (b.xfaces.last() - b.xnodes.last())
                / b.length_x,
                derivative=True,
            )

        @self.Expression(tset, ixnodes, izfaces, comps)
        def dcdz(b, t, ix, iz, i):
            return common._interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=b.znodes,
                faces=b.zfaces,
                phi_func=lambda izf: b.conc_mol_comp[t, ix, izf, i] / b.length_z[None],
                phi_bound_0=0,  # solid wall no flux
                phi_bound_1=0,  # solid wall no flux
                derivative=True,
            )

        @self.Expression(tset, ixfaces, iznodes)
        def dTdx(b, t, ix, iz):
            return common._interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=b.xnodes,
                faces=b.xfaces,
                phi_func=lambda ixf: b.temperature_deviation_x[t, ixf, iz] / b.length_x,
                phi_bound_0=(
                    b.temperature_deviation_x[t, ixnodes.first(), iz]
                    - b.temperature_deviation_x0[t, iz]
                )
                / (b.xnodes.first() - b.xfaces.first())
                / b.length_x,
                phi_bound_1=(
                    b.temperature_deviation_x1[t, iz]
                    - b.temperature_deviation_x[t, ixnodes.last(), iz]
                )
                / (b.xfaces.last() - b.xnodes.last())
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

        @self.Expression(tset, ixfaces, iznodes, comps)
        def diff_eff_coeff_xfaces(b, t, ix, iz, i):
            return common._interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=b.xnodes,
                faces=b.xfaces,
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
            return common._interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=b.znodes,
                faces=b.zfaces,
                phi_func=lambda izf: b.diff_eff_coeff[t, ix, izf, i],
                phi_bound_0=0,  # solid wall no flux
                phi_bound_1=0,  # solid wall no flux
                derivative=False,
            )

        @self.Expression(tset, ixfaces, iznodes)
        def temperature_xfaces(b, t, ix, iz):
            return common._interpolate_2D(
                ic=ix,
                ifaces=ixfaces,
                nodes=b.xnodes,
                faces=b.xfaces,
                phi_func=lambda ixf: b.temperature[t, ixf, iz],
                phi_bound_0=b.temperature_x0[t, iz],
                phi_bound_1=b.temperature_x1[t, iz],
                derivative=False,
            )

        @self.Expression(tset, ixnodes, izfaces)
        def temperature_zfaces(b, t, ix, iz):
            return common._interpolate_2D(
                ic=iz,
                ifaces=izfaces,
                nodes=b.znodes,
                faces=b.zfaces,
                phi_func=lambda izf: b.temperature[t, ix, izf],
                phi_bound_0=b.temperature[t, ix, iznodes.first()],
                phi_bound_1=b.temperature[t, ix, iznodes.last()],
                derivative=False,
            )

        @self.Expression(tset, ixfaces, iznodes, comps)
        def material_flux_x(b, t, ix, iz, i):
            return -b.dcdx[t, ix, iz, i] * b.diff_eff_coeff_xfaces[t, ix, iz, i]

        @self.Constraint(tset, iznodes, comps)
        def material_flux_x0_eqn(b, t, iz, i):
            return (
                b.material_flux_x[t, ixfaces.first(), iz, i]
                == b.material_flux_x0[t, iz, i]
            )

        @self.Constraint(tset, iznodes, comps)
        def material_flux_x1_eqn(b, t, iz, i):
            return (
                b.material_flux_x[t, ixfaces.last(), iz, i]
                == b.material_flux_x1[t, iz, i]
            )

        @self.Expression(tset, ixnodes, izfaces, comps)
        def material_flux_z(b, t, ix, iz, i):
            return -b.dcdz[t, ix, iz, i] * b.diff_eff_coeff_zfaces[t, ix, iz, i]

        @self.Expression(tset, ixfaces, iznodes)
        def heat_flux_x(b, t, ix, iz):
            return -(1 - b.porosity) * b.solid_thermal_conductivity * b.dTdx[t, ix, iz]

        @self.Expression(tset, ixnodes, izfaces)
        def heat_flux_z(b, t, ix, iz):
            return -(1 - b.porosity) * b.solid_thermal_conductivity * b.dTdz[t, ix, iz]

        @self.Constraint(tset, iznodes)
        def heat_flux_x0_eqn(b, t, iz):
            return b.heat_flux_x0[t, iz] == b.heat_flux_x[t, ixfaces.first(), iz]

        @self.Constraint(tset, iznodes)
        def heat_flux_x1_eqn(b, t, iz):
            return b.heat_flux_x1[t, iz] == b.heat_flux_x[t, ixfaces.last(), iz]

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
            return b.porosity * b.node_volume[ix, iz] * b.dconc_mol_compdt[
                t, ix, iz, i
            ] == b.xface_area[iz] * (
                b.material_flux_x[t, ix, iz, i] - b.material_flux_x[t, ix + 1, iz, i]
            ) + b.zface_area[
                ix
            ] * (
                b.material_flux_z[t, ix, iz, i] - b.material_flux_z[t, ix, iz + 1, i]
            )

        if dynamic and self.config.has_gas_holdup:
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
                == b.xface_area[iz]
                * (b.heat_flux_x[t, ix, iz] - b.heat_flux_x[t, ix + 1, iz])
                + b.zface_area[ix]
                * (b.heat_flux_z[t, ix, iz] - b.heat_flux_z[t, ix, iz + 1])
                + b.joule_heating[t, ix, iz]
                # For mass flux heat transfer include exchange with channel
                # probably make little difference, but want to ensure the energy
                # balance closes
                + b.xface_area[iz]
                * sum(
                    b.material_flux_x[t, ix, iz, i]
                    * common._comp_enthalpy_expr(b.temperature_xfaces[t, ix, iz], i)
                    for i in comps
                )
                - b.xface_area[iz]
                * sum(
                    b.material_flux_x[t, ix + 1, iz, i]
                    * common._comp_enthalpy_expr(b.temperature_xfaces[t, ix + 1, iz], i)
                    for i in comps
                )
                + b.zface_area[ix]
                * sum(
                    b.material_flux_z[t, ix, iz, i]
                    * common._comp_enthalpy_expr(b.temperature_zfaces[t, ix, iz], i)
                    for i in comps
                )
                - b.zface_area[ix]
                * sum(
                    b.material_flux_z[t, ix, iz + 1, i]
                    * common._comp_enthalpy_expr(b.temperature_zfaces[t, ix, iz + 1], i)
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

        comps = self.component_list
        # Use conc_x0 and the temperature guess to start filling in initial
        # guess values.
        for t in self.flowsheet().time:
            for iz in self.iznodes:
                for i in comps:
                    _set_if_unfixed(
                        self.conc_mol_comp_deviation_x[t, self.ixnodes.first(), iz, i],
                        self.conc_mol_comp_deviation_x0[t, iz, i],
                    )
                for ix in self.ixnodes:
                    if temperature_guess is not None:
                        _set_if_unfixed(self.temperature_z[t, iz], temperature_guess)
                        _set_if_unfixed(self.temperature_deviation_x[t, ix, iz], 0)
                        _set_if_unfixed(self.temperature_deviation_x0[t, iz], 0)
                        _set_if_unfixed(self.temperature_deviation_x1[t, iz], 0)
                    if pressure_guess is not None:
                        _set_if_unfixed(self.pressure[t, ix, iz], pressure_guess)
                    for i in comps:
                        _set_if_unfixed(self.conc_mol_comp_deviation_x[t, ix, iz, i], 0)
                    mol_dens = pyo.value(
                        sum(self.conc_mol_comp[t, ix, iz, i] for i in comps)
                    )
                    _set_if_unfixed(
                        self.pressure[t, ix, iz],
                        Constants.gas_constant * self.temperature[t, ix, iz] * mol_dens,
                    )
                    for i in comps:
                        _set_if_unfixed(
                            self.mole_frac_comp[t, ix, iz, i],
                            self.conc_mol_comp[t, ix, iz, i] / mol_dens,
                        )
                    _set_if_unfixed(
                        self.enth_mol[t, ix, iz],
                        sum(
                            common._comp_enthalpy_expr(self.temperature[t, ix, iz], i)
                            * self.mole_frac_comp[t, ix, iz, i]
                            for i in comps
                        ),
                    )
                    if self.config.has_gas_holdup:
                        _set_if_unfixed(
                            self.int_energy_mol[t, ix, iz],
                            sum(
                                common._comp_int_energy_expr(
                                    self.temperature[t, ix, iz], i
                                )
                                * self.mole_frac_comp[t, ix, iz, i]
                                for i in comps
                            ),
                        )
                        _set_if_unfixed(
                            self.int_energy_density[t, ix, iz],
                            self.int_energy_mol[t, ix, iz] / self.vol_mol[t, ix, iz],
                        )
                for i in comps:
                    _set_if_unfixed(
                        self.conc_mol_comp_deviation_x1[t, iz, i],
                        self.conc_mol_comp_deviation_x[t, self.ixnodes.last(), iz, i],
                    )

        solver_obj = get_solver(solver, optarg)
        common._init_solve_block(self, solver_obj, solve_log)

    def calculate_scaling_factors(self):
        pass

    def model_check(self, steady_state=True):
        comp_set = set(self.component_list)
        elements_present = set()

        for element in _element_list:
            include_element = False
            for species in _gas_species_list:
                # Floating point equality take warning!
                if species in comp_set and _element_dict[element][species] != 0:
                    include_element = True
            if include_element:
                elements_present.add(element)

        if not steady_state:
            # Mass and energy conservation equations steady state only at present
            return
        for t in self.flowsheet().config.time:
            for element in _element_list:
                if element not in elements_present:
                    continue
                sum_in = 0
                sum_out = 0
                for iz in self.iznodes:
                    sum_in += sum(
                        _element_dict[element][j]
                        * self.material_flux_x0[t, iz, j]
                        * self.xface_area[iz]
                        for j in self.component_list
                    )
                    sum_out += sum(
                        _element_dict[element][j]
                        * self.material_flux_x1[t, iz, j]
                        * self.xface_area[iz]
                        for j in self.component_list
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
                    self.heat_flux_x0[t, iz]
                    + sum(
                        common._comp_enthalpy_expr(self.temperature_x0[t, iz], j)
                        * self.material_flux_x0[t, iz, j]
                        for j in self.component_list
                    )
                )
                enth_out += self.xface_area[iz] * (
                    self.heat_flux_x1[t, iz]
                    + sum(
                        common._comp_enthalpy_expr(self.temperature_x1[t, iz], j)
                        * self.material_flux_x1[t, iz, j]
                        for j in self.component_list
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

        def ssf(c, s):
            iscale.set_scaling_factor(c, s, overwrite=False)

        sgsf = iscale.set_and_get_scaling_factor

        def cst(c, s):
            iscale.constraint_scaling_transform(c, s, overwrite=False)

        sR = 1e-1  # Scaling factor for R
        sD = 1e5  # Heuristic scaling factor for diffusion coefficient
        sy_def = 10  # Mole frac comp scaling
        sh = 1e-2  # Heat xfer coeff
        sH = 1e-4  # Enthalpy/int energy
        sk = 0.1  # Fudge factor to scale temperature_deviation_x
        sLx = sgsf(self.length_x, len(self.ixnodes) / self.length_x.value)
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

                smaterial_flux_x = {}
                for j in self.component_list:
                    # ssf(self.conc_mol_comp_ref[t,iz,j],sy_def*1e-4/(sR*sT))

                    if self.material_flux_x0[t, iz, j].is_reference():
                        smaterial_flux_x0 = gsf(
                            self.material_flux_x0[t, iz, j].referent, default=1e-2
                        )
                    else:
                        smaterial_flux_x0 = sgsf(self.material_flux_x0[t, iz, j], 1e-2)
                    cst(self.material_flux_x0_eqn[t, iz, j], smaterial_flux_x0)
                    if not self.conc_mol_comp_deviation_x0[t, iz, j].is_reference():
                        ssf(
                            self.conc_mol_comp_deviation_x0[t, iz, j],
                            smaterial_flux_x0 * sLx / sD,
                        )

                    if self.material_flux_x1[t, iz, j].is_reference():
                        smaterial_flux_x1 = gsf(
                            self.material_flux_x1[t, iz, j].referent, default=1e-2
                        )
                    else:
                        smaterial_flux_x1 = sgsf(self.material_flux_x1[t, iz, j], 1e-2)
                    cst(self.material_flux_x1_eqn[t, iz, j], smaterial_flux_x1)

                    if not self.conc_mol_comp_deviation_x1[t, iz, j].is_reference():
                        ssf(
                            self.conc_mol_comp_deviation_x1[t, iz, j],
                            smaterial_flux_x1 * sLx / sD,
                        )

                    smaterial_flux_x[j] = min(smaterial_flux_x0, smaterial_flux_x1)

                for ix in self.ixnodes:
                    sP = sgsf(self.pressure[t, ix, iz], 1e-4)
                    sV = sR * sT / sP

                    sDT = sgsf(self.temperature_deviation_x[t, ix, iz], sqx * sLx / sk)
                    sH = sgsf(self.enth_mol[t, ix, iz], sH)
                    cst(self.enth_mol_eqn[t, ix, iz], sH)

                    if self.config.has_holdup:
                        if self.config.has_gas_holdup:
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

                    cst(self.mole_frac_comp_eqn[t, ix, iz], 1)
                    cst(self.energy_balance_solid_eqn[t, ix, iz], sqx * sLy * sLz)

                    for j in self.component_list:
                        sy = sgsf(self.mole_frac_comp[t, ix, iz, j], sy_def)
                        cst(self.conc_mol_comp_eqn[t, ix, iz, j], sy * sP)
                        ssf(
                            self.conc_mol_comp_deviation_x[t, ix, iz, j],
                            smaterial_flux_x[j] * sLx / sD,
                        )

                        cst(
                            self.material_balance_eqn[t, ix, iz, j],
                            smaterial_flux_x[j] * sLy * sLz,
                        )
                        # FIXME come back later to clean up
                        ssf(self.diff_eff_coeff[t, ix, iz, j], 1e5)
                        cst(self.diff_eff_coeff_eqn[t, ix, iz, j], 1e5)
