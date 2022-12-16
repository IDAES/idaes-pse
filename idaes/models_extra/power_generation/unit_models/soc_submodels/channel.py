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
Model for channel in SolidOxideCell. If opposite_flow is false, material flows from
``z=0`` to ``z=1``, otherwise it flows from ``z=1`` to ``z=0``. Pure convection occurs in the z
direction, while mass and heat transfer coefficients govern fluxes out of the ``x=0``
and ``x=1`` sides of the model. The PDEs are discretized in the z direction  by a finite
volume method with upwind differences used to calculate convection fluxes. It is
assumed that the gas is ideal.

Like most submodels, the SocChannel will create terms for temperature and heat flux
at both x=0 and x=1 edges if they are not provided in the config. However,
concentration variables and xflux terms are created only on the side facing the
electrode (``x=1`` if below_electrode is true, ``x=0`` if false). Expressions are created
for the other side, to avoid additional logic with material balances.

Manipulated variables (typically fixed through cell-level ports):
    - ``flow_mol_inlet[t]``: Molar flow into channel
    - ``mole_frac_comp_inlet[t, j]``: Molar composition of flow into channel
    - ``temperature_inlet[t]``: Channel inlet flow rate
    - ``pressure_inlet[t]``: Channel inlet pressure

Boundary variables:
    - ``temperature_deviation_x0[t, iz]``
    - ``temperature_deviation_x1[t, iz]``
    - ``heat_flux_x0[t, iz]``
    - ``heat_flux_x1[t, iz]``

    * If ``below_electrode = False``:
        - ``material_flux_x0[t, iz]``
        - ``conc_mol_comp_deviation_x0[t, iz]``

    * If ``below_electrode = True``:
        - ``material_flux_x1[t, iz]``
        - ``conc_mol_comp_deviation_x1[t, iz]``

Instances of ``Var`` that must be fixed:
    - ``length_x``: Distance between electrode and interconnect
    - ``heat_transfer_coefficient[t, iz]``: Heat transfer coefficient between channel bulk and surfaces.
      Correlations to calculate this are planned, but they will add a large number of equations.
"""

__author__ = "John Eslick, Douglas Allan"

from pyomo.common.config import ConfigValue, In, ListOf, Bool
from pyomo.dae import DerivativeVar
import pyomo.environ as pyo

from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.core.util.constants import Constants
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
from idaes.models_extra.power_generation.unit_models.soc_submodels.common import (
    _set_if_unfixed,
    _gas_species_list,
    _element_list,
    _element_dict,
)
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

from idaes.core.base.var_like_expression import VarLikeExpression
import idaes.logger as idaeslog


@declare_process_block_class("SocChannel")
class SocChannelData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "component_list",
        ConfigValue(
            domain=common._SubsetOf(_gas_species_list), description="List of components"
        ),
    )

    CONFIG.declare(
        "control_volume_zfaces",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates of control volume faces "
            "in z direction. Coordinates must start with zero, be strictly "
            "increasing, and end with one",
        ),
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
    common._thermal_boundary_conditions_config(CONFIG, thin=False)
    CONFIG.declare(
        "opposite_flow",
        ConfigValue(
            domain=Bool,
            default=False,
            description="If True, gas flows from z=1 to z=0. Otherwise "
            "gas flows from z=0 to z=1.",
        ),
    )
    CONFIG.declare(
        "below_electrode",
        ConfigValue(
            domain=Bool,
            description="If True, the top of the channel is the electrode "
            "boundary and material flux terms are created there. If False "
            "the bottom is the electrode boundary and material flux terms "
            "are created there.",
        ),
    )

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        self.component_list = pyo.Set(
            initialize=self.config.component_list,
            ordered=True,
            doc="Set of all components present in submodel",
        )
        # Set up node and face sets and get integer indices for them
        izfaces, iznodes = common._face_initializer(
            self, self.config.control_volume_zfaces, "z"
        )
        # Space saving aliases
        comps = self.component_list
        include_temp_x_thermo = self.config.include_temperature_x_thermo

        common._create_if_none(self, "length_z", idx_set=None, units=pyo.units.m)
        common._create_if_none(self, "length_y", idx_set=None, units=pyo.units.m)

        common._create_thermal_boundary_conditions_if_none(self, thin=False)

        if not self.config.below_electrode:
            self.conc_mol_comp_deviation_x0 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Deviation of concentration at electrode surface from that "
                "in the bulk channel",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3,
            )
            self.material_flux_x0 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Material flux from electrode surface to channel "
                "(positive is in)",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**2 / pyo.units.s,
            )

            self.conc_mol_comp_deviation_x1 = pyo.Param(
                tset,
                iznodes,
                comps,
                doc="Dummy parameter for zero concentration deviation at interconnect wall.",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3,
            )
            self.material_flux_x1 = pyo.Param(
                tset,
                iznodes,
                comps,
                doc="Dummy parameter for zero material flux from channel "
                "to interconnect.",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**2 / pyo.units.s,
            )

        else:
            self.conc_mol_comp_deviation_x1 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Deviation of concentration at electrode surface from that "
                "in the bulk channel",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3,
            )
            self.material_flux_x1 = pyo.Var(
                tset,
                iznodes,
                comps,
                doc="Material flux from channel to electrode surface "
                "(positive is out)",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**2 / pyo.units.s,
            )
            self.conc_mol_comp_deviation_x0 = pyo.Param(
                tset,
                iznodes,
                comps,
                doc="Dummy parameter for zero concentration deviation at interconnect wall.",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3,
            )
            self.material_flux_x0 = pyo.Param(
                tset,
                iznodes,
                comps,
                doc="Dummy parameter for zero material flux from channel "
                "to interconnect.",
                initialize=0,
                units=pyo.units.mol / pyo.units.m**2 / pyo.units.s,
            )

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
        self.conc_mol_comp = pyo.Var(
            tset,
            iznodes,
            comps,
            doc="Component concentration at node centers",
            units=pyo.units.mol / pyo.units.m**3,
            bounds=(0, None),
        )
        self.temperature_deviation_x = pyo.Var(
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
        self.diff_eff_coeff = pyo.Var(
            tset,
            iznodes,
            comps,
            doc="Effective Fick's law diffusion coefficient at node centers",
            initialize=2e-5,
            bounds=(0, None),
            units=pyo.units.m**2 / pyo.units.s,
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
            tset, doc="Inlet temperature", bounds=(300, None), units=pyo.units.K
        )
        self.temperature_outlet = pyo.Var(
            tset, doc="Outlet temperature", bounds=(300, None), units=pyo.units.K
        )
        self.mole_frac_comp_inlet = pyo.Var(
            tset,
            comps,
            doc="Inlet component mole fractions",
            bounds=(0, 1),
            units=pyo.units.dimensionless,
        )

        # Add time derivative variable if steady state use const 0.
        if dynamic:
            self.dconc_mol_compdt = DerivativeVar(
                self.conc_mol_comp,
                wrt=tset,
                initialize=0,
                doc="Component concentration time derivative",
            )
        else:
            self.dconc_mol_compdt = pyo.Param(
                tset,
                iznodes,
                comps,
                initialize=0,
                units=pyo.units.mol / pyo.units.m**3 / pyo.units.s,
            )

        # Add time derivative variable if steady state use const 0.
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
                return b.temperature_z[t, iz] + b.temperature_deviation_x[t, iz]
            else:
                return b.temperature_z[t, iz]

        @self.Expression(tset, iznodes)
        def vol_mol(b, t, iz):
            return Constants.gas_constant * b.temperature[t, iz] / b.pressure[t, iz]

        @self.Expression(tset)
        def vol_mol_inlet(b, t):
            return Constants.gas_constant * b.temperature_inlet[t] / b.pressure_inlet[t]

        # TODO maybe replace with variable-constraint pair?
        @self.Expression(tset)
        def enth_mol_inlet(b, t):
            return sum(
                common._comp_enthalpy_expr(b.temperature_inlet[t], i)
                * b.mole_frac_comp_inlet[t, i]
                for i in comps
            )

        @self.Constraint(tset, iznodes, comps)
        def diff_eff_coeff_eqn(b, t, iz, i):
            T = b.temperature[t, iz]
            P = b.pressure[t, iz]
            x = b.mole_frac_comp
            bfun = common._binary_diffusion_coefficient_expr
            return b.diff_eff_coeff[t, iz, i] == (1.0 - x[t, iz, i]) / sum(
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
            def material_flux_x0_eqn(b, t, iz, i):
                return (
                    -b.material_flux_x0[t, iz, i] / b.mass_transfer_coeff[t, iz, i]
                    == b.conc_mol_comp_deviation_x0[t, iz, i]
                )

        else:

            @self.Constraint(tset, iznodes, comps)
            def material_flux_x1_eqn(b, t, iz, i):
                return (
                    -b.material_flux_x1[t, iz, i] / b.mass_transfer_coeff[t, iz, i]
                    == b.conc_mol_comp_deviation_x1[t, iz, i]
                )

        @self.Constraint(tset, iznodes)
        def flow_mol_eqn(b, t, iz):
            # either way the flow goes, want the flow rate to be positive, but
            # in the opposite flow cases want flux and velocity to be negative
            if self.config.opposite_flow:
                return (
                    b.flow_mol[t, iz]
                    == -b.flow_area * b.velocity[t, iz] / b.vol_mol[t, iz]
                )
            return (
                b.flow_mol[t, iz] == b.flow_area * b.velocity[t, iz] / b.vol_mol[t, iz]
            )

        @self.Constraint(tset, iznodes)
        def constant_pressure_eqn(b, t, iz):
            return b.pressure[t, iz] == b.pressure_inlet[t]

        @self.Constraint(tset, iznodes, comps)
        def conc_mol_comp_eqn(b, t, iz, i):
            return (
                b.conc_mol_comp[t, iz, i]
                * b.temperature[t, iz]
                * Constants.gas_constant
                == b.pressure[t, iz] * b.mole_frac_comp[t, iz, i]
            )

        @self.Constraint(tset, iznodes)
        def enth_mol_eqn(b, t, iz):
            return b.enth_mol[t, iz] == sum(
                common._comp_enthalpy_expr(b.temperature[t, iz], i)
                * b.mole_frac_comp[t, iz, i]
                for i in comps
            )

        if self.config.has_holdup:

            @self.Constraint(tset, iznodes)
            def int_energy_mol_eqn(b, t, iz):
                return b.int_energy_mol[t, iz] == sum(
                    common._comp_int_energy_expr(b.temperature[t, iz], i)
                    * b.mole_frac_comp[t, iz, i]
                    for i in comps
                )

            @self.Constraint(tset, iznodes)
            def int_energy_density_eqn(b, t, iz):
                return (
                    b.int_energy_density[t, iz]
                    == b.int_energy_mol[t, iz] / b.vol_mol[t, iz]
                )

        @self.Constraint(tset, iznodes)
        def mole_frac_comp_eqn(b, t, iz):
            return 1 == sum(b.mole_frac_comp[t, iz, i] for i in comps)

        @self.Expression(tset, comps)
        def flow_mol_comp_inlet(b, t, i):
            return b.flow_mol_inlet[t] * b.mole_frac_comp_inlet[t, i]

        @self.Expression(tset, comps)
        def material_flux_z_inlet(b, t, i):
            # either way the flow goes, want the flow rate to be positive, but
            # in the opposite flow cases want flux and velocity to be negative
            if self.config.opposite_flow:
                return -b.flow_mol_inlet[t] / b.flow_area * b.mole_frac_comp_inlet[t, i]
            return b.flow_mol_inlet[t] / b.flow_area * b.mole_frac_comp_inlet[t, i]

        @self.Expression(tset)
        def enth_flux_z_inlet(b, t):
            # either way the flow goes, want the flow rate to be positive, but
            # in the opposite flow cases want flux and velocity to be negative
            if self.config.opposite_flow:
                return -b.flow_mol_inlet[t] / b.flow_area * b.enth_mol_inlet[t]
            return b.flow_mol_inlet[t] / b.flow_area * b.enth_mol_inlet[t]

        @self.Expression(tset, izfaces, comps)
        def material_flux_z(b, t, iz, i):
            return common._interpolate_channel(
                iz=iz,
                ifaces=izfaces,
                nodes=b.znodes,
                faces=b.zfaces,
                phi_func=lambda iface: b.velocity[t, iface]
                * b.conc_mol_comp[t, iface, i],
                phi_inlet=b.material_flux_z_inlet[t, i],
                opposite_flow=self.config.opposite_flow,
            )

        @self.Expression(tset, izfaces)
        def material_flux_z_enth(b, t, iz):
            return common._interpolate_channel(
                iz=iz,
                ifaces=izfaces,
                nodes=b.znodes,
                faces=b.zfaces,
                phi_func=lambda iface: b.velocity[t, iface]
                / b.vol_mol[t, iface]
                * b.enth_mol[t, iface],
                phi_inlet=b.enth_flux_z_inlet[t],
                opposite_flow=self.config.opposite_flow,
            )

        @self.Expression(tset, izfaces)
        def pressure_face(b, t, iz):
            # Although I'm currently assuming no pressure drop in the channel
            # and don't have a momentum balance, this will let me estimate the
            # outlet pressure in a way that will let me add in a momentum balance
            # later
            return common._interpolate_channel(
                iz=iz,
                ifaces=izfaces,
                nodes=b.znodes,
                faces=b.zfaces,
                phi_func=lambda iface: b.pressure[t, iface],
                phi_inlet=b.pressure_inlet[t],
                opposite_flow=self.config.opposite_flow,
            )

        @self.Constraint(tset, iznodes, comps)
        def material_balance_eqn(b, t, iz, i):
            return b.dconc_mol_compdt[t, iz, i] * b.node_volume[iz] == b.flow_area * (
                b.material_flux_z[t, iz, i] - b.material_flux_z[t, iz + 1, i]
            ) + b.xface_area[iz] * (
                b.material_flux_x0[t, iz, i] - b.material_flux_x1[t, iz, i]
            )

        if dynamic:
            self.material_balance_eqn[tset.first(), :, :].deactivate()

        @self.Constraint(tset, iznodes)
        def energy_balance_eqn(b, t, iz):
            return (
                b.dcedt[t, iz] * b.node_volume[iz]
                == b.flow_area
                * (b.material_flux_z_enth[t, iz] - b.material_flux_z_enth[t, iz + 1])
                + b.xface_area[iz]
                * sum(
                    (b.material_flux_x0[t, iz, i] - b.material_flux_x1[t, iz, i])
                    * common._comp_enthalpy_expr(b.temperature_x1[t, iz], i)
                    for i in comps
                )
                + b.heat_flux_x0[t, iz] * b.xface_area[iz]
                - b.heat_flux_x1[t, iz] * b.xface_area[iz]
            )

        if dynamic:
            self.energy_balance_eqn[tset.first(), :].deactivate()

        @self.Constraint(tset, iznodes)
        def temperature_x0_eqn(b, t, iz):
            return self.heat_flux_x0[t, iz] == self.heat_transfer_coefficient[t, iz] * (
                self.temperature_deviation_x0[t, iz]
                - self.temperature_deviation_x[t, iz]
            )

        @self.Constraint(tset, iznodes)
        def temperature_x1_eqn(b, t, iz):
            return self.heat_flux_x1[t, iz] == self.heat_transfer_coefficient[t, iz] * (
                self.temperature_deviation_x[t, iz]
                - self.temperature_deviation_x1[t, iz]
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
                return -b.material_flux_z[t, izfout, i] * b.flow_area
            else:
                return b.material_flux_z[t, izfout, i] * b.flow_area

        def rule_flow_mol_outlet(b, t):
            return sum(b.flow_mol_comp_outlet[t, i] for i in comps)

        def rule_pressure_outlet(b, t):
            return b.pressure_face[t, izfout]

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
                return (
                    -b.material_flux_z_enth[t, izfout]
                    * b.flow_area
                    / b.flow_mol_outlet[t]
                )
            else:
                return (
                    b.material_flux_z_enth[t, izfout]
                    * b.flow_area
                    / b.flow_mol_outlet[t]
                )

        # know enthalpy need a constraint to back calculate temperature
        @self.Constraint(tset)
        def temperature_outlet_eqn(b, t):
            return b.enth_mol_outlet[t] == sum(
                common._comp_enthalpy_expr(b.temperature_outlet[t], i)
                * b.mole_frac_comp_outlet[t, i]
                for i in comps
            )

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
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
                _set_if_unfixed(self.temperature_deviation_x0[t, iz], 0)
                _set_if_unfixed(self.temperature_deviation_x[t, iz], 0)
                _set_if_unfixed(self.temperature_deviation_x1[t, iz], 0)
                _set_if_unfixed(self.pressure[t, iz], self.pressure_inlet[t])
                _set_if_unfixed(self.flow_mol[t, iz], self.flow_mol_inlet[t])
                for i in self.component_list:
                    _set_if_unfixed(
                        self.mole_frac_comp[t, iz, i], self.mole_frac_comp_inlet[t, i]
                    )
                    _set_if_unfixed(
                        self.conc_mol_comp[t, iz, i],
                        self.pressure[t, iz]
                        * self.mole_frac_comp[t, iz, i]
                        / self.temperature[t, iz]
                        / Constants.gas_constant,
                    )
                _set_if_unfixed(
                    self.velocity[t, iz],
                    self.flow_mol[t, iz] / self.flow_area * self.vol_mol[t, iz],
                )
                _set_if_unfixed(
                    self.enth_mol[t, iz],
                    sum(
                        common._comp_enthalpy_expr(self.temperature[t, iz], i)
                        * self.mole_frac_comp[t, iz, i]
                        for i in self.component_list
                    ),
                )
                if self.config.has_holdup:
                    _set_if_unfixed(
                        self.int_energy_mol[t, iz],
                        sum(
                            common._comp_int_energy_expr(self.temperature[t, iz], i)
                            * self.mole_frac_comp[t, iz, i]
                            for i in self.component_list
                        ),
                    )
                    _set_if_unfixed(
                        self.int_energy_density[t, iz],
                        self.int_energy_mol[t, iz] / self.vol_mol[t, iz],
                    )
        solver_obj = get_solver(solver, optarg)
        common._init_solve_block(self, solver_obj, solve_log)

    def calculate_scaling_factors(self):
        pass

    def model_check(self, steady_state=True):
        if not steady_state:
            # Mass and energy conservation equations steady state only at present
            return

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

        for t in self.flowsheet().config.time:
            for element in _element_list:
                if element not in elements_present:
                    continue
                sum_in = sum(
                    _element_dict[element][j] * self.flow_mol_comp_inlet[t, j]
                    for j in self.component_list
                )
                sum_out = sum(
                    _element_dict[element][j] * self.flow_mol_comp_outlet[t, j]
                    for j in self.component_list
                )
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

        def ssf(c, s):
            iscale.set_scaling_factor(c, s, overwrite=False)

        sgsf = iscale.set_and_get_scaling_factor

        def cst(c, s):
            iscale.constraint_scaling_transform(c, s, overwrite=False)

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
            for j in self.component_list:
                sy_in[j] = sgsf(self.mole_frac_comp_inlet[t, j], sy_def)

            for iz in self.iznodes:
                # These should have been scaled by the cell-level method, so
                # notify the user if they're using a standalone channel
                # and forgot to scale these
                if not self.temperature_z[t, iz].is_reference():
                    gsf(self.temperature_z[t, iz], warning=True)
                gsf(self.heat_flux_x0[t, iz], warning=True)
                gsf(self.heat_flux_x1[t, iz], warning=True)

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

                sq0 = sgsf(self.heat_flux_x0[t, iz], 1e-2)
                cst(self.temperature_x0_eqn[t, iz], sq0)
                sq1 = sgsf(self.heat_flux_x1[t, iz], 1e-2)
                cst(self.temperature_x1_eqn[t, iz], sq1)
                sq = min(sq0, sq1)

                sDT = sgsf(self.temperature_deviation_x[t, iz], sq / sh)
                ssf(self.temperature_deviation_x0[t, iz], sDT)
                ssf(self.temperature_deviation_x1[t, iz], sDT)

                # Pointless other than making a record that this equation
                # is well-scaled by default
                cst(self.mole_frac_comp_eqn[t, iz], 1)

                for j in self.component_list:
                    # These should have been scaled by the cell-level method, so
                    # notify the user if they're using a standalone channel
                    # and forgot to scale these
                    if self.config.below_electrode:
                        gsf(self.material_flux_x1[t, iz, j], warning=True)
                    else:
                        gsf(self.material_flux_x0[t, iz, j], warning=True)

                    sy = sgsf(self.mole_frac_comp[t, iz, j], sy_in[j])
                    cst(self.material_balance_eqn[t, iz, j], s_flow_mol * sy)

                    ssf(self.conc_mol_comp[t, iz, j], sy * sP / (sR * sT))
                    cst(self.conc_mol_comp_eqn[t, iz, j], sy * sP)
                    # FIXME come back later to clean up
                    ssf(self.diff_eff_coeff[t, iz, j], 1e5)
                    cst(self.diff_eff_coeff_eqn[t, iz, j], 1e5)

                    if hasattr(self, "material_flux_x0_eqn"):
                        sXflux = gsf(
                            self.material_flux_x0[t, iz, j], default=1e-1, warning=True
                        )
                        cst(self.material_flux_x0_eqn[t, iz, j], sLx * sXflux / sD)
                        ssf(
                            self.conc_mol_comp_deviation_x0[t, iz, j], sLx * sXflux / sD
                        )
                    if hasattr(self, "material_flux_x1_eqn"):
                        sXflux = gsf(
                            self.material_flux_x1[t, iz, j], default=1e-1, warning=True
                        )
                        cst(self.material_flux_x1_eqn[t, iz, j], sLx * sXflux / sD)
                        ssf(
                            self.conc_mol_comp_deviation_x1[t, iz, j], sLx * sXflux / sD
                        )

            cst(self.temperature_outlet_eqn[t], sH)
