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
Unit model to scale SOC variables to and from cell-scale to module-scale and
translate between the FTPx variables assumed at the cell level to a general
property package. The only module-level variable, besides those contained
in the StateBlocks, is ``number_cells``.

A config block for the SOC unit model must be provided. Further documentation
is provided in the cell model file.

Manipulated variables:
    - ``potential_cell[t]``: Voltage across single cell

Ports:
    * ``fuel_inlet``
    * ``fuel_outlet``
    * ``oxygen_inlet``
    * ``oxygen_outlet``

Instances of ``Var`` that must be fixed:
    - ``number_cells``: Number of cells in module

Expressions:
    - ``electrical_work[t]``: Rate of energy added to module. Greater than zero means energy added to module
      (electrolysis mode) and less than zero means energy removed from module (fuel cell mode)
"""
__author__ = "Douglas Allan"

from pyomo.common.config import ConfigValue, ConfigBlock
import pyomo.environ as pyo

from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.core.util.config import is_physical_parameter_block
import idaes.models_extra.power_generation.unit_models.soc_submodels as soc
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


@declare_process_block_class("SolidOxideModuleSimple")
class SolidOxideModuleSimpleData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    # CONFIG = ConfigBlock()
    CONFIG.declare(
        "solid_oxide_cell_config",
        ConfigBlock(
            implicit=True,
            description="Config block for the SolidOxideCell",
        ),
    )
    CONFIG.declare(
        "fuel_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for fuel-side stream",
            doc="""Property parameter object used to define property
                    calculations for the fuel-side stream,
                    **default** - None.
                    **Valid values:** {
                    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "fuel_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package "
            "of the fuel-side stream",
            doc="""A ConfigBlock with arguments to be passed to the property
                    block associated with the fuel-side stream,
                    **default** - None.
                    **Valid values:** {
                    see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "oxygen_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for oxygen-side stream",
            doc="""Property parameter object used to define property
                        calculations for the oxygen-side stream,
                        **default** - None.
                        **Valid values:** {
                        **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "oxygen_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package "
            "of the oxygen-side stream",
            doc="""A ConfigBlock with arguments to be passed to the property
                        block associated with the oxygen-side stream,
                        **default** - None.
                        **Valid values:** {
                        see property package for documentation.}""",
        ),
    )

    def build(self):
        super().build()
        has_holdup = self.config.has_holdup
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        t0 = tset.first()

        self.number_cells = pyo.Var(
            initialize=1e5,
            units=pyo.units.dimensionless,
            doc="Number of cells in SOC module",
        )

        self.solid_oxide_cell = soc.SolidOxideCell(
            **self.config.solid_oxide_cell_config
        )

        def rule_absent_comp(blk, t, j, props):
            return props[t].flow_mol_comp[j] == 1e-20

        def rule_temperature(blk, t, props, port):
            return props[t].temperature == port.temperature[t]

        def rule_pressure(blk, t, props, port):
            return props[t].pressure == port.pressure[t]

        def rule_flow_mol_comp(blk, t, j, props, port):
            return (
                props[t].flow_mol_comp[j]
                == blk.number_cells * port.flow_mol[t] * port.mole_frac_comp[t, j]
            )

        def rule_mole_frac(blk, t, port, comps):
            return sum(port.mole_frac_comp[t, j] for j in comps) == 1

        for side in ["fuel", "oxygen"]:
            param_block = getattr(self.config, f"{side}_property_package")
            package_args = getattr(self.config, f"{side}_property_package_args")
            setattr(
                self,
                f"{side}_properties_in",
                param_block.build_state_block(
                    self.flowsheet().time,
                    doc=f"Material properties in {side}-side inlet stream",
                    defined_state=True,
                    has_phase_equilibrium=False,
                    **package_args,
                ),
            )
            setattr(
                self,
                f"{side}_properties_out",
                param_block.build_state_block(
                    self.flowsheet().time,
                    doc=f"Material properties in {side}-side outlet stream",
                    defined_state=False,
                    has_phase_equilibrium=False,
                    **package_args,
                ),
            )
            side_comps = getattr(self.solid_oxide_cell, f"{side}_component_list")
            # Support absent components in property package to support legacy property packages that do not
            # permit users to exclude components.
            absent_comp_list = common._set_and_get_attr(
                self,
                f"absent_{side}_component_list",
                pyo.Set(
                    initialize=param_block.component_list - side_comps,
                    ordered=True,
                    doc=f"Components in the {side}-side property package that are not actually present in the stream.",
                ),
            )

            for direction in ["in", "out"]:
                props = getattr(self, f"{side}_properties_{direction}")
                port_name = f"{side}_{direction}let"  # e.g., fuel_inlet, oxygen_outlet
                port = getattr(self.solid_oxide_cell, port_name)
                setattr(
                    self,
                    f"{port_name}_temperature_eqn",
                    pyo.Constraint(
                        tset, rule=lambda blk, t: rule_temperature(blk, t, props, port)
                    ),
                )
                setattr(
                    self,
                    f"{port_name}_pressure_eqn",
                    pyo.Constraint(
                        tset, rule=lambda blk, t: rule_pressure(blk, t, props, port)
                    ),
                )
                setattr(
                    self,
                    f"{port_name}_flow_mol_comp_eqn",
                    pyo.Constraint(
                        tset,
                        side_comps,
                        rule=lambda blk, t, j: rule_flow_mol_comp(
                            blk, t, j, props, port
                        ),
                    ),
                )
                if direction == "in":
                    setattr(
                        self,
                        f"{port_name}_mole_frac_eqn",
                        pyo.Constraint(
                            tset,
                            rule=lambda blk, t: rule_mole_frac(
                                blk, t, port, side_comps
                            ),
                        ),
                    )
                if direction == "out":
                    setattr(
                        self,
                        f"{port_name}_absent_comp_eqn",
                        pyo.Constraint(
                            tset,
                            absent_comp_list,
                            rule=lambda blk, t, j: rule_absent_comp(blk, t, j, props),
                        ),
                    )
                # Add a different port at the module level
                self.add_port(
                    name=port_name,
                    block=props,
                    doc=f"{side.capitalize()}-side {direction.capitalize()}let Port",
                )

        # This is net flow of power *into* the cell. In fuel cell mode, this will
        # be negative
        @self.Expression(tset)
        def electrical_work(b, t):
            return b.solid_oxide_cell.electrical_work[t] * b.number_cells

        self.potential_cell = pyo.Reference(self.solid_oxide_cell.potential)

    def initialize_build(
        self,
        state_args_fuel=None,
        state_args_oxygen=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        current_density_guess=None,
        temperature_guess=None,
    ):
        t0 = self.flowsheet().time.first()

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        tset = self.flowsheet().config.time

        solver_obj = get_solver(solver, optarg)

        number_cells_fixed = self.number_cells.fixed
        potential_cell_fixed = {}
        for t in self.flowsheet().time:
            potential_cell_fixed[t] = self.potential_cell[t].fixed
            self.potential_cell[t].fix()

        self.number_cells.fix()

        flags = {}
        for side, state_args in zip(
            ["fuel", "oxygen"],
            [state_args_fuel, state_args_oxygen],
        ):
            props_in = getattr(self, side + "_properties_in")
            inlet = getattr(self.solid_oxide_cell, side + "_inlet")
            comps = getattr(self.solid_oxide_cell, side + "_component_list")
            flags[side] = props_in.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args,
                hold_state=True,
            )
            for t in self.flowsheet().time:
                inlet.temperature[t].value = pyo.value(props_in[t].temperature)
                inlet.pressure[t].value = pyo.value(props_in[t].pressure)
                inlet.flow_mol[t].value = pyo.value(
                    props_in[t].flow_mol / self.number_cells
                )
                for j in comps:
                    inlet.mole_frac_comp[t, j].value = pyo.value(
                        props_in[t].mole_frac_comp[j]
                    )

        self.solid_oxide_cell.initialize_build(
            outlvl=idaeslog.NOTSET,
            solver=solver,
            optarg=optarg,
            current_density_guess=current_density_guess,
            temperature_guess=temperature_guess,
        )

        # Without knowing the state variables of the property package, the best we can do is just initialize
        # the property blocks and then solve the entire block
        for side, state_args in zip(
            ["fuel", "oxygen"],
            [state_args_fuel, state_args_oxygen],
        ):
            props_out = getattr(self, side + "_properties_out")
            props_out.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args,
            )
        common._init_solve_block(self, solver_obj, solve_log)

        init_log.info(f"{self.name} initialization complete.")

        self.fuel_properties_in.release_state(flags=flags["fuel"], outlvl=outlvl)
        self.oxygen_properties_in.release_state(flags=flags["oxygen"], outlvl=outlvl)
        if not number_cells_fixed:
            self.number_cells.unfix()
        for t in self.flowsheet().time:
            if not potential_cell_fixed[t]:
                self.potential_cell[t].unfix()

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        gsf = iscale.get_scaling_factor

        def ssf(c, s):
            iscale.set_scaling_factor(c, s, overwrite=False)

        def cst(c, s):
            return iscale.constraint_scaling_transform(c, s, overwrite=False)

        s_num_cells = gsf(self.number_cells, default=1 / self.number_cells.value)

        for side in ["fuel", "oxygen"]:
            comps = getattr(self.solid_oxide_cell, f"{side}_component_list")
            channel = getattr(self.solid_oxide_cell, f"{side}_channel")
            for direction in ["in", "out"]:
                props = getattr(self, f"{side}_properties_{direction}")
                port_name = f"{side}_{direction}let"  # e.g., fuel_inlet, oxygen_outlet

                temperature_eqn = getattr(self, f"{port_name}_temperature_eqn")
                pressure_eqn = getattr(self, f"{port_name}_pressure_eqn")
                flow_mol_comp_eqn = getattr(self, f"{port_name}_flow_mol_comp_eqn")
                if direction == "out":
                    absent_comp_eqn = getattr(self, f"{port_name}_absent_comp_eqn")
                    absent_comp_list = getattr(self, f"absent_{side}_component_list")
                for t in self.flowsheet().time:
                    sT = gsf(props[t].temperature)
                    cst(temperature_eqn[t], sT)
                    sP = gsf(props[t].pressure)
                    cst(pressure_eqn[t], sP)
                    sflow = gsf(props[t].flow_mol)

                    if direction == "in":
                        # To avoid Reference issues, use channel variables directly
                        ssf(channel.temperature_inlet[t], sT)
                        ssf(channel.pressure_inlet[t], sP)
                        ssf(channel.flow_mol_inlet[t], sflow / s_num_cells)

                    for j in comps:
                        sx = gsf(props[t].mole_frac_comp[j])
                        cst(flow_mol_comp_eqn[t, j], sflow * sx)
                        if direction == "in":
                            ssf(channel.mole_frac_comp_inlet[t, j], sx)

                    if direction == "out":
                        for j in absent_comp_list:
                            sx = gsf(props[t].mole_frac_comp[j])
                            cst(absent_comp_eqn[t, j], sflow * sx)

        self.solid_oxide_cell.recursive_scaling()

    def model_check(self):
        pass
