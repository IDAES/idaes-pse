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
Model for nonisothermal (single) solid oxide fuel/electrolytic cell. Because
the cell model is reversible, we do not refer to electrodes as "anode" or "cathode",
whose designation would depend on whether the cell was in fuel cell or electrolytic
mode, but rather the but rather as "fuel" or "oxygen" electrodes, being the electrodes
at which fuel or oxygen, respectively, is produced or consumed (depending on the mode).

This model has at least seven submodels:
    - Fuel/oxygen channels
    - Fuel/oxygen electrodes (both ``PorousConductiveSlab`` objects)
    - Fuel/oxygen triple phase boundaries
    - Electrolyte (a ``ConductiveSlab``)
    - (Optional) Contacts between the flow mesh and the interconnect, and between the flow mesh
      and the electrodes (four total contacts)
    - (Optional and not implemented) Interconnect between neighboring cells within the stack
      (also a ``ConductiveSlab``)

This model is two-dimensional.
    - The z-axis is oriented along the direction of flow in the fuel channel, with the fuel
      flow entering at ``z=0`` and exiting at ``z=1``, regardless of whether the cell operates in
      counter-current or co-current mode. If operated in counter-current mode, the oxygen
      enters at ``z=1`` and exists at ``z=0``.
    - The x-axis starts at the fuel channel and ends at the oxygen channel. Each two-dimensional
      submodel has an individual x-axis starting at ``x=0`` at its bottom and ending at ``x=1`` at
      its top. So the fuel electrodes ``x=0`` boundary is located at the fuel channel's ``x=1``
      boundary, the electrolyte's ``x=0`` boundary is located at the fuel electrode's ``x=1``
      boundary, the oxygen electrode's ``x=0`` boundary is located at the electrolyte's ``x=1``
      boundary, and the oxygen channel's ``x=0`` boundary is located at the oxygen electrode's
      ``x=0`` boundary.

Flux terms are positive if a flow occurs from zero to one direction, and negative if a flow
occurs from the one to zero direction. For example, if a positive heat flux term occurs at
the ``x=1`` boundary of the fuel channel, it means heat is flowing from the fuel channel to the
fuel electrode.

The triple phase boundary and the contact resistor terms are "thin". They have only a z-axis,
with no x-axis. However, they do have different flux terms at their ``x=0`` and ``x=1`` sides,
which represent resistive heating and (for the triple phase boundaries) chemical reactions
occurring at those interfaces. Channel models are not "thin", so they have temperatures and
concentrations at their ``x=0`` and ``x=1`` edges different from those in the bulk. However, they
are not discretized in the x direction further. The other submodels have both x and z-axes.

Boundary conditions between submodels are typically handled by adding variables from one
submodel into the Config of another. If a variable is passed to a submodel, a Reference
is created on the other submodel so that it can be used in expressions and constraints.
If no variable is provided in the Config, then a local Var is created automatically.

Because the SOC is extremely thin in the x direction, temperature and concentration variations
in the x direction are small, especially in the electrolytes and electrodes. However,
temperature and concentration gradients in the x direction can be large. To avoid catastrophic
cancellation and ill-conditioning when calculating gradients, temperature variations in the x
direction are expressed as deviations from the average of the fuel and oxygen channel
temperatures at a particular location in the z direction. Concentration values in the electrodes
are expressed as deviations from the concentrations in their respective channels.

This model implements ideal gas properties for certain components, listed in
``common._gas_species_list``. The user must specify which components to include on the fuel and
oxygen sides, using ``fuel_component_list`` and oxygen_component_list. In addition to these gas
phase species, certain solid-phase species are included in  ``common._all_species_list``. These
species can be used for specifying stoichiometry at the fuel and oxygen triple-phase boundaries.
When specifying stoichiometry, an ``"e^-"`` term for electrons produced or consumed by the reaction
must be included.

Submodel PDEs are discretized using a finite volume method. The ``control_volume_zfaces`` config
argument determines the control volume faces in the z direction for all submodels, but
``control_volume_xfaces`` must be supplied separately for the electrodes and electrolyte.

More detailed descriptions of the submodels are found in their respective files. At the cell level,
``temperature_z``, ``current_density``, and ``cell_potential`` variables are created, along with average
temperature equations, voltage drop equations, and boundary conditions at the ``x=0`` end of the fuel
channel and ``x=1`` end of the oxygen channel (no flux conditions until the interconnect is implemented).

Manipulated variables:
    - ``potential[t]``: Voltage across single cell

Ports:
    * ``fuel/oxygen_inlet``: Lifts channel level inlet variables to cell level.
        - ``flow_mol[t]``: Molar flow into channel
        - ``mole_frac_comp[t, j]``: Molar composition of flow into channel
        - ``temperature[t]``: Channel inlet flow rate
        - ``pressure[t]``: Channel inlet pressure
    * ``fuel/oxygen_outlet``: Lifts channel level outlet variables to cell level.
        - ``flow_mol[t]``: Molar flow out of channel
        - ``mole_frac_comp[t, j]``: Molar composition of flow from channel
        - ``temperature[t]``: Channel outlet flow rate
        - ``pressure[t]``: Channel outlet pressure

Instances of ``Var`` that must be fixed:
    - ``length_z``: Length of cell along direction of flow
    - ``length_y``: Width of cell

Expressions:
    - ``electrical_work[t]``: Rate of energy added to cell. Greater than zero means energy added to cell
      (electrolysis mode) and less than zero means energy removed from cell (fuel cell mode)
"""
__author__ = "John Eslick, Douglas Allan"

from pyomo.common.config import ConfigValue, In, Bool, ListOf
import pyomo.environ as pyo
from pyomo.network import Port

from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern
import idaes.models_extra.power_generation.unit_models.soc_submodels as soc
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
from idaes.core.util.constants import Constants
from idaes.models_extra.power_generation.unit_models.soc_submodels.common import (
    _gas_species_list,
    _all_species_list,
    _element_list,
    _element_dict,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


@declare_process_block_class("SolidOxideCell")
class SolidOxideCellData(UnitModelBlockData):
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
        "control_volume_zfaces",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates of control volume faces "
            "in z direction. Coordinates must start with zero, be strictly "
            "increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "control_volume_xfaces_fuel_electrode",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates of control volume faces "
            "in x direction for fuel electrode. Coordinates must start with "
            "zero, be strictly increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "control_volume_xfaces_oxygen_electrode",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates of control volume faces "
            "in x direction for oxygen electrode. Coordinates must start with "
            "zero, be strictly increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "control_volume_xfaces_electrolyte",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates of control volume faces "
            "in x direction for electrolyte. Coordinates must start with "
            "zero, be strictly increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "thin_fuel_electrode",
        ConfigValue(
            domain=Bool,
            default=False,
            description="If True, make 1D fuel electrode",
        ),
    )
    CONFIG.declare(
        "thin_oxygen_electrode",
        ConfigValue(
            domain=Bool,
            default=False,
            description="If True, make 1D oxygen electrode",
        ),
    )
    CONFIG.declare(
        "thin_electrolyte",
        ConfigValue(
            domain=Bool,
            default=False,
            description="If True, make 1D electrolyte",
        ),
    )
    CONFIG.declare(
        "fuel_component_list",
        ConfigValue(
            domain=common._SubsetOf(_gas_species_list),
            default=None,
            description="List of components in fuel stream",
        ),
    )
    CONFIG.declare(
        "oxygen_component_list",
        ConfigValue(
            domain=common._SubsetOf(_gas_species_list),
            default=None,
            description="List of components in the oxygen stream",
        ),
    )
    CONFIG.declare(
        "fuel_triple_phase_boundary_stoich_dict",
        ConfigValue(
            domain=common._SubsetOf(_all_species_list),
            default=None,
            description="Dictionary with species as keys and stoichiometric coefficients as values "
            "for the redox reaction that occurs at the fuel-side triple phase boundary",
        ),
    )
    CONFIG.declare(
        "inert_fuel_species_triple_phase_boundary",
        ConfigValue(
            domain=common._SubsetOf(_gas_species_list),
            default=None,
            description="List of fuel-side species that do not participate in "
            "reactions at the triple phase boundary."
            # But may be involved in reforming
        ),
    )
    CONFIG.declare(
        "oxygen_triple_phase_boundary_stoich_dict",
        ConfigValue(
            domain=common._SubsetOf(_all_species_list),
            default=None,
            description="Dictionary with species as keys and stoichiometric coefficients as values "
            "for the redox reaction that occurs at the oxygen-side triple phase boundary",
        ),
    )
    CONFIG.declare(
        "inert_oxygen_species_triple_phase_boundary",
        ConfigValue(
            domain=common._SubsetOf(_gas_species_list),
            default=None,
            description="List of oxygen-side species that do not participate in "
            "reactions at the triple phase boundary.",
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
            domain=Bool,
            description="If True write periodic constraint "
            "to model flux through interconnect.",
        ),
    )
    CONFIG.declare(
        "thin_interconnect",
        ConfigValue(
            domain=Bool,
            default=False,
            description="If True, make 1D interconnect",
        ),
    )
    CONFIG.declare(
        "control_volume_xfaces_interconnect",
        ConfigValue(
            domain=ListOf(float),
            description="List containing coordinates of control volume faces "
            "in x direction for interconnect. Coordinates must start with "
            "zero, be strictly increasing, and end with one",
        ),
    )
    CONFIG.declare(
        "include_contact_resistance",
        ConfigValue(
            default=False,
            domain=Bool,
            description="If True write periodic constraint "
            "to model flux through interconnect.",
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
        has_gas_holdup = self.config.has_gas_holdup
        dynamic = self.config.dynamic
        tset = self.flowsheet().config.time
        t0 = tset.first()

        if has_gas_holdup and not has_holdup:
            raise ConfigurationError(
                "Creating gas holdup terms while not creating other holdup terms is not supported."
            )

        if self.config.fuel_component_list is None:
            fuel_comp_list = ["H2", "H2O"]
        else:
            fuel_comp_list = self.config.fuel_component_list
        if self.config.oxygen_component_list is None:
            oxygen_comp_list = ["O2"]
        else:
            oxygen_comp_list = self.config.oxygen_component_list
        if self.config.fuel_triple_phase_boundary_stoich_dict is None:
            fuel_tpb_stoich = {"H2": -0.5, "H2O": 0.5, "e^-": 1.0}
        else:
            fuel_tpb_stoich = self.config.fuel_triple_phase_boundary_stoich_dict
        if self.config.oxygen_triple_phase_boundary_stoich_dict is None:
            oxygen_tpb_stoich = {"O2": -0.25, "e^-": -1.0}
        else:
            oxygen_tpb_stoich = self.config.oxygen_triple_phase_boundary_stoich_dict
        self.fuel_component_list = pyo.Set(
            initialize=fuel_comp_list,
            ordered=True,
            doc="Set of all gas-phase components present on fuel side of cell",
        )
        self.oxygen_component_list = pyo.Set(
            initialize=oxygen_comp_list,
            ordered=True,
            doc="Set of all gas-phase components present on oxygen side of cell",
        )
        # Set up node and face sets and get integer indices for them
        izfaces, iznodes = common._face_initializer(
            self, self.config.control_volume_zfaces, "z"
        )
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
            bounds=(300, None),
        )
        self.length_z = pyo.Var(
            doc="Length of cell in direction parallel to channel flow.",
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

        self.fuel_channel = soc.SocChannel(
            has_holdup=has_gas_holdup,
            dynamic=has_gas_holdup and dynamic,
            control_volume_zfaces=self.config.control_volume_zfaces,
            length_z=self.length_z,
            length_y=self.length_y,
            temperature_z=self.temperature_z,
            opposite_flow=False,  # Fuel channel flows from z=0 to z=1 no matter what
            below_electrode=True,
            component_list=self.fuel_component_list,
            include_temperature_x_thermo=include_temp_x_thermo,
        )
        self.oxygen_channel = soc.SocChannel(
            has_holdup=has_gas_holdup,
            dynamic=has_gas_holdup and dynamic,
            control_volume_zfaces=self.config.control_volume_zfaces,
            length_z=self.length_z,
            length_y=self.length_y,
            temperature_z=self.temperature_z,
            below_electrode=False,
            opposite_flow=opposite_flow,  # Oxygen channel flow direction depends on config
            component_list=self.oxygen_component_list,
            include_temperature_x_thermo=include_temp_x_thermo,
        )
        if self.config.include_contact_resistance:
            self.contact_interconnect_fuel_flow_mesh = soc.SocContactResistor(
                has_holdup=False,
                dynamic=False,
                control_volume_zfaces=self.config.control_volume_zfaces,
                length_z=self.length_z,
                length_y=self.length_y,
                temperature_z=self.temperature_z,
                temperature_deviation_x=self.fuel_channel.temperature_deviation_x0,
                heat_flux_x1=self.fuel_channel.heat_flux_x0,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )
            self.contact_flow_mesh_fuel_electrode = soc.SocContactResistor(
                has_holdup=False,
                dynamic=False,
                control_volume_zfaces=self.config.control_volume_zfaces,
                length_z=self.length_z,
                length_y=self.length_y,
                temperature_z=self.temperature_z,
                temperature_deviation_x=self.fuel_channel.temperature_deviation_x1,
                heat_flux_x0=self.fuel_channel.heat_flux_x1,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )
            self.contact_interconnect_oxygen_flow_mesh = soc.SocContactResistor(
                has_holdup=False,
                dynamic=False,
                control_volume_zfaces=self.config.control_volume_zfaces,
                length_z=self.length_z,
                length_y=self.length_y,
                temperature_z=self.temperature_z,
                temperature_deviation_x=self.oxygen_channel.temperature_deviation_x1,
                heat_flux_x0=self.oxygen_channel.heat_flux_x1,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )
            self.contact_flow_mesh_oxygen_electrode = soc.SocContactResistor(
                has_holdup=False,
                dynamic=False,
                control_volume_zfaces=self.config.control_volume_zfaces,
                length_z=self.length_z,
                length_y=self.length_y,
                temperature_z=self.temperature_z,
                temperature_deviation_x=self.oxygen_channel.temperature_deviation_x0,
                heat_flux_x1=self.oxygen_channel.heat_flux_x0,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )
            fuel_electrode_heat_flux_x0 = (
                self.contact_flow_mesh_fuel_electrode.heat_flux_x1
            )
            oxygen_electrode_heat_flux_x1 = (
                self.contact_flow_mesh_oxygen_electrode.heat_flux_x0
            )
            interconnect_heat_flux_x0 = (
                self.contact_interconnect_oxygen_flow_mesh.heat_flux_x1
            )
            interconnect_heat_flux_x1 = (
                self.contact_interconnect_fuel_flow_mesh.heat_flux_x0
            )
        else:
            fuel_electrode_heat_flux_x0 = self.fuel_channel.heat_flux_x1
            oxygen_electrode_heat_flux_x1 = self.oxygen_channel.heat_flux_x0
            interconnect_heat_flux_x0 = self.oxygen_channel.heat_flux_x1
            interconnect_heat_flux_x1 = self.fuel_channel.heat_flux_x0
        if self.config.thin_fuel_electrode:
            if self.config.control_volume_xfaces_fuel_electrode is not None:
                raise ConfigurationError(
                    "User specified thin fuel electrode, but provided xfaces."
                )
            self.fuel_electrode = soc.SocContactResistor(
                has_holdup=False,
                dynamic=False,
                control_volume_zfaces=self.config.control_volume_zfaces,
                length_z=self.length_z,
                length_y=self.length_y,
                temperature_z=self.temperature_z,
                temperature_deviation_x=self.fuel_channel.temperature_deviation_x1,
                heat_flux_x0=fuel_electrode_heat_flux_x0,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )
            # Technically I could use fuel_electrode.temperature_deviation_x, but I want to avoid
            # references to references
            fuel_electrode_temperature_deviation_x1 = (
                self.fuel_channel.temperature_deviation_x1
            )
            fuel_electrode_conc_mol_comp_deviation_x1 = (
                self.fuel_channel.conc_mol_comp_deviation_x1
            )
            fuel_electrode_material_flux_x1 = self.fuel_channel.material_flux_x1
        else:
            self.fuel_electrode = soc.PorousConductiveSlab(
                has_holdup=has_holdup,
                has_gas_holdup=has_gas_holdup,
                dynamic=dynamic,
                control_volume_zfaces=self.config.control_volume_zfaces,
                control_volume_xfaces=self.config.control_volume_xfaces_fuel_electrode,
                component_list=self.fuel_channel.component_list,
                length_z=self.length_z,
                length_y=self.length_y,
                conc_mol_comp_ref=self.fuel_channel.conc_mol_comp,
                dconc_mol_comp_refdt=self.fuel_channel.dconc_mol_compdt,
                conc_mol_comp_deviation_x0=self.fuel_channel.conc_mol_comp_deviation_x1,
                material_flux_x0=self.fuel_channel.material_flux_x1,
                heat_flux_x0=fuel_electrode_heat_flux_x0,
                temperature_z=self.temperature_z,
                temperature_deviation_x0=self.fuel_channel.temperature_deviation_x1,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )
            fuel_electrode_temperature_deviation_x1 = (
                self.fuel_electrode.temperature_deviation_x1
            )
            fuel_electrode_conc_mol_comp_deviation_x1 = (
                self.fuel_electrode.conc_mol_comp_deviation_x1
            )
            fuel_electrode_material_flux_x1 = self.fuel_electrode.material_flux_x1

        if self.config.thin_oxygen_electrode:
            if self.config.control_volume_xfaces_oxygen_electrode is not None:
                raise ConfigurationError(
                    "User specified thin oxygen electrode, but provided xfaces."
                )
            self.oxygen_electrode = soc.SocContactResistor(
                has_holdup=False,
                dynamic=False,
                control_volume_zfaces=self.config.control_volume_zfaces,
                length_z=self.length_z,
                length_y=self.length_y,
                temperature_z=self.temperature_z,
                temperature_deviation_x=self.oxygen_channel.temperature_deviation_x0,
                heat_flux_x1=oxygen_electrode_heat_flux_x1,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )
            # Technically I could use oxygen_electrode.temperature_deviation_x, but I want to avoid
            # references to references
            oxygen_electrode_temperature_deviation_x0 = (
                self.oxygen_channel.temperature_deviation_x0
            )
            oxygen_electrode_conc_mol_comp_deviation_x0 = (
                self.oxygen_channel.conc_mol_comp_deviation_x0
            )
            oxygen_electrode_material_flux_x0 = self.oxygen_channel.material_flux_x0
        else:
            self.oxygen_electrode = soc.PorousConductiveSlab(
                has_holdup=has_holdup,
                has_gas_holdup=has_gas_holdup,
                dynamic=dynamic,
                control_volume_zfaces=self.config.control_volume_zfaces,
                control_volume_xfaces=self.config.control_volume_xfaces_oxygen_electrode,
                component_list=self.oxygen_channel.component_list,
                length_z=self.length_z,
                length_y=self.length_y,
                conc_mol_comp_ref=self.oxygen_channel.conc_mol_comp,
                dconc_mol_comp_refdt=self.oxygen_channel.dconc_mol_compdt,
                conc_mol_comp_deviation_x1=self.oxygen_channel.conc_mol_comp_deviation_x0,
                material_flux_x1=self.oxygen_channel.material_flux_x0,
                heat_flux_x1=oxygen_electrode_heat_flux_x1,
                temperature_z=self.temperature_z,
                temperature_deviation_x1=self.oxygen_channel.temperature_deviation_x0,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )
            oxygen_electrode_temperature_deviation_x0 = (
                self.oxygen_electrode.temperature_deviation_x0
            )
            oxygen_electrode_conc_mol_comp_deviation_x0 = (
                self.oxygen_electrode.conc_mol_comp_deviation_x0
            )
            oxygen_electrode_material_flux_x0 = self.oxygen_electrode.material_flux_x0

        self.fuel_triple_phase_boundary = soc.SocTriplePhaseBoundary(
            has_holdup=False,
            dynamic=False,
            control_volume_zfaces=self.config.control_volume_zfaces,
            length_z=self.length_z,
            length_y=self.length_y,
            component_list=self.fuel_component_list,
            reaction_stoichiometry=fuel_tpb_stoich,
            inert_species=self.config.inert_fuel_species_triple_phase_boundary,
            current_density=self.current_density,
            temperature_z=self.temperature_z,
            temperature_deviation_x=fuel_electrode_temperature_deviation_x1,
            heat_flux_x0=self.fuel_electrode.heat_flux_x1,
            conc_mol_comp_ref=self.fuel_channel.conc_mol_comp,
            conc_mol_comp_deviation_x=fuel_electrode_conc_mol_comp_deviation_x1,
            material_flux_x=fuel_electrode_material_flux_x1,
            include_temperature_x_thermo=include_temp_x_thermo,
            below_electrolyte=True,
        )
        self.oxygen_triple_phase_boundary = soc.SocTriplePhaseBoundary(
            has_holdup=False,
            dynamic=False,
            control_volume_zfaces=self.config.control_volume_zfaces,
            length_z=self.length_z,
            length_y=self.length_y,
            component_list=self.oxygen_component_list,
            reaction_stoichiometry=oxygen_tpb_stoich,
            inert_species=self.config.inert_oxygen_species_triple_phase_boundary,
            current_density=self.current_density,
            temperature_z=self.temperature_z,
            temperature_deviation_x=oxygen_electrode_temperature_deviation_x0,
            heat_flux_x1=self.oxygen_electrode.heat_flux_x0,
            conc_mol_comp_ref=self.oxygen_channel.conc_mol_comp,
            conc_mol_comp_deviation_x=oxygen_electrode_conc_mol_comp_deviation_x0,
            material_flux_x=oxygen_electrode_material_flux_x0,
            include_temperature_x_thermo=include_temp_x_thermo,
            below_electrolyte=False,
        )
        if self.config.thin_electrolyte:
            if self.config.control_volume_xfaces_electrolyte is not None:
                raise ConfigurationError(
                    "User specified thin electrolyte, but provided xfaces."
                )
            self.electrolyte = soc.SocContactResistor(
                has_holdup=False,
                dynamic=False,
                control_volume_zfaces=self.config.control_volume_zfaces,
                length_z=self.length_z,
                length_y=self.length_y,
                temperature_z=self.temperature_z,
                temperature_deviation_x=fuel_electrode_temperature_deviation_x1,
                heat_flux_x0=self.fuel_triple_phase_boundary.heat_flux_x1,
                heat_flux_x1=self.oxygen_triple_phase_boundary.heat_flux_x0,
                current_density=self.current_density,
                include_temperature_x_thermo=include_temp_x_thermo,
            )

            @self.Constraint(tset, iznodes)
            def electrolyte_temperature_continuity_eqn(b, t, iz):
                return (
                    fuel_electrode_temperature_deviation_x1[t, iz]
                    == oxygen_electrode_temperature_deviation_x0[t, iz]
                )

        else:
            self.electrolyte = soc.SocConductiveSlab(
                has_holdup=has_holdup,
                dynamic=dynamic,
                control_volume_zfaces=self.config.control_volume_zfaces,
                control_volume_xfaces=self.config.control_volume_xfaces_electrolyte,
                length_z=self.length_z,
                length_y=self.length_y,
                temperature_z=self.temperature_z,
                current_density=self.current_density,
                temperature_deviation_x0=fuel_electrode_temperature_deviation_x1,
                temperature_deviation_x1=oxygen_electrode_temperature_deviation_x0,
                heat_flux_x0=self.fuel_triple_phase_boundary.heat_flux_x1,
                heat_flux_x1=self.oxygen_triple_phase_boundary.heat_flux_x0,
                include_temperature_x_thermo=include_temp_x_thermo,
            )

        self.state_vars = {"flow_mol", "mole_frac_comp", "temperature", "pressure"}
        for chan, alias in zip(
            [self.fuel_channel, self.oxygen_channel], ["fuel", "oxygen"]
        ):
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

        if self.config.flux_through_interconnect:
            if self.config.thin_interconnect:
                if self.config.control_volume_xfaces_interconnect is not None:
                    raise ConfigurationError(
                        "User specified thin electrolyte, but provided xfaces."
                    )
                self.interconnect = soc.SocContactResistor(
                    has_holdup=False,
                    dynamic=False,
                    control_volume_zfaces=self.config.control_volume_zfaces,
                    length_z=self.length_z,
                    length_y=self.length_y,
                    temperature_z=self.temperature_z,
                    temperature_deviation_x=self.oxygen_channel.temperature_deviation_x1,
                    heat_flux_x0=interconnect_heat_flux_x0,
                    heat_flux_x1=interconnect_heat_flux_x1,
                    current_density=self.current_density,
                    include_temperature_x_thermo=include_temp_x_thermo,
                )

                @self.Constraint(tset, iznodes)
                def interconnect_temperature_continuity_eqn(b, t, iz):
                    return (
                        self.oxygen_channel.temperature_deviation_x1[t, iz]
                        == self.fuel_channel.temperature_deviation_x0[t, iz]
                    )

            else:
                self.interconnect = soc.SocConductiveSlab(
                    has_holdup=has_holdup,
                    dynamic=dynamic,
                    control_volume_zfaces=self.config.control_volume_zfaces,
                    control_volume_xfaces=self.config.control_volume_xfaces_interconnect,
                    length_z=self.length_z,
                    length_y=self.length_y,
                    temperature_z=self.temperature_z,
                    current_density=self.current_density,
                    temperature_deviation_x0=self.oxygen_channel.temperature_deviation_x1,
                    temperature_deviation_x1=self.fuel_channel.temperature_deviation_x0,
                    heat_flux_x0=interconnect_heat_flux_x0,
                    heat_flux_x1=interconnect_heat_flux_x1,
                    include_temperature_x_thermo=include_temp_x_thermo,
                )
        else:
            interconnect_heat_flux_x0.value = 0
            interconnect_heat_flux_x1.value = 0

            @self.Constraint(tset, iznodes)
            def no_heat_flux_fuel_interconnect_eqn(b, t, iz):
                return (
                    0 * pyo.units.W / pyo.units.m**2
                    == interconnect_heat_flux_x1[t, iz]
                )

            @self.Constraint(tset, iznodes)
            def no_heat_flux_oxygen_interconnect_eqn(b, t, iz):
                return (
                    0 * pyo.units.W / pyo.units.m**2
                    == interconnect_heat_flux_x0[t, iz]
                )

        @self.Constraint(tset, iznodes)
        def mean_temperature_eqn(b, t, iz):
            return (
                0 * pyo.units.K
                == b.fuel_channel.temperature_deviation_x[t, iz]
                + b.oxygen_channel.temperature_deviation_x[t, iz]
            )

        if dynamic:
            self.mean_temperature_eqn[tset.first(), :].deactivate()

        @self.Expression(tset, iznodes)
        def voltage_drop_contact(b, t, iz):
            if b.config.include_contact_resistance:
                return (
                    b.contact_interconnect_fuel_flow_mesh.voltage_drop_total[t, iz]
                    + b.contact_flow_mesh_fuel_electrode.voltage_drop_total[t, iz]
                    + b.contact_interconnect_oxygen_flow_mesh.voltage_drop_total[t, iz]
                    + b.contact_flow_mesh_oxygen_electrode.voltage_drop_total[t, iz]
                )
            else:
                return 0 * pyo.units.V

        @self.Expression(tset, iznodes)
        def voltage_drop_interconnect(b, t, iz):
            if b.config.flux_through_interconnect:
                return b.interconnect.voltage_drop_total[t, iz]
            else:
                return 0 * pyo.units.V

        @self.Expression(tset, iznodes)
        def voltage_drop_ohmic(b, t, iz):
            return (
                b.electrolyte.voltage_drop_total[t, iz]
                + b.fuel_electrode.voltage_drop_total[t, iz]
                + b.oxygen_electrode.voltage_drop_total[t, iz]
                + b.voltage_drop_contact[t, iz]
                + b.voltage_drop_interconnect[t, iz]
            )

        @self.Constraint(tset, iznodes)
        def potential_eqn(b, t, iz):
            return b.potential[t] == (
                b.fuel_triple_phase_boundary.potential_nernst[t, iz]
                + b.oxygen_triple_phase_boundary.potential_nernst[t, iz]
                - (
                    b.voltage_drop_ohmic[t, iz]
                    + b.fuel_triple_phase_boundary.voltage_drop_total[t, iz]
                    + b.oxygen_triple_phase_boundary.voltage_drop_total[t, iz]
                )
            )

        @self.Expression(iznodes)
        def dz(b, iz):
            return b.zfaces.at(iz + 1) - b.zfaces.at(iz)

        @self.Expression(iznodes)
        def xface_area(b, iz):
            return b.length_y * b.length_z * b.dz[iz]

        @self.Expression()
        def cell_area(b):
            return b.length_y * b.length_z

        # This is net flow of power *into* the cell. In fuel cell mode, this will
        # be negative
        @self.Expression(tset)
        def electrical_work(b, t):
            return -b.potential[t] * sum(
                b.xface_area[iz] * b.current_density[t, iz]
                for iz in b.electrolyte.iznodes
            )

        @self.Expression(tset)
        def average_current_density(b, t):
            # z is dimensionless and goes from 0 to 1 so sum(dz) = 1 and the
            # cell is a rectangle hence I don't need to worry about width or
            # length or dividing by the sum of dz and the units are right.
            return sum(b.current_density[t, iz] * b.dz[iz] for iz in b.iznodes)

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        current_density_guess=None,
        temperature_guess=None,
    ):
        t0 = self.flowsheet().time.first()
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        if temperature_guess is None:
            temperature_guess = (
                self.fuel_channel.temperature_inlet[t0].value
                + self.oxygen_channel.temperature_inlet[t0].value
            ) / 2
            init_log.warning(
                f"No guess provided for {self.name} average operating temperature, using average of "
                "channel inlet temperatures instead."
            )
        if current_density_guess is None:
            current_density_guess = 0
            init_log.warning(
                f"No guess provided for {self.name} average current density, using initial guess of zero current density."
            )
        tset = self.flowsheet().config.time
        # t0 = tset.first()

        unfix_inlet_var = {}

        # Save fixedness status of inlet variables and fix them for initialization
        for chan, comps in zip(
            ["fuel", "oxygen"], [self.fuel_component_list, self.oxygen_component_list]
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
        solver_obj = get_solver(solver, optarg)
        for t in tset:
            unfix_potential[t] = not self.potential[t].fixed
        self.potential_eqn.deactivate()
        self.current_density.fix(current_density_guess)

        self.temperature_z.fix(temperature_guess)
        self.mean_temperature_eqn.deactivate()
        if self.config.include_contact_resistance:
            self.contact_interconnect_fuel_flow_mesh.initialize_build(
                outlvl=outlvl, solver=solver, optarg=optarg, fix_heat_flux_x0=True
            )
            self.contact_interconnect_oxygen_flow_mesh.initialize_build(
                outlvl=outlvl, solver=solver, optarg=optarg, fix_heat_flux_x0=False
            )

        # Reset the fluxes to zero in case there are stale values
        init_log.info_high("Initializing Fuel Channel")
        self.fuel_channel.material_flux_x1.fix(0)
        self.fuel_channel.heat_flux_x0.fix(0)
        self.fuel_channel.heat_flux_x1.fix(0)
        self.fuel_channel.initialize_build(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )
        self.fuel_channel.material_flux_x1.unfix()
        self.fuel_channel.heat_flux_x0.unfix()
        self.fuel_channel.heat_flux_x1.unfix()

        init_log.info_high("Initializing Oxygen Channel")
        self.oxygen_channel.material_flux_x0.fix(0)
        self.oxygen_channel.heat_flux_x0.fix(0)
        self.oxygen_channel.heat_flux_x1.fix(0)
        self.oxygen_channel.initialize_build(
            outlvl=outlvl,
            solver=solver,
            optarg=optarg,
        )
        self.oxygen_channel.material_flux_x0.unfix()
        self.oxygen_channel.heat_flux_x0.unfix()
        self.oxygen_channel.heat_flux_x1.unfix()
        if self.config.include_contact_resistance:
            init_log.info_high("Initializing Contact Resistors")

            self.contact_flow_mesh_fuel_electrode.initialize_build(
                fix_heat_flux_x0=True
            )
            self.contact_flow_mesh_oxygen_electrode.initialize_build(
                fix_heat_flux_x0=False
            )

        init_log.info_high("Calculating reaction rate at current guess")
        self.fuel_triple_phase_boundary.temperature_deviation_x.fix(0)
        self.fuel_triple_phase_boundary.conc_mol_comp_deviation_x.fix(0)
        self.fuel_triple_phase_boundary.conc_mol_comp_ref.fix()
        self.fuel_triple_phase_boundary.heat_flux_x1.fix(0)
        for t in self.flowsheet().time:
            for iz in self.fuel_triple_phase_boundary.iznodes:
                denom = pyo.value(
                    sum(
                        self.fuel_triple_phase_boundary.conc_mol_comp[t, iz, j]
                        for j in self.fuel_triple_phase_boundary.component_list
                    )
                )
                for j in self.fuel_triple_phase_boundary.component_list:
                    self.fuel_triple_phase_boundary.mole_frac_comp[
                        t, iz, j
                    ].value = pyo.value(
                        self.fuel_triple_phase_boundary.conc_mol_comp[t, iz, j] / denom
                    )
                    if j in self.fuel_triple_phase_boundary.reacting_gas_list:
                        self.fuel_triple_phase_boundary.log_mole_frac_comp[
                            t, iz, j
                        ].value = pyo.value(
                            pyo.log(
                                self.fuel_triple_phase_boundary.mole_frac_comp[t, iz, j]
                            )
                        )

        common._init_solve_block(self.fuel_triple_phase_boundary, solver_obj, solve_log)

        self.fuel_triple_phase_boundary.temperature_deviation_x.unfix()
        self.fuel_triple_phase_boundary.conc_mol_comp_deviation_x.unfix()
        self.fuel_triple_phase_boundary.conc_mol_comp_ref.unfix()
        self.fuel_triple_phase_boundary.heat_flux_x1.unfix()

        self.oxygen_triple_phase_boundary.temperature_deviation_x.fix(0)
        self.oxygen_triple_phase_boundary.conc_mol_comp_deviation_x.fix(0)
        self.oxygen_triple_phase_boundary.conc_mol_comp_ref.fix()
        self.oxygen_triple_phase_boundary.heat_flux_x0.fix(0)
        for t in self.flowsheet().time:
            for iz in self.oxygen_triple_phase_boundary.iznodes:
                denom = pyo.value(
                    sum(
                        self.oxygen_triple_phase_boundary.conc_mol_comp[t, iz, j]
                        for j in self.oxygen_triple_phase_boundary.component_list
                    )
                )
                for j in self.oxygen_triple_phase_boundary.component_list:
                    self.oxygen_triple_phase_boundary.mole_frac_comp[
                        t, iz, j
                    ].value = pyo.value(
                        self.oxygen_triple_phase_boundary.conc_mol_comp[t, iz, j]
                        / denom
                    )
                    if j in self.oxygen_triple_phase_boundary.reacting_gas_list:
                        self.oxygen_triple_phase_boundary.log_mole_frac_comp[
                            t, iz, j
                        ].value = pyo.value(
                            pyo.log(
                                self.oxygen_triple_phase_boundary.mole_frac_comp[
                                    t, iz, j
                                ]
                            )
                        )

        common._init_solve_block(
            self.oxygen_triple_phase_boundary, solver_obj, solve_log
        )

        self.oxygen_triple_phase_boundary.temperature_deviation_x.unfix()
        self.oxygen_triple_phase_boundary.conc_mol_comp_deviation_x.unfix()
        self.oxygen_triple_phase_boundary.conc_mol_comp_ref.unfix()
        self.oxygen_triple_phase_boundary.heat_flux_x0.unfix()

        init_log.info_high("Initializing Fuel Electrode")
        if self.config.thin_fuel_electrode:
            self.fuel_electrode.initialize_build(
                outlvl=outlvl,
                solver=solver,
                optarg=optarg,
                fix_heat_flux_x0=True,
            )
        else:
            self.fuel_electrode.conc_mol_comp_ref.fix()
            self.fuel_electrode.conc_mol_comp_deviation_x0.fix()
            self.fuel_electrode.temperature_deviation_x0.fix()
            self.fuel_electrode.heat_flux_x0.fix()
            self.fuel_electrode.material_flux_x1.fix()

            self.fuel_electrode.initialize_build(
                outlvl=outlvl,
                solver=solver,
                optarg=optarg,
                temperature_guess=temperature_guess,
            )

            self.fuel_electrode.conc_mol_comp_ref.unfix()
            self.fuel_electrode.conc_mol_comp_deviation_x0.unfix()
            self.fuel_electrode.temperature_deviation_x0.unfix()
            self.fuel_electrode.heat_flux_x0.unfix()
            self.fuel_electrode.material_flux_x1.unfix()

        init_log.info_high("Initializing Oxygen Electrode")
        if self.config.thin_oxygen_electrode:
            self.oxygen_electrode.initialize_build(
                outlvl=outlvl,
                solver=solver,
                optarg=optarg,
                fix_heat_flux_x0=False,
            )
        else:
            self.oxygen_electrode.conc_mol_comp_ref.fix()
            self.oxygen_electrode.conc_mol_comp_deviation_x1.fix()
            self.oxygen_electrode.temperature_deviation_x1.fix()
            self.oxygen_electrode.heat_flux_x1.fix()
            self.oxygen_electrode.material_flux_x0.fix()

            self.oxygen_electrode.initialize_build(
                outlvl=outlvl,
                solver=solver,
                optarg=optarg,
                temperature_guess=temperature_guess,
            )

            self.oxygen_electrode.conc_mol_comp_ref.unfix()
            self.oxygen_electrode.conc_mol_comp_deviation_x1.unfix()
            self.oxygen_electrode.temperature_deviation_x1.unfix()
            self.oxygen_electrode.heat_flux_x1.unfix()
            self.oxygen_electrode.material_flux_x0.unfix()

        init_log.info_high("Initializing Triple Phase Boundaries")
        self.fuel_triple_phase_boundary.initialize_build(
            outlvl=outlvl, solver=solver, optarg=optarg, fix_x0=True
        )
        self.oxygen_triple_phase_boundary.initialize_build(
            outlvl=outlvl, solver=solver, optarg=optarg, fix_x0=False
        )

        self.temperature_z.unfix()
        self.mean_temperature_eqn.activate()

        init_log.info_high("Solving cell with fixed current density")
        common._init_solve_block(self, solver_obj, solve_log)

        self.potential_eqn.activate()
        # TODO---does the user ever have a reason to fix current density
        # besides initialization and initial conditions?
        self.current_density.unfix()
        self.potential.fix()

        init_log.info_high("Solving cell with potential equations active")
        common._init_solve_block(self, solver_obj, solve_log)

        # Unfix any inlet variables that were fixed by initialization
        # To be honest, using a state block would probably have been less work
        for channel, comps in zip(
            ["fuel", "oxygen"], [self.fuel_component_list, self.oxygen_component_list]
        ):
            pname = channel + "_inlet"
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
        pass

    def model_check(self, steady_state=True):
        self.fuel_channel.model_check()
        self.oxygen_channel.model_check()
        self.fuel_electrode.model_check()
        self.oxygen_electrode.model_check()
        self.electrolyte.model_check()

        comp_set = set(self.fuel_component_list)
        comp_set = comp_set.union(self.oxygen_component_list)
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
                for chan, comps in zip(
                    [self.fuel_channel, self.oxygen_channel],
                    [self.fuel_component_list, self.oxygen_component_list],
                ):
                    sum_in += sum(
                        _element_dict[element][j] * chan.flow_mol_comp_inlet[t, j]
                        for j in comps
                    )
                    sum_out += sum(
                        _element_dict[element][j] * chan.flow_mol_comp_inlet[t, j]
                        for j in comps
                    )
                fraction_change = pyo.value((sum_out - sum_in) / sum_in)
                if abs(fraction_change) > 1e-5:
                    raise RuntimeError(
                        f"{element} is not being conserved {self.name}; "
                        f"fractional change {fraction_change}"
                    )

            enth_in = (
                self.fuel_channel.enth_mol_inlet[t]
                * self.fuel_channel.flow_mol_inlet[t]
                + self.oxygen_channel.enth_mol_inlet[t]
                * self.oxygen_channel.flow_mol_inlet[t]
            )
            enth_out = (
                self.fuel_channel.enth_mol_outlet[t]
                * self.fuel_channel.flow_mol_outlet[t]
                + self.oxygen_channel.enth_mol_outlet[t]
                * self.oxygen_channel.flow_mol_outlet[t]
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

        def ssf(c, s):
            iscale.set_scaling_factor(c, s, overwrite=False)

        sgsf = iscale.set_and_get_scaling_factor

        def cst(c, s):
            iscale.constraint_scaling_transform(c, s, overwrite=False)

        sdf = common._set_default_factor

        submodels = [
            self.fuel_channel,
            self.fuel_electrode,
            self.fuel_triple_phase_boundary,
            self.oxygen_channel,
            self.oxygen_electrode,
            self.oxygen_triple_phase_boundary,
            self.electrolyte,
        ]
        if self.config.include_contact_resistance:
            submodels.append(self.contact_interconnect_fuel_flow_mesh)
            submodels.append(self.contact_flow_mesh_fuel_electrode)
            submodels.append(self.contact_interconnect_oxygen_flow_mesh)
            submodels.append(self.contact_flow_mesh_oxygen_electrode)
        if self.config.flux_through_interconnect:
            submodels.append(self.interconnect)

        sy_def = 10
        s_inert_flux = 1e4

        sdf(self.current_density, 5e-4)
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

        # Need to scale material_flux_xes by component because inerts have much smaller
        # fluxes than actively reacting species.
        # TODO Revisit when reforming equations are added

        for t in self.flowsheet().time:
            sy_in_fuel = {}
            for j in self.fuel_component_list:
                sy_in_fuel[j] = gsf(self.fuel_inlet.mole_frac_comp[t, j], sy_def)
            sy_in_oxygen = {}
            for j in self.oxygen_component_list:
                sy_in_oxygen[j] = gsf(self.oxygen_inlet.mole_frac_comp[t, j], sy_def)
            for iz in self.iznodes:
                s_react_flux = gsf(self.current_density[t, iz]) * pyo.value(
                    Constants.faraday_constant
                )
                for j in self.fuel_component_list:
                    if j in self.config.inert_fuel_species_triple_phase_boundary:
                        s_flux_j = sy_in_fuel[j] * s_inert_flux
                    else:
                        s_flux_j = sy_in_fuel[j] * s_react_flux
                    ssf(self.fuel_channel.material_flux_x1[t, iz, j], s_flux_j)
                    if not self.config.thin_fuel_electrode:
                        ssf(self.fuel_electrode.material_flux_x1[t, iz, j], s_flux_j)

                    ssf(
                        self.fuel_triple_phase_boundary.mole_frac_comp[t, iz, j],
                        sy_in_fuel[j],
                    )
                for j in self.oxygen_component_list:
                    if j in self.config.inert_oxygen_species_triple_phase_boundary:
                        s_flux_j = sy_in_oxygen[j] * s_inert_flux
                    else:
                        s_flux_j = sy_in_oxygen[j] * s_react_flux
                    # s_flux_j = sy_in_oxygen[j]*s_mat_flux / max(abs(self.oxygen_triple_phase_boundary.tpb_stoich[j]),0.25)
                    ssf(self.oxygen_channel.material_flux_x0[t, iz, j], s_flux_j)
                    if not self.config.thin_oxygen_electrode:
                        ssf(self.oxygen_electrode.material_flux_x0[t, iz, j], s_flux_j)
                    ssf(
                        self.oxygen_triple_phase_boundary.mole_frac_comp[t, iz, j],
                        sy_in_oxygen[j],
                    )

                s_q_flux = s_react_flux * 1e-4  # Chosen heuristically based on TdS_rxn
                # s_q_flux = s_react_flux*1e-6
                for submodel in submodels:
                    ssf(submodel.heat_flux_x0[t, iz], s_q_flux)
                    ssf(submodel.heat_flux_x1[t, iz], s_q_flux)
                if not self.config.flux_through_interconnect:
                    if self.config.include_contact_resistance:
                        sq = gsf(
                            self.contact_interconnect_fuel_flow_mesh.heat_flux_x0[
                                t, iz
                            ],
                            default=s_q_flux,
                        )
                        cst(
                            self.no_heat_flux_fuel_interconnect_eqn[t, iz],
                            sq,
                        )
                        sq = gsf(
                            self.contact_interconnect_oxygen_flow_mesh.heat_flux_x1[
                                t, iz
                            ],
                            default=s_q_flux,
                        )
                        cst(
                            self.no_heat_flux_oxygen_interconnect_eqn[t, iz],
                            sq,
                        )
                    else:
                        sq = gsf(
                            self.fuel_channel.heat_flux_x0[t, iz],
                            default=s_q_flux,
                        )
                        cst(
                            self.no_heat_flux_fuel_interconnect_eqn[t, iz],
                            sq,
                        )
                        sq = gsf(
                            self.oxygen_channel.heat_flux_x1[t, iz],
                            default=s_q_flux,
                        )
                        cst(
                            self.no_heat_flux_oxygen_interconnect_eqn[t, iz],
                            sq,
                        )
        for idx, con in self.mean_temperature_eqn.items():
            cst(con, 1)

        for submodel in submodels:
            submodel.recursive_scaling()

        for t in self.flowsheet().time:
            for iz in self.iznodes:
                if self.config.thin_electrolyte:
                    sT = gsf(
                        self.electrolyte.temperature_deviation_x.referent[t, iz],
                        default=1,
                    )
                    cst(self.electrolyte_temperature_continuity_eqn[t, iz], sT)
                if (
                    self.config.flux_through_interconnect
                    and self.config.thin_interconnect
                ):
                    sT1 = gsf(
                        self.fuel_channel.temperature_deviation_x0[t, iz], default=1
                    )
                    sT2 = gsf(
                        self.oxygen_channel.temperature_deviation_x1[t, iz], default=1
                    )
                    sT = min(sT1, sT2)
                    cst(self.interconnect_temperature_continuity_eqn[t, iz], sT)
