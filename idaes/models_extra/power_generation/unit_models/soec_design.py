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

__author__ = "John Eslick"

import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

import idaes.core.util.constants as iconst
from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core.util.config import is_physical_parameter_block
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)
import idaes.models.unit_models as um  # um = unit models
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    get_rxn,
    EosType,
)
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog
from idaes.core.util import from_json, to_json, StoreSpec


@declare_process_block_class("SoecDesign", doc="Simple SOEC model for process design.")
class SoecDesignData(UnitModelBlockData):
    """Simple 0D SOEC model. This model can be used to develop a flowsheet at the
    design point or for the basis of a surrogate. For off-design applications,
    a more detailed model is required."""

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    CONFIG.declare(
        "oxygen_side_property_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the oxygen side.",
            doc=(
                "Property package for the oxygen side, using "
                "idaes.models_extra.power_generation.properties.natural_gas_PR is "
                "strongly recomended, either Peng-Robinson or Ideal is okay"
            ),
        ),
    )
    CONFIG.declare(
        "hydrogen_side_property_package",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="Property package for the hydrogen side.",
            doc=(
                "Property package for the hydrogen side, using "
                "idaes.models_extra.power_generation.properties.natural_gas_PR is "
                "strongly recomended, either Peng-Robinson or Ideal is okay"
            ),
        ),
    )
    CONFIG.declare(
        "oxygen_side_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Property package arguments for the oxygen side.",
            doc="Property package arguments for the oxygen side.",
        ),
    )
    CONFIG.declare(
        "hydrogen_side_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Property package arguments for the hydrogen side.",
            doc="Property package arguments for the hydrogen side.",
        ),
    )
    CONFIG.declare(
        "reaction_eos",
        ConfigValue(
            default=EosType.PR,
            domain=In(EosType),
            description="Physical properies for electrolysis reactions",
            doc=(
                "Reaction properties equation of state in: "
                "{EosType.PR, EosType.IDEAL}."
            ),
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat transfer term construction flag",
            doc="Indicates whether the SOEC is adiabatic. Default=False (adiabatic)",
        ),
    )

    def build(self):
        """Construct a unit model data block"""
        super().build()

        # Do a little validation.
        # The hydrogen side needs water and hydrogen
        if set(self.config.hydrogen_side_property_package.component_list) != {
            "H2",
            "H2O",
        }:
            raise ConfigurationError(
                "SOEC hydrogen side must contain exactly H2 and H2O"
            )
        # The oxygen side must have oxygen in it.
        if not "O2" in self.config.oxygen_side_property_package.component_list:
            raise ConfigurationError("SOEC oxygen side must contain O2")

        # build the block, use smaller methods to make it easier to follow
        self._add_electrolysis_properties()
        self._add_unit_models()
        self._add_arcs()
        self._add_variables()
        self._add_constraints()
        self._add_o2_translator_constraints()
        self._add_h2_inlet_translator_constraints()
        self._add_h2_outlet_translator_constraints()
        self._add_ports()
        self._scaling_guess()

    def _add_electrolysis_properties(self):
        """Add the electrolysis property package.  This has H2, H2O and O2,
        and is only used internally by the unit for the electrolysis reaction.
        """
        self.electrolysis_prop_params = GenericParameterBlock(
            **get_prop(
                {"H2O", "H2", "O2"}, phases={"Vap"}, eos=self.config.reaction_eos
            ),
            doc="Physical property parameters for the electrolysis reaction",
        )
        # Electrolysis is just the hydrogen combustion reaction backward
        self.electrolysis_rxn_params = GenericReactionParameterBlock(
            **get_rxn(self.electrolysis_prop_params, {"h2_cmb"}),
            doc="Reaction parameters",
        )

    def _add_unit_models(self):
        """Add the unit models that make up the composite model. The main
        things happening are 1) Electrolysis Reaction, 2) O2 Seperation, 3) O2
        mixes with sweep gas, and 4) sweep gas heat transfer.  There are
        translator blocks to switch between property pacakage depending on the
        existance of components in various streams.
        """
        self.h2_inlet_translator = um.Translator(
            doc="Translate hydrogen side properties to electrolysis properties",
            inlet_property_package=self.config.hydrogen_side_property_package,
            inlet_property_package_args=self.config.hydrogen_side_property_package_args,
            outlet_property_package=self.electrolysis_prop_params,
            outlet_state_defined=False,
        )
        self.electrolysis_reactor = um.StoichiometricReactor(
            doc="Electrolysis reator",
            property_package=self.electrolysis_prop_params,
            reaction_package=self.electrolysis_rxn_params,
            has_pressure_change=False,
            has_heat_transfer=True,
        )
        self.o2_seperator = um.Separator(
            property_package=self.electrolysis_prop_params,
            split_basis=um.SplittingType.componentFlow,
            outlet_list=["h2_strm", "o2_strm"],
        )
        self.h2_outlet_translator = um.Translator(
            doc="Translate electrolysis properties to hydrogen properies",
            inlet_property_package=self.electrolysis_prop_params,
            outlet_property_package=self.config.hydrogen_side_property_package,
            outlet_property_package_args=self.config.hydrogen_side_property_package_args,
            outlet_state_defined=False,
        )
        self.o2_translator = um.Translator(
            doc="Translate electrolysis properties to oxygen properies",
            inlet_property_package=self.electrolysis_prop_params,
            outlet_property_package=self.config.oxygen_side_property_package,
            outlet_property_package_args=self.config.oxygen_side_property_package_args,
            outlet_state_defined=False,
        )
        self.o2_mixer = um.Mixer(
            property_package=self.config.oxygen_side_property_package,
            property_package_args=self.config.oxygen_side_property_package_args,
            momentum_mixing_type=um.MomentumMixingType.none,
            inlet_list=["sweep_strm", "o2_strm"],
        )
        self.sweep_heater = um.Heater(
            property_package=self.config.oxygen_side_property_package,
            property_package_args=self.config.oxygen_side_property_package_args,
        )

    def _add_arcs(self):
        """Add streams to connect internal units."""
        self.strm1 = Arc(
            doc="Hydrogen side inlet to electrolysis translator",
            source=self.h2_inlet_translator.outlet,
            destination=self.electrolysis_reactor.inlet,
        )
        self.strm2 = Arc(
            doc="Electrolysis reactor to oxygen seperation",
            source=self.electrolysis_reactor.outlet,
            destination=self.o2_seperator.inlet,
        )
        self.strm3 = Arc(
            doc="Oxygen seperation to hydrogen side outlet translator",
            source=self.o2_seperator.h2_strm,
            destination=self.h2_outlet_translator.inlet,
        )
        self.strm4 = Arc(
            doc="Oxygen from electrolysis to translator to sweep properies",
            source=self.o2_seperator.o2_strm,
            destination=self.o2_translator.inlet,
        )
        self.strm5 = Arc(
            doc="Oxygen translator to O2/sweep mixer",
            source=self.o2_translator.outlet,
            destination=self.o2_mixer.o2_strm,
        )
        self.strm6 = Arc(
            doc="O2/sweep mixer to sweep heater",
            source=self.o2_mixer.outlet,
            destination=self.sweep_heater.inlet,
        )
        # Probably don't need this and may need to call it again for dynamic
        # models, but this will make it easy for a user to have a flowsheet
        # with one steady state model in it.
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _add_variables(self):
        """Add variables for unit model level constraints"""
        self.current = pyo.Var(
            self.flowsheet().time,
            initialize=0,
            units=pyo.units.ampere,
            doc="SOEC module electric currrent",
        )
        self.water_utilization = pyo.Var(
            self.flowsheet().time,
            units=pyo.units.dimensionless,
            doc="Fraction of water used in electrolysis reaction",
        )
        self.cell_potential = pyo.Var(
            self.flowsheet().time,
            initialize=1.3,
            units=pyo.units.volts,
            doc="Electric potential for a single cell",
        )
        self.heat = pyo.Var(
            self.flowsheet().time,
            units=pyo.units.W,
            initialize=0,
            doc="Heat transfered to the module",
        )
        self.hydrogen_side_outlet_temperature = pyo.Reference(
            self.electrolysis_reactor.control_volume.properties_out[:].temperature
        )
        self.oxygen_side_outlet_temperature = pyo.Reference(
            self.sweep_heater.control_volume.properties_out[:].temperature
        )

    def _add_constraints(self):
        """Add unit model constraints"""

        @self.Constraint(self.flowsheet().time)
        def water_utilization_eqn(b, t):
            return b.electrolysis_reactor.control_volume.properties_out[
                t
            ].flow_mol_comp[
                "H2O"
            ] == b.electrolysis_reactor.control_volume.properties_in[
                t
            ].flow_mol_comp[
                "H2O"
            ] * (
                1.0 - b.water_utilization[t]
            )

        @self.Expression(self.flowsheet().time)
        def current_expr(b, t):
            return (
                (
                    (
                        b.electrolysis_reactor.control_volume.properties_in[
                            t
                        ].mole_frac_comp["H2O"]
                        * b.electrolysis_reactor.control_volume.properties_in[
                            t
                        ].flow_mol
                    )
                    - (
                        b.electrolysis_reactor.control_volume.properties_out[
                            t
                        ].mole_frac_comp["H2O"]
                        * b.electrolysis_reactor.control_volume.properties_out[
                            t
                        ].flow_mol
                    )
                )
                * 2
                * iconst.Constants.faraday_constant
            )

        @self.Constraint(self.flowsheet().time)
        def current_eqn(b, t):
            return b.current[t] == b.current_expr[t]

        # Fix the O2 seperator split fractions to just remove O2 from fuel side
        self.o2_seperator.split_fraction[:, "o2_strm", "O2"].fix(1)
        self.o2_seperator.split_fraction[:, "o2_strm", "H2O"].fix(0)
        self.o2_seperator.split_fraction[:, "o2_strm", "H2"].fix(0)

        @self.o2_mixer.Constraint(self.flowsheet().time)
        def pressure_eqn(b, t):
            return b.mixed_state[t].pressure == b.sweep_strm_state[t].pressure

        if not self.config.has_heat_transfer:

            @self.Constraint(self.flowsheet().time)
            def heat_transfer_eqn(b, t):
                return b.heat[t] == 0

        @self.Expression(self.flowsheet().time)
        def delta_enth(b, t):
            return (
                b.o2_seperator.h2_strm_state[t].get_enthalpy_flow_terms("Vap")
                + b.sweep_heater.control_volume.properties_out[
                    t
                ].get_enthalpy_flow_terms("Vap")
                - b.electrolysis_reactor.control_volume.properties_in[
                    t
                ].get_enthalpy_flow_terms("Vap")
                - b.o2_mixer.sweep_strm_state[t].get_enthalpy_flow_terms("Vap")
                + b.heat[t]
            )

        @self.Constraint(self.flowsheet().time)
        def cell_potential_eqn(b, t):
            return b.cell_potential[t] == b.delta_enth[t] / b.current[t]

    @staticmethod
    def _translator_ftp_constraints(translator):
        """Generalize some translator block constraints"""

        @translator.Constraint(translator.flowsheet().time)
        def temperature_eqn(b, t):
            return b.properties_out[t].temperature == b.properties_in[t].temperature

        @translator.Constraint(translator.flowsheet().time)
        def flow_mol_eqn(b, t):
            return b.properties_out[t].flow_mol == b.properties_in[t].flow_mol

        @translator.Constraint(translator.flowsheet().time)
        def pressure_eqn(b, t):
            return b.properties_out[t].pressure == b.properties_in[t].pressure

    def _add_o2_translator_constraints(self):
        """Add constraints to traslate from electrolysis properties to sweep
        properties.  The stream being translated contains only O2.
        """
        t0 = self.flowsheet().time.first()
        comps = set(self.o2_translator.properties_out[t0].mole_frac_comp.keys())
        comps.remove("O2")

        @self.o2_translator.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return b.properties_out[t].mole_frac_comp[c] == 1e-19

        self.o2_translator.properties_out[:].mole_frac_comp["O2"] = 1

        self._translator_ftp_constraints(self.o2_translator)

    def _add_h2_inlet_translator_constraints(self):
        """Translate the inlet hydrogen stream to electrolysis properties. This
        stream contains only H2 and H2O.
        """
        t0 = self.flowsheet().time.first()
        # There is no O2.  So set the value to zero.  Don't need a constraint
        # for this due to the sum = 1 constraint
        self.o2_translator.properties_out[:].mole_frac_comp["O2"] = 1e-19

        # Add constraints for the other componets
        comps = set(self.h2_inlet_translator.properties_in[t0].mole_frac_comp.keys())

        @self.h2_inlet_translator.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return (
                b.properties_out[t].mole_frac_comp[c]
                == b.properties_in[t].mole_frac_comp[c]
            )

        self._translator_ftp_constraints(self.h2_inlet_translator)

    def _add_h2_outlet_translator_constraints(self):
        """Translate the electrolysis properties to hydrogen side properties.
        This stream contains only H2 and H2O.
        """
        t0 = self.flowsheet().time.first()
        # Add constraints for the other componets
        comps = set(self.h2_outlet_translator.properties_out[t0].mole_frac_comp.keys())
        comps.remove("H2")  # Need to remove one component due to sum = 1 constraint

        @self.h2_outlet_translator.Constraint(self.flowsheet().time, comps)
        def mole_frac_comp_eqn(b, t, c):
            return (
                b.properties_out[t].mole_frac_comp[c]
                == b.properties_in[t].mole_frac_comp[c]
            )

        self._translator_ftp_constraints(self.h2_outlet_translator)

    def _add_ports(self):
        """Add unit level ports"""
        self.add_inlet_port(
            name="hydrogen_side_inlet",
            block=self.h2_inlet_translator.properties_in,
            doc="Hydrogen side inlet port",
        )
        self.add_inlet_port(
            name="oxygen_side_inlet",
            block=self.o2_mixer.sweep_strm_state,
            doc="Oxygen side inlet port",
        )
        self.add_outlet_port(
            name="hydrogen_side_outlet",
            block=self.h2_outlet_translator.properties_out,
            doc="Hydrogen side outlet port",
        )
        self.add_outlet_port(
            name="oxygen_side_outlet",
            block=self.sweep_heater.control_volume.properties_out,
            doc="Oxygen side outlet port",
        )

    def set_flow_scale(self, scale=1):
        """Set default flow scaling in the property packages based on the
        expected magnitude of typical flows

        Args:
            scale (float): Expected molar flow variable scale

        Returns:
            None
        """
        self.config.oxygen_side_property_package.set_default_scaling("flow_mol", scale)
        self.config.hydrogen_side_property_package.set_default_scaling(
            "flow_mol", scale
        )
        self.electrolysis_prop_params.set_default_scaling("flow_mol", scale)

        self.config.oxygen_side_property_package.set_default_scaling(
            "flow_mol_phase", scale
        )
        self.config.hydrogen_side_property_package.set_default_scaling(
            "flow_mol_phase", scale
        )
        self.electrolysis_prop_params.set_default_scaling("flow_mol_phase", scale)

        for (
            t,
            i,
        ), v in self.electrolysis_reactor.control_volume.rate_reaction_extent.items():
            iscale.set_scaling_factor(v, scale)
        for (
            t,
            p,
            i,
        ), v in (
            self.electrolysis_reactor.control_volume.rate_reaction_generation.items()
        ):
            iscale.set_scaling_factor(v, scale)

    def set_heat_scale(self, scale=1e-5):
        """Set the heat transfer scale factor roughly based on the size of the
        process.

        Args:
            scale (float): Expected heat transfer scale

        Returns:
            None
        """
        iscale.set_scaling_factor(self.electrolysis_reactor.control_volume.heat, scale)
        iscale.set_scaling_factor(self.sweep_heater.control_volume.heat, scale)
        iscale.set_scaling_factor(self.heat, scale)

    def set_current_scale(self, scale=1e-6):
        """Set the expected electrical current scale

        Args:
            scale (float): Expected electrical current scale

        Returns:
            None
        """
        iscale.set_scaling_factor(self.current, scale)

    def _scaling_guess(self):
        self.set_flow_scale()
        self.set_heat_scale()
        self.set_current_scale()

        # For components that exist we don't expect particularly low
        # concentrations, so change the default to 10, which should be good
        # for mole fractions > 0.01 which is what we expect here.
        self.config.oxygen_side_property_package.set_default_scaling(
            "mole_frac_comp", 10
        )
        self.config.hydrogen_side_property_package.set_default_scaling(
            "mole_frac_comp", 10
        )
        self.electrolysis_prop_params.set_default_scaling("mole_frac_comp", 10)
        self.config.oxygen_side_property_package.set_default_scaling(
            "mole_frac_phase_comp", 10
        )
        self.config.hydrogen_side_property_package.set_default_scaling(
            "mole_frac_phase_comp", 10
        )
        self.electrolysis_prop_params.set_default_scaling("mole_frac_phase_comp", 10)

        # Set some other scale factors that we have a good guess for
        unt = self.electrolysis_reactor
        iscale.set_scaling_factor(
            unt.control_volume.properties_in[0.0].enth_mol_phase["Vap"], 1e-4
        )
        iscale.set_scaling_factor(
            unt.control_volume.properties_out[0.0].enth_mol_phase["Vap"], 1e-4
        )
        unt = self.o2_seperator
        iscale.set_scaling_factor(unt.h2_strm_state[0.0].enth_mol_phase["Vap"], 1e-4)
        unt = self.o2_mixer
        iscale.set_scaling_factor(unt.o2_strm_state[0.0].enth_mol_phase["Vap"], 1e-4)
        iscale.set_scaling_factor(unt.o2_strm_state[0.0].enth_mol_phase["Vap"], 1e-4)
        iscale.set_scaling_factor(unt.mixed_state[0.0].enth_mol_phase["Vap"], 1e-4)
        iscale.set_scaling_factor(unt.sweep_strm_state[0.0].enth_mol_phase["Vap"], 1e-4)
        unt = self.sweep_heater
        iscale.set_scaling_factor(
            unt.control_volume.properties_in[0.0].enth_mol_phase["Vap"], 1e-4
        )
        iscale.set_scaling_factor(
            unt.control_volume.properties_out[0.0].enth_mol_phase["Vap"], 1e-4
        )

    def calculate_scaling_factors(self):
        """Calculate scale factors for the unit model equations"""
        for t in self.flowsheet().time:
            iscale.constraint_scaling_transform(
                self.water_utilization_eqn[t],
                iscale.get_scaling_factor(
                    self.electrolysis_reactor.control_volume.properties_in[t].flow_mol
                ),
            )
            iscale.constraint_scaling_transform(
                self.current_eqn[t], iscale.get_scaling_factor(self.current[t])
            )
            iscale.constraint_scaling_transform(
                self.o2_mixer.pressure_eqn[t],
                iscale.get_scaling_factor(self.o2_mixer.mixed_state[t].pressure),
            )
            iscale.constraint_scaling_transform(
                self.o2_translator.pressure_eqn[t],
                iscale.get_scaling_factor(self.o2_translator.properties_in[t].pressure),
            )
            iscale.constraint_scaling_transform(
                self.o2_translator.temperature_eqn[t],
                iscale.get_scaling_factor(
                    self.o2_translator.properties_in[t].temperature
                ),
            )
            iscale.constraint_scaling_transform(
                self.o2_translator.temperature_eqn[t],
                iscale.get_scaling_factor(
                    self.o2_translator.properties_in[t].temperature
                ),
            )
            iscale.constraint_scaling_transform(
                self.o2_translator.flow_mol_eqn[t],
                iscale.get_scaling_factor(self.o2_translator.properties_in[t].flow_mol),
            )
            iscale.constraint_scaling_transform(
                self.h2_inlet_translator.pressure_eqn[t],
                iscale.get_scaling_factor(
                    self.h2_inlet_translator.properties_in[t].pressure
                ),
            )
            iscale.constraint_scaling_transform(
                self.h2_inlet_translator.temperature_eqn[t],
                iscale.get_scaling_factor(
                    self.h2_inlet_translator.properties_in[t].temperature
                ),
            )
            iscale.constraint_scaling_transform(
                self.h2_inlet_translator.flow_mol_eqn[t],
                iscale.get_scaling_factor(
                    self.h2_inlet_translator.properties_in[t].flow_mol
                ),
            )
            iscale.constraint_scaling_transform(
                self.h2_outlet_translator.pressure_eqn[t],
                iscale.get_scaling_factor(
                    self.h2_outlet_translator.properties_in[t].pressure
                ),
            )
            iscale.constraint_scaling_transform(
                self.h2_outlet_translator.temperature_eqn[t],
                iscale.get_scaling_factor(
                    self.h2_outlet_translator.properties_in[t].temperature
                ),
            )
            iscale.constraint_scaling_transform(
                self.h2_outlet_translator.flow_mol_eqn[t],
                iscale.get_scaling_factor(
                    self.h2_outlet_translator.properties_in[t].flow_mol
                ),
            )
            iscale.constraint_scaling_transform(self.cell_potential_eqn[t], 1)
            if not self.config.has_heat_transfer:
                iscale.constraint_scaling_transform(
                    self.heat_transfer_eqn[t], iscale.get_scaling_factor(self.heat[t])
                )
        for (t, i), c in self.o2_translator.mole_frac_comp_eqn.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.o2_translator.properties_out[t].mole_frac_comp[i]
                ),
            )
        for (t, i), c in self.h2_inlet_translator.mole_frac_comp_eqn.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.h2_inlet_translator.properties_out[t].mole_frac_comp[i]
                ),
            )
        for (t, i), c in self.h2_outlet_translator.mole_frac_comp_eqn.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.h2_outlet_translator.properties_out[t].mole_frac_comp[i]
                ),
            )

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        init_log.info_high("SOEC Initialization Starting")
        solver_obj = get_solver(solver, optarg)

        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        self.oxygen_side_inlet.fix()
        self.hydrogen_side_inlet.fix()
        self.oxygen_side_outlet.unfix()
        self.hydrogen_side_outlet.unfix()

        self.h2_inlet_translator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm1)
        self.electrolysis_reactor.initialize(
            outlvl=outlvl, solver=solver, optarg=optarg
        )
        propagate_state(self.strm2)
        self.o2_seperator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm3)
        self.h2_outlet_translator.initialize(
            outlvl=outlvl, solver=solver, optarg=optarg
        )
        propagate_state(self.strm4)
        self.o2_translator.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm5)
        self.o2_mixer.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        propagate_state(self.strm6)
        self.sweep_heater.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        for t in self.current.keys():
            self.current[t] = pyo.value(self.current_expr[t])

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver_obj.solve(self, tee=slc.tee)
        if not pyo.check_optimal_termination(res):
            raise InitializationError(f"SOEC failed to initialize.")

        from_json(self, sd=istate, wts=sp)
        init_log.info_high("SOEC Initialization Complete")
