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
IDAES 0D Fixed Bed Reactor model.
"""

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Var,
    value,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.dae import DerivativeVar
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    EnergyBalanceType,
    UnitModelBlockData,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.constants import Constants as constants
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

__author__ = "Chinedu Okoli, Andrew Lee"


# Assumptions:
# Only Vap and Sol phases, which are explicitly named as such
# Perfect mixing in Vap phase
# Static solid phase


@declare_process_block_class("FixedBed0D")
class FixedBed0DData(UnitModelBlockData):
    """
    0D Fixed Bed Reactor Model Class
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([True]),
            default=True,
            description="Dynamic model flag - must be True",
            doc="""Indicates whether this model will be dynamic or not,
**default** = True. Fixed beds must be dynamic.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=True,
            domain=In([True]),
            description="Holdup construction flag - must be True",
            doc="""Indicates whether holdup terms should be constructed or not.
**default** - True. Fixed bed reactors must be dynamic, thus this must be
True.""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.enthalpyTotal,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "gas_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for gas phase",
            doc="""Property parameter object used to define property calculations
for the gas phase, **default** - None.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "gas_property_package_args",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Arguments to use for constructing gas phase "
            "property packages",
            doc="""A ConfigBlock with arguments to be passed to a gas phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "solid_property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for solid phase",
            doc="""Property parameter object used to define property calculations
for the solid phase, **default** - None.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "solid_property_package_args",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Arguments to use for constructing solid phase "
            "property packages",
            doc="""A ConfigBlock with arguments to be passed to a solid phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "reaction_package",
        ConfigValue(
            default=None,
            # Removed domain argument due to limitations in ReactionBase for
            # heterogeneous systems
            description="Reaction package to use for unit",
            doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing reaction packages",
            doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
        ),
    )

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()

        # Get units meta data from property packages (only solid needed)
        units_meta_solid = (
            self.config.solid_property_package.get_metadata().get_derived_units
        )

        # Build Gas Phase StateBlock
        # This block only needed so gas conc. to rxn block is calculated
        # Defined state is set to True so that the "sum(mole_frac)=1" eqn in
        # the gas state block is deactivated.
        self.gas = self.config.gas_property_package.state_block_class(
            self.flowsheet().time,
            parameters=self.config.gas_property_package,
            defined_state=True,
        )

        # Build Solid Phase StateBlock
        # As there is no solid flow, only need a single state block
        # Defined state is set to True so that the "sum(mass_frac)=1" eqn in
        # the solid state block is deactivated. This is done here as there is
        # currently no way to deactivate the constraint at the initial time
        # for batch systems (i.e. no inlet or outlet ports).
        # The "sum(mass_frac)=1 for all t neq 0" eqn is written in the
        # unit model instead
        self.solids = self.config.solid_property_package.state_block_class(
            self.flowsheet().time,
            parameters=self.config.solid_property_package,
            defined_state=True,
        )

        tmp_dict = dict(**self.config.reaction_package_args)
        tmp_dict["gas_state_block"] = self.gas
        tmp_dict["solid_state_block"] = self.solids
        tmp_dict["parameters"] = self.config.reaction_package
        self.reactions = self.config.reaction_package.reaction_block_class(
            self.flowsheet().time,
            doc="Reaction properties in control volume",
            **tmp_dict,
        )

        # Volume of reactor
        self.bed_diameter = Var(
            initialize=1, doc="Reactor diameter", units=units_meta_solid("length")
        )
        self.bed_height = Var(
            initialize=1, doc="Bed length", units=units_meta_solid("length")
        )
        self.volume_bed = Var(
            initialize=1.0,
            doc="Volume of the reactor bed",
            units=units_meta_solid("volume"),
        )

        @self.Constraint(doc="Calculating volume of the reactor bed")
        def volume_bed_constraint(b):
            return b.volume_bed == (
                constants.pi * b.bed_height * (0.5 * b.bed_diameter) ** 2
            )

        # Solid phase material balance
        # Volume of solid
        self.volume_solid = Var(
            self.flowsheet().config.time,
            initialize=1.0,
            doc="Solids phase volume including particle pores",
            units=units_meta_solid("volume"),
        )

        @self.Constraint(
            self.flowsheet().config.time, doc="Calculating solid phase volume"
        )
        def volume_solid_constraint(b, t):
            return b.volume_solid[t] == (
                b.volume_bed * (1 - b.solids[t]._params.voidage)
            )

        # Accumulation equal to mass transfer/reaction
        self.solids_material_holdup = Var(
            self.flowsheet().time,
            self.config.solid_property_package.component_list,
            initialize=1,
            doc="Solid phase component holdups",
            units=units_meta_solid("mass"),
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.solid_property_package.component_list,
            doc="Solid phase material holdup constraints",
        )
        def solids_material_holdup_constraints(b, t, j):
            return b.solids_material_holdup[t, j] == (
                b.volume_solid[t] * b.solids[t].get_material_density_terms("Sol", j)
            )

        self.solids_material_accumulation = DerivativeVar(
            self.solids_material_holdup,
            wrt=self.flowsheet().config.time,
            doc="Solids material accumulation",
            units=units_meta_solid("mass"),
        )
        # should be mass/time, to get around DAE limitations var is
        # divided by time units when used

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.solid_property_package.component_list,
            doc="Solid phase material accumulation constraints",
        )
        def solids_material_accumulation_constraints(b, t, j):
            return b.solids_material_accumulation[t, j] / units_meta_solid(
                "time"
            ) == b.volume_solid[t] * b.solids[t]._params.mw_comp[j] * sum(
                b.reactions[t].reaction_rate[r]
                * b.config.reaction_package.rate_reaction_stoichiometry[r, "Sol", j]
                for r in b.config.reaction_package.rate_reaction_idx
            )

        # Add solid mass variable and constraint for TGA tracking
        self.mass_solids = Var(
            self.flowsheet().config.time,
            doc="Total mass of solids",
            units=units_meta_solid("mass"),
        )

        @self.Constraint(
            self.flowsheet().config.time, doc="Calculating total mass of solids"
        )
        def mass_solids_constraint(b, t):
            return b.mass_solids[t] == b.volume_solid[t] * sum(
                b.solids[t].get_material_density_terms("Sol", j)
                for j in b.config.solid_property_package.component_list
            )

        # Sum of mass fractions at all time equals 1
        @self.Constraint(
            self.flowsheet().config.time, doc="Sum of mass fractions at all time"
        )
        def sum_component_constraint(b, t):
            if t == b.flowsheet().config.time.first():
                return Constraint.Skip
            else:
                return 1 == sum(
                    b.solids[t].mass_frac_comp[j]
                    for j in b.solids[t]._params.component_list
                )

        if self.config.energy_balance_type != EnergyBalanceType.none:
            # Solid phase energy balance
            # Accumulation equal to heat transfer
            self.solids_energy_holdup = Var(
                self.flowsheet().time,
                initialize=1,
                doc="Solid phase energy holdup",
                units=units_meta_solid("energy"),
            )

            @self.Constraint(
                self.flowsheet().config.time,
                doc="Solid phase energy holdup constraints",
            )
            def solids_energy_holdup_constraints(b, t):
                return b.solids_energy_holdup[t] == (
                    self.volume_solid[t]
                    * pyunits.convert(
                        b.solids[t].get_energy_density_terms("Sol"),
                        to_units=units_meta_solid("energy")
                        / units_meta_solid("volume"),
                    )
                )

            self.solids_energy_accumulation = DerivativeVar(
                self.solids_energy_holdup,
                wrt=self.flowsheet().config.time,
                doc="Solids energy accumulation",
                units=units_meta_solid("energy"),
            )
            # should be energy/time, to get around DAE limitations var is
            # divided by time units when used

            @self.Constraint(
                self.flowsheet().config.time,
                doc="Solid phase energy accumulation constraints",
            )
            def solids_energy_accumulation_constraints(b, t):
                return b.solids_energy_accumulation[t] / units_meta_solid(
                    "time"
                ) == -sum(
                    b.reactions[t].reaction_rate[r]
                    * b.volume_solid[t]
                    * pyunits.convert(
                        b.reactions[t].dh_rxn[r],
                        to_units=units_meta_solid("energy_mole"),
                    )
                    for r in b.config.reaction_package.rate_reaction_idx
                )

        if self.config.energy_balance_type == EnergyBalanceType.none:
            # Fix solids temperature to initial value for isothermal conditions
            @self.Constraint(
                self.flowsheet().config.time, doc="Isothermal solid phase constraint"
            )
            def isothermal_solid_phase(b, t):
                if t == b.flowsheet().config.time.first():
                    return Constraint.Skip
                else:
                    return b.solids[t].temperature == b.solids[0].temperature

    def initialize_build(
        blk,
        gas_phase_state_args=None,
        solid_phase_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for FB0D unit.

        Keyword Arguments:
            gas_phase_state_args : a dict of arguments to be passed to the
                        property package(s) to provide an initial state for
                        initialization (see documentation of the specific
                        property package) (default = None).
            solid_phase_state_args : a dict of arguments to be passed to the
                        property package(s) to provide an initial state for
                        initialization (see documentation of the specific
                        property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """

        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        opt = get_solver(solver, optarg)  # create solver
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialize control volume block
        init_log.info("Initialize Thermophysical Properties")

        # Initialize gas_phase block
        gas_phase_flags = blk.gas.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=gas_phase_state_args
        )

        # Initialize solid_phase properties block
        solid_phase_flags = blk.solids.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=solid_phase_state_args,
        )

        init_log.info("Initialize reaction properties")
        # Initialize reactions
        blk.reactions.initialize(outlvl=outlvl, optarg=optarg, solver=solver)

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, outlvl) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.gas.release_state(gas_phase_flags, outlvl)
        blk.solids.release_state(solid_phase_flags, outlvl)

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # Get units meta data from property packages (only solid needed)
        units_meta_solid = (
            self.config.solid_property_package.get_metadata().get_derived_units
        )

        # scale some variables
        if hasattr(self, "bed_diameter"):
            if iscale.get_scaling_factor(self.bed_diameter) is None:
                sf = 1 / value(self.bed_diameter)
                iscale.set_scaling_factor(self.bed_diameter, 1e-1 * sf)

        if hasattr(self, "bed_height"):
            if iscale.get_scaling_factor(self.bed_height) is None:
                sf = 1 / value(self.bed_height)
                iscale.set_scaling_factor(self.bed_height, 1e-1 * sf)

        if hasattr(self, "volume_bed"):
            if iscale.get_scaling_factor(self.volume_bed) is None:
                sf = 1 / value(
                    constants.pi * self.bed_height * (0.5 * self.bed_diameter) ** 2
                )
                iscale.set_scaling_factor(self.volume_bed, 1e-1 * sf)

        # Scale some constraints
        if hasattr(self, "volume_bed_constraint"):
            for c in self.volume_bed_constraint.values():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.volume_bed), overwrite=False
                )

        if hasattr(self, "volume_solid_constraint"):
            for t, c in self.volume_solid_constraint.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.volume_bed), overwrite=False
                )

        if hasattr(self, "solids_material_holdup_constraints"):
            for (t, j), c in self.solids_material_holdup_constraints.items():
                sf1 = iscale.get_scaling_factor(self.volume_bed)
                sf2 = iscale.get_scaling_factor(self.solids[t].dens_mass_skeletal)
                iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

        if hasattr(self, "solids_material_accumulation_constraints"):
            for (t, j), c in self.solids_material_accumulation_constraints.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.volume_bed), overwrite=False
                )

        if hasattr(self, "mass_solids_constraint"):
            for t, c in self.mass_solids_constraint.items():
                sf1 = iscale.get_scaling_factor(self.volume_bed)
                sf2 = iscale.get_scaling_factor(self.solids[t].dens_mass_skeletal)
                iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

        if hasattr(self, "sum_component_constraint"):
            # Get a single representative value in component list
            for j in self.config.solid_property_package.component_list:
                break
            for t, c in self.sum_component_constraint.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.solids[t].mass_frac_comp[j]),
                    overwrite=False,
                )

        if self.config.energy_balance_type != EnergyBalanceType.none:

            if hasattr(self, "solids_energy_holdup_constraints"):
                for t, c in self.solids_energy_holdup_constraints.items():
                    sf1 = iscale.get_scaling_factor(self.volume_bed)
                    sf2 = iscale.get_scaling_factor(self.solids[t].enth_mass)
                    iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

            if hasattr(self, "solids_energy_accumulation_constraints"):
                for t, c in self.solids_energy_accumulation_constraints.items():
                    sf1 = iscale.get_scaling_factor(self.volume_bed)
                    sf2 = iscale.get_scaling_factor(self.solids[t].enth_mass)
                    iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

        if self.config.energy_balance_type == EnergyBalanceType.none:

            if hasattr(self, "isothermal_solid_phase"):
                for t, c in self.isothermal_solid_phase.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(self.solids[t].temperature),
                        overwrite=False,
                    )
