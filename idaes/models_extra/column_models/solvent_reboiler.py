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
Reboiler model for solvent columns.

This is a simple model for a reboiler in the case where liquid and vapor phases
have separate proeprty packages, such as the case of solvent columns.

Assumptions:
     * Steady-state only
     * Liquid phase property package has a single phase named Liq
     * Vapor phase property package has a single phase named Vap
     * Liquid and vapor phase properties need not have the same component lists
"""

__author__ = "Andrew Lee, Paul Akula"

# Import Pyomo libraries
from pyomo.environ import (
    check_optimal_termination,
    Constraint,
    Param,
    Reference,
    units as pyunits,
    value,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    MaterialFlowBasis,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError, InitializationError


_log = idaeslog.getIdaesLogger(__name__)


@declare_process_block_class("SolventReboiler")
class SolventReboilerData(UnitModelBlockData):
    """
    Reboiler unit for solvent column models using separate property packages
    for liquid and vpor phases.

    Unit model to reboil the liquid from the bottom of a solvent column.
    """

    CONFIG = ConfigBlock()
    # TOOO: Add dynamics in future
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Equilibrium Reactors do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Equilibrium reactors do not have defined volume, thus
this must be False.""",
        ),
    )
    # TODO : Add boilup ratio back later if needed
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
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
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "liquid_property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for liquid phase",
            doc="""Property parameter object used to define property calculations
for the liquid phase,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "liquid_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing liquid phase properties",
            doc="""A ConfigBlock with arguments to be passed to liquid phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "vapor_property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for vapor phase",
            doc="""Property parameter object used to define property calculations
for the vapor phase,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "vapor_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing vapor phase properties",
            doc="""A ConfigBlock with arguments to be passed to vapor phase
property block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """Build the model.

        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        # Check phase lists match assumptions
        if self.config.vapor_property_package.phase_list != ["Vap"]:
            raise ConfigurationError(
                f"{self.name} SolventReboiler model requires that the vapor "
                f"phase property package have a single phase named 'Vap'"
            )
        if self.config.liquid_property_package.phase_list != ["Liq"]:
            raise ConfigurationError(
                f"{self.name} SolventReboiler model requires that the liquid "
                f"phase property package have a single phase named 'Liq'"
            )

        # Check for at least one common component in component lists
        if not any(
            j in self.config.vapor_property_package.component_list
            for j in self.config.liquid_property_package.component_list
        ):
            raise ConfigurationError(
                f"{self.name} SolventReboiler model requires that the liquid "
                f"and vapor phase property packages have at least one "
                f"common component."
            )

        # ---------------------------------------------------------------------
        # Add Control Volume for the Liquid Phase
        self.liquid_phase = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.liquid_property_package,
            property_package_args=self.config.liquid_property_package_args,
        )

        self.liquid_phase.add_state_blocks(has_phase_equilibrium=True)

        # Separate liquid and vapor phases means that phase equilibrium will
        # be handled at the unit model level, thus has_phase_equilibrium is
        # False, but has_mass_transfer is True.
        self.liquid_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True,
            has_phase_equilibrium=False,
        )

        # Need to include enthalpy transfer term for the mass transfer
        self.liquid_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=True,
            has_enthalpy_transfer=True,
        )

        self.liquid_phase.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # ---------------------------------------------------------------------
        # Add single state block for vapor phase
        tmp_dict = dict(**self.config.vapor_property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = False
        self.vapor_phase = self.config.vapor_property_package.build_state_block(
            self.flowsheet().time, doc="Vapor phase properties", **tmp_dict
        )

        # ---------------------------------------------------------------------
        # Check flow basis is compatable
        # TODO : Could add code to convert flow bases, but not now
        t_init = self.flowsheet().time.first()
        if (
            self.vapor_phase[t_init].get_material_flow_basis()
            != self.liquid_phase.properties_out[t_init].get_material_flow_basis()
        ):
            raise ConfigurationError(
                f"{self.name} vapor and liquid property packages must use the "
                f"same material flow basis."
            )

        # ---------------------------------------------------------------------
        # Add Ports for the reboiler
        self.add_inlet_port(name="inlet", block=self.liquid_phase, doc="Liquid feed")
        self.add_outlet_port(
            name="bottoms", block=self.liquid_phase, doc="Bottoms stream"
        )
        self.add_outlet_port(
            name="vapor_reboil",
            block=self.vapor_phase,
            doc="Vapor stream from reboiler",
        )

        # ---------------------------------------------------------------------
        # Add unit level constraints
        # First, need the union and intersection of component lists
        all_comps = (
            self.vapor_phase.component_list
            | self.liquid_phase.properties_out.component_list
        )
        common_comps = (
            self.vapor_phase.component_list
            & self.liquid_phase.properties_out.component_list
        )

        # Get units for unit conversion
        vunits = self.config.vapor_property_package.get_metadata().get_derived_units
        lunits = self.config.liquid_property_package.get_metadata().get_derived_units
        flow_basis = self.vapor_phase[t_init].get_material_flow_basis()
        if flow_basis == MaterialFlowBasis.molar:
            fb = "flow_mole"
        elif flow_basis == MaterialFlowBasis.mass:
            fb = "flow_mass"
        else:
            raise ConfigurationError(
                f"{self.name} SolventReboiler only supports mass or molar "
                f"basis for MaterialFlowBasis."
            )

        if any(j not in common_comps for j in self.vapor_phase.component_list):
            # We have non-condensable components present, need zero-flow param
            self.zero_flow_param = Param(
                mutable=True, default=1e-8, units=vunits("flow_mole")
            )

        # Material balances
        def rule_material_balance(blk, t, j):
            if j in common_comps:
                # Component is in equilibrium
                # Mass transfer equals vapor flowrate
                return -blk.liquid_phase.mass_transfer_term[
                    t, "Liq", j
                ] == pyunits.convert(
                    blk.vapor_phase[t].get_material_flow_terms("Vap", j),
                    to_units=lunits(fb),
                )
            elif j in self.vapor_phase.component_list:
                # Non-condensable component
                # No mass transfer term
                # Set vapor flowrate to an arbitary small value
                return (
                    blk.vapor_phase[t].get_material_flow_terms("Vap", j)
                    == blk.zero_flow_param
                )
            else:
                # Non-vaporisable comonent
                # Mass transfer term is zero, no vapor flowrate
                return blk.liquid_phase.mass_transfer_term[t, "Liq", j] == 0 * lunits(
                    fb
                )

        self.unit_material_balance = Constraint(
            self.flowsheet().time,
            all_comps,
            rule=rule_material_balance,
            doc="Unit level material balances",
        )

        # Phase equilibrium constraints
        # For all common components, equate fugacity in vapor and liquid
        def rule_phase_equilibrium(blk, t, j):
            return blk.liquid_phase.properties_out[t].fug_phase_comp[
                "Liq", j
            ] == pyunits.convert(
                blk.vapor_phase[t].fug_phase_comp["Vap", j], to_units=lunits("pressure")
            )

        self.unit_phase_equilibrium = Constraint(
            self.flowsheet().time,
            common_comps,
            rule=rule_phase_equilibrium,
            doc="Unit level phase equilibrium constraints",
        )

        # Temperature equality constraint
        def rule_temperature_balance(blk, t):
            return blk.liquid_phase.properties_out[t].temperature == pyunits.convert(
                blk.vapor_phase[t].temperature, to_units=lunits("temperature")
            )

        self.unit_temperature_equality = Constraint(
            self.flowsheet().time,
            rule=rule_temperature_balance,
            doc="Unit level temperature equality",
        )

        # Unit level energy balance
        # Energy leaving in vapor phase must be equal and opposite to enthalpy
        # transfer from liquid phase
        def rule_energy_balance(blk, t):
            return -blk.liquid_phase.enthalpy_transfer[t] == pyunits.convert(
                blk.vapor_phase[t].get_enthalpy_flow_terms("Vap"),
                to_units=lunits("energy") / lunits("time"),
            )

        self.unit_enthalpy_balance = Constraint(
            self.flowsheet().time,
            rule=rule_energy_balance,
            doc="Unit level enthalpy_balance",
        )

        # Pressure balance constraint
        def rule_pressure_balance(blk, t):
            return blk.liquid_phase.properties_out[t].pressure == pyunits.convert(
                blk.vapor_phase[t].pressure, to_units=lunits("pressure")
            )

        self.unit_pressure_balance = Constraint(
            self.flowsheet().time,
            rule=rule_pressure_balance,
            doc="Unit level pressure balance",
        )

        # Set references to balance terms at unit level
        self.heat_duty = Reference(self.liquid_phase.heat[:])

        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != MomentumBalanceType.none
        ):
            self.deltaP = Reference(self.liquid_phase.deltaP[:])

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        common_comps = (
            self.vapor_phase.component_list
            & self.liquid_phase.properties_out.component_list
        )

        for (t, j), v in self.unit_material_balance.items():
            if j in common_comps:
                iscale.constraint_scaling_transform(
                    v,
                    iscale.get_scaling_factor(
                        self.liquid_phase.mass_transfer_term[t, "Liq", j],
                        default=1,
                        warning=True,
                    ),
                )
            elif j in self.vapor_phase.component_list:
                iscale.constraint_scaling_transform(v, value(1 / self.zero_flow_param))
            else:
                pass  # no need to scale this constraint

        for (t, j), v in self.unit_phase_equilibrium.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.liquid_phase.properties_out[t].fug_phase_comp["Liq", j],
                    default=1,
                    warning=True,
                ),
            )

        for t, v in self.unit_temperature_equality.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.liquid_phase.properties_out[t].temperature,
                    default=1,
                    warning=True,
                ),
            )

        for t, v in self.unit_enthalpy_balance.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.liquid_phase.enthalpy_transfer[t], default=1, warning=True
                ),
            )

        for t, v in self.unit_pressure_balance.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.liquid_phase.properties_out[t].pressure,
                    default=1,
                    warning=True,
                ),
            )

    def initialize(
        blk,
        liquid_state_args=None,
        vapor_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for solvent reboiler unit model.

        Keyword Arguments:
            liquid_state_args : a dict of arguments to be passed to the
                liquid property packages to provide an initial state for
                initialization (see documentation of the specific property
                package) (default = none).
            vapor_state_args : a dict of arguments to be passed to the
                vapor property package to provide an initial state for
                initialization (see documentation of the specific property
                package) (default = none).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default IDAES solver)

        Returns:
            None
        """
        if optarg is None:
            optarg = {}

        # Check DOF
        if degrees_of_freedom(blk) != 0:
            raise InitializationError(
                f"{blk.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(blk)}"
            )

        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        solverobj = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize liquid phase control volume block
        flags = blk.liquid_phase.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=liquid_state_args,
            hold_state=True,
        )

        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize vapor phase state block
        if vapor_state_args is None:
            t_init = blk.flowsheet().time.first()
            vapor_state_args = {}
            vap_state_vars = blk.vapor_phase[t_init].define_state_vars()

            liq_state = blk.liquid_phase.properties_out[t_init]

            # Check for unindexed state variables
            for sv in vap_state_vars:
                if "flow" in sv:
                    # Flow varaible, assume 10% vaporization
                    if "phase_comp" in sv:
                        # Flow is indexed by phase and component
                        vapor_state_args[sv] = {}
                        for p, j in vap_state_vars[sv]:
                            if j in liq_state.component_list:
                                vapor_state_args[sv][p, j] = 0.1 * value(
                                    getattr(liq_state, sv)[p, j]
                                )
                            else:
                                vapor_state_args[sv][p, j] = 1e-8
                    elif "comp" in sv:
                        # Flow is indexed by component
                        vapor_state_args[sv] = {}
                        for j in vap_state_vars[sv]:
                            if j in liq_state.component_list:
                                vapor_state_args[sv][j] = 0.1 * value(
                                    getattr(liq_state, sv)[j]
                                )
                            else:
                                vapor_state_args[sv][j] = 1e-8
                    elif "phase" in sv:
                        # Flow is indexed by phase
                        vapor_state_args[sv] = {}
                        for p in vap_state_vars[sv]:
                            vapor_state_args[sv][p] = 0.1 * value(
                                getattr(liq_state, sv)["Liq"]
                            )
                    else:
                        vapor_state_args[sv] = 0.1 * value(getattr(liq_state, sv))
                elif "mole_frac" in sv:
                    vapor_state_args[sv] = {}
                    if "phase" in sv:
                        # Variable is indexed by phase and component
                        for p, j in vap_state_vars[sv].keys():
                            if j in liq_state.component_list:
                                vapor_state_args[sv][p, j] = value(
                                    liq_state.fug_phase_comp["Liq", j]
                                    / liq_state.pressure
                                )
                            else:
                                vapor_state_args[sv][p, j] = 1e-8
                    else:
                        for j in vap_state_vars[sv].keys():
                            if j in liq_state.component_list:
                                vapor_state_args[sv][j] = value(
                                    liq_state.fug_phase_comp["Liq", j]
                                    / liq_state.pressure
                                )
                            else:
                                vapor_state_args[sv][j] = 1e-8
                else:
                    vapor_state_args[sv] = value(getattr(liq_state, sv))

        blk.vapor_phase.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=vapor_state_args,
            hold_state=False,
        )

        init_log.info_high("Initialization Step 2 Complete.")
        # ---------------------------------------------------------------------
        # Solve unit model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solverobj.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 3 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.liquid_phase.release_state(flags, outlvl)

        if not check_optimal_termination(results):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(results)))

    # TODO : performance and stream table methods
