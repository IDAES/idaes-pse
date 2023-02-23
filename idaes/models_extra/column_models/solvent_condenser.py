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
Condenser model for solvent columns.

This is a simple model for a condenser in the case where liquid and vapor
phases have separate proeprty packages, such as the case of solvent columns.

Assumptions:
     * Steady-state only
     * Liquid phase property package has a single phase named Liq
     * Vapor phase property package has a single phase named Vap
     * Liquid and vapor phase proeprtes need not have the same component lists
"""

__author__ = "Andrew Lee, Paul Akula"

# Import Pyomo libraries
from pyomo.environ import (
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


@declare_process_block_class("SolventCondenser")
class SolventCondenserData(UnitModelBlockData):
    """
    Condenser unit for solvent column models using separate property packages
    for liquid and vpor phases.

    Unit model to condense the vapor from the top of a solvent column.
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
                f"{self.name} SolventCondenser model requires that the vapor "
                f"phase property package have a single phase named 'Vap'"
            )
        if self.config.liquid_property_package.phase_list != ["Liq"]:
            raise ConfigurationError(
                f"{self.name} SolventCondenser model requires that the liquid "
                f"phase property package have a single phase named 'Liq'"
            )

        # Check for at least one common component in component lists
        if not any(
            j in self.config.vapor_property_package.component_list
            for j in self.config.liquid_property_package.component_list
        ):
            raise ConfigurationError(
                f"{self.name} SolventCondenser model requires that the liquid "
                f"and vapor phase property packages have at least one "
                f"common component."
            )

        # ---------------------------------------------------------------------
        # Add Control Volume for the Vapor Phase
        self.vapor_phase = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.vapor_property_package,
            property_package_args=self.config.vapor_property_package_args,
        )

        self.vapor_phase.add_state_blocks(has_phase_equilibrium=True)

        # Separate liquid and vapor phases means that phase equilibrium will
        # be handled at the unit model level, thus has_phase_equilibrium is
        # False, but has_mass_transfer is True.
        self.vapor_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True,
            has_phase_equilibrium=False,
        )

        # Need to include enthalpy transfer term for the mass transfer
        self.vapor_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=True,
            has_enthalpy_transfer=True,
        )

        self.vapor_phase.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # ---------------------------------------------------------------------
        # Add single state block for liquid phase
        tmp_dict = dict(**self.config.liquid_property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = False
        self.liquid_phase = self.config.liquid_property_package.build_state_block(
            self.flowsheet().time, doc="Liquid phase properties", **tmp_dict
        )

        # ---------------------------------------------------------------------
        # Check flow basis is compatable
        # TODO : Could add code to convert flow bases, but not now
        t_init = self.flowsheet().time.first()
        if (
            self.liquid_phase[t_init].get_material_flow_basis()
            != self.vapor_phase.properties_out[t_init].get_material_flow_basis()
        ):
            raise ConfigurationError(
                f"{self.name} vapor and liquid property packages must use the "
                f"same material flow basis."
            )

        # ---------------------------------------------------------------------
        # Add Ports for the condenser
        self.add_inlet_port(name="inlet", block=self.vapor_phase, doc="Vapor feed")
        self.add_outlet_port(
            name="vapor_outlet", block=self.vapor_phase, doc="Vapor product"
        )
        self.add_outlet_port(
            name="reflux", block=self.liquid_phase, doc="Liquid reflux from condenser"
        )

        # ---------------------------------------------------------------------
        # Add unit level constraints
        # First, need the union and intersection of component lists
        all_comps = (
            self.liquid_phase.component_list
            | self.vapor_phase.properties_out.component_list
        )
        common_comps = (
            self.liquid_phase.component_list
            & self.vapor_phase.properties_out.component_list
        )

        # Get units for unit conversion
        vunits = self.config.vapor_property_package.get_metadata().get_derived_units
        lunits = self.config.liquid_property_package.get_metadata().get_derived_units
        flow_basis = self.liquid_phase[t_init].get_material_flow_basis()
        if flow_basis == MaterialFlowBasis.molar:
            fb = "flow_mole"
        elif flow_basis == MaterialFlowBasis.molar:
            fb = "flow_mass"
        else:
            raise ConfigurationError(
                f"{self.name} SolventCondenser only supports mass or molar "
                f"basis for MaterialFlowBasis."
            )

        if any(j not in common_comps for j in self.liquid_phase.component_list):
            # We have non-volatile components present, need zero-flow param
            self.zero_flow_param = Param(
                mutable=True, default=1e-8, units=lunits("flow_mole")
            )

        # Material balances
        def rule_material_balance(blk, t, j):
            if j in common_comps:
                # Component is in equilibrium
                # Mass transfer equals liquid flowrate
                return -blk.vapor_phase.mass_transfer_term[
                    t, "Vap", j
                ] == pyunits.convert(
                    blk.liquid_phase[t].get_material_flow_terms("Liq", j),
                    to_units=vunits(fb),
                )
            elif j in self.liquid_phase.component_list:
                # Non-volatile component
                # No mass transfer term
                # Set liquid flowrate to an arbitary small value
                return (
                    blk.liquid_phase[t].get_material_flow_terms("Liq", j)
                    == blk.zero_flow_param
                )
            else:
                # Non-condensable comonent
                # Mass transfer term is zero, no vapor flowrate
                return blk.vapor_phase.mass_transfer_term[t, "Vap", j] == 0 * vunits(fb)

        self.unit_material_balance = Constraint(
            self.flowsheet().time,
            all_comps,
            rule=rule_material_balance,
            doc="Unit level material balances",
        )

        # Phase equilibrium constraints
        # For all common components, equate fugacity in vapor and liquid
        def rule_phase_equilibrium(blk, t, j):
            return blk.vapor_phase.properties_out[t].fug_phase_comp[
                "Vap", j
            ] == pyunits.convert(
                blk.liquid_phase[t].fug_phase_comp["Liq", j],
                to_units=vunits("pressure"),
            )

        self.unit_phase_equilibrium = Constraint(
            self.flowsheet().time,
            common_comps,
            rule=rule_phase_equilibrium,
            doc="Unit level phase equilibrium constraints",
        )

        # Temperature equality constraint
        def rule_temperature_balance(blk, t):
            return blk.vapor_phase.properties_out[t].temperature == pyunits.convert(
                blk.liquid_phase[t].temperature, to_units=vunits("temperature")
            )

        self.unit_temperature_equality = Constraint(
            self.flowsheet().time,
            rule=rule_temperature_balance,
            doc="Unit level temperature equality",
        )

        # Unit level energy balance
        # Energy leaving in liquid phase must be equal and opposite to enthalpy
        # transfer from vapor phase
        def rule_energy_balance(blk, t):
            return -blk.vapor_phase.enthalpy_transfer[t] == pyunits.convert(
                blk.liquid_phase[t].get_enthalpy_flow_terms("Liq"),
                to_units=vunits("energy") / lunits("time"),
            )

        self.unit_enthalpy_balance = Constraint(
            self.flowsheet().time,
            rule=rule_energy_balance,
            doc="Unit level enthalpy_balance",
        )

        # Pressure balance constraint
        def rule_pressure_balance(blk, t):
            return blk.vapor_phase.properties_out[t].pressure == pyunits.convert(
                blk.liquid_phase[t].pressure, to_units=vunits("pressure")
            )

        self.unit_pressure_balance = Constraint(
            self.flowsheet().time,
            rule=rule_pressure_balance,
            doc="Unit level pressure balance",
        )

        # Set references to balance terms at unit level
        self.heat_duty = Reference(self.vapor_phase.heat[:])

        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != MomentumBalanceType.none
        ):
            self.deltaP = Reference(self.vapor_phase.deltaP[:])

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        common_comps = (
            self.liquid_phase.component_list
            & self.vapor_phase.properties_out.component_list
        )

        for (t, j), v in self.unit_material_balance.items():
            if j in common_comps:
                iscale.constraint_scaling_transform(
                    v,
                    iscale.get_scaling_factor(
                        self.vapor_phase.mass_transfer_term[t, "Vap", j],
                        default=1,
                        warning=True,
                    ),
                )
            elif j in self.liquid_phase.component_list:
                iscale.constraint_scaling_transform(v, value(1 / self.zero_flow_param))
            else:
                pass  # no need to scale this constraint

        for (t, j), v in self.unit_phase_equilibrium.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.properties_out[t].fug_phase_comp["Vap", j],
                    default=1,
                    warning=True,
                ),
            )

        for t, v in self.unit_temperature_equality.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.properties_out[t].temperature,
                    default=1,
                    warning=True,
                ),
            )

        for t, v in self.unit_enthalpy_balance.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.enthalpy_transfer[t], default=1, warning=True
                ),
            )

        for t, v in self.unit_pressure_balance.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.properties_out[t].pressure, default=1, warning=True
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
        Initialization routine for solvent condenser unit model.

        Keyword Arguments:
            liquid_state_args : a dict of arguments to be passed to the
                liquid property package to provide an initial state for
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
        flags = blk.vapor_phase.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=vapor_state_args,
            hold_state=True,
        )

        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize liquid phase state block
        if liquid_state_args is None:
            t_init = blk.flowsheet().time.first()
            liquid_state_args = {}
            liq_state_vars = blk.liquid_phase[t_init].define_state_vars()

            vap_state = blk.vapor_phase.properties_out[t_init]

            # Check for unindexed state variables
            for sv in liq_state_vars:
                if "flow" in sv:
                    # Flow varaible, assume 10% condensation
                    if "phase_comp" in sv:
                        # Flow is indexed by phase and component
                        liquid_state_args[sv] = {}
                        for p, j in liq_state_vars[sv]:
                            if j in vap_state.component_list:
                                liquid_state_args[sv][p, j] = 0.1 * value(
                                    getattr(vap_state, sv)[p, j]
                                )
                            else:
                                liquid_state_args[sv][p, j] = 1e-8
                    elif "comp" in sv:
                        # Flow is indexed by component
                        liquid_state_args[sv] = {}
                        for j in liq_state_vars[sv]:
                            if j in vap_state.component_list:
                                liquid_state_args[sv][j] = 0.1 * value(
                                    getattr(vap_state, sv)[j]
                                )
                            else:
                                liquid_state_args[sv][j] = 1e-8
                    elif "phase" in sv:
                        # Flow is indexed by phase
                        liquid_state_args[sv] = {}
                        for p in liq_state_vars[sv]:
                            liquid_state_args[sv][p] = 0.1 * value(
                                getattr(vap_state, sv)["Vap"]
                            )
                    else:
                        liquid_state_args[sv] = 0.1 * value(getattr(vap_state, sv))
                elif "mole_frac" in sv:
                    liquid_state_args[sv] = {}
                    if "phase" in sv:
                        # Variable is indexed by phase and component
                        for p, j in liq_state_vars[sv].keys():
                            if j in vap_state.component_list:
                                liquid_state_args[sv][p, j] = value(
                                    vap_state.fug_phase_comp["Vap", j]
                                    / vap_state.pressure
                                )
                            else:
                                liquid_state_args[sv][p, j] = 1e-8
                    else:
                        for j in liq_state_vars[sv].keys():
                            if j in vap_state.component_list:
                                liquid_state_args[sv][j] = value(
                                    vap_state.fug_phase_comp["Vap", j]
                                    / vap_state.pressure
                                )
                            else:
                                liquid_state_args[sv][j] = 1e-8
                else:
                    liquid_state_args[sv] = value(getattr(vap_state, sv))

        blk.liquid_phase.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=liquid_state_args,
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
        blk.vapor_phase.release_state(flags, outlvl)

        # TODO : This fails in the current model
        # if not check_optimal_termination(results):
        #     raise InitializationError(
        #         f"{blk.name} failed to initialize successfully. Please check "
        #         f"the output logs for more information.")

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(results)))

    # TODO : performance and stream table methods
