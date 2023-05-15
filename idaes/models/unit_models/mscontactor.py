#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
IDAES model for a generic multiple-stream contactor unit.
"""
# Import Pyomo libraries
from pyomo.environ import Block, Constraint, RangeSet, Reals, Set, units, Var
from pyomo.common.config import ConfigDict, ConfigValue, Bool, In
from pyomo.contrib.incidence_analysis import solve_strongly_connected_components

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    FlowDirection,
    MaterialFlowBasis,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    BurntToast,
    PropertyNotSupportedError,
)
from idaes.core.initialization import ModularInitializerBase
from idaes.core.initialization.initializer_base import StoreState
from idaes.core.solvers import get_solver
from idaes.core.util.model_serializer import to_json, from_json
import idaes.logger as idaeslog

__author__ = "Andrew Lee"

# TODO: Initializer object
# TODO: Could look at using Pyomo DAE for the length domain, but this would make
# it harder to do side feeds.


class MSContactorInitializer(ModularInitializerBase):
    """
    This is a general purpose sequential-modular Initializer object for
    multi-stream contactor unit models.

    This routine starts by deactivating any constraints that are not
    part of the base model and fixing all inter-stream transfer variables.
    The model is then solved using the Pyomo ssc_solver function
    to initialize each stream separately.

    The inter-stream transfer variables are then unfixed and the additional
    constraints reactivated, and the full model solved using the user-specified
    solver.

    """

    CONFIG = ModularInitializerBase.CONFIG()

    CONFIG.declare(
        "ssc_solver_options",
        ConfigDict(
            implicit=True,
            description="Dict of arguments for solver calls by ssc_solver",
        ),
    )
    CONFIG.declare(
        "calculate_variable_options",
        ConfigDict(
            implicit=True,
            description="Dict of options to pass to 1x1 block solver",
            doc="Dict of options to pass to calc_var_kwds argument in "
            "scc_solver method.",
        ),
    )

    def initialization_routine(
        self,
        model: Block,
    ):
        """
        Initialization routine for MSContactor Blocks.

        Args:
            model: model to be initialized

        Returns:
            None
        """
        # Get loggers
        init_log = idaeslog.getInitLogger(
            model.name, self.get_output_level(), tag="unit"
        )
        solve_log = idaeslog.getSolveLogger(
            model.name, self.get_output_level(), tag="unit"
        )

        # Get current model state
        initial_state = to_json(model, wts=StoreState, return_dict=True)

        # Isolate streams by fixing inter-stream variables
        model.material_transfer_term.fix()

        # Deactivate any constraints that are not part of the base model
        # First, build list of names for known constraints
        const_names = []
        for s in model.config.streams.keys():
            const_names.append(s + "_rate_reaction_constraint")
            const_names.append(s + "_equilibrium_reaction_constraint")
            const_names.append(s + "_inherent_reaction_constraint")
            const_names.append(s + "_material_balance")
            const_names.append(s + "_energy_balance")
            const_names.append(s + "_pressure_balance")
            const_names.append(s + "_side_stream_pressure_balance")

        # Iterate through all constraints attached to model - do not search sub-blocks
        for c in model.component_objects(Constraint, descend_into=False):
            # Deactivate constraint if its name is not in list of known names
            if c.name not in const_names:
                c.deactivate()

        # Call css_solver
        solver = get_solver(self.config.solver, options=self.config.solver_options)
        solve_strongly_connected_components(
            model,
            solver=solver,
            solve_kwds=self.config.ssc_solver_options,
            calc_var_kwds=self.config.calculate_variable_options,
        )
        init_log.info("Stream Initialization Completed.")

        # Revert state
        from_json(model, sd=initial_state, wts=StoreState)

        # Solve full model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(model, tee=slc.tee)
        init_log.info(f"Initialization Completed, {idaeslog.condition(res)}")

        return res


# Config dict to contain information on each stream
STREAM_CONFIG = ConfigDict()
STREAM_CONFIG.declare(
    "property_package",
    ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for given stream",
        doc="""Property parameter object used to define property calculations for given stream,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
    ),
)
STREAM_CONFIG.declare(
    "property_package_args",
    ConfigDict(
        implicit=True,
        description="Dict of arguments to use for constructing property package",
        doc="""A ConfigDict with arguments to be passed to property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
    ),
)
STREAM_CONFIG.declare(
    "reaction_package",
    ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package to use for stream",
        doc="""Reaction parameter object used to define reaction calculations
for stream. **default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
    ),
)
STREAM_CONFIG.declare(
    "reaction_package_args",
    ConfigDict(
        implicit=True,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
    ),
)
STREAM_CONFIG.declare(
    "flow_direction",
    ConfigValue(
        default=FlowDirection.forward,
        domain=In(FlowDirection),
        doc="Direction of flow for stream",
        description="FlowDirection Enum indicating direction of "
        "flow for given stream. Default=FlowDirection.forward.",
    ),
)
STREAM_CONFIG.declare(
    "has_feed",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether stream has a feed.",
        description="Bool indicating whether stream has a feed Port and inlet "
        "state, or if all flow is provided via mass transfer. Default=True.",
    ),
)
STREAM_CONFIG.declare(
    "has_rate_reactions",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether rate-based reactions occur in stream.",
    ),
)
STREAM_CONFIG.declare(
    "has_equilibrium_reactions",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether equilibrium-based reactions occur in stream.",
    ),
)
STREAM_CONFIG.declare(
    "has_energy_balance",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether to include energy balance for stream. Default=True.",
    ),
)
STREAM_CONFIG.declare(
    "has_heat_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether to include external heat transfer terms in energy "
        "balance for stream. Default=False.",
    ),
)
STREAM_CONFIG.declare(
    "has_heat_of_reaction",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether heat of reaction terms should be included in energy "
        "balance for stream (required reactions). Default=False.",
    ),
)
STREAM_CONFIG.declare(
    "has_pressure_balance",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether to include pressure balance for stream. Default=True.",
    ),
)
STREAM_CONFIG.declare(
    "has_pressure_change",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether to include deltaP terms in pressure balance for "
        "stream. Default=False.",
    ),
)
STREAM_CONFIG.declare(
    "side_streams",
    ConfigValue(
        default=None,
        domain=list,
        doc="List of finite elements at which a side stream should be included.",
    ),
)


@declare_process_block_class("MSContactor")
class MSContactorData(UnitModelBlockData):
    """
    Multi-Stream Contactor Unit Model Class
    """

    # Set default initializer
    default_initializer = MSContactorInitializer

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "streams",
        ConfigDict(
            implicit=True,
            implicit_domain=STREAM_CONFIG,
            description="Dict of streams and associated property packages",
            doc="ConfigDict with keys indicating names for each stream in system and "
            "values indicating property package and associated arguments.",
        ),
    )
    CONFIG.declare(
        "number_of_finite_elements",
        ConfigValue(domain=int, description="Number of finite elements to use"),
    )
    CONFIG.declare(
        "interacting_streams",
        ConfigValue(
            domain=list,
            doc="List of interacting stream pairs as 2-tuples ('stream1', 'stream2').",
        ),
    )

    def build(self):
        """
        Begin building model.

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build
        super().build()

        self._verify_inputs()
        flow_basis, uom = self._build_state_blocks()
        self._build_material_balance_constraints(flow_basis, uom)
        self._build_energy_balance_constraints(uom)
        self._build_pressure_balance_constraints(uom)
        self._build_ports()

    def _verify_inputs(self):
        # Check that at least two streams were declared
        if len(self.config.streams) < 2:
            raise ConfigurationError(
                f"MSContactor models must define at least two streams; received "
                f"{list(self.config.streams.keys())}"
            )

        if self.config.dynamic:
            raise NotImplementedError(
                "MSContactor model does not support dynamics yet."
            )

        # Build indexing sets
        self.elements = RangeSet(
            1,
            self.config.number_of_finite_elements,
            doc="Set of finite elements in cascade (1 to number of elements)",
        )

        interacting_streams = self.config.interacting_streams
        # If user did not provide interacting streams list, assume all steams interact
        if interacting_streams is None:
            interacting_streams = []
            for s1 in self.config.streams:
                for s2 in self.config.streams:
                    if not s1 == s2 and (s2, s1) not in interacting_streams:
                        interacting_streams.append((s1, s2))
        self.stream_interactions = Set(
            initialize=interacting_streams, doc="Set of interacting streams."
        )

        self.stream_component_interactions = Set(
            doc="Set of interacting components between streams."
        )
        for (stream1, stream2) in self.stream_interactions:
            for j in self.config.streams[stream1].property_package.component_list:
                if (
                    j in self.config.streams[stream2].property_package.component_list
                    and (stream2, stream1, j) not in self.stream_component_interactions
                ):
                    # Common component, assume interaction
                    self.stream_component_interactions.add((stream1, stream2, j))
        if len(self.stream_component_interactions) == 0:
            raise ConfigurationError(
                "No common components found in property packages. MSContactor model assumes "
                "mass transfer occurs between components with the same name in different streams."
            )

        # Check that reaction block was provided if reactions requested
        for s, sconfig in self.config.streams.items():
            if (
                sconfig.has_rate_reactions or sconfig.has_equilibrium_reactions
            ) and sconfig.reaction_package is None:
                raise ConfigurationError(
                    f"Stream {s} was set to include reactions, but no reaction package was provided."
                )

    def _build_state_blocks(self):
        # Build state blocks
        # Placeholders for things we will get from first StateBlock
        flow_basis = None
        uom = None

        for stream, pconfig in self.config.streams.items():
            ppack = pconfig.property_package

            arg_dict1 = dict(**pconfig.property_package_args)
            arg_dict1["defined_state"] = False

            state = ppack.build_state_block(
                self.flowsheet().time,
                self.elements,
                doc=f"States for stream {stream} in each stage.",
                **arg_dict1,
            )
            self.add_component(stream, state)

            # Add feed state if required
            if pconfig.has_feed:
                arg_dict0 = dict(**pconfig.property_package_args)
                arg_dict0["defined_state"] = True

                inlet_state = ppack.build_state_block(
                    self.flowsheet().time,
                    doc=f"Inlet states for stream {stream}.",
                    **arg_dict0,
                )
                self.add_component(stream + "_inlet_state", inlet_state)

            # Add side streams if required
            if pconfig.side_streams is not None:
                # First, verify side stream set is a sub-set of elements
                for ss in pconfig.side_streams:
                    if ss not in self.elements:
                        raise ConfigurationError(
                            f"side_streams must be a sub-set of the set of elements. "
                            f"Found {ss} in side_streams which is not in elements."
                        )
                # Add indexing Set for side streams
                ss_set = Set(initialize=pconfig.side_streams)
                self.add_component(stream + "_side_stream_set", ss_set)

                ss_state = ppack.build_state_block(
                    self.flowsheet().time,
                    ss_set,
                    doc=f"States for {stream} side streams in each stage.",
                    **arg_dict1,
                )
                self.add_component(stream + "_side_stream_state", ss_state)

            tref = self.flowsheet().time.first()
            sref = self.elements.first()
            if flow_basis is None:
                # Set unit level flow basis and units from first stream

                flow_basis = state[tref, sref].get_material_flow_basis()
                uom = state[tref, sref].params.get_metadata().derived_units
            else:
                # Check that flow bases are consistent
                if not state[tref, sref].get_material_flow_basis() == flow_basis:
                    raise ConfigurationError(
                        f"Property packages use different flow bases: ExtractionCascade "
                        f"requires all property packages to use the same basis. "
                        f"{stream} uses {state[tref, sref].get_material_flow_basis()}, "
                        f"whilst first stream uses {flow_basis}."
                    )

            # Build reactions blocks if provided
            if pconfig.reaction_package is not None:
                tmp_dict = dict(**pconfig.reaction_package_args)
                tmp_dict["state_block"] = state
                tmp_dict["has_equilibrium"] = pconfig.has_equilibrium_reactions

                reactions = pconfig.reaction_package.build_reaction_block(
                    self.flowsheet().time,
                    self.elements,
                    doc=f"Reaction properties for stream {stream}",
                    **tmp_dict,
                )
                self.add_component(stream + "_reactions", reactions)

        return flow_basis, uom

    def _build_material_balance_constraints(self, flow_basis, uom):
        # Get units for transfer terms
        if flow_basis is MaterialFlowBasis.molar:
            mb_units = uom.FLOW_MOLE
        elif flow_basis is MaterialFlowBasis.mass:
            mb_units = uom.FLOW_MASS
        else:
            # Flow type other, so cannot determine units
            mb_units = None

        # Material transfer terms are indexed by stream pairs and components.
        # Convention is that a positive material flow term indicates flow into
        # the first stream from the second stream for a given component.
        self.material_transfer_term = Var(
            self.flowsheet().time,
            self.elements,
            self.stream_component_interactions,
            initialize=0,
            units=mb_units,
            doc="Inter-stream mass transfer term",
        )

        # Build balance equations
        for stream, sconfig in self.config.streams.items():
            state_block = getattr(self, stream)
            component_list = state_block.component_list
            phase_list = state_block.phase_list
            pc_set = state_block.phase_component_set

            # Get reaction block if present
            if hasattr(self, stream + "_reactions"):
                reaction_block = getattr(self, stream + "_reactions")

            # Add equilibrium reaction terms (if required)
            if sconfig.has_rate_reactions:
                if not hasattr(sconfig.reaction_package, "rate_reaction_idx"):
                    raise PropertyNotSupportedError(
                        f"Reaction package for {stream} does not contain a list of "
                        "rate reactions (rate_reaction_idx), thus "
                        "does not support rate-based reactions."
                    )
                # Add extents of reaction and stoichiometric constraints
                # We will assume the user will define how extent will be calculated
                rate_reaction_extent = Var(
                    self.flowsheet().time,
                    self.elements,
                    sconfig.reaction_package.rate_reaction_idx,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Extent of rate-based reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_rate_reaction_extent",
                    rate_reaction_extent,
                )

                rate_reaction_generation = Var(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Generation due to rate-based reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_rate_reaction_generation",
                    rate_reaction_generation,
                )

                def rate_reaction_rule(b, t, s, p, j):
                    if (p, j) in pc_set:
                        return rate_reaction_generation[t, s, p, j] == (
                            sum(
                                reaction_block[t, s].params.rate_reaction_stoichiometry[
                                    r, p, j
                                ]
                                * rate_reaction_extent[t, s, r]
                                for r in sconfig.reaction_package.rate_reaction_idx
                            )
                        )
                    return Constraint.Skip

                rate_reaction_constraint = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    doc=f"Rate-based reaction stoichiometry for stream {stream}",
                    rule=rate_reaction_rule,
                )
                self.add_component(
                    stream + "_rate_reaction_constraint",
                    rate_reaction_constraint,
                )

            # Add equilibrium reaction terms (if required)
            if sconfig.has_equilibrium_reactions:
                if not hasattr(sconfig.reaction_package, "equilibrium_reaction_idx"):
                    raise PropertyNotSupportedError(
                        f"Reaction package for {stream} does not contain a list of "
                        "equilibrium reactions (equilibrium_reaction_idx), thus "
                        "does not support equilibrium-based reactions."
                    )
                # Add extents of reaction and stoichiometric constraints
                equilibrium_reaction_extent = Var(
                    self.flowsheet().time,
                    self.elements,
                    sconfig.reaction_package.equilibrium_reaction_idx,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Extent of equilibrium reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_equilibrium_reaction_extent",
                    equilibrium_reaction_extent,
                )

                equilibrium_reaction_generation = Var(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Generation due to equilibrium reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_equilibrium_reaction_generation",
                    equilibrium_reaction_generation,
                )

                def equilibrium_reaction_rule(b, t, s, p, j):
                    if (p, j) in pc_set:
                        return equilibrium_reaction_generation[t, s, p, j] == (
                            sum(
                                reaction_block[
                                    t, s
                                ].params.equilibrium_reaction_stoichiometry[r, p, j]
                                * equilibrium_reaction_extent[t, s, r]
                                for r in sconfig.reaction_package.equilibrium_reaction_idx
                            )
                        )
                    return Constraint.Skip

                equilibrium_reaction_constraint = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    doc=f"Equilibrium reaction stoichiometry for stream {stream}",
                    rule=equilibrium_reaction_rule,
                )
                self.add_component(
                    stream + "_equilibrium_reaction_constraint",
                    equilibrium_reaction_constraint,
                )

            # Inherent reaction terms (if required)
            if state_block.include_inherent_reactions:
                # Add extents of reaction and stoichiometric constraints
                inherent_reaction_extent = Var(
                    self.flowsheet().time,
                    self.elements,
                    sconfig.property_package.inherent_reaction_idx,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Extent of inherent reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_inherent_reaction_extent",
                    inherent_reaction_extent,
                )

                inherent_reaction_generation = Var(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Generation due to inherent reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_inherent_reaction_generation",
                    inherent_reaction_generation,
                )

                def inherent_reaction_rule(b, t, s, p, j):
                    if (p, j) in pc_set:
                        return inherent_reaction_generation[t, s, p, j] == (
                            sum(
                                state_block[
                                    t, s
                                ].params.inherent_reaction_stoichiometry[r, p, j]
                                * inherent_reaction_extent[t, s, r]
                                for r in sconfig.property_package.inherent_reaction_idx
                            )
                        )
                    return Constraint.Skip

                inherent_reaction_constraint = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    doc=f"Inherent reaction stoichiometry for stream {stream}",
                    rule=inherent_reaction_rule,
                )
                self.add_component(
                    stream + "_inherent_reaction_constraint",
                    inherent_reaction_constraint,
                )

            # Material balance for stream
            def material_balance_rule(b, t, s, j):
                in_state, out_state, side_state = _get_state_blocks(b, t, s, stream)

                if in_state is not None:
                    rhs = sum(
                        in_state.get_material_flow_terms(p, j)
                        for p in phase_list
                        if (p, j) in pc_set
                    ) - sum(
                        out_state.get_material_flow_terms(p, j)
                        for p in phase_list
                        if (p, j) in pc_set
                    )
                else:
                    rhs = -sum(
                        out_state.get_material_flow_terms(p, j)
                        for p in phase_list
                        if (p, j) in pc_set
                    )

                # Add side stream energy flow if required
                if side_state is not None:
                    rhs += sum(
                        side_state.get_material_flow_terms(p, j)
                        for p in phase_list
                        if (p, j) in pc_set
                    )

                # As overall units come from the first StateBlock constructed, we
                # cannot guarantee that any units are consistent, so convert all flow terms
                if mb_units is not None:
                    rhs = units.convert(rhs, mb_units)

                # Add mass transfer terms
                for k in b.stream_component_interactions:
                    if k[0] == stream and k[2] == j:
                        # Positive mass transfer
                        rhs += b.material_transfer_term[t, s, k]
                    elif k[1] == stream and k[2] == j:
                        # Negative mass transfer
                        rhs += -b.material_transfer_term[t, s, k]

                # Add rate reactions (if required)
                if sconfig.has_rate_reactions:
                    rhs += sum(rate_reaction_generation[t, s, p, j] for p in phase_list)

                # Add equilibrium reactions (if required)
                if sconfig.has_equilibrium_reactions:
                    rhs += sum(
                        equilibrium_reaction_generation[t, s, p, j] for p in phase_list
                    )

                # Add inherent reactions (if required)
                if state_block.include_inherent_reactions:
                    rhs += sum(
                        inherent_reaction_generation[t, s, p, j] for p in phase_list
                    )

                return 0 == rhs

            mbal = Constraint(
                self.flowsheet().time,
                self.elements,
                component_list,
                rule=material_balance_rule,
            )
            self.add_component(stream + "_material_balance", mbal)

    def _build_energy_balance_constraints(self, uom):
        # Energy Balances

        # Assume that if energy balances are enabled that energy transfer
        # occurs between all interacting phases.
        # # For now, we will not distinguish different types of energy transfer.
        # Convention is that a positive material flow term indicates flow into
        # the first stream from the second stream.
        self.energy_transfer_term = Var(
            self.flowsheet().time,
            self.elements,
            self.stream_interactions,
            initialize=0,
            units=uom.POWER,
            doc="Inter-stream energy transfer term",
        )

        for stream, pconfig in self.config.streams.items():
            if pconfig.has_energy_balance:
                state_block = getattr(self, stream)
                phase_list = state_block.phase_list

                if pconfig.has_heat_transfer:
                    heat = Var(
                        self.flowsheet().time,
                        self.elements,
                        initialize=0,
                        units=uom.POWER,
                        doc=f"External heat transfer term for stream {stream}",
                    )
                    self.add_component(stream + "_heat", heat)

                def energy_balance_rule(b, t, s):
                    in_state, out_state, side_state = _get_state_blocks(b, t, s, stream)

                    if in_state is not None:
                        rhs = sum(
                            in_state.get_enthalpy_flow_terms(p) for p in phase_list
                        ) - sum(
                            out_state.get_enthalpy_flow_terms(p) for p in phase_list
                        )
                    else:
                        rhs = -sum(
                            out_state.get_enthalpy_flow_terms(p) for p in phase_list
                        )

                    # Add side stream energy flow if required
                    if side_state is not None:
                        rhs += sum(
                            side_state.get_enthalpy_flow_terms(p) for p in phase_list
                        )

                    # As overall units come from the first StateBlock constructed, we
                    # cannot guarantee that any units are consistent, so convert all flow terms
                    rhs = units.convert(rhs, uom.POWER)

                    # Add interstream transfer terms
                    for k in b.stream_interactions:
                        if k[0] == stream:
                            # Positive mass transfer
                            rhs += b.energy_transfer_term[t, s, k]
                        elif k[1] == stream:
                            # Negative mass transfer
                            rhs += -b.energy_transfer_term[t, s, k]

                    # Add external heat term if required
                    if pconfig.has_heat_transfer:
                        rhs += heat[t, s]

                    # Add heat of reaction terms if required
                    if pconfig.has_heat_of_reaction:
                        if not (
                            hasattr(b, stream + "_rate_reaction_extent")
                            or hasattr(b, stream + "_equilibrium_reaction_extent")
                        ):
                            raise ConfigurationError(
                                f"Stream {stream} was set to include heats of reaction, "
                                "but no extent of reactions terms could be found. "
                                "Please ensure that you defined a reaction package for this "
                                "stream and that the material balances were set to include "
                                "reactions."
                            )
                        reactions = getattr(b, stream + "_reactions")

                        if hasattr(b, stream + "_rate_reaction_extent"):
                            rate_extent = getattr(b, stream + "_rate_reaction_extent")
                            rhs += -sum(
                                rate_extent[t, s, r] * reactions[t, s].dh_rxn[r]
                                for r in pconfig.reaction_package.rate_reaction_idx
                            )

                        if hasattr(b, stream + "_equilibrium_reaction_extent"):
                            equil_extent = getattr(
                                b, stream + "_equilibrium_reaction_extent"
                            )
                            rhs += -sum(
                                equil_extent[t, s, e] * reactions[t, s].dh_rxn[e]
                                for e in pconfig.reaction_package.equilibrium_reaction_idx
                            )

                    return 0 == rhs

                ebal = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    rule=energy_balance_rule,
                )
                self.add_component(stream + "_energy_balance", ebal)

    def _build_pressure_balance_constraints(self, uom):
        # Pressure Balances
        for stream, pconfig in self.config.streams.items():
            if pconfig.has_pressure_balance:

                if pconfig.has_pressure_change:
                    deltaP = Var(
                        self.flowsheet().time,
                        self.elements,
                        initialize=0,
                        units=uom.PRESSURE,
                        doc=f"DeltaP term for stream {stream}",
                    )
                    self.add_component(stream + "_deltaP", deltaP)

                def pressure_balance_rule(b, t, s):
                    in_state, out_state, _ = _get_state_blocks(b, t, s, stream)

                    if in_state is None:
                        # If there is no feed, then there is no need for a pressure balance
                        return Constraint.Skip

                    rhs = in_state.pressure - out_state.pressure

                    # As overall units come from the first StateBlock constructed, we
                    # cannot guarantee that any units are consistent, so convert all flow terms
                    rhs = units.convert(rhs, uom.PRESSURE)

                    # Add deltaP term if required
                    if pconfig.has_pressure_change:
                        rhs += deltaP[t, s]

                    return 0 == rhs

                pbal = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    rule=pressure_balance_rule,
                )
                self.add_component(stream + "_pressure_balance", pbal)

                # Add side stream pressure equality if required
                if hasattr(self, stream + "_side_stream_state"):
                    side_set = getattr(self, stream + "_side_stream_set")

                    def ss_presure_rule(b, t, s):
                        stage_state = getattr(b, stream)[t, s]
                        side_state = getattr(b, stream + "_side_stream_state")[t, s]

                        return stage_state.pressure == side_state.pressure

                    side_pbal = Constraint(
                        self.flowsheet().time,
                        side_set,
                        rule=ss_presure_rule,
                    )
                    self.add_component(
                        stream + "_side_stream_pressure_balance", side_pbal
                    )

    def _build_ports(self):
        # Add Ports
        for stream, pconfig in self.config.streams.items():
            sblock = getattr(self, stream)
            flow_dir = pconfig.flow_direction

            if pconfig.has_feed:
                inlet_state = getattr(self, stream + "_inlet_state")
                in_port, _ = inlet_state.build_port(
                    f"{stream} Inlet", slice_index=(slice(None))
                )
                self.add_component(stream + "_inlet", in_port)

            if flow_dir == FlowDirection.forward:
                outlet = self.elements.last()
            elif flow_dir == FlowDirection.backward:
                outlet = self.elements.first()
            else:
                raise BurntToast("If/else overrun when constructing Ports")

            out_port, _ = sblock.build_port(
                f"{stream} Outlet", slice_index=(slice(None), outlet)
            )
            self.add_component(stream + "_outlet", out_port)

    # TODO: Initialization - use the new framework
    def initialize(self, **kwargs):
        raise NotImplementedError(
            "The MSContactor unit model does not support the old initialization API. "
            "Please use the new API (InitializerObjects) instead."
        )

    def _get_performance_contents(self, time_point=0):
        raise NotImplementedError()


def _get_state_blocks(b, t, s, stream):
    """
    Utility method for collecting states representing flows into and out of
    a stage for a given stream.
    """
    state_block = getattr(b, stream)

    if b.config.streams[stream].flow_direction == FlowDirection.forward:
        if s == b.elements.first():
            if not b.config.streams[stream].has_feed:
                in_state = None
            else:
                in_state = getattr(b, stream + "_inlet_state")[t]
        else:
            in_state = state_block[t, b.elements.prev(s)]
    elif b.config.streams[stream].flow_direction == FlowDirection.backward:
        if s == b.elements.last():
            if not b.config.streams[stream].has_feed:
                in_state = None
            else:
                in_state = getattr(b, stream + "_inlet_state")[t]
        else:
            in_state = state_block[t, b.elements.next(s)]
    else:
        raise BurntToast("If/else overrun when constructing balances")

    out_state = state_block[t, s]

    # Look for side state
    side_state = None
    if hasattr(b, stream + "_side_stream_state"):
        try:
            side_state = getattr(b, stream + "_side_stream_state")[t, s]
        except KeyError:
            pass

    return in_state, out_state, side_state
