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
IDAES model for a generic multiple-stream, multi-stage extractor cascade.
"""
# Import Pyomo libraries
from pyomo.environ import Constraint, RangeSet, Set, units, Var
from pyomo.common.config import ConfigDict, ConfigValue, Bool, In

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
)
from idaes.core.util.exceptions import ConfigurationError, BurntToast

__author__ = "Andrew Lee"

# TODO: Initializer object
# TODO: Consider renaming stages to elements
# TODO: Consider possibility of using Pyomo DAE for elements (makes side feeds harder to implement?)
# TODO: Ordering of interaction terms, and n-1 issue/verification
# TODO: Add heat transfer, enthalpy transfer, pressure change, reactions, heat of reaction, phase change


@declare_process_block_class("ExtractorCascade")
class ExtractorCascadeData(UnitModelBlockData):
    """
    Standard Extractor Cascade Unit Model Class
    """

    CONFIG = UnitModelBlockData.CONFIG()

    # Config dict to contain information on each stream
    _stream_config = ConfigDict()
    _stream_config.declare(
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
    _stream_config.declare(
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
    _stream_config.declare(
        "flow_direction",
        ConfigValue(
            default=FlowDirection.forward,
            domain=In(FlowDirection),
            doc="Direction of flow in stages",
            description="FlowDirection Enum indicating direction of "
            "flow for given stream. Default=FlowDirection.forward.",
        ),
    )
    _stream_config.declare(
        "has_feed",
        ConfigValue(
            default=True,
            domain=Bool,
            doc="Bool indicating whether stream has a feed.",
            description="Bool indicating whether stream has a feed Port and inlet "
            "state, or if all flow is provided via mass transfer. Default=True.",
        ),
    )
    _stream_config.declare(
        "has_energy_balance",
        ConfigValue(
            default=True,
            domain=Bool,
            doc="Bool indicating whether to include energy balance for stream. Default=True.",
        ),
    )
    _stream_config.declare(
        "has_pressure_balance",
        ConfigValue(
            default=True,
            domain=Bool,
            doc="Bool indicating whether to include pressure balance for stream. Default=True.",
        ),
    )
    _stream_config.declare(
        "side_streams",
        ConfigValue(
            default=None,
            domain=list,
            doc="List of stages at which a side stream should be included.",
        ),
    )

    CONFIG.declare(
        "streams",
        ConfigDict(
            implicit=True,
            implicit_domain=_stream_config,
            description="Dict of streams and associated property packages",
            doc="ConfigDict with keys indicating names for each stream in system and "
            "values indicating property package and associated arguments.",
        ),
    )
    CONFIG.declare(
        "number_of_stages",
        ConfigValue(domain=int, description="Number of stages in cascade"),
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
                f"ExtractorCascade models must define at least two streams; received "
                f"{list(self.config.streams.keys())}"
            )

        if self.config.dynamic:
            raise NotImplementedError(
                "ExtractorCascade model does not support dynamics yet."
            )

        # Build indexing sets
        self.stages = RangeSet(
            1,
            self.config.number_of_stages,
            doc="Set of stages in cascade (1 to number of stages)",
        )

        self.stream_component_interactions = Set(
            doc="Set of interacting components between streams."
        )
        for stream1 in self.config.streams:
            for stream2 in self.config.streams:
                if stream1 == stream2:
                    continue
                if self.config.interacting_streams is not None and not (
                    (stream1, stream2) in self.config.interacting_streams
                    or (stream2, stream1) in self.config.interacting_streams
                ):
                    # Not an interacting stream pair, or one we already caught
                    continue
                # Interacting stream pair
                for j in self.config.streams[stream1].property_package.component_list:
                    if (
                        j
                        in self.config.streams[stream2].property_package.component_list
                        and (stream2, stream1, j)
                        not in self.stream_component_interactions
                    ):
                        # Common component, assume interaction
                        self.stream_component_interactions.add((stream1, stream2, j))
        if len(self.stream_component_interactions) == 0:
            raise ConfigurationError(
                "No common components found in property packages. Extractor model assumes "
                "mass transfer occurs between components with the same name in different streams."
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
                self.stages,
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
                # First, verify side stream set is a sub-set of stages
                for ss in pconfig.side_streams:
                    if ss not in self.stages:
                        raise ConfigurationError(
                            f"side_streams must be a sub-set of the set of stages. "
                            f"Found {ss} in side_streams which is not in stages."
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
            sref = self.stages.first()
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
            self.stages,
            self.stream_component_interactions,
            initialize=0,
            units=mb_units,
            doc="Inter-stream mass transfer term",
        )

        # Build balance equations
        for stream in self.config.streams:
            state_block = getattr(self, stream)
            component_list = state_block.component_list
            phase_list = state_block.phase_list
            pc_set = state_block.phase_component_set

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

                for k in b.stream_component_interactions:
                    if k[0] == stream and k[2] == j:
                        # Positive mass transfer
                        rhs += b.material_transfer_term[t, s, k]
                    elif k[1] == stream and k[2] == j:
                        # Negative mass transfer
                        rhs += -b.material_transfer_term[t, s, k]

                return 0 == rhs

            mbal = Constraint(
                self.flowsheet().time,
                self.stages,
                component_list,
                rule=material_balance_rule,
            )
            self.add_component(stream + "_material_balance", mbal)

    def _build_energy_balance_constraints(self, uom):
        # Energy Balances
        for stream, pconfig in self.config.streams.items():
            if pconfig.has_energy_balance:
                state_block = getattr(self, stream)
                phase_list = state_block.phase_list

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

                    # TODO: Add transfer terms

                    return 0 == rhs

                ebal = Constraint(
                    self.flowsheet().time,
                    self.stages,
                    rule=energy_balance_rule,
                )
                self.add_component(stream + "_energy_balance", ebal)

    def _build_pressure_balance_constraints(self, uom):
        # Pressure Balances
        for stream, pconfig in self.config.streams.items():
            if pconfig.has_pressure_balance:

                def pressure_balance_rule(b, t, s):
                    in_state, out_state, _ = _get_state_blocks(b, t, s, stream)

                    if in_state is None:
                        # If there is no feed, then there is no need for a pressure balance
                        return Constraint.Skip

                    rhs = in_state.pressure - out_state.pressure

                    # As overall units come from the first StateBlock constructed, we
                    # cannot guarantee that any units are consistent, so convert all flow terms
                    rhs = units.convert(rhs, uom.PRESSURE)

                    # TODO: Add deltaP terms

                    return 0 == rhs

                pbal = Constraint(
                    self.flowsheet().time,
                    self.stages,
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
                outlet = self.stages.last()
            elif flow_dir == FlowDirection.backward:
                outlet = self.stages.first()
            else:
                raise BurntToast("If/else overrun when constructing Ports")

            out_port, _ = sblock.build_port(
                f"{stream} Outlet", slice_index=(slice(None), outlet)
            )
            self.add_component(stream + "_outlet", out_port)

    # TODO: Initialization - use the new framework
    def initialize(self, **kwargs):
        raise NotImplementedError(
            "The ExtractorCascade unit model does not support the old initialization API. "
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
        if s == b.stages.first():
            if not b.config.streams[stream].has_feed:
                in_state = None
            else:
                in_state = getattr(b, stream + "_inlet_state")[t]
        else:
            in_state = state_block[t, b.stages.prev(s)]
    elif b.config.streams[stream].flow_direction == FlowDirection.backward:
        if s == b.stages.last():
            if not b.config.streams[stream].has_feed:
                in_state = None
            else:
                in_state = getattr(b, stream + "_inlet_state")[t]
        else:
            in_state = state_block[t, b.stages.next(s)]
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