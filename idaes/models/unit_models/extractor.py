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

# TODO: Add reactions

# Import Pyomo libraries
from pyomo.environ import Constraint, RangeSet, Set, units, Var
from pyomo.common.config import ConfigDict, ConfigList, ConfigValue, Bool, In

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


# TODO: Energy and momentum balances
# TODO: Add heat transfer, enthalpy transfer, pressure change, reactions, heat of reaction, phase change
# TODO: Options for balance types?
# TODO: Side feed/draw
# TODO: No feed streams


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
            "state, or if all flow is provided via mass transfer. DEfault=True.",
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

        # Check that at least two streams were declared
        if len(self.config.streams) < 2:
            raise ConfigurationError(
                f"Extractor models must define at least two streams; received "
                f"{self.config.streams}"
            )

        if self.config.dynamic:
            raise NotImplementedError("Extractor model does not support dynamics yet.")

        # Build indexing sets
        self.stages = RangeSet(
            1,
            self.config.number_of_stages,
            doc="Set of stages in cascade (1 to number of stages)",
        )

        self.stream_component_interactions = Set(
            doc="Set of interacting components between streams."
        )
        for p1 in self.config.streams:
            for p2 in self.config.streams:
                if p1 == p2:
                    continue
                if self.config.interacting_streams is not None and not (
                    (p1, p2) in self.config.interacting_streams
                    or (p2, p1) in self.config.interacting_streams
                ):
                    # Not an interacting stream pair
                    continue
                # Interacting stream pair
                for j in self.config.streams[p1].property_package.component_list:
                    if (
                        j in self.config.streams[p2].property_package.component_list
                        and (p2, p1, j) not in self.stream_component_interactions
                    ):
                        # Common component, assume interaction
                        self.stream_component_interactions.add((p1, p2, j))
        if len(self.stream_component_interactions) == 0:
            raise ConfigurationError(
                "No common components found in property packages. Extractor model assumes "
                "mass transfer occurs between components with the same name in different streams."
            )

        # Build state blocks
        # Placeholders for things we will get from first StateBlock
        flow_basis = None
        uom = None

        for stream, pconfig in self.config.streams.items():
            ppack = pconfig.property_package
            flow_dir = pconfig.flow_direction

            arg_dict1 = dict(**pconfig.property_package_args)
            arg_dict1["defined_state"] = False
            # d0 = dict(**pconfig.property_package_args)
            # d0["defined_state"] = True

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

            if flow_basis is None:
                # Set unit level flow basis and units from first stream
                t = self.flowsheet().time.first()
                s = self.stages.first()
                flow_basis = state[t, s].get_material_flow_basis()
                uom = state[t, s].params.get_metadata().derived_units
            else:
                # Check that flow basis is consistent
                t = self.flowsheet().time.first()
                s = self.stages.first()
                if not state[t, s].get_material_flow_basis() == flow_basis:
                    raise ConfigurationError(
                        f"Property packages use different flow bases: ExtractionCascade "
                        f"requires all property packages to use the same basis. "
                        f"{stream} uses {state[t, s].get_material_flow_basis()}, "
                        f"whilst first stream uses {flow_basis}."
                    )

        # Build material balances
        if flow_basis is MaterialFlowBasis.molar:
            mb_units = uom.FLOW_MOLE
        elif flow_basis is MaterialFlowBasis.mass:
            mb_units = uom.FLOW_MASS
        else:
            # Flow type other, so cannot determine units
            mb_units = None

        self.material_transfer_term = Var(
            self.flowsheet().time,
            self.stages,
            self.stream_component_interactions,
            initialize=0,
            units=mb_units,
            doc="Inter-stream mass transfer term",
        )

        # Build balance equations
        for stream, pconfig in self.config.streams.items():
            # Material balance for stream
            state_block = getattr(self, stream)
            component_list = state_block.component_list
            phase_list = state_block.phase_list
            pc_set = state_block.phase_component_set
            flow_dir = pconfig.flow_direction

            def material_balance_rule(b, t, s, j):
                if flow_dir == FlowDirection.forward:
                    if s == b.stages.first():
                        in_state = getattr(self, stream + "_inlet_state")[t]
                    else:
                        in_state = state_block[t, b.stages.prev(s)]
                elif flow_dir == FlowDirection.backward:
                    if s == b.stages.last():
                        in_state = getattr(self, stream + "_inlet_state")[t]
                    else:
                        in_state = state_block[t, b.stages.next(s)]
                else:
                    raise BurntToast("If/else overrun when constructing balances")
                out_state = state_block[t, s]

                # As overall units come from the first StateBlock constructed, we
                # cannot guarantee that any units are consistent, so convert all flow terms
                rhs = sum(
                    in_state.get_material_flow_terms(p, j)
                    for p in phase_list
                    if (p, j) in pc_set
                ) - sum(
                    out_state.get_material_flow_terms(p, j)
                    for p in phase_list
                    if (p, j) in pc_set
                )
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

            # TODO: Optional energy and momentum balances

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
        raise NotImplemented(
            "The ExtractorCascade unit model does not support the old initialization API. "
            "Please use the new API (InitializerObjects) instead."
        )

    def _get_performance_contents(self, time_point=0):
        raise NotImplementedError()
