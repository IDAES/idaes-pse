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
IDAES model for a generic multiple-phase, multi-stage extractor cascade.
"""

# TODO: Add reactions

# Import Pyomo libraries
from pyomo.environ import Constraint, RangeSet, Set, units, Var
from pyomo.common.config import ConfigDict, ConfigValue, In

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


@declare_process_block_class("ExtractorCascade")
class ExtractorCascadeData(UnitModelBlockData):
    """
    Standard Extractor Cascade Unit Model Class
    """

    CONFIG = UnitModelBlockData.CONFIG()

    # Config dict to contain information on each phase
    _phase_config = ConfigDict()
    _phase_config.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for given phase",
            doc="""Property parameter object used to define property calculations for given phase,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    _phase_config.declare(
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
    _phase_config.declare(
        "flow_direction",
        ConfigValue(
            default=FlowDirection.forward,
            domain=In(FlowDirection),
            doc="Direction of flow in stages",
            description="FlowDirection Enum indicating direction of "
            "flow for given phase. Default=FlowDirection.forward.",
        ),
    )

    CONFIG.declare(
        "phases",
        ConfigDict(
            implicit=True,
            implicit_domain=_phase_config,
            description="Dict of phases and associated property packages",
            doc="ConfigDict with keys indicating names for each phase in system and "
            "values indicating property package and associated arguments.",
        ),
    )
    CONFIG.declare(
        "number_of_stages",
        ConfigValue(domain=int, description="Number of stages in cascade"),
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

        # Check that at least two phases were declared
        if len(self.config.phases) < 2:
            raise ConfigurationError(
                f"Extractor models must define at least two phases; received "
                f"{self.config.phases}"
            )

        if self.config.dynamic:
            raise NotImplementedError("Extractor model does not support dynamics yet.")

        # Build indexing sets
        self.stages = RangeSet(1, self.config.number_of_stages)
        self.states = RangeSet(0, self.config.number_of_stages)

        self.phase_component_interactions = Set()
        for p1 in self.config.phases:
            for p2 in self.config.phases:
                if p1 == p2:
                    continue
                # TODO: Add check for non-interacting phases
                for j in self.config.phases[p1].property_package.component_list:
                    if (
                        j in self.config.phases[p2].property_package.component_list
                        and (p2, p1, j) not in self.phase_component_interactions
                    ):
                        # Common component, assume interaction
                        self.phase_component_interactions.add((p1, p2, j))
        if len(self.phase_component_interactions) == 0:
            raise ConfigurationError(
                "No common components found in property packages. Extractor model assumes "
                "mass transfer occurs between components with the same name in different phases."
            )

        # Build state blocks
        # Placeholders for things we will get from first StateBlock
        flow_basis = None
        uom = None

        for phase, pconfig in self.config.phases.items():
            ppack = pconfig.property_package
            flow_dir = pconfig.flow_direction

            d1 = dict(**pconfig.property_package_args)
            d1["defined_state"] = False
            d0 = dict(**pconfig.property_package_args)
            d0["defined_state"] = True

            def idx_map(i):  # i = (t, s)
                if flow_dir == FlowDirection.forward and i[1] == self.states.first():
                    return 0
                if flow_dir == FlowDirection.backward and i[1] == self.states.last():
                    return 0
                return 1

            state = ppack.build_state_block(
                self.flowsheet().time,
                self.states,
                doc=f"States for phase {phase}",
                initialize={0: d0, 1: d1},
                idx_map=idx_map,
            )
            self.add_component(phase, state)

            if flow_basis is None:
                # Set unit level flow basis and units from first phase
                t = self.flowsheet().time.first()
                s = self.states.first()
                flow_basis = state[t, s].get_material_flow_basis()
                uom = state[t, s].params.get_metadata().derived_units
            else:
                # Check that flow basis is consistent
                t = self.flowsheet().time.first()
                s = self.states.first()
                if not state[t, s].get_material_flow_basis() == flow_basis:
                    raise ConfigurationError(
                        f"Property packages use different flow bases: ExtractionCascade "
                        f"requires all property packages to use the same basis. "
                        f"{phase} uses {state[t, s].get_material_flow_basis()}, "
                        f"whilst first phase uses {flow_basis}."
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
            self.phase_component_interactions,
            initialize=0,
            units=mb_units,
            doc="Interphase mass transfer term",
        )

        for phase, pconfig in self.config.phases.items():
            # Material balance for phase
            state_block = getattr(self, phase)
            component_list = state_block.component_list
            phase_list = state_block.phase_list
            pc_set = state_block.phase_component_set
            flow_dir = pconfig.flow_direction

            def material_balance_rule(b, t, s, j):
                if flow_dir == FlowDirection.forward:
                    in_state = state_block[t, b.states.prev(s)]
                    out_state = state_block[t, s]
                elif flow_dir == FlowDirection.backward:
                    in_state = state_block[t, s]
                    out_state = state_block[t, b.states.prev(s)]
                else:
                    raise BurntToast(
                        "If/else overrun when constructing material balances"
                    )

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

                for k in b.phase_component_interactions:
                    if k[0] == phase and k[2] == j:
                        # Positive mass transfer
                        rhs += b.material_transfer_term[t, s, k]
                    elif k[1] == phase and k[2] == j:
                        # Negative mass transfer
                        rhs += -b.material_transfer_term[t, s, k]

                return 0 == rhs

            mbal = Constraint(
                self.flowsheet().time,
                self.stages,
                component_list,
                rule=material_balance_rule,
            )
            self.add_component(phase + "_material_balance", mbal)

            # TODO: Optional energy and momentum balances

        # Add Ports
        for phase, pconfig in self.config.phases.items():
            sblock = getattr(self, phase)
            flow_dir = pconfig.flow_direction

            if flow_dir == FlowDirection.forward:
                inlet = self.states.first()
                outlet = self.states.last()
            elif flow_dir == FlowDirection.backward:
                inlet = self.states.last()
                outlet = self.states.first()
            else:
                raise BurntToast("If/else overrun when constructing Ports")

            in_port, _ = sblock.build_port(
                f"{phase} Inlet", slice_index=(slice(None), inlet)
            )
            self.add_component(phase + "_inlet", in_port)

            out_port, _ = sblock.build_port(
                f"{phase} Outlet", slice_index=(slice(None), outlet)
            )
            self.add_component(phase + "_outlet", out_port)

    # TODO: Initialization - use the new framework
    def initialize(self, **kwargs):
        raise NotImplemented(
            "The ExtractorCascade unit model does not support the old initialization API. "
            "Please use the new API (InitializerObjects) instead."
        )

    def _get_performance_contents(self, time_point=0):
        raise NotImplementedError()
