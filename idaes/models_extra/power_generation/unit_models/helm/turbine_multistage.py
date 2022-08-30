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
Multistage steam turbine for power generation.

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136, November
"""
import copy

import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.common.config import ConfigBlock, ConfigValue, ConfigList, In

from idaes.core import declare_process_block_class, UnitModelBlockData, useDefault
from idaes.models_extra.power_generation.unit_models.helm import (
    HelmSplitter,
    HelmMixer,
    MomentumMixingType,
    HelmTurbineInletStage,
    HelmTurbineStage,
    HelmTurbineOutletStage,
    ValveFunctionType,
)
from idaes.models_extra.power_generation.unit_models.helm import HelmValve as SteamValve

from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def _define_turbine_multistage_config(config):
    config.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag",
            doc="Only False, in a dynamic flowsheet this is psuedo-steady-state.",
        ),
    )
    config.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag",
            doc="Only False, in a dynamic flowsheet this is psuedo-steady-state.",
        ),
    )
    config.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    config.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    config.declare(
        "num_parallel_inlet_stages",
        ConfigValue(
            default=4,
            domain=int,
            description="Number of parallel inlet stages to simulate partial arc "
            "admission.  Default=4",
        ),
    )
    config.declare(
        "throttle_valve_function",
        ConfigValue(
            default=ValveFunctionType.linear,
            domain=In(ValveFunctionType),
            description="Valve function type, if custom provide an expression rule",
            doc="""The type of valve function, if custom provide an expression rule
with the valve_function_rule argument.
**default** - ValveFunctionType.linear
**Valid values** - {
ValveFunctionType.linear,
ValveFunctionType.quick_opening,
ValveFunctionType.equal_percentage,
ValveFunctionType.custom}""",
        ),
    )
    config.declare(
        "throttle_valve_function_callback",
        ConfigValue(
            default=None,
            description="A callback to add a custom valve function to the "
            "throttle valves or None.  If a callback is provided, it should "
            "take the valve block data as an argument and add a "
            "valve_function expressions to it. Default=None",
        ),
    )
    config.declare(
        "num_hp",
        ConfigValue(
            default=2,
            domain=int,
            description="Number of high pressure stages not including inlet stage",
            doc="Number of high pressure stages not including inlet stage",
        ),
    )
    config.declare(
        "num_ip",
        ConfigValue(
            default=10,
            domain=int,
            description="Number of intermediate pressure stages",
            doc="Number of intermediate pressure stages",
        ),
    )
    config.declare(
        "num_lp",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of low pressure stages not including outlet stage",
            doc="Number of low pressure stages not including outlet stage",
        ),
    )
    config.declare(
        "hp_split_locations",
        ConfigList(
            default=[],
            domain=int,
            description="Locations of splitters in HP section",
            doc="A list of index locations of splitters in the HP section. The "
            "indexes indicate after which stage to include splitters.  0 is "
            "between the inlet stage and the first regular HP stage.",
        ),
    )
    config.declare(
        "ip_split_locations",
        ConfigList(
            default=[],
            domain=int,
            description="Locations of splitters in IP section",
            doc="A list of index locations of splitters in the IP section. The "
            "indexes indicate after which stage to include splitters.",
        ),
    )
    config.declare(
        "lp_split_locations",
        ConfigList(
            default=[],
            domain=int,
            description="Locations of splitter in LP section",
            doc="A list of index locations of splitters in the LP section. The "
            "indexes indicate after which stage to include splitters.",
        ),
    )
    config.declare(
        "hp_disconnect",
        ConfigList(
            default=[],
            domain=int,
            description="HP Turbine stages to not connect to next with an arc.",
            doc="HP Turbine stages to not connect to next with an arc. This is "
            "usually used to insert additional units between stages on a "
            "flowsheet, such as a reheater",
        ),
    )
    config.declare(
        "ip_disconnect",
        ConfigList(
            default=[],
            domain=int,
            description="IP Turbine stages to not connect to next with an arc.",
            doc="IP Turbine stages to not connect to next with an arc. This is "
            "usually used to insert additional units between stages on a "
            "flowsheet, such as a reheater",
        ),
    )
    config.declare(
        "lp_disconnect",
        ConfigList(
            default=[],
            domain=int,
            description="LP Turbine stages to not connect to next with an arc.",
            doc="LP Turbine stages to not connect to next with an arc. This is "
            "usually used to insert additional units between stages on a "
            "flowsheet, such as a reheater",
        ),
    )
    config.declare(
        "hp_split_num_outlets",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict, hp split index: number of splitter outlets, if not 2",
        ),
    )
    config.declare(
        "ip_split_num_outlets",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict, ip split index: number of splitter outlets, if not 2",
        ),
    )
    config.declare(
        "lp_split_num_outlets",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict, lp split index: number of splitter outlets, if not 2",
        ),
    )


@declare_process_block_class(
    "HelmTurbineMultistage",
    doc="Multistage steam turbine with optional reheat and extraction",
)
class HelmTurbineMultistageData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    _define_turbine_multistage_config(CONFIG)

    def build(self):
        super().build()
        config = self.config
        unit_cfg = {  # general unit model config
            "dynamic": config.dynamic,
            "has_holdup": config.has_holdup,
            "property_package": config.property_package,
            "property_package_args": config.property_package_args,
        }
        ni = self.config.num_parallel_inlet_stages
        inlet_idx = self.inlet_stage_idx = pyo.RangeSet(ni)

        thrtl_cfg = unit_cfg.copy()
        thrtl_cfg["valve_function"] = self.config.throttle_valve_function
        thrtl_cfg[
            "valve_function_callback"
        ] = self.config.throttle_valve_function_callback

        # Adding unit models
        # ------------------------

        # Splitter to inlet that splits main flow into parallel flows for
        # paritial arc admission to the turbine
        self.inlet_split = HelmSplitter(**self._split_cfg(unit_cfg, ni))
        self.throttle_valve = SteamValve(inlet_idx, **thrtl_cfg)
        self.inlet_stage = HelmTurbineInletStage(inlet_idx, **unit_cfg)
        # mixer to combine the parallel flows back together
        self.inlet_mix = HelmMixer(**self._mix_cfg(unit_cfg, ni))
        # add turbine sections.
        # inlet stage -> hp stages -> ip stages -> lp stages -> outlet stage
        self.hp_stages = HelmTurbineStage(pyo.RangeSet(config.num_hp), **unit_cfg)
        self.ip_stages = HelmTurbineStage(pyo.RangeSet(config.num_ip), **unit_cfg)
        self.lp_stages = HelmTurbineStage(pyo.RangeSet(config.num_lp), **unit_cfg)
        self.outlet_stage = HelmTurbineOutletStage(**unit_cfg)

        for i in self.hp_stages:
            self.hp_stages[i].ratioP.fix()
            self.hp_stages[i].efficiency_isentropic.fix()
        for i in self.ip_stages:
            self.ip_stages[i].ratioP.fix()
            self.ip_stages[i].efficiency_isentropic.fix()
        for i in self.lp_stages:
            self.lp_stages[i].ratioP.fix()
            self.lp_stages[i].efficiency_isentropic.fix()

        # Then make splitter config.  If number of outlets is specified
        # make a specific config, otherwise use default with 2 outlets
        s_sfg_default = self._split_cfg(unit_cfg, 2)
        hp_splt_cfg = {}
        ip_splt_cfg = {}
        lp_splt_cfg = {}
        # Now to finish up if there are more than two outlets, set that
        for i, v in config.hp_split_num_outlets.items():
            hp_splt_cfg[i] = self._split_cfg(unit_cfg, v)
        for i, v in config.ip_split_num_outlets.items():
            ip_splt_cfg[i] = self._split_cfg(unit_cfg, v)
        for i, v in config.lp_split_num_outlets.items():
            lp_splt_cfg[i] = self._split_cfg(unit_cfg, v)
        # put in splitters for turbine steam extractions
        if config.hp_split_locations:
            self.hp_split = HelmSplitter(
                config.hp_split_locations, **s_sfg_default, initialize=hp_splt_cfg
            )
        else:
            self.hp_split = {}
        if config.ip_split_locations:
            self.ip_split = HelmSplitter(
                config.ip_split_locations, **s_sfg_default, initialize=ip_splt_cfg
            )
        else:
            self.ip_split = {}
        if config.lp_split_locations:
            self.lp_split = HelmSplitter(
                config.lp_split_locations, **s_sfg_default, initialize=lp_splt_cfg
            )
        else:
            self.lp_split = {}
        # Done with unit models.  Adding Arcs (streams).
        # ------------------------------------------------

        # First up add streams in the inlet section
        def _split_to_rule(b, i):
            return {
                "source": getattr(self.inlet_split, "outlet_{}".format(i)),
                "destination": self.throttle_valve[i].inlet,
            }

        def _valve_to_rule(b, i):
            return {
                "source": self.throttle_valve[i].outlet,
                "destination": self.inlet_stage[i].inlet,
            }

        def _inlet_to_rule(b, i):
            return {
                "source": self.inlet_stage[i].outlet,
                "destination": getattr(self.inlet_mix, "inlet_{}".format(i)),
            }

        self.stream_throttle_inlet = Arc(inlet_idx, rule=_split_to_rule)
        self.stream_throttle_outlet = Arc(inlet_idx, rule=_valve_to_rule)
        self.stream_inlet_mix_inlet = Arc(inlet_idx, rule=_inlet_to_rule)

        # There are three sections HP, IP, and LP which all have the same sort
        # of internal connctions, so the functions below provide some generic
        # capcbilities for adding the internal Arcs (streams).
        def _arc_indexes(nstages, index_set, discon, splits):
            """
            This takes the index set of all possible streams in a turbine
            section and throws out arc indexes for stages that are disconnected
            and arc indexes that are not needed because there is no splitter
            after a stage.

            Args:
                nstages (int): Number of stages in section
                index_set (Set): Index set for arcs in the section
                discon (list): Disconnected stages in the section
                splits (list): Spliter locations
            """
            sr = set()  # set of things to remove from the Arc index set
            for i in index_set:
                if (i[0] in discon or i[0] == nstages) and i[0] in splits:
                    # don't connect stage i to next remove stream after split
                    sr.add((i[0], 2))
                elif (i[0] in discon or i[0] == nstages) and i[0] not in splits:
                    # no splitter and disconnect so remove both streams
                    sr.add((i[0], 1))
                    sr.add((i[0], 2))
                elif i[0] not in splits:
                    # no splitter and not disconnected so just second stream
                    sr.add((i[0], 2))
                else:
                    # has splitter so need both streams don't remove anything
                    pass
            for i in sr:  # remove the unneeded Arc indexes
                index_set.remove(i)

        def _arc_rule(turbines, splitters):
            """
            This creates a rule function for arcs in a turbine section. When
            this is used, the indexes for nonexistant stream will have already
            been removed, so any indexes the rule will get should have a stream
            associated.

            Args:
                turbines (TurbineStage): Indexed block with turbine section stages
                splitters (Separator): Indexed block of splitters
            """

            def _rule(b, i, j):
                if i in splitters and j == 1:  # stage to splitter
                    return {
                        "source": turbines[i].outlet,
                        "destination": splitters[i].inlet,
                    }
                elif j == 2:  # splitter to next stage
                    return {
                        "source": splitters[i].outlet_1,
                        "destination": turbines[i + 1].inlet,
                    }
                else:  # no splitter, stage to next stage
                    return {
                        "source": turbines[i].outlet,
                        "destination": turbines[i + 1].inlet,
                    }

            return _rule

        # Create initial arcs index sets with all possible streams
        self.hp_stream_idx = pyo.Set(initialize=self.hp_stages.index_set() * [1, 2])
        self.ip_stream_idx = pyo.Set(initialize=self.ip_stages.index_set() * [1, 2])
        self.lp_stream_idx = pyo.Set(initialize=self.lp_stages.index_set() * [1, 2])

        # Throw out unneeded streams for disconnected stages or no splitter
        _arc_indexes(
            config.num_hp,
            self.hp_stream_idx,
            config.hp_disconnect,
            config.hp_split_locations,
        )
        _arc_indexes(
            config.num_ip,
            self.ip_stream_idx,
            config.ip_disconnect,
            config.ip_split_locations,
        )
        _arc_indexes(
            config.num_lp,
            self.lp_stream_idx,
            config.lp_disconnect,
            config.lp_split_locations,
        )

        # Create connections internal to each turbine section (hp, ip, and lp)
        self.hp_stream = Arc(
            self.hp_stream_idx, rule=_arc_rule(self.hp_stages, self.hp_split)
        )
        self.ip_stream = Arc(
            self.ip_stream_idx, rule=_arc_rule(self.ip_stages, self.ip_split)
        )
        self.lp_stream = Arc(
            self.lp_stream_idx, rule=_arc_rule(self.lp_stages, self.lp_split)
        )

        # Connect hp section to ip section unless its a disconnect location
        last_hp = config.num_hp
        if 0 not in config.ip_disconnect and last_hp not in config.hp_disconnect:
            # Not disconnected stage so add stream, depending on splitter existance
            if last_hp in config.hp_split_locations:  # connect splitter to ip
                self.hp_to_ip_stream = Arc(
                    source=self.hp_split[last_hp].outlet_1,
                    destination=self.ip_stages[1].inlet,
                )
            else:  # connect last hp to ip
                self.hp_to_ip_stream = Arc(
                    source=self.hp_stages[last_hp].outlet,
                    destination=self.ip_stages[1].inlet,
                )
        # Connect ip section to lp section unless its a disconnect location
        last_ip = config.num_ip
        if 0 not in config.lp_disconnect and last_ip not in config.ip_disconnect:
            if last_ip in config.ip_split_locations:  # connect splitter to ip
                self.ip_to_lp_stream = Arc(
                    source=self.ip_split[last_ip].outlet_1,
                    destination=self.lp_stages[1].inlet,
                )
            else:  # connect last hp to ip
                self.ip_to_lp_stream = Arc(
                    source=self.ip_stages[last_ip].outlet,
                    destination=self.lp_stages[1].inlet,
                )
        # Connect inlet stage to hp section
        #   not allowing disconnection of inlet and first regular hp stage
        if 0 in config.hp_split_locations:
            # connect inlet mix to splitter and splitter to hp section
            self.inlet_to_splitter_stream = Arc(
                source=self.inlet_mix.outlet, destination=self.hp_split[0].inlet
            )
            self.splitter_to_hp_stream = Arc(
                source=self.hp_split[0].outlet_1, destination=self.hp_stages[1].inlet
            )
        else:  # connect mixer to first hp turbine stage
            self.inlet_to_hp_stream = Arc(
                source=self.inlet_mix.outlet, destination=self.hp_stages[1].inlet
            )

        self.power = pyo.Var(
            self.flowsheet().time,
            initialize=-1e8,
            doc="total turbine power",
            units=pyo.units.W,
        )

        @self.Constraint(self.flowsheet().time)
        def power_eqn(b, t):
            return b.power[t] == b.outlet_stage.control_volume.work[
                t
            ] * b.outlet_stage.efficiency_mech + sum(
                b.inlet_stage[i].control_volume.work[t]
                * b.inlet_stage[i].efficiency_mech
                for i in b.inlet_stage
            ) + sum(
                b.hp_stages[i].control_volume.work[t] * b.hp_stages[i].efficiency_mech
                for i in b.hp_stages
            ) + sum(
                b.ip_stages[i].control_volume.work[t] * b.ip_stages[i].efficiency_mech
                for i in b.ip_stages
            ) + sum(
                b.lp_stages[i].control_volume.work[t] * b.lp_stages[i].efficiency_mech
                for i in b.lp_stages
            )

        # Connect lp section to outlet stage, not allowing outlet stage to be
        # disconnected
        last_lp = config.num_lp
        if last_lp in config.lp_split_locations:  # connect splitter to outlet
            self.lp_to_outlet_stream = Arc(
                source=self.lp_split[last_lp].outlet_1,
                destination=self.outlet_stage.inlet,
            )
        else:  # connect last lpstage to outlet
            self.lp_to_outlet_stream = Arc(
                source=self.lp_stages[last_lp].outlet,
                destination=self.outlet_stage.inlet,
            )
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

    def _split_cfg(self, unit_cfg, no=2):
        """
        This creates a configuration dictionary for a splitter.

        Args:
            unit_cfg: The base unit config dict.
            no: Number of outlets, default=2
        """
        # Create a dict for splitter config args
        cfg = copy.copy(unit_cfg)
        cfg.update(num_outlets=no)
        return cfg

    def _mix_cfg(self, unit_cfg, ni=2):
        """
        This creates a configuration dictionary for a mixer.

        Args:
            unit_cfg: The base unit config dict.
            ni: Number of inlets, default=2
        """
        cfg = copy.copy(unit_cfg)
        cfg.update(
            num_inlets=ni, momentum_mixing_type=MomentumMixingType.minimize_and_equality
        )
        return cfg

    def throttle_cv_fix(self, value):
        """
        Fix the thottle valve coefficients.  These are generally the same for
        each of the parallel stages so this provides a convenient way to set
        them.

        Args:
            value: The value to fix the turbine inlet flow coefficients at
        """
        for i in self.throttle_valve:
            self.throttle_valve[i].Cv.fix(value)

    def turbine_inlet_cf_fix(self, value):
        """
        Fix the inlet turbine stage flow coefficient.  These are
        generally the same for each of the parallel stages so this provides
        a convenient way to set them.

        Args:
            value: The value to fix the turbine inlet flow coefficients at
        """
        for i in self.inlet_stage:
            self.inlet_stage[i].flow_coeff.fix(value)

    def _init_section(
        self,
        stages,
        splits,
        disconnects,
        prev_port,
        outlvl,
        solver,
        optarg,
        copy_disconneted_flow,
        copy_disconneted_pressure,
    ):
        """Reuse the initializtion for HP, IP and, LP sections."""
        if 0 in splits:
            propagate_state(splits[0].inlet, prev_port)
            splits[0].initialize(outlvl=outlvl, solver=solver, optarg=optarg)
            prev_port = splits[0].outlet_1
        for i in stages:
            if i - 1 not in disconnects:
                propagate_state(stages[i].inlet, prev_port)
            else:
                if copy_disconneted_flow:
                    for t in stages[i].inlet.flow_mol:
                        stages[i].inlet.flow_mol[t] = pyo.value(prev_port.flow_mol[t])
                if copy_disconneted_pressure:
                    for t in stages[i].inlet.pressure:
                        stages[i].inlet.pressure[t] = pyo.value(prev_port.pressure[t])
            stages[i].initialize(outlvl=outlvl, solver=solver, optarg=optarg)
            prev_port = stages[i].outlet
            if i in splits:
                propagate_state(splits[i].inlet, prev_port)
                splits[i].initialize(outlvl=outlvl, solver=solver, optarg=optarg)
                prev_port = splits[i].outlet_1
        return prev_port

    def turbine_outlet_cf_fix(self, value):
        """
        Fix the inlet turbine stage flow coefficient.  These are
        generally the same for each of the parallel stages so this provides
        a convenient way to set them.

        Args:
            value: The value to fix the turbine inlet flow coefficients at
        """
        self.outlet_stage.flow_coeff.fix(value)

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        flow_iterate=2,
        optarg=None,
        copy_disconneted_flow=True,
        copy_disconneted_pressure=True,
        calculate_outlet_cf=False,
        calculate_inlet_cf=False,
    ):
        """
        Initialize

        Args:
            outlvl: logging level default is NOTSET, which inherits from the
                parent logger
            solver: the NL solver
            flow_iterate: If not calculating flow coefficients, this is the
                number of times to update the flow and repeat initialization
                (1 to 5 where 1 does not update the flow guess)
            optarg: solver arguments, default is None
            copy_disconneted_flow: Copy the flow through the disconnected stages
                default is True
            copy_disconneted_pressure: Copy the pressure through the disconnected
                stages default is True
            calculate_outlet_cf: Use the flow initial flow guess to calculate
                the outlet stage flow coefficient, default is False,
            calculate_inlet_cf: Use the inlet stage ratioP to calculate the flow
                coefficent for the inlet stage default is False

        Returns:
            None
        """
        # Setup loggers
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        # Store initial model specs, restored at the end of initializtion, so
        # the problem is not altered.  This can restore fixed/free vars,
        # active/inactive constraints, and fixed variable values.
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # Assume the flow into the turbine is a reasonable guess for
        # initializtion
        flow_guess = self.inlet_split.inlet.flow_mol[0].value

        for it_count in range(flow_iterate):
            self.inlet_split.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

            # Initialize valves
            for i in self.inlet_stage_idx:
                u = self.throttle_valve[i]
                propagate_state(
                    u.inlet, getattr(self.inlet_split, "outlet_{}".format(i))
                )
                u.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

            # Initialize turbine
            for i in self.inlet_stage_idx:
                u = self.inlet_stage[i]
                propagate_state(u.inlet, self.throttle_valve[i].outlet)
                u.initialize(
                    outlvl=outlvl,
                    solver=solver,
                    optarg=optarg,
                    calculate_cf=calculate_inlet_cf,
                )

            # Initialize Mixer
            self.inlet_mix.use_minimum_inlet_pressure_constraint()
            for i in self.inlet_stage_idx:
                propagate_state(
                    getattr(self.inlet_mix, "inlet_{}".format(i)),
                    self.inlet_stage[i].outlet,
                )
                getattr(self.inlet_mix, "inlet_{}".format(i)).fix()
            self.inlet_mix.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
            for i in self.inlet_stage_idx:
                getattr(self.inlet_mix, "inlet_{}".format(i)).unfix()
            self.inlet_mix.use_equal_pressure_constraint()

            prev_port = self.inlet_mix.outlet
            prev_port = self._init_section(
                self.hp_stages,
                self.hp_split,
                self.config.hp_disconnect,
                prev_port,
                outlvl,
                solver,
                optarg,
                copy_disconneted_flow=copy_disconneted_flow,
                copy_disconneted_pressure=copy_disconneted_pressure,
            )
            if len(self.hp_stages) in self.config.hp_disconnect:
                self.config.ip_disconnect.append(0)
            prev_port = self._init_section(
                self.ip_stages,
                self.ip_split,
                self.config.ip_disconnect,
                prev_port,
                outlvl,
                solver,
                optarg,
                copy_disconneted_flow=copy_disconneted_flow,
                copy_disconneted_pressure=copy_disconneted_pressure,
            )
            if len(self.ip_stages) in self.config.ip_disconnect:
                self.config.lp_disconnect.append(0)
            prev_port = self._init_section(
                self.lp_stages,
                self.lp_split,
                self.config.lp_disconnect,
                prev_port,
                outlvl,
                solver,
                optarg,
                copy_disconneted_flow=copy_disconneted_flow,
                copy_disconneted_pressure=copy_disconneted_pressure,
            )

            propagate_state(self.outlet_stage.inlet, prev_port)
            self.outlet_stage.initialize(
                outlvl=outlvl,
                solver=solver,
                optarg=optarg,
                calculate_cf=calculate_outlet_cf,
            )
            if calculate_outlet_cf:
                break
            if it_count < flow_iterate - 1:
                for t in self.inlet_split.inlet.flow_mol:
                    self.inlet_split.inlet.flow_mol[
                        t
                    ].value = self.outlet_stage.inlet.flow_mol[t].value

                    for s in self.hp_split.values():
                        for i, o in enumerate(s.outlet_list):
                            if i == 0:
                                continue
                            o = getattr(s, o)
                            self.inlet_split.inlet.flow_mol[t].value += o.flow_mol[
                                t
                            ].value
                    for s in self.ip_split.values():
                        for i, o in enumerate(s.outlet_list):
                            if i == 0:
                                continue
                            o = getattr(s, o)
                            self.inlet_split.inlet.flow_mol[t].value += o.flow_mol[
                                t
                            ].value
                    for s in self.lp_split.values():
                        for i, o in enumerate(s.outlet_list):
                            if i == 0:
                                continue
                            o = getattr(s, o)
                            self.inlet_split.inlet.flow_mol[t].value += o.flow_mol[
                                t
                            ].value

        if calculate_inlet_cf:
            # cf was probably fixed, so will have to set the value agian here
            # if you ask for it to be calculated.
            icf = {}
            for i in self.inlet_stage:
                for t in self.inlet_stage[i].flow_coeff:
                    icf[i, t] = pyo.value(self.inlet_stage[i].flow_coeff[t])
        if calculate_outlet_cf:
            ocf = pyo.value(self.outlet_stage.flow_coeff)

        from_json(self, sd=istate, wts=sp)

        if calculate_inlet_cf:
            # cf was probably fixed, so will have to set the value agian here
            # if you ask for it to be calculated.
            for t in self.inlet_stage[i].flow_coeff:
                for i in self.inlet_stage:
                    self.inlet_stage[i].flow_coeff[t] = icf[i, t]
        if calculate_outlet_cf:
            self.outlet_stage.flow_coeff = ocf

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # Add a default power scale
        # pretty safe to say power is around 100 to 1000 MW

        for t in self.power:
            if iscale.get_scaling_factor(self.power[t]) is None:
                iscale.set_scaling_factor(self.power[t], 1e-8)

        for t, c in self.power_eqn.items():
            power_scale = iscale.get_scaling_factor(
                self.power[t], default=1, warning=True
            )
            # Set power equation scale factor
            iscale.constraint_scaling_transform(c, power_scale, overwrite=False)
