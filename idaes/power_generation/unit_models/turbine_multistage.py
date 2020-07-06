##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Multistage steam turbine for power generation.

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136, November
"""
import copy

from pyomo.environ import (
    RangeSet,
    Set,
    TransformationFactory,
    Var,
    value,
    SolverFactory,
)
from pyomo.network import Arc
from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    useDefault,
)
from idaes.generic_models.unit_models import (
    Separator,
    Mixer,
    SplittingType,
    EnergySplittingType,
    MomentumMixingType,
)
from idaes.power_generation.unit_models import (
    TurbineInletStage,
    TurbineStage,
    TurbineOutletStage,
    SteamValve,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.util.misc import copy_port_values as _set_port
from pyomo.common.config import ConfigBlock, ConfigValue, In, ConfigList
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def _define_turbine_multistage_config(config):
    config.declare(
        "dynamic",
        ConfigValue(
            domain=In([True, False]),
            default=False,
            description="Dynamic model flag",
            doc="Indicates whether the model is dynamic.",
        ),
    )
    config.declare(
        "has_holdup",
        ConfigValue(
            default=useDefault,
            domain=In([useDefault, True, False]),
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )
    config.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Calculate phase equilibrium in mixed stream",
            doc="""Argument indicating whether phase equilibrium should be
calculated for the resulting mixed stream,
**default** - False.
**Valid values:** {
**True** - calculate phase equilibrium in mixed stream,
**False** - do not calculate equilibrium in mixed stream.}""",
        ),
    )
    config.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.componentTotal,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentTotal`.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
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
            "usually used to insert addtional units between stages on a "
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
            "usually used to insert addtional units between stages on a "
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
            "usually used to insert addtional units between stages on a "
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
    "TurbineMultistage",
    doc="Multistage steam turbine with optional reheat and extraction",
)
class TurbineMultistageData(UnitModelBlockData):

    CONFIG = ConfigBlock()
    _define_turbine_multistage_config(CONFIG)

    def build(self):
        super(TurbineMultistageData, self).build()
        config = self.config
        unit_cfg = {  # general unit model config
            "dynamic": config.dynamic,
            "has_holdup": config.has_holdup,
            "has_phase_equilibrium": config.has_phase_equilibrium,
            "material_balance_type": config.material_balance_type,
            "property_package": config.property_package,
            "property_package_args": config.property_package_args,
        }
        ni = self.config.num_parallel_inlet_stages
        inlet_idx = self.inlet_stage_idx = RangeSet(ni)

        # Adding unit models
        # ------------------------

        # Splitter to inlet that splits main flow into parallel flows for
        # paritial arc admission to the turbine
        self.inlet_split = Separator(default=self._split_cfg(unit_cfg, ni))
        self.throttle_valve = SteamValve(inlet_idx, default=unit_cfg)
        self.inlet_stage = TurbineInletStage(inlet_idx, default=unit_cfg)
        # mixer to combine the parallel flows back together
        self.inlet_mix = Mixer(default=self._mix_cfg(unit_cfg, ni))
        # add turbine sections.
        # inlet stage -> hp stages -> ip stages -> lp stages -> outlet stage
        self.hp_stages = TurbineStage(RangeSet(config.num_hp), default=unit_cfg)
        self.ip_stages = TurbineStage(RangeSet(config.num_ip), default=unit_cfg)
        self.lp_stages = TurbineStage(RangeSet(config.num_lp), default=unit_cfg)
        self.outlet_stage = TurbineOutletStage(default=unit_cfg)

        for i in self.hp_stages:
            self.hp_stages[i].ratioP.fix()
            self.hp_stages[i].efficiency_isentropic[:].fix()
        for i in self.ip_stages:
            self.ip_stages[i].ratioP.fix()
            self.ip_stages[i].efficiency_isentropic[:].fix()
        for i in self.lp_stages:
            self.lp_stages[i].ratioP.fix()
            self.lp_stages[i].efficiency_isentropic[:].fix()

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
            self.hp_split = Separator(
                config.hp_split_locations, default=s_sfg_default, initialize=hp_splt_cfg
            )
        if config.ip_split_locations:
            self.ip_split = Separator(
                config.ip_split_locations, default=s_sfg_default, initialize=ip_splt_cfg
            )
        if config.lp_split_locations:
            self.lp_split = Separator(
                config.lp_split_locations, default=s_sfg_default, initialize=lp_splt_cfg
            )

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

        self.split_to_valve_stream = Arc(inlet_idx, rule=_split_to_rule)
        self.valve_to_inlet_stage_stream = Arc(inlet_idx, rule=_valve_to_rule)
        self.inlet_stage_to_mix = Arc(inlet_idx, rule=_inlet_to_rule)

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
            this is used the indexes for nonexistant stream will have already
            been removed, so any indexes the rule will get should have a stream
            associated.

            Args:
                turbines (TurbineStage): Indexed block with turbine section stages
                splitters (Separator): Indexed block of splitters
            """

            def _rule(b, i, j):
                if i in splitters and j == 1:
                    return {
                        "source": turbines[i].outlet,
                        "destination": splitters[i].inlet,
                    }
                elif j == 2:
                    return {
                        "source": splitters[i].outlet_1,
                        "destination": turbines[i + 1].inlet,
                    }
                else:
                    return {
                        "source": turbines[i].outlet,
                        "destination": turbines[i + 1].inlet,
                    }

            return _rule

        # Create initial arcs index sets with all possible streams
        self.hp_stream_idx = Set(initialize=self.hp_stages.index_set() * [1, 2])
        self.ip_stream_idx = Set(initialize=self.ip_stages.index_set() * [1, 2])
        self.lp_stream_idx = Set(initialize=self.lp_stages.index_set() * [1, 2])

        # Throw out unneeded streams
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

        @self.Expression(self.flowsheet().config.time)
        def power(b, t):
            return (
                sum(b.inlet_stage[i].power_shaft[t] for i in b.inlet_stage)
                + b.outlet_stage.power_shaft[t]
                + sum(b.hp_stages[i].power_shaft[t] for i in b.hp_stages)
                + sum(b.ip_stages[i].power_shaft[t] for i in b.ip_stages)
                + sum(b.lp_stages[i].power_shaft[t] for i in b.lp_stages)
            )

        # Connect inlet stage to hp section
        #   not allowing disconnection of inlet and first regular hp stage
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
        TransformationFactory("network.expand_arcs").apply_to(self)

    def _split_cfg(self, unit_cfg, no=2):
        """
        This creates a configuration dictionary for a splitter.

        Args:
            unit_cfg: The base unit config dict.
            no: Number of outlets, default=2
        """
        # Create a dict for splitter config args
        s_cfg = copy.copy(unit_cfg)  # splitter config based on unit_cfg
        s_cfg.update(
            split_basis=SplittingType.totalFlow,
            ideal_separation=False,
            num_outlets=no,
            energy_split_basis=EnergySplittingType.equal_molar_enthalpy,
        )
        return s_cfg

    def _mix_cfg(self, unit_cfg, ni=2):
        """
        This creates a configuration dictionary for a mixer.

        Args:
            unit_cfg: The base unit config dict.
            ni: Number of inlets, default=2
        """
        m_cfg = copy.copy(unit_cfg)  # splitter config based on unit_cfg
        m_cfg.update(
            num_inlets=ni, momentum_mixing_type=MomentumMixingType.minimize_and_equality
        )
        return m_cfg

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

    def turbine_outlet_cf_fix(self, value):
        """
        Fix the inlet turbine stage flow coefficient.  These are
        generally the same for each of the parallel stages so this provides
        a convenient way to set them.

        Args:
            value: The value to fix the turbine inlet flow coefficients at
        """
        self.outlet_stage.flow_coeff.fix(value)

    def initialize(
        self,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6, "max_iter": 35},
        copy_disconneted_flow=True,
    ):
        """
        Initialize
        """
        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but initializes
        #   the values of free vars
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)
        ni = self.config.num_parallel_inlet_stages
        flow_guess = self.inlet_split.inlet.flow_mol[0].value

        def init_section(stages, splits, disconnects, prev_port):
            if 0 in splits:
                _set_port(splits[0].inlet, prev_port)
                splits[0].initialize(outlvl=outlvl, solver=solver, optarg=optarg)
                prev_port = splits[0].outlet_1
            for i in stages:
                if i - 1 not in disconnects:
                    _set_port(stages[i].inlet, prev_port)
                else:
                    if copy_disconneted_flow:
                        for t in stages[i].stages[i].inlet.flow_mol[t]:
                            stages[i].inlet.flow_mol[t] = prev_port.flow_mol[t]
                stages[i].initialize(outlvl=outlvl, solver=solver, optarg=optarg)
                prev_port = stages[i].outlet
                if i in splits:
                    _set_port(splits[i].inlet, prev_port)
                    splits[i].initialize(outlvl=outlvl, solver=solver, optarg=optarg)
                    prev_port = splits[i].outlet_1
            return prev_port

        for k in [1, 2]:
            # Initialize Splitter
            # Fix n - 1 split fractions
            self.inlet_split.split_fraction[0, "outlet_1"].value = 1.0 / ni
            for i in self.inlet_stage_idx:
                if i == 1:  # fix rest of splits at leaving first one free
                    continue
                self.inlet_split.split_fraction[0, "outlet_{}".format(i)].fix(1.0 / ni)
            # fix inlet and free outlet
            self.inlet_split.inlet.fix()
            for i in self.inlet_stage_idx:
                ol = getattr(self.inlet_split, "outlet_{}".format(i))
                ol.unfix()
            self.inlet_split.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
            # free split fractions
            for i in self.inlet_stage_idx:
                self.inlet_split.split_fraction[0, "outlet_{}".format(i)].unfix()

            # Initialize valves
            for i in self.inlet_stage_idx:
                _set_port(
                    self.throttle_valve[i].inlet,
                    getattr(self.inlet_split, "outlet_{}".format(i)),
                )
                self.throttle_valve[i].initialize(
                    outlvl=outlvl, solver=solver, optarg=optarg
                )

            # Initialize turbine
            for i in self.inlet_stage_idx:
                _set_port(self.inlet_stage[i].inlet, self.throttle_valve[i].outlet)
                self.inlet_stage[i].initialize(
                    outlvl=outlvl, solver=solver, optarg=optarg
                )

            # Initialize Mixer
            self.inlet_mix.use_minimum_inlet_pressure_constraint()
            for i in self.inlet_stage_idx:
                _set_port(
                    getattr(self.inlet_mix, "inlet_{}".format(i)),
                    self.inlet_stage[i].outlet,
                )
                getattr(self.inlet_mix, "inlet_{}".format(i)).fix()
            self.inlet_mix.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
            for i in self.inlet_stage_idx:
                getattr(self.inlet_mix, "inlet_{}".format(i)).unfix()
            self.inlet_mix.use_equal_pressure_constraint()

            prev_port = self.inlet_mix.outlet
            prev_port = init_section(
                self.hp_stages, self.hp_split, self.config.hp_disconnect, prev_port
            )
            if len(self.hp_stages) in self.config.hp_disconnect:
                prev_port = self.ip_stages[1].inlet
            prev_port = init_section(
                self.ip_stages, self.ip_split, self.config.ip_disconnect, prev_port
            )
            if len(self.ip_stages) in self.config.ip_disconnect:
                prev_port = self.lp_stages[1].inlet
            prev_port = init_section(
                self.lp_stages, self.lp_split, self.config.lp_disconnect, prev_port
            )

            _set_port(self.outlet_stage.inlet, prev_port)
            self.outlet_stage.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
            for t in self.flowsheet().time:
                self.inlet_split.inlet.flow_mol[
                    t
                ].value = self.outlet_stage.inlet.flow_mol[t].value

        from_json(self, sd=istate, wts=sp)
