##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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

from pyomo.environ import RangeSet, Set, TransformationFactory, Var, value
from pyomo.network import Arc
from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.unit_models import (Separator, Mixer, SplittingType,
    EnergySplittingType, MomentumMixingType)
from idaes.unit_models.power_generation import (
    TurbineInletStage, TurbineStage, TurbineOutletStage, SteamValve)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util import from_json, to_json, StoreSpec
from .turbine_multistage_config import _define_turbine_multistage_config
from idaes.ui.report import degrees_of_freedom

def _set_port(p1, p2):
    for k, v in p1.vars.items():
        if isinstance(v, Var):
            for i in v:
                v[i].value = value(p2.vars[k][i])

@declare_process_block_class("TurbineMultistage",
    doc="Multistage steam turbine with optional reheat and extraction")
class TurbineMultistageData(UnitModelBlockData):

    CONFIG = ConfigBlock()
    _define_turbine_multistage_config(CONFIG)

    def build(self):
        super(TurbineMultistageData, self).build()
        config = self.config
        unit_cfg = { # general unit model config
            "dynamic":config.dynamic,
            "has_holdup":config.has_holdup,
            "has_phase_equilibrium":config.has_phase_equilibrium,
            "property_package":config.property_package,
            "property_package_args":config.property_package_args,
        }
        ni = self.config.num_parallel_inlet_stages
        inlet_idx = self.inlet_stage_idx = RangeSet(ni)

        # Adding unit models
        #------------------------

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
            self.hp_split = Separator(config.hp_split_locations,
                default=s_sfg_default, initialize=hp_splt_cfg)
        if config.ip_split_locations:
            self.ip_split = Separator(config.ip_split_locations,
                default=s_sfg_default, initialize=ip_splt_cfg)
        if config.lp_split_locations:
            self.lp_split = Separator(config.lp_split_locations,
                default=s_sfg_default, initialize=lp_splt_cfg)

        # Done with unit models.  Adding Arcs (streams).
        #------------------------------------------------

        # First up add streams in the inlet section
        def _split_to_rule(b, i):
            return {"source":getattr(self.inlet_split, "outlet_{}".format(i)),
                    "destination":self.throttle_valve[i].inlet}
        def _valve_to_rule(b, i):
            return {"source":self.throttle_valve[i].outlet,
                    "destination":self.inlet_stage[i].inlet}
        def _inlet_to_rule(b, i):
            return {"source":self.inlet_stage[i].outlet,
                    "destination":getattr(self.inlet_mix, "inlet_{}".format(i))}

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
            sr = set() # set of things to remove from the Arc index set
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
            for i in sr: # remove the unneeded Arc indexes
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
                    return {"source":turbines[i].outlet,
                            "destination":splitters[i].inlet}
                elif j == 2:
                    return {"source":splitters[i].outlet_1,
                            "destination":turbines[i+1].inlet}
                else:
                    return {"source":turbines[i].outlet,
                            "destination":turbines[i+1].inlet}
            return _rule

        # Create initial arcs index sets with all possible streams
        self.hp_stream_idx = Set(initialize=self.hp_stages.index_set()*[1,2])
        self.ip_stream_idx = Set(initialize=self.ip_stages.index_set()*[1,2])
        self.lp_stream_idx = Set(initialize=self.lp_stages.index_set()*[1,2])

        # Throw out unneeded streams
        _arc_indexes(config.num_hp,
            self.hp_stream_idx, config.hp_disconnect, config.hp_split_locations)
        _arc_indexes(config.num_ip,
            self.ip_stream_idx, config.ip_disconnect, config.ip_split_locations)
        _arc_indexes(config.num_lp,
            self.lp_stream_idx, config.lp_disconnect, config.lp_split_locations)

        # Create connections internal to each turbine section (hp, ip, and lp)
        self.hp_stream = Arc(self.hp_stream_idx,
            rule=_arc_rule(self.hp_stages, self.hp_split))
        self.ip_stream = Arc(self.ip_stream_idx,
            rule=_arc_rule(self.ip_stages, self.ip_split))
        self.lp_stream = Arc(self.lp_stream_idx,
            rule=_arc_rule(self.lp_stages, self.lp_split))

        # Connect hp section to ip section unless its a disconnect location
        last_hp = config.num_hp
        if 0 not in config.ip_disconnect and last_hp not in config.hp_disconnect:
            if last_hp in config.hp_split_locations: # connect splitter to ip
                self.hp_to_ip_stream = Arc(
                    source=self.hp_split[last_hp].outlet_1,
                    destination=self.ip_stages[1].inlet)
            else: # connect last hp to ip
                self.hp_to_ip_stream = Arc(
                    source=self.hp_stages[last_hp].outlet,
                    destination=self.ip_stages[1].inlet)
        # Connect ip section to lp section unless its a disconnect location
        last_ip = config.num_ip
        if 0 not in config.lp_disconnect and last_ip not in config.ip_disconnect:
            if last_ip in config.ip_split_locations: # connect splitter to ip
                self.ip_to_lp_stream = Arc(
                    source=self.ip_split[last_ip].outlet_1,
                    destination=self.lp_stages[1].inlet)
            else: # connect last hp to ip
                self.ip_to_lp_stream = Arc(
                    source=self.ip_stages[last_ip].outlet,
                    destination=self.lp_stages[1].inlet)
        # Connect inlet stage to hp section
        #   not allowing disconnection of inlet and first regular hp stage
        if 0 in config.hp_split_locations:
            # connect inlet mix to splitter and splitter to hp section
            self.inlet_to_splitter_stream = Arc(
                source=self.inlet_mix.outlet,
                destination=self.hp_split[0].inlet)
            self.splitter_to_hp_stream = Arc(
                source=self.hp_split[0].outlet_1,
                destination=self.hp_stages[1].inlet)
        else: # connect mixer to first hp turbine stage
            self.inlet_to_hp_stream = Arc(
                source=self.inlet_mix.outlet,
                destination=self.hp_stages[1].inlet)

        # Connect inlet stage to hp section
        #   not allowing disconnection of inlet and first regular hp stage
        last_lp = config.num_lp
        if last_lp in config.lp_split_locations: # connect splitter to outlet
            self.lp_to_outlet_stream = Arc(
                source=self.lp_split[last_lp].outlet_1,
                destination=self.outlet_stage.inlet)
        else: # connect last lpstage to outlet
            self.lp_to_outlet_stream = Arc(
                source=self.lp_stages[last_lp].outlet,
                destination=self.outlet_stage.inlet)
        TransformationFactory("network.expand_arcs").apply_to(self)

    def _split_cfg(self, unit_cfg, no=2):
        """
        This creates a configuration dictionary for a splitter.

        Args:
            unit_cfg: The base unit config dict.
            no: Number of outlets, default=2
        """
        # Create a dict for splitter config args
        s_cfg = copy.copy(unit_cfg) # splitter config based on unit_cfg
        s_cfg.update(
            split_basis=SplittingType.totalFlow,
            ideal_separation=False,
            num_outlets=no,
            energy_split_basis=EnergySplittingType.equal_molar_enthalpy)
        del s_cfg["has_phase_equilibrium"]
        return s_cfg

    def _mix_cfg(self, unit_cfg, ni=2):
        """
        This creates a configuration dictionary for a mixer.

        Args:
            unit_cfg: The base unit config dict.
            ni: Number of inlets, default=2
        """
        m_cfg = copy.copy(unit_cfg) # splitter config based on unit_cfg
        m_cfg.update(
            num_inlets=ni,
            momentum_mixing_type=MomentumMixingType.minimize_and_equality)
        del m_cfg["has_phase_equilibrium"]
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
            self.throttle_valve.Cv.fix(value)

    def turbine_inlet_cf_fix(self, value):
        """
        Fix the inlet turbine stage flow coefficient.  These are
        generally the same for each of the parallel stages so this provides
        a convenient way to set them.

        Args:
            value: The value to fix the turbine inlet flow coefficients at
        """
        for i in self.throttle_valve:
            self.inlet_stage.flow_coeff.fix(value)

    def initialize(self, outlvl=0, solver='ipopt',
        optarg={'tol':1e-6, 'max_iter':35}):
        """
        Initialize
        """
        stee = True if outlvl >= 3 else False
        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but initializes
        #   the values of free vars
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)


        ni = self.config.num_parallel_inlet_stages

        # Initialize Splitter
        # Fix n - 1 split fractions
        self.inlet_split.split_fraction[0,"outlet_1"].value = 1.0/ni
        for i in self.inlet_stage_idx:
            if i == 1: #fix rest of splits at leaving first one free
                continue
            self.inlet_split.split_fraction[0, "outlet_{}".format(i)].fix(1.0/ni)
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
            _set_port(self.throttle_valve[i].inlet,
                      getattr(self.inlet_split, "outlet_{}".format(i)))
            self.throttle_valve[i].initialize(outlvl=outlvl, solver=solver,
                optarg=optarg)

        # Initialize turbine
        for i in self.inlet_stage_idx:
            _set_port(self.inlet_stage[i].inlet, self.throttle_valve[i].outlet)
            self.inlet_stage[i].initialize(outlvl=outlvl, solver=solver,
                optarg=optarg)

        # Initialize Mixer
        self.inlet_mix.use_minimum_inlet_pressure_constraint()
        for i in self.inlet_stage_idx:
            _set_port(getattr(self.inlet_mix, "inlet_{}".format(i)),
                      self.inlet_stage[i].outlet)
            getattr(self.inlet_mix, "inlet_{}".format(i)).fix()
        self.inlet_mix.initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        for i in self.inlet_stage_idx:
            getattr(self.inlet_mix, "inlet_{}".format(i)).unfix()
        self.inlet_mix.use_equal_pressure_constraint()

        def init_section(stages, splits, disconnects, prev_port):
            if 0 in splits:
                _set_port(splits[0].inlet, prev_port)
                splits[0].initialize(outlvl=outlvl, solver=solver, optarg=optarg)
                prev_port = splits[0].outlet_1
            for i in stages:
                if i - 1 not in disconnects:
                    _set_port(stages[i].inlet, prev_port)
                stages[i].initialize(
                    outlvl=outlvl, solver=solver, optarg=optarg)
                prev_port = stages[i].outlet
                if i in splits:
                    _set_port(splits[i].inlet, prev_port)
                    splits[i].initialize(
                        outlvl=outlvl, solver=solver, optarg=optarg)
                    prev_port = splits[i].outlet_1
            return prev_port

        prev_port = self.inlet_mix.outlet
        prev_port = init_section(
            self.hp_stages, self.hp_split, self.config.hp_disconnect, prev_port)
        if len(self.hp_stages) in self.config.hp_disconnect:
            prev_port = self.ip_stages[1].inlet
        prev_port = init_section(
            self.ip_stages, self.ip_split, self.config.ip_disconnect, prev_port)
        if len(self.ip_stages) in self.config.ip_disconnect:
            prev_port = self.lp_stages[1].inlet
        prev_port = init_section(
            self.lp_stages, self.lp_split, self.config.lp_disconnect, prev_port)

        _set_port(self.outlet_stage.inlet, prev_port)
        self.outlet_stage.initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        from_json(self, sd=istate, wts=sp)
