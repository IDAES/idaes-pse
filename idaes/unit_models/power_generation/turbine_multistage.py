##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Multistage steam turbine for power generation.

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136, November
"""
import copy

from pyomo.environ import RangeSet, Set, TransformationFactory
from pyomo.network import Arc
from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.unit_models import Separator, Mixer, SplittingType
from idaes.unit_models.power_generation import (
    TurbineInletStage, TurbineStage, TurbineOutletStage)
from idaes.core.util.config import is_physical_parameter_block
from .turbine_multistage_config import _define_turbine_mutlistage_config


@declare_process_block_class("TurbineMultistage",
    doc="Multistage steam turbine with optional reheat and extraction")
class TurbineMultistageData(UnitModelBlockData):

    CONFIG = ConfigBlock()
    _define_turbine_mutlistage_config(CONFIG)

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
        # add inlet stage
        self._add_inlet_stage(unit_cfg)
        # add turbine stages.
        # inlet stage -> hp stages -> ip stages -> lp stages -> outlet stage

        self.hp_stages = TurbineStage(RangeSet(config.num_hp), default=unit_cfg)
        self.ip_stages = TurbineStage(RangeSet(config.num_ip), default=unit_cfg)
        self.lp_stages = TurbineStage(RangeSet(config.num_lp), default=unit_cfg)
        self.outlet_stage = TurbineOutletStage(default=unit_cfg)

        for i in self.hp_stages:
            self.hp_stages[i].ratioP[:].fix()
            self.hp_stages[i].efficiency_isentropic[:].fix()
        for i in self.ip_stages:
            self.ip_stages[i].ratioP[:].fix()
            self.ip_stages[i].efficiency_isentropic[:].fix()
        for i in self.lp_stages:
            self.lp_stages[i].ratioP[:].fix()
            self.lp_stages[i].efficiency_isentropic[:].fix()

        # Make the config args for each splitter
        # first the generic splitter config
        s_cfg = copy.copy(unit_cfg) # splitter config based on unit_cfg
        s_cfg.update(split_basis=SplittingType.totalFlow, ideal_separation=False)
        del s_cfg["has_holdup"]
        del s_cfg["has_phase_equilibrium"]
        s_cfg["num_outlets"] = 2
        # Then make one for each splitter
        hp_splt_cfg = {}
        for i in config.hp_split_locations:
            hp_splt_cfg[i] = copy.copy(s_cfg)
        ip_splt_cfg = {}
        for i in config.ip_split_locations:
            ip_splt_cfg[i] = copy.copy(s_cfg)
        lp_splt_cfg = {}
        for i in config.lp_split_locations:
            lp_splt_cfg[i] = copy.copy(s_cfg)
        # Now to finish up if there are more than two outlets, set that
        for i, v in config.hp_split_num_outlets.items():
            hp_splt_cfg[i]["num_outlets"] = v
        for i, v in config.ip_split_num_outlets.items():
            ip_splt_cfg[i]["num_outlets"] = v
        for i, v in config.lp_split_num_outlets.items():
            lp_splt_cfg[i]["num_outlets"] = v
        # put in splitters for turbine steam extractions
        if config.hp_split_locations:
            self.hp_split = Separator(
                config.hp_split_locations, initialize=hp_splt_cfg)
        if config.ip_split_locations:
            self.ip_split = Separator(
                config.ip_split_locations, initialize=ip_splt_cfg)
        if config.lp_split_locations:
            self.lp_split = Separator(
                config.lp_split_locations, initialize=lp_splt_cfg)

        # All the unit models are in place next step is to connect them
        #   arcs come up in this section and Pyomo arcs for connecting ports
        #   not related for turbine steam admission.

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
            sr = set()
            for i in index_set:
                if (i[0] in discon or i[0] == nstages) and i[0] in splits:
                    # don't connect stage i to next remove stream after split
                    sr.add((i[0], 2))
                elif (i[0] in discon or i[0] == nstages) and i[0] not in splits:
                    # no splitter and disconnect so remove both streams
                    sr.add((i[0], 1))
                    sr.add((i[0], 2))
                elif i[0] not in splits:
                    # no splitter and not disconnected so just second steam
                    sr.add((i[0], 2))
                else:
                    # has splitter so need both streams don't remove anything
                    pass
            for i in sr:
                index_set.remove(i)

        def _arc_rule(turbines, splitters):
            """
            This creates a rule function for arcs in a tubrine section. when
            this is used the indexes for nonexitant stream will have already
            been removed, so any indexes the rule will get should have a stream
            accossiated.

            Args:
                turbines (TurbinStage): Indexed block with turbine section stages
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

        # Connect hp section to ip section unless its a disconnect locations
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
        # Connect ip section to lp section unless its a disconnect locations
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

    def _add_inlet_stage(self, unit_cfg):
        """
        Add parallel turbine inlet stage models to simulate partial arc
        addmission. This includes a splitter mixer and a control valve before
        each stage. (TODO add throttel valves once a valve model is available
        from the basic framework models.)

        Args:
            unit_cfg: The base unit config dict.
        """
        ni = self.config.num_parallel_inlet_stages

        # Splitter config
        s_cfg = copy.copy(unit_cfg) # splitter config based on unit_cfg
        s_cfg.update(split_basis=SplittingType.totalFlow,
            ideal_separation=False, num_outlets=ni)
        del s_cfg["has_holdup"]
        del s_cfg["has_phase_equilibrium"]
        # add splitter
        self.inlet_split = Separator(default=s_cfg)

        # Mixer config
        m_cfg = copy.copy(unit_cfg) # splitter config based on unit_cfg
        m_cfg.update(num_inlets=ni)
        del m_cfg["has_holdup"]
        del m_cfg["has_phase_equilibrium"]
        # Add mixer
        self.inlet_mix = Mixer(default=m_cfg)

        #TODO<jce> Add throttel valves for each inlet stage
        # Add turbine stages
        self.inlet_stage = TurbineInletStage(RangeSet(ni), default=unit_cfg)

        # Add connections
        # TODO<jce> when framework is updated fix indexing
        def _split_to_rule(b, i):
            return {"source":getattr(self.inlet_split, "outlet_{}".format(i)),
                "destination":self.inlet_stage[i].inlet}
        def _inlet_to_rule(b, i):
            return {"source":self.inlet_stage[i].outlet,
                "destination":getattr(self.inlet_mix, "inlet_{}".format(i))}

        self.split_to_inlet_stage_stream = Arc(RangeSet(ni), rule=_split_to_rule)
        self.inlet_stage_to_mix = Arc(RangeSet(ni), rule=_inlet_to_rule)
