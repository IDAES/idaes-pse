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

from pyomo.environ import RangeSet, Set
from pyomo.network import Arc
from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.unit_models import Separator, Mixer, SplittingType
from idaes.unit_models.power_generation import (
    TurbineInletStage, TurbineStage, TurbineOutletStage)
from idaes.core import UnitBlockData
from idaes.core.util.config import is_physical_parameter_block
from .turbine_multistage_config import (
    _define_turbine_mutlistage_config, _ReheatType)

from idaes.core import declare_process_block_class, UnitBlockData


@declare_process_block_class("TurbineMultistage",
    doc="Multistage steam turbine with optional reheat and extraction")
class TurbineMultistageData(UnitBlockData):

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

        # add turbine stages.
        # inlet stage -> hp stages -> ip stages -> lp stages -> outlet stage
        self.inlet_stage = TurbineInletStage(default=unit_cfg)
        self.hp_stages = TurbineStage(RangeSet(config.num_hp), default=unit_cfg)
        self.ip_stages = TurbineStage(RangeSet(config.num_ip), default=unit_cfg)
        self.lp_stages = TurbineStage(RangeSet(config.num_lp), default=unit_cfg)
        self.outlet_stage = TurbineOutletStage(default=unit_cfg)

        # put in splitters for feedwater heater, reheat, ...
        s_cfg = copy.copy(unit_cfg) # splitter config based on unit_cfg
        s_cfg.update(split_basis=SplittingType.totalFlow, ideal_separation=False)
        del s_cfg["has_holdup"]
        del s_cfg["has_phase_equilibrium"]
        s_cfg["num_outlets"] = 2
        if config.hp_split_locations:
            self.hp_split = Separator(config.hp_split_locations, default=s_cfg)
        if config.ip_split_locations:
            self.ip_split = Separator(config.ip_split_locations, default=s_cfg)
        if config.lp_split_locations:
            self.lp_split = Separator(config.lp_split_locations, default=s_cfg)

        self.hp_stream_idx = Set(initialize=self.hp_stages.index_set()*[1,2])
        self.ip_stream_idx = Set(initialize=self.ip_stages.index_set()*[1,2])
        self.lp_stream_idx = Set(initialize=self.lp_stages.index_set()*[1,2])

        # okay I'm gonna talk about arcs here. They are the Pyomo components
        # used to make the material streams, this has nothing to do with the
        # turbine, an is totally unrealated to partial arc admission
        def _arc_indexes(nstages, index_set, discon, splits):
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
                    # no splitter and not disconnected so just first steam
                    sr.add((i[0], 2))
                else:
                    # has splitter so need both streams don't remove anything
                    pass
            for i in sr:
                index_set.remove(i)

        def _arc_rule(split_loc, turbines, splitters):
            # already removed indexes for unneeded or disconnected streams, so
            # I don't need to worry about that here
            #TODO<jce> Need to track the port indexing which is going to change
            #          dynamics are broken so just hooking time zero, after
            #          dynamics are fixed in framwwork the port index should go
            #          away.
            def _rule(b, i, j):
                if i in splitters and j == 1:
                    return {"source":turbines[i].outlet[0],
                            "destination":splitters[i].inlet[0]}
                elif j == 2:
                    return {"source":splitters[i].outlet_1[0],
                            "destination":turbines[i+1].inlet[0]}
                else:
                    return {"source":turbines[i].outlet[0],
                            "destination":turbines[i+1].inlet[0]}
            return _rule

        _arc_indexes(config.num_hp,
            self.hp_stream_idx, config.hp_disconnect, config.hp_split_locations)
        _arc_indexes(config.num_ip,
            self.ip_stream_idx, config.ip_disconnect, config.ip_split_locations)
        _arc_indexes(config.num_lp,
            self.lp_stream_idx, config.lp_disconnect, config.lp_split_locations)

        self.hp_stream = Arc(self.hp_stream_idx, rule=
            _arc_rule(config.hp_split_locations, self.hp_stages, self.hp_split))
        self.ip_stream = Arc(self.ip_stream_idx, rule=
            _arc_rule(config.hp_split_locations, self.ip_stages, self.ip_split))
        self.lp_stream = Arc(self.lp_stream_idx, rule=
            _arc_rule(config.hp_split_locations, self.lp_stages, self.lp_split))

        if 0 not in conf.ip_disconnect and conf.num_hp not in conf.hp_disconnect:
            #connect hp section to the lp section
            last_hp = self.hp_stages[-1]
            if last_hp in config.kp_split_locations: # connect splitter to ip
                self.hp_to_ip_stream = Arc(
                    
                )
            else: # connect last hp to ip
                self.hp_to_ip_stream = Arc()
        if 0 not in conf.lp_disconnect and conf.num_ip not in conf.ip_disconnect:
            #connect ip section to the lp section
            last_ip = self.ip_stages[-1]

            self.ip_to_lp_stream =
            self.ip_to_lp_stream =
        # you are not allowed to stick stuff between the inlet stage and hp
        # or the last lp stage and the outlet stage.  So connect those
        last_lp = self.lp_stages[-1]
        if last_lp in config.lp_split_locations: # connect splitter to outlet
            self.lp_to_outlet_stream = Arc(
                source=self.lp_split[last_lp].outlet_1[0],
                destination=self.outlet_stage.inlet[0])
        else: # connect last lpstage to outlet
            self.lp_to_outlet_stream = Arc(
                source=self.lp_stages[last_lp].outlet[0],
                destination=self.outlet_stage.inlet[0])
