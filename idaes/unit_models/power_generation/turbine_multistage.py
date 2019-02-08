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

from pyomo.environ import RangeSet
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
        unit_cfg = {
            "dynamic":config.dynamic,
            "has_holdup":config.has_holdup,
            "has_phase_equilibrium":config.has_phase_equilibrium,
            "property_package":config.property_package,
            "property_package_args":config.property_package_args,
        }

        self.inlet_stage = TurbineInletStage(default=unit_cfg)
        self.hp_stages = TurbineStage(RangeSet(config.num_hp), default=unit_cfg)
        self.ip_stages = TurbineStage(RangeSet(config.num_ip), default=unit_cfg)
        self.lp_stages = TurbineStage(RangeSet(config.num_lp), default=unit_cfg)
        self.outlet_stage = TurbineOutletStage(default=unit_cfg)

        s_cfg = copy.copy(unit_cfg)
        s_cfg.update(split_basis=SplittingType.totalFlow, ideal_separation=False)
        del s_cfg["has_holdup"]
        del s_cfg["has_phase_equilibrium"]
        if config.hp_split_locations:
            self.hp_split = Separator(config.hp_split_locations, default=s_cfg)
        if config.ip_split_locations:
            self.ip_split = Separator(config.ip_split_locations, default=s_cfg)
        if config.lp_split_locations:
            self.lp_split = Separator(config.lp_split_locations, default=s_cfg)

        m_cfg = copy.copy(unit_cfg)
        del m_cfg["has_holdup"]
        del m_cfg["has_phase_equilibrium"]
        if config.hp_mix_locations:
            self.hp_mix = Mixer(config.hp_mix_locations, default=m_cfg)
        if config.ip_mix_locations:
            self.ip_mix = Mixer(config.ip_mix_locations, default=m_cfg)
        if config.lp_mix_locations:
            self.lp_mix = Mixer(config.lp_mix_locations, default=m_cfg)
