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
This file contains 0D feedwater heater models. These models are suitable for
steady state calculations. For dynamic modeling 1D models are required.
"""

__author__ = "John Eslick"

import logging
from pyomo.environ import SolverFactory, TransformationFactory, Var, value
from pyomo.opt import TerminationCondition
from pyomo.network import Arc

from idaes.core import declare_process_block_class, UnitModelBlockData
from idaes.unit_models.heat_exchanger import HeatExchangerData
from idaes.unit_models import Mixer, MomentumMixingType, HeatExchanger
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.ui.report import degrees_of_freedom
from .feedwater_heater_0D_config import _define_feedwater_heater_0D_config

_log = logging.getLogger(__name__)

def _set_port(p1, p2):
    for k, v in p1.vars.items():
        if isinstance(v, Var):
            for i in v:
                v[i].value = value(p2.vars[k][i])


@declare_process_block_class("FWHCondensing0D",
    doc="Feedwater heater condensing section")
class FWHCondensing0DData(HeatExchangerData):
    def build(self):
        super().build()
        @self.Constraint(self.time_ref,
            doc="Calculate steam extraction rate such that all steam condenses")
        def extraction_rate_constraint(b, t):
            return b.side_1.properties_out[t].enth_mol_sat_phase["Liq"] == \
                   b.side_1.properties_out[t].enth_mol

    def initialize(self, *args, **kwargs):
        """
        Use the regular heat exchanger initilization, with the extraction rate
        constraint deactivated; then it activates the constraint and calculates
        a steam inlet flow rate.
        """
        self.extraction_rate_constraint.deactivate()
        super().initialize(*args, **kwargs)
        self.extraction_rate_constraint.activate()

        solver = kwargs.get("solver", "ipopt")
        optarg = kwargs.get("oparg", {})
        outlvl = kwargs.get("outlvl", 0)

        opt = SolverFactory(solver)
        opt.options = optarg
        tee = True if outlvl >= 3 else False
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        self.area.fix()
        self.overall_heat_transfer_coefficient.fix()
        self.inlet_1.fix()
        self.inlet_2.fix()
        self.outlet_1.unfix()
        self.outlet_2.unfix()
        self.inlet_1.flow_mol.unfix()
        results = opt.solve(self, tee=tee)

        if results.solver.termination_condition == TerminationCondition.optimal:
            if outlvl >= 2:
                _log.info('{} Initialization Failed.'.format(self.name))
        else:
            _log.warning('{} Initialization Failed.'.format(self.name))

        from_json(self, sd=istate, wts=sp)


@declare_process_block_class("FWH0D", doc="Feedwater heater")
class FWH0DData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    _define_feedwater_heater_0D_config(CONFIG)

    def build(self):
        super().build()
        config = self.config

        # All feedwater heaters have a condensing section
        self.condense = FWHCondensing0D(default=config.condense)

        # Add a mixer to add the drain stream from another feedwater heater
        if config.has_drain_mixer:
            mix_cfg = { # general unit model config
                "dynamic":config.dynamic,
                "has_holdup":config.has_holdup,
                "property_package":config.property_package,
                "property_package_args":config.property_package_args,
                "momentum_mixing_type":MomentumMixingType.none,
                "inlet_list":["steam", "drain"]}
            self.drain_mix = Mixer(default=mix_cfg)
            @self.drain_mix.Constraint(self.drain_mix.time_ref)
            def mixer_pressure_constraint(b, t):
                return b.steam_state[t].pressure == b.mixed_state[t].pressure
            self.mix_to_cond = Arc(source=self.drain_mix.outlet,
                                   destination=self.condense.inlet_1)

        # Add a desuperheat section before the condensing section
        if config.has_desuperheat:
            self.desuperheat = HeatExchanger(default=config.desuperheat)
            self.desuperheat.area.value = 10
            if config.has_drain_mixer:
                self.desuperheat_drain_arc = Arc(
                    source=self.desuperheat.outlet_1,
                    destination=self.drain_mix.steam)
            else:
                self.desuperheat_drain_arc = Arc(
                    source=self.desuperheat.outlet_1,
                    destination=self.condense.inlet_1)
            self.desuperheat_fw_arc = Arc(
                source=self.desuperheat.outlet_2,
                destination=self.condense.inlet_2)

        # Add a drain cooling section after the condensing section
        if config.has_drain_cooling:
            self.cooling = HeatExchanger(default=config.cooling)
            self.cooling.area.value = 10

        TransformationFactory("network.expand_arcs").apply_to(self)

    def initialize(self, *args, **kwargs):
        config = self.config
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)
        # initialize desuperheat if include
        if config.has_desuperheat:
            self.desuperheat.initialize(*args, **kwargs)
            self.desuperheat.inlet_1.flow_mol.unfix()
            if config.has_drain_mixer:
                _set_port(self.drain_mix.steam, self.desuperheat.outlet_1)
            else:
                _set_port(self.condense.inlet_1, self.desuperheat.outlet_1)
            _set_port(self.condense.inlet_2, self.desuperheat.outlet_2)
        # initialize mixer if included
        if config.has_drain_mixer:
            self.drain_mix.steam.fix()
            self.drain_mix.drain.fix()
            self.drain_mix.outlet.unfix()
            self.drain_mix.initialize(*args, **kwargs)
            _set_port(self.condense.inlet_1, self.drain_mix.outlet)
            if config.has_desuperheat:
                self.drain_mix.steam.unfix()
            else:
                self.drain_mix.steam.flow_mol.unfix()
        # Initialize condense section
        self.condense.initialize(*args, **kwargs)
        # Initialize drain cooling if included
        if config.has_drain_cooling:
            pass

        # Solve all together
        outlvl = kwargs.get("outlvl", 0)
        opt = SolverFactory(kwargs.get("solver", "ipopt"))
        opt.options = kwargs.get("oparg", {})
        tee = True if outlvl >= 3 else False
        assert(degrees_of_freedom(self)==0)
        results = opt.solve(self, tee=tee)
        if results.solver.termination_condition == TerminationCondition.optimal:
            if outlvl >= 2:
                _log.info('{} Initialization Failed.'.format(self.name))
        else:
            _log.warning('{} Initialization Failed.'.format(self.name))

        from_json(self, sd=istate, wts=sp)
