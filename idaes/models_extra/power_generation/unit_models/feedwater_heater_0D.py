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
This file contains 0D feedwater heater models. These models are suitable for
steady state calculations. For dynamic modeling 1D models are required. There
are two models included here.

1) FWHCondensing0D: this is a regular 0D heat exchanger model with a constraint
   added to ensure all the steam fed to the feedwater heater is condensed at
   the outlet. At the shell outlet the molar enthalpy is equal to the the
   saturated liquid molar enthalpy.
2) FWH0D is a feedwater heater model with three sections and a mixer for
   combining another feedwater heater's drain outlet with steam extracted from
   the turbine.  The drain mixer, desuperheat, and drain cooling sections are
   optional.  Only the condensing section is required.
"""

__author__ = "John Eslick"

from pyomo.common.config import ConfigValue, ConfigBlock, Bool
from pyomo.environ import TransformationFactory, Var, value
from pyomo.network import Arc

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    MaterialBalanceType,
)
from idaes.models.unit_models.heat_exchanger import HeatExchangerData
from idaes.models.unit_models import Mixer, MomentumMixingType, HeatExchanger
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import useDefault
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def _define_feedwater_heater_0D_config(config):
    config.declare(
        "has_drain_mixer",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Add a mixer to the inlet of the condensing section",
            doc="""Add a mixer to the inlet of the condensing section to add
water from the drain of another feedwaterheater to the steam, if True""",
        ),
    )
    config.declare(
        "has_desuperheat",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Add a mixer desuperheat section to the heat exchanger",
            doc="Add a mixer desuperheat section to the heat exchanger",
        ),
    )
    config.declare(
        "has_drain_cooling",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Add a section after condensing section cool condensate.",
            doc="Add a section after condensing section to cool condensate.",
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
    config.declare("condense", HeatExchangerData.CONFIG())
    config.declare("desuperheat", HeatExchangerData.CONFIG())
    config.declare("cooling", HeatExchangerData.CONFIG())


def _set_port(p1, p2):
    """
    Copy the values from port p2 to port p1.

    Args:
        p1: port to copy values to
        p2: port to compy values from
    """
    for k, v in p1.vars.items():
        if isinstance(v, Var):
            for i in v:
                v[i].value = value(p2.vars[k][i])


def _set_prop_pack(hxcfg, fwhcfg):
    """
    Set the property package and property pacakge args to the values given for
    the overall feedwater heater model if not otherwise specified.

    Args:
        hxcfg: Heat exchanger subblock config block
        fwhcfg: Overall feedwater heater config block
    """
    # this sets the property pack for the hot and cold side, but if the user
    # provides a specific property package using the tube and shell names it
    # will override this.  I think this behavior is fine, and what we'd want.
    if hxcfg.hot_side.property_package == useDefault:
        hxcfg.hot_side.property_package = fwhcfg.property_package
        hxcfg.hot_side.property_package_args = fwhcfg.property_package_args
    if hxcfg.cold_side.property_package == useDefault:
        hxcfg.cold_side.property_package = fwhcfg.property_package
        hxcfg.cold_side.property_package_args = fwhcfg.property_package_args


@declare_process_block_class(
    "FWHCondensing0D",
    doc="""Feedwater Heater Condensing Section
The feedwater heater condensing section model is a normal 0D heat exchanger
model with an added constraint to calculate the steam flow such that the outlet
of shell is a saturated liquid.""",
)
class FWHCondensing0DData(HeatExchangerData):
    def build(self):
        super().build()
        units_meta = self.hot_side.config.property_package.get_metadata()
        self.enth_sub = Var(
            self.flowsheet().time,
            initialize=0,
            units=units_meta.get_derived_units("energy_mole"),
        )
        self.enth_sub.fix()

        @self.Constraint(
            self.flowsheet().time,
            doc="Calculate steam extraction rate such that all steam condenses",
        )
        def extraction_rate_constraint(b, t):
            return (
                b.hot_side.properties_out[t].enth_mol - b.enth_sub[t]
                == b.hot_side.properties_out[t].enth_mol_sat_phase["Liq"]
            )

    def initialize_build(self, *args, **kwargs):
        """
        Use the regular heat exchanger initialization, with the extraction rate
        constraint deactivated; then it activates the constraint and calculates
        a steam inlet flow rate.
        """
        solver = kwargs.get("solver", None)
        optarg = kwargs.get("oparg", {})
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        self.extraction_rate_constraint.deactivate()
        self.area.fix()
        self.overall_heat_transfer_coefficient.fix()
        self.hot_side_inlet.fix()
        self.cold_side_inlet.fix()
        self.hot_side_outlet.unfix()
        self.cold_side_outlet.unfix()

        # Do regular heat exchanger intialization
        super().initialize_build(*args, **kwargs)
        self.extraction_rate_constraint.activate()
        self.hot_side_inlet.flow_mol.unfix()

        # Create solver
        opt = get_solver(solver, optarg)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            "Initialization Complete (w/ extraction calc): {}".format(
                idaeslog.condition(res)
            )
        )

        from_json(self, sd=istate, wts=sp)


@declare_process_block_class(
    "FWH0D",
    doc="""Feedwater Heater Model
This is a 0D feedwater heater model.  The model may contain three 0D heat
exchanger models representing the desuperheat, condensing and drain cooling
sections of the feedwater heater. Only the condensing section must be included.
A drain mixer can also be optionally included, which mixes the drain outlet of
another feedwater heater with the steam fed into the condensing section.
""",
)
class FWH0DData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    _define_feedwater_heater_0D_config(CONFIG)

    def build(self):
        super().build()
        config = self.config  # sorter ref to config for less line splitting

        # All feedwater heaters have a condensing section
        _set_prop_pack(config.condense, config)
        self.condense = FWHCondensing0D(**config.condense)

        # Add a mixer to add the drain stream from another feedwater heater
        if config.has_drain_mixer:
            mix_cfg = {  # general unit model config
                "dynamic": config.dynamic,
                "has_holdup": config.has_holdup,
                "property_package": config.property_package,
                "property_package_args": config.property_package_args,
                "momentum_mixing_type": MomentumMixingType.none,
                "material_balance_type": MaterialBalanceType.componentTotal,
                "inlet_list": ["steam", "drain"],
            }
            self.drain_mix = Mixer(**mix_cfg)

            @self.drain_mix.Constraint(self.drain_mix.flowsheet().time)
            def mixer_pressure_constraint(b, t):
                """
                Constraint to set the drain mixer pressure to the pressure of
                the steam extracted from the turbine. The drain inlet should
                always be a higher pressure than the steam inlet.
                """
                return b.steam_state[t].pressure == b.mixed_state[t].pressure

            # Connect the mixer to the condensing section inlet
            self.SMX = Arc(
                source=self.drain_mix.outlet, destination=self.condense.hot_side_inlet
            )

        # Add a desuperheat section before the condensing section
        if config.has_desuperheat:
            _set_prop_pack(config.desuperheat, config)
            self.desuperheat = HeatExchanger(**config.desuperheat)
            # set default area less than condensing section area, this will
            # almost always be overridden by the user fixing an area later
            self.desuperheat.area.value = 10
            if config.has_drain_mixer:
                self.SDS = Arc(
                    source=self.desuperheat.hot_side_outlet,
                    destination=self.drain_mix.steam,
                )
            else:
                self.SDS = Arc(
                    source=self.desuperheat.hot_side_outlet,
                    destination=self.condense.hot_side_inlet,
                )
            self.FW2 = Arc(
                source=self.condense.cold_side_outlet,
                destination=self.desuperheat.cold_side_inlet,
            )

        # Add a drain cooling section after the condensing section
        if config.has_drain_cooling:
            _set_prop_pack(config.cooling, config)
            self.cooling = HeatExchanger(**config.cooling)
            # set default area less than condensing section area, this will
            # almost always be overridden by the user fixing an area later
            self.cooling.area.value = 10
            self.FW1 = Arc(
                source=self.cooling.cold_side_outlet,
                destination=self.condense.cold_side_inlet,
            )
            self.SC = Arc(
                source=self.condense.hot_side_outlet,
                destination=self.cooling.hot_side_inlet,
            )

        TransformationFactory("network.expand_arcs").apply_to(self)

    def initialize_build(self, *args, **kwargs):
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        config = self.config  # shorter ref to config for less line splitting
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # the initialization here isn't straight forward since the heat
        # exchanger may have 3 stages and they are countercurrent.  For
        # simplicity each stage in initialized with the same cooling water
        # inlet conditions then the whole feedwater heater is solved together.
        # There are more robust approaches which can be implimented if the
        # need arises.

        # initialize desuperheat if include
        if config.has_desuperheat:
            if config.has_drain_cooling:
                _set_port(
                    self.desuperheat.cold_side_inlet, self.cooling.cold_side_inlet
                )
            else:
                _set_port(
                    self.desuperheat.cold_side_inlet, self.condense.cold_side_inlet
                )
            self.desuperheat.initialize(*args, **kwargs)
            self.desuperheat.hot_side_inlet.flow_mol.unfix()
            if config.has_drain_mixer:
                _set_port(self.drain_mix.steam, self.desuperheat.hot_side_outlet)
            else:
                _set_port(
                    self.condense.hot_side_inlet, self.desuperheat.hot_side_outlet
                )
            # fix the steam and fwh inlet for init
            self.desuperheat.hot_side_inlet.fix()
            self.desuperheat.hot_side_inlet.flow_mol.unfix()  # unfix for extract calc
        # initialize mixer if included
        if config.has_drain_mixer:
            self.drain_mix.steam.fix()
            self.drain_mix.drain.fix()
            self.drain_mix.outlet.unfix()
            self.drain_mix.initialize(*args, **kwargs)
            _set_port(self.condense.hot_side_inlet, self.drain_mix.outlet)
            if config.has_desuperheat:
                self.drain_mix.steam.unfix()
            else:
                self.drain_mix.steam.flow_mol.unfix()
        # Initialize condense section
        if config.has_drain_cooling:
            _set_port(self.condense.cold_side_inlet, self.cooling.cold_side_inlet)
            self.cooling.cold_side_inlet.fix()
        else:
            self.condense.cold_side_inlet.fix()
        if not config.has_drain_mixer and not config.has_desuperheat:
            self.condense.hot_side_inlet.fix()
            self.condense.hot_side_inlet.flow_mol.unfix()

        tempsat = value(self.condense.hot_side.properties_in[0].temperature_sat)
        temp = value(self.condense.cold_side.properties_in[0].temperature)
        if tempsat - temp < 30:
            init_log.warning(
                "The steam sat. temperature ({}) is near the feedwater"
                " inlet temperature ({})".format(tempsat, temp)
            )

        self.condense.initialize(*args, **kwargs)
        # Initialize drain cooling if included
        if config.has_drain_cooling:
            _set_port(self.cooling.hot_side_inlet, self.condense.hot_side_outlet)
            self.cooling.initialize(*args, **kwargs)

        # Solve all together
        opt = get_solver(kwargs.get("solver"), kwargs.get("oparg", {}))

        assert degrees_of_freedom(self) == 0
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            "Condensing hot side inlet delta T = {}".format(
                value(self.condense.delta_temperature_in[0])
            )
        )
        init_log.info(
            "Condensing hot side outlet delta T = {}".format(
                value(self.condense.delta_temperature_out[0])
            )
        )
        init_log.info(
            "Steam Flow = {}".format(value(self.condense.hot_side_inlet.flow_mol[0]))
        )
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        from_json(self, sd=istate, wts=sp)
