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
   staurated liquid molar enthalpy.
2) FWH0D is a feedwater heater model with three sections and a mixer for
   combining another feedwater heater's drain outlet with steam extracted from
   the turbine.  The drain mixer, desuperheat, and drain cooling sections are
   optional. Only the condensing section is required.
"""

__author__ = "John Eslick, Jinliang Ma"
from pyomo.common.config import ConfigValue, In, ConfigBlock
from pyomo.environ import (
    SolverFactory,
    TransformationFactory,
    Var,
    value,
    asin,
    cos
    )
from pyomo.network import Arc

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
)
from idaes.generic_models.unit_models.heat_exchanger import HeatExchangerData
from idaes.generic_models.unit_models import MomentumMixingType, HeatExchanger
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import useDefault
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog
from idaes.power_generation.unit_models.helm import \
    HelmNtuCondenserData as CondenserData
from idaes.core.util.constants import Constants as const
import idaes.core.util.scaling as iscale
from idaes.power_generation.unit_models.helm import HelmMixer as Mixer
from idaes.core.util import get_solver, copy_port_values as _set_port

_log = idaeslog.getLogger(__name__)


def _define_feedwater_heater_0D_config(config):
    config.declare(
        "has_drain_mixer",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Add a mixer to the inlet of the condensing section",
            doc="""Add a mixer to the inlet of the condensing section to add
water from the drain of another feedwaterheater to the steam, if True""",
        ),
    )
    config.declare(
        "has_desuperheat",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Add a desuperheat section to the heat exchanger",
            doc="Add a mixer desuperheat section to the heat exchanger",
        ),
    )
    config.declare(
        "has_drain_cooling",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Add a drain cooler section to the heat exchanger.",
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
    config.declare("condense", CondenserData.CONFIG())
    config.declare("desuperheat", HeatExchangerData.CONFIG())
    config.declare("cooling", HeatExchangerData.CONFIG())


def _set_prop_pack(hxcfg, fwhcfg):
    """
    Set the property package and property pacakge args to the values given for
    the overall feedwater heater model if not otherwise specified.

    Args:
        hxcfg: Heat exchanger subblock config block
        fwhcfg: Overall feedwater heater config block
    """
    # this method sets the property pack for the hot and cold side,
    # but if the user provides a specific property package using
    # the tube and shell names it will override this.
    if hxcfg.hot_side_config.property_package == useDefault:
        hxcfg.hot_side_config.property_package = fwhcfg.property_package
        hxcfg.hot_side_config.property_package_args = fwhcfg.\
            property_package_args
    if hxcfg.cold_side_config.property_package == useDefault:
        hxcfg.cold_side_config.property_package = fwhcfg.property_package
        hxcfg.cold_side_config.property_package_args = fwhcfg.\
            property_package_args


@declare_process_block_class(
    "FWHCondensing0D",
    doc="""Feedwater Heater Condensing Section
The feedwater heater condensing section model is a normal 0D heat exchanger
model with an added constraint to calculate the steam flow such that the outlet
of shell is a saturated liquid.""",
)
class FWHCondensing0DData(CondenserData):
    def build(self):
        super().build()
        self.vol_frac_shell = Var(initialize=0.7,
                                  doc="Volume fraction of shell side fluid")
        self.heater_diameter = Var(initialize=1,
                                   doc="Inside diameter of FWH tank")
        self.cond_sect_length = Var(initialize=10,
                                    doc="Length of condensing section")
        self.level = Var(self.flowsheet().time,
                         initialize=1,
                         doc="Water level in condensing section")

        # Radius of feed water heater tank
        @self.Expression(doc="Radius of feed water heater tank")
        def heater_radius(b):
            return b.heater_diameter/2

        # Expressure for the angle from the tank cylinder center to
        # the circumference point at water level between -pi/2 and pi/2
        @self.Expression(self.flowsheet().time,
                         doc="Angle of water level")
        def alpha(b, t):
            return asin((b.level[t]-b.heater_radius)/b.heater_radius)

        # Constraint to calculate the shell side liquid volume
        @self.Constraint(self.flowsheet().time,
                         doc="Calculate the volume of "
                         "shell side based on water level")
        def shell_volume_eqn(b, t):
            return b.shell.volume[t] == ((
                (b.alpha[t]+0.5*const.pi)*b.heater_radius**2
                + b.heater_radius*cos(b.alpha[t])
                * (b.level[t]-b.heater_radius))
                * b.cond_sect_length*b.vol_frac_shell)

        # Total pressure change equation
        @self.Constraint(self.flowsheet().time,
                         doc="Pressure drop")
        def pressure_change_total_eqn(b, t):
            return b.shell.deltaP[t] == const.acceleration_gravity*b.level[t]\
                * b.shell.properties_out[t].dens_mass_phase["Liq"]

    def initialize(self, *args, **kwargs):
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

        self.area.fix()
        self.overall_heat_transfer_coefficient.fix()
        self.inlet_1.fix()
        self.inlet_2.fix()
        self.outlet_1.unfix()
        self.outlet_2.unfix()

        # Do condenser initialization
        self.inlet_1.flow_mol.unfix()
        # fix volume and pressure drop since
        # the condenser initialization dosen't require them
        self.side_1.volume.fix(10)
        self.side_1.deltaP.fix(0)
        self.shell_volume_eqn.deactivate()
        self.pressure_change_total_eqn.deactivate()
        super().initialize()
        self.side_1.volume.unfix()
        self.side_1.deltaP.unfix()
        self.shell_volume_eqn.activate()
        self.pressure_change_total_eqn.activate()

        # Create solver
        opt = get_solver(solver, optarg)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            "Initialization Complete (w/ steam flow calc): {}".format(
                idaeslog.condition(res)
            )
        )

        from_json(self, sd=istate, wts=sp)


@declare_process_block_class(
    "FWH0DDynamic",
    doc="""Feedwater Heater Model
This is a 0D feedwater heater model.  The model may contain three 0D heat
exchanger models representing the desuperheat, condensing and drain cooling
sections of the feedwater heater. Only the condensing section must be included.
A drain mixer can also be optionally included, which mixes the drain outlet of
another feedwater heater with the steam fed into the condensing section.
""",
)
class FWH0DDynamicData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    _define_feedwater_heater_0D_config(CONFIG)

    def build(self):
        super().build()
        config = self.config  # sorter ref to config for less line splitting

        # All feedwater heaters have a condensing section
        _set_prop_pack(config.condense, config)
        self.condense = FWHCondensing0D(default=config.condense)

        # Add a mixer to add the drain stream from another feedwater heater
        if config.has_drain_mixer:
            mix_cfg = {"dynamic": False,
                       "inlet_list": ["steam", "drain"],
                       "property_package": config.property_package,
                       "momentum_mixing_type": MomentumMixingType.none,
                       }
            self.drain_mix = Mixer(default=mix_cfg)

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
                source=self.drain_mix.outlet, destination=self.condense.inlet_1
            )

        # Add a desuperheat section before the condensing section
        if config.has_desuperheat:
            _set_prop_pack(config.desuperheat, config)
            self.desuperheat = HeatExchanger(default=config.desuperheat)
            # set default area less than condensing section area, this will
            # almost always be overridden by the user fixing an area later
            self.desuperheat.area.value = 10
            if config.has_drain_mixer:
                self.SDS = Arc(
                    source=self.desuperheat.outlet_1,
                    destination=self.drain_mix.steam
                )
            else:
                self.SDS = Arc(
                    source=self.desuperheat.outlet_1,
                    destination=self.condense.inlet_1
                )
            self.FW2 = Arc(
                source=self.condense.outlet_2,
                destination=self.desuperheat.inlet_2
            )

        # Add a drain cooling section after the condensing section
        if config.has_drain_cooling:
            _set_prop_pack(config.cooling, config)
            self.cooling = HeatExchanger(default=config.cooling)
            # set default area less than condensing section area, this will
            # almost always be overridden by the user fixing an area later
            self.cooling.area.value = 10
            self.FW1 = Arc(
                source=self.cooling.outlet_2, destination=self.condense.inlet_2
            )
            self.SC = Arc(
                source=self.condense.outlet_1, destination=self.cooling.inlet_1
            )

        TransformationFactory("network.expand_arcs").apply_to(self)

    def set_initial_condition(self):
        # currently assume steady-state for desuperheater
        if self.config.dynamic is True:
            if self.condense.config.has_holdup is True:
                self.condense.tube.material_accumulation[:, :, :].value = 0
                self.condense.tube.energy_accumulation[:, :].value = 0
                self.condense.tube.material_accumulation[0, :, :].fix(0)
                self.condense.tube.energy_accumulation[0, :].fix(0)
                self.condense.shell.material_accumulation[:, :, :].value = 0
                self.condense.shell.energy_accumulation[:, :].value = 0
                self.condense.shell.material_accumulation[0, :, :].fix(0)
                self.condense.shell.energy_accumulation[0, :].fix(0)
            if self.config.has_drain_cooling is True:
                if self.cooling.config.has_holdup is True:
                    self.cooling.tube.material_accumulation[:, :, :].value = 0
                    self.cooling.tube.energy_accumulation[:, :].value = 0
                    self.cooling.tube.material_accumulation[0, :, :].fix(0)
                    self.cooling.tube.energy_accumulation[0, :].fix(0)
                    self.cooling.shell.material_accumulation[:, :, :].value = 0
                    self.cooling.shell.energy_accumulation[:, :].value = 0
                    self.cooling.shell.material_accumulation[0, :, :].fix(0)
                    self.cooling.shell.energy_accumulation[0, :].fix(0)

    def initialize(self, *args, **kwargs):
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        config = self.config  # shorter ref to config for less line splitting
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # the initialization here isn't straight forward since
        # the heat exchanger may have 3 stages and they are countercurrent.
        # For simplicity each stage in initialized with the same cooling water
        # inlet conditions then the whole feedwater heater is solved together.
        # There are more robust aproaches which can be implimented if needed.

        # initialize desuperheat if any
        if config.has_desuperheat:
            if config.has_drain_cooling:
                _set_port(self.desuperheat.inlet_2, self.cooling.inlet_2)
            else:
                _set_port(self.desuperheat.inlet_2, self.condense.inlet_2)
            self.desuperheat.initialize(*args, **kwargs)
            self.desuperheat.inlet_1.flow_mol.unfix()
            if config.has_drain_mixer:
                _set_port(self.drain_mix.steam, self.desuperheat.outlet_1)
            else:
                _set_port(self.condense.inlet_1, self.desuperheat.outlet_1)
            # fix the steam and fwh inlet for init
            self.desuperheat.inlet_1.fix()
            self.desuperheat.inlet_1.flow_mol.unfix()  # unfix for extract calc
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
        if config.has_drain_cooling:
            _set_port(self.condense.inlet_2, self.cooling.inlet_2)
            self.cooling.inlet_2.fix()
        else:
            self.condense.inlet_2.fix()
        if not config.has_drain_mixer and not config.has_desuperheat:
            self.condense.inlet_1.fix()
            self.condense.inlet_1.flow_mol.unfix()

        tempsat = value(self.condense.shell.properties_in[0].temperature_sat)
        temp = value(self.condense.tube.properties_in[0].temperature)
        if tempsat - temp < 30:
            init_log.warning(
                "The steam sat. temperature ({}) is near the feedwater"
                " inlet temperature ({})".format(tempsat, temp)
            )

        self.condense.initialize()
        # Initialize drain cooling if included
        if config.has_drain_cooling:
            _set_port(self.cooling.inlet_1, self.condense.outlet_1)
            self.cooling.initialize(*args, **kwargs)

        # Create solver
        opt = get_solver(kwargs.get("solver"), kwargs.get("oparg", {}))

        assert degrees_of_freedom(self) == 0
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            "Condensing shell inlet delta T = {}".format(
                value(self.condense.delta_temperature_in[0])
            )
        )
        init_log.info(
            "Condensing Shell outlet delta T = {}".format(
                value(self.condense.delta_temperature_out[0])
            )
        )
        init_log.info(
            "Steam Flow = {}".format(value(self.condense.inlet_1.flow_mol[0]))
        )
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

        from_json(self, sd=istate, wts=sp)

    def calculate_scaling_factors(self):
        if hasattr(self, "mixer_pressure_constraint"):
            for t, c in self.mixer_pressure_constraint.items():
                sf = iscale.get_scaling_factor(
                    self.steam_state[t].pressure, default=1, warning=True)
                iscale.constraint_scaling_transform(c, sf, overwrite=False)
        if hasattr(self, "pressure_change_total_eqn"):
            for t, c in self.pressure_change_total_eqn.items():
                sf = iscale.get_scaling_factor(
                    self.steam_state[t].pressure, default=1, warning=True)
                iscale.constraint_scaling_transform(c, sf, overwrite=False)
