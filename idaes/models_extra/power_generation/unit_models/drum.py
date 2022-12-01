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
Drum model
The drum model consists of three main unit operations
1) a flash model (only for water and steam)
2) a mixer model
3) a tank drum level/accumulation model

Inlet Ports:

* water/steam mixture from water wall
* feedwater_inlet: feedwater from economizer/pipe

Outlet Ports:

* liquid_outlet: liquid to downcomer
* steam_outlet: steam leaving the drum

The drum model receives saturated water from the water wall, and this is
separated at the flash unit, the steam leaves the drum, while the water/liquid
mixes with the feedwater from the economizer (or water pipe model).
Finally, the mixed state liquid (stream leaving the mixer), is used as the
control volume of the drum level model, which computes velocity and pressure
drop. The exit of the drum model is the liquid outlet.

main assumptions:
1.- Heat loss is a variable given by the user (zero heat loss can
be specified if adiabatic)
2.- Calculate pressure change due to gravity based on water level
and contraction to downcomer
3.- Water level is either fixed for steady-state model or calculated for
dynamic model
4.- Assume enthalpy_in == enthalpy_out + heat loss
5.- Subcooled water from economizer and saturated water from waterwall
are mixed before entering drum

Created: August 19 2020
"""
# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
import pyomo.environ as pyo

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
import idaes.logger as idaeslog

from idaes.core.util.config import is_physical_parameter_block, DefaultBool

# Additional import for the unit operation
from pyomo.environ import Var, asin, cos
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from pyomo.network import Port
import idaes.core.util.scaling as iscale
from pyomo.network import Arc

from idaes.models_extra.power_generation.unit_models.helm.phase_separator import (
    HelmPhaseSeparator,
)
from idaes.models_extra.power_generation.unit_models.helm.mixer import HelmMixer
from idaes.core.util.constants import Constants as const
from idaes.core.solvers import get_solver

__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"
__version__ = "2.0.0"


@declare_process_block_class("Drum")
class DrumData(UnitModelBlockData):
    """
    Boiler Drum Unit Operation Class
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=DefaultBool,
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.componentPhase,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of material balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.enthalpyTotal,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat transfer term construction flag",
            doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
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

    def build(self):
        """
        Begin building model (pre-DAE transformation).
        Args:
            None
        Returns:
            None
        """

        # Call UnitModel.build to setup dynamics
        super().build()

        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_geometry()

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer,
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
        )

        self.flash = HelmPhaseSeparator(
            dynamic=False, property_package=self.config.property_package
        )

        self.mixer = HelmMixer(
            dynamic=False,
            property_package=self.config.property_package,
            inlet_list=["FeedWater", "SaturatedWater"],
        )

        # Inlet Ports
        # FeedWater to Drum (from Pipe or Economizer)
        self.feedwater_inlet = Port(extends=self.mixer.FeedWater)
        # Sat water from water wall
        self.water_steam_inlet = Port(extends=self.flash.inlet)

        # Exit Ports
        # Liquid to Downcomer
        # self.liquid_outlet = Port(extends=self.mixer.outlet)
        self.add_outlet_port("liquid_outlet", self.control_volume)
        # Steam to superheaters
        self.steam_outlet = Port(extends=self.flash.vap_outlet)

        # constraint to make pressures of two inlets of drum mixer the same
        @self.Constraint(self.flowsheet().time, doc="Mixter pressure identical")
        def mixer_pressure_eqn(b, t):
            return (
                b.mixer.SaturatedWater.pressure[t] * 1e-6
                == b.mixer.FeedWater.pressure[t] * 1e-6
            )

        self.stream_flash_out = Arc(
            source=self.flash.liq_outlet, destination=self.mixer.SaturatedWater
        )

        # Pyomo arc connect flash liq_outlet with mixer SaturatedWater inlet
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

        # connect internal units (Mixer to Water Tank Model)
        # Mixer Outlet (mixed_state) to unit control volume.properties_in
        @self.Constraint(self.flowsheet().time)
        def connection_material_balance(b, t):
            return (
                1e-4 * b.mixer.mixed_state[t].flow_mol
                == b.control_volume.properties_in[t].flow_mol * 1e-4
            )

        @self.Constraint(self.flowsheet().time)
        def connection_enthalpy_balance(b, t):
            return (
                b.mixer.mixed_state[t].enth_mol * 1e-4
                == b.control_volume.properties_in[t].enth_mol * 1e-4
            )

        @self.Constraint(self.flowsheet().time)
        def connection_pressure_balance(b, t):
            return (
                b.mixer.mixed_state[t].pressure * 1e-6
                == b.control_volume.properties_in[t].pressure * 1e-6
            )

        # Add object references
        self.volume = pyo.Reference(self.control_volume.volume)

        # Set references to balance terms at unit level
        if (
            self.config.has_heat_transfer is True
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            self.heat_duty = pyo.Reference(self.control_volume.heat)

        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != "none"
        ):
            self.deltaP = pyo.Reference(self.control_volume.deltaP)

        # Set Unit Geometry and Holdup Volume
        self._set_geometry()

        # Construct performance equations
        self._make_performance()

    def _set_geometry(self):
        """
        Define the geometry of the unit as necessary
        """
        units_meta = self.config.property_package.get_metadata()

        # Inside diameter of drum
        self.drum_diameter = Var(
            initialize=1.0,
            doc="Inside diameter of drum",
            units=units_meta.get_derived_units("length"),
        )
        # Length of drum
        self.drum_length = Var(
            initialize=10,
            doc="Horizontal length of drum",
            units=units_meta.get_derived_units("length"),
        )
        # Number of downcomers connected at the bottom of drum,
        # used to calculate contrac
        self.number_downcomers = Var(
            initialize=4, doc="Number of downcomers connected to drum"
        )
        # Inside diameter of downcomer
        self.downcomer_diameter = Var(
            initialize=0.6,
            doc="Inside diameter of downcomer",
            units=units_meta.get_derived_units("length"),
        )

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.
        """
        units_meta = self.config.property_package.get_metadata()

        # Add performance variables
        self.drum_level = Var(
            self.flowsheet().time,
            within=pyo.PositiveReals,
            initialize=1.0,
            doc="Water level from the bottom of the drum",
            units=units_meta.get_derived_units("length"),
        )

        # Velocity of fluid inside downcomer pipe
        self.downcomer_velocity = Var(
            self.flowsheet().time,
            initialize=10.0,
            doc="Liquid water velocity at the top of downcomer",
            units=units_meta.get_derived_units("velocity"),
        )

        # Pressure change due to contraction
        self.deltaP_contraction = Var(
            self.flowsheet().time,
            initialize=-1.0,
            doc="Pressure change due to contraction",
            units=units_meta.get_derived_units("pressure"),
        )

        # Pressure change due to gravity
        self.deltaP_gravity = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Pressure change due to gravity",
            units=units_meta.get_derived_units("pressure"),
        )

        # Radius expression
        @self.Expression(doc="Radius of drum")
        def drum_radius(b):
            return 0.5 * b.drum_diameter

        # Expression for the angle from the drum center
        # to the circumference point at water level
        @self.Expression(self.flowsheet().time, doc="angle of water level")
        def alpha_drum(b, t):
            return (
                asin((b.drum_level[t] - b.drum_radius) / b.drum_radius) / pyo.units.rad
            )

        # Constraint for volume liquid in drum
        @self.Constraint(self.flowsheet().time, doc="volume of liquid in drum")
        def volume_eqn(b, t):
            return (
                b.volume[t]
                == (
                    (b.alpha_drum[t] + 0.5 * const.pi) * b.drum_radius**2
                    + b.drum_radius
                    * cos(b.alpha_drum[t])
                    * (b.drum_level[t] - b.drum_radius)
                )
                * b.drum_length
            )

        # Equation for velocity at the entrance of downcomer
        @self.Constraint(self.flowsheet().time, doc="Velocity at entrance of downcomer")
        def velocity_eqn(b, t):
            return (
                b.downcomer_velocity[t]
                * 0.25
                * const.pi
                * b.downcomer_diameter**2
                * b.number_downcomers
                == b.control_volume.properties_out[t].flow_vol
            )

        # Pressure change equation for contraction
        # (acceleration pressure change)
        @self.Constraint(
            self.flowsheet().time, doc="pressure change due to contraction"
        )
        def pressure_change_contraction_eqn(b, t):
            return (
                1e-3 * b.deltaP_contraction[t]
                == -1e-3
                * 0.75
                * b.control_volume.properties_out[t].dens_mass_phase["Liq"]
                * b.downcomer_velocity[t] ** 2
            )

        # Pressure change equation for gravity, density*gravity*height
        g_units = units_meta.get_derived_units("acceleration")

        @self.Constraint(self.flowsheet().time, doc="pressure change due to gravity")
        def pressure_change_gravity_eqn(b, t):
            return (
                1e-3 * b.deltaP_gravity[t]
                == 1e-3
                * b.control_volume.properties_out[t].dens_mass_phase["Liq"]
                * pyo.units.convert(const.acceleration_gravity, to_units=g_units)
                * b.drum_level[t]
            )

        # Total pressure change equation
        @self.Constraint(self.flowsheet().time, doc="pressure drop")
        def pressure_change_total_eqn(b, t):
            return 1e-3 * b.deltaP[t] == 1e-3 * (
                b.deltaP_contraction[t] + b.deltaP_gravity[t]
            )

    def set_initial_condition(self):
        if self.config.dynamic is True:
            self.control_volume.material_accumulation[:, :, :].value = 0
            self.control_volume.energy_accumulation[:, :].value = 0
            self.control_volume.material_accumulation[0, :, :].fix(0)
            self.control_volume.energy_accumulation[0, :].fix(0)

    def initialize_build(
        blk,
        state_args_feedwater=None,
        state_args_water_steam=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Drum initialization routine.

        Keyword Arguments:
            state_args_feedwater : a dict of arguments to be passed to the
                           property package(s) for the control_volume of the
                           model to provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            state_args_steam : a dict of arguments to be passed to the property
                           package(s) for the control_volume of the model to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = return solver state for each step in subroutines
                     * 3 = include solver output infomation (tee=True)

            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        init_log.info_low("Starting initialization...")
        # fix FeedWater Inlet
        flags_fw = fix_state_vars(blk.mixer.FeedWater_state, state_args_feedwater)

        # expecting 2 DOF due to pressure driven constraint
        if degrees_of_freedom(blk) != 2:
            raise Exception(degrees_of_freedom(blk))

        blk.flash.initialize(
            state_args_water_steam=state_args_water_steam,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )
        init_log.info("Initialization Step 1 Complete.")

        blk.mixer.SaturatedWater.flow_mol[:].fix(blk.flash.liq_outlet.flow_mol[0].value)
        blk.mixer.SaturatedWater.pressure[:].fix(blk.flash.liq_outlet.pressure[0].value)
        blk.mixer.SaturatedWater.enth_mol[:].fix(blk.flash.liq_outlet.enth_mol[0].value)
        blk.mixer.initialize(outlvl=outlvl, optarg=optarg, solver=solver)
        init_log.info("Initialization Step 2 Complete.")

        blk.control_volume.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, hold_state=False
        )
        init_log.info("Initialization Step 3 Complete.")

        # fix flash Inlet
        flags_steam = fix_state_vars(blk.flash.mixed_state, state_args_water_steam)
        # unfix inlets (connected with arc)
        blk.mixer.SaturatedWater.flow_mol[:].unfix()
        blk.mixer.SaturatedWater.enth_mol[:].unfix()
        blk.mixer.SaturatedWater.pressure[:].unfix()
        blk.mixer.FeedWater.pressure[0].unfix()

        # solve model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))
        revert_state_vars(blk.mixer.FeedWater_state, flags_fw)
        revert_state_vars(blk.flash.mixed_state, flags_steam)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        for v in self.deltaP_gravity.values():
            if iscale.get_scaling_factor(v, warning=True) is None:
                iscale.set_scaling_factor(v, 1e-3)

        for v in self.deltaP_contraction.values():
            if iscale.get_scaling_factor(v, warning=True) is None:
                iscale.set_scaling_factor(v, 1e-3)

        for t, c in self.pressure_change_contraction_eqn.items():
            sf = iscale.get_scaling_factor(
                self.deltaP_contraction[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for t, c in self.pressure_change_gravity_eqn.items():
            sf = iscale.get_scaling_factor(
                self.deltaP_gravity[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for t, c in self.pressure_change_total_eqn.items():
            sf = iscale.get_scaling_factor(self.deltaP[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf, overwrite=False)
