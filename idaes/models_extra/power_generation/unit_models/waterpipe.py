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
Pipe model for water or steam, Liq or Vap phase must be provided (mixed phase
not supported).

main equations:

* Heat is given (zero if adiabatic)
* Calculate pressure change due to friction and gravity
* Assume enthalpy_in = enthalpy_out + heat

Created: April 2019 by Jinliang Ma
"""
# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

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

from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.constants import Constants as const
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog


# Additional import for the unit operation
from pyomo.environ import value, Var, Reference

__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"
__version__ = "2.0.0"


@declare_process_block_class("WaterPipe")
class WaterPipeData(UnitModelBlockData):
    """
    Water or steam pipe Unit Class
    """

    CONFIG = UnitModelBlockData.CONFIG()
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
            default=True,
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
    CONFIG.declare(
        "contraction_expansion_at_end",
        ConfigValue(
            default="None",
            domain=In(["None", "contraction", "expansion"]),
            description="Any contraction or expansion at the end",
            doc="Define if pressure drop due to contraction"
            " or expansion needs to be considered",
        ),
    )
    CONFIG.declare(
        "water_phase",
        ConfigValue(
            default="Liq",
            domain=In(["Liq", "Vap"]),
            description="Water phase",
            doc="""Define water phase for property calls,
mixed phase not supported""",
        ),
    )

    def build(self):
        """
        Begin building model
        """
        # Call UnitModel.build to setup dynamics
        super(WaterPipeData, self).build()

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

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Add object references
        self.volume = Reference(self.control_volume.volume)

        # Set references to balance terms at unit level
        if (
            self.config.has_heat_transfer is True
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            self.heat_duty = Reference(self.control_volume.heat)

        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != "none"
        ):
            self.deltaP = Reference(self.control_volume.deltaP)

        # Set Unit Geometry and Volume
        self._set_geometry()

        # Construct performance equations
        self._make_performance()

    def _set_geometry(self):
        """
        Define the geometry of the unit as necessary
        """
        # Number of pipe
        self.number_of_pipes = Var(initialize=4, doc="Number of pipes")
        # Length of pipe
        self.length = Var(initialize=10.0, doc="Total length of the straight pipe")
        # Elevation change of pipe, elevation of outlet - elevation of inlet
        self.elevation_change = Var(
            initialize=10.0, doc="Change of elevation from inlet to outlet"
        )
        # Inside diameter of pipe
        self.diameter = Var(initialize=0.6, doc="Inside diameter of pipe")

        if not self.config.contraction_expansion_at_end == "None":
            self.area_ratio = Var(
                initialize=1, doc="Cross section area ratio of exit end to pipe"
            )

        # Volume constraint
        @self.Constraint(self.flowsheet().time, doc="Total volume of all pipes")
        def volume_eqn(b, t):
            return (
                b.volume[t]
                == 0.25 * const.pi * b.diameter**2 * b.length * b.number_of_pipes
            )

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.
        """
        phase = self.config.water_phase

        # Add performance variables
        # Velocity of fluid inside pipe
        self.velocity = Var(
            self.flowsheet().time, initialize=10.0, doc="Fluid velocity inside pipe"
        )

        # Reynolds number
        self.N_Re = Var(
            self.flowsheet().time, initialize=10000.0, doc="Reynolds number"
        )

        # Darcy friction factor
        self.friction_factor_darcy = Var(
            self.flowsheet().time, initialize=0.005, doc="Darcy friction factor"
        )

        # Correction factor for pressure drop due to friction
        self.fcorrection_dp = Var(
            initialize=1.0, doc="Correction factor for pressure drop"
        )

        # Pressure change due to friction
        self.deltaP_friction = Var(
            self.flowsheet().time,
            initialize=-1.0,
            doc="Pressure change due to friction",
        )

        # Pressure change due to gravity
        self.deltaP_gravity = Var(
            self.flowsheet().time, initialize=0.0, doc="Pressure change due to gravity"
        )

        # Pressure change due to area change at end
        self.deltaP_area_change = Var(
            self.flowsheet().time,
            initialize=0.0,
            doc="Pressure change due to area change",
        )

        # Equation for calculating velocity
        @self.Constraint(self.flowsheet().time, doc="Velocity of fluid inside pipe")
        def velocity_eqn(b, t):
            return (
                b.velocity[t] * 0.25 * const.pi * b.diameter**2 * b.number_of_pipes
                == b.control_volume.properties_in[t].flow_vol
            )

        # Equation for calculating Reynolds number
        @self.Constraint(self.flowsheet().time, doc="Reynolds number")
        def Reynolds_number_eqn(b, t):
            return (
                b.N_Re[t] * b.control_volume.properties_in[t].visc_d_phase[phase]
                == b.diameter
                * b.velocity[t]
                * b.control_volume.properties_in[t].dens_mass_phase[phase]
            )

        # Friction factor expression depending on laminar or turbulent flow
        @self.Constraint(
            self.flowsheet().time,
            doc="Darcy friction factor as" " a function of Reynolds number",
        )
        def friction_factor_darcy_eqn(b, t):
            return (
                b.friction_factor_darcy[t] * b.N_Re[t] ** (0.25)
                == 0.3164 * b.fcorrection_dp
            )

        # Pressure change equation for friction
        @self.Constraint(self.flowsheet().time, doc="Pressure change due to friction")
        def pressure_change_friction_eqn(b, t):
            return (
                b.deltaP_friction[t] * b.diameter
                == -0.5
                * b.control_volume.properties_in[t].dens_mass_phase[phase]
                * b.velocity[t] ** 2
                * b.friction_factor_darcy[t]
                * b.length
            )

        # Pressure change equation for gravity
        @self.Constraint(self.flowsheet().time, doc="Pressure change due to gravity")
        def pressure_change_gravity_eqn(b, t):
            return (
                b.deltaP_gravity[t]
                == -const.acceleration_gravity
                * b.control_volume.properties_in[t].dens_mass_phase[phase]
                * b.elevation_change
            )

        # Pressure change equation for contraction or expansion
        @self.Constraint(self.flowsheet().time, doc="Pressure change due to gravity")
        def pressure_change_area_change_eqn(b, t):
            if self.config.contraction_expansion_at_end == "contraction":
                return (
                    b.deltaP_area_change[t]
                    == -(0.1602 * b.area_ratio**2 - 0.646 * b.area_ratio + 1.4858)
                    * 0.5
                    * b.control_volume.properties_out[t].dens_mass_phase[phase]
                    * (b.velocity[t] / b.area_ratio) ** 2
                    + 0.5
                    * b.control_volume.properties_out[t].dens_mass_phase[phase]
                    * b.velocity[t] ** 2
                )
            elif self.config.contraction_expansion_at_end == "expansion":
                return (
                    b.deltaP_area_change[t]
                    == b.control_volume.properties_out[t].dens_mass_phase[phase]
                    * b.velocity[t] ** 2
                    * (b.area_ratio - 1)
                    / b.area_ratio**2
                )
            else:
                return b.deltaP_area_change[t] == 0

        # Total pressure change equation
        @self.Constraint(self.flowsheet().time, doc="Pressure drop")
        def pressure_change_total_eqn(b, t):
            return b.deltaP[t] == (
                b.deltaP_friction[t] + b.deltaP_gravity[t] + b.deltaP_area_change[t]
            )

    def set_initial_condition(self):
        if self.config.dynamic is True:
            self.control_volume.material_accumulation[:, :, :].value = 0
            self.control_volume.energy_accumulation[:, :].value = 0
            self.control_volume.material_accumulation[0, :, :].fix(0)
            self.control_volume.energy_accumulation[0, :].fix(0)

    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        WaterPipe initialization routine.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) for the control_volume of the model to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            outlvl : sets output level of initialisation routine
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

        flags = blk.control_volume.initialize(
            outlvl=outlvl + 1, optarg=optarg, solver=solver, state_args=state_args
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # Fix outlet enthalpy and pressure
        for t in blk.flowsheet().time:
            blk.control_volume.properties_out[t].pressure.fix(
                value(blk.control_volume.properties_in[t].pressure)
            )
        blk.pressure_change_total_eqn.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Unfix outlet enthalpy and pressure
        for t in blk.flowsheet().time:
            blk.control_volume.properties_out[t].pressure.unfix()
        blk.pressure_change_total_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        blk.control_volume.release_state(flags, outlvl)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        pass
