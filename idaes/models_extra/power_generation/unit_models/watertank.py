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
Watertank model
The water tank has only one inlet and one outlet

main assumptions:

1.- Heat loss is a variable given by the user (zero heat loss can be
specified if adiabatic)
2.- Calculate pressure change due to gravity based on water level
and contraction to downcomer
3.- Water level is either fixed for steady-state model or calculated for
dynamic model
4.- Assume enthalpy_in = enthalpy_out + heat loss

5.- Subcooled water from economizer and saturated water from waterwall
are mixed before entering the tank

Created: November 04 2020
"""
# Import Pyomo libraries
from pyomo.environ import value, Var, Reference, acos
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
from idaes.core.solvers import get_solver
from idaes.core.util.constants import Constants as const

import idaes.logger as idaeslog

__author__ = "Boiler Subsystem Team (J. Ma, D. Caballero, M. Zamarripa)"
__version__ = "2.0.0"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("WaterTank")
class WaterTankData(UnitModelBlockData):
    """
    Water Tank Unit Operation Class
    """

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "tank_type",
        ConfigValue(
            default="simple_tank",
            domain=In(
                [
                    "simple_tank",
                    "rectangular_tank",
                    "vertical_cylindrical_tank",
                    "horizontal_cylindrical_tank",
                ]
            ),
            description="Flag indicating the tank type",
            doc="""Flag indicating the type of tank to be modeled, and
then calculate the volume of the filled level consequently,
**default** - simple_tank.
**Valid values:** {
**simple_tank** - use a general tank and provide the area,
**rectangular_tank** - use a rectangular tank and provide the width and length,
**vertical_cylindrical_tank** - use a vertical cylindrical tank
and provide the diameter,
**horizontal_cylindrical_tank** - use a horizontal cylindrical tank and
provide the length and diameter.}""",
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
            balance_type=self.config.material_balance_type
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer,
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
        )

        # Add Inlet and Outlet Ports
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

        # Set Unit Geometry and Holdup Volume
        self._set_geometry()

        # Construct performance equations
        self._make_performance()

    def _set_geometry(self):
        """
        Define the geometry of the unit as necessary
        """
        if self.config.tank_type == "simple_tank":
            # Declare a variable for cross sectional area
            self.tank_cross_sect_area = Var(
                initialize=1.0, doc="Cross-sectional" " area of the tank"
            )

        elif self.config.tank_type == "rectangular_tank":
            # Declare variables for width and length
            self.tank_width = Var(initialize=1.0, doc="Width of the tank")
            self.tank_length = Var(initialize=1.0, doc="Length of the tank")

        elif (
            self.config.tank_type == "horizontal_cylindrical_tank"
            or "vertical_cylindrical_tank"
        ):
            # Declare a variable for diameter of the tank
            self.tank_diameter = Var(initialize=0.5, doc="Inside diameter of the tank")
            if self.config.tank_type == "horizontal_cylindrical_tank":
                # Declare a variable for length of the tank
                self.tank_length = Var(initialize=1, doc="Length of the tank")

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model
        """

        # Add performance variables
        self.tank_level = Var(
            self.flowsheet().time, initialize=1.0, doc="Water level from in the tank"
        )

        # Auxiliar expressions for volume
        # Rectangular tank
        if self.config.tank_type == "rectangular_tank":
            # Calculation of cross-sectional area of the rectangle
            @self.Expression(doc="Cross-sectional area of the tank")
            def tank_cross_sect_area(b):
                return b.tank_width * b.tank_length

        # Vertical cylindrical tank
        elif self.config.tank_type == "vertical_cylindrical_tank":

            @self.Expression(doc="Radius of the tank")
            def tank_radius(b):
                return b.tank_diameter / 2

            # Calculation of cross-sectional area of the vertical cylinder
            @self.Expression(doc="Cross-sectional area of the tank")
            def tank_cross_sect_area(b):
                return const.pi * b.tank_radius**2

        # Horizontal cylindrical tank
        elif self.config.tank_type == "horizontal_cylindrical_tank":
            # Calculation of area covered by the liquid level
            # at one end of the tank
            @self.Expression(doc="Radius of the tank")
            def tank_radius(b):
                return b.tank_diameter / 2

            # Angle of the circular sector used to calculate the area

            @self.Expression(
                self.flowsheet().time,
                doc="Angle of the circular" " sector of liquid level",
            )
            def alpha_tank(b, t):
                return 2 * acos((b.tank_radius - b.tank_level[t]) / b.tank_radius)

            @self.Expression(
                self.flowsheet().time,
                doc="Area covered by the liquid level" " at one end of the tank",
            )
            def tank_area(b, t):
                return (
                    0.5 * b.alpha_tank[t] * b.tank_radius**2
                    - (b.tank_radius - b.tank_level[t])
                    * (2 * b.tank_radius * b.tank_level[t] - b.tank_level[t] ** 2)
                    ** 0.5
                )

        # Constraint for volume of the liquid in tank
        @self.Constraint(self.flowsheet().time, doc="volume of liquid in the tank")
        def volume_eqn(b, t):
            if self.config.tank_type == "horizontal_cylindrical_tank":
                return b.volume[t] == b.tank_length * b.tank_area[t]
            else:
                return b.volume[t] == b.tank_level[t] * b.tank_cross_sect_area

        # Pressure change equation due gravity
        @self.Constraint(self.flowsheet().time, doc="pressure drop")
        def pressure_change_eqn(b, t):
            return (
                b.deltaP[t]
                == b.control_volume.properties_in[t].dens_mass_phase["Liq"]
                * const.acceleration_gravity
                * b.tank_level[t]
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
        Water tank initialization routine.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
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

        flags = blk.control_volume.initialize(
            state_args=state_args, outlvl=outlvl, optarg=optarg, solver=solver
        )
        init_log.info_high("Initialization Step 1 Complete.")

        # Fix outlet pressure
        for t in blk.flowsheet().time:
            blk.control_volume.properties_out[t].pressure.fix(
                value(blk.control_volume.properties_in[t].pressure)
            )
        blk.pressure_change_eqn.deactivate()

        # solve model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Unfix outlet pressure
        for t in blk.flowsheet().time:
            blk.control_volume.properties_out[t].pressure.unfix()
        blk.pressure_change_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        blk.control_volume.release_state(flags, outlvl)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        pass
