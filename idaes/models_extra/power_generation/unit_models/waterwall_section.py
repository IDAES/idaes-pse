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
Waterwall section model

main equations:

* Heat is given by fire-side boiler model
* Calculate pressure change due to friction and gravity
* Calculate slag layer wall temperature
* Consider a layer of metal and a layer of slag


"""
# Import Pyomo libraries
from pyomo.environ import (
    value,
    Var,
    Param,
    asin,
    cos,
    sqrt,
    log10,
    PositiveReals,
    Reference,
    units as pyunits,
)
from pyomo.dae import DerivativeVar
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.core.expr.current import Expr_if

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

from idaes.core.util.config import is_physical_parameter_block, DefaultBool
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog


# Additional import for the unit operation
from idaes.core.util.constants import Constants as const


__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"
__version__ = "2.0.0"


@declare_process_block_class("WaterwallSection")
class WaterwallSectionData(UnitModelBlockData):
    """
    WaterwallSection Unit Class
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
    CONFIG.declare(
        "rigorous_boiling",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat of reaction term construction flag",
            doc="""Indicates whether terms for heat of reaction terms should be
constructed,
**default** - False.
**Valid values:** {
**True** - include heat of reaction terms,
**False** - exclude heat of reaction terms.}""",
        ),
    )

    def build(self):
        """
        Begin building model (pre-DAE transformation)


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
        # This model requires the IAPWS95 property package with the mixed phase
        # option, therefore, phase equilibrium calculations are handled by
        # the property package.
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

        # Set Unit Geometry and Holdup Volume
        self._set_geometry()

        # Construct performance equations
        self._make_performance()

    def _set_geometry(self):
        """
        Define the geometry of the unit as necessary, and link to holdup volume
        """
        units_meta = self.config.property_package.get_metadata()
        # Total projected wall area of waterwall section from fire-side model
        self.projected_area = Var(
            initialize=10,
            doc="Total projected wall area of waterwall section",
            units=units_meta.get_derived_units("area"),
        )
        # Number of waterwall tubes
        self.number_tubes = Var(initialize=4, doc="Number of waterwall tubes")
        # Height of waterwall section, given by boiler model
        self.height = Var(
            initialize=5.0,
            doc="Height of waterwall section",
            units=units_meta.get_derived_units("length"),
        )
        # Length of waterwall tubes calculated based on given total area
        # and perimeter of waterwall
        self.tube_length = Var(
            initialize=5.0,
            doc="length waterwall tube",
            units=units_meta.get_derived_units("length"),
        )
        # Inner diameter of waterwall tubes
        self.tube_diameter = Var(
            initialize=0.05,
            doc="Inside diameter of waterwall tube",
            units=units_meta.get_derived_units("length"),
        )
        # Inside radius of waterwall tube
        @self.Expression(doc="Inside radius of waterwall tube")
        def radius_inner(b):
            return 0.5 * b.tube_diameter

        # Total cross section area of fluid flow
        @self.Expression(doc="Cross section area of fluid")
        def area_cross_fluid_total(b):
            return 0.25 * const.pi * b.tube_diameter**2 * b.number_tubes

        # Tube thickness
        self.tube_thickness = Var(
            initialize=0.005,
            doc="Thickness of waterwall tube",
            units=units_meta.get_derived_units("length"),
        )
        # Outside radius of waterwall tube
        @self.Expression(doc="Outside radius of waterwall tube")
        def radius_outer(b):
            return b.radius_inner + b.tube_thickness

        # Thickness of waterwall fin
        self.fin_thickness = Var(
            initialize=0.004,
            doc="Thickness of waterwall Fin",
            units=units_meta.get_derived_units("length"),
        )
        # Half of the waterwall fin thickness
        @self.Expression(doc="Half of the waterwall fin thickness")
        def fin_thickness_half(b):
            return 0.5 * b.fin_thickness

        # Length of waterwall fin
        self.fin_length = Var(
            initialize=0.005,
            doc="Length of waterwall fin",
            units=units_meta.get_derived_units("length"),
        )
        # Thickness of slag layer
        self.slag_thickness = Var(
            self.flowsheet().time,
            bounds=(0.0001, 0.009),
            initialize=0.001,
            doc="Thickness of slag layer",
            units=units_meta.get_derived_units("length"),
        )

        @self.Expression(doc="Pitch of two neighboring tubes")
        def pitch(b):
            return b.fin_length + b.radius_outer * 2.0

        # Equivalent tube length (not neccesarily equal to height)
        @self.Constraint(doc="Equivalent length of tube")
        def tube_length_eqn(b):
            return b.tube_length * b.pitch * b.number_tubes == b.projected_area

        @self.Expression(doc="Angle at joint of tube and fin")
        def alpha_tube(b):
            return asin(b.fin_thickness_half / b.radius_outer)

        @self.Expression(
            self.flowsheet().time,
            doc="Angle at joint of tube and " "fin at outside slag layer",
        )
        def alpha_slag(b, t):
            return asin(
                (b.fin_thickness_half + b.slag_thickness[t])
                / (b.radius_outer + b.slag_thickness[t])
            )

        @self.Expression(doc="Perimeter of interface between slag and tube")
        def perimeter_interface(b):
            return (
                (const.pi - 2 * b.alpha_tube) * b.radius_outer
                + b.pitch
                - 2 * b.radius_outer * cos(b.alpha_tube)
            )

        @self.Expression(doc="Perimeter on the inner tube side")
        def perimeter_ts(b):
            return const.pi * b.tube_diameter

        @self.Expression(self.flowsheet().time, doc="Perimeter on the outer slag side")
        def perimeter_ss(b, t):
            return (
                (const.pi - 2 * b.alpha_slag[t])
                * (b.radius_outer + b.slag_thickness[t])
                + b.pitch
                - 2 * (b.radius_outer + b.slag_thickness[t]) * cos(b.alpha_slag[t])
            )

        # Cross section area of tube and fin metal
        @self.Expression(doc="Cross section area of tube and fin metal")
        def area_cross_metal(b):
            return (
                const.pi * (b.radius_outer**2 - b.radius_inner**2)
                + b.fin_thickness * b.fin_length
            )

        # Cross section area of slag layer
        @self.Expression(
            self.flowsheet().time, doc="Cross section area of slag layer per tube"
        )
        def area_cross_slag(b, t):
            return b.perimeter_interface * b.slag_thickness[t]

        # Volume constraint
        @self.Constraint(
            self.flowsheet().time, doc="waterwall fluid volume of all tubes"
        )
        def volume_eqn(b, t):
            return (
                b.volume[t]
                == 0.25
                * const.pi
                * b.tube_diameter**2
                * b.tube_length
                * b.number_tubes
            )

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.
        """
        units_meta = self.config.property_package.get_metadata()
        self.fcorrection_dp = Var(
            initialize=1.2,
            within=PositiveReals,
            doc="correction factor for pressure drop due to acceleration"
            "and unsmooth tube applied to friction term",
        )
        # Thermal conductivity of metal
        self.therm_cond_metal = Param(
            initialize=43.0,
            mutable=True,
            doc="Thermal conductivity of tube metal",
            units=units_meta.get_derived_units("thermal_conductivity"),
        )
        # Thermal conductivity of slag
        self.therm_cond_slag = Param(
            initialize=1.3,
            mutable=True,
            doc="Thermal conductivity of slag",
            units=units_meta.get_derived_units("thermal_conductivity"),
        )
        # Heat capacity of metal
        self.cp_metal = Param(
            initialize=500.0,
            mutable=True,
            doc="Heat capacity of tube metal",
            units=units_meta.get_derived_units("heat_capacity_mass"),
        )
        # Heat Capacity of slag
        self.cp_slag = Param(
            initialize=250,
            mutable=True,
            doc="Heat capacity of slag",
            units=units_meta.get_derived_units("heat_capacity_mass"),
        )
        # Density of metal
        self.dens_metal = Param(
            initialize=7800.0,
            mutable=True,
            doc="Density of tube metal",
            units=units_meta.get_derived_units("density_mass"),
        )
        # Density of slag
        self.dens_slag = Param(
            initialize=2550,
            mutable=True,
            doc="Density of slag",
            units=units_meta.get_derived_units("density_mass"),
        )
        # Shape factor of tube metal conduction resistance
        # based on projected area
        self.fshape_metal = Param(
            initialize=0.7718, mutable=True, doc="Shape factor of tube metal conduction"
        )
        # Shape factor of slag conduction resistance
        # based on projected area
        self.fshape_slag = Param(
            initialize=0.6858, mutable=True, doc="Shape factor of slag conduction"
        )
        # Shape factor of convection on tube side resistance
        # based on projected area
        self.fshape_conv = Param(
            initialize=0.8496, mutable=True, doc="Shape factor of convection"
        )

        # Heat conduction resistance of half of metal wall thickness
        # based on interface perimeter
        @self.Expression(doc="half metal layer conduction resistance")
        def half_resistance_metal(b):
            return (
                b.tube_thickness
                / 2
                / b.therm_cond_metal
                * b.fshape_metal
                * b.perimeter_interface
                / b.pitch
            )

        # Heat conduction resistance of half of slag thickness
        # based on mid slag layer perimeter
        @self.Expression(
            self.flowsheet().time, doc="half slag layer conduction resistance"
        )
        def half_resistance_slag(b, t):
            return (
                b.slag_thickness[t]
                / 2
                / b.therm_cond_slag
                * b.fshape_slag
                * (b.perimeter_ss[t] + b.perimeter_interface)
                / 2
                / b.pitch
            )

        # Add performance variables
        # Heat from fire side boiler model
        self.heat_fireside = Var(
            self.flowsheet().time,
            initialize=1e7,
            doc="total heat from fire side model for the section",
            units=units_meta.get_derived_units("power"),
        )
        # Tube boundary wall temperature
        self.temp_tube_boundary = Var(
            self.flowsheet().time,
            initialize=400.0,
            doc="Temperature of tube boundary wall",
            units=units_meta.get_derived_units("temperature"),
        )
        # Tube center point wall temperature
        self.temp_tube_center = Var(
            self.flowsheet().time,
            initialize=450.0,
            doc="Temperature of tube center wall",
            units=units_meta.get_derived_units("temperature"),
        )
        # Slag boundary wall temperature
        self.temp_slag_boundary = Var(
            self.flowsheet().time,
            initialize=600.0,
            doc="Temperature of slag boundary wall",
            units=units_meta.get_derived_units("temperature"),
        )
        # Slag center point wall temperature
        self.temp_slag_center = Var(
            self.flowsheet().time,
            initialize=500.0,
            doc="Temperature of slag layer center point",
            units=units_meta.get_derived_units("temperature"),
        )

        # Energy holdup for slag layer
        self.energy_holdup_slag = Var(
            self.flowsheet().time,
            initialize=1e4,
            doc="Energy holdup of slag layer",
            units=(
                units_meta.get_derived_units("energy")
                * units_meta.get_derived_units("length") ** -1
            ),
        )

        # Energy holdup for metal (tube + fin)
        self.energy_holdup_metal = Var(
            self.flowsheet().time,
            initialize=1e6,
            doc="Energy holdup of metal",
            units=(
                units_meta.get_derived_units("energy")
                * units_meta.get_derived_units("length") ** -1
            ),
        )

        # Energy accumulation for slag and metal
        if self.config.dynamic is True:
            self.energy_accumulation_slag = DerivativeVar(
                self.energy_holdup_slag,
                wrt=self.flowsheet().time,
                doc="Energy accumulation of slag layer",
            )
            self.energy_accumulation_metal = DerivativeVar(
                self.energy_holdup_metal,
                wrt=self.flowsheet().time,
                doc="Energy accumulation of tube and fin metal",
            )

        def energy_accumulation_term_slag(b, t):
            return b.energy_accumulation_slag[t] if b.config.dynamic else 0

        def energy_accumulation_term_metal(b, t):
            return b.energy_accumulation_metal[t] if b.config.dynamic else 0

        # Velocity of liquid only
        self.velocity_liquid = Var(
            self.flowsheet().time,
            initialize=3.0,
            doc="Velocity of liquid only",
            units=units_meta.get_derived_units("velocity"),
        )

        # Reynolds number based on liquid only flow
        self.N_Re = Var(self.flowsheet().time, initialize=1.0e6, doc="Reynolds number")

        # Prandtl number of liquid phase
        self.N_Pr = Var(self.flowsheet().time, initialize=2.0, doc="Reynolds number")

        # Darcy friction factor
        self.friction_factor_darcy = Var(
            self.flowsheet().time, initialize=0.01, doc="Darcy friction factor"
        )

        # Vapor fraction at inlet
        self.vapor_fraction = Var(
            self.flowsheet().time,
            initialize=0.0,
            doc="Vapor fractoin of vapor-liquid mixture",
        )

        # Liquid fraction at inlet
        self.liquid_fraction = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Liquid fractoin of vapor-liquid mixture",
        )

        # Density ratio of liquid to vapor
        self.ratio_density = Var(
            self.flowsheet().time, initialize=1.0, doc="Liquid to vapor density ratio"
        )

        # Void fraction at inlet
        self.void_fraction = Var(
            self.flowsheet().time, initialize=0.0, doc="void fraction at inlet"
        )

        # Exponent n for gamma at inlet, typical range in (0.75,0.8294)
        self.n_exp = Var(
            self.flowsheet().time, initialize=0.82, doc="exponent for gamma at inlet"
        )

        # Gamma for velocity slip at inlet
        self.gamma = Var(
            self.flowsheet().time,
            initialize=3.0,
            doc="gamma for velocity slip at inlet",
        )

        # Mass flux
        self.mass_flux = Var(
            self.flowsheet().time,
            initialize=1000.0,
            doc="mass flux",
            units=units_meta.get_derived_units("flux_mass"),
        )

        # Reduced pressure
        self.reduced_pressure = Var(
            self.flowsheet().time, initialize=0.85, doc="reduced pressure"
        )

        # Two-phase correction factor
        self.phi_correction = Var(
            self.flowsheet().time,
            initialize=1.01,
            doc="Two-phase flow correction factor",
        )

        # Convective heat transfer coefficient on tube side,
        # typically in range (1000, 5e5)
        self.hconv = Var(
            self.flowsheet().time,
            initialize=30000.0,
            doc="Convective heat transfer coefficient",
            units=units_meta.get_derived_units("heat_transfer_coefficient"),
        )

        # Convective heat transfer coefficient for liquid only,
        # typically in range (1000.0, 1e5)
        self.hconv_liquid = Var(
            self.flowsheet().time,
            initialize=20000.0,
            doc="Convective heat transfer coefficient of liquid only",
            units=units_meta.get_derived_units("heat_transfer_coefficient"),
        )

        # Pool boiling heat transfer coefficient, typically in range (1e4, 5e5)
        self.hpool = Var(
            self.flowsheet().time,
            initialize=1e5,
            doc="Pool boiling heat transfer coefficient",
            units=units_meta.get_derived_units("heat_transfer_coefficient"),
        )

        # Boiling number, typical range in (1e-7, 5e-4) in original formula.
        # we define here as boiling_number_scaled == 1e6*boiling_number
        self.boiling_number_scaled = Var(
            self.flowsheet().time, initialize=1, doc="Scaled boiling number"
        )

        # Enhancement factor, typical range in (1.0, 3.0)
        self.enhancement_factor = Var(
            self.flowsheet().time, initialize=1.3, doc="Enhancement factor"
        )

        # Suppression factor, typical range in (0.005, 1.0)
        self.suppression_factor = Var(
            self.flowsheet().time, initialize=0.03, doc="Suppression factor"
        )

        # Convective heat flux to fluid
        self.heat_flux_conv = Var(
            self.flowsheet().time,
            initialize=7e4,
            doc="Convective heat flux to fluid",
            units=units_meta.get_derived_units("flux_energy"),
        )

        # Slag-tube interface heat flux
        self.heat_flux_interface = Var(
            self.flowsheet().time,
            initialize=100000.0,
            doc="Slag-tube interface heat flux",
            units=units_meta.get_derived_units("flux_energy"),
        )

        # Pressure change due to friction
        self.deltaP_friction = Var(
            self.flowsheet().time,
            initialize=-1000.0,
            doc="Pressure change due to friction",
            units=units_meta.get_derived_units("pressure"),
        )

        # Pressure change due to gravity
        self.deltaP_gravity = Var(
            self.flowsheet().time,
            initialize=-1000.0,
            doc="Pressure change due to gravity",
            units=units_meta.get_derived_units("pressure"),
        )

        # Equation to calculate heat flux to slag boundary
        @self.Expression(self.flowsheet().time, doc="heat flux at slag outer layer")
        def heat_flux_fireside(b, t):
            return b.heat_fireside[t] * b.pitch / (b.projected_area * b.perimeter_ss[t])

        # Equation to calculate slag layer boundary temperature
        @self.Constraint(self.flowsheet().time, doc="slag layer boundary temperature")
        def slag_layer_boundary_temperature_eqn(b, t):
            return b.heat_flux_fireside[t] * b.half_resistance_slag[t] == (
                b.temp_slag_boundary[t] - b.temp_slag_center[t]
            )

        # Equation to calculate heat flux at the slag-metal interface
        @self.Constraint(self.flowsheet().time, doc="heat flux at slag-tube interface")
        def heat_flux_interface_eqn(b, t):
            return b.heat_flux_interface[t] * (
                b.half_resistance_metal + b.half_resistance_slag[t]
            ) == (b.temp_slag_center[t] - b.temp_tube_center[t])

        # Equation to calculate heat flux at tube boundary
        @self.Constraint(
            self.flowsheet().time, doc="convective heat flux at tube boundary"
        )
        def heat_flux_conv_eqn(b, t):
            return b.heat_flux_conv[
                t
            ] * b.fshape_conv * b.perimeter_ts == b.pitch * b.hconv[t] * (
                b.temp_tube_boundary[t] - b.control_volume.properties_in[t].temperature
            )

        # Equation to calculate tube boundary wall temperature
        @self.Constraint(self.flowsheet().time, doc="tube bounary wall temperature")
        def temperature_tube_boundary_eqn(b, t):
            return b.heat_flux_conv[
                t
            ] * b.half_resistance_metal * b.perimeter_ts == b.perimeter_interface * (
                b.temp_tube_center[t] - b.temp_tube_boundary[t]
            )

        # Equation to calculate energy holdup for slag layer per tube length
        @self.Constraint(self.flowsheet().time, doc="energy holdup for slag layer")
        def energy_holdup_slag_eqn(b, t):
            return (
                b.energy_holdup_slag[t]
                == b.temp_slag_center[t]
                * b.cp_slag
                * b.dens_slag
                * b.area_cross_slag[t]
            )

        # Equation to calculate energy holdup for metal
        # (tube + fin) per tube length
        @self.Constraint(self.flowsheet().time, doc="energy holdup for metal")
        def energy_holdup_metal_eqn(b, t):
            return (
                b.energy_holdup_metal[t]
                == b.temp_tube_center[t]
                * b.cp_metal
                * b.dens_metal
                * b.area_cross_metal
            )

        # Energy balance for slag layer
        @self.Constraint(self.flowsheet().time, doc="energy balance for slag layer")
        def energy_balance_slag_eqn(b, t):
            return (
                energy_accumulation_term_slag(b, t)
                == b.heat_flux_fireside[t] * b.perimeter_ss[t]
                - b.heat_flux_interface[t] * b.perimeter_interface
            )

        # Energy balance for metal
        @self.Constraint(self.flowsheet().time, doc="energy balance for metal")
        def energy_balance_metal_eqn(b, t):
            return (
                energy_accumulation_term_metal(b, t)
                == b.heat_flux_interface[t] * b.perimeter_interface
                - b.heat_flux_conv[t] * b.perimeter_ts
            )

        # Expression to calculate slag/tube metal interface wall temperature
        @self.Expression(
            self.flowsheet().time, doc="Slag tube interface wall temperature"
        )
        def temp_interface(b, t):
            return (
                b.temp_tube_center[t]
                + b.heat_flux_interface[t] * b.half_resistance_metal
            )

        # Equations for calculate pressure drop
        # and convective heat transfer coefficient for 2-phase flow
        # Equation to calculate liquid to vapor density ratio
        @self.Constraint(self.flowsheet().time, doc="liquid to vapor density ratio")
        def ratio_density_eqn(b, t):
            return (
                0.01
                * b.ratio_density[t]
                * b.control_volume.properties_in[t].dens_mol_phase["Vap"]
                == 0.01 * b.control_volume.properties_in[t].dens_mol_phase["Liq"]
            )

        # Equation for calculating velocity if the flow is liquid only
        @self.Constraint(self.flowsheet().time, doc="Vecolity of fluid if liquid only")
        def velocity_lo_eqn(b, t):
            return (
                1e-4
                * b.velocity_liquid[t]
                * b.area_cross_fluid_total
                * b.control_volume.properties_in[t].dens_mol_phase["Liq"]
                == 1e-4 * b.control_volume.properties_in[t].flow_mol
            )

        # Equation for calculating Reynolds number if liquid only
        @self.Constraint(self.flowsheet().time, doc="Reynolds number if liquid only")
        def Reynolds_number_eqn(b, t):
            return (
                b.N_Re[t] * b.control_volume.properties_in[t].visc_d_phase["Liq"]
                == b.tube_diameter
                * b.velocity_liquid[t]
                * b.control_volume.properties_in[t].dens_mass_phase["Liq"]
            )

        # Friction factor depending on laminar or turbulent flow,
        # usually always turbulent (>1187.385)
        @self.Constraint(self.flowsheet().time, doc="Darcy friction factor")
        def friction_factor_darcy_eqn(b, t):
            return b.friction_factor_darcy[t] * b.N_Re[t] ** 0.25 / 0.3164 == 1.0

        # Vapor fractoin equation at inlet,
        # add 1e-5 such that vapor fraction is always positive
        @self.Constraint(self.flowsheet().time, doc="Vapor fractoin at inlet")
        def vapor_fraction_eqn(b, t):
            return 100 * b.vapor_fraction[t] == 100 * (
                b.control_volume.properties_in[t].vapor_frac + 1e-5
            )

        # n-exponent equation for inlet
        @self.Constraint(self.flowsheet().time, doc="n-exponent")
        def n_exp_eqn(b, t):
            return 0.001 * (0.8294 - b.n_exp[t]) * b.control_volume.properties_in[
                t
            ].pressure == 8.0478 * units_meta.get_derived_units("pressure")

        # Gamma equation for inlet
        @self.Constraint(self.flowsheet().time, doc="Gamma at inlet")
        def gamma_eqn(b, t):
            return b.gamma[t] == b.ratio_density[t] ** b.n_exp[t]

        # void faction at inlet equation
        @self.Constraint(self.flowsheet().time, doc="Void fractoin at inlet")
        def void_fraction_eqn(b, t):
            return (
                b.void_fraction[t] * (1.0 + b.vapor_fraction[t] * (b.gamma[t] - 1.0))
                == b.vapor_fraction[t] * b.gamma[t]
            )

        # Two-phase flow correction factor equation
        @self.Constraint(self.flowsheet().time, doc="Correction factor")
        def correction_factor_eqn(b, t):
            return (b.phi_correction[t] - 0.027 * b.liquid_fraction[t]) ** 2 == (
                0.973 * b.liquid_fraction[t] + b.vapor_fraction[t] * b.ratio_density[t]
            ) * (0.973 * b.liquid_fraction[t] + b.vapor_fraction[t])

        # Pressure change equation due to friction,
        # -1/2*density*velocity^2*fD/diameter*length*phi^2
        @self.Constraint(self.flowsheet().time, doc="pressure change due to friction")
        def pressure_change_friction_eqn(b, t):
            return (
                0.01 * b.deltaP_friction[t] * b.tube_diameter
                == -0.01
                * b.fcorrection_dp
                * 0.5
                * b.control_volume.properties_in[t].dens_mass_phase["Liq"]
                * b.velocity_liquid[t] ** 2
                * b.friction_factor_darcy[t]
                * b.tube_length
                * b.phi_correction[t] ** 2
            )

        # Pressure change equation due to gravity,
        # density_mixture*gravity*height
        @self.Constraint(self.flowsheet().time, doc="pressure change due to gravity")
        def pressure_change_gravity_eqn(b, t):
            return 1e-3 * b.deltaP_gravity[
                t
            ] == -1e-3 * const.acceleration_gravity * b.height * (
                b.control_volume.properties_in[t].dens_mass_phase["Vap"]
                * b.void_fraction[t]
                + b.control_volume.properties_in[t].dens_mass_phase["Liq"]
                * (1.0 - b.void_fraction[t])
            )

        # Mass flux of vapor-liquid mixture
        # (density*velocity or mass_flow/area)
        @self.Constraint(self.flowsheet().time, doc="mass flux")
        def mass_flux_eqn(b, t):
            return (
                b.mass_flux[t] * b.area_cross_fluid_total
                == b.control_volume.properties_in[t].flow_mol
                * b.control_volume.properties_in[0].mw
            )

        # Liquid fraction at inlet
        @self.Constraint(self.flowsheet().time, doc="liquid fraction")
        def liquid_fraction_eqn(b, t):
            return b.liquid_fraction[t] + b.vapor_fraction[t] == 1.0

        # Total pressure change equation
        @self.Constraint(self.flowsheet().time, doc="pressure drop")
        def pressure_change_total_eqn(b, t):
            return b.deltaP[t] == b.deltaP_friction[t] + b.deltaP_gravity[t]

        # Total heat added to control_volume
        @self.Constraint(
            self.flowsheet().time, doc="total heat added to fluid control_volume"
        )
        def heat_eqn(b, t):
            return (
                b.heat_duty[t]
                == b.number_tubes * b.heat_flux_conv[t] * b.tube_length * b.perimeter_ts
            )

        # Reduced pressure
        @self.Constraint(self.flowsheet().time, doc="reduced pressure")
        def reduced_pressure_eqn(b, t):
            return (
                b.reduced_pressure[t] * self.config.property_package.pressure_crit
                == b.control_volume.properties_in[t].pressure
            )

        # Prandtl number of liquid
        @self.Constraint(self.flowsheet().time, doc="liquid Prandtl number")
        def N_Pr_eqn(b, t):
            return (
                b.N_Pr[t]
                * b.control_volume.properties_in[t].therm_cond_phase["Liq"]
                * b.control_volume.properties_in[0].mw
                == b.control_volume.properties_in[t].cp_mol_phase["Liq"]
                * b.control_volume.properties_in[t].visc_d_phase["Liq"]
            )

        # Forced convection heat transfer coefficient for liquid only
        @self.Constraint(
            self.flowsheet().time,
            doc="forced convection heat transfer " "coefficient for liquid only",
        )
        def hconv_lo_eqn(b, t):
            return (
                b.hconv_liquid[t] * b.tube_diameter
                == 0.023
                * b.N_Re[t] ** 0.8
                * b.N_Pr[t] ** 0.4
                * b.control_volume.properties_in[t].therm_cond_phase["Liq"]
            )

        # Pool boiling heat transfer coefficient
        @self.Constraint(
            self.flowsheet().time, doc="pool boiling heat transfer coefficient"
        )
        def hpool_eqn(b, t):
            return (
                1e-4
                * b.hpool[t]
                * sqrt(
                    pyunits.convert(
                        b.control_volume.properties_in[0].mw,
                        to_units=pyunits.g / pyunits.mol,
                    )
                )
                * (-log10(b.reduced_pressure[t])) ** (0.55)
                == 1e-4
                * 55.0
                * b.reduced_pressure[t] ** 0.12
                * b.heat_flux_conv[t] ** 0.67
            )

        # Boiling number scaled by a factor of 1e6
        @self.Constraint(self.flowsheet().time, doc="boiling number")
        def boiling_number_eqn(b, t):
            return (
                1e-10
                * b.boiling_number_scaled[t]
                * b.control_volume.properties_in[t].dh_vap_mol
                * b.mass_flux[t]
                == b.heat_flux_conv[t] * b.control_volume.properties_in[0].mw * 1e-4
            )

        if self.config.rigorous_boiling is True:
            """
            Due to low contribution to the enhancement factor and highly nonlinear
            constraint, the Martinelli paramter has been removed from the model.
            if required user needs to declare the variable and constraint, and
            update constraint in line 909 to add the Martinelli parameter factor
            """
            # Reciprocal of Martinelli parameter to the power of 0.86,
            # typical range in (1e-3, 1.0)
            self.martinelli_reciprocal_p86 = Var(
                self.flowsheet().time,
                initialize=0.2,
                doc="Reciprocal of Martinelli parameter " "to the power of 0.86",
            )

            # Reciprocal of Martinelli parameter to the power of 0.86
            @self.Constraint(
                self.flowsheet().time,
                doc="Reciprocal of Martineli parameter " "to the power of 0.86",
            )
            def martinelli_reciprocal_p86_eqn(b, t):
                return b.martinelli_reciprocal_p86[t] * b.liquid_fraction[
                    t
                ] ** 0.774 * b.control_volume.properties_in[t].visc_d_phase[
                    "Liq"
                ] ** 0.086 == Expr_if(
                    b.vapor_fraction[t] > 0,
                    (b.vapor_fraction[t]) ** 0.774
                    * b.ratio_density[t] ** 0.43
                    * b.control_volume.properties_in[t].visc_d_phase["Vap"] ** 0.086,
                    0.0,
                )

        # Enhancement factor
        # due to low contribution the reciprocal Martinalli parameter,
        # it can be removed from the enhancement factor equation
        @self.Constraint(
            self.flowsheet().time, doc="Forced convection enhancement factor"
        )
        def enhancement_factor_eqn(b, t):
            if self.config.rigorous_boiling is True:
                return b.enhancement_factor[t] == 1.0 + 24000.0 * (
                    b.boiling_number_scaled[t] / 1e6
                ) ** 1.16 + Expr_if(
                    b.vapor_fraction[t] > 0, 1.37 * b.martinelli_reciprocal_p86[t], 0.0
                )
            else:
                return (
                    b.enhancement_factor[t]
                    == 1.0 + 24000.0 * (b.boiling_number_scaled[t] / 1e6) ** 1.16
                )

        # Suppression factor
        @self.Constraint(self.flowsheet().time, doc="Pool boiler suppression factor")
        def suppression_factor_eqn(b, t):
            return (
                b.suppression_factor[t]
                * (1.0 + 1.15e-6 * b.enhancement_factor[t] ** 2 * b.N_Re[t] ** 1.17)
                == 1.0
            )

        @self.Constraint(
            self.flowsheet().time, doc="convective heat transfer coefficient"
        )
        def hconv_eqn(b, t):
            return (
                1e-3 * b.hconv[t]
                == 1e-3 * b.hconv_liquid[t] * b.enhancement_factor[t]
                + 1e-3 * b.hpool[t] * b.suppression_factor[t]
            )

    def set_initial_condition(self):
        """Initialization of dynamic accumulation terms"""

        if self.config.dynamic is True:
            self.control_volume.material_accumulation[:, :, :].value = 0
            self.control_volume.energy_accumulation[:, :].value = 0
            self.control_volume.material_accumulation[0, :, :].fix(0)
            self.control_volume.energy_accumulation[0, :].fix(0)
            self.energy_accumulation_slag[:].value = 0
            self.energy_accumulation_metal[:].value = 0
            self.energy_accumulation_slag[0].fix(0)
            self.energy_accumulation_metal[0].fix(0)

    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        Waterwall section initialization routine.

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
            blk.control_volume.properties_out[t].enth_mol.fix(
                value(blk.control_volume.properties_in[t].enth_mol)
                + value(blk.heat_fireside[t])
                / value(blk.control_volume.properties_in[t].flow_mol)
            )
            blk.control_volume.properties_out[t].pressure.fix(
                value(blk.control_volume.properties_in[t].pressure) - 1.0
            )
        blk.heat_eqn.deactivate()
        blk.pressure_change_total_eqn.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Unfix outlet enthalpy and pressure
        for t in blk.flowsheet().time:
            blk.control_volume.properties_out[t].enth_mol.unfix()
            blk.control_volume.properties_out[t].pressure.unfix()
        blk.heat_eqn.activate()
        blk.pressure_change_total_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        blk.control_volume.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for t, c in self.energy_holdup_slag_eqn.items():
            s = iscale.get_scaling_factor(
                self.energy_holdup_slag[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.friction_factor_darcy_eqn.items():
            s = iscale.get_scaling_factor(self.N_Re[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.volume_eqn.items():
            s = iscale.get_scaling_factor(self.volume[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.heat_flux_conv_eqn.items():
            s = iscale.get_scaling_factor(
                self.heat_flux_conv[t], default=1, warning=True
            )
            s *= iscale.get_scaling_factor(self.tube_diameter, default=1, warning=True)
            iscale.constraint_scaling_transform(c, s / 10.0, overwrite=False)
        for t, c in self.energy_holdup_slag_eqn.items():
            s = iscale.get_scaling_factor(
                self.energy_holdup_slag[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.energy_holdup_metal_eqn.items():
            s = iscale.get_scaling_factor(
                self.energy_holdup_metal[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, s, overwrite=False)
