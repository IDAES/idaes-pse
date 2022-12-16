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
The 1D drum model consists of three main unit operations
1) a flash model (only for water and steam)
2) a mixer model
3) a horizontal cylindrical water tank model
with multiple downcomers as water outlet

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
control volume of the cylindrical water tank model,
which computes velocity and pressure
drop. The exit of the drum model is the liquid outlet.

main assumptions:

1.- Heat loss to ambient is calculated based on natural heat convection
2.- Calculate pressure change due to gravity based on water level
and contraction to downcomer
3.- Water level is either fixed for steady-state model or calculated for
dynamic model
4.- Assume enthalpy_in == enthalpy_out + heat loss
5.- Subcooled water from economizer and saturated water from waterwall
are mixed before entering drum
6.- Heat conduction through drum wall thickness is modeled
by 1D heat conduction

Created: October 27 2020
"""
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

# Import Pyomo libraries
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from idaes.core.util.config import is_physical_parameter_block

# Additional import for the unit operation
from pyomo.environ import (
    value,
    Var,
    Param,
    Reference,
    exp,
    asin,
    cos,
    log10,
    log,
    sqrt,
    Constraint,
    TransformationFactory,
)

# from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from pyomo.network import Port
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from pyomo.network import Arc

from idaes.models_extra.power_generation.unit_models.helm.phase_separator import (
    HelmPhaseSeparator,
)
from idaes.models_extra.power_generation.unit_models.helm import (
    HelmMixer,
    MomentumMixingType,
)
from idaes.core.util.constants import Constants as const

__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"
__version__ = "2.0.0"


@declare_process_block_class("Drum1D")
class Drum1DData(UnitModelBlockData):
    """
    1-D Boiler Drum Unit Operation Class
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
            default=True,
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
        "finite_elements",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of finite elements in length (drum radius) domain",
            doc="""Number of finite elements to use when discretizing length
domain (default=5). Should set to the number of tube rows""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=3,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)""",
        ),
    )
    CONFIG.declare(
        "drum_inner_diameter",
        ConfigValue(
            default=1.0,
            description="inside diameter of drum",
            doc="define inside diameter of drum",
        ),
    )

    CONFIG.declare(
        "drum_thickness",
        ConfigValue(
            default=0.1,
            description="drum wall thickness",
            doc="define drum wall thickness",
        ),
    )

    def build(self):
        """
        Begin building model
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
            momentum_mixing_type=MomentumMixingType.equality,
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

        self.stream_flash_out = Arc(
            source=self.flash.liq_outlet, destination=self.mixer.SaturatedWater
        )

        # Pyomo arc connect flash liq_outlet with mixer SaturatedWater inlet
        TransformationFactory("network.expand_arcs").apply_to(self)

        # connect internal units (Mixer to Water Tank Model)
        # Mixer Outlet (mixed_state) to unit control volume.properties_in
        @self.Constraint(self.flowsheet().time)
        def connection_material_balance(b, t):
            return (
                b.mixer.mixed_state[t].flow_mol
                == b.control_volume.properties_in[t].flow_mol
            )

        @self.Constraint(self.flowsheet().time)
        def connection_enthalpy_balance(b, t):
            return (
                b.mixer.mixed_state[t].enth_mol
                == b.control_volume.properties_in[t].enth_mol
            )

        @self.Constraint(self.flowsheet().time)
        def connection_pressure_balance(b, t):
            return (
                b.mixer.mixed_state[t].pressure
                == b.control_volume.properties_in[t].pressure
            )

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
        Define the Geometry of the Unit
        """
        # Inside diameter of drum
        self.drum_diameter = Param(
            initialize=self.config.drum_inner_diameter, doc="Inside Diameter of Drum"
        )

        # Thickness of drum wall
        self.drum_thickness = Param(
            initialize=self.config.drum_thickness, doc="Wall Thickness of Drum"
        )

        # Thickness of insulation layer
        self.insulation_thickness = Var(
            initialize=0.2, doc="Insulation Layer Thickness"
        )

        # Length of drum
        self.drum_length = Var(initialize=10, doc="Horizontal Length of Drum")

        # Number of downcomers connected at the bottom of drum
        self.number_downcomer = Var(
            initialize=4, doc="Number of Downcomers Connected to Drum"
        )

        # Inside diameter of downcomer
        self.downcomer_diameter = Var(
            initialize=0.6, doc="Inside Diameter of Downcomer"
        )

        # Inside Radius expression
        @self.Expression(doc="Inside Radius of Drum")
        def drum_ri(b):
            return 0.5 * b.drum_diameter

        # Outside Radius expression
        @self.Expression(doc="Outside radius of drum")
        def drum_ro(b):
            return b.drum_ri + b.drum_thickness

        # Outside diameter expression
        @self.Expression(doc="Outside Radius of Drum")
        def drum_do(b):
            return b.drum_diameter + 2 * b.drum_thickness

        # Inner surface area (ignore two hemispheres at ends)
        @self.Expression(doc="Inner Surface Area")
        def drum_area(b):
            return const.pi * b.drum_diameter * b.drum_length

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.
        """

        # Thermal conductivity of drum
        self.therm_cond_metal = Param(initialize=40, mutable=True)

        # Thermal conductivity of insulation
        self.therm_cond_insulation = Param(initialize=0.08, mutable=True)

        # Thermal conductivity of air
        self.therm_cond_air = Param(initialize=0.03, mutable=True)

        # Density of drum material
        self.dens_metal = Param(initialize=7753, mutable=True)

        # Heat capacity of drum material
        self.cp_metal = Param(initialize=486, mutable=True)

        # Young modulus
        self.Young_modulus = Param(initialize=2.07e5, mutable=True)

        # Poisson's ratio
        self.Poisson_ratio = Param(initialize=0.292, mutable=True)

        # Coefficient of thermal expansion
        self.coefficient_therm_expansion = Param(initialize=1.4e-5, mutable=True)

        # Constant related to Rayleigh number for free convection
        # (gravity*expansion_coefficient*density
        # /viscosity/thermal_diffusivity)^(1/6)
        # Use properties at 50 C and 1 atm, 6.84e7^(1/6)=20.223
        self.const_Ra_root6 = Param(
            initialize=20.223, mutable=True, doc="Rayleigh number for free convection"
        )

        # Constant related to Nu for free convection,
        # 0.387/(1+0.721*Pr^(-9/16))^(8/27)
        # Use properties at 50 C and 1 atm
        self.const_Nu = Param(
            initialize=0.322,
            mutable=True,
            doc="constant related to Nu number for free" "convection",
        )

        # Ambient pressure
        self.pressure_amb = Param(initialize=1.0e5, doc="Ambient Pressure")

        # Ambient temperature
        self.temperature_ambient = Var(
            self.flowsheet().time, initialize=298.15, doc="Ambient Temperature"
        )

        # Inside heat transfer coefficient
        self.heat_transfer_in = Var(
            self.flowsheet().time, initialize=1, doc="Inside Heat Transfer Coefficient"
        )

        # Outside heat transfer coefficient
        self.heat_transfer_out = Var(
            self.flowsheet().time, initialize=1, doc="Outside Heat Transfer Coefficient"
        )

        # Insulation free convection heat transfer coefficient
        self.heat_transfer_free_conv = Var(
            self.flowsheet().time,
            initialize=1,
            doc="Insulation Free Convection" "Heat Transfer Coefficient",
        )

        # Ra number of free convection
        self.N_Ra_root6 = Var(
            self.flowsheet().time,
            initialize=80,
            doc="1/6 Power of Ra" "Number of Free Convection of Air",
        )

        # Nu number  of free convection
        self.N_Nu = Var(
            self.flowsheet().time,
            initialize=1,
            doc="Nu Number of Free Convection of Air",
        )

        # Define the continuous domains for model
        self.radial_domain = ContinuousSet(bounds=(self.drum_ri, self.drum_ro))

        # Temperature across wall thickness
        self.drum_wall_temperature = Var(
            self.flowsheet().time, self.radial_domain, bounds=(280, 800), initialize=550
        )

        # Declare derivatives in the model
        if self.config.dynamic is True:
            self.dTdt = DerivativeVar(
                self.drum_wall_temperature, wrt=self.flowsheet().time
            )
        self.dTdr = DerivativeVar(self.drum_wall_temperature, wrt=self.radial_domain)
        self.d2Tdr2 = DerivativeVar(
            self.drum_wall_temperature, wrt=(self.radial_domain, self.radial_domain)
        )

        discretizer = TransformationFactory("dae.finite_difference")
        discretizer.apply_to(
            self,
            nfe=self.config.finite_elements,
            wrt=self.radial_domain,
            scheme="CENTRAL",
        )

        # Add performance variables
        self.level = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Water Level from the Bottom of the Drum",
        )

        # Velocity of fluid inside downcomer pipe
        self.velocity_downcomer = Var(
            self.flowsheet().time,
            initialize=10.0,
            doc="Liquid Water Velocity at the Top of Downcomer",
        )

        # Pressure change due to contraction
        self.deltaP_contraction = Var(
            self.flowsheet().time,
            initialize=-1.0,
            doc="Pressure Change due to Contraction",
        )

        # Pressure change due to gravity
        self.deltaP_gravity = Var(
            self.flowsheet().time, initialize=1.0, doc="Pressure Change due to Gravity"
        )

        # Thermal diffusivity of drum
        @self.Expression(doc="Thermal Diffusivity of Drum Material")
        def diff_therm_metal(b):
            return b.therm_cond_metal / (b.dens_metal * b.cp_metal)

        # Expressure for the angle from the drum center
        # to the circumference point at water level
        @self.Expression(self.flowsheet().time, doc="Angle of Water Level")
        def alpha_drum(b, t):
            return asin((b.level[t] - b.drum_ri) / b.drum_ri)

        # Expressure for the fraction of wet area
        @self.Expression(self.flowsheet().time, doc="Fraction of Wet Area")
        def frac_wet_area(b, t):
            return (b.alpha_drum[t] + const.pi / 2) / const.pi

        # Constraint for volume liquid in drum
        @self.Constraint(self.flowsheet().time, doc="Volume of Liquid in Drum")
        def volume_eqn(b, t):
            return (
                b.volume[t]
                == (
                    (b.alpha_drum[t] + 0.5 * const.pi) * b.drum_ri**2
                    + b.drum_ri * cos(b.alpha_drum[t]) * (b.level[t] - b.drum_ri)
                )
                * b.drum_length
            )

        # Equation for velocity at the entrance of downcomer
        @self.Constraint(self.flowsheet().time, doc="Velocity at Entrance of Downcomer")
        def velocity_eqn(b, t):
            return (
                b.velocity_downcomer[t]
                * 0.25
                * const.pi
                * b.downcomer_diameter**2
                * b.number_downcomer
                == b.control_volume.properties_out[t].flow_vol
            )

        # Pressure change equation for contraction,
        # -0.5*1/2*density*velocity^2 for stagnation head loss
        # plus 1/2*density*velocity^2 dynamic head
        # (acceleration pressure change)
        @self.Constraint(
            self.flowsheet().time, doc="Pressure Change due To Contraction"
        )
        def pressure_change_contraction_eqn(b, t):
            return (
                b.deltaP_contraction[t]
                == -0.75
                * b.control_volume.properties_out[t].dens_mass_phase["Liq"]
                * b.velocity_downcomer[t] ** 2
            )

        # Pressure change equation for gravity, density*gravity*height
        @self.Constraint(self.flowsheet().time, doc="Pressure Change due to Gravity")
        def pressure_change_gravity_eqn(b, t):
            return (
                b.deltaP_gravity[t]
                == b.control_volume.properties_out[t].dens_mass_phase["Liq"]
                * const.acceleration_gravity
                * b.level[t]
            )

        # Total pressure change equation
        @self.Constraint(self.flowsheet().time, doc="Pressure Drop")
        def pressure_change_total_eqn(b, t):
            return b.deltaP[t] == b.deltaP_contraction[t] + b.deltaP_gravity[t]

        # Constraint for heat conduction equation
        @self.Constraint(
            self.flowsheet().time,
            self.radial_domain,
            doc="1-D Heat Conduction Equation Through Radius",
        )
        def heat_conduction_eqn(b, t, r):
            if r == b.radial_domain.first() or r == b.radial_domain.last():
                return Constraint.Skip
            if self.config.dynamic is True:
                return (
                    b.dTdt[t, r]
                    == b.diff_therm_metal * b.d2Tdr2[t, r]
                    + b.diff_therm_metal * (1 / r) * b.dTdr[t, r]
                )
            else:
                return (
                    0
                    == b.diff_therm_metal * b.d2Tdr2[t, r]
                    + b.diff_therm_metal * (1 / r) * b.dTdr[t, r]
                )

        @self.Constraint(self.flowsheet().time, doc="Inner Wall Boundary")
        def inner_wall_bc_eqn(b, t):
            return (
                b.heat_transfer_in[t]
                * (
                    b.control_volume.properties_out[t].temperature
                    - b.drum_wall_temperature[t, b.radial_domain.first()]
                )
                == -b.dTdr[t, b.radial_domain.first()] * b.therm_cond_metal
            )

        @self.Constraint(self.flowsheet().time, doc="Outer Wall Boundary")
        def outer_wall_bc_eqn(b, t):
            return (
                b.heat_transfer_out[t]
                * (
                    b.drum_wall_temperature[t, b.radial_domain.last()]
                    - b.temperature_ambient[t]
                )
                == -b.dTdr[t, b.radial_domain.last()] * b.therm_cond_metal
            )

        # Inner wall BC for dTdt
        @self.Constraint(
            self.flowsheet().time, doc="Extra Inner Wall Temperature Derivative"
        )
        def extra_at_inner_wall_eqn(b, t):
            if self.config.dynamic is True:
                term = b.dTdt[t, b.radial_domain.first()]
            else:
                term = 0
            return term == 4 * b.diff_therm_metal * (
                b.radial_domain.first() + b.radial_domain.at(2)
            ) / (b.radial_domain.at(2) - b.radial_domain.first()) ** 2 / (
                3 * b.radial_domain.first() + b.radial_domain.at(2)
            ) * (
                b.drum_wall_temperature[t, b.radial_domain.at(2)]
                - b.drum_wall_temperature[t, b.radial_domain.first()]
            ) + 8 * b.diff_therm_metal / b.therm_cond_metal * b.heat_transfer_in[
                t
            ] * b.radial_domain.first() / (
                b.radial_domain.at(2) - b.radial_domain.first()
            ) / (
                3 * b.radial_domain.first() + b.radial_domain.at(2)
            ) * (
                b.control_volume.properties_out[t].temperature
                - b.drum_wall_temperature[t, b.radial_domain.first()]
            )

        @self.Constraint(
            self.flowsheet().time, doc="Extra Outer Wall Temperature Derivative"
        )
        def extra_at_outer_wall_eqn(b, t):
            if self.config.dynamic is True:
                term = b.dTdt[t, b.radial_domain.last()]
            else:
                term = 0
            return term == 4 * b.diff_therm_metal * (
                b.radial_domain.last() + b.radial_domain.at(-2)
            ) / (b.radial_domain.last() - b.radial_domain.at(-2)) ** 2 / (
                3 * b.radial_domain.last() + b.radial_domain.at(-2)
            ) * (
                b.drum_wall_temperature[t, b.radial_domain.at(-2)]
                - b.drum_wall_temperature[t, b.radial_domain.last()]
            ) + 8 * b.diff_therm_metal / b.therm_cond_metal * b.heat_transfer_out[
                t
            ] * b.radial_domain.last() / (
                b.radial_domain.last() - b.radial_domain.at(-2)
            ) / (
                3 * b.radial_domain.last() + b.radial_domain.at(-2)
            ) * (
                b.temperature_ambient[t]
                - b.drum_wall_temperature[t, b.radial_domain.last()]
            )

        # Reduced pressure expression
        @self.Expression(self.flowsheet().time, doc="Reduced Pressure")
        def pres_reduced(b, t):
            return b.control_volume.properties_out[t].pressure / 2.2048e7

        # Calculate inner side heat transfer coefficient
        # with minimum temperature difference set to sqrt(0.1)
        # multipling wet area fraction to convert it
        # to the value based on total circumference
        @self.Constraint(
            self.flowsheet().time, doc="Inner Side Heat Transfer Coefficient"
        )
        def h_in_eqn(b, t):
            return (
                b.heat_transfer_in[t]
                == 2178.6
                * (b.control_volume.properties_out[t].pressure / 2.2048e7) ** 0.36
                / (-log10(b.control_volume.properties_out[t].pressure / 2.2048e7))
                ** 1.65
                * (
                    0.1
                    + (
                        b.control_volume.properties_out[t].temperature
                        - b.drum_wall_temperature[t, b.radial_domain.first()]
                    )
                    ** 2
                )
                * b.frac_wet_area[t]
            )

        # Expressure for insulation heat transfer (conduction)
        # resistance based on drum metal outside diameter
        @self.Expression(doc="Heat Transfer Resistance of Insulation Layer")
        def resistance_insulation(b):
            return (
                b.drum_ro
                * log((b.drum_ro + b.insulation_thickness) / b.drum_ro)
                / b.therm_cond_insulation
            )

        # heat_transfer_out equation considering conduction through insulation
        # and free convection between insulation and ambient
        @self.Constraint(
            self.flowsheet().time, doc="Outer Side Heat Transfer Coefficient"
        )
        def h_out_eqn(b, t):
            return (
                b.heat_transfer_out[t]
                * (b.resistance_insulation + 1 / b.heat_transfer_free_conv[t])
                == 1.0
            )

        # Expressure for outside insulation wall temperature (skin temperature)
        @self.Expression(
            self.flowsheet().time, doc="Outside Insulation Wall Temperature"
        )
        def temp_insulation_outside(b, t):
            return (
                b.temperature_ambient[t]
                + (
                    b.drum_wall_temperature[t, b.radial_domain.last()]
                    - b.temperature_ambient[t]
                )
                * b.heat_transfer_out[t]
                / b.heat_transfer_free_conv[t]
            )

        # Ra number equation
        @self.Constraint(self.flowsheet().time, doc="Ra Number of Free Convection")
        def Ra_number_eqn(b, t):
            return (
                b.N_Ra_root6[t]
                == b.const_Ra_root6
                * sqrt(b.drum_do + 2 * b.insulation_thickness)
                * (b.temp_insulation_outside[t] - b.temperature_ambient[t]) ** 0.166667
            )

        # Nu number equation
        @self.Constraint(self.flowsheet().time, doc="Nu Number of Free Convection")
        def Nu_number_eqn(b, t):
            return b.N_Nu[t] == (0.6 + b.const_Nu * b.N_Ra_root6[t]) ** 2

        # Free convection coefficient based on the drum metal outside diameter
        @self.Constraint(
            self.flowsheet().time,
            doc="Free Convection Heat Transfer Coefficient"
            "between Insulation Wall and Ambient",
        )
        def h_free_conv_eqn(b, t):
            return (
                b.heat_transfer_free_conv[t] == b.N_Nu[t] * b.therm_cond_air / b.drum_do
            )

        @self.Constraint(self.flowsheet().time, doc="Heat Loss of Water")
        def heat_loss_eqn(b, t):
            return b.heat_duty[t] == b.drum_area * b.heat_transfer_in[t] * (
                b.drum_wall_temperature[t, b.radial_domain.first()]
                - b.control_volume.properties_out[t].temperature
            )

        # Calculate mechanical and thermal stresses based on
        # EN 13445 Standard
        # Integer indexing for radius domain
        self.rindex = Param(
            self.radial_domain,
            initialize=1,
            mutable=True,
            doc="Integer Indexing for Radius Domain",
        )

        # calculate integral point for mean temperature in the wall
        @self.Expression(self.flowsheet().time, doc="Mean Temperature across the Wall")
        def mean_temperature(b, t):
            return (
                2
                * (b.radial_domain.at(2) - b.radial_domain.at(1))
                / (b.drum_ro**2 - b.drum_ri**2)
                * (
                    sum(
                        0.5
                        * (
                            b.radial_domain.at(i - 1)
                            * b.drum_wall_temperature[t, b.radial_domain.at(i - 1)]
                            + b.radial_domain.at(i)
                            * b.drum_wall_temperature[t, b.radial_domain.at(i)]
                        )
                        for i in range(2, len(b.radial_domain) + 1)
                    )
                )
            )

        for index_r, value_r in enumerate(self.radial_domain, 1):
            self.rindex[value_r] = index_r

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Discrete Point Mean Temperature",
        )
        def discrete_mean_temperature(b, t, r):
            if b.rindex[r].value == 1:
                return b.drum_wall_temperature[t, b.radial_domain.first()]
            else:
                return (
                    2
                    * (b.radial_domain.at(2) - b.radial_domain.at(1))
                    / (b.radial_domain.at(b.rindex[r].value) ** 2 - b.drum_ri**2)
                    * (
                        sum(
                            0.5
                            * (
                                b.radial_domain.at(j - 1)
                                * b.drum_wall_temperature[t, b.radial_domain.at(j - 1)]
                                + b.radial_domain.at(j)
                                * b.drum_wall_temperature[t, b.radial_domain.at(j)]
                            )
                            for j in range(2, b.rindex[r].value + 1)
                        )
                    )
                )

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Thermal Stress at Radial Direction for Drum",
        )
        def therm_sigma_r(b, t, r):
            if r == b.radial_domain.first() or r == b.radial_domain.last():
                return 0
            else:
                return (
                    0.5
                    * b.Young_modulus
                    * b.coefficient_therm_expansion
                    / (1 - b.Poisson_ratio)
                    * (
                        (1 - b.drum_ri**2 / r**2)
                        * (b.mean_temperature[t] - b.discrete_mean_temperature[t, r])
                    )
                )

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Thermal Stress at Circumferential Direction" "for Drum",
        )
        def therm_sigma_theta(b, t, r):
            return (
                0.5
                * b.Young_modulus
                * b.coefficient_therm_expansion
                / (1 - b.Poisson_ratio)
                * (
                    (1 + b.drum_ri**2 / r**2) * b.mean_temperature[t]
                    + (1 - b.drum_ri**2 / r**2) * b.discrete_mean_temperature[t, r]
                    - 2 * b.drum_wall_temperature[t, r]
                )
            )

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Thermal Stress at Axial Direction for Drum",
        )
        def therm_sigma_z(b, t, r):
            return (
                b.Young_modulus
                * b.coefficient_therm_expansion
                / (1 - b.Poisson_ratio)
                * (b.mean_temperature[t] - b.drum_wall_temperature[t, r])
            )

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Mechanical Stress at Radial Direction for Drum",
        )
        def mech_sigma_r(b, t, r):
            if r == b.radial_domain.first():
                return 1e-6 * (-b.control_volume.properties_out[t].pressure)
            elif r == b.radial_domain.last():
                return 1e-6 * (-b.pressure_amb)
            else:
                return 0.1 * (
                    1e-5
                    * (
                        b.control_volume.properties_out[t].pressure * b.drum_ri**2
                        - b.pressure_amb * b.drum_ro**2
                    )
                    / (b.drum_ro**2 - b.drum_ri**2)
                    + (
                        1e-5
                        * (b.pressure_amb - b.control_volume.properties_out[t].pressure)
                        * b.drum_ri**2
                        * b.drum_ro**2
                        / (r**2 * (b.drum_ro**2 - b.drum_ri**2))
                    )
                )

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Mechanical Stress at Circumferential Direction" "for Drum",
        )
        def mech_sigma_theta(b, t, r):
            return 0.1 * (
                1e-5
                * (
                    b.control_volume.properties_out[t].pressure * b.drum_ri**2
                    - b.pressure_amb * b.drum_ro**2
                )
                / (b.drum_ro**2 - b.drum_ri**2)
                - (
                    1e-5
                    * (b.pressure_amb - b.control_volume.properties_out[t].pressure)
                    * b.drum_ri**2
                    * b.drum_ro**2
                    / (r**2 * (b.drum_ro**2 - b.drum_ri**2))
                )
            )

        @self.Expression(
            self.flowsheet().time, doc="Mechanical Stress at Axial Direction" "for Drum"
        )
        def mech_sigma_z(b, t):
            return 0.1 * (
                1e-5
                * (
                    b.control_volume.properties_out[t].pressure * b.drum_ri**2
                    - b.pressure_amb * b.drum_ro**2
                )
                / (b.drum_ro**2 - b.drum_ri**2)
            )

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Principal Structural Stress" "at Radial Direction for Drum",
        )
        def sigma_r(b, t, r):
            if r == b.radial_domain.first():
                return 1e-6 * (-b.control_volume.properties_out[t].pressure)
            elif r == b.radial_domain.last():
                return 1e-6 * (-b.pressure_amb)
            else:
                return b.mech_sigma_r[t, r] + b.therm_sigma_r[t, r]

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Principal Structural Stress" "at Circumferential Direction for Drum",
        )
        def sigma_theta(b, t, r):
            return b.mech_sigma_theta[t, r] + b.therm_sigma_theta[t, r]

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Principal Structural Stress" "at Axial Direction for Drum",
        )
        def sigma_z(b, t, r):
            return b.mech_sigma_z[t] + b.therm_sigma_z[t, r]

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Equivalent von Mises Stress for Drum",
        )
        def sigma_von_Mises(b, t, r):
            return sqrt(
                b.sigma_r[t, r] ** 2
                + b.sigma_theta[t, r] ** 2
                + b.sigma_z[t, r] ** 2
                - (
                    b.sigma_r[t, r] * b.sigma_theta[t, r]
                    + b.sigma_r[t, r] * b.sigma_z[t, r]
                    + b.sigma_theta[t, r] * b.sigma_z[t, r]
                )
            )

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Variation Principal Stress"
            "between Radial-Circumferential Directions for Drum",
        )
        def delta_sigma_r_theta(b, t, r):
            return abs(b.sigma_r[t, r] - b.sigma_theta[t, r])

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Variation Principal Stress"
            "between Circumferential-Axial Directions for Drum",
        )
        def delta_sigma_theta_z(b, t, r):
            return abs(b.sigma_theta[t, r] - b.sigma_z[t, r])

        @self.Expression(
            self.flowsheet().time,
            self.radial_domain,
            doc="Variation Principal Stress" "between Axial-Radial Directions for Drum",
        )
        def delta_sigma_z_r(b, t, r):
            return abs(b.sigma_z[t, r] - b.sigma_r[t, r])

        # Calculate stress based on EN12952 standard
        # -----------------------------------------------------------------
        # mechanical stress of circumferential direction
        # mean radius
        r_ms_drum = self.drum_ri + self.drum_thickness / 2

        # thickness of downcomer //m
        self.downcomer_thickness = Param(
            initialize=0.0254, mutable=True, doc="Thickness of Downcomer"
        )

        # mean diameter of downcomer
        # pipe_d = pipe_d(inner) + self.downcomer_thickness ==
        # pipe_d(outer)-self.downcomer_thickness
        pipe_d = self.downcomer_diameter + self.downcomer_thickness

        # mechanical coefficients
        k_m_A = (
            -1.14 * (self.downcomer_thickness / self.drum_thickness) ** 2
            - 0.89 * (self.downcomer_thickness / self.drum_thickness)
            + 1.43
        )
        k_m_B = (
            0.326 * (self.downcomer_thickness / self.drum_thickness) ** 2
            - 0.59 * (self.downcomer_thickness / self.drum_thickness)
            + 1.08
        )
        k_m_C = (
            pipe_d / (2 * r_ms_drum) * sqrt((2 * r_ms_drum) / (2 * self.drum_thickness))
        )
        k_m = 2.2 + exp(k_m_A) * k_m_C**k_m_B

        # thermal stress concentration factor
        k_t_A = pipe_d / (2 * r_ms_drum)
        # heat transfer coefficient:3000 for water , 1000 for steam
        k_t_B = 1000
        k_t = sqrt(
            (
                2
                - (k_t_B + 2700) / (k_t_B + 1700) * k_t_A
                + k_t_B / (k_t_B + 1700) * (exp(-7 * k_t_A) - 1)
            )
            ** 2
            + 0.81 * k_t_A**2
        )

        # mechanical stress at circumferential direction
        @self.Expression(
            self.flowsheet().time,
            doc="Mechanical Stress "
            "at Circumferential Direction"
            "for Drum (EN 12952-3)",
        )
        def sigma_p(b, t):
            return (
                0.1
                * 1e-5
                * b.control_volume.properties_out[t].pressure
                * r_ms_drum
                / self.drum_thickness
            )

        # thermal stress at circumferential direction
        @self.Expression(
            self.flowsheet().time,
            doc="Thermal Stress at Circumferential" "Direction for Drum (EN 12952-3)",
        )
        def sigma_t(b, t):
            delta_T = (
                b.mean_temperature[t]
                - b.drum_wall_temperature[t, b.radial_domain.first()]
            )
            return (
                b.coefficient_therm_expansion
                * b.Young_modulus
                / (1 - b.Poisson_ratio)
                * delta_T
            )

        # calculate stress at 2 locations at the hole
        # stress at crotch corner P1 and location P2

        # mechanical stress by pressure at crotch corner
        @self.Expression(
            self.flowsheet().time, doc="Mechanical Stress at Crotch Corner" "for Drum"
        )
        def sigma_p_P1(b, t):
            return b.sigma_p[t] * k_m

        # mechanical stress at location P2
        @self.Expression(
            self.flowsheet().time,
            doc="Mechanical Stress at" "Critical Point P2 for Drum",
        )
        def sigma_p_P2(b, t):
            return b.sigma_p[t] * k_m / 5

        # thermal stress at crotch corner
        @self.Expression(
            self.flowsheet().time, doc="Thermal Stress at Crotch Corner for Drum"
        )
        def sigma_t_P1(b, t):
            return b.sigma_t[t] * k_t

        # thermal stress at location P2
        @self.Expression(
            self.flowsheet().time, doc="Thermal Stress" "at Critical Point P2 for Drum"
        )
        def sigma_t_P2(b, t):
            return b.sigma_t[t] * k_t

        # total circumferential stress with notch effect
        # crotch corner P1
        @self.Expression(
            self.flowsheet().time,
            doc="Circumferential Stress" "at Crotch Corner for Drum",
        )
        def sigma_theta_P1(b, t):
            return b.sigma_p_P1[t] + b.sigma_t_P1[t]

        # location P2
        @self.Expression(
            self.flowsheet().time,
            doc="Circumferential Stress" "at Critical Point P2 for Drum",
        )
        def sigma_theta_P2(b, t):
            return b.sigma_p_P2[t] + b.sigma_t_P2[t]

        # total stress with notch effect // f1 - f2
        # crotch corner
        @self.Expression(
            self.flowsheet().time, doc="Total Stress at Crotch Corner for Drum"
        )
        def sigma_notch_P1(b, t):
            return (
                b.sigma_theta_P1[t] + 1e-6 * b.control_volume.properties_out[t].pressure
            )

        # location P2
        @self.Expression(
            self.flowsheet().time, doc="Total Stress at Critial Point P2 for Header"
        )
        def sigma_notch_P2(b, t):
            return (
                b.sigma_theta_P2[t] + 1e-6 * b.control_volume.properties_out[t].pressure
            )

        # Von Mises equivalent stress
        # VM stress at crotch corner
        @self.Expression(
            self.flowsheet().time,
            doc="Equivalent von Mises Stress" "at Crotch Corner for Drum",
        )
        def sigma_eff_P1(b, t):
            sigma_comer_P1 = b.sigma_theta_P1[t]
            p_in = b.control_volume.properties_out[t].pressure
            return sqrt(
                sigma_comer_P1**2
                + (-1e-6 * p_in) ** 2
                + (-1e-6 * p_in) ** 2
                - (
                    (-1e-6 * p_in) * sigma_comer_P1
                    + (-1e-6 * p_in) * sigma_comer_P1
                    + (-1e-6 * p_in) * (-1e-6 * p_in)
                )
            )

        # VM stress at location P2
        @self.Expression(
            self.flowsheet().time,
            doc="Equivalent von Mises Stress" "at Critical Point P2 for Drum",
        )
        def sigma_eff_P2(b, t):
            sigma_comer_P2 = b.sigma_theta_P2[t]
            p_in = b.control_volume.properties_out[t].pressure
            return sqrt(
                sigma_comer_P2**2
                + (-1e-6 * p_in) ** 2
                + (-1e-6 * p_in) ** 2
                - (
                    (-1e-6 * p_in) * sigma_comer_P2
                    + (-1e-6 * p_in) * sigma_comer_P2
                    + (-1e-6 * p_in) * (-1e-6 * p_in)
                )
            )

    def set_initial_condition(self):
        if self.config.dynamic is True:
            self.control_volume.material_accumulation[:, :, :].value = 0
            self.control_volume.energy_accumulation[:, :].value = 0
            self.dTdt[:, :].value = 0
            self.control_volume.material_accumulation[0, :, :].fix(0)
            self.control_volume.energy_accumulation[0, :].fix(0)
            self.dTdt[0, :].fix(0)

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
        state_args_feedwater : a dict of arguments to be passed to the property
        package(s) for the control_volume of the model to
        provide an initial state for initialization
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

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        init_log.info_low("Starting Initialization...")
        # fix FeedWater Inlet
        flags_fw = fix_state_vars(blk.mixer.FeedWater_state, state_args_feedwater)

        blk.flash.initialize(
            state_args_water_steam=state_args_water_steam,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )
        init_log.info("Initialization Step 1 Complete. Flash Model Solved")

        # Copy the state values from flash model liquid outlet to
        # mixer model saturated water inlet
        blk.mixer.SaturatedWater.flow_mol[:].value = blk.flash.liq_outlet.flow_mol[
            0
        ].value
        blk.mixer.SaturatedWater.pressure[:].value = blk.flash.liq_outlet.pressure[
            0
        ].value
        blk.mixer.SaturatedWater.enth_mol[:].value = blk.flash.liq_outlet.enth_mol[
            0
        ].value

        # Since the mixer model requires the same inlet pressure,
        # unfix feed water pressure
        # use the water_steam inlet pressure as mixer outlet pressure
        blk.mixer.FeedWater.pressure.unfix()

        blk.mixer.initialize(outlvl=outlvl, optarg=optarg, solver=solver)
        init_log.info("Initialization Step 2 Complete. Mixer Model Solved")

        # Copy the state values from mixer model outlet to control volume inlet
        blk.control_volume.properties_in[:].flow_mol.value = blk.mixer.outlet.flow_mol[
            0
        ].value
        blk.control_volume.properties_in[:].enth_mol.value = blk.mixer.outlet.enth_mol[
            0
        ].value
        blk.control_volume.properties_in[:].pressure.value = blk.mixer.outlet.pressure[
            0
        ].value

        blk.control_volume.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, hold_state=False
        )
        init_log.info("Initialization Step 3 Complete." "Control Volume Initialized")

        # fix flash model inlet
        flags_steam = fix_state_vars(blk.flash.mixed_state, state_args_water_steam)

        # set initial values for T
        r_mid = value((blk.radial_domain.first() + blk.radial_domain.last()) / 2)
        # assume outside wall temperature is 1 K lower than fluid temperature
        T_out = value(blk.control_volume.properties_in[0].temperature - 1)
        T_mid = value((T_out + blk.control_volume.properties_in[0].temperature) / 2)
        slope = value(
            (T_out - blk.control_volume.properties_in[0].temperature)
            / (blk.radial_domain.last() - blk.radial_domain.first())
            / 3
        )
        for x in blk.radial_domain:
            blk.drum_wall_temperature[:, x].fix(T_mid + slope * (x - r_mid))
        blk.drum_wall_temperature[:, :].unfix()

        # Fix outlet enthalpy and pressure
        for t in blk.flowsheet().time:
            blk.control_volume.properties_out[t].pressure.fix(
                value(blk.control_volume.properties_in[0].pressure) - 5000.0
            )
            blk.control_volume.properties_out[t].enth_mol.fix(
                value(blk.control_volume.properties_in[0].enth_mol)
            )
        blk.pressure_change_total_eqn.deactivate()
        blk.heat_loss_eqn.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
            "Initialization Step 4 Complete {}.".format(idaeslog.condition(res))
        )

        # Unfix outlet enthalpy and pressure
        for t in blk.flowsheet().time:
            blk.control_volume.properties_out[t].pressure.unfix()
            blk.control_volume.properties_out[t].enth_mol.unfix()
        blk.pressure_change_total_eqn.activate()
        blk.heat_loss_eqn.activate()

        # Finally solve model with all constraints activated
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 5 {}.".format(idaeslog.condition(res)))
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

        for t, c in self.connection_material_balance.items():
            sf = iscale.get_scaling_factor(1e-4, default=1)
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for t, c in self.connection_enthalpy_balance.items():
            sf = iscale.get_scaling_factor(1e-4, default=1)
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for t, c in self.connection_pressure_balance.items():
            sf = iscale.get_scaling_factor(1e-6, default=1)
            iscale.constraint_scaling_transform(c, sf, overwrite=False)
