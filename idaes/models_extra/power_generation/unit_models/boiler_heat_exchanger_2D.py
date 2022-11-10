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
The boiler 2D heat exchanger model consist of a cross flow shell and tube
heat exchanger. 1-D Cross Flow Heat Exchanger Model with wall temperatures,
discretization based on tube rows


The model includes shell and tube rigorous heat transfer calculations and
pressure drop calculations for shell side. Note that this model assumes no
phase transitions (if user requires phase transitions, they need a general
model)

"""
# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    Constraint,
    TransformationFactory,
    Reference,
    value,
    exp,
    sqrt,
    log,
    log10,
    sin,
    cos,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    ControlVolume1DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
    UnitModelBlockData,
)
from pyomo.dae import ContinuousSet, DerivativeVar
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.constants import Constants as const
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.math import smooth_max
import idaes.logger as idaeslog

__author__ = "Jinliang Ma, Q. M. Le, M. Zamarripa "

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("HeatExchangerCrossFlow2D_Header")
class HeatExchangerCrossFlow2D_HeaderData(UnitModelBlockData):
    """Standard Heat Exchanger Cross Flow Unit Model Class."""

    CONFIG = UnitModelBlockData.CONFIG()

    # Template for config arguments for shell and tube side
    _SideTemplate = ConfigBlock()
    _SideTemplate.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.componentTotal,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentTotal.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    _SideTemplate.declare(
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
    _SideTemplate.declare(
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
    _SideTemplate.declare(
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
    _SideTemplate.declare(
        "property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object""",
        ),
    )
    _SideTemplate.declare(
        "property_package_args",
        ConfigValue(
            default={},
            description="Arguments for constructing shell property package",
            doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)""",
        ),
    )

    # Create individual config blocks for shell and tube side
    CONFIG.declare("shell_side", _SideTemplate(doc="shell side config arguments"))
    CONFIG.declare("tube_side", _SideTemplate(doc="tube side config arguments"))

    # Common config args for both sides
    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default="dae.finite_difference",
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See Pyomo
documentation for supported transformations.""",
        ),
    )
    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default="BACKWARD",
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transformating domain.
See Pyomo documentation for supported schemes.""",
        ),
    )
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of finite elements length domain",
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
        "flow_type",
        ConfigValue(
            default="co_current",
            domain=In(["co_current", "counter_current"]),
            description="Flow configuration of heat exchanger",
            doc="""Flow configuration of heat exchanger
co_current: shell and tube flows from 0 to 1
counter_current: shell side flows from 0 to 1
tube side flows from 1 to 0""",
        ),
    )
    CONFIG.declare(
        "tube_arrangement",
        ConfigValue(
            default="in-line",
            domain=In(["in-line", "staggered"]),
            description="tube configuration",
            doc="Tube arrangement could be in-line or staggered",
        ),
    )
    CONFIG.declare(
        "tube_side_water_phase",
        ConfigValue(
            default="Liq",
            domain=In(["Liq", "Vap"]),
            description="tube side water phase",
            doc="Define water phase for property calls",
        ),
    )
    CONFIG.declare(
        "has_radiation",
        ConfigValue(
            default=False,
            domain=In([False, True]),
            description="Has side 2 gas radiation",
            doc="Define if shell side gas radiation is to be considered",
        ),
    )
    CONFIG.declare(
        "tube_inner_diameter",
        ConfigValue(
            default=None,
            description="Inner diameter of tube",
            doc="User must define inner diameter of tube",
        ),
    )
    CONFIG.declare(
        "tube_thickness",
        ConfigValue(
            default=None,
            description="Tube wall thickness",
            doc="User must define tube wall thickness",
        ),
    )
    CONFIG.declare(
        "radial_elements",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of finite elements in radius domain",
            doc="""Number of finite elements to use when discretizing radius
        domain (default=5).""",
        ),
    )
    CONFIG.declare(
        "header_inner_diameter",
        ConfigValue(
            default=None,
            description="Inner diameter of header",
            doc="User must define inner diameter of header",
        ),
    )
    CONFIG.declare(
        "header_wall_thickness",
        ConfigValue(
            default=None,
            description="Header wall thickness",
            doc="User must define header wall thickness",
        ),
    )
    CONFIG.declare(
        "header_radial_elements",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of finite elements in radius domain",
            doc="""Number of finite elements to use when discretizing radius
        domain (default=5).""",
        ),
    )
    CONFIG.declare(
        "has_header",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Flag to include tube header",
            doc="""If has_header is True, user must provide header thickness and
        inner diameter.""",
        ),
    )

    def build(self):
        """
        Begin building model.

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(HeatExchangerCrossFlow2D_HeaderData, self).build()

        # Set flow directions for the control volume blocks and specify
        # dicretisation if not specified.
        if self.config.flow_type == "co_current":
            set_direction_shell = FlowDirection.forward
            set_direction_tube = FlowDirection.forward
        else:
            set_direction_shell = FlowDirection.forward
            set_direction_tube = FlowDirection.backward

        # Control volume 1D for shell and tube, set to steady-state
        # for fluid on both sides
        self.shell = ControlVolume1DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.shell_side.property_package,
            property_package_args=self.config.shell_side.property_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )

        self.tube = ControlVolume1DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.tube_side.property_package,
            property_package_args=self.config.tube_side.property_package_args,
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )

        self.shell.add_geometry(flow_direction=set_direction_shell)
        self.tube.add_geometry(flow_direction=set_direction_tube)

        self.shell.add_state_blocks(
            information_flow=set_direction_shell, has_phase_equilibrium=False
        )
        self.tube.add_state_blocks(
            information_flow=set_direction_tube, has_phase_equilibrium=False
        )

        # Populate shell
        self.shell.add_material_balances(
            balance_type=self.config.shell_side.material_balance_type,
            has_phase_equilibrium=False,
        )

        self.shell.add_energy_balances(
            balance_type=self.config.shell_side.energy_balance_type,
            has_heat_transfer=True,
        )

        self.shell.add_momentum_balances(
            balance_type=self.config.shell_side.momentum_balance_type,
            has_pressure_change=self.config.shell_side.has_pressure_change,
        )

        self.shell.apply_transformation()

        # Populate tube
        self.tube.add_material_balances(
            balance_type=self.config.tube_side.material_balance_type,
            has_phase_equilibrium=False,
        )

        self.tube.add_energy_balances(
            balance_type=self.config.tube_side.energy_balance_type,
            has_heat_transfer=True,
        )

        self.tube.add_momentum_balances(
            balance_type=self.config.tube_side.momentum_balance_type,
            has_pressure_change=self.config.tube_side.has_pressure_change,
        )

        self.tube.apply_transformation()

        # Add Ports for shell side
        self.add_inlet_port(name="shell_inlet", block=self.shell)
        self.add_outlet_port(name="shell_outlet", block=self.shell)

        # Add Ports for tube side
        self.add_inlet_port(name="tube_inlet", block=self.tube)
        self.add_outlet_port(name="tube_outlet", block=self.tube)

        # Check input arguments
        # tube inputs:
        if self.config.tube_inner_diameter is None:
            raise ConfigurationError(
                "User must provide a value for " "tube_inner_diameter"
            )
        if self.config.tube_thickness is None:
            raise ConfigurationError("User must provide a value for " "tube_thickness")
        # header inputs:
        if (
            self.config.has_header is False
            and self.config.header_inner_diameter is True
        ):
            _log.info_high(
                "User set has_header to False " "and provided header_inner_diameter"
            )

        if self.config.has_header and self.config.header_inner_diameter is None:
            raise ConfigurationError(
                "If has_heder is True, user must " "provide header_inner_diameter"
            )

        if (
            self.config.has_header is False
            and self.config.header_wall_thickness is True
        ):
            _log.info_high(
                "User set has_header to False " "and provided header_wall_thickness"
            )

        if self.config.has_header and self.config.header_wall_thickness is None:
            raise ConfigurationError(
                "If has_heder is True, user must " "provide header_wall_thickness"
            )

        self._make_geometry()

        self._make_performance()

    def _make_geometry(self):
        """
        Constraints for Unit Model.

        Args:
            None

        Returns:
            None
        """
        # Add object reference to control volume geometry
        add_object_reference(self, "area_flow_shell", self.shell.area)
        add_object_reference(self, "length_flow_shell", self.shell.length)
        add_object_reference(self, "area_flow_tube", self.tube.area)
        add_object_reference(self, "length_flow_tube", self.tube.length)

        # Elevation difference (outlet - inlet) for static pressure calculation
        self.delta_elevation = Var(
            initialize=0.0,
            doc="Elevation Increase Used for" "Static Pressure Calculation",
        )

        # Number of tube columns in the cross section plane perpendicular
        # to shell side fluid flow (y direction)
        self.tube_ncol = Var(initialize=10.0, doc="Number of Tube Columns")

        # Number of segments of tube bundles
        self.tube_nseg = Var(initialize=10.0, doc="Number of Tube Segments")

        # Number of inlet tube rows
        self.tube_inlet_nrow = Var(initialize=1, doc="Number of Inlet Tube Rows")

        # DAE discretization scaling (tube_radius_scaling)
        self.ri_scaling = Param(
            initialize=0.01,
            mutable=True,
            doc="Discretization scaling Tube inner radius",
        )

        # Inner diameter of tubes
        self.tube_di = Param(
            initialize=self.config.tube_inner_diameter,
            mutable=False,
            doc="Inner diameter of tube",
        )

        # Thickness of header
        self.tube_thickness = Param(
            initialize=self.config.tube_thickness, mutable=False, doc="Tube thickness"
        )
        if self.config.has_header is True:
            self.head_di = Param(
                initialize=self.config.header_inner_diameter,
                mutable=False,
                doc="Inner Diameter of Tube",
            )

            # Thickness of header
            self.head_thickness = Param(
                initialize=self.config.header_wall_thickness,
                mutable=False,
                doc="Header wall thickness",
            )

            self.head_ri_scaling = Param(
                initialize=0.1,
                mutable=True,
                doc="Discretization scaling" "of the header inner radius",
            )

        # Pitch of tubes between two neighboring columns (in y direction).
        # Always greater than tube outer diameter
        self.pitch_y = Var(initialize=0.1, doc="Pitch between Two Neighboring Columns")

        # Pitch of tubes between two neighboring rows (in x direction).
        # Always greater than tube outer diameter
        self.pitch_x = Var(initialize=0.1, doc="Pitch between Two Neighboring Rows")

        # Length of tube per segment in z direction
        self.tube_length_seg = Var(initialize=1.0, doc="Length of Tube per Segment")

        # Minimum cross section area on shell side
        self.area_flow_shell_min = Var(
            initialize=1.0, doc="Minimum Flow Area on Shell Side"
        )

        # total number of tube rows
        @self.Expression(doc="Total Number of Tube Rows")
        def tube_nrow(b):
            return b.tube_nseg * b.tube_inlet_nrow

        # Tube outside diameter
        @self.Expression(doc="Outside Diameter of Tube")
        def tube_do(b):
            return b.tube_di + b.tube_thickness * 2.0

        # Mean beam length for radiation
        if self.config.has_radiation is True:

            @self.Expression(doc="Mean Beam Length")
            def mbl(b):
                return 3.6 * (
                    b.pitch_x * b.pitch_y / const.pi / b.tube_do - b.tube_do / 4.0
                )

            # Mean beam length for radiation divided by sqrt(2)
            @self.Expression(doc="Sqrt(1/2) of Mean Beam Length")
            def mbl_div2(b):
                return b.mbl / sqrt(2.0)

            # Mean beam length for radiation multiplied by sqrt(2)
            @self.Expression(doc="Sqrt(2) of Mean Beam Length")
            def mbl_mul2(b):
                return b.mbl * sqrt(2.0)

        # Ratio of pitch_x/tube_do
        @self.Expression(doc="Ratio of Pitch in x Direction" "to Tube Outside Diameter")
        def pitch_x_to_do(b):
            return b.pitch_x / b.tube_do

        # Ratio of pitch_y/tube_do
        @self.Expression(
            doc="Ratio of Pitch in y Direction " "to Tube Outside Diameter"
        )
        def pitch_y_to_do(b):
            return b.pitch_y / b.tube_do

        # Total cross section area of tube metal per segment
        @self.Expression(doc="Total Cross Section Area of" "Tube Metal Per Segment")
        def area_wall_seg(b):
            return (
                0.25
                * const.pi
                * (b.tube_do**2 - b.tube_di**2)
                * b.tube_ncol
                * b.tube_inlet_nrow
            )

        # Length of shell side flow
        @self.Constraint(doc="Length of Shell Side Flow")
        def length_flow_shell_eqn(b):
            return b.length_flow_shell == b.tube_nrow * b.pitch_x

        # Length of tube side flow
        @self.Constraint(doc="Length of Tube Side Flow")
        def length_flow_tube_eqn(b):
            return b.length_flow_tube == b.tube_nseg * b.tube_length_seg

        # Total flow area on tube side
        @self.Constraint(doc="Total Area of Tube Flow")
        def area_flow_tube_eqn(b):
            return (
                b.area_flow_tube
                == 0.25 * const.pi * b.tube_di**2.0 * b.tube_ncol * b.tube_inlet_nrow
            )

        # Average flow area on shell side
        @self.Constraint(doc="Average Cross Section Area of Shell Side Flow")
        def area_flow_shell_eqn(b):
            return (
                b.length_flow_shell * b.area_flow_shell
                == b.tube_length_seg * b.length_flow_shell * b.pitch_y * b.tube_ncol
                - b.tube_ncol
                * b.tube_nrow
                * 0.25
                * const.pi
                * b.tube_do**2
                * b.tube_length_seg
            )

        # Minimum flow area on shell side
        @self.Constraint(doc="Minimum Flow Area on Shell Side")
        def area_flow_shell_min_eqn(b):
            return (
                b.area_flow_shell_min
                == b.tube_length_seg * (b.pitch_y - b.tube_do) * b.tube_ncol
            )

        # Note that the volumes of both sides are
        # calculated by the ControlVolume1D
        @self.Expression(doc="Inner Radius of Tube")
        def tube_ri(b):
            return b.tube_di / 2.0

        @self.Expression(doc="Outside Radius of Tube")
        def tube_ro(b):
            return b.tube_ri + b.tube_thickness

        if self.config.has_header is True:
            # Header outside diameter
            @self.Expression(doc="Outside Diameter of Header")
            def head_do(b):
                return b.head_di + b.head_thickness * 2.0

            @self.Expression(doc="Inner Radius of Header")
            def head_ri(b):
                return b.head_di / 2.0

            @self.Expression(doc="Outside Radius of Header")
            def head_ro(b):
                return b.head_ri + b.head_thickness

            self.head_r = ContinuousSet(
                bounds=(
                    value(self.head_ri / self.head_ri_scaling),
                    value(self.head_ro / self.head_ri_scaling),
                )
            )

        # Define the continuous domains for model
        self.r = ContinuousSet(
            bounds=(
                value(self.tube_ri / self.ri_scaling),
                value(self.tube_ro / self.ri_scaling),
            )
        )

    def _make_performance(self):
        """
        Constraints for Unit Model.

        Args:
            None

        Returns:
            None
        """
        # add object Reference
        self.tube_heat = Reference(self.tube.heat)
        self.shell_heat = Reference(self.shell.heat)
        phase_s = self.config.tube_side_water_phase

        # Parameters
        if self.config.has_radiation is True:
            # tube wall emissivity, converted from parameter to variable
            self.emissivity_wall = Var(initialize=0.7, doc="Shell Side Wall Emissivity")

        # Wall thermal conductivity
        self.therm_cond_wall = Param(
            initialize=43.0,
            mutable=True,
            doc="Thermal Conductivity of" "Tube Wall Material",
        )

        # Wall heat capacity
        self.cp_wall = Param(
            initialize=502.4, mutable=True, doc="Tube Wall Heat Capacity"
        )

        # Wall density
        self.dens_wall = Param(initialize=7800.0, mutable=True, doc="Tube Wall Density")

        # Young modulus
        self.Young_modulus = Param(
            initialize=1.90e5, mutable=True, doc="Tube Wall Young Modulus"
        )

        # Poisson's ratio
        self.Poisson_ratio = Param(
            initialize=0.29, mutable=True, doc="Tube Wall Poisson Ratio"
        )

        # Coefficient of thermal expansion
        self.coefficient_thermal_expansion = Param(
            initialize=1.2e-5,
            mutable=True,
            doc="Tube Wall Coefficient" "of Thermal Expansion",
        )

        # thermal diffusivity of wall
        @self.Expression(doc="Thermal Diffusivity of Tube Wall Material")
        def diff_therm_wall(b):
            return b.therm_cond_wall / (b.dens_wall * b.cp_wall)

        # Loss coefficient for a 180 degree bend (u-turn),
        # usually related to radius to inside diameter ratio
        self.kloss_uturn = Param(
            initialize=0.5, mutable=True, doc="Loss Coefficient of a Tube U-Turn"
        )

        # Heat transfer resistance due to the fouling on tube side
        self.tube_r_fouling = Param(
            initialize=0.0, mutable=True, doc="Fouling resistance on tube side"
        )

        # Heat transfer resistance due to the fouling on shell side
        self.shell_r_fouling = Param(
            initialize=0.0001, mutable=True, doc="Fouling Resistance on Tube Side"
        )

        # Correction factor for convective heat transfer
        # coefficient on shell side
        self.fcorrection_htc_shell = Var(
            initialize=1.0, doc="Correction Factor for" "Convective HTC on Shell"
        )

        # Correction factor for convective heat transfer
        # coefficient on tube side
        self.fcorrection_htc_tube = Var(
            initialize=1.0, doc="Correction Factor for Convective" "HTC on Tube Side"
        )

        # Correction factor for tube side pressure drop due to friction
        self.fcorrection_dp_tube = Var(
            initialize=1.0, doc="Correction Factor for Tube Side" "Pressure Drop"
        )

        # Correction factor for shell side pressure drop due to friction
        self.fcorrection_dp_shell = Var(
            initialize=1.0, doc="Correction Factor for Shell Side" "Pressure Drop"
        )

        # Performance variables
        if self.config.has_radiation is True:
            # Gas emissivity at mbl
            self.gas_emissivity = Var(
                self.flowsheet().time,
                self.shell.length_domain,
                initialize=0.5,
                doc="Emissivity at Given" "Mean Beam Length",
            )

            # Gas emissivity at mbl/sqrt(2)
            self.gas_emissivity_div2 = Var(
                self.flowsheet().time,
                self.shell.length_domain,
                initialize=0.4,
                doc="Emissivity at Mean Beam Length" "Divided by Sqrt of 2",
            )

            # Gas emissivity at mbl*sqrt(2)
            self.gas_emissivity_mul2 = Var(
                self.flowsheet().time,
                self.shell.length_domain,
                initialize=0.6,
                doc="Emissivity at Mean Beam" "Length Multiplied by Sqrt Of 2",
            )

            # Gray fraction of gas in entire spectrum
            self.gas_gray_fraction = Var(
                self.flowsheet().time,
                self.shell.length_domain,
                initialize=0.5,
                doc="Gray Fraction of Gas" "in Entire Spectrum",
            )

            # Gas-surface radiation exchange factor for shell side wall
            self.frad_gas_shell = Var(
                self.flowsheet().time,
                self.shell.length_domain,
                initialize=0.5,
                doc="Gas-Surface Radiation Exchange" "Factor for Shell Side Wall",
            )

            # Shell side equivalent convective heat transfer coefficient
            # due to radiation
            self.hconv_shell_rad = Var(
                self.flowsheet().time,
                self.shell.length_domain,
                initialize=100.0,
                doc="Shell Side Convective Heat"
                "Transfer Coefficient due to Radiation",
            )

        # Tube side convective heat transfer coefficient
        self.hconv_tube = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=100.0,
            doc="Tube Side Convective" "Heat Transfer Coefficient",
        )

        # Tube side convective heat transfer coefficient combined with fouling
        self.hconv_tube_foul = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=100.0,
            doc="Tube Side Convective Heat Transfer"
            "Coefficient Combined with Fouling",
        )

        # Shell side convective heat transfer coefficient
        # due to convection only
        self.hconv_shell_conv = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=100.0,
            doc="Shell Side Convective Heat Transfer" "Coefficient due to Convection",
        )

        # Total shell side convective heat transfer coefficient
        # including convection and radiation
        self.hconv_shell_total = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=150.0,
            doc="Total Shell Side Convective" "Heat Transfer Coefficient",
        )

        # Total shell side convective heat transfer coefficient
        # combined with fouling
        self.hconv_shell_foul = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=150.0,
            doc="Shell Side Convective Heat Transfer"
            "Coefficient Combined with Fouling",
        )

        # Constraint for hconv_tube_foul
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Tube Side Convective Heat Transfer" "Coefficient with Fouling",
        )
        def hconv_tube_foul_eqn(b, t, x):
            return (
                0.01
                * b.hconv_tube_foul[t, x]
                * (1 + b.hconv_tube[t, x] * b.tube_r_fouling)
                == 0.01 * b.hconv_tube[t, x]
            )

        # Constraint for hconv_shell_foul
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Shell Side Convective Heat" "Transfer Coefficient with Fouling",
        )
        def hconv_shell_foul_eqn(b, t, x):
            return (
                0.1
                * b.hconv_shell_foul[t, x]
                * (1 + b.hconv_shell_total[t, x] * b.shell_r_fouling)
                == 0.1 * b.hconv_shell_total[t, x]
            )

        # Tube metal wall temperature profile across radius
        self.tube_wall_temperature = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            initialize=500,
            doc="Tube Wall Temperature",
        )

        self.shell_wall_temperature = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=500,
            doc="Shell Side Fouling Wall" "Surface Temperature",
        )

        # Fouling wall surface temperature on shell side
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Fouling Wall Surface" "Temperature on Shell Side",
        )
        def temp_wall_shell_eqn(b, t, x):
            return b.shell_wall_temperature[t, x] == b.tube_wall_temperature[
                t, x, b.r.last()
            ] + b.shell_r_fouling * b.hconv_shell_foul[t, x] * (
                b.shell.properties[t, x].temperature
                - b.tube_wall_temperature[t, x, b.r.last()]
            )

        # Declare derivatives in the model
        if self.config.dynamic is True:
            self.dTdt = DerivativeVar(
                self.tube_wall_temperature, wrt=self.flowsheet().time
            )
        self.dTdr = DerivativeVar(self.tube_wall_temperature, wrt=self.r)
        self.d2Tdr2 = DerivativeVar(self.tube_wall_temperature, wrt=(self.r, self.r))

        discretizer = TransformationFactory("dae.finite_difference")
        discretizer.apply_to(
            self, nfe=self.config.radial_elements, wrt=self.r, scheme="CENTRAL"
        )

        # Constraint for heat conduction equation
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="1-D Heat Conduction Equation Through Radius",
        )
        def heat_conduction_eqn(b, t, x, r):
            if r == b.r.first() or r == b.r.last():
                return Constraint.Skip
            if self.config.dynamic is True:
                return b.dTdt[t, x, r] == b.diff_therm_wall / b.ri_scaling**2 * (
                    b.d2Tdr2[t, x, r] + b.dTdr[t, x, r] / r
                )
            else:
                return 0 == b.diff_therm_wall / b.ri_scaling**2 * (
                    b.d2Tdr2[t, x, r] + b.dTdr[t, x, r] / r
                )

        @self.Constraint(
            self.flowsheet().time, self.tube.length_domain, doc="Inner Wall Boundary"
        )
        def inner_wall_bc_eqn(b, t, x):
            return (
                0.01
                * b.hconv_tube_foul[t, x]
                * (
                    b.tube.properties[t, x].temperature
                    - b.tube_wall_temperature[t, x, b.r.first()]
                )
                == -0.01 * b.dTdr[t, x, b.r.first()] / b.ri_scaling * b.therm_cond_wall
            )

        @self.Constraint(
            self.flowsheet().time, self.shell.length_domain, doc="Outer Wall Boundary"
        )
        def outer_wall_bc_eqn(b, t, x):
            return (
                0.01
                * b.hconv_shell_foul[t, x]
                * (
                    b.tube_wall_temperature[t, x, b.r.last()]
                    - b.shell.properties[t, x].temperature
                )
                == -0.01 * b.dTdr[t, x, b.r.last()] / b.ri_scaling * b.therm_cond_wall
            )

        # Inner wall BC for dTdt
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Extra Inner Wall Temperature Derivative",
        )
        def extra_at_inner_wall_eqn(b, t, x):
            if self.config.dynamic is True:
                term = b.dTdt[t, x, b.r.first()]
            else:
                term = 0
            return term == 4 * b.diff_therm_wall * (b.r.first() + b.r.at(2)) / (
                b.r.at(2) - b.r.first()
            ) ** 2 / (3 * b.r.first() + b.r.at(2)) / b.ri_scaling**2 * (
                b.tube_wall_temperature[t, x, b.r.at(2)]
                - b.tube_wall_temperature[t, x, b.r.first()]
            ) + 8 * b.diff_therm_wall / b.therm_cond_wall * b.hconv_tube_foul[
                t, x
            ] * b.r.first() / (
                b.r.at(2) - b.r.first()
            ) / (
                3 * b.r.first() + b.r.at(2)
            ) / b.ri_scaling * (
                b.tube.properties[t, x].temperature
                - b.tube_wall_temperature[t, x, b.r.first()]
            )

        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Extra Outer Wall Temperature Derivative",
        )
        def extra_at_outer_wall_eqn(b, t, x):
            if self.config.dynamic is True:
                term = b.dTdt[t, x, b.r.last()]
            else:
                term = 0
            return term == 4 * b.diff_therm_wall * (b.r.last() + b.r.at(-2)) / (
                b.r.last() - b.r.at(-2)
            ) ** 2 / (3 * b.r.last() + b.r.at(-2)) / b.ri_scaling**2 * (
                b.tube_wall_temperature[t, x, b.r.at(-2)]
                - b.tube_wall_temperature[t, x, b.r.last()]
            ) + 8 * b.diff_therm_wall / b.therm_cond_wall * b.hconv_shell_foul[
                t, x
            ] * b.r.last() / (
                b.r.last() - b.r.at(-2)
            ) / (
                3 * b.r.last() + b.r.at(-2)
            ) / b.ri_scaling * (
                b.shell.properties[t, x].temperature
                - b.tube_wall_temperature[t, x, b.r.last()]
            )

        if self.config.has_radiation is True:
            # Constraints for gas emissivity
            @self.Constraint(
                self.flowsheet().time, self.shell.length_domain, doc="Gas Emissivity"
            )
            def gas_emissivity_eqn(b, t, x):
                # This is a surrogate model, so need to do units manually
                X1 = b.shell.properties[t, x].temperature
                X2 = b.mbl
                X3 = b.shell.properties[t, x].pressure
                X4 = b.shell.properties[t, x].mole_frac_comp["CO2"]
                X5 = b.shell.properties[t, x].mole_frac_comp["H2O"]
                X6 = b.shell.properties[t, x].mole_frac_comp["O2"]
                # Surrogate model fitted using rigorous calc. - 500 samples
                # Wide operating range:
                #       X1: 700 – 1500    (Gas Temperature)
                #       X2: 0.2 – 1       (Mean beam length)
                #       X3: 79000-102000  (pressure in Pa)
                #       X4: 0.12-0.16     (mol frac CO2)
                #       X5: 0.075-0.15    (mol frac H2O)
                #       X6: 0.01-0.07     (mol frac O2)
                return (
                    b.gas_emissivity[t, x]
                    == -0.000116906 * X1
                    + 1.02113 * X2
                    + 4.81687e-07 * X3
                    + 0.922679 * X4
                    - 0.0708822 * X5
                    - 0.0368321 * X6
                    + 0.121843 * log(X1)
                    + 0.0353343 * log(X2)
                    + 0.0346181 * log(X3)
                    + 0.0180859 * log(X5)
                    - 0.256274 * exp(X2)
                    - 0.674791 * exp(X4)
                    - 0.724802 * sin(X2)
                    - 0.0206726 * cos(X2)
                    - 9.01012e-05 * cos(X3)
                    - 3.09283e-05 * X1 * X2
                    - 5.44339e-10 * X1 * X3
                    - 0.000196134 * X1 * X5
                    + 4.54838e-05 * X1 * X6
                    + 7.57411e-07 * X2 * X3
                    + 0.0395456 * X2 * X4
                    + 0.726625 * X2 * X5
                    - 0.034842 * X2 * X6
                    + 4.00056e-06 * X3 * X5
                    + 5.71519e-09 * (X1 * X2) ** 2
                    - 1.27853 * (X2 * X5) ** 2
                )

            # Constraints for gas emissivity at mbl/sqrt(2)
            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="Gas Emissivity at a Lower Mean Beam Length",
            )
            def gas_emissivity_div2_eqn(b, t, x):
                X1 = b.shell.properties[t, x].temperature
                X2 = b.mbl_div2
                X3 = b.shell.properties[t, x].pressure
                X4 = b.shell.properties[t, x].mole_frac_comp["CO2"]
                X5 = b.shell.properties[t, x].mole_frac_comp["H2O"]
                X6 = b.shell.properties[t, x].mole_frac_comp["O2"]
                # Surrogate model fitted using rigorous calc. - 500 samples
                # Wide operating range:
                #       X1: 700 – 1500    (Gas Temperature)
                #       X2: 0.2 – 1       (Mean beam length)
                #       X3: 79000-102000  (pressure in Pa)
                #       X4: 0.12-0.16     (mol frac CO2)
                #       X5: 0.075-0.15    (mol frac H2O)
                #       X6: 0.01-0.07     (mol frac O2)
                return (
                    b.gas_emissivity_div2[t, x]
                    == -0.000116906 * X1
                    + 1.02113 * X2
                    + 4.81687e-07 * X3
                    + 0.922679 * X4
                    - 0.0708822 * X5
                    - 0.0368321 * X6
                    + 0.121843 * log(X1)
                    + 0.0353343 * log(X2)
                    + 0.0346181 * log(X3)
                    + 0.0180859 * log(X5)
                    - 0.256274 * exp(X2)
                    - 0.674791 * exp(X4)
                    - 0.724802 * sin(X2)
                    - 0.0206726 * cos(X2)
                    - 9.01012e-05 * cos(X3)
                    - 3.09283e-05 * X1 * X2
                    - 5.44339e-10 * X1 * X3
                    - 0.000196134 * X1 * X5
                    + 4.54838e-05 * X1 * X6
                    + 7.57411e-07 * X2 * X3
                    + 0.0395456 * X2 * X4
                    + 0.726625 * X2 * X5
                    - 0.034842 * X2 * X6
                    + 4.00056e-06 * X3 * X5
                    + 5.71519e-09 * (X1 * X2) ** 2
                    - 1.27853 * (X2 * X5) ** 2
                )

            # Constraints for gas emissivity at mbl*sqrt(2)
            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="Gas Emissivity at a Higher Mean Beam Length",
            )
            def gas_emissivity_mul2_eqn(b, t, x):
                X1 = b.shell.properties[t, x].temperature
                X2 = b.mbl_mul2
                X3 = b.shell.properties[t, x].pressure
                X4 = b.shell.properties[t, x].mole_frac_comp["CO2"]
                X5 = b.shell.properties[t, x].mole_frac_comp["H2O"]
                X6 = b.shell.properties[t, x].mole_frac_comp["O2"]
                # Surrogate model fitted using rigorous calc. - 500 samples
                # Wide operating range:
                #       X1: 700 – 1500    (Gas Temperature)
                #       X2: 0.2 – 1       (Mean beam length)
                #       X3: 79000-102000  (pressure in Pa)
                #       X4: 0.12-0.16     (mol frac CO2)
                #       X5: 0.075-0.15    (mol frac H2O)
                #       X6: 0.01-0.07     (mol frac O2)
                return (
                    b.gas_emissivity_mul2[t, x]
                    == -0.000116906 * X1
                    + 1.02113 * X2
                    + 4.81687e-07 * X3
                    + 0.922679 * X4
                    - 0.0708822 * X5
                    - 0.0368321 * X6
                    + 0.121843 * log(X1)
                    + 0.0353343 * log(X2)
                    + 0.0346181 * log(X3)
                    + 0.0180859 * log(X5)
                    - 0.256274 * exp(X2)
                    - 0.674791 * exp(X4)
                    - 0.724802 * sin(X2)
                    - 0.0206726 * cos(X2)
                    - 9.01012e-05 * cos(X3)
                    - 3.09283e-05 * X1 * X2
                    - 5.44339e-10 * X1 * X3
                    - 0.000196134 * X1 * X5
                    + 4.54838e-05 * X1 * X6
                    + 7.57411e-07 * X2 * X3
                    + 0.0395456 * X2 * X4
                    + 0.726625 * X2 * X5
                    - 0.034842 * X2 * X6
                    + 4.00056e-06 * X3 * X5
                    + 5.71519e-09 * (X1 * X2) ** 2
                    - 1.27853 * (X2 * X5) ** 2
                )

            # fraction of gray gas spectrum
            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="Fraction of Gray Gas Spectrum",
            )
            def gas_gray_fraction_eqn(b, t, x):
                return (
                    b.gas_gray_fraction[t, x]
                    * (2 * b.gas_emissivity_div2[t, x] - b.gas_emissivity_mul2[t, x])
                    == b.gas_emissivity_div2[t, x] ** 2
                )

            # Gas-surface radiation exchange factor between
            # gas and shell side wall
            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="Gas-Surface Radiation Exchange"
                "Factor between Gas and Shell Side Wall",
            )
            def frad_gas_shell_eqn(b, t, x):
                return (
                    b.frad_gas_shell[t, x]
                    * (
                        (1 / b.emissivity_wall - 1) * b.gas_emissivity[t, x]
                        + b.gas_gray_fraction[t, x]
                    )
                    == b.gas_gray_fraction[t, x] * b.gas_emissivity[t, x]
                )

            # equivalent convective heat transfer coefficent due to radiation
            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="Equivalent Convective Heat Transfer"
                "Coefficient due to Radiation",
            )
            def hconv_shell_rad_eqn(b, t, x):
                return b.hconv_shell_rad[
                    t, x
                ] == const.stefan_constant * b.frad_gas_shell[t, x] * (
                    b.shell.properties[t, x].temperature
                    + b.shell_wall_temperature[t, x]
                ) * (
                    b.shell.properties[t, x].temperature ** 2
                    + b.shell_wall_temperature[t, x] ** 2
                )

        # Tube side heat transfer coefficient and pressure drop
        # -----------------------------------------------------
        # Velocity on tube side
        self.tube_velocity = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=1.0,
            doc="Velocity on Tube Side",
        )

        # Reynalds number on tube side
        self.tube_N_Re = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=10000.0,
            doc="Reynolds Number on Tube Side",
        )

        # Friction factor on tube side
        self.tube_friction_factor = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=1.0,
            doc="Friction Factor on Tube Side",
        )

        # Pressure drop due to friction on tube side
        self.deltaP_tube_friction = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=-10.0,
            doc="Pressure Drop due to" "Friction on Tube Side",
        )

        # Pressure drop due to 180 degree turn on tube side
        self.deltaP_tube_uturn = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=-10.0,
            doc="Pressure Drop due to U-Turn on" "Tube Side",
        )

        # Prandtl number on tube side
        self.tube_N_Pr = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=1.0,
            doc="Prandtl Number on Tube Side",
        )

        # Nusselt number on tube side
        self.tube_N_Nu = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=1,
            doc="Nusselts Number on Tube Side",
        )

        # Velocity equation
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Tube Side Velocity Equation",
        )
        def v_tube_eqn(b, t, x):
            return (
                0.001
                * b.tube_velocity[t, x]
                * b.area_flow_tube
                * b.tube.properties[t, x].dens_mol_phase[phase_s]
                == 0.001 * b.tube.properties[t, x].flow_mol
            )

        # Reynolds number
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Reynolds Number Equation on Tube Side",
        )
        def N_Re_tube_eqn(b, t, x):
            return (
                b.tube_N_Re[t, x] * b.tube.properties[t, x].visc_d_phase[phase_s]
                == b.tube_di
                * b.tube_velocity[t, x]
                * b.tube.properties[t, x].dens_mass_phase[phase_s]
            )

        # Friction factor
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Darcy Friction Factor on Tube Side",
        )
        def friction_factor_tube_eqn(b, t, x):
            return (
                b.tube_friction_factor[t, x]
                * smooth_max(b.tube_N_Re[t, x], 1, 1e-5) ** 0.25
                == 0.3164 * b.fcorrection_dp_tube
            )

        # Pressure drop due to friction per tube length
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Pressure Drop due to Friction per Tube Length",
        )
        def deltaP_tube_friction_eqn(b, t, x):
            return (
                b.deltaP_tube_friction[t, x] * b.tube_di
                == -0.5
                * b.tube.properties[t, x].dens_mass_phase[phase_s]
                * b.tube_velocity[t, x] ** 2
                * b.tube_friction_factor[t, x]
            )

        # Pressure drop due to u-turn
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Pressure Drop due to U-Turn on Tube Side",
        )
        def deltaP_tube_uturn_eqn(b, t, x):
            return (
                b.deltaP_tube_uturn[t, x] * b.tube_length_seg
                == -0.5
                * b.tube.properties[t, x].dens_mass_phase[phase_s]
                * b.tube_velocity[t, x] ** 2
                * b.kloss_uturn
            )

        # Total pressure drop on tube side
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Total Pressure Drop on Tube Side",
        )
        def deltaP_tube_eqn(b, t, x):
            return b.tube.deltaP[t, x] == (
                b.deltaP_tube_friction[t, x]
                + b.deltaP_tube_uturn[t, x]
                - b.delta_elevation
                / b.tube_nseg
                * const.acceleration_gravity
                * b.tube.properties[t, x].dens_mass_phase[phase_s]
                / b.tube_length_seg
            )

        # Prandtl number
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Prandtl Number Equation on Tube Side",
        )
        def N_Pr_tube_eqn(b, t, x):
            return (
                b.tube_N_Pr[t, x]
                * b.tube.properties[t, x].therm_cond_phase[phase_s]
                * b.tube.properties[t, x].mw
                == b.tube.properties[t, x].cp_mol_phase[phase_s]
                * b.tube.properties[t, x].visc_d_phase[phase_s]
            )

        # Nusselts number
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Nusselts Number Equation on Tube Side",
        )
        def N_Nu_tube_eqn(b, t, x):
            return (
                b.tube_N_Nu[t, x]
                == 0.023
                * smooth_max(b.tube_N_Re[t, x], 1, 1e-5) ** 0.8
                * b.tube_N_Pr[t, x] ** 0.4
            )

        # Heat transfer coefficient
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Convective Heat Transfer Coefficient" "Equation on Tube Side",
        )
        def hconv_tube_eqn(b, t, x):
            return (
                b.hconv_tube[t, x] * self.tube_di
                == b.tube_N_Nu[t, x]
                * b.tube.properties[t, x].therm_cond_phase[phase_s]
                * b.fcorrection_htc_tube
            )

        # Pressure drop and heat transfer coefficient on shell side
        # ----------------------------------------------------------
        # Tube arrangement factor
        if self.config.tube_arrangement == "in-line":
            self.f_arrangement = Param(
                initialize=0.788, doc="In-Line Tube Arrangement Factor"
            )
        elif self.config.tube_arrangement == "staggered":
            self.f_arrangement = Param(
                initialize=1.0, doc="Staggered Tube Arrangement Factor"
            )
        else:
            raise Exception("Tube Arrangement Not Supported")

        # Velocity on shell side
        self.shell_velocity = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=1.0,
            doc="Velocity on Shell Side",
        )

        # Reynalds number on shell side
        self.shell_N_Re = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=10000.0,
            doc="Reynolds Number on Shell Side",
        )

        # Friction factor on shell side
        self.shell_friction_factor = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=1.0,
            doc="Friction Factor on Shell Side",
        )

        # Prandtl number on shell side
        self.shell_N_Pr = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=1,
            doc="Prandtl Number on Shell Side",
        )

        # Nusselt number on shell side
        self.shell_N_Nu = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=1,
            doc="Nusselts Number on Shell Side",
        )

        # Velocity equation on shell side, using inlet molar flow rate
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Velocity on Shell Side",
        )
        def v_shell_eqn(b, t, x):
            return (
                b.shell_velocity[t, x]
                * b.shell.properties[t, x].dens_mol_phase["Vap"]
                * b.area_flow_shell_min
                == b.shell.properties[t, 0].flow_mol
            )

        # Reynolds number
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Reynolds Number Equation on Shell Side",
        )
        def N_Re_shell_eqn(b, t, x):
            return (
                b.shell_N_Re[t, x] * b.shell.properties[t, x].visc_d
                == b.tube_do
                * b.shell_velocity[t, x]
                * b.shell.properties[t, x].dens_mol_phase["Vap"]
                * b.shell.properties[t, x].mw
            )

        # Friction factor on shell side
        if self.config.tube_arrangement == "in-line":

            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="In-Line Friction Factor on Shell Side",
            )
            def friction_factor_shell_eqn(b, t, x):
                return (
                    b.shell_friction_factor[t, x]
                    * smooth_max(b.shell_N_Re[t, x], 1, 1e-5) ** 0.15
                    == (
                        0.044
                        + 0.08
                        * b.pitch_x_to_do
                        / (b.pitch_y_to_do - 1.0) ** (0.43 + 1.13 / b.pitch_x_to_do)
                    )
                    * b.fcorrection_dp_shell
                )

        elif self.config.tube_arrangement == "staggered":

            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="Staggered Friction Factor on Shell Side",
            )
            def friction_factor_shell_eqn(b, t, x):
                return (
                    b.shell_friction_factor[t, x]
                    * smooth_max(b.shell_N_Re[t, x], 1, 1e-5) ** 0.16
                    == (0.25 + 0.118 / (b.pitch_y_to_do - 1.0) ** 1.08)
                    * b.fcorrection_dp_shell
                )

        else:
            raise Exception("Tube Arrangement Not Supported")

        # Pressure drop on shell side
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Pressure Change on Shell Side",
        )
        def deltaP_shell_eqn(b, t, x):
            return (
                b.shell.deltaP[t, x] * b.pitch_x
                == -1.4
                * b.shell_friction_factor[t, x]
                * b.shell.properties[t, x].dens_mol_phase["Vap"]
                * b.shell.properties[t, x].mw
                * b.shell_velocity[t, x] ** 2
            )

        # Prandtl number
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Prandtl Number Equation on Shell Side",
        )
        def N_Pr_shell_eqn(b, t, x):
            return (
                b.shell_N_Pr[t, x]
                * b.shell.properties[t, x].therm_cond
                * b.shell.properties[t, x].mw
                == b.shell.properties[t, x].cp_mol * b.shell.properties[t, x].visc_d
            )

        # Nusselt number, currently assume Re > 300
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Nusselts Number Equation on Shell Side",
        )
        def N_Nu_shell_eqn(b, t, x):
            return (
                b.shell_N_Nu[t, x]
                == b.f_arrangement
                * 0.33
                * smooth_max(b.shell_N_Re[t, x], 1, 1e-5) ** 0.6
                * b.shell_N_Pr[t, x] ** 0.333333
            )

        # Convective heat transfer coefficient on shell side
        # due to convection only
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Convective Heat Transfer Coefficient Equation"
            "on Shell Side due to Convection",
        )
        def hconv_shell_conv_eqn(b, t, x):
            return (
                b.hconv_shell_conv[t, x] * b.tube_do
                == b.shell_N_Nu[t, x]
                * b.shell.properties[t, x].therm_cond
                * b.fcorrection_htc_shell
            )

        # Total convective heat transfer coefficient on shell side
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Total Convective Heat Transfer Coefficient" "Equation on Shell Side",
        )
        def hconv_shell_total_eqn(b, t, x):
            if self.config.has_radiation is True:
                return (
                    b.hconv_shell_total[t, x]
                    == b.hconv_shell_conv[t, x] + b.hconv_shell_rad[t, x]
                )
            else:
                return b.hconv_shell_total[t, x] == b.hconv_shell_conv[t, x]

        # Energy balance with tube wall
        # ------------------------------------
        # Heat to wall per length on tube side
        @self.Constraint(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Heat per Length on Tube Side",
        )
        def heat_tube_eqn(b, t, x):
            return b.tube_heat[t, x] == b.hconv_tube_foul[
                t, x
            ] * const.pi * b.tube_di * b.tube_inlet_nrow * b.tube_ncol * (
                b.tube_wall_temperature[t, x, b.r.first()]
                - b.tube.properties[t, x].temperature
            )

        # Heat to wall per length on shell side
        @self.Constraint(
            self.flowsheet().time,
            self.shell.length_domain,
            doc="Heat per Length on Shell Side",
        )
        def heat_shell_eqn(b, t, x):
            return b.shell_heat[
                t, x
            ] * b.length_flow_shell == b.length_flow_tube * b.hconv_shell_foul[
                t, x
            ] * const.pi * b.tube_do * b.tube_inlet_nrow * b.tube_ncol * (
                b.tube_wall_temperature[t, x, b.r.last()]
                - b.shell.properties[t, x].temperature
            )

        # Calculate mechanical and thermal stresses based on
        # EN 13445 for SH thick walled component
        # ---------------------------------------------------------------------
        # Integer indexing for radius domain
        self.rindex = Param(
            self.r, initialize=1, mutable=True, doc="Integer Indexing for Radius Domain"
        )

        # Calculate integral point for mean temperature in the wall
        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Mean Temperature across the Wall",
        )
        def mean_temperature(b, t, x):
            return (
                2
                * (b.r.at(2) - b.r.at(1))
                * b.ri_scaling**2
                / (b.tube_ro**2 - b.tube_ri**2)
                * (
                    sum(
                        0.5
                        * (
                            b.r.at(i - 1) * b.tube_wall_temperature[t, x, b.r.at(i - 1)]
                            + b.r.at(i) * b.tube_wall_temperature[t, x, b.r.at(i)]
                        )
                        for i in range(2, len(b.r) + 1)
                    )
                )
            )

        for index_r, value_r in enumerate(self.r, 1):
            self.rindex[value_r] = index_r

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Discrete Point Mean Temperature",
        )
        def discrete_mean_temperature(b, t, x, r):
            if b.rindex[r].value == 1:
                return b.tube_wall_temperature[t, x, b.r.first()]
            else:
                return (
                    2
                    * (b.r.at(2) - b.r.at(1))
                    * b.ri_scaling**2
                    / ((b.r.at(b.rindex[r].value) * b.ri_scaling) ** 2 - b.tube_ri**2)
                    * (
                        sum(
                            0.5
                            * (
                                b.r.at(j - 1)
                                * b.tube_wall_temperature[t, x, b.r.at(j - 1)]
                                + b.r.at(j) * b.tube_wall_temperature[t, x, b.r.at(j)]
                            )
                            for j in range(2, b.rindex[r].value + 1)
                        )
                    )
                )

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Thermal Stress at Radial Direction for Tube",
        )
        def therm_sigma_r(b, t, x, r):
            if r == b.r.first() or r == b.r.last():
                return 0
            else:
                return (
                    0.5
                    * b.Young_modulus
                    * b.coefficient_thermal_expansion
                    / (1 - b.Poisson_ratio)
                    * (
                        (1 - b.tube_ri**2 / (r * b.ri_scaling) ** 2)
                        * (
                            b.mean_temperature[t, x]
                            - b.discrete_mean_temperature[t, x, r]
                        )
                    )
                )

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Thermal Stress at " "Circumferential Direction for Tube",
        )
        def therm_sigma_theta(b, t, x, r):
            r_2 = (r * b.ri_scaling) ** 2
            return (
                0.5
                * b.Young_modulus
                * b.coefficient_thermal_expansion
                / (1 - b.Poisson_ratio)
                * (
                    (1 + b.tube_ri**2 / r_2) * b.mean_temperature[t, x]
                    + (1 - b.tube_ri**2 / r_2) * b.discrete_mean_temperature[t, x, r]
                    - 2 * b.tube_wall_temperature[t, x, r]
                )
            )

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Thermal Stress at Axial Direction for Tube",
        )
        def therm_sigma_z(b, t, x, r):
            return (
                b.Young_modulus
                * b.coefficient_thermal_expansion
                / (1 - b.Poisson_ratio)
                * (b.mean_temperature[t, x] - b.tube_wall_temperature[t, x, r])
            )

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Mechanical Stress at Radial Direction for Tube",
        )
        def mech_sigma_r(b, t, x, r):
            if r == b.r.first():
                return 1e-6 * (-b.tube.properties[t, x].pressure)
            elif r == b.r.last():
                return 1e-6 * (-b.shell.properties[t, x].pressure)
            else:
                return 0.1 * (
                    1e-5
                    * (
                        b.tube.properties[t, x].pressure * b.tube_ri**2
                        - b.shell.properties[t, x].pressure * b.tube_ro**2
                    )
                    / (b.tube_ro**2 - b.tube_ri**2)
                    + (
                        1e-5
                        * (
                            b.shell.properties[t, x].pressure
                            - b.tube.properties[t, x].pressure
                        )
                        * b.tube_ri**2
                        * b.tube_ro**2
                        / ((r * b.ri_scaling) ** 2 * (b.tube_ro**2 - b.tube_ri**2))
                    )
                )

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Mechanical Stress at" "Circumferential Direction for Tube",
        )
        def mech_sigma_theta(b, t, x, r):
            return 0.1 * (
                1e-5
                * (
                    b.tube.properties[t, x].pressure * b.tube_ri**2
                    - b.shell.properties[t, x].pressure * b.tube_ro**2
                )
                / (b.tube_ro**2 - b.tube_ri**2)
                - (
                    1e-5
                    * (
                        b.shell.properties[t, x].pressure
                        - b.tube.properties[t, x].pressure
                    )
                    * b.tube_ri**2
                    * b.tube_ro**2
                    / ((r * b.ri_scaling) ** 2 * (b.tube_ro**2 - b.tube_ri**2))
                )
            )

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            doc="Mechanical Stress at Axial Direction for Tube",
        )
        def mech_sigma_z(b, t, x):
            return 0.1 * (
                1e-5
                * (
                    b.tube.properties[t, x].pressure * b.tube_ri**2
                    - b.shell.properties[t, x].pressure * b.tube_ro**2
                )
                / (b.tube_ro**2 - b.tube_ri**2)
            )

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Principal Structural Stress" "at Radial Direction for Tube",
        )
        def sigma_r(b, t, x, r):
            if r == b.r.first():
                return 1e-6 * (-b.tube.properties[t, x].pressure)
            elif r == b.r.last():
                return 1e-6 * (-b.shell.properties[t, x].pressure)
            else:
                return b.mech_sigma_r[t, x, r] + b.therm_sigma_r[t, x, r]

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Principal Structural Stress" "at Circumferential Direction for Tube",
        )
        def sigma_theta(b, t, x, r):
            return b.mech_sigma_theta[t, x, r] + b.therm_sigma_theta[t, x, r]

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Principal Structural Stress" "at Axial Direction for Tube",
        )
        def sigma_z(b, t, x, r):
            return b.mech_sigma_z[t, x] + b.therm_sigma_z[t, x, r]

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Variation Principal Stress"
            "between Radial-Circumferential Directions for Tube",
        )
        def delta_sigma_r_theta(b, t, x, r):
            return abs(b.sigma_r[t, x, r] - b.sigma_theta[t, x, r])

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Variation Principal Stress"
            "between Circumferential-Axial Directions for Tube",
        )
        def delta_sigma_theta_z(b, t, x, r):
            return abs(b.sigma_theta[t, x, r] - b.sigma_z[t, x, r])

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Variation Principal Stress" "between Axial-Radial Directions for Tube",
        )
        def delta_sigma_z_r(b, t, x, r):
            return abs(b.sigma_z[t, x, r] - b.sigma_r[t, x, r])

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Equivalent von Mises Stress for Tube",
        )
        def sigma_von_Mises(b, t, x, r):
            return sqrt(
                b.sigma_r[t, x, r] ** 2
                + b.sigma_theta[t, x, r] ** 2
                + b.sigma_z[t, x, r] ** 2
                - (
                    b.sigma_r[t, x, r] * b.sigma_theta[t, x, r]
                    + b.sigma_r[t, x, r] * b.sigma_z[t, x, r]
                    + b.sigma_theta[t, x, r] * b.sigma_z[t, x, r]
                )
            )

        # Calculate creep // Using property of Steel SA 209 T1
        # data obtained from NIMS databank, date: 11/20/2020
        # https://mits.nims.go.jp/en/
        # data were assessed using Minimum-Commitment parameter
        self.creep_a = Param(initialize=168.15, mutable=True)
        self.creep_b = Param(initialize=-262.34, mutable=True)
        self.creep_c = Param(initialize=123.72, mutable=True)
        self.creep_d = Param(initialize=-19.62, mutable=True)
        self.creep_e = Param(initialize=348582.80, mutable=True)

        @self.Expression(
            self.flowsheet().time,
            self.tube.length_domain,
            self.r,
            doc="Rupture Time for Tube",
        )
        def rupture_time(b, t, x, r):
            return 10 ** (
                b.creep_a
                + b.creep_b * log10(b.sigma_von_Mises[t, x, r])
                + b.creep_c * (log10(b.sigma_von_Mises[t, x, r])) ** 2
                + b.creep_d * (log10(b.sigma_von_Mises[t, x, r])) ** 3
                + b.creep_e / (19.1425 * b.tube_wall_temperature[t, x, r])
            )

        if self.config.has_header is True:
            # -----------------------------------------------------------------
            # calculate temperature distribution and stress for headers
            # @ outlet header (Because the counter current flow type,
            # be careful to define the Tgas, Tst, hgas and hst)
            # material density of header
            # thermal conductivity of header
            # Calculate outlet header of primary superheater
            # Material is SA 387 Grade 12
            # (other name EN: 13CrMo45 and DIN 13CrMo44)

            self.therm_cond_header = Param(initialize=44, mutable=True)

            # density of header
            self.dens_header = Param(initialize=7800, mutable=True)

            # heat capacity of header
            self.cp_header = Param(initialize=470, mutable=True)

            # Young modulus
            self.Young_modulus_header = Param(initialize=2.00e5, mutable=True)

            # Poisson's ratio
            self.Poisson_ratio_header = Param(initialize=0.29, mutable=True)

            # Coefficient of thermal expansion
            self.coefficient_therm_expansion_header = Param(
                initialize=1.2e-5, mutable=True
            )

            # Parameters for heat transfer of outlet header of
            # counter-current flow with insulation
            # Thermal conductivity of air
            self.therm_cond_air = Param(initialize=0.03, mutable=True)

            # Header insulation thermal conductivity
            self.head_therm_cond_insulation = Param(initialize=0.08, mutable=True)

            # Constant related to Rayleigh number for free convection
            # (gravity*expansion_coefficient*density
            # /viscosity/thermal_diffusivity)^(1/6)
            # Use properties at 50 C and 1 atm, 6.84e7^(1/6)=20.223
            self.const_head_Ra_root6 = Param(
                initialize=20.223,
                mutable=True,
                doc="Rayleigh number for free convection",
            )

            # Constant related to Nu for free convection,
            # 0.387/(1+0.721*Pr^(-9/16))^(8/27)
            # Use properties at 50 C and 1 atm
            self.const_head_Nu = Param(
                initialize=0.322,
                mutable=True,
                doc="constant related to Nu number for free convection",
            )

            # Variables for heat transfer of header with insulation
            # Ambient temperature
            self.temperature_ambient = Var(
                self.flowsheet().time, initialize=298.15, doc="Ambient Temperature"
            )

            # Header insulation thickness
            self.head_insulation_thickness = Var(
                initialize=0.05, doc="Header insulation thickness"
            )

            # Header inside heat transfer coefficient
            self.head_hin = Var(
                self.flowsheet().time,
                initialize=1,
                doc="Inside Heat Transfer Coefficient of Header",
            )

            # Header outside heat transfer coefficient
            self.head_hout = Var(
                self.flowsheet().time,
                initialize=1,
                doc="Outside Heat Transfer Coefficient of Header",
            )

            # Insulation free convection heat transfer coefficient
            self.head_h_free_conv = Var(
                self.flowsheet().time,
                initialize=1,
                doc="Insulation Free Convection Heat Transfer Coefficient",
            )

            # Ra number of free convection
            self.head_N_Ra_root6 = Var(
                self.flowsheet().time,
                initialize=80,
                doc="1/6 Power of Ra Number of Free Convection of Air",
            )

            # Nu number of free convection for header shell side
            self.head_N_Nu_shell = Var(
                self.flowsheet().time,
                initialize=1,
                doc="Nu Number of Free Convection of Air",
            )

            # header tube side velocity
            self.head_velocity_tube = Var(
                self.flowsheet().time, initialize=1, doc="header velocity on tube side"
            )

            # Re number on header tube side
            self.head_N_Re_tube = Var(
                self.flowsheet().time, initialize=1, doc="Nu Number on header tube side"
            )

            # Nu number of free convection on header tube side
            self.head_N_Nu_tube = Var(
                self.flowsheet().time, initialize=1, doc="Nu Number on header tube side"
            )

            # Temperature across header wall thickness
            self.header_wall_temperature = Var(
                self.flowsheet().time, self.head_r, bounds=(500, 900), initialize=680
            )

            # Declare derivatives in the model
            if self.config.dynamic is True:
                self.head_dTdt = DerivativeVar(
                    self.header_wall_temperature, wrt=self.flowsheet().config.time
                )
            self.head_dTdr = DerivativeVar(
                self.header_wall_temperature, wrt=self.head_r
            )
            self.head_d2Tdr2 = DerivativeVar(
                self.header_wall_temperature, wrt=(self.head_r, self.head_r)
            )

            discretizer_2 = TransformationFactory("dae.finite_difference")
            discretizer_2.apply_to(
                self,
                nfe=self.config.header_radial_elements,
                wrt=self.head_r,
                scheme="CENTRAL",
            )

            # thermal diffusivity of header
            @self.Expression(doc="Thermal Diffusivity of Header Material")
            def therm_diffus_header(b):
                return b.therm_cond_header / (b.dens_header * b.cp_header)

            # Constraint for heat conduction equation
            @self.Constraint(
                self.flowsheet().time,
                self.head_r,
                doc="1-D PDE Heat Conduction for Header",
            )
            def head_heat_conduction_eqn(b, t, r):
                if r == b.head_r.first() or r == b.head_r.last():
                    return Constraint.Skip
                if self.config.dynamic is True:
                    return b.head_dTdt[
                        t, r
                    ] == b.therm_diffus_header / b.head_ri_scaling**2 * (
                        b.head_d2Tdr2[t, r] + b.head_dTdr[t, r] / r
                    )
                else:
                    return 0 == b.therm_diffus_header / b.head_ri_scaling**2 * (
                        b.head_d2Tdr2[t, r] + b.head_dTdr[t, r] / r
                    )

            @self.Constraint(self.flowsheet().time, doc="Inner Wall Boundary")
            def head_inner_wall_bc_eqn(b, t):
                return (
                    0.01
                    * b.head_hin[t]
                    * (
                        b.tube.properties[t, b.tube.length_domain.first()].temperature
                        - b.header_wall_temperature[t, b.head_r.first()]
                    )
                    == -0.01
                    * b.head_dTdr[t, b.head_r.first()]
                    / b.head_ri_scaling
                    * b.therm_cond_header
                )

            @self.Constraint(self.flowsheet().time, doc="Outer Wall Boundary")
            def head_outer_wall_bc_eqn(b, t):
                return (
                    0.01
                    * b.head_hout[t]
                    * (
                        b.header_wall_temperature[t, b.head_r.last()]
                        - b.temperature_ambient[t]
                    )
                    == -0.01
                    * b.head_dTdr[t, b.head_r.last()]
                    / b.head_ri_scaling
                    * b.therm_cond_header
                )

            # Inner wall BC for dTdt
            @self.Constraint(
                self.flowsheet().time,
                doc="Extra Boundary at Inner Wall" "Temperature Derivative",
            )
            def head_extra_at_inner_wall_eqn(b, t):
                if self.config.dynamic is True:
                    term = b.head_dTdt[t, b.head_r.first()]
                else:
                    term = 0
                return term == (
                    4
                    * b.therm_diffus_header
                    * (b.head_r.first() + b.head_r.at(2))
                    / (b.head_r.at(2) - b.head_r.first()) ** 2
                    / (3 * b.head_r.first() + b.head_r.at(2))
                    / b.head_ri_scaling**2
                    * (
                        b.header_wall_temperature[t, b.head_r.at(2)]
                        - b.header_wall_temperature[t, b.head_r.first()]
                    )
                    + 8
                    * b.therm_diffus_header
                    / b.therm_cond_header
                    * b.head_hin[t]
                    * b.head_r.first()
                    / (b.head_r.at(2) - b.head_r.first())
                    / (3 * b.head_r.first() + b.head_r.at(2))
                    / b.head_ri_scaling
                    * (
                        b.tube.properties[t, b.tube.length_domain.first()].temperature
                        - b.header_wall_temperature[t, b.head_r.first()]
                    )
                )

            @self.Constraint(
                self.flowsheet().time,
                doc="Extra Boundary at Outer Wall " "Temperature Derivative",
            )
            def head_extra_at_outer_wall_eqn(b, t):
                if self.config.dynamic is True:
                    term = b.head_dTdt[t, b.head_r.last()]
                else:
                    term = 0
                return term == (
                    4
                    * b.therm_diffus_header
                    * (b.head_r.last() + b.head_r.at(-2))
                    / (b.head_r.last() - b.head_r.at(-2)) ** 2
                    / (3 * b.head_r.last() + b.head_r.at(-2))
                    / b.head_ri_scaling**2
                    * (
                        b.header_wall_temperature[t, b.head_r.at(-2)]
                        - b.header_wall_temperature[t, b.head_r.last()]
                    )
                    + 8
                    * b.therm_diffus_header
                    / b.therm_cond_header
                    * b.head_hout[t]
                    * b.head_r.last()
                    / (b.head_r.last() - b.head_r.at(-2))
                    / (3 * b.head_r.last() + b.head_r.at(-2))
                    / b.head_ri_scaling
                    * (
                        b.temperature_ambient[t]
                        - b.header_wall_temperature[t, b.head_r.last()]
                    )
                )

            # Calculate header outside and inside heat transfer coefficient
            # Expressure for insulation heat transfer (conduction)
            # resistance based on drum metal outside diameter
            @self.Expression(doc="Heat Transfer Resistance of Insulation")
            def resistance_insulation(b):
                return (
                    b.head_ro
                    * log((b.head_ro + b.head_insulation_thickness) / b.head_ro)
                    / b.head_therm_cond_insulation
                )

            # Equation considering conduction through insulation
            # and free convection between insulation and ambient
            @self.Constraint(
                self.flowsheet().time, doc="Outer Side Heat Transfer Coefficient"
            )
            def head_hout_eqn(b, t):
                return (
                    b.head_hout[t]
                    * (b.resistance_insulation + 1 / b.head_h_free_conv[t])
                    == 1.0
                )

            # Expressure for outside insulation wall temperature
            @self.Expression(
                self.flowsheet().time, doc="Outside Insulation Wall Temperature"
            )
            def head_temp_insulation_outside(b, t):
                return (
                    b.temperature_ambient[t]
                    + (
                        b.header_wall_temperature[t, b.head_r.last()]
                        - b.temperature_ambient[t]
                    )
                    * b.head_hout[t]
                    / b.head_h_free_conv[t]
                )

            # Ra number equation
            @self.Constraint(self.flowsheet().time, doc="Ra Number of Free Convection")
            def head_Ra_number_eqn(b, t):
                return (
                    b.head_N_Ra_root6[t]
                    == b.const_head_Ra_root6
                    * sqrt(b.head_do + 2 * b.head_insulation_thickness)
                    * (b.head_temp_insulation_outside[t] - b.temperature_ambient[t])
                    ** 0.166667
                )

            # Nu number equation
            @self.Constraint(self.flowsheet().time, doc="Nu Number of Free Convection")
            def head_Nu_number_shell_eqn(b, t):
                return (
                    b.head_N_Nu_shell[t]
                    == (0.6 + b.const_head_Nu * b.head_N_Ra_root6[t]) ** 2
                )

            # Free convection coefficient based on the drum outside diameter
            @self.Constraint(
                self.flowsheet().time,
                doc="Free Convection Heat Transfer Coefficient"
                "between Insulation Wall and Ambient",
            )
            def head_h_free_conv_eqn(b, t):
                return (
                    b.head_h_free_conv[t]
                    == b.head_N_Nu_shell[t] * b.therm_cond_air / b.head_do
                )

            # Calculate header inside heat transfer coefficient
            @self.Constraint(self.flowsheet().time, doc="velocity inside header")
            def head_velocity_tube_eqn(b, t):
                return (
                    0.0001
                    * b.head_velocity_tube[t]
                    * const.pi
                    * b.head_ri**2
                    * b.tube.properties[t, b.tube.length_domain.first()].dens_mol_phase[
                        phase_s
                    ]
                    == 0.0001
                    * b.tube.properties[t, b.tube.length_domain.first()].flow_mol
                )

            # Calculate header Re number on tube side
            @self.Constraint(self.flowsheet().time, doc="Reynolds number inside header")
            def head_N_Re_tube_eqn(b, t):
                return (
                    b.head_N_Re_tube[t]
                    * b.tube.properties[t, b.tube.length_domain.first()].visc_d_phase[
                        phase_s
                    ]
                    == b.head_di
                    * b.head_velocity_tube[t]
                    * b.tube.properties[
                        t, b.tube.length_domain.first()
                    ].dens_mass_phase[phase_s]
                )

            # Calculate header Nu number
            @self.Constraint(self.flowsheet().time, doc="Nuseldt number inside header")
            def head_N_Nu_tube_eqn(b, t):
                return (
                    b.head_N_Nu_tube[t]
                    == 0.023
                    * b.head_N_Re_tube[t] ** 0.8
                    * b.tube_N_Pr[t, b.tube.length_domain.first()] ** 0.4
                )

            # Calculate header hin
            @self.Constraint(
                self.flowsheet().time, doc="inside heat transfer coefficient of header"
            )
            def head_hin_eqn(b, t):
                return (
                    b.head_hin[t] * b.head_di
                    == b.head_N_Nu_tube[t]
                    * b.tube.properties[
                        t, b.tube.length_domain.first()
                    ].therm_cond_phase[phase_s]
                )

            # Calculate mechanical and thermal stresses based on EN 13445
            # for thick walled component
            # ------------------------------------------------------------
            # Integer indexing for radius domain
            self.head_rindex = Param(
                self.head_r,
                initialize=1,
                mutable=True,
                doc="Integer Indexing for Radius Domain",
            )

            # calculate integral point for mean temperature in the wall
            @self.Expression(self.flowsheet().time, doc="Mean Temperature for Header")
            def mean_temperature_header(b, t):
                return (
                    2
                    * (b.head_r.at(2) - b.head_r.at(1))
                    * b.head_ri_scaling**2
                    / (b.head_ro**2 - b.head_ri**2)
                    * (
                        sum(
                            0.5
                            * (
                                b.head_r.at(i - 1)
                                * b.header_wall_temperature[t, b.head_r.at(i - 1)]
                                + b.head_r.at(i)
                                * b.header_wall_temperature[t, b.head_r.at(i)]
                            )
                            for i in range(2, len(b.head_r) + 1)
                        )
                    )
                )

            for head_index_r, head_value_r in enumerate(self.head_r, 1):
                self.head_rindex[head_value_r] = head_index_r

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Discrete Point Mean Temperature for Header",
            )
            def discrete_mean_temperature_header(b, t, r):
                if b.head_rindex[r].value == 1:
                    return b.header_wall_temperature[t, b.head_r.first()]
                else:
                    return (
                        2
                        * (b.head_r.at(2) - b.head_r.at(1))
                        * b.head_ri_scaling**2
                        / (
                            (b.head_r.at(b.head_rindex[r].value) * b.head_ri_scaling)
                            ** 2
                            - b.head_ri**2
                        )
                        * (
                            sum(
                                0.5
                                * (
                                    b.head_r.at(j - 1)
                                    * b.header_wall_temperature[t, b.head_r.at(j - 1)]
                                    + b.head_r.at(j)
                                    * b.header_wall_temperature[t, b.head_r.at(j)]
                                )
                                for j in range(2, b.head_rindex[r].value + 1)
                            )
                        )
                    )

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Thermal Stress at" "Radial Direction for Header",
            )
            def therm_sigma_r_header(b, t, r):
                if r == b.head_r.first() or r == b.head_r.last():
                    return 0
                else:
                    return (
                        0.5
                        * b.Young_modulus_header
                        * b.coefficient_therm_expansion_header
                        / (1 - b.Poisson_ratio_header)
                        * (
                            (1 - b.head_ri**2 / (r * b.head_ri_scaling) ** 2)
                            * (
                                b.mean_temperature_header[t]
                                - b.discrete_mean_temperature_header[t, r]
                            )
                        )
                    )

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Thermal Stress at Circumferential" "Direction for Header",
            )
            def therm_sigma_theta_header(b, t, r):
                r_2 = (r * b.head_ri_scaling) ** 2
                return (
                    0.5
                    * b.Young_modulus_header
                    * b.coefficient_therm_expansion_header
                    / (1 - b.Poisson_ratio_header)
                    * (
                        (1 + b.head_ri**2 / r_2) * b.mean_temperature_header[t]
                        + (1 - b.head_ri**2 / r_2)
                        * b.discrete_mean_temperature_header[t, r]
                        - 2 * b.header_wall_temperature[t, r]
                    )
                )

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Thermal Stress at " " Axial Direction for Header",
            )
            def therm_sigma_z_header(b, t, r):
                return (
                    b.Young_modulus_header
                    * b.coefficient_therm_expansion_header
                    / (1 - b.Poisson_ratio_header)
                    * (b.mean_temperature_header[t] - b.header_wall_temperature[t, r])
                )

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Mechanical Stress " "at Radial Direction for Header",
            )
            def mech_sigma_r_header(b, t, r):
                if r == b.head_r.first():
                    return 1e-6 * (
                        -b.tube.properties[t, b.tube.length_domain.first()].pressure
                    )
                elif r == b.head_r.last():
                    return 1e-6 * (
                        -b.shell.properties[t, b.shell.length_domain.first()].pressure
                    )
                else:
                    return 0.1 * (
                        1e-5
                        * (
                            b.tube.properties[t, b.tube.length_domain.first()].pressure
                            * b.head_ri**2
                            - b.shell.properties[
                                t, b.shell.length_domain.first()
                            ].pressure
                            * b.head_ro**2
                        )
                        / (b.head_ro**2 - b.head_ri**2)
                        + (
                            1e-5
                            * (
                                b.shell.properties[
                                    t, b.shell.length_domain.first()
                                ].pressure
                                - b.tube.properties[
                                    t, b.tube.length_domain.first()
                                ].pressure
                            )
                            * b.head_ri**2
                            * b.head_ro**2
                            / (
                                (r * b.head_ri_scaling) ** 2
                                * (b.head_ro**2 - b.head_ri**2)
                            )
                        )
                    )

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Mechanical Stress" "at Circumferential Direction for Header",
            )
            def mech_sigma_theta_header(b, t, r):
                return 0.1 * (
                    1e-5
                    * (
                        b.tube.properties[t, b.tube.length_domain.first()].pressure
                        * b.head_ri**2
                        - b.shell.properties[t, b.shell.length_domain.first()].pressure
                        * b.head_ro**2
                    )
                    / (b.head_ro**2 - b.head_ri**2)
                    - (
                        1e-5
                        * (
                            b.shell.properties[
                                t, b.shell.length_domain.first()
                            ].pressure
                            - b.tube.properties[
                                t, b.tube.length_domain.first()
                            ].pressure
                        )
                        * b.head_ri**2
                        * b.head_ro**2
                        / (
                            (r * b.head_ri_scaling) ** 2
                            * (b.head_ro**2 - b.head_ri**2)
                        )
                    )
                )

            @self.Expression(
                self.flowsheet().time,
                doc="Mechanical Stress" "at Axial Direction for Header",
            )
            def mech_sigma_z_header(b, t):
                return 0.1 * (
                    1e-5
                    * (
                        b.tube.properties[t, b.tube.length_domain.first()].pressure
                        * b.head_ri**2
                        - b.shell.properties[t, b.shell.length_domain.first()].pressure
                        * b.head_ro**2
                    )
                    / (b.head_ro**2 - b.head_ri**2)
                )

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Principal Structural Stress " "at Radial Direction for Header",
            )
            def sigma_r_header(b, t, r):
                if r == b.head_r.first():
                    return 1e-6 * (
                        -b.tube.properties[t, b.tube.length_domain.first()].pressure
                    )
                elif r == b.head_r.last():
                    return 1e-6 * (
                        -b.shell.properties[t, b.shell.length_domain.first()].pressure
                    )
                else:
                    return b.mech_sigma_r_header[t, r] + b.therm_sigma_r_header[t, r]

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Principal Structural Stress"
                "at Circumferential Direction for Header",
            )
            def sigma_theta_header(b, t, r):
                return (
                    b.mech_sigma_theta_header[t, r] + b.therm_sigma_theta_header[t, r]
                )

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Principal Structural Stress" "at Axial Direction for Header",
            )
            def sigma_z_header(b, t, r):
                return b.mech_sigma_z_header[t] + b.therm_sigma_z_header[t, r]

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Variation Principal Stress"
                "between Radial - Circumferential Directions"
                "for Header",
            )
            def delta_sigma_r_theta_header(b, t, r):
                return abs(b.sigma_r_header[t, r] - b.sigma_theta_header[t, r])

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Variation Principal Stress"
                "between Circumferential-Axial Directions"
                "for Header",
            )
            def delta_sigma_theta_z_header(b, t, r):
                return abs(b.sigma_theta_header[t, r] - b.sigma_z_header[t, r])

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Variation Principal Stress"
                "between Axial-Radial Directions for Header",
            )
            def delta_sigma_z_r_header(b, t, r):
                return abs(b.sigma_z_header[t, r] - b.sigma_r_header[t, r])

            @self.Expression(
                self.flowsheet().time,
                self.head_r,
                doc="Equivalent von Mises Stress for Header",
            )
            def sigma_von_Mises_header(b, t, r):
                return sqrt(
                    b.sigma_r_header[t, r] ** 2
                    + b.sigma_theta_header[t, r] ** 2
                    + b.sigma_z_header[t, r] ** 2
                    - (
                        b.sigma_r_header[t, r] * b.sigma_theta_header[t, r]
                        + b.sigma_r_header[t, r] * b.sigma_z_header[t, r]
                        + b.sigma_theta_header[t, r] * b.sigma_z_header[t, r]
                    )
                )

            # Calculate creep
            # Using property of Steel SA 387 Grade 12
            # (Steel 13CrMo45 or UNS K11757))
            # Data obtained from EN 13445 textbook page 812/838
            # Using Manson-Haferd model
            # log(tr) = f(sigma)*(T-T0) + beta // f(sigma) = b0 + b1*log(sigma)
            # + b2*log(sigma)^2 + b3*log(sigma)^3

            self.creep_a_header = Param(initialize=0.066684094, mutable=True)
            self.creep_b_header = Param(initialize=-0.143434107, mutable=True)
            self.creep_c_header = Param(initialize=0.073764831, mutable=True)
            self.creep_d_header = Param(initialize=-0.013083912, mutable=True)
            self.creep_e_header = Param(initialize=20.32884026, mutable=True)
            self.creep_f_header = Param(initialize=280, mutable=True)

            @self.Expression(
                self.flowsheet().time, self.head_r, doc="Rupture Time for Header"
            )
            def rupture_time_header(b, t, r):
                return 10 ** (
                    (
                        b.creep_a_header
                        + b.creep_b_header * log10(b.sigma_von_Mises_header[t, r])
                        + b.creep_c_header
                        * (log10(b.sigma_von_Mises_header[t, r])) ** 2
                        + b.creep_d_header
                        * (log10(b.sigma_von_Mises_header[t, r])) ** 3
                    )
                    * (b.header_wall_temperature[t, r] - b.creep_f_header)
                    + b.creep_e_header
                )

            # Calculate stress based on EN12952 standard
            # -----------------------------------------------------------------
            # mechanical stress of circumferential direction / thin walled drum
            header_thickness = self.head_thickness
            # mean radius
            r_ms_head = self.head_ri + header_thickness / 2

            # mechanical stress concentration factor
            # thickness of pipe //m
            pipe_th = self.tube_thickness

            # mean diameter of pipe //
            pipe_d = self.tube_do - self.tube_thickness

            # mechanical coefficients
            k_m_A = (
                -1.14 * (pipe_th / header_thickness) ** 2
                - 0.89 * (pipe_th / header_thickness)
                + 1.43
            )
            k_m_B = (
                0.326 * (pipe_th / header_thickness) ** 2
                - 0.59 * (pipe_th / header_thickness)
                + 1.08
            )
            k_m_C = (
                pipe_d
                / (2 * r_ms_head)
                * sqrt((2 * r_ms_head) / (2 * header_thickness))
            )
            k_m_header = 2.2 + exp(k_m_A) * k_m_C**k_m_B

            # thermal stress concentration factor
            k_t_A = pipe_d / (2 * r_ms_head)
            # heat transfer coefficient:3000 for water , 1000 for steam
            k_t_B = 1000
            k_t_header = sqrt(
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
                "for Header (EN 12952-3)",
            )
            def sigma_p(b, t):
                return (
                    0.1
                    * 1e-5
                    * (
                        b.tube.properties[t, b.tube.length_domain.first()].pressure
                        - b.shell.properties[t, b.shell.length_domain.first()].pressure
                    )
                    * r_ms_head
                    / header_thickness
                )

            # thermal stress at circumferential direction
            @self.Expression(
                self.flowsheet().time,
                doc="Thermal Stress at Circumferential"
                "Direction for Header (EN 12952-3)",
            )
            def sigma_t(b, t):
                delta_T = (
                    b.mean_temperature_header[t]
                    - b.header_wall_temperature[t, b.head_r.first()]
                )
                return (
                    b.coefficient_therm_expansion_header
                    * b.Young_modulus_header
                    / (1 - b.Poisson_ratio_header)
                    * delta_T
                )

            # calculate stress at 2 locations at the hole
            # stress at crotch corner P1 and location P2

            # mechanical stress by pressure at crotch corner
            @self.Expression(
                self.flowsheet().time,
                doc="Mechanical Stress at Crotch Corner" "for Header",
            )
            def sigma_p_P1(b, t):
                return b.sigma_p[t] * k_m_header

            # mechanical stress at location P2
            @self.Expression(
                self.flowsheet().time,
                doc="Mechanical Stress at" "Critical Point P2 for Header",
            )
            def sigma_p_P2(b, t):
                return b.sigma_p[t] * k_m_header / 5

            # thermal stress at crotch corner
            @self.Expression(
                self.flowsheet().time, doc="Thermal Stress at Crotch Corner for Header"
            )
            def sigma_t_P1(b, t):
                return b.sigma_t[t] * k_t_header

            # thermal stress at location P2
            @self.Expression(
                self.flowsheet().time,
                doc="Thermal Stress" "at Critical Point P2 for Header",
            )
            def sigma_t_P2(b, t):
                return b.sigma_t[t] * k_t_header

            # total circumferential stress with notch effect
            # crotch corner P1
            @self.Expression(
                self.flowsheet().time,
                doc="Circumferential Stress" "at Crotch Corner for Header",
            )
            def sigma_theta_P1(b, t):
                return b.sigma_p_P1[t] + b.sigma_t_P1[t]

            # location P2
            @self.Expression(
                self.flowsheet().time,
                doc="Circumferential Stress" "at Critical Point P2 for Header",
            )
            def sigma_theta_P2(b, t):
                return b.sigma_p_P2[t] + b.sigma_t_P2[t]

            # total stress with notch effect // f1 - f2
            # crotch corner
            @self.Expression(
                self.flowsheet().time, doc="Total Stress at Crotch Corner for Header"
            )
            def sigma_notch_P1(b, t):
                return (
                    b.sigma_theta_P1[t]
                    + 1e-6 * b.tube.properties[t, b.tube.length_domain.first()].pressure
                )

            # location P2
            @self.Expression(
                self.flowsheet().time, doc="Total Stress at Critial Point P2 for Header"
            )
            def sigma_notch_P2(b, t):
                return (
                    b.sigma_theta_P2[t]
                    + 1e-6 * b.tube.properties[t, b.tube.length_domain.first()].pressure
                )

            # Von Mises equivalent stress
            # VM stress at crotch corner
            @self.Expression(
                self.flowsheet().time,
                doc="Equivalent von Mises Stress" "at Crotch Corner for Header",
            )
            def sigma_eff_P1(b, t):
                sigma_comer_P1 = b.sigma_theta_P1[t]
                p_in = b.tube.properties[t, b.tube.length_domain.first()].pressure
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
                doc="Equivalent von Mises Stress" "at Critical Point P2 for Header",
            )
            def sigma_eff_P2(b, t):
                sigma_comer_P2 = b.sigma_theta_P2[t]
                p_in = b.tube.properties[t, b.tube.length_domain.first()].pressure
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

            # rupture time calculation at crotch corner
            @self.Expression(
                self.flowsheet().time, doc="Rupture Tme at Crotch Corner for Header"
            )
            def rupture_time_crotch_corner(b, t):
                if value(b.sigma_eff_P1[t]) > 10:  # MPa
                    return 10 ** (
                        (
                            b.creep_a_header
                            + b.creep_b_header * log10(b.sigma_eff_P1[t])
                            + b.creep_c_header * (log10(b.sigma_eff_P1[t])) ** 2
                            + b.creep_d_header * (log10(b.sigma_eff_P1[t])) ** 3
                        )
                        * (
                            b.header_wall_temperature[t, b.head_r.first()]
                            - b.creep_f_header
                        )
                        + b.creep_e_header
                    )
                else:
                    return 10**13

            # rupture time calculation at location P2
            @self.Expression(
                self.flowsheet().time,
                doc="Rupture Time at Critical Point P2" "for Header",
            )
            def rupture_time_P2(b, t):
                if value(b.sigma_eff_P2[t]) > 10:  # MPa
                    return 10 ** (
                        (
                            b.creep_a_header
                            + b.creep_b_header * log10(b.sigma_eff_P2[t])
                            + b.creep_c_header * (log10(b.sigma_eff_P2[t])) ** 2
                            + b.creep_d_header * (log10(b.sigma_eff_P2[t])) ** 3
                        )
                        * (
                            b.header_wall_temperature[t, b.head_r.first()]
                            - b.creep_f_header
                        )
                        + b.creep_e_header
                    )
                else:
                    return 10**13

        # ---------------------------------------------------------------------
        # total heat released by shell side fluid assuming even discretization.
        # shell side always in forward direction and the first point is skiped

        @self.Expression(
            self.flowsheet().time, doc="Total Heat Released from Shell Side"
        )
        def total_heat(b, t):
            return (
                -(
                    sum(b.shell_heat[t, x] for x in b.shell.length_domain)
                    - b.shell_heat[t, b.shell.length_domain.first()]
                )
                * b.length_flow_shell
                / b.config.finite_elements
            )

    def set_initial_condition(self):
        if self.config.dynamic is True:
            t0 = self.flowsheet().time.first()
            self.dTdt[:, :, :].value = 0
            self.dTdt[t0, :, :].fix(0)
            if self.config.has_header is True:
                self.head_dTdt[:, :].value = 0
                self.head_dTdt[t0, :].fix(0)
            # no accumulation terms for tube and shell side fluids
            # since currently the fluid flows are modeled as steady-state only

    def initialize_build(
        blk,
        shell_state_args=None,
        tube_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        HeatExchangerCrossFlow1D initialisation routine

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = None).
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

        # ---------------------------------------------------------------------
        # Initialize shell block

        flags_tube = blk.tube.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=tube_state_args
        )
        flags_shell = blk.shell.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=shell_state_args
        )

        init_log.info_high("Initialization Step 1 Complete.")

        phase_s = blk.config.tube_side_water_phase

        # In Step 2, fix tube metal temperatures
        # fix fluid state variables (enthalpy/temperature and pressure)
        # calculate maximum heat duty assuming infinite area and
        # use half of the maximum duty as initial guess
        # to calculate outlet temperature
        if blk.config.flow_type == "co_current":
            mcp_shell = value(
                blk.shell.properties[0, 0].flow_mol * blk.shell.properties[0, 0].cp_mol
            )
            mcp_tube = value(
                blk.tube_inlet.flow_mol[0]
                * blk.tube.properties[0, 0].cp_mol_phase[phase_s]
            )
            tout_max = (
                mcp_tube * value(blk.tube.properties[0, 0].temperature)
                + mcp_shell * value(blk.shell.properties[0, 0].temperature)
            ) / (mcp_tube + mcp_shell)
            q_guess = (
                mcp_tube
                * value(tout_max - value(blk.tube.properties[0, 0].temperature))
                / 2
            )
            temp_out_tube_guess = (
                value(blk.tube.properties[0, 0].temperature) + q_guess / mcp_tube
            )
            temp_out_shell_guess = (
                value(blk.shell.properties[0, 0].temperature) - q_guess / mcp_shell
            )
            if phase_s == "Liq" and temp_out_tube_guess > value(
                blk.tube.properties[0, 0].temperature_sat
            ):
                init_log.info(
                    "Estimated Outlet Liquid Water Temperature"
                    "Exceeds the Saturation Temperature."
                )
                init_log.info(
                    "Estimated Outlet"
                    " Liquid Water Temperature = {}.".format(temp_out_tube_guess)
                )
                init_log.info(
                    "Saturation Temperature at"
                    " Inlet Pressure = {}.".format(
                        value(blk.tube.properties[0, 0].temperature_sat)
                    )
                )
                temp_out_tube_guess = value(
                    0.9 * blk.tube.properties[0, 0].temperature_sat
                    + 0.1 * blk.tube.properties[0, 0].temperature
                )
                init_log.info(
                    "Reset Estimated Outlet Liquid "
                    "Water Ttemperature = {}.".format(temp_out_tube_guess)
                )
        else:
            mcp_shell = value(
                blk.shell.properties[0, 0].flow_mol * blk.shell.properties[0, 0].cp_mol
            )
            mcp_tube = value(
                blk.tube_inlet.flow_mol[0]
                * blk.tube.properties[0, 1].cp_mol_phase[phase_s]
            )
            if mcp_tube < mcp_shell:
                q_guess = (
                    mcp_tube
                    * value(
                        blk.shell.properties[0, 0].temperature
                        - blk.tube.properties[0, 1].temperature
                    )
                    / 2
                )
            else:
                q_guess = (
                    mcp_shell
                    * value(
                        blk.shell.properties[0, 0].temperature
                        - blk.tube.properties[0, 1].temperature
                    )
                    / 2
                )
            temp_out_tube_guess = (
                value(blk.tube.properties[0, 1].temperature) + q_guess / mcp_tube
            )
            temp_out_shell_guess = (
                value(blk.shell.properties[0, 0].temperature) - q_guess / mcp_shell
            )
            if phase_s == "Liq" and temp_out_tube_guess > value(
                blk.tube.properties[0, 1].temperature_sat
            ):
                init_log.info(
                    "Estimated Outlet Liquid Water Temperature"
                    " Exceeds the Saturation temperature."
                )
                init_log.info(
                    "Estimated Outlet Liquid Water "
                    "Temperature = {}.".format(temp_out_tube_guess)
                )
                init_log.info(
                    "Saturation Temperature at Inlet"
                    " Pressure = {}.".format(
                        value(blk.tube.properties[0, 1].temperature_sat)
                    )
                )
                temp_out_tube_guess = value(
                    0.9 * blk.tube.properties[0, 1].temperature_sat
                    + 0.1 * blk.tube.properties[0, 1].temperature
                )
                init_log.info(
                    "Reset Estimated Outlet Liquid Water"
                    " Temperature = {}.".format(temp_out_tube_guess)
                )

        for t in blk.flowsheet().time:
            for z in blk.tube.length_domain:
                if blk.config.flow_type == "co_current":
                    blk.tube_wall_temperature[t, z, :].fix(
                        value(
                            0.05
                            * (
                                (1 - z) * blk.shell.properties[0, 0].temperature
                                + z * temp_out_shell_guess
                            )
                            + 0.95
                            * (
                                (1 - z) * blk.tube.properties[0, 0].temperature
                                + z * temp_out_tube_guess
                            )
                        )
                    )
                else:
                    blk.tube_wall_temperature[t, z, :].fix(
                        value(
                            0.05
                            * (
                                (1 - z) * blk.shell.properties[0, 0].temperature
                                + z * temp_out_shell_guess
                            )
                            + 0.95
                            * (
                                (1 - z) * temp_out_tube_guess
                                + z * blk.tube.properties[0, 1].temperature
                            )
                        )
                    )

        for t in blk.flowsheet().time:
            for z in blk.tube.length_domain:
                blk.tube.properties[t, z].enth_mol.fix(
                    value(blk.tube_inlet.enth_mol[0])
                )
                blk.tube.properties[t, z].pressure.fix(
                    value(blk.tube_inlet.pressure[0])
                )

        for t in blk.flowsheet().time:
            for z in blk.shell.length_domain:
                blk.shell.properties[t, z].temperature.fix(
                    value(blk.shell.properties[0, 0].temperature)
                )
                blk.shell.properties[t, z].pressure.fix(
                    value(blk.shell.properties[0, 0].pressure)
                )

        blk.heat_conduction_eqn.deactivate()
        blk.inner_wall_bc_eqn.deactivate()
        blk.outer_wall_bc_eqn.deactivate()
        blk.extra_at_inner_wall_eqn.deactivate()
        blk.extra_at_outer_wall_eqn.deactivate()
        blk.deltaP_tube_eqn.deactivate()
        blk.deltaP_shell_eqn.deactivate()
        blk.heat_tube_eqn.deactivate()
        blk.heat_shell_eqn.deactivate()

        if blk.config.has_header is True:
            blk.head_heat_conduction_eqn.deactivate()
            blk.head_inner_wall_bc_eqn.deactivate()
            blk.head_outer_wall_bc_eqn.deactivate()
            blk.head_extra_at_inner_wall_eqn.deactivate()
            blk.head_extra_at_outer_wall_eqn.deactivate()
            blk.head_hout_eqn.deactivate()
            blk.head_Ra_number_eqn.deactivate()
            blk.head_Nu_number_shell_eqn.deactivate()
            blk.head_h_free_conv_eqn.deactivate()
            blk.head_velocity_tube_eqn.deactivate()
            blk.head_N_Re_tube_eqn.deactivate()
            blk.head_N_Nu_tube_eqn.deactivate()
            blk.head_hin_eqn.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # In Step 3, unfix fluid state variables
        # (enthalpy/temperature and pressure)
        # keep the inlet state variables fixed,
        # otherwise, the degree of freedom > 0
        for t in blk.flowsheet().time:
            for z in blk.tube.length_domain:
                blk.tube.properties[t, z].enth_mol.unfix()
                blk.tube.properties[t, z].pressure.unfix()
            if blk.config.flow_type == "co_current":
                blk.tube.properties[t, 0].enth_mol.fix(
                    value(blk.tube_inlet.enth_mol[0])
                )
                blk.tube.properties[t, 0].pressure.fix(
                    value(blk.tube_inlet.pressure[0])
                )
            else:
                blk.tube.properties[t, 1].enth_mol.fix(
                    value(blk.tube_inlet.enth_mol[0])
                )
                blk.tube.properties[t, 1].pressure.fix(
                    value(blk.tube_inlet.pressure[0])
                )

        for t in blk.flowsheet().time:
            for z in blk.shell.length_domain:
                blk.shell.properties[t, z].temperature.unfix()
                blk.shell.properties[t, z].pressure.unfix()
            blk.shell.properties[t, 0].temperature.fix(
                value(blk.shell_inlet.temperature[0])
            )
            blk.shell.properties[t, 0].pressure.fix(value(blk.shell_inlet.pressure[0]))

        blk.deltaP_tube_eqn.activate()
        blk.deltaP_shell_eqn.activate()
        blk.heat_tube_eqn.activate()
        blk.heat_shell_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        blk.tube_wall_temperature[:, :, :].unfix()
        blk.heat_conduction_eqn.activate()
        blk.inner_wall_bc_eqn.activate()
        blk.outer_wall_bc_eqn.activate()
        blk.extra_at_inner_wall_eqn.activate()
        blk.extra_at_outer_wall_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))

        if blk.config.has_header is True:
            blk.head_heat_conduction_eqn.activate()
            blk.head_inner_wall_bc_eqn.activate()
            blk.head_outer_wall_bc_eqn.activate()
            blk.head_extra_at_inner_wall_eqn.activate()
            blk.head_extra_at_outer_wall_eqn.activate()
            blk.head_hout_eqn.activate()
            blk.head_Ra_number_eqn.activate()
            blk.head_Nu_number_shell_eqn.activate()
            blk.head_h_free_conv_eqn.activate()
            blk.head_velocity_tube_eqn.activate()
            blk.head_N_Re_tube_eqn.activate()
            blk.head_N_Nu_tube_eqn.activate()
            blk.head_hin_eqn.activate()
            temp_head_in = (
                value(
                    blk.tube.properties[0, blk.tube.length_domain.first()].temperature
                )
                - 0.1
            )
            temp_head_out = temp_head_in - 0.5
            for x in blk.head_r:
                blk.header_wall_temperature[:, x].value = temp_head_in + (
                    temp_head_out - temp_head_in
                ) / (blk.head_r.last() - blk.head_r.first()) * (x - blk.head_r.first())

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
            init_log.info_high(
                "Initialization Step 5 {}.".format(idaeslog.condition(res))
            )

        blk.tube.release_state(flags_tube)
        blk.shell.release_state(flags_shell)
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        for i, c in self.heat_tube_eqn.items():
            sf = iscale.get_scaling_factor(self.tube_heat[i], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for i, c in self.heat_shell_eqn.items():
            sf = iscale.get_scaling_factor(self.shell_heat[i], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf, overwrite=False)
