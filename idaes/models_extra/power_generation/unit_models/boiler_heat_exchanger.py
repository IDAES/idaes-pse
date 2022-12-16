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
Power Plant IDAES heat exchanger model.

The boiler heat exchanger model consist of a cross flow shell and tube hx
that can be used for any of the boiler components, such as economizer,
reheater or superheaters (primary, secondary, etc).
The model includes shell and tube rigorous heat transfer calculations and
pressure drop calculations for shell side. Note that this model assumes no
phase transitions (if user requires phase transitions, they need a general
model)

The main config arguments:
    - delta T method: counter-current or co-current
    - tube_arrangement: in-line or staggered
    - has radiation: True if model is used as a reheater or superheater unit
        Gas emissivity calculated (Gas temperature above 700 K)

General assumtpions:
    - SI units (consistent with prop pack)
    - heat transfer calc U = f(Nu, Re, Pr)
    - Pressure drop tube and shell side (friction factor calc.)

"""

__author__ = "Boiler subsystem team (J Ma, M Zamarripa)"
__version__ = "1.0.0"

# Import Python libraries
from enum import Enum

# Import Pyomo libraries
from pyomo.common.config import ConfigValue, In

# Additional import for the unit operation
from pyomo.environ import (
    value,
    Var,
    Param,
    exp,
    sqrt,
    log,
    PositiveReals,
    NonNegativeReals,
    Reference,
    units as pyunits,
)

# Import IDAES cores
from idaes.core import declare_process_block_class

from idaes.core.util.misc import add_object_reference
from idaes.core.util.constants import Constants as c
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.exceptions import ConfigurationError

from idaes.models.unit_models.heat_exchanger import (
    HeatExchangerData,
    delta_temperature_lmtd_callback,
    delta_temperature_lmtd2_callback,
    delta_temperature_lmtd3_callback,
    delta_temperature_amtd_callback,
    delta_temperature_underwood_callback,
    HeatExchangerFlowPattern,
)


# Set up logger
_log = idaeslog.getLogger(__name__)


class TubeArrangement(Enum):
    inLine = 0
    staggered = 1


@declare_process_block_class("BoilerHeatExchanger")
class BoilerHeatExchangerData(HeatExchangerData):
    CONFIG = HeatExchangerData.CONFIG(implicit=True)

    CONFIG.declare(
        "tube_arrangement",
        ConfigValue(
            default=TubeArrangement.inLine,
            domain=In(TubeArrangement),
            description="tube configuration",
            doc="Tube arrangement could be in-line and staggered",
        ),
    )
    CONFIG.declare(
        "cold_side_water_phase",
        ConfigValue(
            default="Liq",
            domain=In(["Liq", "Vap"]),
            description="side 1 water phase",
            doc="Define water phase for property calls",
        ),
    )
    CONFIG.declare(
        "has_radiation",
        ConfigValue(
            default=False,
            domain=In([False, True]),
            description="Has hot side gas radiation",
            doc="Define if hot side gas radiation is to be considered",
        ),
    )

    def _process_config(self):
        """
        Check for boiler specific config arguments.
        """
        config = self.config

        if config.flow_pattern == HeatExchangerFlowPattern.crossflow:
            raise ConfigurationError("Boiler heat exchanger does not support crossflow")

    def build(self):
        """
        Build method for Boiler heat exchanger model

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()
        self._process_config()
        self.deltaT_1 = Reference(self.delta_temperature_in)
        self.deltaT_2 = Reference(self.delta_temperature_out)

        self._set_geometry()
        self.cold_side_fluid_phase = self.config.cold_side_water_phase
        # Construct performance equations
        self._make_performance()

    def _set_geometry(self):
        """
        Define the geometry of the unit as necessary, and link to holdup volume

        Args:
            None

        Returns:
            None
        """
        # Elevation difference (outlet - inlet) for static pressure calculation
        self.delta_elevation = Var(
            initialize=0,
            within=NonNegativeReals,
            doc="Elevation increase used for static pressure calculation",
            units=pyunits.m,
        )

        # Number of tube columns in the cross section plane
        # perpendicular to shell side fluid flow (y direction)
        self.tube_ncol = Var(
            initialize=10.0, within=PositiveReals, doc="Number of tube columns"
        )

        # Number of tube rows in the direction of shell side
        # fluid flow (x direction)
        self.tube_nrow = Var(
            initialize=10.0, within=PositiveReals, doc="Number of tube rows"
        )

        # Number of inlet tube rows
        self.nrow_inlet = Var(
            initialize=1, within=PositiveReals, doc="Number of inlet tube rows"
        )

        # Length of a tube in z direction for each path
        self.tube_length = Var(
            initialize=5.0, within=PositiveReals, doc="Tube length", units=pyunits.m
        )

        # Inner diameter of tubes
        self.tube_di = Var(
            initialize=0.05,
            within=PositiveReals,
            doc="Inner diameter of tube",
            units=pyunits.m,
        )

        # Thickness of tube
        self.tube_thickness = Var(
            initialize=0.005,
            within=PositiveReals,
            doc="Tube thickness",
            units=pyunits.m,
        )

        # Pitch of tubes between two neighboring columns (in y direction).
        # Always greater than tube outside diameter
        self.pitch_y = Var(
            initialize=0.1,
            within=PositiveReals,
            doc="Pitch between two neighboring columns",
            units=pyunits.m,
        )

        # Pitch of tubes between two neighboring rows (in x direction).
        # Always greater than tube outside diameter
        self.pitch_x = Var(
            initialize=0.1,
            within=PositiveReals,
            doc="Pitch between two neighboring rows",
            units=pyunits.m,
        )

        # Tube outside diameter
        @self.Expression(doc="Outside diameter of tube")
        def do_tube(b):
            return b.tube_di + b.tube_thickness * 2.0

        if self.config.has_radiation is True:
            # Mean beam length for radiation
            @self.Expression(doc="Mean beam length")
            def mbl(b):
                return 3.6 * (
                    b.pitch_x * b.pitch_y / c.pi / b.do_tube - b.do_tube / 4.0
                )

            # Mean beam length for radiation divided by sqrt(2)
            @self.Expression(doc="Mean beam length")
            def mbl_div2(b):
                return b.mbl / sqrt(2.0)

            # Mean beam length for radiation multiplied by sqrt(2)
            @self.Expression(doc="Mean beam length")
            def mbl_mul2(b):
                return b.mbl * sqrt(2.0)

        # Number of 180 degree bends for the tube
        @self.Expression(doc="Nbend_tube")
        def nbend_tube(b):
            return b.tube_nrow / b.nrow_inlet

        # Total flow area on tube side
        @self.Expression(doc="Total flow area on tube side")
        def area_flow_tube(b):
            return 0.25 * c.pi * b.tube_di**2.0 * b.tube_ncol * b.nrow_inlet

        # Total flow area on shell side
        @self.Expression(doc="Total flow area on shell side")
        def area_flow_shell(b):
            return b.tube_length * (b.pitch_y - b.do_tube) * b.tube_ncol

        # Total heat transfer area based on outside diameter
        @self.Constraint(doc="Total heat transfer area based on tube outside diameter")
        def area_eqn(b):
            return (
                b.area == c.pi * b.do_tube * b.tube_length * b.tube_ncol * b.tube_nrow
            )

        # Ratio of pitch_x/do_tube
        @self.Expression(
            doc="Ratio of pitch in x " "direction to tube outside diameter"
        )
        def pitch_x_to_do(b):
            return b.pitch_x / b.do_tube

        # Ratio of pitch_y/do_tube
        @self.Expression(
            doc="Ratio of pitch in y " "direction to tube outside diameter"
        )
        def pitch_y_to_do(b):
            return b.pitch_y / b.do_tube

        if self.config.has_holdup is True:
            add_object_reference(self, "volume_cold_side", self.cold_side.volume)
            add_object_reference(self, "volume_hot_side", self.hot_side.volume)
            # Total tube side valume
            self.Constraint(doc="Total tube side volume")

            def volume_cold_side_eqn(b):
                return b.volume_cold_side == (
                    0.25
                    * c.pi
                    * b.tube_di**2.0
                    * b.tube_length
                    * b.tube_ncol
                    * b.tube_nrow
                )

            # Total shell side valume
            self.Constraint(doc="Total shell side volume")

            def volume_hot_side_eqn(b):
                return (
                    b.volume_hot_side
                    == b.tube_ncol * b.pitch_y * b.tube_length * b.tube_nrow * b.pitch_x
                    - 0.25
                    * c.pi
                    * b.do_tube**2.0
                    * b.tube_length
                    * b.tube_ncol
                    * b.tube_nrow
                )

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.

        Args:
            None

        Returns:
            None
        """
        # Set references to balance terms at unit level
        add_object_reference(self, "heat_duty", self.cold_side.heat)
        if self.config.cold_side.has_pressure_change is True:
            add_object_reference(self, "deltaP_tube", self.cold_side.deltaP)
        if self.config.hot_side.has_pressure_change is True:
            add_object_reference(self, "deltaP_shell", self.hot_side.deltaP)

        # Performance parameters and variables
        # Wall thermal conductivity
        self.therm_cond_wall = Param(
            initialize=43.0,
            within=PositiveReals,
            doc="Thermal conductivity of the wall",
            units=pyunits.W / pyunits.m / pyunits.K,
        )

        # Loss coefficient for a 180 degree bend (u-turn),
        # usually related to radius to inside diameter ratio
        self.k_loss_uturn = Param(
            initialize=0.5,
            within=PositiveReals,
            mutable=True,
            doc="Loss coefficient of a tube u-turn",
        )

        # Heat transfer resistance due to the fouling on tube side
        # (typical boiler hx)
        self.tube_r_fouling = Param(
            initialize=0.00017,
            within=NonNegativeReals,
            mutable=True,
            doc="Fouling resistance on tube side",
            units=pyunits.K * pyunits.m**2 * pyunits.W**-1,
        )

        # Heat transfer resistance due to the fouling on shell side
        self.shell_r_fouling = Param(
            initialize=0.0008,
            within=NonNegativeReals,
            mutable=True,
            doc="Fouling resistance on tube side",
            units=pyunits.K * pyunits.m**2 * pyunits.W**-1,
        )

        # Correction factor for overall heat transfer coefficient
        self.fcorrection_htc = Var(
            initialize=1.0, within=NonNegativeReals, doc="Correction factor for HTC"
        )

        # Correction factor for tube side pressure drop due to friction
        self.fcorrection_dp_tube = Var(
            initialize=1.0, doc="Correction factor for tube side pressure drop"
        )

        # Correction factor for shell side pressure drop due to friction
        self.fcorrection_dp_shell = Var(
            initialize=1.0, doc="Correction factor for shell side pressure drop"
        )

        if self.config.has_radiation is True:
            # Shell side wall emissivity, converted from parameter to variable
            self.emissivity_wall = Var(initialize=0.7, doc="Shell side wall emissivity")
            # Gas emissivity at mbl
            self.gas_emissivity = Var(
                self.flowsheet().time,
                initialize=0.5,
                doc="Emissivity at given mean beam length",
            )

            # Gas emissivity at mbl/sqrt(2)
            self.gas_emissivity_div2 = Var(
                self.flowsheet().time,
                initialize=0.4,
                doc="Emissivity at mean beam length divided by sqrt of 2",
            )

            # Gas emissivity at mbl*sqrt(2)
            self.gas_emissivity_mul2 = Var(
                self.flowsheet().time,
                initialize=0.6,
                doc="Emissivity at mean beam length multiplied by sqrt of 2",
            )

            # Gray fraction of gas in entire spectrum
            self.gas_gray_fraction = Var(
                self.flowsheet().time,
                initialize=0.5,
                doc="Gray fraction of gas in entire spectrum",
            )

            # Gas-surface radiation exchange factor for shell side wall
            self.frad_gas_shell = Var(
                self.flowsheet().time,
                initialize=0.5,
                doc="Gas-surface radiation exchange " "factor for shell side wall",
            )

            # Shell side equivalent convective heat transfer coefficient
            # due to radiation
            self.hconv_shell_rad = Var(
                self.flowsheet().time,
                initialize=100.0,
                doc="Shell convective heat transfer coefficient due to radiation",
                units=pyunits.W / pyunits.m**2 / pyunits.K,
            )

        # Tube side convective heat transfer coefficient
        self.hconv_tube = Var(
            self.flowsheet().time,
            initialize=100.0,
            doc="Tube side convective heat transfer coefficient",
            units=pyunits.W / pyunits.m**2 / pyunits.K,
        )

        # Shell side convective heat transfer coefficient due to convection
        self.hconv_shell_conv = Var(
            self.flowsheet().time,
            initialize=100.0,
            doc="Shell side convective heat transfer coefficient due to convection",
            units=pyunits.W / pyunits.m**2 / pyunits.K,
        )

        # Total shell side convective heat transfer coefficient
        # including convection and radiation
        self.hconv_shell_total = Var(
            self.flowsheet().time,
            initialize=150.0,
            doc="Total shell side convective heat transfer coefficient",
            units=pyunits.W / pyunits.m**2 / pyunits.K,
        )

        # Heat conduction resistance of tube wall
        self.rcond_wall = Var(
            initialize=1.0,
            doc="Heat conduction resistance of wall",
            units=pyunits.m**2 * pyunits.K / pyunits.W,
        )

        if self.config.has_radiation is True:
            # Constraints for gas emissivity
            @self.Constraint(self.flowsheet().time, doc="Gas emissivity")
            def gas_emissivity_eqn(b, t):
                # This is a surrogate model, so need to do units manually
                X1 = (
                    (
                        b.hot_side.properties_in[t].temperature
                        + b.hot_side.properties_out[t].temperature
                    )
                    / 2
                    / pyunits.K
                )
                X2 = b.mbl / pyunits.m
                X3 = b.hot_side.properties_in[t].pressure / pyunits.Pa
                X4 = b.hot_side.properties_in[t].mole_frac_comp["CO2"]
                X5 = b.hot_side.properties_in[t].mole_frac_comp["H2O"]
                X6 = b.hot_side.properties_in[t].mole_frac_comp["O2"]

                # Surrogate model fitted using rigorous calc. - 500 samples
                # Wide operating range:
                #      X1: 700 – 1500    (Gas Temperature)
                #      X2: 0.2 – 1       (Mean beam length)
                #      X3: 79000-102000  (pressure in Pa)
                #      X4: 0.12-0.16     (mol frac CO2)
                #      X5: 0.075-0.15    (mol frac H2O)
                #      X6: 0.01-0.07     (mol frac O2)

                return b.gas_emissivity[t] == (
                    -0.116916606892e-003 * X1
                    - 0.29111124038936179309056e-001 * X2
                    + 0.50509651230704191577346e-006 * X3
                    + 1.1844222822155641150488 * X4
                    - 0.64720757767102773949652e-001 * X5
                    - 0.35853593221454795048064e-001 * X6
                    + 0.12227919099126832724878 * log(X1)
                    + 0.45102118316418124410738e-001 * log(X2)
                    + 0.33111863480179408447679e-001 * log(X3)
                    + 0.17674928397780117345084e-001 * log(X5)
                    - 0.12541139396423576016226e-001 * exp(X2)
                    - 0.90251708836308952577099 * exp(X4)
                    + 0.32447078857791738538963e-002 * X2**2
                    - 0.31332075610864829615706e-004 * X1 * X2
                    - 0.54639645449809960433102e-009 * X1 * X3
                    - 0.19721467902854980460033e-003 * X1 * X5
                    + 0.45275517692290622763507e-004 * X1 * X6
                    + 0.75458754990630776904396e-006 * X2 * X3
                    + 0.39691751689931338564765e-001 * X2 * X4
                    + 0.73169514231974708273754 * X2 * X5
                    - 0.35852614507684822664491e-001 * X2 * X6
                    + 0.39743672195685803976177e-005 * X3 * X5
                    + 0.58802879141883679897383e-008 * (X1 * X2) ** 2
                    - 1.2994610452829884472692 * (X2 * X5) ** 2
                )

            # Constraints for gas emissivity at mbl/sqrt(2)
            @self.Constraint(
                self.flowsheet().time, doc="Gas emissivity at a lower mean beam length"
            )
            def gas_emissivity_div2_eqn(b, t):
                # This is a surrogate model, so need to do units manually
                X1 = (
                    (
                        b.hot_side.properties_in[t].temperature
                        + b.hot_side.properties_out[t].temperature
                    )
                    / 2
                    / pyunits.K
                )
                X2 = b.mbl_div2 / pyunits.m
                X3 = b.hot_side.properties_in[t].pressure / pyunits.Pa
                X4 = b.hot_side.properties_in[t].mole_frac_comp["CO2"]
                X5 = b.hot_side.properties_in[t].mole_frac_comp["H2O"]
                X6 = b.hot_side.properties_in[t].mole_frac_comp["O2"]

                # Surrogate model fitted using rigorous calc. - 500 samples
                # Wide operating range:
                #       X1: 700 – 1500    (Gas Temperature)
                #       X2: 0.2 – 1       (Mean beam length)
                #       X3: 79000-102000  (pressure in Pa)
                #       X4: 0.12-0.16     (mol frac CO2)
                #       X5: 0.075-0.15    (mol frac H2O)
                #       X6: 0.01-0.07     (mol frac O2)
                return b.gas_emissivity_div2[t] == (
                    -0.116916606892e-003 * X1
                    - 0.29111124038936179309056e-001 * X2
                    + 0.50509651230704191577346e-006 * X3
                    + 1.1844222822155641150488 * X4
                    - 0.64720757767102773949652e-001 * X5
                    - 0.35853593221454795048064e-001 * X6
                    + 0.12227919099126832724878 * log(X1)
                    + 0.45102118316418124410738e-001 * log(X2)
                    + 0.33111863480179408447679e-001 * log(X3)
                    + 0.17674928397780117345084e-001 * log(X5)
                    - 0.12541139396423576016226e-001 * exp(X2)
                    - 0.90251708836308952577099 * exp(X4)
                    + 0.32447078857791738538963e-002 * X2**2
                    - 0.31332075610864829615706e-004 * X1 * X2
                    - 0.54639645449809960433102e-009 * X1 * X3
                    - 0.19721467902854980460033e-003 * X1 * X5
                    + 0.45275517692290622763507e-004 * X1 * X6
                    + 0.75458754990630776904396e-006 * X2 * X3
                    + 0.39691751689931338564765e-001 * X2 * X4
                    + 0.73169514231974708273754 * X2 * X5
                    - 0.35852614507684822664491e-001 * X2 * X6
                    + 0.39743672195685803976177e-005 * X3 * X5
                    + 0.58802879141883679897383e-008 * (X1 * X2) ** 2
                    - 1.2994610452829884472692 * (X2 * X5) ** 2
                )

            # Constraints for gas emissivity at mbl*sqrt(2)
            @self.Constraint(
                self.flowsheet().time, doc="Gas emissivity at a higher mean beam length"
            )
            def gas_emissivity_mul2_eqn(b, t):
                # This is a surrogate model, so need to do units manually
                X1 = (
                    (
                        b.hot_side.properties_in[t].temperature
                        + b.hot_side.properties_out[t].temperature
                    )
                    / 2
                    / pyunits.K
                )
                X2 = b.mbl_mul2 / pyunits.m
                X3 = b.hot_side.properties_in[t].pressure / pyunits.Pa
                X4 = b.hot_side.properties_in[t].mole_frac_comp["CO2"]
                X5 = b.hot_side.properties_in[t].mole_frac_comp["H2O"]
                X6 = b.hot_side.properties_in[t].mole_frac_comp["O2"]

                # Surrogate model fitted using rigorous calc. 500 samples
                # Wide operating range:
                #       X1: 700 – 1500    (Gas Temperature)
                #       X2: 0.2 – 1       (Mean beam length)
                #       X3: 79000-102000  (pressure in Pa)
                #       X4: 0.12-0.16     (mol frac CO2)
                #       X5: 0.075-0.15    (mol frac H2O)
                #       X6: 0.01-0.07     (mol frac O2)
                return b.gas_emissivity_mul2[t] == (
                    -0.116916606892e-003 * X1
                    - 0.29111124038936179309056e-001 * X2
                    + 0.50509651230704191577346e-006 * X3
                    + 1.1844222822155641150488 * X4
                    - 0.64720757767102773949652e-001 * X5
                    - 0.35853593221454795048064e-001 * X6
                    + 0.12227919099126832724878 * log(X1)
                    + 0.45102118316418124410738e-001 * log(X2)
                    + 0.33111863480179408447679e-001 * log(X3)
                    + 0.17674928397780117345084e-001 * log(X5)
                    - 0.12541139396423576016226e-001 * exp(X2)
                    - 0.90251708836308952577099 * exp(X4)
                    + 0.32447078857791738538963e-002 * X2**2
                    - 0.31332075610864829615706e-004 * X1 * X2
                    - 0.54639645449809960433102e-009 * X1 * X3
                    - 0.19721467902854980460033e-003 * X1 * X5
                    + 0.45275517692290622763507e-004 * X1 * X6
                    + 0.75458754990630776904396e-006 * X2 * X3
                    + 0.39691751689931338564765e-001 * X2 * X4
                    + 0.73169514231974708273754 * X2 * X5
                    - 0.35852614507684822664491e-001 * X2 * X6
                    + 0.39743672195685803976177e-005 * X3 * X5
                    + 0.58802879141883679897383e-008 * (X1 * X2) ** 2
                    - 1.2994610452829884472692 * (X2 * X5) ** 2
                )

            # fraction of gray gas spectrum
            @self.Constraint(self.flowsheet().time, doc="Fraction of gray gas spectrum")
            def gas_gray_fraction_eqn(b, t):
                return (
                    b.gas_gray_fraction[t]
                    * (2 * b.gas_emissivity_div2[t] - b.gas_emissivity_mul2[t])
                    == b.gas_emissivity_div2[t] ** 2
                )

            # gas-surface radiation exchange factor
            # between gas and shell side wall
            @self.Constraint(
                self.flowsheet().time,
                doc="Gas-surface radiation exchange "
                "factor between gas and shell side wall",
            )
            def frad_gas_shell_eqn(b, t):
                return (
                    b.frad_gas_shell[t]
                    * (
                        (1 / b.emissivity_wall - 1) * b.gas_emissivity[t]
                        + b.gas_gray_fraction[t]
                    )
                    == b.gas_gray_fraction[t] * b.gas_emissivity[t]
                )

            # equivalent convective heat transfer coefficent due to radiation
            @self.Constraint(
                self.flowsheet().time,
                doc="Equivalent convective heat transfer "
                "coefficent due to radiation",
            )
            def hconv_shell_rad_eqn(b, t):
                return b.hconv_shell_rad[t] == c.stefan_constant * b.frad_gas_shell[
                    t
                ] * (
                    (
                        b.hot_side.properties_in[t].temperature
                        + b.hot_side.properties_out[t].temperature
                    )
                    / 2
                    + b.cold_side.properties_in[t].temperature
                ) * (
                    (
                        (
                            b.hot_side.properties_in[t].temperature
                            + b.hot_side.properties_out[t].temperature
                        )
                        / 2
                    )
                    ** 2
                    + b.cold_side.properties_in[t].temperature ** 2
                )

        # Tube side heat transfer coefficient and pressure drop
        # -----------------------------------------------------
        # Velocity on tube side
        self.v_tube = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Velocity on tube side",
            units=pyunits.m / pyunits.s,
        )

        # Reynalds number on tube side
        self.N_Re_tube = Var(
            self.flowsheet().time,
            initialize=10000.0,
            doc="Reynolds number on tube side",
        )
        if self.config.cold_side.has_pressure_change is True:
            # Friction factor on tube side
            self.friction_factor_tube = Var(
                self.flowsheet().time,
                initialize=1.0,
                doc="Friction factor on tube side",
            )

            # Pressure drop due to friction on tube side
            self.deltaP_tube_friction = Var(
                self.flowsheet().time,
                initialize=-10.0,
                doc="Pressure drop due to friction on tube side",
                units=pyunits.Pa,
            )

            # Pressure drop due to 180 degree turn on tube side
            self.deltaP_tube_uturn = Var(
                self.flowsheet().time,
                initialize=-10.0,
                doc="Pressure drop due to u-turn on tube side",
                units=pyunits.Pa,
            )

        # Prandtl number on tube side
        self.N_Pr_tube = Var(
            self.flowsheet().time, initialize=1, doc="Prandtl number on tube side"
        )

        # Nusselt number on tube side
        self.N_Nu_tube = Var(
            self.flowsheet().time, initialize=1, doc="Nusselts number on tube side"
        )

        # Velocity equation
        @self.Constraint(self.flowsheet().time, doc="Tube side velocity equation")
        def v_tube_eqn(b, t):
            return (
                b.v_tube[t]
                * b.area_flow_tube
                * b.cold_side.properties_in[t].dens_mol_phase[
                    self.cold_side_fluid_phase
                ]
                == b.cold_side.properties_in[t].flow_mol
            )

        # Reynolds number
        @self.Constraint(
            self.flowsheet().time, doc="Reynolds number equation on tube side"
        )
        def N_Re_tube_eqn(b, t):
            return (
                b.N_Re_tube[t]
                * b.cold_side.properties_in[t].visc_d_phase[self.cold_side_fluid_phase]
                == b.tube_di
                * b.v_tube[t]
                * b.cold_side.properties_in[t].dens_mass_phase[
                    self.cold_side_fluid_phase
                ]
            )

        if self.config.cold_side.has_pressure_change is True:
            # Friction factor
            @self.Constraint(
                self.flowsheet().time, doc="Darcy friction factor on tube side"
            )
            def friction_factor_tube_eqn(b, t):
                return (
                    b.friction_factor_tube[t] * b.N_Re_tube[t] ** 0.25
                    == 0.3164 * b.fcorrection_dp_tube
                )

            # Pressure drop due to friction
            @self.Constraint(
                self.flowsheet().time, doc="Pressure drop due to friction on tube side"
            )
            def deltaP_tube_friction_eqn(b, t):
                return (
                    b.deltaP_tube_friction[t]
                    == -0.5
                    * b.cold_side.properties_in[t].dens_mass_phase[
                        self.cold_side_fluid_phase
                    ]
                    * b.v_tube[t] ** 2
                    * b.friction_factor_tube[t]
                    * b.tube_length
                    * b.tube_nrow
                    / b.tube_di
                    / b.nrow_inlet
                )

            # Pressure drop due to u-turn
            @self.Constraint(
                self.flowsheet().time, doc="Pressure drop due to u-turn on tube side"
            )
            def deltaP_tube_uturn_eqn(b, t):
                return (
                    b.deltaP_tube_uturn[t]
                    == -0.5
                    * b.cold_side.properties_in[t].dens_mass_phase[
                        self.cold_side_fluid_phase
                    ]
                    * b.v_tube[t] ** 2
                    * b.k_loss_uturn
                )

            # Total pressure drop on tube side
            @self.Constraint(
                self.flowsheet().time, doc="Total pressure drop on tube side"
            )
            def deltaP_tube_eqn(b, t):
                return (
                    b.deltaP_tube[t]
                    == b.deltaP_tube_friction[t]
                    + b.deltaP_tube_uturn[t]
                    - b.delta_elevation
                    * c.acceleration_gravity
                    * (
                        b.cold_side.properties_in[t].dens_mass_phase[
                            self.cold_side_fluid_phase
                        ]
                        + b.cold_side.properties_out[t].dens_mass_phase[
                            self.cold_side_fluid_phase
                        ]
                    )
                    / 2.0
                )

        # Prandtl number
        @self.Constraint(
            self.flowsheet().time, doc="Prandtl number equation on tube side"
        )
        def N_Pr_tube_eqn(b, t):
            return (
                b.N_Pr_tube[t]
                * b.cold_side.properties_in[t].therm_cond_phase[
                    self.cold_side_fluid_phase
                ]
                * b.cold_side.properties_in[t].mw
                == b.cold_side.properties_in[t].cp_mol_phase[self.cold_side_fluid_phase]
                * b.cold_side.properties_in[t].visc_d_phase[self.cold_side_fluid_phase]
            )

        # Nusselts number
        @self.Constraint(
            self.flowsheet().time, doc="Nusselts number equation on tube side"
        )
        def N_Nu_tube_eqn(b, t):
            return (
                b.N_Nu_tube[t]
                == 0.023 * b.N_Re_tube[t] ** 0.8 * abs(b.N_Pr_tube[t]) ** 0.4
            )

        # Heat transfer coefficient
        @self.Constraint(
            self.flowsheet().time,
            doc="Convective heat transfer " "coefficient equation on tube side",
        )
        def hconv_tube_eqn(b, t):
            return (
                b.hconv_tube[t] * self.tube_di / 1000
                == b.N_Nu_tube[t]
                * b.cold_side.properties_in[t].therm_cond_phase[
                    self.cold_side_fluid_phase
                ]
                / 1000
            )

        # Pressure drop and heat transfer coefficient on shell side
        # ----------------------------------------------------------
        # Tube arrangement factor
        if self.config.tube_arrangement == TubeArrangement.inLine:
            self.f_arrangement = Param(
                initialize=0.788, doc="In-line tube arrangement factor"
            )
        elif self.config.tube_arrangement == TubeArrangement.staggered:
            self.f_arrangement = Param(
                initialize=1.0, doc="Staggered tube arrangement factor"
            )
        else:
            raise Exception("tube arrangement type not supported")
        # Velocity on shell side
        self.v_shell = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Velocity on shell side",
            units=pyunits.m / pyunits.s,
        )

        # Reynalds number on shell side
        self.N_Re_shell = Var(
            self.flowsheet().time,
            initialize=10000.0,
            doc="Reynolds number on shell side",
        )

        # Friction factor on shell side
        self.friction_factor_shell = Var(
            self.flowsheet().time, initialize=1.0, doc="Friction factor on shell side"
        )

        # Prandtl number on shell side
        self.N_Pr_shell = Var(
            self.flowsheet().time, initialize=1, doc="Prandtl number on shell side"
        )

        # Nusselt number on shell side
        self.N_Nu_shell = Var(
            self.flowsheet().time, initialize=1, doc="Nusselts number on shell side"
        )

        # Velocity equation on shell side
        @self.Constraint(self.flowsheet().time, doc="Velocity on shell side")
        def v_shell_eqn(b, t):
            return b.v_shell[t] * b.hot_side.properties_in[t].dens_mol_phase[
                "Vap"
            ] * b.area_flow_shell == sum(
                b.hot_side.properties_in[t].flow_mol_comp[j]
                for j in b.hot_side.properties_in[t].params.component_list
            )

        # Reynolds number
        @self.Constraint(
            self.flowsheet().time, doc="Reynolds number equation on shell side"
        )
        def N_Re_shell_eqn(b, t):
            return b.N_Re_shell[t] * b.hot_side.properties_in[
                t
            ].visc_d == b.do_tube * b.v_shell[t] * b.hot_side.properties_in[
                t
            ].dens_mol_phase[
                "Vap"
            ] * sum(
                b.hot_side.properties_in[t].mw_comp[c]
                * b.hot_side.properties_in[t].mole_frac_comp[c]
                for c in b.hot_side.properties_in[t].params.component_list
            )

        if self.config.hot_side.has_pressure_change is True:
            # Friction factor on shell side
            if self.config.tube_arrangement == TubeArrangement.inLine:

                @self.Constraint(
                    self.flowsheet().time, doc="In-line friction factor on shell side"
                )
                def friction_factor_shell_eqn(b, t):
                    return (
                        b.friction_factor_shell[t] * b.N_Re_shell[t] ** 0.15
                        == (
                            0.044
                            + 0.08
                            * b.pitch_x_to_do
                            / (b.pitch_y_to_do - 1.0) ** (0.43 + 1.13 / b.pitch_x_to_do)
                        )
                        * b.fcorrection_dp_shell
                    )

            elif self.config.tube_arrangement == TubeArrangement.staggered:

                @self.Constraint(
                    self.flowsheet().time, doc="Staggered friction factor on shell side"
                )
                def friction_factor_shell_eqn(b, t):
                    return (
                        b.friction_factor_shell[t] * b.N_Re_shell[t] ** 0.16
                        == (0.25 + 0.118 / (b.pitch_y_to_do - 1.0) ** 1.08)
                        * b.fcorrection_dp_shell
                    )

            else:
                raise Exception("tube arrangement type not supported")

            # Pressure drop on shell side
            @self.Constraint(self.flowsheet().time, doc="Pressure change on shell side")
            def deltaP_shell_eqn(b, t):
                return (
                    b.deltaP_shell[t]
                    == -1.4
                    * b.friction_factor_shell[t]
                    * b.tube_nrow
                    * b.hot_side.properties_in[t].dens_mol_phase["Vap"]
                    * sum(
                        b.hot_side.properties_in[t].mw_comp[c]
                        * b.hot_side.properties_in[t].mole_frac_comp[c]
                        for c in b.hot_side.properties_in[t].params.component_list
                    )
                    * b.v_shell[t] ** 2
                )

        # Prandtl number
        @self.Constraint(
            self.flowsheet().time, doc="Prandtl number equation on shell side"
        )
        def N_Pr_shell_eqn(b, t):
            return (
                b.N_Pr_shell[t]
                * b.hot_side.properties_in[t].therm_cond
                * sum(
                    b.hot_side.properties_in[t].mw_comp[c]
                    * b.hot_side.properties_in[t].mole_frac_comp[c]
                    for c in b.hot_side.properties_in[t].params.component_list
                )
                == b.hot_side.properties_in[t].cp_mol
                * b.hot_side.properties_in[t].visc_d
            )

        # Nusselt number, currently assume Re>300
        @self.Constraint(
            self.flowsheet().time, doc="Nusselts number equation on shell side"
        )
        def N_Nu_shell_eqn(b, t):
            return (
                b.N_Nu_shell[t]
                == b.f_arrangement
                * 0.33
                * b.N_Re_shell[t] ** 0.6
                * b.N_Pr_shell[t] ** 0.333333
            )

        # Convective heat transfer coefficient on shell side due to convection
        @self.Constraint(
            self.flowsheet().time,
            doc="Convective heat transfer coefficient equation"
            "on shell side due to convection",
        )
        def hconv_shell_conv_eqn(b, t):
            return (
                b.hconv_shell_conv[t] * b.do_tube / 1000
                == b.N_Nu_shell[t] * b.hot_side.properties_in[t].therm_cond / 1000
            )

        # Total convective heat transfer coefficient on shell side
        @self.Constraint(
            self.flowsheet().time,
            doc="Total convective heat transfer " "coefficient equation on shell side",
        )
        def hconv_shell_total_eqn(b, t):
            if self.config.has_radiation is True:
                return (
                    b.hconv_shell_total[t]
                    == b.hconv_shell_conv[t] + b.hconv_shell_rad[t]
                )
            else:
                return b.hconv_shell_total[t] == b.hconv_shell_conv[t]

        # Wall conduction heat transfer resistance
        # based on outside surface area
        @self.Constraint(doc="Wall conduction heat transfer resistance")
        def rcond_wall_eqn(b):
            return b.rcond_wall * b.therm_cond_wall == 0.5 * b.do_tube * log(
                b.do_tube / b.tube_di
            )

        # Overall heat transfer coefficient
        @self.Constraint(
            self.flowsheet().time, doc="Wall conduction heat transfer resistance"
        )
        def overall_heat_transfer_coefficient_eqn(b, t):
            return (
                b.overall_heat_transfer_coefficient[t]
                * (
                    b.rcond_wall
                    + b.tube_r_fouling
                    + b.shell_r_fouling
                    + 1.0 / b.hconv_shell_total[t]
                    + b.do_tube / b.hconv_tube[t] / b.tube_di
                )
                == b.fcorrection_htc
            )

    def model_check(blk):
        """
        Model checks for unit - calls model checks for both control volume
        Blocks.

        Args:
            None

        Returns:
            None
        """
        # Run control volume block model checks
        blk.cold_side.model_check()
        blk.hot_side.model_check()

    def initialize_build(
        blk,
        state_args_1=None,
        state_args_2=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General Heat Exchanger initialisation routine.

        Keyword Arguments:
            state_args_1 : a dict of arguments to be passed to the property
                           package(s) for side 1 of the heat exchanger to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            state_args_2 : a dict of arguments to be passed to the property
                           package(s) for side 2 of the heat exchanger to
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
        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize inlet property blocks
        flags1 = blk.cold_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_1
        )

        flags2 = blk.hot_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_2
        )
        init_log.info("{} Initialisation Step 1 Complete.".format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize temperature differentials
        p1_flags = {}
        p2_flags = {}
        h1_flags = {}
        t2_flags = {}
        for t in blk.flowsheet().time:
            p1_flags[t] = blk.cold_side.properties_out[t].pressure.fixed
            if (
                not blk.cold_side.properties_out[t].pressure.fixed
                and blk.config.cold_side.has_pressure_change
            ):
                blk.cold_side.properties_out[t].pressure.fix(
                    value(blk.cold_side.properties_in[t].pressure)
                )

            p2_flags[t] = blk.hot_side.properties_out[t].pressure.fixed
            if (
                not blk.hot_side.properties_out[t].pressure.fixed
                and blk.config.hot_side.has_pressure_change
            ):
                blk.hot_side.properties_out[t].pressure.fix(
                    value(blk.hot_side.properties_in[t].pressure)
                )

            h1_flags[t] = blk.cold_side.properties_out[t].enth_mol.fixed
            if not blk.cold_side.properties_out[t].enth_mol.fixed:
                blk.cold_side.properties_out[t].enth_mol.fix(
                    value(blk.cold_side.properties_in[t].enth_mol) + 100.0
                )

            t2_flags[t] = blk.hot_side.properties_out[t].temperature.fixed
            if not blk.hot_side.properties_out[t].temperature.fixed:
                blk.hot_side.properties_out[t].temperature.fix(
                    value(blk.hot_side.properties_in[t].temperature) - 5.0
                )
                #                                assuming Delta T min approach
        # Deactivate Constraints
        blk.heat_transfer_equation.deactivate()
        blk.unit_heat_balance.deactivate()
        if blk.config.cold_side.has_pressure_change:
            blk.deltaP_tube_eqn.deactivate()
        if blk.config.hot_side.has_pressure_change:
            blk.deltaP_shell_eqn.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Activate energy balance and driving force
        for t in blk.flowsheet().time:
            if not p1_flags[t]:
                blk.cold_side.properties_out[t].pressure.unfix()
            if not p2_flags[t]:
                blk.hot_side.properties_out[t].pressure.unfix()
            if not h1_flags[t]:
                blk.cold_side.properties_out[t].enth_mol.unfix()
            if not t2_flags[t]:
                blk.hot_side.properties_out[t].temperature.unfix()
        blk.heat_transfer_equation.activate()
        blk.unit_heat_balance.activate()

        if blk.config.cold_side.has_pressure_change:
            blk.deltaP_tube_eqn.activate()
        if blk.config.hot_side.has_pressure_change:
            blk.deltaP_shell_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.cold_side.release_state(flags1, outlvl)
        blk.hot_side.release_state(flags2, outlvl)

        init_log.info("{} Initialisation Complete.".format(blk.name))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # We have a pretty good idea that the delta Ts will be between about
        # 1 and 100 regardless of process of temperature units, so a default
        # should be fine, so don't warn.  Guessing a typical delta t around 10
        # the default scaling factor is set to 0.1
        sf_dT1 = dict(
            zip(
                self.delta_temperature_in.keys(),
                [
                    iscale.get_scaling_factor(v, default=0.1)
                    for v in self.delta_temperature_in.values()
                ],
            )
        )
        sf_dT2 = dict(
            zip(
                self.delta_temperature_out.keys(),
                [
                    iscale.get_scaling_factor(v, default=0.1)
                    for v in self.delta_temperature_out.values()
                ],
            )
        )

        # U depends a lot on the process and units of measure so user should set
        # this one.
        sf_u = dict(
            zip(
                self.overall_heat_transfer_coefficient.keys(),
                [
                    iscale.get_scaling_factor(v, default=0.01, warning=True)
                    for v in self.overall_heat_transfer_coefficient.values()
                ],
            )
        )

        # Since this depends on the process size this is another scaling factor
        # the user should always set.
        sf_a = iscale.get_scaling_factor(self.area, default=1e-4, warning=True)

        iscale.constraint_scaling_transform(self.area_eqn, sf_a)

        for t, c in self.v_shell_eqn.items():
            s = iscale.min_scaling_factor(
                self.hot_side.properties_in[t].flow_mol_comp,
                default=0,
                warning=False,
                hint=None,
            )
            if s == 0:
                s = iscale.get_scaling_factor(
                    self.hot_side.properties_in[t].flow_mol, default=1, warning=True
                )
            iscale.constraint_scaling_transform(c, s, overwrite=False)

        for t, c in self.v_tube_eqn.items():
            s = iscale.min_scaling_factor(
                self.cold_side.properties_in[t].flow_mol_comp,
                default=0,
                warning=False,
                hint=None,
            )
            if s == 0:
                s = iscale.get_scaling_factor(
                    self.cold_side.properties_in[t].flow_mol, default=1, warning=True
                )
            iscale.constraint_scaling_transform(c, s, overwrite=False)

        for t, c in self.N_Nu_tube_eqn.items():
            s = iscale.get_scaling_factor(self.N_Nu_tube[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)

        for t, c in self.N_Nu_shell_eqn.items():
            s = iscale.get_scaling_factor(self.N_Nu_shell[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, s, overwrite=False)

        for t, c in self.hconv_shell_total_eqn.items():
            s = iscale.get_scaling_factor(
                self.hconv_shell_total[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, s, overwrite=False)

        if hasattr(self, "deltaP_shell_eqn"):
            for t, c in self.deltaP_shell_eqn.items():
                s = iscale.get_scaling_factor(
                    self.deltaP_shell[t], default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, s, overwrite=False)

        if hasattr(self, "deltaP_tube_eqn"):
            for t, c in self.deltaP_tube_eqn.items():
                s = iscale.get_scaling_factor(
                    self.deltaP_tube[t], default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, s, overwrite=False)

        if hasattr(self, "deltaP_tube_friction_eqn"):
            for t, c in self.deltaP_tube_friction_eqn.items():
                s = iscale.get_scaling_factor(
                    self.deltaP_tube_friction[t], default=1, warning=True
                )
                iscale.constraint_scaling_transform(c, s, overwrite=False)
