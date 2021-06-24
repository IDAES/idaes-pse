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
# Import Python libraries
import logging
from enum import Enum

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
# Additional import for the unit operation
from pyomo.environ import SolverFactory, value, Var, Param, exp, sqrt,\
    log, PositiveReals, NonNegativeReals, units as pyunits
from pyomo.opt import TerminationCondition

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)

from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
from idaes.core.util.constants import Constants as c
from idaes.core.util import get_solver

import idaes.logger as idaeslog


__author__ = "Boiler subsystem team (J Ma, M Zamarripa)"
__version__ = "1.0.0"

# Set up logger
_log = logging.getLogger(__name__)


class TubeArrangement(Enum):
    inLine = 0
    staggered = 1


class DeltaTMethod(Enum):
    counterCurrent = 0
    coCurrent = 1


@declare_process_block_class("BoilerHeatExchanger")
class BoilerHeatExchangerData(UnitModelBlockData):
    """
    Standard Heat Exchanger Unit Model Class
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([useDefault, True, False]),
        default=useDefault,
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}"""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))
    CONFIG.declare("side_1_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("side_1_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("side_2_property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("side_2_property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.useDefault,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of material balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}"""))
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.useDefault,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    CONFIG.declare("momentum_balance_type", ConfigValue(
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
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    CONFIG.declare("delta_T_method", ConfigValue(
        default=DeltaTMethod.counterCurrent,
        domain=In(DeltaTMethod),
        description="Flow configuration in unit to compute delta T",
        doc="""Flag indicating type of flow arrangement to use for delta
**default** - DeltaTMethod.counterCurrent
**Valid values:** {
**DeltaTMethod.counterCurrent**}"""))
    CONFIG.declare("tube_arrangement", ConfigValue(
        default=TubeArrangement.inLine,
        domain=In(TubeArrangement),
        description='tube configuration',
        doc='Tube arrangement could be in-line and staggered'))
    CONFIG.declare("side_1_water_phase", ConfigValue(
        default='Liq',
        domain=In(['Liq', 'Vap']),
        description='side 1 water phase',
        doc='Define water phase for property calls'))
    CONFIG.declare("has_radiation", ConfigValue(
        default=False,
        domain=In([False, True]),
        description='Has side 2 gas radiation',
        doc='Define if side 2 gas radiation is to be considered'))

    def build(self):
        """
        Build method for Boiler heat exchanger model

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(BoilerHeatExchangerData, self).build()

        # Build ControlVolume Block
        self.side_1 = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.side_1_property_package,
            "property_package_args": self.config.side_1_property_package_args})

        self.side_2 = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.side_2_property_package,
            "property_package_args": self.config.side_2_property_package_args})

        # Add Geometry
        self.side_1.add_geometry()
        self.side_2.add_geometry()

        # Add state block
        self.side_1.add_state_blocks(has_phase_equilibrium=False)

        # Add material balance
        self.side_1.add_material_balances(
            balance_type=self.config.material_balance_type)
        # add energy balance
        self.side_1.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=True)
        # add momentum balance
        self.side_1.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # Add state block
        self.side_2.add_state_blocks(has_phase_equilibrium=False)

        # Add material balance
        self.side_2.add_material_balances(
            balance_type=self.config.material_balance_type)
        # add energy balance
        self.side_2.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=True)
        # add momentum balance
        self.side_2.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # Set Unit Geometry and control volume
        self._set_geometry()

        self.side_1_fluid_phase = self.config.side_1_water_phase

        # Construct performance equations
        self._make_performance()

        # Construct performance equations
        if self.config.delta_T_method == DeltaTMethod.counterCurrent:
            self._make_counter_current()
        else:
            self._make_co_current()

        self.add_inlet_port(name="side_1_inlet", block=self.side_1)
        self.add_inlet_port(name="side_2_inlet", block=self.side_2)
        self.add_outlet_port(name="side_1_outlet", block=self.side_1)
        self.add_outlet_port(name="side_2_outlet", block=self.side_2)

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
            doc='Elevation increase used for static pressure calculation - m',
            units=pyunits.m)

        # Number of tube columns in the cross section plane
        # perpendicular to shell side fluid flow (y direction)
        self.tube_ncol = Var(initialize=10.0,
                             within=PositiveReals,
                             doc='Number of tube columns')

        # Number of tube rows in the direction of shell side
        # fluid flow (x direction)
        self.tube_nrow = Var(initialize=10.0,
                             within=PositiveReals,
                             doc='Number of tube rows')

        # Number of inlet tube rows
        self.nrow_inlet = Var(initialize=1,
                              within=PositiveReals,
                              doc='Number of inlet tube rows')

        # Length of a tube in z direction for each path
        self.tube_length = Var(initialize=5.0,
                               within=PositiveReals,
                               doc='Tube length - m',
                               units=pyunits.m)

        # Inner diameter of tubes
        self.tube_di = Var(initialize=0.05,
                           within=PositiveReals,
                           doc='Inner diameter of tube - m',
                           units=pyunits.m)

        # Thickness of tube
        self.tube_thickness = Var(initialize=0.005,
                                  within=PositiveReals,
                                  doc='Tube thickness - m',
                                  units=pyunits.m)

        # Pitch of tubes between two neighboring columns (in y direction).
        # Always greater than tube outside diameter
        self.pitch_y = Var(initialize=0.1,
                           within=PositiveReals,
                           doc='Pitch between two neighboring columns - m',
                           units=pyunits.m)

        # Pitch of tubes between two neighboring rows (in x direction).
        # Always greater than tube outside diameter
        self.pitch_x = Var(initialize=0.1,
                           within=PositiveReals,
                           doc='Pitch between two neighboring rows - m',
                           units=pyunits.m)

        # Tube outside diameter
        @self.Expression(doc="Outside diameter of tube - m")
        def do_tube(b):
            return b.tube_di + b.tube_thickness * 2.0

        if self.config.has_radiation is True:
            # Mean beam length for radiation
            @self.Expression(doc="Mean beam length - m")
            def mbl(b):
                return 3.6*(b.pitch_x*b.pitch_y/c.pi/b.do_tube - b.do_tube/4.0)

            # Mean beam length for radiation divided by sqrt(2)
            @self.Expression(doc="Mean beam length - m")
            def mbl_div2(b):
                return b.mbl/sqrt(2.0)

            # Mean beam length for radiation multiplied by sqrt(2)
            @self.Expression(doc="Mean beam length - m")
            def mbl_mul2(b):
                return b.mbl*sqrt(2.0)

        # Number of 180 degree bends for the tube
        @self.Expression(doc="Nbend_tube")
        def nbend_tube(b):
            return b.tube_nrow / b.nrow_inlet

        # Total flow area on tube side
        @self.Expression(doc="Total flow area on tube side - m2")
        def area_flow_tube(b):
            return 0.25 * c.pi * b.tube_di**2.0 * b.tube_ncol * b.nrow_inlet

        # Total flow area on shell side
        @self.Expression(doc="Total flow area on shell side - m2")
        def area_flow_shell(b):
            return b.tube_length * (b.pitch_y - b.do_tube) * b.tube_ncol

        # Total heat transfer area based on outside diameter
        @self.Expression(doc="Total heat transfer "
                         "area based on tube outside diamer - m2")
        def area_heat_transfer(b):
            return c.pi * b.do_tube * b.tube_length * b.tube_ncol * b.tube_nrow

        # Ratio of pitch_x/do_tube
        @self.Expression(doc="Ratio of pitch in x "
                         "direction to tube outside diamer")
        def pitch_x_to_do(b):
            return b.pitch_x / b.do_tube

        # Ratio of pitch_y/do_tube
        @self.Expression(doc="Ratio of pitch in y "
                         "direction to tube outside diamer")
        def pitch_y_to_do(b):
            return b.pitch_y / b.do_tube

        if self.config.has_holdup is True:
            add_object_reference(self, "volume_side_1", self.side_1.volume)
            add_object_reference(self, "volume_side_2", self.side_2.volume)
            # Total tube side valume
            self.Constraint(doc="Total tube side volume")

            def volume_side_1_eqn(b):
                return b.volumne_side_1 == (
                    0.25 * c.pi * b.tube_di**2.0 * b.tube_length
                    * b.tube_ncol * b.tube_nrow)
            # Total shell side valume
            self.Constraint(doc="Total shell side volume")

            def volume_side_2_eqn(b):
                return b.volumne_side_2 == \
                    b.tube_ncol * b.pitch_y * b.tube_length \
                    * b.tube_nrow * b.pitch_x - 0.25 * c.pi * b.do_tube**2.0 \
                    * b.tube_length * b.tube_ncol * b.tube_nrow

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.

        Args:
            None

        Returns:
            None
        """
        # Set references to balance terms at unit level
        add_object_reference(self, "heat_duty", self.side_1.heat)
        if self.config.has_pressure_change is True:
            add_object_reference(self, "deltaP_tube", self.side_1.deltaP)
            add_object_reference(self, "deltaP_shell", self.side_2.deltaP)

        # Performance parameters and variables
        # Wall thermal conductivity
        self.therm_cond_wall = Param(
            initialize=43.0,
            within=PositiveReals,
            doc="Thermal conductivity of the wall - W/(m K)",
            units=pyunits.W/pyunits.m/pyunits.K)

        # Loss coefficient for a 180 degree bend (u-turn),
        # usually related to radius to inside diameter ratio
        self.k_loss_uturn = Param(initialize=0.5,
                                  within=PositiveReals,
                                  mutable=True,
                                  doc='Loss coefficient of a tube u-turn')

        # Heat transfer resistance due to the fouling on tube side
        # (typical boiler hx)
        self.tube_r_fouling = Param(
            initialize=0.00017,
            within=NonNegativeReals,
            mutable=True,
            doc="Fouling resistance on tube side - K m2 / W",
            units=pyunits.K*pyunits.m**2*pyunits.W**-1)

        # Heat transfer resistance due to the fouling on shell side
        self.shell_r_fouling = Param(
            initialize=0.0008,
            within=NonNegativeReals,
            mutable=True,
            doc="Fouling resistance on tube side - K m2 / W",
            units=pyunits.K*pyunits.m**2*pyunits.W**-1)

        # Correction factor for overall heat transfer coefficient
        self.fcorrection_htc = Var(initialize=1.0,
                                   within=NonNegativeReals,
                                   doc="Correction factor for HTC")

        # Correction factor for tube side pressure drop due to friction
        self.fcorrection_dp_tube = Var(
            initialize=1.0,
            doc="Correction factor for tube side pressure drop")

        # Correction factor for shell side pressure drop due to friction
        self.fcorrection_dp_shell = Var(
            initialize=1.0,
            doc="Correction factor for shell side pressure drop")

        # Temperature driving force
        self.temperature_driving_force = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Mean driving force for heat exchange - K",
            units=pyunits.K)

        if self.config.has_radiation is True:
            # Shell side wall emissivity, converted from parameter to variable
            self.emissivity_wall = Var(initialize=0.7,
                                       doc='Shell side wall emissivity')
            # Gas emissivity at mbl
            self.gas_emissivity = Var(
                self.flowsheet().time,
                initialize=0.5,
                doc="Emissivity at given mean beam length")

            # Gas emissivity at mbl/sqrt(2)
            self.gas_emissivity_div2 = Var(
                self.flowsheet().time,
                initialize=0.4,
                doc="Emissivity at mean beam length divided by sqrt of 2")

            # Gas emissivity at mbl*sqrt(2)
            self.gas_emissivity_mul2 = Var(
                self.flowsheet().time,
                initialize=0.6,
                doc="Emissivity at mean beam length multiplied by sqrt of 2")

            # Gray fraction of gas in entire spectrum
            self.gas_gray_fraction = Var(
                self.flowsheet().time,
                initialize=0.5,
                doc="Gray fraction of gas in entire spectrum")

            # Gas-surface radiation exchange factor for shell side wall
            self.frad_gas_shell = Var(self.flowsheet().time,
                                      initialize=0.5,
                                      doc="Gas-surface radiation exchange "
                                      "factor for shell side wall")

            # Shell side equivalent convective heat transfer coefficient
            # due to radiation
            self.hconv_shell_rad = Var(
                self.flowsheet().time,
                initialize=100.0,
                doc="Shell convective heat transfer coefficient due to radiation",
                units=pyunits.W/pyunits.m**2/pyunits.K)

        # Temperature difference at side 1 inlet
        self.deltaT_1 = Var(self.flowsheet().time,
                            initialize=1.0,
                            doc="Temperature difference at side 1 inlet - K",
                            units=pyunits.K)

        # Temperature difference at side 1 outlet
        self.deltaT_2 = Var(self.flowsheet().time,
                            initialize=1.0,
                            doc="Temperature difference at side 1 outlet - K",
                            units=pyunits.K)

        # Overall heat transfer coefficient
        self.overall_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            initialize=1.0,
            units=pyunits.W/pyunits.m**2/pyunits.K)

        # Tube side convective heat transfer coefficient
        self.hconv_tube = Var(
            self.flowsheet().time,
            initialize=100.0,
            doc="Tube side convective heat transfer coefficient - W / (m2 K)",
            units=pyunits.W/pyunits.m**2/pyunits.K)

        # Shell side convective heat transfer coefficient due to convection
        self.hconv_shell_conv = Var(
            self.flowsheet().time,
            initialize=100.0,
            doc="Shell side convective heat transfer coefficient due to convection",
            units=pyunits.W/pyunits.m**2/pyunits.K)

        # Total shell side convective heat transfer coefficient
        # including convection and radiation
        self.hconv_shell_total = Var(
            self.flowsheet().time,
            initialize=150.0,
            doc="Total shell side convective heat transfer coefficient",
            units=pyunits.W/pyunits.m**2/pyunits.K)

        # Heat conduction resistance of tube wall
        self.rcond_wall = Var(
            initialize=1.0,
            doc="Heat conduction resistance of wall - K m2 / W",
            units=pyunits.m**2*pyunits.K/pyunits.W)

        if self.config.has_radiation is True:
            # Constraints for gas emissivity
            @self.Constraint(self.flowsheet().time, doc="Gas emissivity")
            def gas_emissivity_eqn(b, t):
                # This is a surrogate model, so need to do units manually
                X1 = (b.side_2.properties_in[t].temperature
                      + b.side_2.properties_out[t].temperature)/2/pyunits.K
                X2 = b.mbl/pyunits.m
                X3 = b.side_2.properties_in[t].pressure/pyunits.Pa
                X4 = b.side_2.properties_in[t].mole_frac_comp['CO2']
                X5 = b.side_2.properties_in[t].mole_frac_comp['H2O']
                X6 = b.side_2.properties_in[t].mole_frac_comp['O2']

                # Surrogate model fitted using rigorous calc. - 500 samples
                # Wide operating range:
                #      X1: 700 – 1500    (Gas Temperature)
                #      X2: 0.2 – 1       (Mean beam length)
                #      X3: 79000-102000  (pressure in Pa)
                #      X4: 0.12-0.16     (mol frac CO2)
                #      X5: 0.075-0.15    (mol frac H2O)
                #      X6: 0.01-0.07     (mol frac O2)

                return b.gas_emissivity[t] == \
                    (- 0.116916606892E-003 * X1
                     - 0.29111124038936179309056E-001 * X2
                     + 0.50509651230704191577346E-006 * X3
                     + 1.1844222822155641150488 * X4
                     - 0.64720757767102773949652E-001 * X5
                     - 0.35853593221454795048064E-001 * X6
                     + 0.12227919099126832724878 * log(X1)
                     + 0.45102118316418124410738E-001 * log(X2)
                     + 0.33111863480179408447679E-001 * log(X3)
                     + 0.17674928397780117345084E-001 * log(X5)
                     - 0.12541139396423576016226E-001 * exp(X2)
                     - 0.90251708836308952577099 * exp(X4)
                     + 0.32447078857791738538963E-002 * X2**2
                     - 0.31332075610864829615706E-004 * X1*X2
                     - 0.54639645449809960433102E-009 * X1*X3
                     - 0.19721467902854980460033E-003 * X1*X5
                     + 0.45275517692290622763507E-004 * X1*X6
                     + 0.75458754990630776904396E-006 * X2*X3
                     + 0.39691751689931338564765E-001 * X2*X4
                     + 0.73169514231974708273754 * X2*X5
                     - 0.35852614507684822664491E-001 * X2*X6
                     + 0.39743672195685803976177E-005 * X3*X5
                     + 0.58802879141883679897383E-008 * (X1*X2)**2
                     - 1.2994610452829884472692 * (X2*X5)**2)

            # Constraints for gas emissivity at mbl/sqrt(2)
            @self.Constraint(self.flowsheet().time,
                             doc="Gas emissivity at a lower mean beam length")
            def gas_emissivity_div2_eqn(b, t):
                # This is a surrogate model, so need to do units manually
                X1 = (b.side_2.properties_in[t].temperature
                      + b.side_2.properties_out[t].temperature)/2/pyunits.K
                X2 = b.mbl_div2/pyunits.m
                X3 = b.side_2.properties_in[t].pressure/pyunits.Pa
                X4 = b.side_2.properties_in[t].mole_frac_comp['CO2']
                X5 = b.side_2.properties_in[t].mole_frac_comp['H2O']
                X6 = b.side_2.properties_in[t].mole_frac_comp['O2']

                # Surrogate model fitted using rigorous calc. - 500 samples
                # Wide operating range:
                #       X1: 700 – 1500    (Gas Temperature)
                #       X2: 0.2 – 1       (Mean beam length)
                #       X3: 79000-102000  (pressure in Pa)
                #       X4: 0.12-0.16     (mol frac CO2)
                #       X5: 0.075-0.15    (mol frac H2O)
                #       X6: 0.01-0.07     (mol frac O2)
                return b.gas_emissivity_div2[t] == \
                    (- 0.116916606892E-003 * X1
                     - 0.29111124038936179309056E-001 * X2
                     + 0.50509651230704191577346E-006 * X3
                     + 1.1844222822155641150488 * X4
                     - 0.64720757767102773949652E-001 * X5
                     - 0.35853593221454795048064E-001 * X6
                     + 0.12227919099126832724878 * log(X1)
                     + 0.45102118316418124410738E-001 * log(X2)
                     + 0.33111863480179408447679E-001 * log(X3)
                     + 0.17674928397780117345084E-001 * log(X5)
                     - 0.12541139396423576016226E-001 * exp(X2)
                     - 0.90251708836308952577099 * exp(X4)
                     + 0.32447078857791738538963E-002 * X2**2
                     - 0.31332075610864829615706E-004 * X1*X2
                     - 0.54639645449809960433102E-009 * X1*X3
                     - 0.19721467902854980460033E-003 * X1*X5
                     + 0.45275517692290622763507E-004 * X1*X6
                     + 0.75458754990630776904396E-006 * X2*X3
                     + 0.39691751689931338564765E-001 * X2*X4
                     + 0.73169514231974708273754 * X2*X5
                     - 0.35852614507684822664491E-001 * X2*X6
                     + 0.39743672195685803976177E-005 * X3*X5
                     + 0.58802879141883679897383E-008 * (X1*X2)**2
                     - 1.2994610452829884472692 * (X2*X5)**2)

            # Constraints for gas emissivity at mbl*sqrt(2)
            @self.Constraint(self.flowsheet().time,
                             doc="Gas emissivity at a higher mean beam length")
            def gas_emissivity_mul2_eqn(b, t):
                # This is a surrogate model, so need to do units manually
                X1 = (b.side_2.properties_in[t].temperature
                      + b.side_2.properties_out[t].temperature)/2/pyunits.K
                X2 = b.mbl_mul2/pyunits.m
                X3 = b.side_2.properties_in[t].pressure/pyunits.Pa
                X4 = b.side_2.properties_in[t].mole_frac_comp['CO2']
                X5 = b.side_2.properties_in[t].mole_frac_comp['H2O']
                X6 = b.side_2.properties_in[t].mole_frac_comp['O2']

                # Surrogate model fitted using rigorous calc. 500 samples
                # Wide operating range:
                #       X1: 700 – 1500    (Gas Temperature)
                #       X2: 0.2 – 1       (Mean beam length)
                #       X3: 79000-102000  (pressure in Pa)
                #       X4: 0.12-0.16     (mol frac CO2)
                #       X5: 0.075-0.15    (mol frac H2O)
                #       X6: 0.01-0.07     (mol frac O2)
                return b.gas_emissivity_mul2[t] == \
                    (- 0.116916606892E-003 * X1
                     - 0.29111124038936179309056E-001 * X2
                     + 0.50509651230704191577346E-006 * X3
                     + 1.1844222822155641150488 * X4
                     - 0.64720757767102773949652E-001 * X5
                     - 0.35853593221454795048064E-001 * X6
                     + 0.12227919099126832724878 * log(X1)
                     + 0.45102118316418124410738E-001 * log(X2)
                     + 0.33111863480179408447679E-001 * log(X3)
                     + 0.17674928397780117345084E-001 * log(X5)
                     - 0.12541139396423576016226E-001 * exp(X2)
                     - 0.90251708836308952577099 * exp(X4)
                     + 0.32447078857791738538963E-002 * X2**2
                     - 0.31332075610864829615706E-004 * X1*X2
                     - 0.54639645449809960433102E-009 * X1*X3
                     - 0.19721467902854980460033E-003 * X1*X5
                     + 0.45275517692290622763507E-004 * X1*X6
                     + 0.75458754990630776904396E-006 * X2*X3
                     + 0.39691751689931338564765E-001 * X2*X4
                     + 0.73169514231974708273754 * X2*X5
                     - 0.35852614507684822664491E-001 * X2*X6
                     + 0.39743672195685803976177E-005 * X3*X5
                     + 0.58802879141883679897383E-008 * (X1*X2)**2
                     - 1.2994610452829884472692 * (X2*X5)**2)

            # fraction of gray gas spectrum
            @self.Constraint(self.flowsheet().time,
                             doc="Fraction of gray gas spectrum")
            def gas_gray_fraction_eqn(b, t):
                return (b.gas_gray_fraction[t]*(2*b.gas_emissivity_div2[t] -
                        b.gas_emissivity_mul2[t]) ==
                        b.gas_emissivity_div2[t]**2)

            # gas-surface radiation exchange factor
            # between gas and shell side wall
            @self.Constraint(self.flowsheet().time,
                             doc="Gas-surface radiation exchange "
                             "factor between gas and shell side wall")
            def frad_gas_shell_eqn(b, t):
                return (b.frad_gas_shell[t] *
                        ((1/b.emissivity_wall-1)*b.gas_emissivity[t] +
                         b.gas_gray_fraction[t]) ==
                        b.gas_gray_fraction[t]*b.gas_emissivity[t])

            # equivalent convective heat transfer coefficent due to radiation
            @self.Constraint(self.flowsheet().time,
                             doc="Equivalent convective heat transfer "
                             "coefficent due to radiation")
            def hconv_shell_rad_eqn(b, t):
                return b.hconv_shell_rad[t] == \
                    c.stefan_constant * b.frad_gas_shell[t] * \
                    ((b.side_2.properties_in[t].temperature +
                      b.side_2.properties_out[t].temperature)/2
                     + b.side_1.properties_in[t].temperature) * \
                    (((b.side_2.properties_in[t].temperature
                       + b.side_2.properties_out[t].temperature)/2)**2 +
                     b.side_1.properties_in[t].temperature**2)

        # Energy balance equation
        @self.Constraint(self.flowsheet().time,
                         doc="Energy balance between two sides")
        def energy_balance(b, t):
            return b.side_1.heat[t] / 1e6 == -b.side_2.heat[t] / 1e6

        # Heat transfer correlation
        @self.Constraint(self.flowsheet().time,
                         doc="Heat transfer correlation")
        def heat_transfer_correlation(b, t):
            return b.heat_duty[t] / 1e6 == \
                (b.overall_heat_transfer_coefficient[t] *
                 b.area_heat_transfer *
                 b.temperature_driving_force[t]) / 1e6

        # Driving force
        @self.Constraint(self.flowsheet().time,
                         doc="Simplified Log mean temperature "
                         "difference calculation")
        def LMTD(b, t):
            return b.temperature_driving_force[t] == \
                ((b.deltaT_1[t]**0.3241 +
                  b.deltaT_2[t]**0.3241)/1.99996)**(1/0.3241)

        # Tube side heat transfer coefficient and pressure drop
        # -----------------------------------------------------
        # Velocity on tube side
        self.v_tube = Var(self.flowsheet().time,
                          initialize=1.0,
                          doc="Velocity on tube side - m/s",
                          units=pyunits.m/pyunits.s)

        # Reynalds number on tube side
        self.N_Re_tube = Var(self.flowsheet().time,
                             initialize=10000.0,
                             doc="Reynolds number on tube side")
        if self.config.has_pressure_change is True:
            # Friction factor on tube side
            self.friction_factor_tube = Var(self.flowsheet().time,
                                            initialize=1.0,
                                            doc='Friction factor on tube side')

            # Pressure drop due to friction on tube side
            self.deltaP_tube_friction = Var(
                self.flowsheet().time,
                initialize=-10.0,
                doc="Pressure drop due to friction on tube side - Pa",
                units=pyunits.Pa)

            # Pressure drop due to 180 degree turn on tube side
            self.deltaP_tube_uturn = Var(
                self.flowsheet().time,
                initialize=-10.0,
                doc="Pressure drop due to u-turn on tube side - Pa",
                units=pyunits.Pa)

        # Prandtl number on tube side
        self.N_Pr_tube = Var(self.flowsheet().time, initialize=1,
                             doc="Prandtl number on tube side")

        # Nusselt number on tube side
        self.N_Nu_tube = Var(self.flowsheet().time, initialize=1,
                             doc="Nusselts number on tube side")

        # Velocity equation
        @self.Constraint(self.flowsheet().time,
                         doc="Tube side velocity equation - m/s")
        def v_tube_eqn(b, t):
            return (b.v_tube[t] * b.area_flow_tube *
                    b.side_1.properties_in[t].dens_mol_phase[
                        self.side_1_fluid_phase] ==
                    b.side_1.properties_in[t].flow_mol)

        # Reynolds number
        @self.Constraint(self.flowsheet().time,
                         doc="Reynolds number equation on tube side")
        def N_Re_tube_eqn(b, t):
            return (b.N_Re_tube[t] *
                    b.side_1.properties_in[t].visc_d_phase[
                        self.side_1_fluid_phase] ==
                    b.tube_di * b.v_tube[t] *
                    b.side_1.properties_in[t].dens_mass_phase[
                        self.side_1_fluid_phase])

        if self.config.has_pressure_change is True:
            # Friction factor
            @self.Constraint(self.flowsheet().time,
                             doc="Darcy friction factor on tube side")
            def friction_factor_tube_eqn(b, t):
                return b.friction_factor_tube[t]*b.N_Re_tube[t]**0.25 == \
                    0.3164*b.fcorrection_dp_tube

            # Pressure drop due to friction
            @self.Constraint(self.flowsheet().time,
                             doc="Pressure drop due to friction on tube side")
            def deltaP_tube_friction_eqn(b, t):
                return (b.deltaP_tube_friction[t]*b.tube_di*b.nrow_inlet ==
                        -0.5 * b.side_1.properties_in[t].dens_mass_phase[
                            self.side_1_fluid_phase] *
                        b.v_tube[t]**2 * b.friction_factor_tube[t] *
                        b.tube_length * b.tube_nrow)

            # Pressure drop due to u-turn
            @self.Constraint(self.flowsheet().time,
                             doc="Pressure drop due to u-turn on tube side")
            def deltaP_tube_uturn_eqn(b, t):
                return (b.deltaP_tube_uturn[t] ==
                        -0.5 * b.side_1.properties_in[t].dens_mass_phase[
                            self.side_1_fluid_phase] *
                        b.v_tube[t]**2 * b.k_loss_uturn)

            # Total pressure drop on tube side
            @self.Constraint(self.flowsheet().time,
                             doc="Total pressure drop on tube side")
            def deltaP_tube_eqn(b, t):
                return (b.deltaP_tube[t] ==
                        b.deltaP_tube_friction[t] + b.deltaP_tube_uturn[t] -
                        b.delta_elevation * c.acceleration_gravity *
                        (b.side_1.properties_in[t].dens_mass_phase[
                            self.side_1_fluid_phase] +
                         b.side_1.properties_out[t].dens_mass_phase[
                            self.side_1_fluid_phase]) / 2.0)

        # Prandtl number
        @self.Constraint(self.flowsheet().time,
                         doc="Prandtl number equation on tube side")
        def N_Pr_tube_eqn(b, t):
            return (b.N_Pr_tube[t] *
                    b.side_1.properties_in[t].therm_cond_phase[
                        self.side_1_fluid_phase] *
                    b.side_1.properties_in[t].mw ==
                    b.side_1.properties_in[t].cp_mol_phase[
                        self.side_1_fluid_phase] *
                    b.side_1.properties_in[t].visc_d_phase[
                        self.side_1_fluid_phase])

        # Nusselts number
        @self.Constraint(self.flowsheet().time,
                         doc="Nusselts number equation on tube side")
        def N_Nu_tube_eqn(b, t):
            return b.N_Nu_tube[t] == \
                0.023 * b.N_Re_tube[t]**0.8 * b.N_Pr_tube[t]**0.4

        # Heat transfer coefficient
        @self.Constraint(self.flowsheet().time,
                         doc="Convective heat transfer "
                         "coefficient equation on tube side")
        def hconv_tube_eqn(b, t):
            return (b.hconv_tube[t]*self.tube_di/1000 ==
                    b.N_Nu_tube[t] *
                    b.side_1.properties_in[t].therm_cond_phase[
                        self.side_1_fluid_phase]/1000)

        # Pressure drop and heat transfer coefficient on shell side
        # ----------------------------------------------------------
        # Tube arrangement factor
        if self.config.tube_arrangement == TubeArrangement.inLine:
            self.f_arrangement = Param(initialize=0.788,
                                       doc="In-line tube arrangement factor")
        elif self.config.tube_arrangement == TubeArrangement.staggered:
            self.f_arrangement = Param(initialize=1.0,
                                       doc="Staggered tube arrangement factor")
        else:
            raise Exception('tube arrangement type not supported')
        # Velocity on shell side
        self.v_shell = Var(self.flowsheet().time,
                           initialize=1.0,
                           doc="Velocity on shell side - m/s",
                           units=pyunits.m/pyunits.s)

        # Reynalds number on shell side
        self.N_Re_shell = Var(self.flowsheet().time,
                              initialize=10000.0,
                              doc="Reynolds number on shell side")

        # Friction factor on shell side
        self.friction_factor_shell = Var(self.flowsheet().time,
                                         initialize=1.0,
                                         doc='Friction factor on shell side')

        # Prandtl number on shell side
        self.N_Pr_shell = Var(self.flowsheet().time,
                              initialize=1,
                              doc="Prandtl number on shell side")

        # Nusselt number on shell side
        self.N_Nu_shell = Var(self.flowsheet().time,
                              initialize=1,
                              doc="Nusselts number on shell side")

        # Velocity equation on shell side
        @self.Constraint(self.flowsheet().time, doc="Velocity on shell side")
        def v_shell_eqn(b, t):
            return b.v_shell[t] * \
                b.side_2.properties_in[t].dens_mol_phase["Vap"] * \
                b.area_flow_shell == \
                sum(b.side_2.properties_in[t].flow_mol_comp[j]
                    for j in b.side_2.properties_in[t].params.component_list)

        # Reynolds number
        @self.Constraint(self.flowsheet().time,
                         doc="Reynolds number equation on shell side")
        def N_Re_shell_eqn(b, t):
            return b.N_Re_shell[t] * b.side_2.properties_in[t].visc_d == \
                b.do_tube * b.v_shell[t] \
                * b.side_2.properties_in[t].dens_mol_phase["Vap"] *\
                   sum(b.side_2.properties_in[t].mw_comp[c]
                       * b.side_2.properties_in[t].mole_frac_comp[c]
                       for c in b.side_2.properties_in[t].
                       params.component_list)

        if self.config.has_pressure_change is True:
            # Friction factor on shell side
            if self.config.tube_arrangement == TubeArrangement.inLine:
                @self.Constraint(self.flowsheet().time,
                                 doc="In-line friction factor on shell side")
                def friction_factor_shell_eqn(b, t):
                    return b.friction_factor_shell[t] \
                        * b.N_Re_shell[t]**0.15 == \
                        (0.044 + 0.08 * b.pitch_x_to_do
                         / (b.pitch_y_to_do - 1.0)**(0.43 + 1.13
                                                     / b.pitch_x_to_do)
                         ) * b.fcorrection_dp_shell

            elif self.config.tube_arrangement == TubeArrangement.staggered:
                @self.Constraint(self.flowsheet().time,
                                 doc="Staggered friction factor on shell side")
                def friction_factor_shell_eqn(b, t):
                    return b.friction_factor_shell[t] \
                        * b.N_Re_shell[t]**0.16 == \
                        (0.25 + 0.118 / (b.pitch_y_to_do - 1.0)**1.08) \
                        * b.fcorrection_dp_shell
            else:
                raise Exception('tube arrangement type not supported')

            # Pressure drop on shell side
            @self.Constraint(self.flowsheet().time,
                             doc="Pressure change on shell side")
            def deltaP_shell_eqn(b, t):
                return (
                    b.deltaP_shell[t] ==
                    -1.4 * b.friction_factor_shell[t] * b.tube_nrow *
                    b.side_2.properties_in[t].dens_mol_phase["Vap"] *
                    sum(b.side_2.properties_in[t].mw_comp[c] *
                        b.side_2.properties_in[t].mole_frac_comp[c] for c
                        in b.side_2.properties_in[t].params.component_list) *
                    b.v_shell[t]**2)

        # Prandtl number
        @self.Constraint(self.flowsheet().time,
                         doc="Prandtl number equation on shell side")
        def N_Pr_shell_eqn(b, t):
            return b.N_Pr_shell[t] * b.side_2.properties_in[t].therm_cond \
                * sum(b.side_2.properties_in[t].mw_comp[c]
                      * b.side_2.properties_in[t].mole_frac_comp[c]
                      for c in b.side_2.properties_in[t].
                      params.component_list) == \
                b.side_2.properties_in[t].cp_mol * \
                b.side_2.properties_in[t].visc_d

        # Nusselt number, currently assume Re>300
        @self.Constraint(self.flowsheet().time,
                         doc="Nusselts number equation on shell side")
        def N_Nu_shell_eqn(b, t):
            return b.N_Nu_shell[t] == b.f_arrangement * 0.33 \
                * b.N_Re_shell[t]**0.6 * b.N_Pr_shell[t]**0.333333

        # Convective heat transfer coefficient on shell side due to convection
        @self.Constraint(self.flowsheet().time,
                         doc="Convective heat transfer coefficient equation"
                         "on shell side due to convection")
        def hconv_shell_conv_eqn(b, t):
            return b.hconv_shell_conv[t] * b.do_tube / 1000 == \
                b.N_Nu_shell[t] * b.side_2.properties_in[t].therm_cond\
                / 1000

        # Total convective heat transfer coefficient on shell side
        @self.Constraint(self.flowsheet().time,
                         doc="Total convective heat transfer "
                         "coefficient equation on shell side")
        def hconv_shell_total_eqn(b, t):
            if self.config.has_radiation is True:
                return b.hconv_shell_total[t] == \
                    b.hconv_shell_conv[t] + b.hconv_shell_rad[t]
            else:
                return b.hconv_shell_total[t] == b.hconv_shell_conv[t]

        # Wall conduction heat transfer resistance
        # based on outside surface area
        @self.Constraint(doc="Wall conduction heat transfer resistance")
        def rcond_wall_eqn(b):
            return b.rcond_wall * b.therm_cond_wall == \
                0.5 * b.do_tube * log(b.do_tube / b.tube_di)

        # Overall heat transfer coefficient
        @self.Constraint(self.flowsheet().time,
                         doc="Wall conduction heat transfer resistance")
        def overall_heat_transfer_coefficient_eqn(b, t):
            return b.overall_heat_transfer_coefficient[t] \
                * (b.rcond_wall + b.tube_r_fouling + b.shell_r_fouling +
                   1.0 / b.hconv_shell_total[t]
                   + b.do_tube / b.hconv_tube[t] / b.tube_di) == \
                b.fcorrection_htc

    def _make_co_current(self):
        """
        Add temperature driving force Constraints for co-current flow.

        Args:
            None

        Returns:
            None
        """
        # Temperature Differences
        @self.Constraint(self.flowsheet().time,
                         doc="Side 1 inlet temperature difference")
        def temperature_difference_1(b, t):
            return b.deltaT_1[t] == (
                       b.side_2.properties_in[t].temperature -
                       b.side_1.properties_in[t].temperature)

        @self.Constraint(self.flowsheet().time,
                         doc="Side 1 outlet temperature difference")
        def temperature_difference_2(b, t):
            return b.deltaT_2[t] == (
                       b.side_2.properties_out[t].temperature -
                       b.side_1.properties_out[t].temperature)

    def _make_counter_current(self):
        """
        Add temperature driving force Constraints for counter-current flow.

        Args:
            None

        Returns:
            None
        """
        # Temperature Differences
        @self.Constraint(self.flowsheet().time,
                         doc="Side 1 inlet temperature difference")
        def temperature_difference_1(b, t):
            return b.deltaT_1[t] == (
                       b.side_2.properties_out[t].temperature -
                       b.side_1.properties_in[t].temperature)

        @self.Constraint(self.flowsheet().time,
                         doc="Side 1 outlet temperature difference")
        def temperature_difference_2(b, t):
            return b.deltaT_2[t] == (
                       b.side_2.properties_in[t].temperature -
                       b.side_1.properties_out[t].temperature)

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
        blk.side_1.model_check()
        blk.side_2.model_check()

    def initialize(blk, state_args_1=None, state_args_2=None,
                   outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        '''
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
        '''
        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize inlet property blocks
        flags1 = blk.side_1.initialize(outlvl=outlvl,
                                       optarg=optarg,
                                       solver=solver,
                                       state_args=state_args_1)

        flags2 = blk.side_2.initialize(outlvl=outlvl,
                                       optarg=optarg,
                                       solver=solver,
                                       state_args=state_args_2)
        init_log.info('{} Initialisation Step 1 Complete.'.format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize temperature differentials
        p1_flags = {}
        p2_flags = {}
        h1_flags = {}
        t2_flags = {}
        for t in blk.flowsheet().time:
            p1_flags[t] = blk.side_1.properties_out[t].pressure.fixed
            if not blk.side_1.properties_out[t].pressure.fixed \
                    and blk.config.has_pressure_change:
                blk.side_1.properties_out[t].pressure.fix(
                        value(blk.side_1.properties_in[t].pressure))

            p2_flags[t] = blk.side_2.properties_out[t].pressure.fixed
            if not blk.side_2.properties_out[t].pressure.fixed \
                    and blk.config.has_pressure_change:
                blk.side_2.properties_out[t].pressure.fix(
                        value(blk.side_2.properties_in[t].pressure))

            h1_flags[t] = blk.side_1.properties_out[t].enth_mol.fixed
            if not blk.side_1.properties_out[t].enth_mol.fixed:
                blk.side_1.properties_out[t].enth_mol.fix(
                        value(blk.side_1.properties_in[t].enth_mol)+100.0)

            t2_flags[t] = blk.side_2.properties_out[t].temperature.fixed
            if not blk.side_2.properties_out[t].temperature.fixed:
                blk.side_2.properties_out[t].temperature.fix(
                        value(blk.side_2.properties_in[t].temperature)-5.0)
                #                                assuming Delta T min approach
        # Deactivate Constraints
        blk.heat_transfer_correlation.deactivate()
        blk.LMTD.deactivate()
        blk.energy_balance.deactivate()
        if blk.config.has_pressure_change:
            blk.deltaP_tube_eqn.deactivate()
            blk.deltaP_shell_eqn.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Activate energy balance and driving force
        for t in blk.flowsheet().time:
            if not p1_flags[t]:
                blk.side_1.properties_out[t].pressure.unfix()
            if not p2_flags[t]:
                blk.side_2.properties_out[t].pressure.unfix()
            if not h1_flags[t]:
                blk.side_1.properties_out[t].enth_mol.unfix()
            if not t2_flags[t]:
                blk.side_2.properties_out[t].temperature.unfix()
        blk.heat_transfer_correlation.activate()
        blk.LMTD.activate()
        blk.energy_balance.activate()

        if blk.config.has_pressure_change:
            blk.deltaP_tube_eqn.activate()
            blk.deltaP_shell_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.side_1.release_state(flags1, outlvl)
        blk.side_2.release_state(flags2, outlvl)

        init_log.info('{} Initialisation Complete.'.format(blk.name))
