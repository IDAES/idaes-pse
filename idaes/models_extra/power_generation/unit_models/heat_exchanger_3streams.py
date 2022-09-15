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
3 stream IDAES heat exchanger model with given UA.
side 1 is hot stream, side 2 and 3 are cold streams
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
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


# Additional import for the unit operation
from pyomo.environ import Var, Reference

__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"
__version__ = "1.0.0"


@declare_process_block_class("HeatExchangerWith3Streams")
class HeatExchangerWith3StreamsData(UnitModelBlockData):
    """
    Standard Heat Exchanger Unit Model Class
    """

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "side_1_property_package",
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
        "side_1_property_package_args",
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
        "side_2_property_package",
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
        "side_2_property_package_args",
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
        "side_3_property_package",
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
        "side_3_property_package_args",
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
        "flow_type_side_2",
        ConfigValue(
            default="counter-current",
            domain=In(["counter-current", "co-current"]),
            description="Flow configuration in unit",
            doc="""Flag indicating type of flow arrangement to use for heat
exchanger, **default** 'counter-current' counter-current flow arrangement""",
        ),
    )
    CONFIG.declare(
        "flow_type_side_3",
        ConfigValue(
            default="counter-current",
            domain=In(["counter-current", "co-current"]),
            description="Flow configuration in unit",
            doc="""Flag indicating type of flow arrangement to use for heat
exchanger (default = 'counter-current' - counter-current flow arrangement""",
        ),
    )

    def build(self):
        """
        Begin building model
        """
        # Call UnitModel.build to setup dynamics
        super(HeatExchangerWith3StreamsData, self).build()

        # Build Holdup Block
        self.side_1 = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.side_1_property_package,
            property_package_args=self.config.side_1_property_package_args,
        )

        self.side_2 = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.side_2_property_package,
            property_package_args=self.config.side_2_property_package_args,
        )

        self.side_3 = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.side_3_property_package,
            property_package_args=self.config.side_3_property_package_args,
        )

        # Add Geometry
        self.side_1.add_geometry()
        self.side_2.add_geometry()
        self.side_3.add_geometry()

        # Add state block
        self.side_1.add_state_blocks(has_phase_equilibrium=False)

        # Add material balance
        self.side_1.add_material_balances(
            balance_type=self.config.material_balance_type
        )
        # add energy balance
        self.side_1.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer,
        )
        # add momentum balance
        self.side_1.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Add state block
        self.side_2.add_state_blocks(has_phase_equilibrium=False)

        # Add material balance
        self.side_2.add_material_balances(
            balance_type=self.config.material_balance_type
        )
        # add energy balance
        self.side_2.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer,
        )
        # add momentum balance
        self.side_2.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Add state block
        self.side_3.add_state_blocks(has_phase_equilibrium=False)

        # Add material balance
        self.side_3.add_material_balances(
            balance_type=self.config.material_balance_type
        )
        # add energy balance
        self.side_3.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer,
        )
        # add momentum balance
        self.side_3.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        self._set_geometry()

        # Construct performance equations
        self._make_performance()

        # Construct performance equations
        if self.config.flow_type_side_2 == "counter-current":
            self._make_counter_current_side_2()
        else:
            self._make_co_current_side_2()

        # Construct performance equations
        if self.config.flow_type_side_3 == "counter-current":
            self._make_counter_current_side_3()
        else:
            self._make_co_current_side_3()

        self.add_inlet_port(name="side_1_inlet", block=self.side_1)
        self.add_inlet_port(name="side_2_inlet", block=self.side_2)
        self.add_inlet_port(name="side_3_inlet", block=self.side_3)
        self.add_outlet_port(name="side_1_outlet", block=self.side_1)
        self.add_outlet_port(name="side_2_outlet", block=self.side_2)
        self.add_outlet_port(name="side_3_outlet", block=self.side_3)

    def _set_geometry(self):
        """
        Define the geometry of the unit as necessary, and link to holdup volume

        """

        # UA (product of overall heat transfer coefficient and area)
        # between side 1 and side 2
        self.ua_side_2 = Var(
            self.flowsheet().time, initialize=10.0, doc="UA between side 1 and side 2"
        )

        # UA (product of overall heat transfer coefficient and area)
        # between side 1 and side 3
        self.ua_side_3 = Var(
            self.flowsheet().time, initialize=10.0, doc="UA between side 1 and side 3"
        )

        # fraction of heat from hot stream as heat loss to ambient
        self.frac_heatloss = Var(
            initialize=0.05, doc="Fraction of heat loss to ambient"
        )

        if self.config.has_holdup is True:
            self.volume_side_1 = Reference(self.side_1.volume)
            self.volume_side_2 = Reference(self.side_2.volume)
            self.volume_side_3 = Reference(self.side_3.volume)

    def _make_performance(self):
        """
        Define constraints which describe the behaviour of the unit model.

        Args:
            None

        Returns:
            None
        """
        # Set references to balance terms at unit level
        self.heat_duty_side_1 = Reference(self.side_1.heat)
        self.heat_duty_side_2 = Reference(self.side_2.heat)
        self.heat_duty_side_3 = Reference(self.side_3.heat)

        if self.config.has_pressure_change is True:
            self.deltaP_side_1 = Reference(self.side_1.deltaP)
            self.deltaP_side_2 = Reference(self.side_2.deltaP)
            self.deltaP_side_3 = Reference(self.side_3.deltaP)

        # Performance parameters and variables
        # Temperature driving force
        self.temperature_driving_force_side_2 = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Mean driving force " "for heat exchange",
        )

        # Temperature driving force
        self.temperature_driving_force_side_3 = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Mean driving force " "for heat exchange",
        )

        # Temperature difference at side 2 inlet
        self.side_2_inlet_dT = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Temperature difference " "at side 2 inlet",
        )

        # Temperature difference at side 2 outlet
        self.side_2_outlet_dT = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Temperature difference " "at side 2 outlet",
        )

        # Temperature difference at side 3 inlet
        self.side_3_inlet_dT = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Temperature difference" " at side 3 inlet",
        )

        # Temperature difference at side 3 outlet
        self.side_3_outlet_dT = Var(
            self.flowsheet().time,
            initialize=1.0,
            doc="Temperature difference " "at side 3 outlet",
        )

        # Driving force side 2 (Underwood approximation)
        @self.Constraint(
            self.flowsheet().time,
            doc="Log mean temperature difference calculation "
            "using Underwood approximation",
        )
        def LMTD_side_2(b, t):
            return b.temperature_driving_force_side_2[t] == (
                (b.side_2_inlet_dT[t] ** (1 / 3) + b.side_2_outlet_dT[t] ** (1 / 3)) / 2
            ) ** (3)

        # Driving force side 3 (Underwood approximation)
        @self.Constraint(
            self.flowsheet().time,
            doc="Log mean temperature difference calculation "
            "using Underwood approximation",
        )
        def LMTD_side_3(b, t):
            return b.temperature_driving_force_side_3[t] == (
                (b.side_3_inlet_dT[t] ** (1 / 3) + b.side_3_outlet_dT[t] ** (1 / 3)) / 2
            ) ** (3)

        # Heat duty side 2
        @self.Constraint(self.flowsheet().time, doc="Heat transfer rate")
        def heat_duty_side_2_eqn(b, t):
            return b.heat_duty_side_2[t] == (
                b.ua_side_2[t] * b.temperature_driving_force_side_2[t]
            )

        # Heat duty side 3
        @self.Constraint(self.flowsheet().time, doc="Heat transfer rate")
        def heat_duty_side_3_eqn(b, t):
            return b.heat_duty_side_3[t] == (
                b.ua_side_3[t] * b.temperature_driving_force_side_3[t]
            )

        # Energy balance equation
        @self.Constraint(self.flowsheet().time, doc="Energy balance between two sides")
        def heat_duty_side_1_eqn(b, t):
            return -b.heat_duty_side_1[t] * (1 - b.frac_heatloss) == (
                b.heat_duty_side_2[t] + b.heat_duty_side_3[t]
            )

    def _make_co_current_side_2(self):
        """
        Add temperature driving force Constraints for co-current flow.


        """
        # Temperature Differences
        @self.Constraint(
            self.flowsheet().time, doc="Side 2 inlet temperature difference"
        )
        def side_2_inlet_dT_eqn(b, t):
            return b.side_2_inlet_dT[t] == (
                b.side_1.properties_in[t].temperature
                - b.side_2.properties_in[t].temperature
            )

        @self.Constraint(
            self.flowsheet().time, doc="Side 2 outlet temperature difference"
        )
        def side_2_outlet_dT_eqn(b, t):
            return b.side_2_outlet_dT[t] == (
                b.side_1.properties_out[t].temperature
                - b.side_2.properties_out[t].temperature
            )

    def _make_counter_current_side_2(self):
        """
        Add temperature driving force Constraints for counter-current flow.
        """
        # Temperature Differences
        @self.Constraint(
            self.flowsheet().time, doc="Side 2 inlet temperature difference"
        )
        def side_2_inlet_dT_eqn(b, t):
            return b.side_2_inlet_dT[t] == (
                b.side_1.properties_out[t].temperature
                - b.side_2.properties_in[t].temperature
            )

        @self.Constraint(
            self.flowsheet().time, doc="Side 2 outlet temperature difference"
        )
        def side_2_outlet_dT_eqn(b, t):
            return b.side_2_outlet_dT[t] == (
                b.side_1.properties_in[t].temperature
                - b.side_2.properties_out[t].temperature
            )

    def _make_co_current_side_3(self):
        """
        Add temperature driving force Constraints for co-current flow.
        """
        # Temperature Differences
        @self.Constraint(
            self.flowsheet().time, doc="Side 3 inlet temperature difference"
        )
        def side_3_inlet_dT_eqn(b, t):
            return b.side_3_inlet_dT[t] == (
                b.side_1.properties_in[t].temperature
                - b.side_3.properties_in[t].temperature
            )

        @self.Constraint(
            self.flowsheet().time, doc="Side 3 outlet temperature difference"
        )
        def side_3_outlet_dT_eqn(b, t):
            return b.side_3_outlet_dT[t] == (
                b.side_1.properties_out[t].temperature
                - b.side_3.properties_out[t].temperature
            )

    def _make_counter_current_side_3(self):
        """
        Add temperature driving force Constraints for counter-current flow.
        """
        # Temperature Differences
        @self.Constraint(
            self.flowsheet().time, doc="Side 3 inlet temperature difference"
        )
        def side_3_inlet_dT_eqn(b, t):
            return b.side_3_inlet_dT[t] == (
                b.side_1.properties_out[t].temperature
                - b.side_3.properties_in[t].temperature
            )

        @self.Constraint(
            self.flowsheet().time, doc="Side 3 outlet temperature difference"
        )
        def side_3_outlet_dT_eqn(b, t):
            return b.side_3_outlet_dT[t] == (
                b.side_1.properties_in[t].temperature
                - b.side_3.properties_out[t].temperature
            )

    def initialize_build(
        blk,
        state_args_1=None,
        state_args_2=None,
        state_args_3=None,
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
            state_args_3 : a dict of arguments to be passed to the property
                           package(s) for side 3 of the heat exchanger to
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

        # ---------------------------------------------------------------------
        # Initialize inlet property blocks
        flags1 = blk.side_1.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_1
        )

        flags2 = blk.side_2.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_2
        )

        flags3 = blk.side_3.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args_3
        )

        init_log.info("Initialisation Step 1 Complete.")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info(
            "Initialization Step 2 Complete: {}".format(idaeslog.condition(res))
        )
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.side_1.release_state(flags1, outlvl)
        blk.side_2.release_state(flags2, outlvl)
        blk.side_3.release_state(flags3, outlvl)

        init_log.info_low("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        for t, c in self.heat_duty_side_1_eqn.items():
            sf = iscale.get_scaling_factor(
                self.heat_duty_side_1[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for t, c in self.heat_duty_side_2_eqn.items():
            sf = iscale.get_scaling_factor(
                self.heat_duty_side_2[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, sf, overwrite=False)

        for t, c in self.heat_duty_side_3_eqn.items():
            sf = iscale.get_scaling_factor(
                self.heat_duty_side_3[t], default=1, warning=True
            )
            iscale.constraint_scaling_transform(c, sf, overwrite=False)
