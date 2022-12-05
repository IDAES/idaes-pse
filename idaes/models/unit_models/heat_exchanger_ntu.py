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
IDAES heat exchanger model using effectiveness-NTU method

Assumptions:
    * No phase equilibrium or reactions occur within unit

"""

# Import Pyomo libraries
from pyomo.environ import (
    check_optimal_termination,
    Constraint,
    Expression,
    Param,
    PositiveReals,
    Reference,
    units as pyunits,
    Var,
)
from pyomo.common.config import Bool, ConfigBlock, ConfigValue, In

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
from idaes.models.unit_models.heat_exchanger import hx_process_config, add_hx_references
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.math import smooth_min, smooth_max
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

__author__ = "Paul Akula, Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("HeatExchangerNTU")
class HeatExchangerNTUData(UnitModelBlockData):
    """Heat Exchanger Unit Model using NTU method."""

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)

    # Configuration template for fluid specific  arguments
    _SideCONFIG = ConfigBlock()

    _SideCONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    _SideCONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    _SideCONFIG.declare(
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
    _SideCONFIG.declare(
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
    _SideCONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use ",
            doc="""Property parameter object used to define property calculations
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    _SideCONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property package",
            doc="""A ConfigBlock with arguments to be passed to
        property block(s) and used when constructing these,
        **default** - None.
        **Valid values:** {
        see property package for documentation.}""",
        ),
    )

    # Create individual config blocks for hot and cold sides
    CONFIG.declare("hot_side", _SideCONFIG(doc="Hot fluid config arguments"))
    CONFIG.declare("cold_side", _SideCONFIG(doc="Cold fluid config arguments"))
    CONFIG.declare(
        "hot_side_name",
        ConfigValue(
            default=None,
            domain=str,
            doc="Hot side name, sets control volume and inlet and outlet names",
        ),
    )
    CONFIG.declare(
        "cold_side_name",
        ConfigValue(
            default=None,
            domain=str,
            doc="Cold side name, sets control volume and inlet and outlet names",
        ),
    )

    def build(self):
        # Call UnitModel.build to setup model
        super().build()
        hx_process_config(self)

        # ---------------------------------------------------------------------
        # Build hot-side control volume
        self.hot_side = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.hot_side.property_package,
            property_package_args=self.config.hot_side.property_package_args,
        )

        # TODO : Add support for phase equilibrium?
        self.hot_side.add_state_blocks(has_phase_equilibrium=False)

        self.hot_side.add_material_balances(
            balance_type=self.config.hot_side.material_balance_type,
            has_phase_equilibrium=False,
        )

        self.hot_side.add_energy_balances(
            balance_type=self.config.hot_side.energy_balance_type,
            has_heat_transfer=True,
        )

        self.hot_side.add_momentum_balances(
            balance_type=self.config.hot_side.momentum_balance_type,
            has_pressure_change=self.config.hot_side.has_pressure_change,
        )

        # ---------------------------------------------------------------------
        # Build cold-side control volume
        self.cold_side = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.cold_side.property_package,
            property_package_args=self.config.cold_side.property_package_args,
        )

        self.cold_side.add_state_blocks(has_phase_equilibrium=False)

        self.cold_side.add_material_balances(
            balance_type=self.config.cold_side.material_balance_type,
            has_phase_equilibrium=False,
        )

        self.cold_side.add_energy_balances(
            balance_type=self.config.cold_side.energy_balance_type,
            has_heat_transfer=True,
        )

        self.cold_side.add_momentum_balances(
            balance_type=self.config.cold_side.momentum_balance_type,
            has_pressure_change=self.config.cold_side.has_pressure_change,
        )

        # ---------------------------------------------------------------------
        # Add Ports to control volumes
        self.add_inlet_port(
            name="hot_side_inlet", block=self.hot_side, doc="Hot side inlet port"
        )
        self.add_outlet_port(
            name="hot_side_outlet", block=self.hot_side, doc="Hot side outlet port"
        )

        self.add_inlet_port(
            name="cold_side_inlet", block=self.cold_side, doc="Cold side inlet port"
        )
        self.add_outlet_port(
            name="cold_side_outlet", block=self.cold_side, doc="Cold side outlet port"
        )

        # ---------------------------------------------------------------------
        # Add unit level References
        # Set references to balance terms at unit level
        self.heat_duty = Reference(self.cold_side.heat[:])

        # Add references to the user provided aliases (if applicable).
        add_hx_references(self)

        # ---------------------------------------------------------------------
        # Add performance equations
        # All units of measurement will be based on hot side
        hunits = self.config.hot_side.property_package.get_metadata().get_derived_units

        # Common heat exchanger variables
        self.area = Var(
            initialize=1,
            units=hunits("area"),
            domain=PositiveReals,
            doc="Heat transfer area",
        )

        self.heat_transfer_coefficient = Var(
            self.flowsheet().time,
            initialize=1,
            units=hunits("heat_transfer_coefficient"),
            domain=PositiveReals,
            doc="Overall heat transfer coefficient",
        )

        # Overall energy balance
        def rule_energy_balance(blk, t):
            return blk.hot_side.heat[t] == -pyunits.convert(
                blk.cold_side.heat[t], to_units=hunits("power")
            )

        self.energy_balance_constraint = Constraint(
            self.flowsheet().time, rule=rule_energy_balance
        )

        # Add e-NTU variables
        self.effectiveness = Var(
            self.flowsheet().time,
            initialize=1,
            units=pyunits.dimensionless,
            domain=PositiveReals,
            doc="Effectiveness factor for NTU method",
        )

        # Minimum heat capacitance ratio for e-NTU method
        self.eps_cmin = Param(
            initialize=1e-3,
            mutable=True,
            units=hunits("power") / hunits("temperature"),
            doc="Epsilon parameter for smooth Cmin and Cmax",
        )

        # TODO : Support both mass and mole based flows
        def rule_Cmin(blk, t):
            caph = (
                blk.hot_side.properties_in[t].flow_mol
                * blk.hot_side.properties_in[t].cp_mol
            )
            capc = pyunits.convert(
                blk.cold_side.properties_in[t].flow_mol
                * blk.cold_side.properties_in[t].cp_mol,
                to_units=hunits("power") / hunits("temperature"),
            )
            return smooth_min(caph, capc, eps=blk.eps_cmin)

        self.Cmin = Expression(
            self.flowsheet().time, rule=rule_Cmin, doc="Minimum heat capacitance rate"
        )

        def rule_Cmax(blk, t):
            caph = (
                blk.hot_side.properties_in[t].flow_mol
                * blk.hot_side.properties_in[t].cp_mol
            )
            capc = pyunits.convert(
                blk.cold_side.properties_in[t].flow_mol
                * blk.cold_side.properties_in[t].cp_mol,
                to_units=hunits("power") / hunits("temperature"),
            )
            return smooth_max(caph, capc, eps=blk.eps_cmin)

        self.Cmax = Expression(
            self.flowsheet().time, rule=rule_Cmax, doc="Maximum heat capacitance rate"
        )

        # Heat capacitance ratio
        def rule_Cratio(blk, t):
            return blk.Cmin[t] / blk.Cmax[t]

        self.Cratio = Expression(
            self.flowsheet().time, rule=rule_Cratio, doc="Heat capacitance ratio"
        )

        def rule_NTU(blk, t):
            return blk.heat_transfer_coefficient[t] * blk.area / blk.Cmin[t]

        self.NTU = Expression(
            self.flowsheet().time, rule=rule_NTU, doc="Number of heat transfer units"
        )

        # Heat transfer by e-NTU method
        def rule_entu(blk, t):
            return blk.hot_side.heat[t] == -(
                blk.effectiveness[t]
                * blk.Cmin[t]
                * (
                    blk.hot_side.properties_in[t].temperature
                    - pyunits.convert(
                        blk.cold_side.properties_in[t].temperature,
                        to_units=hunits("temperature"),
                    )
                )
            )

        self.heat_duty_constraint = Constraint(self.flowsheet().time, rule=rule_entu)

    # TODO : Add scaling methods

    def initialize_build(
        self,
        hot_side_state_args=None,
        cold_side_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        duty=None,
    ):
        """
        Heat exchanger initialization method.

        Args:
            hot_side_state_args : a dict of arguments to be passed to the
                property initialization for the hot side (see documentation of
                the specific property package) (default = None).
            cold_side_state_args : a dict of arguments to be passed to the
                property initialization for the cold side (see documentation of
                the specific property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            duty : an initial guess for the amount of heat transfered. This
                should be a tuple in the form (value, units), (default
                = (1000 J/s))

        Returns:
            None

        """
        # Set solver options
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        hot_side = self.hot_side
        cold_side = self.cold_side

        # Create solver
        opt = get_solver(solver, optarg)

        flags1 = hot_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=hot_side_state_args
        )

        init_log.info_high("Initialization Step 1a (hot side) Complete.")

        flags2 = cold_side.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=cold_side_state_args
        )

        init_log.info_high("Initialization Step 1b (cold side) Complete.")

        # ---------------------------------------------------------------------
        # Solve unit without heat transfer equation
        self.energy_balance_constraint.deactivate()

        # Get side 1 and side 2 heat units, and convert duty as needed
        s1_units = hot_side.heat.get_units()
        s2_units = cold_side.heat.get_units()

        if duty is None:
            # Assume 1000 J/s and check for unitless properties
            if s1_units is None and s2_units is None:
                # Backwards compatability for unitless properties
                s1_duty = -1000
                s2_duty = 1000
            else:
                s1_duty = pyunits.convert_value(
                    -1000, from_units=pyunits.W, to_units=s1_units
                )
                s2_duty = pyunits.convert_value(
                    1000, from_units=pyunits.W, to_units=s2_units
                )
        else:
            # Duty provided with explicit units
            s1_duty = -pyunits.convert_value(
                duty[0], from_units=duty[1], to_units=s1_units
            )
            s2_duty = pyunits.convert_value(
                duty[0], from_units=duty[1], to_units=s2_units
            )

        cold_side.heat.fix(s2_duty)
        for i in hot_side.heat:
            hot_side.heat[i].value = s1_duty

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        cold_side.heat.unfix()
        self.energy_balance_constraint.activate()

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        hot_side.release_state(flags1, outlvl=outlvl)
        cold_side.release_state(flags2, outlvl=outlvl)

        init_log.info("Initialization Completed, {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Hot Inlet": self.hot_side_inlet,
                "Hot Outlet": self.hot_side_outlet,
                "Cold Inlet": self.cold_side_inlet,
                "Cold Outlet": self.cold_side_outlet,
            },
            time_point=time_point,
        )
