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
Generic IDAES 1D Heat Exchanger Model with overall area and heat transfer coefficient
"""
# Import Pyomo libraries
from pyomo.environ import (
    Var,
    check_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    ControlVolume1DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
    useDefault,
)
from idaes.models.unit_models.heat_exchanger import (
    HeatExchangerFlowPattern,
    hx_process_config,
    add_hx_references,
)
from idaes.core.util.config import is_physical_parameter_block, DefaultBool
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


__author__ = "Jaffer Ghouse, Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("HeatExchanger1D")
class HeatExchanger1DData(UnitModelBlockData):
    """Standard Heat Exchanger 1D Unit Model Class."""

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    # Template for config arguments for hot and cold side
    _SideTemplate = ConfigBlock()
    _SideTemplate.declare(
        "dynamic",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}""",
        ),
    )
    _SideTemplate.declare(
        "has_holdup",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )
    _SideTemplate.declare(
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
    _SideTemplate.declare(
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
    _SideTemplate.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should
be constructed,
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
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Phase equilibrium term construction flag",
            doc="""Argument to enable phase equilibrium.
- True - include phase equilibrium term
- False - do not include phase equilibrium term""",
        ),
    )
    _SideTemplate.declare(
        "property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property
calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object""",
        ),
    )
    _SideTemplate.declare(
        "property_package_args",
        ConfigValue(
            default={},
            description="Arguments for constructing property package",
            doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)""",
        ),
    )
    # TODO : We should probably think about adding a consistency check for the
    # TODO : discretisation methods as well.
    _SideTemplate.declare(
        "transformation_method",
        ConfigValue(
            default=useDefault,
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See
Pyomo documentation for supported transformations.""",
        ),
    )
    _SideTemplate.declare(
        "transformation_scheme",
        ConfigValue(
            default=useDefault,
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transformating domain. See
Pyomo documentation for supported schemes.""",
        ),
    )

    # Create individual config blocks for hot and cold side
    CONFIG.declare("hot_side", _SideTemplate(doc="hot side config arguments"))
    CONFIG.declare("cold_side", _SideTemplate(doc="cold side config arguments"))

    # Common config args for both sides
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=20,
            domain=int,
            description="Number of finite elements length domain",
            doc="""Number of finite elements to use when discretizing length
domain (default=20)""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)""",
        ),
    )
    CONFIG.declare(
        "flow_type",
        ConfigValue(
            default=HeatExchangerFlowPattern.cocurrent,
            domain=In(HeatExchangerFlowPattern),
            description="Flow configuration of heat exchanger",
            doc="""Flow configuration of heat exchanger
- HeatExchangerFlowPattern.cocurrent: hot and cold flows from 0 to 1
(default)
- HeatExchangerFlowPattern.countercurrent: hot side flows from 0 to 1
cold side flows from 1 to 0""",
        ),
    )
    CONFIG.declare(
        "hot_side_name",
        ConfigValue(
            default=None,
            domain=str,
            doc="Hot side name, sets control volume and inlet and outlet names. "
            "Default = None.",
        ),
    )
    CONFIG.declare(
        "cold_side_name",
        ConfigValue(
            default=None,
            domain=str,
            doc="Cold side name, sets control volume and inlet and outlet names. "
            "Default = None.",
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
        self._process_config()

        # Set flow directions for the control volume blocks and specify
        # dicretisation if not specified.
        if self.config.flow_type == HeatExchangerFlowPattern.cocurrent:
            set_direction_hot = FlowDirection.forward
            set_direction_cold = FlowDirection.forward
            if (
                self.config.hot_side.transformation_method
                != self.config.cold_side.transformation_method
            ) or (
                self.config.hot_side.transformation_scheme
                != self.config.cold_side.transformation_scheme
            ):
                raise ConfigurationError(
                    "HeatExchanger1D only supports similar transformation "
                    "schemes on the hot and cold side domains for "
                    "both cocurrent and countercurrent flow patterns."
                )
            if self.config.hot_side.transformation_method is useDefault:
                _log.warning(
                    "Discretization method was "
                    "not specified for the hot side of the "
                    "co-current heat exchanger. "
                    "Defaulting to finite "
                    "difference method on the hot side."
                )
                self.config.hot_side.transformation_method = "dae.finite_difference"
            if self.config.cold_side.transformation_method is useDefault:
                _log.warning(
                    "Discretization method was "
                    "not specified for the cold side of the "
                    "co-current heat exchanger. "
                    "Defaulting to finite "
                    "difference method on the cold side."
                )
                self.config.cold_side.transformation_method = "dae.finite_difference"
            if self.config.hot_side.transformation_scheme is useDefault:
                _log.warning(
                    "Discretization scheme was "
                    "not specified for the hot side of the "
                    "co-current heat exchanger. "
                    "Defaulting to backward finite "
                    "difference on the hot side."
                )
                self.config.hot_side.transformation_scheme = "BACKWARD"
            if self.config.cold_side.transformation_scheme is useDefault:
                _log.warning(
                    "Discretization scheme was "
                    "not specified for the cold side of the "
                    "co-current heat exchanger. "
                    "Defaulting to backward finite "
                    "difference on the cold side."
                )
                self.config.cold_side.transformation_scheme = "BACKWARD"
        elif self.config.flow_type == HeatExchangerFlowPattern.countercurrent:
            set_direction_hot = FlowDirection.forward
            set_direction_cold = FlowDirection.backward
            if self.config.hot_side.transformation_method is useDefault:
                _log.warning(
                    "Discretization method was "
                    "not specified for the hot side of the "
                    "counter-current heat exchanger. "
                    "Defaulting to finite "
                    "difference method on the hot side."
                )
                self.config.hot_side.transformation_method = "dae.finite_difference"
            if self.config.cold_side.transformation_method is useDefault:
                _log.warning(
                    "Discretization method was "
                    "not specified for the cold side of the "
                    "counter-current heat exchanger. "
                    "Defaulting to finite "
                    "difference method on the cold side."
                )
                self.config.cold_side.transformation_method = "dae.finite_difference"
            if self.config.hot_side.transformation_scheme is useDefault:
                _log.warning(
                    "Discretization scheme was "
                    "not specified for the hot side of the "
                    "counter-current heat exchanger. "
                    "Defaulting to backward finite "
                    "difference on the hot side."
                )
                self.config.hot_side.transformation_scheme = "BACKWARD"
            if self.config.cold_side.transformation_scheme is useDefault:
                _log.warning(
                    "Discretization scheme was "
                    "not specified for the cold side of the "
                    "counter-current heat exchanger. "
                    "Defaulting to forward finite "
                    "difference on the cold side."
                )
                self.config.cold_side.transformation_scheme = "BACKWARD"
        else:
            raise ConfigurationError(
                "{} HeatExchanger1D only supports cocurrent and "
                "countercurrent flow patterns, but flow_type configuration"
                " argument was set to {}.".format(self.name, self.config.flow_type)
            )

        # Control volume 1D for hot
        self.hot_side = ControlVolume1DBlock(
            dynamic=self.config.hot_side.dynamic,
            has_holdup=self.config.hot_side.has_holdup,
            property_package=self.config.hot_side.property_package,
            property_package_args=self.config.hot_side.property_package_args,
            transformation_method=self.config.hot_side.transformation_method,
            transformation_scheme=self.config.hot_side.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )

        self.cold_side = ControlVolume1DBlock(
            dynamic=self.config.cold_side.dynamic,
            has_holdup=self.config.cold_side.has_holdup,
            property_package=self.config.cold_side.property_package,
            property_package_args=self.config.cold_side.property_package_args,
            transformation_method=self.config.cold_side.transformation_method,
            transformation_scheme=self.config.cold_side.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
        )

        self.hot_side.add_geometry(flow_direction=set_direction_hot)
        self.cold_side.add_geometry(flow_direction=set_direction_cold)

        self.hot_side.add_state_blocks(
            information_flow=set_direction_hot,
            has_phase_equilibrium=self.config.hot_side.has_phase_equilibrium,
        )
        self.cold_side.add_state_blocks(
            information_flow=set_direction_cold,
            has_phase_equilibrium=self.config.cold_side.has_phase_equilibrium,
        )

        # Populate hot
        self.hot_side.add_material_balances(
            balance_type=self.config.hot_side.material_balance_type,
            has_phase_equilibrium=self.config.hot_side.has_phase_equilibrium,
        )

        self.hot_side.add_energy_balances(
            balance_type=self.config.hot_side.energy_balance_type,
            has_heat_transfer=True,
        )

        self.hot_side.add_momentum_balances(
            balance_type=self.config.hot_side.momentum_balance_type,
            has_pressure_change=self.config.hot_side.has_pressure_change,
        )

        self.hot_side.apply_transformation()

        # Populate cold side
        self.cold_side.add_material_balances(
            balance_type=self.config.cold_side.material_balance_type,
            has_phase_equilibrium=self.config.cold_side.has_phase_equilibrium,
        )

        self.cold_side.add_energy_balances(
            balance_type=self.config.cold_side.energy_balance_type,
            has_heat_transfer=True,
        )

        self.cold_side.add_momentum_balances(
            balance_type=self.config.cold_side.momentum_balance_type,
            has_pressure_change=self.config.cold_side.has_pressure_change,
        )

        self.cold_side.apply_transformation()

        # Add Ports for hot side
        self.add_inlet_port(name="hot_side_inlet", block=self.hot_side)
        self.add_outlet_port(name="hot_side_outlet", block=self.hot_side)

        # Add Ports for cold side
        self.add_inlet_port(name="cold_side_inlet", block=self.cold_side)
        self.add_outlet_port(name="cold_side_outlet", block=self.cold_side)

        # Add references to the user provided aliases if applicable
        add_hx_references(self)

        self._make_geometry()
        self._make_performance()

        hot_side_units = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )
        q_units = hot_side_units("power") / hot_side_units("length")

        @self.Constraint(
            self.flowsheet().time,
            self.hot_side.length_domain,
            doc="Heat conservation equality",
        )
        def heat_conservation(self, t, x):
            return pyunits.convert(self.cold_side.heat[t, x], to_units=q_units) == -(
                self.hot_side.heat[t, x]
            )

    def _process_config(self):
        hx_process_config(self)

    def _make_geometry(self):
        hot_side_units = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )

        self.area = Var(
            initialize=1, units=hot_side_units("area"), doc="Heat transfer area"
        )

        # Add reference to control volume geometry
        add_object_reference(self, "length", self.hot_side.length)

        # Equate hot and cold side lengths
        @self.Constraint(doc="Equating hot and cold side lengths")
        def length_equality(self):
            return (
                pyunits.convert(
                    self.cold_side.length, to_units=hot_side_units("length")
                )
                == self.hot_side.length
            )

    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """
        hot_side_units = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )

        # Performance variables
        self.heat_transfer_coefficient = Var(
            self.flowsheet().time,
            self.hot_side.length_domain,
            initialize=50,
            doc="Average heat transfer coefficient",
            units=hot_side_units("heat_transfer_coefficient"),
        )

        @self.Constraint(
            self.flowsheet().time,
            self.hot_side.length_domain,
            doc="Heat transfer between hot_side and cold_side",
        )
        def heat_transfer_eq(self, t, x):
            return self.hot_side.heat[t, x] == -(
                self.heat_transfer_coefficient[t, x]
                * self.area
                / self.length
                * (
                    self.hot_side.properties[t, x].temperature
                    - pyunits.convert(
                        self.cold_side.properties[t, x].temperature,
                        to_units=hot_side_units("temperature"),
                    )
                )
            )

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
        Initialization routine for the unit.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            duty : an initial guess for the amount of heat transferred. This
                should be a tuple in the form (value, units). A default value
                is calculated based on stream temperatures, the overall
                heat transfer coefficient, and exchanger area

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Get length values
        if self.length.fixed:
            # Most likely case
            self.cold_side.length.set_value(self.length)
        elif self.cold_side.length.fixed:
            # This would be unusual, but check
            self.length.set_value(self.cold_side.length)
        else:
            # No fixed value for length - we will assume the user knows what they are doing
            pass

        # Initialize control volumes blocks
        Lfix = self.hot_side.length.fixed
        self.hot_side.length.fix()
        flags_hot_side = self.hot_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=hot_side_state_args,
        )
        if not Lfix:
            self.hot_side.length.unfix()

        Lfix = self.cold_side.length.fixed
        self.cold_side.length.fix()
        flags_cold_side = self.cold_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=cold_side_state_args,
        )
        if not Lfix:
            self.cold_side.length.unfix()

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit with fixed heat duty
        # Guess heat duty based on 1/4 of maximum driving force
        hot_side_units = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )
        cold_side_units = (
            self.config.cold_side.property_package.get_metadata().get_derived_units
        )
        if duty is None:
            duty = value(
                0.25
                * self.heat_transfer_coefficient[0, 0]
                * self.area
                * (
                    self.hot_side.properties[0, 0].temperature
                    - pyunits.convert(
                        self.cold_side.properties[0, 0].temperature,
                        to_units=hot_side_units("temperature"),
                    )
                )
            )
        else:
            duty = pyunits.convert_value(
                duty[0], from_units=duty[1], to_units=hot_side_units("power")
            )
        duty_per_length = value(duty / self.length)
        # Fix heat duties
        for v in self.hot_side.heat.values():
            v.fix(-duty_per_length)
        for v in self.cold_side.heat.values():
            v.fix(
                pyunits.convert_value(
                    duty_per_length,
                    to_units=cold_side_units("power") / cold_side_units("length"),
                    from_units=hot_side_units("power") / hot_side_units("length"),
                )
            )

        # Deactivate heat duty constraints
        self.heat_transfer_eq.deactivate()
        self.heat_conservation.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Unfix heat duty and reactivate constraints
        for v in self.hot_side.heat.values():
            v.unfix()
        for v in self.cold_side.heat.values():
            v.unfix()
        self.heat_transfer_eq.activate()
        self.heat_conservation.activate()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        self.hot_side.release_state(flags_hot_side)
        self.cold_side.release_state(flags_cold_side)

        if res is not None and not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info("Initialization Complete.")

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Area"] = self.area
        var_dict["Length"] = self.length

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        # Get names for hot and cold sides
        hot_name = self.config.hot_side_name
        if hot_name is None:
            hot_name = "Hot Side"
        cold_name = self.config.cold_side_name
        if cold_name is None:
            cold_name = "Cold Side"
        return create_stream_table_dataframe(
            {
                f"{hot_name} Inlet": self.hot_side_inlet,
                f"{hot_name} Outlet": self.hot_side_outlet,
                f"{cold_name} Inlet": self.cold_side_inlet,
                f"{cold_name} Outlet": self.cold_side_outlet,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for i, c in self.heat_transfer_eq.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.hot_side.heat[i], default=1, warning=True
                ),
                overwrite=False,
            )

        for i, c in self.heat_conservation.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.hot_side.heat[i], default=1, warning=True
                ),
                overwrite=False,
            )
