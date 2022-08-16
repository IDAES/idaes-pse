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
Basic IDAES 1D Heat Exchanger Model.

1D Single pass shell and tube HX model with 0D wall conduction model
"""
# Import Python libraries
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    check_optimal_termination,
    Constraint,
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
from idaes.core.util.constants import Constants as c
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


__author__ = "Jaffer Ghouse"

# Set up logger
_log = idaeslog.getLogger(__name__)


class WallConductionType(Enum):
    zero_dimensional = 0
    one_dimensional = 1
    two_dimensional = 2


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
        "has_wall_conduction",
        ConfigValue(
            default=WallConductionType.zero_dimensional,
            domain=In(WallConductionType),
            description="Conduction model for cold side wall",
            doc="""Argument to enable type of wall heat conduction model.
- WallConductionType.zero_dimensional - 0D wall model (default),
- WallConductionType.one_dimensional - 1D wall model along the thickness of the
tube,
- WallConductionType.two_dimensional - 2D wall model along the lenghth and
thickness of the tube""",
        ),
    )
    CONFIG.declare(
        "hot_side_name",
        ConfigValue(
            default="shell",
            domain=str,
            doc="Hot side name, sets control volume and inlet and outlet names. "
            "Default = 'shell'.",
        ),
    )
    CONFIG.declare(
        "cold_side_name",
        ConfigValue(
            default="tube",
            domain=str,
            doc="Cold side name, sets control volume and inlet and outlet names. "
            "Default = 'tube'.",
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
        hx_process_config(self)

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
            default={
                "dynamic": self.config.hot_side.dynamic,
                "has_holdup": self.config.hot_side.has_holdup,
                "property_package": self.config.hot_side.property_package,
                "property_package_args": self.config.hot_side.property_package_args,
                "transformation_method": self.config.hot_side.transformation_method,
                "transformation_scheme": self.config.hot_side.transformation_scheme,
                "finite_elements": self.config.finite_elements,
                "collocation_points": self.config.collocation_points,
            }
        )

        self.cold_side = ControlVolume1DBlock(
            default={
                "dynamic": self.config.cold_side.dynamic,
                "has_holdup": self.config.cold_side.has_holdup,
                "property_package": self.config.cold_side.property_package,
                "property_package_args": self.config.cold_side.property_package_args,
                "transformation_method": self.config.cold_side.transformation_method,
                "transformation_scheme": self.config.cold_side.transformation_scheme,
                "finite_elements": self.config.finite_elements,
                "collocation_points": self.config.collocation_points,
            }
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

        # Add reference to control volume geometry
        add_object_reference(self, "hot_side_area", self.hot_side.area)
        add_object_reference(self, "hot_side_length", self.hot_side.length)
        add_object_reference(self, "cold_side_area", self.cold_side.area)
        add_object_reference(self, "cold_side_length", self.cold_side.length)

        # Add references to the user provided aliases if applicable
        add_hx_references(self)

        self._make_performance()

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
        cold_side_units = (
            self.config.cold_side.property_package.get_metadata().get_derived_units
        )

        # Unit model variables
        # HX dimensions
        self.d_hot_side = Var(
            initialize=1, doc="Diameter of hot side", units=hot_side_units("length")
        )
        self.d_cold_side_outer = Var(
            initialize=0.011,
            doc="Outer diameter of cold side",
            units=hot_side_units("length"),
        )
        self.d_cold_side_inner = Var(
            initialize=0.010,
            doc="Inner diameter of cold side",
            units=hot_side_units("length"),
        )
        self.N_tubes = Var(
            initialize=1, doc="Number of tubes", units=pyunits.dimensionless
        )

        # Note: In addition to the above variables, "hot_side_length" and
        # "cold_side_length" need to be fixed at the flowsheet level

        # Performance variables
        self.hot_side_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            self.hot_side.length_domain,
            initialize=50,
            doc="Heat transfer coefficient",
            units=hot_side_units("heat_transfer_coefficient"),
        )
        self.cold_side_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            self.cold_side.length_domain,
            initialize=50,
            doc="Heat transfer coefficient",
            units=cold_side_units("heat_transfer_coefficient"),
        )

        # Wall 0D model (Q_hot_side = Q_cold_side*N_tubes)
        if self.config.has_wall_conduction == WallConductionType.zero_dimensional:
            self.temperature_wall = Var(
                self.flowsheet().time,
                self.cold_side.length_domain,
                initialize=298.15,
                units=hot_side_units("temperature"),
            )

            # Performance equations
            # Energy transfer between hot_side and cold_side wall

            @self.Constraint(
                self.flowsheet().time,
                self.hot_side.length_domain,
                doc="Heat transfer between hot_side and cold_side",
            )
            def hot_side_heat_transfer_eq(self, t, x):
                return self.hot_side.heat[t, x] == -self.N_tubes * (
                    self.hot_side_heat_transfer_coefficient[t, x]
                    * c.pi
                    * self.d_cold_side_outer
                    * (
                        self.hot_side.properties[t, x].temperature
                        - self.temperature_wall[t, x]
                    )
                )

            # Energy transfer between cold_side wall and cold_side
            @self.Constraint(
                self.flowsheet().time,
                self.cold_side.length_domain,
                doc="Convective heat transfer",
            )
            def cold_side_heat_transfer_eq(self, t, x):
                return self.cold_side.heat[
                    t, x
                ] == self.cold_side_heat_transfer_coefficient[
                    t, x
                ] * c.pi * pyunits.convert(
                    self.d_cold_side_inner, to_units=cold_side_units("length")
                ) * (
                    pyunits.convert(
                        self.temperature_wall[t, x],
                        to_units=cold_side_units("temperature"),
                    )
                    - self.cold_side.properties[t, x].temperature
                )

            if hot_side_units("length") is None:
                # Backwards compatability check
                q_units = None
            else:
                q_units = hot_side_units("power") / hot_side_units("length")
            # Wall 0D model
            @self.Constraint(
                self.flowsheet().time,
                self.hot_side.length_domain,
                doc="wall 0D model",
            )
            def wall_0D_model(self, t, x):
                return pyunits.convert(
                    self.cold_side.heat[t, x], to_units=q_units
                ) == -(self.hot_side.heat[t, x] / self.N_tubes)

        else:
            raise NotImplementedError(
                "{} HeatExchanger1D has not yet implemented support for "
                "wall conduction models."
            )

        # Define cold_side area in terms of tube diameter
        self.area_calc_cold_side = Constraint(
            expr=4 * self.cold_side_area
            == c.pi
            * pyunits.convert(
                self.d_cold_side_inner, to_units=cold_side_units("length")
            )
            ** 2
        )

        # Define hot_side area in terms of hot_side and tube diameter
        self.area_calc_hot_side = Constraint(
            expr=4 * self.hot_side_area
            == c.pi
            * (self.d_hot_side**2 - self.N_tubes * self.d_cold_side_outer**2)
        )

    def initialize_build(
        self,
        hot_side_state_args=None,
        cold_side_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
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

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize hot_side block
        flags_hot_side = self.hot_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=hot_side_state_args,
        )

        flags_cold_side = self.cold_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=cold_side_state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        # Wall 0D
        if self.config.has_wall_conduction == WallConductionType.zero_dimensional:
            hot_side_units = (
                self.config.hot_side.property_package.get_metadata().get_derived_units
            )
            for t in self.flowsheet().time:
                for z in self.hot_side.length_domain:
                    self.temperature_wall[t, z].fix(
                        value(
                            0.5
                            * (
                                self.hot_side.properties[t, 0].temperature
                                + pyunits.convert(
                                    self.cold_side.properties[t, 0].temperature,
                                    to_units=hot_side_units("temperature"),
                                )
                            )
                        )
                    )

            self.cold_side.deactivate()
            self.cold_side_heat_transfer_eq.deactivate()
            self.wall_0D_model.deactivate()

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(self, tee=slc.tee)
            init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(res))
            )

            self.cold_side.activate()
            self.cold_side_heat_transfer_eq.activate()

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(self, tee=slc.tee)
            init_log.info_high(
                "Initialization Step 3 {}.".format(idaeslog.condition(res))
            )

            self.wall_0D_model.activate()
            self.temperature_wall.unfix()

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(self, tee=slc.tee)
            init_log.info_high(
                "Initialization Step 4 {}.".format(idaeslog.condition(res))
            )
        else:
            res = None

        self.hot_side.release_state(flags_hot_side)
        self.cold_side.release_state(flags_cold_side)

        if res is not None and not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info("Initialization Complete.")

    def _get_performance_contents(self, time_point=0):
        # TODO: Set this up to use user names if available
        var_dict = {}
        var_dict["Hot Side Area"] = self.hot_side.area
        var_dict["Hot Side Diameter"] = self.d_hot_side
        var_dict["Hot Side Length"] = self.hot_side.length
        var_dict["Cold Side Area"] = self.cold_side.area
        var_dict["Cold Side Outer Diameter"] = self.d_cold_side_outer
        var_dict["Cold Side Inner Diameter"] = self.d_cold_side_inner
        var_dict["Cold Side Length"] = self.cold_side.length
        var_dict["Number of Tubes"] = self.N_tubes

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        # TODO : Set this up to use user provided names if available
        return create_stream_table_dataframe(
            {
                "Hot Side Inlet": self.hot_side_inlet,
                "Hot Side Outlet": self.hot_side_outlet,
                "Cold Side Inlet": self.cold_side_inlet,
                "Cold Side Outlet": self.cold_side_outlet,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for i, c in self.hot_side_heat_transfer_eq.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.hot_side.heat[i], default=1, warning=True
                ),
                overwrite=False,
            )

        for i, c in self.cold_side_heat_transfer_eq.items():
            iscale.constraint_scaling_transform(
                c,
                iscale.get_scaling_factor(
                    self.cold_side.heat[i], default=1, warning=True
                ),
                overwrite=False,
            )
