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
    SolverFactory,
    Var,
    Constraint,
    value,
    units as pyunits
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

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
from idaes.generic_models.unit_models.heat_exchanger \
    import HeatExchangerFlowPattern
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants as c
from idaes.core.util import get_solver, scaling as iscale

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

    CONFIG = UnitModelBlockData.CONFIG()
    # Template for config arguments for shell and tube side
    _SideTemplate = ConfigBlock()
    _SideTemplate.declare(
        "dynamic",
        ConfigValue(
            default=useDefault,
            domain=In([useDefault, True, False]),
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
            domain=In([useDefault, True, False]),
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
            domain=In([True, False]),
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
            domain=In([True, False]),
            description="Phase equilibrium term construction flag",
            doc="""Argument to enable phase equilibrium on the shell side.
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
            description="Arguments for constructing shell property package",
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

    # Create individual config blocks for shell and tube side
    CONFIG.declare(
        "shell_side", _SideTemplate(doc="shell side config arguments"))
    CONFIG.declare(
        "tube_side", _SideTemplate(doc="tube side config arguments"))

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
- HeatExchangerFlowPattern.cocurrent: shell and tube flows from 0 to 1
(default)
- HeatExchangerFlowPattern.countercurrent: shell side flows from 0 to 1
tube side flows from 1 to 0""",
        ),
    )
    CONFIG.declare(
        "has_wall_conduction",
        ConfigValue(
            default=WallConductionType.zero_dimensional,
            domain=In(WallConductionType),
            description="Conduction model for tube wall",
            doc="""Argument to enable type of wall heat conduction model.
- WallConductionType.zero_dimensional - 0D wall model (default),
- WallConductionType.one_dimensional - 1D wall model along the thickness of the
tube,
- WallConductionType.two_dimensional - 2D wall model along the lenghth and
thickness of the tube""",
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
        super(HeatExchanger1DData, self).build()

        # Set flow directions for the control volume blocks and specify
        # dicretisation if not specified.
        if self.config.flow_type == HeatExchangerFlowPattern.cocurrent:
            set_direction_shell = FlowDirection.forward
            set_direction_tube = FlowDirection.forward
            if (
                self.config.shell_side.transformation_method
                != self.config.tube_side.transformation_method
            ) or (
                self.config.shell_side.transformation_scheme
                != self.config.tube_side.transformation_scheme
            ):
                raise ConfigurationError(
                    "HeatExchanger1D only supports similar transformation "
                    "schemes on the shell side and tube side domains for "
                    "both cocurrent and countercurrent flow patterns."
                )
            if self.config.shell_side.transformation_method is useDefault:
                _log.warning(
                    "Discretization method was "
                    "not specified for the shell side of the "
                    "co-current heat exchanger. "
                    "Defaulting to finite "
                    "difference method on the shell side."
                )
                self.config.shell_side.transformation_method = \
                    "dae.finite_difference"
            if self.config.tube_side.transformation_method is useDefault:
                _log.warning(
                    "Discretization method was "
                    "not specified for the tube side of the "
                    "co-current heat exchanger. "
                    "Defaulting to finite "
                    "difference method on the tube side."
                )
                self.config.tube_side.transformation_method = \
                    "dae.finite_difference"
            if self.config.shell_side.transformation_scheme is useDefault:
                _log.warning(
                    "Discretization scheme was "
                    "not specified for the shell side of the "
                    "co-current heat exchanger. "
                    "Defaulting to backward finite "
                    "difference on the shell side."
                )
                self.config.shell_side.transformation_scheme = "BACKWARD"
            if self.config.tube_side.transformation_scheme is useDefault:
                _log.warning(
                    "Discretization scheme was "
                    "not specified for the tube side of the "
                    "co-current heat exchanger. "
                    "Defaulting to backward finite "
                    "difference on the tube side."
                )
                self.config.tube_side.transformation_scheme = "BACKWARD"
        elif self.config.flow_type == HeatExchangerFlowPattern.countercurrent:
            set_direction_shell = FlowDirection.forward
            set_direction_tube = FlowDirection.backward
            if self.config.shell_side.transformation_method is useDefault:
                _log.warning(
                    "Discretization method was "
                    "not specified for the shell side of the "
                    "counter-current heat exchanger. "
                    "Defaulting to finite "
                    "difference method on the shell side."
                )
                self.config.shell_side.transformation_method = \
                    "dae.finite_difference"
            if self.config.tube_side.transformation_method is useDefault:
                _log.warning(
                    "Discretization method was "
                    "not specified for the tube side of the "
                    "counter-current heat exchanger. "
                    "Defaulting to finite "
                    "difference method on the tube side."
                )
                self.config.tube_side.transformation_method = \
                    "dae.finite_difference"
            if self.config.shell_side.transformation_scheme is useDefault:
                _log.warning(
                    "Discretization scheme was "
                    "not specified for the shell side of the "
                    "counter-current heat exchanger. "
                    "Defaulting to backward finite "
                    "difference on the shell side."
                )
                self.config.shell_side.transformation_scheme = "BACKWARD"
            if self.config.tube_side.transformation_scheme is useDefault:
                _log.warning(
                    "Discretization scheme was "
                    "not specified for the tube side of the "
                    "counter-current heat exchanger. "
                    "Defaulting to forward finite "
                    "difference on the tube side."
                )
                self.config.tube_side.transformation_scheme = "BACKWARD"
        else:
            raise ConfigurationError(
                "{} HeatExchanger1D only supports cocurrent and "
                "countercurrent flow patterns, but flow_type configuration"
                " argument was set to {}."
                .format(self.name, self.config.flow_type)
            )

        # Control volume 1D for shell
        self.shell = ControlVolume1DBlock(
            default={
                "dynamic": self.config.shell_side.dynamic,
                "has_holdup": self.config.shell_side.has_holdup,
                "property_package": self.config.shell_side.property_package,
                "property_package_args":
                self.config.shell_side.property_package_args,
                "transformation_method":
                self.config.shell_side.transformation_method,
                "transformation_scheme":
                self.config.shell_side.transformation_scheme,
                "finite_elements": self.config.finite_elements,
                "collocation_points": self.config.collocation_points,
            }
        )

        self.tube = ControlVolume1DBlock(
            default={
                "dynamic": self.config.tube_side.dynamic,
                "has_holdup": self.config.tube_side.has_holdup,
                "property_package": self.config.tube_side.property_package,
                "property_package_args":
                self.config.tube_side.property_package_args,
                "transformation_method":
                self.config.tube_side.transformation_method,
                "transformation_scheme":
                self.config.tube_side.transformation_scheme,
                "finite_elements": self.config.finite_elements,
                "collocation_points": self.config.collocation_points,
            }
        )

        self.shell.add_geometry(flow_direction=set_direction_shell)
        self.tube.add_geometry(flow_direction=set_direction_tube)

        self.shell.add_state_blocks(
            information_flow=set_direction_shell,
            has_phase_equilibrium=self.config.shell_side.has_phase_equilibrium,
        )
        self.tube.add_state_blocks(
            information_flow=set_direction_tube,
            has_phase_equilibrium=self.config.tube_side.has_phase_equilibrium,
        )

        # Populate shell
        self.shell.add_material_balances(
            balance_type=self.config.shell_side.material_balance_type,
            has_phase_equilibrium=self.config.shell_side.has_phase_equilibrium,
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
            has_phase_equilibrium=self.config.tube_side.has_phase_equilibrium,
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

        # Add reference to control volume geometry
        add_object_reference(self, "shell_area", self.shell.area)
        add_object_reference(self, "shell_length", self.shell.length)
        add_object_reference(self, "tube_area", self.tube.area)
        add_object_reference(self, "tube_length", self.tube.length)

        self._make_performance()

    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """
        shell_units = self.config.shell_side.property_package.\
            get_metadata().get_derived_units
        tube_units = self.config.tube_side.property_package.\
            get_metadata().get_derived_units

        # Unit model variables
        # HX dimensions
        self.d_shell = Var(initialize=1,
                           doc="Diameter of shell",
                           units=shell_units("length"))
        self.d_tube_outer = Var(initialize=0.011,
                                doc="Outer diameter of tube",
                           units=shell_units("length"))
        self.d_tube_inner = Var(initialize=0.010,
                                doc="Inner diameter of tube",
                           units=shell_units("length"))
        self.N_tubes = Var(initialize=1,
                           doc="Number of tubes",
                           units=pyunits.dimensionless)

        # Note: In addition to the above variables, "shell_length" and
        # "tube_length" need to be fixed at the flowsheet level

        # Performance variables
        self.shell_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            self.shell.length_domain,
            initialize=50,
            doc="Heat transfer coefficient",
            units=shell_units("heat_transfer_coefficient")
        )
        self.tube_heat_transfer_coefficient = Var(
            self.flowsheet().time,
            self.tube.length_domain,
            initialize=50,
            doc="Heat transfer coefficient",
            units=tube_units("heat_transfer_coefficient")
        )

        # Wall 0D model (Q_shell = Q_tube*N_tubes)
        if self.config.has_wall_conduction == \
                WallConductionType.zero_dimensional:
            self.temperature_wall = Var(
                self.flowsheet().time,
                self.tube.length_domain,
                initialize=298.15,
                units=shell_units("temperature")
            )

            # Performance equations
            # Energy transfer between shell and tube wall

            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="Heat transfer between shell and tube",
            )
            def shell_heat_transfer_eq(self, t, x):
                return self.shell.heat[t, x] == -self.N_tubes * (
                    self.shell_heat_transfer_coefficient[t, x]
                    * c.pi
                    * self.d_tube_outer
                    * (
                        self.shell.properties[t, x].temperature
                        - self.temperature_wall[t, x]
                    )
                )

            # Energy transfer between tube wall and tube
            @self.Constraint(
                self.flowsheet().time,
                self.tube.length_domain,
                doc="Convective heat transfer",
            )
            def tube_heat_transfer_eq(self, t, x):
                return self.tube.heat[t, x] == \
                    self.tube_heat_transfer_coefficient[
                    t, x
                ] * c.pi * pyunits.convert(self.d_tube_inner,
                                           to_units=tube_units("length")) * (
                    pyunits.convert(self.temperature_wall[t, x],
                                    to_units=tube_units('temperature')) -
                    self.tube.properties[t, x].temperature
                )

            if shell_units("length") is None:
                # Backwards compatability check
                q_units = None
            else:
                q_units = shell_units("power")/shell_units("length")
            # Wall 0D model
            @self.Constraint(
                self.flowsheet().time,
                self.shell.length_domain,
                doc="wall 0D model",
            )
            def wall_0D_model(self, t, x):
                return pyunits.convert(self.tube.heat[t, x],
                                       to_units=q_units) == -(
                            self.shell.heat[t, x] / self.N_tubes)

        else:
            raise NotImplementedError(
                "{} HeatExchanger1D has not yet implemented support for "
                "wall conduction models."
            )

        # Define tube area in terms of tube diameter
        self.area_calc_tube = Constraint(
            expr=4 * self.tube_area == c.pi * pyunits.convert(
                self.d_tube_inner, to_units=tube_units("length"))**2
        )

        # Define shell area in terms of shell and tube diameter
        self.area_calc_shell = Constraint(
            expr=4 * self.shell_area
            == c.pi * (self.d_shell**2 - self.N_tubes*self.d_tube_outer**2)
        )

    def initialize(
        self,
        shell_state_args=None,
        tube_state_args=None,
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
        # Initialize shell block
        flags_shell = self.shell.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=shell_state_args,
        )

        flags_tube = self.tube.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=tube_state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        # Wall 0D
        if self.config.has_wall_conduction == \
            WallConductionType.zero_dimensional:
            shell_units = self.config.shell_side.property_package.\
                get_metadata().get_derived_units
            for t in self.flowsheet().time:
                for z in self.shell.length_domain:
                    self.temperature_wall[t, z].fix(
                        value(
                            0.5
                            * (
                                self.shell.properties[t, 0].temperature
                                + pyunits.convert(
                                    self.tube.properties[t, 0].temperature,
                                    to_units=shell_units('temperature'))
                            )
                        )
                    )

            self.tube.deactivate()
            self.tube_heat_transfer_eq.deactivate()
            self.wall_0D_model.deactivate()

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(self, tee=slc.tee)
            init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(res))
            )

            self.tube.activate()
            self.tube_heat_transfer_eq.activate()

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

        self.shell.release_state(flags_shell)
        self.tube.release_state(flags_tube)

        init_log.info("Initialization Complete.")

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Shell Area"] = self.shell.area
        var_dict["Shell Diameter"] = self.d_shell
        var_dict["Shell Length"] = self.shell.length
        var_dict["Tube Area"] = self.tube.area
        var_dict["Tube Outer Diameter"] = self.d_tube_outer
        var_dict["Tube Inner Diameter"] = self.d_tube_inner
        var_dict["Tube Length"] = self.tube.length
        var_dict["Number of Tubes"] = self.N_tubes

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Shell Inlet": self.shell_inlet,
                "Shell Outlet": self.shell_outlet,
                "Tube Inlet": self.tube_inlet,
                "Tube Outlet": self.tube_outlet,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for i, c in self.shell_heat_transfer_eq.items():
            iscale.constraint_scaling_transform(
                c, iscale.get_scaling_factor(
                    self.shell.heat[i], default=1, warning=True),
                overwrite=False)

        for i, c in self.tube_heat_transfer_eq.items():
            iscale.constraint_scaling_transform(c, iscale.get_scaling_factor(
                    self.tube.heat[i], default=1, warning=True),
                overwrite=False)
