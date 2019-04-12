##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Basic IDAES 1D Heat Exchanger Model.

1D Single pass shell and tube HX model with 0D wall conduction model
"""
from __future__ import division

# Import Python libraries
import math
import logging
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (SolverFactory, Var, Param, Constraint,
                           value, TerminationCondition)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume1DBlock,
                        UnitModelBlockData,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        FlowDirection,
                        useDefault)
from idaes.unit_models.heat_exchanger import HeatExchangerFlowPattern
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import ConfigurationError

__author__ = "Jaffer Ghouse"

# Set up logger
_log = logging.getLogger(__name__)


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
    _SideTemplate.declare("dynamic", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}"""))
    _SideTemplate.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))
    _SideTemplate.declare("material_balance_type", ConfigValue(
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
**MaterialBalanceType.total** - use total material balance.}"""))
    _SideTemplate.declare("energy_balance_type", ConfigValue(
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
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    _SideTemplate.declare("momentum_balance_type", ConfigValue(
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
    _SideTemplate.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    _SideTemplate.declare("has_phase_equilibrium", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Phase equilibrium term construction flag",
        doc="""Argument to enable phase equilibrium on the shell side.
- True - include phase equilibrium term
- False - do not include phase equilinrium term"""))
    _SideTemplate.declare("property_package", ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object"""))
    _SideTemplate.declare("property_package_args", ConfigValue(
        default={},
        description="Arguments for constructing shell property package",
        doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)"""))
    # TODO : We should probably think about adding a consistency check for the
    # TODO : discretisation methods as well.
    _SideTemplate.declare("transformation_method", ConfigValue(
        default=useDefault,
        description="Discretization method to use for DAE transformation",
        doc="""Discretization method to use for DAE transformation. See Pyomo
documentation for supported transformations."""))
    _SideTemplate.declare("transformation_scheme", ConfigValue(
        default=useDefault,
        description="Discretization scheme to use for DAE transformation",
        doc="""Discretization scheme to use when transformating domain. See
Pyomo documentation for supported schemes."""))

    # Create individual config blocks for shell and tube side
    CONFIG.declare("shell_side",
                   _SideTemplate(doc="shell side config arguments"))
    CONFIG.declare("tube_side",
                   _SideTemplate(doc="tube side config arguments"))

    # Common config args for both sides
    CONFIG.declare("finite_elements", ConfigValue(
        default=20,
        domain=int,
        description="Number of finite elements length domain",
        doc="""Number of finite elements to use when discretizing length
domain (default=20)"""))
    CONFIG.declare("collocation_points", ConfigValue(
        default=5,
        domain=int,
        description="Number of collocation points per finite element",
        doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)"""))
    CONFIG.declare("flow_type", ConfigValue(
        default=HeatExchangerFlowPattern.cocurrent,
        domain=In(HeatExchangerFlowPattern),
        description="Flow configuration of heat exchanger",
        doc="""Flow configuration of heat exchanger
- HeatExchangerFlowPattern.cocurrent: shell and tube flows from 0 to 1
- HeatExchangerFlowPattern.countercurrent: shell side flows from 0 to 1
tube side flows from 1 to 0"""))
    CONFIG.declare("has_wall_conduction", ConfigValue(
        default=WallConductionType.zero_dimensional,
        domain=In(WallConductionType),
        description="Conduction model for tube wall",
        doc="""Argument to enable type of wall heat conduction model.
- WallConductionType.zero_dimensional - 0D wall model,
- WallConductionType.one_dimensional - 1D wall model along the thickness of the
tube,
- WallConductionType.two_dimensional - 2D wall model along the lenghth and
thickness of the tube"""))

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
            if self.config.shell_side.transformation_method is useDefault:
                _log.warning("Discretization method was "
                             "not specified for the shell side of the "
                             "co-current heat exchanger. "
                             "Defaulting to finite "
                             "difference method on the shell side.")
                self.config.shell_side.transformation_method = \
                    "dae.finite_difference"
            if self.config.tube_side.transformation_method is useDefault:
                _log.warning("Discretization method was "
                             "not specified for the tube side of the "
                             "co-current heat exchanger. "
                             "Defaulting to finite "
                             "difference method on the tube side.")
                self.config.tube_side.transformation_method = \
                    "dae.finite_difference"
            if self.config.shell_side.transformation_scheme is useDefault:
                _log.warning("Discretization scheme was "
                             "not specified for the shell side of the "
                             "co-current heat exchanger. "
                             "Defaulting to backward finite "
                             "difference on the shell side.")
                self.config.shell_side.transformation_scheme = "BACKWARD"
            if self.config.tube_side.transformation_scheme is useDefault:
                _log.warning("Discretization scheme was "
                             "not specified for the tube side of the "
                             "co-current heat exchanger. "
                             "Defaulting to backward finite "
                             "difference on the tube side.")
                self.config.tube_side.transformation_scheme = "BACKWARD"
        elif self.config.flow_type == HeatExchangerFlowPattern.countercurrent:
            set_direction_shell = FlowDirection.forward
            set_direction_tube = FlowDirection.backward
            if self.config.shell_side.transformation_method is useDefault:
                _log.warning("Discretization method was "
                             "not specified for the shell side of the "
                             "counter-current heat exchanger. "
                             "Defaulting to finite "
                             "difference method on the shell side.")
                self.config.shell_side.transformation_method = \
                    "dae.finite_difference"
            if self.config.tube_side.transformation_method is useDefault:
                _log.warning("Discretization method was "
                             "not specified for the tube side of the "
                             "counter-current heat exchanger. "
                             "Defaulting to finite "
                             "difference method on the tube side.")
                self.config.tube_side.transformation_method = \
                    "dae.finite_difference"
            if self.config.shell_side.transformation_scheme is useDefault:
                _log.warning("Discretization scheme was "
                             "not specified for the shell side of the "
                             "counter-current heat exchanger. "
                             "Defaulting to backward finite "
                             "difference on the shell side.")
                self.config.shell_side.transformation_scheme = "BACKWARD"
            if self.config.tube_side.transformation_scheme is useDefault:
                _log.warning("Discretization scheme was "
                             "not specified for the tube side of the "
                             "counter-current heat exchanger. "
                             "Defaulting to forward finite "
                             "difference on the tube side.")
                self.config.tube_side.transformation_scheme = "FORWARD"
        else:
            raise ConfigurationError(
                    "{} HeatExchanger1D only supports cocurrent and "
                    "countercurrent flow patterns, but flow_type configuration"
                    " argument was set to {}."
                    .format(self.name,
                            self.config.flow_type))

        # Control volume 1D for shell
        self.shell = ControlVolume1DBlock(default={
            "dynamic": self.config.shell_side.dynamic,
            "has_holdup": self.config.shell_side.has_holdup,
            "property_package": self.config.shell_side.property_package,
            "property_package_args":
                self.config.shell_side.property_package_args})

        self.tube = ControlVolume1DBlock(default={
            "dynamic": self.config.tube_side.dynamic,
            "has_holdup": self.config.tube_side.has_holdup,
            "property_package": self.config.tube_side.property_package,
            "property_package_args":
                self.config.tube_side.property_package_args})

        self.shell.add_geometry(flow_direction=set_direction_shell)
        self.tube.add_geometry(flow_direction=set_direction_tube)

        self.shell.add_state_blocks(
            information_flow=set_direction_shell,
            has_phase_equilibrium=self.config.shell_side.has_phase_equilibrium)
        self.tube.add_state_blocks(
            information_flow=set_direction_tube,
            has_phase_equilibrium=self.config.tube_side.has_phase_equilibrium)

        # Populate shell
        self.shell.add_material_balances(
            balance_type=self.config.shell_side.material_balance_type,
            has_phase_equilibrium=self.config.shell_side.has_phase_equilibrium)

        self.shell.add_energy_balances(
            balance_type=self.config.shell_side.energy_balance_type,
            has_heat_transfer=True)

        self.shell.add_momentum_balances(
            balance_type=self.config.shell_side.momentum_balance_type,
            has_pressure_change=self.config.shell_side.has_pressure_change)

        self.shell.apply_transformation(
            transformation_method=self.config.shell_side.transformation_method,
            transformation_scheme=self.config.shell_side.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points)

        # Populate tube
        self.tube.add_material_balances(
            balance_type=self.config.tube_side.material_balance_type,
            has_phase_equilibrium=self.config.tube_side.has_phase_equilibrium)

        self.tube.add_energy_balances(
            balance_type=self.config.tube_side.energy_balance_type,
            has_heat_transfer=True)

        self.tube.add_momentum_balances(
            balance_type=self.config.tube_side.momentum_balance_type,
            has_pressure_change=self.config.tube_side.has_pressure_change)

        self.tube.apply_transformation(
            transformation_method=self.config.tube_side.transformation_method,
            transformation_scheme=self.config.tube_side.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points)

        # Add Ports for shell side
        self.add_inlet_port(name="shell_inlet", block=self.shell)
        self.add_outlet_port(name="shell_outlet", block=self.shell)

        # Add Ports for tube side
        self.add_inlet_port(name="tube_inlet", block=self.tube)
        self.add_outlet_port(name="tube_outlet", block=self.tube)

        # Add reference to control volume geometry
        add_object_reference(self,
                             "shell_area",
                             self.shell.area)
        add_object_reference(self,
                             "shell_length",
                             self.shell.length)
        add_object_reference(self,
                             "tube_area",
                             self.tube.area)
        add_object_reference(self,
                             "tube_length",
                             self.tube.length)

        self._make_performance()

    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """
        # Unit model parameters
        self.pi = Param(initialize=math.pi, doc="pi")

        # Unit model variables
        # HX dimensions
        self.d_shell = Var(initialize=1, doc="Diameter of shell")
        self.d_tube_outer = Var(initialize=0.011, doc="Outer diameter of tube")
        self.d_tube_inner = Var(initialize=0.010, doc="Inner diameter of tube")
        self.N_tubes = Var(initialize=1, doc="Number of tubes")

        # Note: In addition to the above variables, "shell_length" and
        # "tube_length" need to be fixed at the flowsheet level

        # Performance variables
        self.shell_heat_transfer_coefficient = Var(self.time_ref,
                                                   self.shell.length_domain,
                                                   initialize=50,
                                                   doc="Heat transfer "
                                                   "coefficient")
        self.tube_heat_transfer_coefficient = Var(self.time_ref,
                                                  self.tube.length_domain,
                                                  initialize=50,
                                                  doc="Heat transfer "
                                                  "coefficient")

        # Wall 0D model (Q_shell = Q_tube*N_tubes)
        if (self.config.has_wall_conduction ==
                WallConductionType.zero_dimensional):
            self.temperature_wall = Var(self.time_ref, self.tube.length_domain,
                                        initialize=298.15)

            # Performance equations
            # Energy transfer between shell and tube wall

            @self.Constraint(self.time_ref, self.shell.length_domain,
                             doc="Heat transfer between shell and tube")
            def shell_heat_transfer_eq(self, t, x):
                return self.shell.heat[t, x] == - self.N_tubes *\
                    (self.shell_heat_transfer_coefficient[t, x] *
                     self.pi * self.d_tube_outer *
                     (self.shell.properties[t, x].temperature -
                      self.temperature_wall[t, x]))

            # Energy transfer between tube wall and tube
            @self.Constraint(self.time_ref, self.tube.length_domain,
                             doc="Convective heat transfer")
            def tube_heat_transfer_eq(self, t, x):
                return self.tube.heat[t, x] == \
                    self.tube_heat_transfer_coefficient[t, x] *\
                    self.pi * self.d_tube_inner * \
                    (self.temperature_wall[t, x] -
                     self.tube.properties[t, x].temperature)

            # Wall 0D model
            @self.Constraint(self.time_ref, self.shell.length_domain,
                             doc="wall 0D model")
            def wall_0D_model(self, t, x):
                return self.tube.heat[t, x] == -(self.shell.heat[t, x] /
                                                 self.N_tubes)

        else:
            raise NotImplementedError(
                    "{} HeatExchanger1D has not yet implemented support for "
                    "wall conduction models.")

        # Define tube area in terms of tube diameter
        self.area_calc_tube = Constraint(expr=4 * self.tube_area == self.pi *
                                         self.d_tube_inner**2)

        # Define shell area in terms of shell and tube diameter
        self.area_calc_shell = Constraint(expr=4 * self.shell_area ==
                                          self.pi * (self.d_shell**2 -
                                                     self.N_tubes *
                                                     self.d_tube_outer**2))

    def initialize(blk, shell_state_args=None, tube_state_args=None, outlvl=1,
                   solver='ipopt', optarg={'tol': 1e-6}):
        """
        Initialisation routine for the unit (default solver ipopt).

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = return solver state for each step in subroutines
                     * 3 = include solver output infomation (tee=True)

            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        """
        # Set solver options
        if outlvl > 3:
            stee = True
        else:
            stee = False

        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialize shell block
        flags_shell = blk.shell.initialize(outlvl=outlvl - 1,
                                           optarg=optarg,
                                           solver=solver,
                                           state_args=shell_state_args)

        flags_tube = blk.tube.initialize(outlvl=outlvl - 1,
                                         optarg=optarg,
                                         solver=solver,
                                         state_args=tube_state_args)

        if outlvl > 0:
            _log.info('{} Initialisation Step 1 Complete.'.format(blk.name))

        # ---------------------------------------------------------------------
        # Solve unit
        # Wall 0D
        if blk.config.has_wall_conduction == \
                WallConductionType.zero_dimensional:
            for t in blk.time_ref:
                for z in blk.shell.length_domain:
                    blk.temperature_wall[t, z].fix(value(
                        0.5 * (blk.shell.properties[t, 0].temperature +
                               blk.tube.properties[t, 0].temperature)))

            blk.tube.deactivate()
            blk.tube_heat_transfer_eq.deactivate()
            blk.wall_0D_model.deactivate()

            results = opt.solve(blk, tee=stee)

            if outlvl > 0:
                if results.solver.termination_condition \
                        == TerminationCondition.optimal:
                    _log.info('{} Initialisation Step 2 Complete.'
                              .format(blk.name))
                else:
                    _log.warning('{} Initialisation Step 2 Failed.'
                                 .format(blk.name))

            blk.tube.activate()
            blk.tube_heat_transfer_eq.activate()

            results = opt.solve(blk, tee=stee)
            if outlvl > 0:
                if results.solver.termination_condition \
                        == TerminationCondition.optimal:
                    _log.info('{} Initialisation Step 3 Complete.'
                              .format(blk.name))
                else:
                    _log.warning('{} Initialisation Step 3 Failed.'
                                 .format(blk.name))
            blk.wall_0D_model.activate()
            blk.temperature_wall.unfix()

            results = opt.solve(blk, tee=stee)
            if outlvl > 0:
                if results.solver.termination_condition \
                        == TerminationCondition.optimal:
                    _log.info('{} Initialisation Step 4 Complete.'
                              .format(blk.name))
                else:
                    _log.warning('{} Initialisation Step 4 Failed.'
                                 .format(blk.name))

        blk.shell.release_state(flags_shell)
        blk.tube.release_state(flags_tube)

        if outlvl > 0:
            _log.info('{} Initialisation Complete.'.format(blk.name))
