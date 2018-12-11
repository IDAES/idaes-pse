##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Basic IDAES 1D Heat Exchanger Model.

1D Single pass shell and tube HX model with 0D/1D/2D wall conduction model
"""
from __future__ import division

# Import Pyomo libraries
from pyomo.environ import (SolverFactory, Var, Param, Constraint, Expression,
                           value, TransformationFactory)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume1D,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitBlockData,
                        useDefault)
from idaes.core.util.config import (is_physical_parameter_block,
                                    list_of_floats)
from idaes.core.util.misc import add_object_reference

__author__ = "Jaffer Ghouse"


@declare_process_block_class("HeatExchanger1D")
class HeatExchanger1DData(UnitBlockData):
    """Standard Heat Exchanger 1D Unit Model Class."""

    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic or not,
**default** = useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - set as a dynamic model,
**False** - set as a steady-state model.}"""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))
    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.componentPhase,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}"""))
    CONFIG.declare("energy_balance_type", ConfigValue(
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
    CONFIG.declare("shell_has_phase_equilibrium", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Phase equilibrium for shell side",
        doc="""Argument to enable phase equilibrium on the shell side.
- True - include phase equilibrium term
- False - do not include phase equilinrium term"""))
    CONFIG.declare("tube_has_phase_equilibrium", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Phase equilibrium for tube side",
        doc="""Argument to enable phase equilibrium on the tube side.
- True - include phase equilibrium term
- False - do not include phase equilinrium term"""))
    CONFIG.declare("shell_has_mass_diffusion", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Mass diffusion flag",
        doc="""Flag indicating  whether mass diffusion/dispersion should be
included in material balance equations (default=False)"""))
    CONFIG.declare("shell_has_energy_diffusion", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Energy diffusion flag",
        doc="""Flag indicating  whether energy diffusion/dispersion should be
included in energy balance equations (default=False)"""))
    CONFIG.declare("tube_has_mass_diffusion", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Mass diffusion flag",
        doc="""Flag indicating  whether mass diffusion/dispersion should be
included in material balance equations (default=False)"""))
    CONFIG.declare("tube_has_energy_diffusion", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Energy diffusion flag",
        doc="""Flag indicating  whether energy diffusion/dispersion should be
included in energy balance equations (default=False)"""))
    CONFIG.declare("shell_property_package", ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for shell control volume",
        doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
    - 'use_parent_value' - get package from parent (default = None)
    - a ParameterBlock object"""))
    CONFIG.declare("shell_property_package_args", ConfigValue(
        default={},
        description="Arguments for constructing shell property package",
        doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
    - 'use_parent_value' - get package from parent (default = None)
    - a dict (see property package for documentation)"""))
    CONFIG.declare("tube_property_package", ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for tube control volume",
        doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
    - 'use_parent_value' - get package from parent (default = None)
    - a ParameterBlock object"""))
    CONFIG.declare("tube_property_package_args", ConfigValue(
        default={},
        description="Arguments for constructing tube property package",
        doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
    - 'use_parent_value' - get package from parent (default = None)
    - a dict (see property package for documentation)"""))
    CONFIG.declare("shell_inlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of inlet names for shell side of heat exchanger",
        doc="""A list containing names of inlets for shell of heat exchanger
(default = None)
    - None - default single inlet
    - list - a list of names for inlets"""))
    CONFIG.declare("shell_num_inlets", ConfigValue(
        domain=int,
        description="Number of inlets to shell side",
        doc="""Argument indication number (int) of inlets to construct
(default = None). Not used if inlet_list arg is provided.
    - None - use inlet_list arg instead
    - int - Inlets will be named with sequential numbers from 1
            to num_inlets"""))
    CONFIG.declare("shell_outlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of outlet names for shell side of heat exchanger",
        doc="""A list containing names of outlets for shell of heat exchanger
(default = None)
    - None - default single inlet
    - list - a list of names for inlets"""))
    CONFIG.declare("shell_num_outlets", ConfigValue(
        domain=int,
        description="Number of outlets to shell side",
        doc="""Argument indication number (int) of outlets to construct
(default = None). Not used if outlet_list arg is provided.
    - None - use outlet_list arg instead
    - int - Outlets will be named with sequential numbers from 1
            to num_outlets"""))
    CONFIG.declare("tube_inlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of inlet names for tube side of heat exchanger",
        doc="""A list containing names of inlets for tube of heat exchanger
(default = None)
    - None - default single inlet
    - list - a list of names for inlets"""))
    CONFIG.declare("tube_num_inlets", ConfigValue(
        domain=int,
        description="Number of inlets to tube side",
        doc="""Argument indication number (int) of inlets to construct
(default = None). Not used if inlet_list arg is provided.
    - None - use inlet_list arg instead
    - int - Inlets will be named with sequential numbers from 1
            to num_inlets"""))
    CONFIG.declare("tube_outlet_list", ConfigValue(
        domain=list_of_strings,
        description="List of outlet names for tube side of heat exchanger",
        doc="""A list containing names of outlets for tube of heat exchanger
(default = None)
    - None - default single inlet
    - list - a list of names for inlets"""))
    CONFIG.declare("tube_num_outlets", ConfigValue(
        domain=int,
        description="Number of outlets to tube side",
        doc="""Argument indication number (int) of outlets to construct
(default = None). Not used if outlet_list arg is provided.
    - None - use outlet_list arg instead
    - int - Outlets will be named with sequential numbers from 1
            to num_outlets"""))
    # TODO : We should probably think about adding a consistency check for the
    # TODO : discretisation methdos as well.
    CONFIG.declare("shell_discretization_method", ConfigValue(
        default="OCLR",
        domain=In(['OCLR', 'OCLL', 'BFD', 'FFD']),
        description="Discretization method to apply for length domain",
        doc="""Method to be used by DAE transformation when discretizing length
domain on shell side (default = `OCLR`).
- 'OCLR' - orthogonal collocation (Radau roots)
- 'OCLL' - orthogonal collocation (Legendre roots)
- 'BFD' - backwards finite difference (1st order)
- 'FFD' - forwards finite difference (1st order)"""))
    CONFIG.declare("tube_discretization_method", ConfigValue(
        default="OCLR",
        domain=In(['OCLR', 'OCLL', 'BFD', 'FFD']),
        description="Discretization method to apply for length domain",
        doc="""Method to be used by DAE transformation when discretizing
length domain on tube side (default = `OCLR`).
- 'OCLR' - orthogonal collocation (Radau roots)
- 'OCLL' - orthogonal collocation (Legendre roots)
- 'BFD' - backwards finite difference (1st order)
- 'FFD' - forwards finite difference (1st order)"""))
    CONFIG.declare("finite_elements", ConfigValue(
        default=20,
        domain=int,
        description="Number of finite elements length domain",
        doc="""Number of finitie elements to use when discretizing length
domain (default=20)"""))
    CONFIG.declare("collocation_points", ConfigValue(
        default=5,
        domain=int,
        description="Number of collocation points per finite element",
        doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)"""))
    CONFIG.declare("flow_type", ConfigValue(
        default="co_current",
        domain=In(['co_current', 'counter_current']),
        description="Flow configuration of heat exchanger",
        doc="""Flow configuration of heat exchanger
- co_current: shell and tube flows from 0 to 1
- counter_current: shell side flows from 0 to 1
tube side flows from 1 to 0"""))
    CONFIG.declare("has_wall_conduction", ConfigValue(
        default="none",
        domain=In(["none", "1D", "2D"]),
        description="Conduction model for tube wall",
        doc="""Argument to enable type of wall heat conduction model.
- none - 0D wall model
- 1D - 1D wall model along the thickness of the tube
- 2D - 2D wall model along the lenghth and thickness of the tube"""))


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

        # Set flow directions for the control volume blocks
        if self.config.flow_type == "co_current":
            set_direction_shell = "forward"
            set_direction_tube = "forward"
        else:
            set_direction_shell = "forward"
            set_direction_tube = "backward"

        # Control volume 1D for shell
        self.shell = ControlVolume1D(default={
            "has_phase_equilibrium": self.config.shell_has_phase_equilibrium,
            "has_mass_diffusion": self.config.shell_has_mass_diffusion,
            "has_energy_diffusion": self.config.shell_has_energy_diffusion,
            "property_package": self.config.shell_property_package,
            "property_package_args": self.config.shell_property_package_args,
            "flow_direction": set_direction_shell,
            "discretization_method": self.config.shell_discretization_method,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points})

        self.tube = ControlVolume1D(default={
            "has_phase_equilibrium": self.config.tube_has_phase_equilibrium,
            "has_mass_diffusion": self.config.tube_has_mass_diffusion,
            "has_energy_diffusion": self.config.tube_has_energy_diffusion,
            "property_package": self.config.tube_property_package,
            "property_package_args": self.config.tube_property_package_args,
            "flow_direction": set_direction_tube,
            "discretization_method": self.config.tube_discretization_method,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points})

        self.shell.add_geometry(
            length_domain_set=self.config.length_domain_set)
        self.tube.add_geometry(
            length_domain_set=self.config.length_domain_set)

        self.shell.add_state_blocks()
        self.tube.add_state_blocks()

        self.shell.add_material_balances(
            balance_type=self.config.material_balance_type)

        self.shell.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer)

        self.shell.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        self.shell.apply_transformation(
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points)

        self.tube.add_material_balances(
            balance_type=self.config.material_balance_type)

        self.tube.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer)

        self.tube.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        self.tube.apply_transformation(
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points)

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

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
        self.shell_heat_transfer_coefficient = Var(self.time,
                                                   self.shell.ldomain,
                                                   initialize=50,
                                                   doc="Heat transfer "
                                                   "coefficient")
        self.tube_heat_transfer_coefficient = Var(self.time,
                                                  self.tube.ldomain,
                                                  initialize=50,
                                                  doc="Heat transfer "
                                                  "coefficient")

        # Wall 0D model (Q_shell = Q_tube*N_tubes)
        if self.config.has_wall_conduction == "none":
            self.temperature_wall = Var(self.time, self.tube.ldomain,
                                        initialize=298.15)

            # Performance equations
            # Energy transfer between shell and tube wall

            @self.Constraint(self.time, self.shell.ldomain,
                             doc="Heat transfer between shell and tube")
            def shell_heat_transfer_eq(self, t, x):
                return self.shell.heat[t, x] == - self.N_tubes *\
                    (self.shell_heat_transfer_coefficient[t, x] *
                     self.pi * self.d_tube_inner *
                     (self.shell.properties[t, x].temperature -
                      self.temperature_wall[t, x]))

            # Energy transfer between tube wall and tube
            @self.Constraint(self.time, self.tube.ldomain,
                             doc="Convective heat transfer")
            def tube_heat_transfer_eq(self, t, x):
                return self.tube.heat[t, x] == \
                    self.tube_heat_transfer_coefficient[t, x] *\
                    self.pi * self.d_tube_inner * \
                    (self.temperature_wall[t, x] -
                     self.tube.properties[t, x].temperature)

            # Wall 0D model
            @self.Constraint(self.time, self.shell.ldomain,
                             doc="wall 0D model")
            def wall_0D_model(self, t, x):
                return self.tube.heat[t, x] == -(self.shell.heat[t, x] /
                                                 self.N_tubes)
