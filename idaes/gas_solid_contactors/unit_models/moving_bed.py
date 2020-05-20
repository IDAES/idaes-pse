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
IDAES Moving Bed Model.

The moving bed model is a 1D axially discretized model with a gas and
solid phase and a counter-current flow direction. The model captures
the gas-solid interaction between both phases through reaction, mass
and heat transfer.

Equations written in this model were derived from:
A. Ostace, A. Lee, C.O. Okoli, A.P. Burgard, D.C. Miller, D. Bhattacharyya,
Mathematical modeling of a moving-bed reactor for chemical looping combustion
of methane, in: M.R. Eden, M. Ierapetritou, G.P. Towler (Eds.),13th Int. Symp.
Process Syst. Eng. (PSE 2018), Computer-Aided Chemical Engineering 2018,
pp. 325–330 , San Diego, CA.

Assumptions:
Property package contains temperature and pressure variables.
Property package contains minimum fluidization velocity.

"""
from __future__ import division

# Import Python libraries
import math
import matplotlib.pyplot as plt

# Import Pyomo libraries
from pyomo.environ import (SolverFactory, Var, Param, Reals, value,
                           TransformationFactory, Constraint,
                           TerminationCondition)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.dae import ContinuousSet

# Import IDAES cores
from idaes.core import (ControlVolume1DBlock, UnitModelBlockData,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        FlowDirection,
                        useDefault)
from idaes.core.util.config import (is_physical_parameter_block,
                                    is_reaction_parameter_block)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.control_volume1d import DistributedVars
import idaes.logger as idaeslog


__author__ = "Chinedu Okoli", "Anca Ostace"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("MBR")
class MBRData(UnitModelBlockData):
    """Standard Moving Bed Unit Model Class."""

    # Create template for unit level config arguments
    CONFIG = ConfigBlock()

    # Create template for phase/side specific config arguments
    _SideTemplate = ConfigBlock()

    # Unit level config arguments
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
    CONFIG.declare("finite_elements", ConfigValue(
        default=10,
        domain=int,
        description="Number of finite elements length domain",
        doc="""Number of finite elements to use when discretizing length
domain (default=20)"""))
    CONFIG.declare("length_domain_set", ConfigValue(
        default=[0.0, 1.0],
        domain=list,
        description="Number of finite elements length domain",
        doc="""length_domain_set - (optional) list of point to use to
initialize a new ContinuousSet if length_domain is not
provided (default = [0.0, 1.0]).
domain (default = [0.0, 1.0])"""))
    CONFIG.declare("transformation_method", ConfigValue(
        default="dae.finite_difference",
        description="Method to use for DAE transformation",
        doc="""Method to use to transform domain. Must be a method recognised
by the Pyomo TransformationFactory,
**default** - "dae.finite_difference".
**Valid values:** {
**"dae.finite_difference"** - Use a finite difference transformation method,
**"dae.collocation"** - use a collocation transformation method}"""))
    CONFIG.declare("transformation_scheme", ConfigValue(
        default=None,
        domain=In([None, "BACKWARD", "FORWARD", "LAGRANGE-RADAU"]),
        description="Scheme to use for DAE transformation",
        doc="""Scheme to use when transforming domain. See Pyomo
documentation for supported schemes,
**default** - None.
**Valid values:** {
**None** - defaults to "BACKWARD" for finite difference transformation method,
and to "LAGRANGE-RADAU" for collocation transformation method,
**"BACKWARD"** - Use a finite difference transformation method,
**"FORWARD""** - use a finite difference transformation method,
**"LAGRANGE-RADAU""** - use a collocation transformation method}"""))
    CONFIG.declare("collocation_points", ConfigValue(
        default=3,
        domain=int,
        description="Number of collocation points per finite element",
        doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)"""))
    CONFIG.declare("flow_type", ConfigValue(
        default="counter_current",
        domain=In(['counter_current']),
        description="Flow configuration of Moving Bed",
        doc="""Flow configuration of Moving Bed
- counter_current: gas side flows from 0 to 1
solid side flows from 1 to 0"""))
    CONFIG.declare("material_balance_type", ConfigValue(
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
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    CONFIG.declare("area_definition", ConfigValue(
        default=DistributedVars.variant,
        domain=In(DistributedVars),
        description="Argument for defining form of area variable",
        doc="""Argument defining whether area variable should be spatially
        variant or not. **default** - DistributedVars.uniform.
        **Valid values:** {
        DistributedVars.uniform - area does not vary across spatial domian,
        DistributedVars.variant - area can vary over the domain and is indexed
        by time and space.}"""))

    # Populate the side template to default values
    _SideTemplate.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.none,
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
    _SideTemplate.declare("pressure_drop_type", ConfigValue(
        default="simple_correlation",
        domain=In(["simple_correlation", "ergun_correlation"]),
        description="Construction flag for type of pressure drop",
        doc="""Indicates what type of pressure drop correlation should be used,
**default** - "simple_correlation".
**Valid values:** {
**"simple_correlation"** - Use a simplified pressure drop correlation,
**"ergun_correlation"** - Use the ergun equation.}"""))
    _SideTemplate.declare("has_equilibrium_reactions", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Equilibrium reaction construction flag",
        doc="""Indicates whether terms for equilibrium controlled reactions
should be constructed,
**default** - True.
**Valid values:** {
**True** - include equilibrium reaction terms,
**False** - exclude equilibrium reaction terms.}"""))
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
        description="Arguments for constructing gas property package",
        doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)"""))
    _SideTemplate.declare("reaction_package", ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package to use for control volume",
        doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}"""))
    _SideTemplate.declare("reaction_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}"""))

    # Create individual config blocks for gas and solid sides
    CONFIG.declare("gas_phase_config",
                   _SideTemplate(doc="gas side config arguments"))
    CONFIG.declare("solid_phase_config",
                   _SideTemplate(doc="solid side config arguments"))

    # Set gas side momentum balance to pressureTotal as default
    CONFIG.gas_phase_config.momentum_balance_type = \
        MomentumBalanceType.pressureTotal

    # =========================================================================
    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to build default attributes
        super(MBRData, self).build()

    # =========================================================================
        """ Set argument values for gas and solid sides"""
        # Consistency check for transformation method and transformation scheme
        if (self.config.transformation_method == "dae.finite_difference" and
                self.config.transformation_scheme is None):
            self.config.transformation_scheme = "BACKWARD"
        elif (self.config.transformation_method == "dae.collocation" and
                self.config.transformation_scheme is None):
            self.config.transformation_scheme = "LAGRANGE-RADAU"
        elif (self.config.transformation_method == "dae.finite_difference" and
              self.config.transformation_scheme == "LAGRANGE-RADAU"):
            raise ConfigurationError("{} invalid value for"
                                     "transformation_scheme argument."
                                     "Must be ""BACKWARD"" or ""FORWARD"" "
                                     "if transformation_method is"
                                     " ""dae.finite_difference""."
                                     .format(self.name))
        elif (self.config.transformation_method == "dae.collocation" and
                self.config.transformation_scheme != "LAGRANGE-RADAU"):
            raise ConfigurationError("{} invalid value for "
                                     "transformation_scheme argument."
                                     "Must be ""LAGRANGE-RADAU"" if "
                                     "transformation_method is"
                                     " ""dae.collocation""."
                                     .format(self.name))

        # Create a unit model length domain
        self.length_domain = ContinuousSet(
                                bounds=(0.0, 1.0),
                                initialize=self.config.length_domain_set,
                                doc="Normalized length domain")

        # Set flow directions for the control volume blocks
        # Gas flows from 0 to 1, solid flows from 1 to 0
        if self.config.flow_type == "counter_current":
            set_direction_gas = FlowDirection.forward
            set_direction_solid = FlowDirection.backward

        # Set arguments for gas sides if homoogeneous reaction block
        has_rate_reaction_gas_phase = False
        if self.config.gas_phase_config.reaction_package is not None:
            has_rate_reaction_gas_phase = True

        # Set arguments for gas and solid sides if heterogeneous reaction block
        has_rate_reaction_solid_phase = False
        has_mass_transfer_gas_phase = False
        if self.config.solid_phase_config.reaction_package is not None:
            has_rate_reaction_solid_phase = True
            has_mass_transfer_gas_phase = True

        # Set heat transfer terms
        has_heat_transfer = False
        if self.config.energy_balance_type != EnergyBalanceType.none:
            has_heat_transfer = True

        # Set heat of reaction terms
        has_heat_of_reaction_gas_phase = False
        has_heat_of_reaction_solid_phase = False
        if (self.config.energy_balance_type != EnergyBalanceType.none
                and self.config.gas_phase_config.reaction_package is not None):
            has_heat_of_reaction_gas_phase = True

        if (self.config.energy_balance_type != EnergyBalanceType.none
                and self.config.solid_phase_config.
                reaction_package is not None):
            has_heat_of_reaction_solid_phase = True

    # =========================================================================
        """ Build Control volume 1D for gas phase and
            populate gas control volume"""

        self.gas_phase = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": self.config.area_definition,
            "property_package": self.config.gas_phase_config.property_package,
            "property_package_args":
                self.config.gas_phase_config.property_package_args,
            "reaction_package": self.config.gas_phase_config.reaction_package,
            "reaction_package_args":
                self.config.gas_phase_config.reaction_package_args})

        self.gas_phase.add_geometry(
                length_domain=self.length_domain,
                length_domain_set=self.config.length_domain_set,
                flow_direction=set_direction_gas)

        self.gas_phase.add_state_blocks(
            information_flow=set_direction_gas,
            has_phase_equilibrium=False)

        if self.config.gas_phase_config.reaction_package is not None:
            self.gas_phase.add_reaction_blocks(
                    has_equilibrium=self.config.gas_phase_config.
                    has_equilibrium_reactions)

        self.gas_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
            has_mass_transfer=has_mass_transfer_gas_phase,
            has_rate_reactions=has_rate_reaction_gas_phase)

        self.gas_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=has_heat_transfer,
            has_heat_of_reaction=has_heat_of_reaction_gas_phase)

        self.gas_phase.add_momentum_balances(
            balance_type=self.config.gas_phase_config.momentum_balance_type,
            has_pressure_change=self.config.gas_phase_config.
            has_pressure_change)

    # =========================================================================
        """ Build Control volume 1D for solid phase and
            populate solid control volume"""

        # Set argument for heterogeneous reaction block
        self.solid_phase = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": self.config.area_definition,
            "property_package":
                self.config.solid_phase_config.property_package,
            "property_package_args":
                self.config.solid_phase_config.property_package_args,
            "reaction_package":
                self.config.solid_phase_config.reaction_package,
            "reaction_package_args":
                self.config.solid_phase_config.reaction_package_args})

        self.solid_phase.add_geometry(
                length_domain=self.length_domain,
                length_domain_set=self.config.length_domain_set,
                flow_direction=set_direction_solid)

        self.solid_phase.add_state_blocks(
            information_flow=set_direction_solid,
            has_phase_equilibrium=False)

        if self.config.solid_phase_config.reaction_package is not None:
            # TODO - a generalization of the heterogeneous reaction block
            # The heterogeneous reaction block does not use the
            # add_reaction_blocks in control volumes as control volumes are
            # currently setup to handle only homogeneous reaction properties.
            # Thus appending the heterogeneous reaction block to the
            # solid state block is currently hard coded here.

            tmp_dict = dict(**self.config.solid_phase_config.
                            reaction_package_args)
            tmp_dict["gas_state_block"] = self.gas_phase.properties
            tmp_dict["solid_state_block"] = self.solid_phase.properties
            tmp_dict["has_equilibrium"] = (self.config.solid_phase_config.
                                           has_equilibrium_reactions)
            tmp_dict["parameters"] = (self.config.solid_phase_config.
                                      reaction_package)
            self.solid_phase.reactions = (
                    self.config.solid_phase_config.reaction_package.
                    reaction_block_class(
                        self.flowsheet().config.time,
                        self.length_domain,
                        doc="Reaction properties in control volume",
                        default=tmp_dict))

        self.solid_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
            has_mass_transfer=False,
            has_rate_reactions=has_rate_reaction_solid_phase)

        self.solid_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=has_heat_transfer,
            has_heat_of_reaction=has_heat_of_reaction_solid_phase)

        self.solid_phase.add_momentum_balances(
            balance_type=self.config.solid_phase_config.momentum_balance_type,
            has_pressure_change=self.config.solid_phase_config.
            has_pressure_change)

    # =========================================================================
        """ Add ports"""
        # Add Ports for gas side
        self.add_inlet_port(name="gas_inlet", block=self.gas_phase)
        self.add_outlet_port(name="gas_outlet", block=self.gas_phase)

        # Add Ports for solid side
        self.add_inlet_port(name="solid_inlet", block=self.solid_phase)
        self.add_outlet_port(name="solid_outlet", block=self.solid_phase)

    # =========================================================================
        """ Add performace equation method"""

        self._apply_transformation()
        self._make_performance()

    # =========================================================================
    def _apply_transformation(self):
        """
        Method to apply DAE transformation to the Control Volume length domain.
        Transformation applied will be based on the Control Volume
        configuration arguments.
        """
        if self.config.finite_elements is None:
            raise ConfigurationError(
                    "{} was not provided a value for the finite_elements"
                    " configuration argument. Please provide a valid value."
                    .format(self.name))

        if self.config.transformation_method == "dae.finite_difference":
            self.discretizer = TransformationFactory(
                                    self.config.transformation_method)
            self.discretizer.apply_to(self,
                                      wrt=self.length_domain,
                                      nfe=self.config.finite_elements,
                                      scheme=self.config.transformation_scheme)
        elif self.config.transformation_method == "dae.collocation":
            self.discretizer = TransformationFactory(
                                    self.config.transformation_method)
            self.discretizer.apply_to(
                self,
                wrt=self.length_domain,
                nfe=self.config.finite_elements,
                ncp=self.config.collocation_points,
                scheme=self.config.transformation_scheme)

    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """

    # =========================================================================
        # Add object references - Sets
        add_object_reference(self,
                             "gas_component_list_ref",
                             self.config.gas_phase_config.property_package.
                             component_list)
        add_object_reference(self,
                             "solid_component_list_ref",
                             self.config.solid_phase_config.property_package.
                             component_list)
        add_object_reference(self,
                             "gas_phase_list_ref",
                             self.config.gas_phase_config.property_package.
                             phase_list)
        add_object_reference(self,
                             "solid_phase_list_ref",
                             self.config.solid_phase_config.property_package.
                             phase_list)

        if self.config.gas_phase_config.reaction_package is not None:
            add_object_reference(self,
                                 "homogeneous_rate_reaction_idx_ref",
                                 self.config.gas_phase_config.
                                 reaction_package.rate_reaction_idx)

        if self.config.solid_phase_config.reaction_package is not None:
            add_object_reference(self,
                                 "heterogeneous_rate_reaction_idx_ref",
                                 self.config.solid_phase_config.
                                 reaction_package.rate_reaction_idx)

        # Model parameters
        self.pi = Param(initialize=math.pi, doc="pi")
        self.bed_voidage = Param(initialize=0.4, doc="Bed voidage [-]")

        self.smooth_abs_param = Param(default=1e-8, mutable=True)

        # Unit Model variables
        self.bed_diameter = Var(domain=Reals, initialize=1,
                                doc='Reactor diameter [m]')
        self.bed_area = Var(domain=Reals, initialize=1,
                            doc='Reactor cross-sectional area [m2]')
        self.bed_height = Var(domain=Reals, initialize=1,
                              doc='Bed length [m]')

        # Phase specific variables
        self.velocity_superficial_gas = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals, initialize=0.05,
                doc='Gas superficial velocity [m/s]')
        self.velocity_superficial_solid = Var(
                self.flowsheet().config.time,
                domain=Reals, initialize=0.005,
                doc='Solid superficial velocity [m/s]')

        # Dimensionless numbers, mass and heat transfer coefficients
        self.Re_particle = Var(self.flowsheet().config.time,
                               self.length_domain,
                               domain=Reals, initialize=1.0,
                               doc='Particle Reynolds number [-]')

        self.Pr = Var(self.flowsheet().config.time,
                      self.length_domain,
                      domain=Reals, initialize=1.0,
                      doc='Prandtl number of gas in bed [-]')

        self.Nu_particle = Var(self.flowsheet().config.time,
                               self.length_domain,
                               domain=Reals, initialize=1.0,
                               doc='Particle Nusselt number [-]')
        self.gas_solid_htc = Var(self.flowsheet().config.time,
                                 self.length_domain,
                                 domain=Reals, initialize=1.0,
                                 doc='Gas-solid heat transfer coefficient'
                                 '[kJ/(m2Ks)]')

    # =========================================================================
        # Add performance equations

        # ---------------------------------------------------------------------
        # Geometry contraints

        # Bed area
        @self.Constraint(doc="Bed area")
        def bed_area_eqn(b):
            return b.bed_area == (
                    b.pi*(0.5*b.bed_diameter)**2)

        # Area of gas side, and solid side
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Gas side area")
        def gas_phase_config_area(b, t, x):
            return (b.gas_phase.area[t, x] ==
                    b.bed_area*b.bed_voidage)

        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Solid side area")
        def solid_phase_config_area(b, t, x):
            return (b.solid_phase.area[t, x] ==
                    b.bed_area*(1-b.bed_voidage))

        # Length of gas side, and solid side
        @self.Constraint(doc="Gas side length")
        def gas_phase_config_length(b):
            return (b.gas_phase.length == b.bed_height)

        @self.Constraint(doc="Solid side length")
        def solid_phase_config_length(b):
            return (b.solid_phase.length == b.bed_height)

        # ---------------------------------------------------------------------
        # Hydrodynamic contraints

        # Gas superficial velocity
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Gas superficial velocity")
        def gas_super_vel(b, t, x):
            return (b.velocity_superficial_gas[t, x] * b.bed_area *
                    b.gas_phase.properties[t, x].dens_mole_vap ==
                    b.gas_phase.properties[t, x].flow_mol)

        # Solid superficial velocity
        @self.Constraint(self.flowsheet().config.time,
                         doc="Solid superficial velocity")
        def solid_super_vel(b, t):
            return (b.velocity_superficial_solid[t] * b.bed_area *
                    b.solid_phase.properties[t, 1].dens_mass_sol ==
                    b.solid_phase.properties[t, 1].flow_mass)

        # Gas side pressure drop calculation
        if (self.config.gas_phase_config.has_pressure_change and
            self.config.gas_phase_config.pressure_drop_type ==
                "simple_correlation"):
            # Simplified pressure drop
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Gas side pressure drop calculation -"
                                 "simplified pressure drop")
            def gas_phase_config_pressure_drop(b, t, x):
                return b.gas_phase.deltaP[t, x]*1e5 == -0.2*(
                        b.velocity_superficial_gas[t, x] *
                        (b.solid_phase.properties[t, x].dens_mass_sol -
                         b.gas_phase.properties[t, x].dens_mass_vap))
        elif (self.config.gas_phase_config.has_pressure_change and
              self.config.gas_phase_config.pressure_drop_type ==
              "ergun_correlation"):
            # Ergun equation
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Gas side pressure drop calculation -"
                                 "ergun equation")
            def gas_phase_config_pressure_drop(b, t, x):
                return (1e2*-b.gas_phase.deltaP[t, x]*1e5 ==
                        1e2*(
                        150*(1 - b.bed_voidage) ** 2 *
                        b.gas_phase.properties[t, x].visc_d *
                        (b.velocity_superficial_gas[t, x] +
                         b.velocity_superficial_solid[t]) /
                        (b.solid_phase.properties[t, x].
                         _params.particle_dia ** 2 * b.bed_voidage ** 3)) +
                        1e2*(
                        1.75*b.gas_phase.properties[t, x].dens_mass_vap *
                        (1 - b.bed_voidage) *
                        (b.velocity_superficial_gas[t, x] +
                         b.velocity_superficial_solid[t]) ** 2 /
                        (b.solid_phase.properties[t, x]._params.particle_dia *
                         b.bed_voidage**3)))
            # The above expression has no absolute values - assumes:
            # (velocity_superficial_gas + velocity_superficial_solid) > 0

        # ---------------------------------------------------------------------
        # Reaction contraints

        # Build homogeneous reaction constraints
        if self.config.gas_phase_config.reaction_package is not None:
            # Gas side rate reaction extent
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             self.homogeneous_rate_reaction_idx_ref,
                             doc="Gas side rate reaction extent")
            def gas_phase_config_rxn_ext(b, t, x, r):
                return 1e3*b.gas_phase.rate_reaction_extent[t, x, r] == 1e3*(
                        b.gas_phase.reactions[t, x].reaction_rate[r] *
                        b.gas_phase.area[t, x])

        # Build hetereogeneous reaction constraints
        if self.config.solid_phase_config.reaction_package is not None:
            # Solid side rate reaction extent
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             self.heterogeneous_rate_reaction_idx_ref,
                             doc="Solid side rate reaction extent")
            def solid_phase_config_rxn_ext(b, t, x, r):
                return 1e3*b.solid_phase.rate_reaction_extent[t, x, r] == 1e3*(
                        b.solid_phase.reactions[t, x].reaction_rate[r] *
                        b.solid_phase.area[t, x])

            # Gas side heterogeneous rate reaction generation
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             self.gas_phase_list_ref,
                             self.gas_component_list_ref,
                             doc='Gas side heterogeneous'
                             'rate reaction generation')
            def gas_comp_hetero_rxn(b, t, x, p, j):
                return 1e3*b.gas_phase.mass_transfer_term[t, x, p, j] == 1e3*(
                        sum(b.solid_phase.reactions[t, x].
                            rate_reaction_stoichiometry[r, p, j] *
                            b.solid_phase.reactions[t, x].reaction_rate[r]
                            for r in b.heterogeneous_rate_reaction_idx_ref) *
                        b.solid_phase.area[t, x])

        # ---------------------------------------------------------------------
        # Dimensionless numbers, mass and heat transfer coefficients

        # Particle Reynolds number
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Particle Reynolds number")
        def reynolds_number_particle(b, t, x):
            return (b.Re_particle[t, x] *
                    b.gas_phase.properties[t, x].visc_d ==
                    b.velocity_superficial_gas[t, x] *
                    b.solid_phase.properties[t, x]._params.particle_dia *
                    b.gas_phase.properties[t, x].dens_mass_vap)

        # Prandtl number
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Prandtl number of gas in bed")
        def prandtl_number(b, t, x):
            return (b.Pr[t, x] *
                    b.gas_phase.properties[t, x].therm_cond ==
                    b.solid_phase.properties[t, x].cp_mass *
                    b.gas_phase.properties[t, x].visc_d)

        # Particle Nusselt number
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Particle Nusselt number")
        def nusselt_number_particle(b, t, x):
            return (b.Nu_particle[t, x] ** 3 ==
                    ((2.0 + 1.1 * ((b.Re_particle[t, x] ** 2 +
                                    b.smooth_abs_param) ** 0.5) ** 0.6) ** 3) *
                    b.Pr[t, x])

        # Gas-solid heat transfer coefficient
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Gas-solid heat transfer coefficient")
        def gas_solid_htc_eqn(b, t, x):
            return (1e-3*b.gas_solid_htc[t, x] *
                    b.solid_phase.properties[t, x].
                    _params.particle_dia == 1e-3 *
                    b.Nu_particle[t, x] *
                    b.gas_phase.properties[t, x].therm_cond)

        if self.config.energy_balance_type != EnergyBalanceType.none:
            # Gas phase - gas to solid heat transfer
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Gas phase - gas to solid heat transfer")
            def gas_phase_heat_transfer(b, t, x):
                return (b.gas_phase.heat[t, x] *
                        b.solid_phase.properties[t, x]._params.particle_dia ==
                        -6 * b.gas_solid_htc[t, x] *
                        (b.gas_phase.properties[t, x].temperature -
                         b.solid_phase.properties[t, x].temperature) *
                        b.solid_phase.area[t, x])

            # Solid phase - gas to solid heat transfer
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Solid phase - gas to solid heat transfer")
            def solid_phase_heat_transfer(b, t, x):
                return (b.solid_phase.heat[t, x] *
                        b.solid_phase.properties[t, x]._params.particle_dia ==
                        6 * b.gas_solid_htc[t, x] *
                        (b.gas_phase.properties[t, x].temperature -
                         b.solid_phase.properties[t, x].temperature) *
                        b.solid_phase.area[t, x])

    # =========================================================================
    # Model initialization routine

    def initialize(blk, gas_phase_state_args={}, solid_phase_state_args={},
                   outlvl=idaeslog.NOTSET,
                   solver='ipopt', optarg={'tol': 1e-6}):
        """
        Initialisation routine for MB unit (default solver ipopt).

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialisation routine
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        """
        # Set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Set solver options
        opt = SolverFactory(solver)
        opt.options = optarg

        # Keep all unit model geometry constraints, derivative_var constraints,
        # and property block constraints active
        # Exceptions:
        # In control volumes - keep conservation linking constraints and
        # holdup calculation (for dynamic flowsheets) constraints active

        geometry_constraints_list = ["bed_area_eqn", "solid_phase_config_area",
                                     "gas_phase_config_area",
                                     "gas_phase_config_length",
                                     "solid_phase_config_length"]

        def constraints_generator(obj):
            for c in obj.component_objects(Constraint, descend_into=True):
                if not c.local_name.endswith("_disc_eq") \
                    and \
                    "_linking_constraint" not in c.local_name \
                    and not \
                    c.parent_block().local_name.startswith("properties") \
                    and \
                        "_holdup_calculation" not in c.local_name:
                    yield c

        for c in constraints_generator(blk):
            if (c.local_name not in geometry_constraints_list):
                c.deactivate()

        # ---------------------------------------------------------------------

        print()
        print('Initialize thermophysical properties')
        # Initialize gas_phase block
        gas_phase_flags = blk.gas_phase.properties.initialize(
                                state_args=gas_phase_state_args,
                                outlvl=outlvl,
                                optarg=optarg,
                                solver=solver,
                                hold_state=True)

        # Initialize solid_phase properties block
        solid_phase_flags = blk.solid_phase.properties.initialize(
                                state_args=solid_phase_state_args,
                                outlvl=outlvl,
                                optarg=optarg,
                                solver=solver,
                                hold_state=True)

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Initialize hydrodynamics (velocities)
        for t in blk.flowsheet().config.time:
            calculate_variable_from_constraint(
                blk.velocity_superficial_solid[t],
                blk.solid_super_vel[t])
            for x in blk.length_domain:
                calculate_variable_from_constraint(
                    blk.velocity_superficial_gas[t, x],
                    blk.gas_super_vel[t, x])

        blk.gas_super_vel.activate()
        blk.solid_super_vel.activate()

        print()
        print('Initialize hydrodynamics')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
        if results.solver.termination_condition \
                == TerminationCondition.optimal:
            init_log.info_high(
                "Initialization Step 2 {}.".format(
                        idaeslog.condition(results))
                        )
        else:
            _log.warning('{} Initialisation Step 2 Failed.'
                         .format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize mass balance - no reaction and no pressure drop

        # Unfix material balance state variables but keep other states fixed
        blk.gas_phase.properties.release_state(
                gas_phase_flags)
        blk.solid_phase.properties.release_state(
                solid_phase_flags)
        for t in blk.flowsheet().config.time:
            for x in blk.length_domain:
                blk.gas_phase.properties[t, x].pressure.fix()
                blk.gas_phase.properties[t, x].temperature.fix()
                blk.solid_phase.properties[t, x].temperature.fix()

        blk.gas_phase.material_balances.activate()

        if blk.config.gas_phase_config.reaction_package is not None:
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    for p in blk.gas_phase_list_ref:
                        for j in blk.gas_component_list_ref:
                            (blk.gas_phase.
                             rate_reaction_generation[t, x, p, j].fix(0.0))

        blk.solid_phase.material_balances.activate()

        if blk.config.solid_phase_config.reaction_package is not None:
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    for p in blk.solid_phase_list_ref:
                        for j in blk.solid_component_list_ref:
                            (blk.solid_phase.
                             rate_reaction_generation[t, x, p, j].fix(0.0))

                # Gas side heterogeneous rate reaction generation
                for x in blk.length_domain:
                    for p in blk.gas_phase_list_ref:
                        for j in blk.gas_component_list_ref:
                            (blk.gas_phase.
                             mass_transfer_term[t, x, p, j].fix(0.0))

        print()
        print('Mass balance - no reaction and no pressure drop')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
        if results.solver.termination_condition \
                == TerminationCondition.optimal:
            init_log.info_high(
                "Initialization Step 3a {}.".format(
                        idaeslog.condition(results))
                        )
        else:
            _log.warning('{} Initialisation Step 3a Failed.'
                         .format(blk.name))

        # Initialize mass balance - with reaction and no pressure drop
        if blk.config.gas_phase_config.reaction_package is not None:
            blk.gas_phase.rate_reaction_stoichiometry_constraint.activate()
            blk.gas_phase_config_rxn_ext.activate()

            # Initialize reaction property package
            blk.gas_phase.reactions.activate()
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    obj = blk.gas_phase.reactions[t, x]
                    for c in obj.component_objects(
                            Constraint, descend_into=False):
                        c.activate()

            blk.gas_phase.reactions.initialize(outlvl=outlvl,
                                               optarg=optarg,
                                               solver=solver)

            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    for r in blk.homogeneous_rate_reaction_idx_ref:
                        calculate_variable_from_constraint(
                                blk.gas_phase.rate_reaction_extent[t, x, r],
                                blk.gas_phase_config_rxn_ext[t, x, r])
                    for p in blk.gas_phase_list_ref:
                        for j in blk.gas_component_list_ref:
                            (blk.gas_phase.
                             rate_reaction_generation[t, x, p, j].unfix())
                            calculate_variable_from_constraint(
                                blk.gas_phase.
                                rate_reaction_generation[t, x, p, j],
                                blk.gas_phase.
                                rate_reaction_stoichiometry_constraint
                                [t, x, p, j])

        if blk.config.solid_phase_config.reaction_package is not None:
            blk.gas_comp_hetero_rxn.activate()
            blk.solid_phase.rate_reaction_stoichiometry_constraint.activate()
            blk.solid_phase_config_rxn_ext.activate()

            # Initialize reaction property package
            blk.solid_phase.reactions.activate()
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    obj = blk.solid_phase.reactions[t, x]
                    for c in obj.component_objects(
                            Constraint, descend_into=False):
                        c.activate()

            blk.solid_phase.reactions.initialize(outlvl=outlvl,
                                                 optarg=optarg,
                                                 solver=solver)

            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    for p in blk.gas_phase_list_ref:
                        for j in blk.gas_component_list_ref:
                            (blk.gas_phase.
                             mass_transfer_term[t, x, p, j].unfix())
                            calculate_variable_from_constraint(
                                blk.gas_phase.mass_transfer_term[t, x, p, j],
                                blk.gas_comp_hetero_rxn[t, x, p, j])
                for x in blk.length_domain:
                    for r in blk.heterogeneous_rate_reaction_idx_ref:
                        calculate_variable_from_constraint(
                            blk.solid_phase.rate_reaction_extent[t, x, r],
                            blk.solid_phase_config_rxn_ext[t, x, r])
                    for p in blk.solid_phase_list_ref:
                        for j in blk.solid_component_list_ref:
                            (blk.solid_phase.
                             rate_reaction_generation[t, x, p, j].unfix())
                            calculate_variable_from_constraint(
                                blk.solid_phase.
                                rate_reaction_generation[t, x, p, j],
                                blk.solid_phase.
                                rate_reaction_stoichiometry_constraint
                                [t, x, p, j])

        if (blk.config.gas_phase_config.reaction_package is not None or
                blk.config.solid_phase_config.reaction_package is not None):
            print()
            print('Mass balance - with reaction and no pressure drop')
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                init_log.info_high(
                    "Initialization Step 3b {}.".format(
                            idaeslog.condition(results))
                            )
            else:
                _log.warning('{} Initialisation Step 3b Failed.'
                             .format(blk.name))

        # Initialize mass balance - with pressure drop
        for t in blk.flowsheet().config.time:
            for x in blk.length_domain:
                # Unfix all pressure variables except at the inlet
                if (blk.gas_phase.properties[t, x].
                        config.defined_state) is False:
                    blk.gas_phase.properties[t, x].pressure.unfix()

        blk.gas_phase.pressure_balance.activate()

        # Set scaling factors for pressure balance equation
        blk.gas_phase.scaling_factor_pressure = 1e2

        if blk.config.gas_phase_config.has_pressure_change:
            blk.gas_phase_config_pressure_drop.activate()

            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    calculate_variable_from_constraint(
                        blk.gas_phase.deltaP[t, x],
                        blk.gas_phase_config_pressure_drop[t, x])

            print()
            print('Mass balance - with pressure drop')
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                init_log.info_high(
                    "Initialization Step 3c {}.".format(
                            idaeslog.condition(results))
                            )
            else:
                _log.warning('{} Initialisation Step 3c Failed.'
                             .format(blk.name))
        # ---------------------------------------------------------------------
        # Initialize energy balance
        if blk.config.energy_balance_type != EnergyBalanceType.none:
            # Initialize dimensionless numbers,
            # mass and heat transfer coefficients
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    calculate_variable_from_constraint(
                        blk.Re_particle[t, x],
                        blk.reynolds_number_particle[t, x])
                    calculate_variable_from_constraint(
                        blk.Pr[t, x],
                        blk.prandtl_number[t, x])
                    calculate_variable_from_constraint(
                        blk.Nu_particle[t, x],
                        blk.nusselt_number_particle[t, x])
                    calculate_variable_from_constraint(
                            blk.gas_solid_htc[t, x],
                            blk.gas_solid_htc_eqn[t, x])
                    calculate_variable_from_constraint(
                        blk.gas_phase.heat[t, x],
                        blk.gas_phase_heat_transfer[t, x])
                    calculate_variable_from_constraint(
                        blk.solid_phase.heat[t, x],
                        blk.solid_phase_heat_transfer[t, x])

            # Unfix temperatures
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    # Unfix all gas temperature variables except at the inlet
                    if (blk.gas_phase.properties[t, x].
                            config.defined_state) is False:
                        blk.gas_phase.properties[t, x].temperature.unfix()
                for x in blk.length_domain:
                    # Unfix all solid temperature variables except at the inlet
                    if (blk.solid_phase.properties[t, x].
                            config.defined_state) is False:
                        blk.solid_phase.properties[t, x].temperature.unfix()

            blk.reynolds_number_particle.activate()
            blk.prandtl_number.activate()
            blk.nusselt_number_particle.activate()
            blk.gas_solid_htc_eqn.activate()

            # Activate gas phase energy balance equations
            blk.gas_phase_heat_transfer.activate()
            blk.gas_phase.enthalpy_balances.activate()

            # Activate solid phase energy balance equations
            blk.solid_phase_heat_transfer.activate()
            blk.solid_phase.enthalpy_balances.activate()

            print()
            print('Energy balance')
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                init_log.info_high(
                    "Initialization Step 4 {}.".format(
                            idaeslog.condition(results))
                            )
            else:
                _log.warning('{} Initialisation Step 4 Failed.'
                             .format(blk.name))

    def _results_plot(blk):
        '''
        Plot method for common moving bed variables

        Variables plotted:
            Tg : Temperature in gas phase
            Ts : Temperature in solid phase
            vg : Superficial gas velocity
            P : Pressure in gas phase
            Ftotal : Total molar flowrate of gas
            Mtotal : Total mass flowrate of solid
            Cg : Concentration of gas components in the gas phase
            y_frac : Mole fraction of gas components in the gas phase
            x_frac : Mass fraction of solid components in the solid phase
        '''
        print()
        print('================================= Reactor plots ==============='
              '==================')
        Tg = []
        Ts = []
        P = []
        Ftotal = []
        Mtotal = []
        vg = []

        for t in blk.flowsheet().config.time:
            for x in blk.gas_phase.length_domain:
                vg.append(value(blk.velocity_superficial_gas[t, x]))

        fig_vg = plt.figure(1)
        plt.plot(blk.gas_phase.length_domain, vg, label='vg')
        plt.legend(loc=0, ncol=2)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Superficial gas velocity [m/s]")
        fig_vg.savefig('superficial_vel.png')

        # Pressure profile
        for t in blk.flowsheet().config.time:
            for x in blk.gas_phase.length_domain:
                P.append(blk.gas_phase.properties[t, x].pressure.value)

        fig_P = plt.figure(2)
        plt.plot(blk.gas_phase.length_domain, P, label='P')
        plt.legend(loc=0, ncol=2)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Total Pressure [bar]")
        fig_P.savefig('Pressure.png')

        # Temperature profile
        for t in blk.flowsheet().config.time:
            for x in blk.gas_phase.length_domain:
                Tg.append(blk.gas_phase.properties[t, x].temperature.value)
            for x in blk.solid_phase.length_domain:
                Ts.append(blk.solid_phase.properties[t, x].temperature.value)
        fig_T = plt.figure(3)
        plt.plot(blk.gas_phase.length_domain, Tg, label='Tg')
        plt.plot(blk.solid_phase.length_domain, Ts, label='Ts')
        plt.legend(loc=0, ncol=2)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Temperature [K]")
        fig_T.savefig('Temperature.png')

        # Bulk gas phase total molar flow rate
        for t in blk.flowsheet().config.time:
            for x in blk.gas_phase.length_domain:
                Ftotal.append(blk.gas_phase.properties[t, x].flow_mol.value)
        fig_Ftotal = plt.figure(4)
        plt.plot(blk.gas_phase.length_domain, Ftotal)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Total molar gas flow rate [mol/s]")
        fig_Ftotal.savefig('Total_gas_flow.png')

        # Bulk solid phase total mass flow rate
        for t in blk.flowsheet().config.time:
            for x in blk.solid_phase.length_domain:
                Mtotal.append(blk.solid_phase.properties[t, x].flow_mass.value)
        fig_Mtotal = plt.figure(5)
        plt.plot(blk.solid_phase.length_domain, Mtotal)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Solid total mass flow rate [kg/s]")
        fig_Mtotal.savefig('Total_solid_flow.png')

        # Gas phase mole fractions
        for t in blk.flowsheet().config.time:
            for j in blk.gas_component_list_ref:
                y_frac = []
                for x in blk.gas_phase.length_domain:
                    y_frac.append(
                        blk.gas_phase.properties[t, x].mole_frac[j].value)
                fig_y = plt.figure(6)
                plt.plot(blk.gas_phase.length_domain, y_frac, label=j)
        plt.legend(loc=0, ncol=len(blk.gas_component_list_ref))
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("y_frac [-]")
        fig_y.savefig('Gas_mole_fractions.png')

        # Solid phase mass fractions
        for t in blk.flowsheet().config.time:
            for j in blk.solid_component_list_ref:
                x_frac = []
                for x in blk.solid_phase.length_domain:
                    x_frac.append(
                        blk.solid_phase.properties[t, x].mass_frac[j].value)
                fig_x = plt.figure(7)
                plt.plot(blk.solid_phase.length_domain, x_frac, label=j)
        plt.legend(loc=0, ncol=len(blk.solid_component_list_ref))
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("x_frac [-]")
        fig_x.savefig('Solid_mass_fractions.png')

        # Gas phase concentrations
        for t in blk.flowsheet().config.time:
            for j in blk.gas_component_list_ref:
                Cg = []
                for x in blk.gas_phase.length_domain:
                    Cg.append(
                        blk.gas_phase.properties[t, x].
                        dens_mole_comp_vap[j].value)
                fig_Cg = plt.figure(8)
                plt.plot(blk.gas_phase.length_domain, Cg, label=j)
        plt.legend(loc=0, ncol=len(blk.gas_component_list_ref))
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Concentration [mol/m3]")
        fig_Cg.savefig('Gas_concentration.png')

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Gas Inlet": self.gas_inlet,
                "Gas Outlet": self.gas_outlet,
                "Solid Inlet": self.solid_inlet,
                "Solid Outlet": self.solid_outlet,
            },
            time_point=time_point,
            )
