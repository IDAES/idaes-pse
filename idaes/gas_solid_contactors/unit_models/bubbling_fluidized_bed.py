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
IDAES Bubbling Fluidized Bed Model.
The 2-region bubbling fluidized bed model is a 1D axially discretized model
with two phases (gas and solid), and two regions (bubble and emulsion)
resulting in 3 control volume_1D blocks (bubble, gas_emulsion solid_emulsion).
The model captures the gas-solid interaction between both phases and regions
through reaction, mass and heat transfer.
Equations written in this model were derived from:
A. Lee, D.C. Miller. A one-dimensional (1-D) three-region model for a bubbling
fluidized-bed Adsorber, Ind. Eng. Chem. Res. 52 (2013) 469–484.
Assumptions:
Property package contains temperature and pressure variables
Property package contains minimum fluidization velocity and voidage parameters
Gas emulsion is at minimum fluidization conditions
Gas feeds into emulsion region before the excess enters into the bubble region
Solid superficial velocity is constant throughout the bed
"""

# Import Python libraries
import matplotlib.pyplot as plt

# Import Pyomo libraries
from pyomo.environ import (Var, Param, Reals, SolverFactory,
                           TerminationCondition, Constraint,
                           TransformationFactory, sqrt, value)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.dae import ContinuousSet, DerivativeVar

# Import IDAES cores
from idaes.core import (ControlVolume1DBlock, UnitModelBlockData,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        FlowDirection)
from idaes.core.util.config import (is_physical_parameter_block,
                                    is_reaction_parameter_block)
from idaes.core.util.exceptions import (ConfigurationError)
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.control_volume1d import DistributedVars
from idaes.core.util.constants import Constants as constants
from idaes.core.util.math import smooth_min, smooth_max
import idaes.logger as idaeslog
from idaes.core.util import get_solver

__author__ = "Chinedu Okoli"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("BubblingFluidizedBed")
class BubblingFluidizedBedData(UnitModelBlockData):
    """
    Standard Bubbling Fluidized Bed Unit Model Class
    """

    # Create template for unit level config arguments
    CONFIG = UnitModelBlockData.CONFIG()

    # Unit level config arguments
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
provided (default = [0.0, 1.0])."""))
    CONFIG.declare("transformation_method", ConfigValue(
        default="dae.finite_difference",
        domain=In(["dae.finite_difference", "dae.collocation"]),
        description="Method to use for DAE transformation",
        doc="""Method to use to transform domain. Must be a method recognised
by the Pyomo TransformationFactory,
**default** - "dae.finite_difference".
**Valid values:** {
**"dae.finite_difference"** - Use a finite difference transformation scheme,
**"dae.collocation"** - use a collocation transformation scheme}"""))
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
        default="co_current",
        domain=In(['co_current', 'counter_current']),
        description="Flow configuration of Bubbling Fluidized Bed",
        doc="""Flow configuration of Bubbling Fluidized Bed
**default** - "co_current".
**Valid values:** {
**"co_current"** - gas flows from 0 to 1, solid flows from 0 to 1,
**"counter_current"** -  gas flows from 0 to 1, solid flows from 1 to 0.}"""))
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
    CONFIG.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.none.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))

    # Create template for phase specific config arguments
    _PhaseTemplate = UnitModelBlockData.CONFIG()
    _PhaseTemplate.declare("has_equilibrium_reactions", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Equilibrium reaction construction flag",
        doc="""Indicates whether terms for equilibrium controlled reactions
should be constructed,
**default** - True.
**Valid values:** {
**True** - include equilibrium reaction terms,
**False** - exclude equilibrium reaction terms.}"""))
    _PhaseTemplate.declare("property_package", ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object"""))
    _PhaseTemplate.declare("property_package_args", ConfigValue(
        default={},
        domain=dict,
        description="Arguments for constructing gas property package",
        doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)"""))
    _PhaseTemplate.declare("reaction_package", ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package to use for control volume",
        doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}"""))
    _PhaseTemplate.declare("reaction_package_args", ConfigBlock(
        implicit=True,
        implicit_domain=ConfigBlock,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}"""))

    # Create individual config blocks for the gas and solid phases
    CONFIG.declare("gas_phase_config",
                   _PhaseTemplate(doc="gas phase config arguments"
                                  ))
    CONFIG.declare("solid_phase_config",
                   _PhaseTemplate(doc="solid phase config arguments"
                                  ))

    # =========================================================================
    def build(self):
        """
        Begin building model
        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to build default attributes
        super(BubblingFluidizedBedData, self).build()

        # Consistency check for transformation method and transformation scheme
        if (self.config.transformation_method == "dae.finite_difference" and
                self.config.transformation_scheme is None):
            self.config.transformation_scheme = "BACKWARD"
        elif (self.config.transformation_method == "dae.collocation" and
                self.config.transformation_scheme is None):
            self.config.transformation_scheme = "LAGRANGE-RADAU"
        elif (self.config.transformation_method == "dae.finite_difference" and
              self.config.transformation_scheme != "BACKWARD" and
              self.config.transformation_scheme != "FORWARD"):
            raise ConfigurationError("{} invalid value for "
                                     "transformation_scheme argument. "
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

        # Set flow directions for the control volume blocks
        # Gas flows from 0 to 1, solid flows from 0 to 1
        if self.config.flow_type == "co_current":
            set_direction_gas = FlowDirection.forward
            set_direction_solid = FlowDirection.forward
            # Gas flows from 0 to 1, solid flows from 1 to 0
        if self.config.flow_type == "counter_current":
            set_direction_gas = FlowDirection.forward
            set_direction_solid = FlowDirection.backward

        # Set arguments for gas phase if homogeneous reaction block
        if self.config.gas_phase_config.reaction_package is not None:
            has_rate_reaction_gas = True
        else:
            has_rate_reaction_gas = False

        # Set arguments for emulsion region if heterogeneous reaction block
        if self.config.solid_phase_config.reaction_package is not None:
            has_rate_reaction_solid = True
        else:
            has_rate_reaction_solid = False

        # Set heat transfer terms
        if self.config.energy_balance_type != EnergyBalanceType.none:
            has_heat_transfer = True
        else:
            has_heat_transfer = False

        # Set heat of reaction terms
        if (self.config.energy_balance_type != EnergyBalanceType.none
                and self.config.gas_phase_config.reaction_package
                is not None):
            has_heat_of_reaction_gas = True
        else:
            has_heat_of_reaction_gas = False

        if (self.config.energy_balance_type != EnergyBalanceType.none
                and self.config.solid_phase_config.reaction_package
                is not None):
            has_heat_of_reaction_solid = True
        else:
            has_heat_of_reaction_solid = False

        # Create a unit model length domain
        self.length_domain = ContinuousSet(
                                bounds=(0.0, 1.0),
                                initialize=self.config.length_domain_set,
                                doc="Normalized length domain")

    # =========================================================================
        """ Build Control volume 1D for the bubble region and
            populate its control volume"""

        self.bubble = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": DistributedVars.variant,
            "property_package": self.config.gas_phase_config.property_package,
            "property_package_args":
                self.config.gas_phase_config.property_package_args,
            "reaction_package": self.config.gas_phase_config.reaction_package,
            "reaction_package_args":
                self.config.gas_phase_config.reaction_package_args})

        self.bubble.add_geometry(
                length_domain=self.length_domain,
                length_domain_set=self.config.length_domain_set,
                flow_direction=set_direction_gas)

        self.bubble.add_state_blocks(
            information_flow=set_direction_gas,
            has_phase_equilibrium=False)

        if self.config.gas_phase_config.reaction_package is not None:
            self.bubble.add_reaction_blocks(
                    has_equilibrium=(
                        self.config.gas_phase_config.has_equilibrium_reactions)
                    )

        self.bubble.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
            has_mass_transfer=True,
            has_rate_reactions=has_rate_reaction_gas)

        self.bubble.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=has_heat_transfer,
            has_heat_of_reaction=has_heat_of_reaction_gas)

        self.bubble.add_momentum_balances(
            balance_type=MomentumBalanceType.none,
            has_pressure_change=False)

    # =========================================================================
        """ Build Control volume 1D for the gas_emulsion region and
            populate its control volume"""

        self.gas_emulsion = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": DistributedVars.variant,
            "property_package": self.config.gas_phase_config.property_package,
            "property_package_args":
                self.config.gas_phase_config.property_package_args,
            "reaction_package":
                self.config.gas_phase_config.reaction_package,
            "reaction_package_args":
                self.config.gas_phase_config.reaction_package_args})

        self.gas_emulsion.add_geometry(
                length_domain=self.length_domain,
                length_domain_set=self.config.
                length_domain_set,
                flow_direction=set_direction_gas)

        self.gas_emulsion.add_state_blocks(
            information_flow=set_direction_gas,
            has_phase_equilibrium=False)

        if self.config.gas_phase_config.reaction_package is not None:
            self.gas_emulsion.add_reaction_blocks(
                    has_equilibrium=(
                        self.config.gas_phase_config.has_equilibrium_reactions)
                    )

        self.gas_emulsion.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
            has_mass_transfer=True,
            has_rate_reactions=has_rate_reaction_gas)

        self.gas_emulsion.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=has_heat_transfer,
            has_heat_of_reaction=has_heat_of_reaction_gas)

        self.gas_emulsion.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

    # =========================================================================
        """ Build Control volume 1D for solid emulsion region and
            populate solid control volume"""

        self.solid_emulsion = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": DistributedVars.variant,
            "property_package":
                self.config.solid_phase_config.property_package,
            "property_package_args":
                self.config.solid_phase_config.property_package_args,
            "reaction_package":
                self.config.solid_phase_config.reaction_package,
            "reaction_package_args":
                self.config.solid_phase_config.reaction_package_args})

        self.solid_emulsion.add_geometry(
            length_domain=self.length_domain,
            length_domain_set=self.config.
            length_domain_set,
            flow_direction=set_direction_solid)

        self.solid_emulsion.add_state_blocks(
            information_flow=set_direction_solid,
            has_phase_equilibrium=False)

        if self.config.solid_phase_config.reaction_package is not None:
            # TODO - a generalization of the heterogeneous reaction block
            # The heterogeneous reaction block does not use the
            # add_reaction_blocks in control volumes as control volumes are
            # currently setup to handle only homogeneous reaction properties.
            # Thus appending the heterogeneous reaction block to the
            # solid state block is currently hard coded here.

            tmp_dict = dict(
                **self.config.solid_phase_config.reaction_package_args)
            tmp_dict["gas_state_block"] = self.gas_emulsion.properties
            tmp_dict["solid_state_block"] = (
                self.solid_emulsion.properties)
            tmp_dict["has_equilibrium"] = (
                self.config.solid_phase_config.has_equilibrium_reactions)
            tmp_dict["parameters"] = (
                self.config.solid_phase_config.reaction_package)
            self.solid_emulsion.reactions = (
                self.config.solid_phase_config.reaction_package.
                reaction_block_class(
                    self.flowsheet().config.time,
                    self.length_domain,
                    doc="Reaction properties in control volume",
                    default=tmp_dict))

        self.solid_emulsion.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
            has_mass_transfer=False,
            has_rate_reactions=has_rate_reaction_solid)

        self.solid_emulsion.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=has_heat_transfer,
            has_heat_of_reaction=has_heat_of_reaction_solid)

        self.solid_emulsion.add_momentum_balances(
            balance_type=MomentumBalanceType.none,
            has_pressure_change=False)

    # =========================================================================
        """ Add Ports for gas and solid inlets and outlets"""
        # Build inlet and outlet state blocks for model to be attached to ports
        # Extra state blocks are included at the inlet and outlets to model
        # the mass and energy balances (splitting and mixing) that takes place
        # between the gas feed/product from the reactor, and the bubble and
        # gas_emulsion regions at the reactor boundaries. It is also used to
        # capture the different solid phase boundaries resulting from the
        # varying solid flow directions (co_current and counter_current)
        self.gas_inlet_block = (
            self.config.gas_phase_config.property_package.build_state_block(
                self.flowsheet().config.time,
                default={"defined_state": True}))

        self.gas_outlet_block = (
            self.config.gas_phase_config.property_package.build_state_block(
                self.flowsheet().config.time,
                default={"defined_state": False}))

        self.solid_inlet_block = (
            self.config.solid_phase_config.property_package.build_state_block(
                self.flowsheet().config.time,
                default={"defined_state": True}))

        self.solid_outlet_block = (
            self.config.solid_phase_config.property_package.build_state_block(
                self.flowsheet().config.time,
                default={"defined_state": False}))

        # Add Ports for gas side
        self.add_inlet_port(name="gas_inlet",
                            block=self.gas_inlet_block)
        self.add_outlet_port(name="gas_outlet",
                             block=self.gas_outlet_block)

        # Add Ports for solid side
        self.add_inlet_port(name="solid_inlet",
                            block=self.solid_inlet_block)
        self.add_outlet_port(name="solid_outlet",
                             block=self.solid_outlet_block)

    # =========================================================================
        """ Apply transformation and add performace equation method"""
        self._make_vars_params()
        self._apply_transformation()
        self._make_performance()

    # =========================================================================
    def _make_vars_params(self):
        """
        Make model variables and parameters.
        Args:
            None
        Returns:
            None
        """

        # Declare Mutable Parameters
        self.eps = Param(mutable=True,
                         default=1e-8,
                         doc='Smoothing Factor for Smooth IF Statements')

        # Vessel dimensions
        self.bed_diameter = Var(domain=Reals,
                                initialize=1,
                                doc='Reactor Diameter [m]')
        self.bed_area = Var(domain=Reals,
                            initialize=1,
                            doc='Reactor Cross-sectional Area [m2]')
        self.bed_height = Var(domain=Reals,
                              initialize=1,
                              doc='Bed Height [m]')

        # Distributor Design
        self.area_orifice = Var(
                domain=Reals,
                initialize=1,
                doc='Distributor Plate Area per Orifice [m^2/orifice]')
        self.number_orifice = Var(
                domain=Reals,
                doc='Distributor Plate Orifices per Area [orifices/m^2]')

        # Velocities
        self.velocity_superficial_gas = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Gas Superficial Velocity [m/s]')
        self.velocity_bubble = Var(
                self.flowsheet().config.time,
                self.length_domain, domain=Reals,
                doc='Bubble Velocity [m/s]')
        self.velocity_emulsion_gas = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Emulsion Region Superficial Gas Velocity [m/s]')
        self.velocity_superficial_solid = Var(
                self.flowsheet().config.time,
                domain=Reals, initialize=0.01,
                doc='Solid superficial velocity [m/s]')

        # Bubble Dimensions and Hydrodynamics
        self.bubble_diameter = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                initialize=1,
                doc='Average Bubble Diameter [m]')
        self.delta = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                initialize=1,
                doc='Volume Fraction Occupied by Bubble Region [m^3/m^3]')
        self.delta_e = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                initialize=1,
                doc='Volume Fraction Occupied by Emulsion Region [m^3/m^3]')
        self.voidage_average = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                initialize=1,
                doc='Cross-Sectional Average Voidage Fraction [m^3/m^3]')
        self.voidage_emulsion = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                initialize=1,
                doc='Emulsion Region Voidage Fraction [m^3/m^3]')
        self.bubble_growth_coeff = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Bubble Growth Coefficient [-]')
        self.bubble_diameter_max = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Maximum Theoretical Bubble Diameter [m]')
        self.velocity_bubble_rise = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Bubble Rise Velocity [m/s]')
        self.average_gas_density = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='average gas density [mol/m3]')

        # Gas emulsion heterogeneous reaction variable
        self.gas_emulsion_hetero_rxn = Var(
                self.flowsheet().config.time,
                self.length_domain,
                self.config.gas_phase_config.property_package.component_list,
                domain=Reals,
                initialize=0.0,
                doc='Heterogeneous Rate Reaction'
                    'Generation in the Gas Emulsion')

        # Mass transfer coefficients
        self.Kbe = Var(
                self.flowsheet().config.time,
                self.length_domain,
                self.config.gas_phase_config.property_package.component_list,
                domain=Reals,
                initialize=1,
                doc='Bubble to Emulsion Gas Mass Transfer Coefficient [1/s]')
        self.Kgbulk_c = Var(
                self.flowsheet().config.time,
                self.length_domain,
                self.config.gas_phase_config.property_package.component_list,
                domain=Reals,
                initialize=1,
                doc='Gas Phase Component Bulk Transfer Rate [mol/m.s]')

        # Heat transfer coefficients
        self.Hbe = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Bubble to Emulsion Gas Heat Transfer Coefficient'
                    '[kJ/m^3.K.s]')
        self.Hgbulk = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Gas Phase Bulk Enthalpy Transfer Rate [kJ/m.s]')
        self.htc_conv = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Gas to Solid Energy Convective Heat Transfer'
                    'Coefficient [kJ/m^2.K.s]')

        # Heat transfer terms
        self.ht_conv = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Gas to Solid Convective Enthalpy Transfer in'
                    'Emulsion Region [kJ/m^2.K.s]')

        # Reformulation variables
        self._reform_var_1 = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                doc='Reformulation Variable in Bubble'
                    'Diameter Equation [reform var 1]')
        self._reform_var_2 = Var(
                self.flowsheet().config.time,
                self.length_domain,
                self.config.gas_phase_config.property_package.component_list,
                domain=Reals,
                initialize=1,
                doc='Bubble to Emulsion Gas Mass Transfer'
                    'Coefficient Reformulation Variable [reform var 2]')
        self._reform_var_3 = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                initialize=1,
                doc='Bubble to Emulsion Gas Mass Transfer'
                    'Coefficient Reformulation Variable [reform var 3]')
        self._reform_var_4 = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                initialize=1,
                doc='Bubble to Emulsion Gas Heat Transfer'
                    'Coefficient Reformulation Variable [reform var 4]')
        self._reform_var_5 = Var(
                self.flowsheet().config.time,
                self.length_domain,
                domain=Reals,
                initialize=1,
                doc='Convective Heat Transfer'
                    'Coefficient Reformulation Variable [reform var 5]')

        # Derivative variables
        self.ddia_bubbledx = DerivativeVar(
                self.bubble_diameter,
                wrt=self.length_domain,
                doc='Derivative of Bubble Diameter with Respect to Bed Height')

        # Fixed variables (these are parameters that can be estimated)
        self.Kd = Var(domain=Reals,
                      initialize=1,
                      doc='Bulk Gas Permeation Coefficient [m/s]')
        self.Kd.fix()
        self.deltaP_orifice = Var(domain=Reals,
                                  initialize=3.400,
                                  doc='Pressure Drop Across Orifice [bar]')
        self.deltaP_orifice.fix()

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

    # =========================================================================
    def _make_performance(self):
        # local aliases used to shorten object names
        gas_phase = self.config.gas_phase_config
        solid_phase = self.config.solid_phase_config

        # Add performance equations

        # ---------------------------------------------------------------------
        # Geometry contraints

        # Distributor design - Area of orifice
        @self.Constraint(doc="Area of Orifice")
        def orifice_area(b):
            return 1e1 * b.number_orifice*b.area_orifice == 1e1

        # Bed area
        @self.Constraint(doc="Bed Area")
        def bed_area_eqn(b):
            return b.bed_area == (
                    constants.pi*(0.5*b.bed_diameter)**2)

        # Area of bubble, gas_emulsion, solid_emulsion
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Cross-sectional Area Occupied by Bubbles")
        def bubble_area(b, t, x):
            return (b.bubble.area[t, x] ==
                    b.bed_area*b.delta[t, x])

        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Cross-sectional Area Occupied by Gas Emulsion")
        def gas_emulsion_area(b, t, x):
            return (b.gas_emulsion.area[t, x] ==
                    b.bed_area*b.delta_e[t, x]*b.voidage_emulsion[t, x])

        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Cross-sectional Area Occupied by Solid Emulsion")
        def solid_emulsion_area(b, t, x):
            return (b.solid_emulsion.area[t, x] ==
                    b.bed_area*b.delta_e[t, x]*(1-b.voidage_emulsion[t, x]))

        # Length of bubble, gas_emulsion, solid_emulsion
        @self.Constraint(doc="Bubble Region Length")
        def bubble_length(b):
            return (b.bubble.length == b.bed_height)

        @self.Constraint(doc="Gas Emulsion Region Length")
        def gas_emulsion_length(b):
            return (b.gas_emulsion.length == b.bed_height)

        @self.Constraint(doc="Solid Emulsion Region Length")
        def solid_emulsion_length(b):
            return (b.solid_emulsion.length == b.bed_height)

        # ---------------------------------------------------------------------
        # Hydrodynamic contraints

        # Emulsion region volume fraction
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Emulsion Region Volume Fraction")
        def emulsion_vol_frac(b, t, x):
            return (b.delta_e[t, x] == 1 - b.delta[t, x])

        # Average cross-sectional voidage
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Average Cross-sectional Voidage")
        def average_voidage(b, t, x):
            return (1 - b.voidage_average[t, x] == (1 - b.delta[t, x]) *
                    (1 - b.voidage_emulsion[t, x]))

        # Bubble growth coefficient
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Bubble Growth Coefficient")
        def bubble_growth_coefficient(b, t, x):
            return (1e2*(b.bubble_growth_coeff[t, x] *
                         b.solid_emulsion.properties[t, x]._params.velocity_mf
                         )**2 == 1e2 *
                    (2.56e-2**2) * (
                        b.bed_diameter/constants.acceleration_gravity))

        # Maximum bubble diameter
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Maximum Bubble Diameter")
        def bubble_diameter_maximum(b, t, x):
            return ((b.bubble_diameter_max[t, x]**5) *
                    constants.acceleration_gravity ==
                    (2.59**5)*((b.velocity_superficial_gas[t, x] -
                                b.velocity_emulsion_gas[t, x]) *
                               b.bed_area)**2)

        # Bubble diameter reformulation equation
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Bubble Diameter Reformulation")
        def _reformulation_eqn_1(b, t, x):
            return (b._reform_var_1[t, x]**2 ==
                    b.bed_diameter*b.bubble_diameter[t, x])

        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Bubble Diameter")
        def bubble_diameter_eqn(b, t, x):
            if x == b.length_domain.first():
                return (1e3*b.bubble_diameter[t, x] ** 5 == 1e3 * (1.38 ** 5) *
                        (constants.acceleration_gravity ** -1) *
                        ((b.velocity_superficial_gas[t, x]
                            - b.velocity_emulsion_gas[t, x]) *
                        b.area_orifice) ** 2)
            else:
                return (1e2 * b.ddia_bubbledx[t, x] * b.bed_diameter ==
                        1e2 * b.bubble.length *
                        0.3 * (b.bubble_diameter_max[t, x]
                               - b.bubble_diameter[t, x] -
                        b.bubble_growth_coeff[t, x]*b._reform_var_1[t, x]))

        # Bubble rise velocity
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Bubble Rise Velocity")
        def bubble_velocity_rise(b, t, x):
            return (b.velocity_bubble_rise[t, x]**2 == (0.711**2) *
                    constants.acceleration_gravity *
                    b.bubble_diameter[t, x])

        # Emulsion region voidage - Davidson model
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Emulsion Region Voidage")
        def emulsion_voidage(b, t, x):
            return (b.voidage_emulsion[t, x] ==
                    b.solid_emulsion.properties[t, x]._params.voidage_mf)

        # Bubble velocity - Davidson model
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Bubble Velocity")
        def bubble_velocity(b, t, x):
            return (
                b.velocity_bubble[t, x] ==
                b.velocity_superficial_gas[t, x] -
                b.solid_emulsion.properties[t, x]._params.velocity_mf +
                b.velocity_bubble_rise[t, x])

        # Average gas density (moles)
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Average gas density")
        def average_gas_density_eqn(b, t, x):
            return (
                b.average_gas_density[t, x] *
                (b.bubble.properties[t, x].flow_mol +
                 b.gas_emulsion.properties[t, x].flow_mol) ==
                ((b.bubble.properties[t, x].flow_mol *
                  b.bubble.properties[t, x].dens_mol) +
                 (b.gas_emulsion.properties[t, x].flow_mol *
                  b.gas_emulsion.properties[t, x].dens_mol))
                )

        # Superficial gas velocity
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Superficial Gas Velocity")
        def velocity_gas_superficial(b, t, x):
            return (
                b.bubble.properties[t, x].flow_mol +
                b.gas_emulsion.properties[t, x].flow_mol
                ==
                b.velocity_superficial_gas[t, x] * b.bed_area *
                b.average_gas_density[t, x])

        # Bubble volume fraction
        # This computes delta for emulsion at minimum fluidization velocity
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="volume fraction occupied by bubbles")
        def bubble_vol_frac_eqn(b, t, x):
            return (b.velocity_superficial_gas[t, x] ==
                    b.velocity_bubble[t, x] * b.delta[t, x] +
                    b.solid_emulsion.properties[t, x]._params.velocity_mf)

        # Solid superficial velocity
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Solid superficial velocity")
        # This equation uses inlet values to compute the constant solid
        # superficial velocity, and then computes the
        # solid particle density or flowrate through the rest of the bed.
        def solid_super_vel(b, t, x):
            return (b.velocity_superficial_solid[t] * b.bed_area *
                    b.solid_emulsion.properties[t, x].dens_mass_particle ==
                    b.solid_emulsion.properties[t, x].flow_mass)

        # Gas_emulsion pressure drop calculation
        if self.config.has_pressure_change:
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Gas Emulsion Pressure Drop Calculation")
            def gas_emulsion_pressure_drop(b, t, x):
                # 1e5 = pressure unit conversion factor from Pa to bar
                return (1e-2*(b.gas_emulsion.deltaP[t, x] *
                              1e5) ==
                        1e-2*(- constants.acceleration_gravity *
                              (1 - b.voidage_average[t, x]) *
                        b.solid_emulsion.properties[t, x].dens_mass_particle)
                        )

        elif self.config.has_pressure_change is False:
            # If pressure change is false set pressure drop to zero
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Isobaric Gas emulsion")
            def isobaric_gas_emulsion(b, t, x):
                return (1e2*b.gas_emulsion.properties[t, x].pressure ==
                        1e2*b.gas_inlet.pressure[0])

        # ---------------------------------------------------------------------
        # Mass transfer constraints
        @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                gas_phase.property_package.component_list,
                doc="Bubble to Emulsion Gas Mass Transfer"
                    "Coefficient Reformulation [reform eqn 2]")
        def _reformulation_eqn_2(b, t, x, j):
            # 1e-4 = diffusion unit conversion factor from cm2/s to m2/s
            return (1e2*b._reform_var_2[t, x, j]**2 == 1e2 *
                    34.2225 * (
                    1e-4 *
                    b.gas_emulsion.properties[t, x].diffusion_comp[j]) *
                    constants.acceleration_gravity ** 0.5)

        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Bubble to Emulsion Gas Mass Transfer"
                             "Coefficient Reformulation [reform eqn 3]")
        def _reformulation_eqn_3(b, t, x):
            return (1e2*b._reform_var_3[t, x]**4 ==
                    1e2*b.bubble_diameter[t, x])

        @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                gas_phase.property_package.component_list,
                doc="Bubble to Emulsion Gas Mass Transfer Coefficient")
        def bubble_cloud_mass_trans_coeff(b, t, x, j):
            return (
                b.Kbe[t, x, j] * b._reform_var_3[t, x]**5 ==
                0.36 * 4.5 * b._reform_var_3[t, x] *
                b.solid_emulsion.properties[t, x]._params.velocity_mf +
                b._reform_var_2[t, x, j])

        # Bulk gas mass transfer
        @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                gas_phase.property_package.component_list,
                doc="Bulk Gas Mass Transfer Between Bubble and Emulsion")
        def bubble_cloud_bulk_mass_trans(b, t, x, j):
            conc_diff = (b.gas_emulsion.properties[t, x].dens_mol -
                         b.bubble.properties[t, x].dens_mol)
            return (b.Kgbulk_c[t, x, j] * b.bubble_diameter[t, x] ==
                    6 * b.Kd * b.delta[t, x] * b.bed_area *
                    (b.gas_emulsion.properties[t, x].mole_frac_comp[j] *
                     smooth_max(conc_diff, 0, b.eps) +
                     b.bubble.properties[t, x].mole_frac_comp[j] *
                     smooth_min(conc_diff, 0, b.eps)))
        # ---------------------------------------------------------------------
        # Heat transfer constraints

        if self.config.energy_balance_type != EnergyBalanceType.none:
            # Bubble to emulsion gas heat transfer coefficient
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Bubble to Emulsion Gas Heat Transfer"
                                 "Coeff. Reformulation Eqn [reform eqn 4]")
            def _reformulation_eqn_4(b, t, x):
                return (b._reform_var_4[t, x]**2 ==
                        34.2225 * b.bubble.properties[t, x].therm_cond *
                        b.bubble.properties[t, x].enth_mol *
                        b.bubble.properties[t, x].dens_mol *
                        (constants.acceleration_gravity ** 0.5))

            # Convective heat transfer coefficient
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Convective Heat Transfer"
                                 "Coeff. Reformulation Eqn [reform eqn 5]")
            def _reformulation_eqn_5(b, t, x):
                return (
                    1e2 * b._reform_var_5[t, x] *
                    b.gas_emulsion.properties[t, x].visc_d ==
                    1e2 * b.velocity_emulsion_gas[t, x] *
                    b.solid_emulsion.properties[t, x]._params.particle_dia *
                    b.gas_emulsion.properties[t, x].dens_mol)

            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Bubble to Emulsion Gas Heat Transfer"
                                 "Coefficient")
            def bubble_cloud_heat_trans_coeff(b, t, x):
                return (
                    b.Hbe[t, x] * b._reform_var_3[t, x] ** 5 == 4.5 *
                    b.solid_emulsion.properties[t, x]._params.velocity_mf *
                    b.bubble.properties[t, x].enth_mol *
                    b.bubble.properties[t, x].dens_mol *
                    b._reform_var_3[t, x] +
                    b._reform_var_4[t, x])

            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Convective Heat Transfer Coefficient")
            def convective_heat_trans_coeff(b, t, x):
                return (
                    1e6*b.htc_conv[t, x] *
                    b.solid_emulsion.properties[t, x]._params.particle_dia ==
                    1e6 * 0.03 * b.gas_emulsion.properties[t, x].therm_cond *
                    ((b._reform_var_5[t, x]**2 + b.eps)**0.5) ** 1.3)

            # Gas to solid convective heat transfer # replaced "ap" with
            # "6/(dp*rho_sol)"
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Gas to Solid Convective Enthalpy Transfer"
                                 "in Emulsion Region")
            def convective_heat_transfer(b, t, x):
                return (
                    b.ht_conv[t, x] *
                    b.solid_emulsion.properties[t, x]._params.particle_dia ==
                    6 * b.delta_e[t, x] * (1 - b.voidage_emulsion[t, x]) *
                    b.htc_conv[t, x] *
                    (b.gas_emulsion.properties[t, x].temperature -
                     b.solid_emulsion.properties[t, x].temperature))

            # Bulk gas heat transfer
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Bulk Gas Heat Transfer Between"
                                 "Bubble and Emulsion")
            def bubble_cloud_bulk_heat_trans(b, t, x):
                conc_diff = (
                        b.gas_emulsion.properties[t, x].dens_mol -
                        b.bubble.properties[t, x].dens_mol)
                return (
                    b.Hgbulk[t, x] * b.bubble_diameter[t, x] ==
                    6 * b.Kd * b.delta[t, x] * b.bed_area *
                    (b.gas_emulsion.properties[t, x].enth_mol *
                     smooth_max(conc_diff, 0, b.eps) +
                     b.bubble.properties[t, x].enth_mol *
                     smooth_min(conc_diff, 0, b.eps)))
        # ---------------------------------------------------------------------
        # Mass and heat transfer terms in control volumes

        # Bubble mass transfer
        @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                gas_phase.property_package.component_list,
                doc="Bubble Mass Transfer")
        def bubble_mass_transfer(b, t, x, j):
            comp_conc_diff = (
                    b.bubble.properties[t, x].dens_mol_comp[j] -
                    b.gas_emulsion.properties[t, x].dens_mol_comp[j])
            return (1e3*b.bubble.mass_transfer_term[t, x, 'Vap', j] ==
                    1e3*(b.Kgbulk_c[t, x, j] -
                    b.bubble.area[t, x] *
                         b.Kbe[t, x, j] * comp_conc_diff))

        # Gas_emulsion mass transfer
        def gas_emulsion_hetero_rxn_term(b, t, x, j):
            return (b.gas_emulsion_hetero_rxn[t, x, j]
                    if solid_phase.reaction_package else 0)

        @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                gas_phase.property_package.component_list,
                doc="Gas Emulsion Mass Transfer")
        def gas_emulsion_mass_transfer(b, t, x, j):
            comp_conc_diff = (
                    b.bubble.properties[t, x].dens_mol_comp[j] -
                    b.gas_emulsion.properties[t, x].dens_mol_comp[j])
            return (1e3*b.gas_emulsion.mass_transfer_term[t, x, 'Vap', j] ==
                    1e3*(- b.Kgbulk_c[t, x, j] +
                    b.bubble.area[t, x] *
                         b.Kbe[t, x, j] * comp_conc_diff +
                         gas_emulsion_hetero_rxn_term(b, t, x, j)))

        if self.config.energy_balance_type != EnergyBalanceType.none:
            # Bubble - heat transfer
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Bubble - Heat Transfer")
            def bubble_heat_transfer(b, t, x):
                return (b.bubble.heat[t, x] == b.Hgbulk[t, x] -
                        b.Hbe[t, x] *
                        (b.bubble.properties[t, x].temperature -
                         b.gas_emulsion.properties[t, x].temperature) *
                        b.bubble.area[t, x])

            # Gas emulsion - heat transfer
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Gas Emulsion - Heat Transfer")
            def gas_emulsion_heat_transfer(b, t, x):
                return (b.gas_emulsion.heat[t, x] ==
                        b.Hbe[t, x] *
                        (b.bubble.properties[t, x].temperature -
                         b.gas_emulsion.properties[t, x].temperature) *
                        b.bubble.area[t, x] - b.Hgbulk[t, x] -
                        b.ht_conv[t, x] * b.bed_area)

            # Solid emulsion - heat transfer
            @self.Constraint(self.flowsheet().config.time,
                             self.length_domain,
                             doc="Solid Emulsion - Heat Transfer")
            def solid_emulsion_heat_transfer(b, t, x):
                return (b.solid_emulsion.heat[t, x] ==
                        b.ht_conv[t, x] * b.bed_area)

        # ---------------------------------------------------------------------
        # Reaction  contraints

        # Build homogeneous reaction constraints

        if gas_phase.reaction_package is not None:
            # Bubble rate reaction extent
            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                gas_phase.reaction_package.rate_reaction_idx,
                doc="Bubble Rate Reaction Extent")
            def bubble_rxn_ext_constraint(b, t, x, r):
                return b.bubble.rate_reaction_extent[t, x, r] == (
                        b.bubble.reactions[t, x].reaction_rate[r] *
                        b.bubble.area[t, x])

        if gas_phase.reaction_package is not None:
            # Gas emulsion rate reaction extent
            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                gas_phase.reaction_package.rate_reaction_idx,
                doc="Gas Emulsion Rate Reaction Extent")
            def gas_emulsion_rxn_ext_constraint(b, t, x, r):
                return b.gas_emulsion.rate_reaction_extent[t, x, r] == (
                        b.gas_emulsion.reactions[t, x].reaction_rate[r] *
                        b.gas_emulsion.area[t, x])

        # Build hetereogeneous reaction constraints
        if solid_phase.reaction_package is not None:
            # Solid side rate reaction extent
            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                solid_phase.reaction_package.rate_reaction_idx,
                doc="Solid Emulsion Rate Reaction Extent")
            def solid_emulsion_rxn_ext_constraint(b, t, x, r):
                return (
                    b.solid_emulsion.rate_reaction_extent[t, x, r] ==
                    b.solid_emulsion.reactions[t, x].reaction_rate[r] *
                    b.solid_emulsion.area[t, x])

            # Gas side heterogeneous rate reaction generation
            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                gas_phase.property_package.component_list,
                doc="Gas Emulsion Heterogeneous Rate Reaction Generation")
            def gas_emulsion_hetero_rxn_eqn(b, t, x, j):
                return (
                    b.gas_emulsion_hetero_rxn[t, x, j] ==
                    sum(
                        b.solid_emulsion.reactions[t, x].
                        rate_reaction_stoichiometry[r, 'Vap', j] *
                        b.solid_emulsion.reactions[t, x].reaction_rate[r]
                        for r in solid_phase.reaction_package.rate_reaction_idx
                        ) * b.solid_emulsion.area[t, x])

        # ---------------------------------------------------------------------
        # Flowrate constraints

        # Bubble gas flowrate - this eqn calcs bubble conc. (dens_mol)
        # which is then used to calculate the bubble pressure
        @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Bubble Gas Flowrate Constraint")
        def bubble_gas_flowrate(b, t, x):
            return (b.bubble.properties[t, x].flow_mol ==
                    b.bed_area * b.delta[t, x] * b.velocity_bubble[t, x] *
                    b.bubble.properties[t, x].dens_mol)

        # Emulsion gas flowrate - this eqn indirectly calcs bubble voidage
        # (delta) in conjuction with the mass conservation eqns in the
        # gas_emulsion CV1D.
        # This eqn arises because the emulsion region gas is assumed to be at
        # minimum fluidization conditions (delta varies to account for this)
        @self.Constraint(self.flowsheet().config.time,
                         self.length_domain,
                         doc="Emulsion Gas Flowrate Constraint")
        def emulsion_gas_flowrate(b, t, x):
            return (b.gas_emulsion.properties[t, x].flow_mol ==
                    b.bed_area * b.velocity_emulsion_gas[t, x] *
                    b.gas_emulsion.properties[t, x].dens_mol)

        # ---------------------------------------------------------------------
        # Inlet boundary Conditions

        # Gas_emulsion pressure at inlet
        if self.config.has_pressure_change:
            @self.Constraint(self.flowsheet().config.time,
                             doc="Gas Emulsion Pressure at Inlet")
            def gas_emulsion_pressure_in(b, t):
                return (1e2*b.gas_emulsion.properties[t, 0].pressure ==
                        1e2*b.gas_inlet_block[t].pressure - b.deltaP_orifice)

        # Total gas balance at inlet
        @self.Constraint(self.flowsheet().config.time,
                         doc="Total Gas Balance at Inlet")
        def gas_mole_flow_in(b, t):
            return (b.gas_inlet_block[t].flow_mol ==
                    b.bubble.properties[t, 0].flow_mol +
                    b.gas_emulsion.properties[t, 0].flow_mol)

        # Particle porosity at inlet
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constant particle porosity")
        def particle_porosity_in(b, t):
            x = b.length_domain.first()
            return (
                b.solid_emulsion.properties[t, x].particle_porosity ==
                b.solid_inlet_block[t].particle_porosity)

        # Emulsion region gas velocity at inlet - Davidson model
        @self.Constraint(self.flowsheet().config.time,
                         doc="Emulsion Region Superficial Gas Velocity")
        def emulsion_gas_velocity_in(b, t):
            x = self.length_domain.first()
            return (
                b.velocity_emulsion_gas[t, x] ==
                b.solid_emulsion.properties[t, x]._params.velocity_mf)

        # Bubble mole frac at inlet
        @self.Constraint(
                self.flowsheet().config.time,
                gas_phase.property_package.component_list,
                doc="Bubble Mole Fraction at Inlet")
        def bubble_mole_frac_in(b, t, j):
            return (1e2*b.gas_inlet_block[t].mole_frac_comp[j] ==
                    1e2*b.bubble.properties[t, 0].mole_frac_comp[j])

        # Gas_emulsion mole frac at inlet
        @self.Constraint(
                self.flowsheet().config.time,
                gas_phase.property_package.component_list,
                doc="Gas Emulsion Mole Fraction at Inlet")
        def gas_emulsion_mole_frac_in(b, t, j):
            return (1e2*b.gas_inlet_block[t].mole_frac_comp[j] ==
                    1e2*b.gas_emulsion.properties[t, 0].mole_frac_comp[j])

        # Solid emulsion mass flow at inlet
        @self.Constraint(self.flowsheet().config.time,
                         doc="Solid Emulsion Mass Flow at Inlet")
        def solid_emulsion_mass_flow_in(b, t):
            if (self.config.flow_type == "co_current"):
                return (b.solid_inlet_block[t].flow_mass ==
                        b.solid_emulsion.properties[t, 0].flow_mass)
            elif (self.config.flow_type == "counter_current"):
                return (b.solid_inlet_block[t].flow_mass ==
                        b.solid_emulsion.properties[t, 1].flow_mass)

        # Solid emulsion mass frac at inlet
        @self.Constraint(
                self.flowsheet().config.time,
                solid_phase.property_package.component_list,
                doc="Solid Emulsion Mass Fraction at Inlet")
        def solid_emulsion_mass_frac_in(b, t, j):
            if (self.config.flow_type == "co_current"):
                return (1e2 * b.solid_inlet_block[t].mass_frac_comp[j] ==
                        1e2 *
                        b.solid_emulsion.properties[t, 0].mass_frac_comp[j])
            elif (self.config.flow_type == "counter_current"):
                return (1e2 * b.solid_inlet_block[t].mass_frac_comp[j] ==
                        1e2 *
                        b.solid_emulsion.properties[t, 1].mass_frac_comp[j])

        if self.config.energy_balance_type != EnergyBalanceType.none:
            @self.Constraint(
                    self.flowsheet().config.time,
                    gas_phase.property_package.phase_list,
                    doc="Gas Inlet Energy Balance")
            def gas_energy_balance_in(b, t, p):
                return (b.gas_inlet_block[t].get_enthalpy_flow_terms(p) ==
                        b.bubble._enthalpy_flow[t, 0, p] +
                        b.gas_emulsion._enthalpy_flow[t, 0, p])

            # Gas emulsion temperature at inlet
            @self.Constraint(
                    self.flowsheet().config.time,
                    doc="Gas Emulsion Temperature at Inlet")
            def gas_emulsion_temperature_in(b, t):
                return (b.gas_inlet_block[t].temperature ==
                        b.gas_emulsion.properties[t, 0].temperature)

            # Solid inlet energy balance
            @self.Constraint(
                    self.flowsheet().config.time,
                    doc="Solid Inlet Energy Balance")
            def solid_energy_balance_in(b, t):
                if (self.config.flow_type == "co_current"):
                    return (b.solid_inlet_block[t].
                            get_enthalpy_flow_terms('Sol') ==
                            b.solid_emulsion._enthalpy_flow[t, 0, 'Sol'])
                elif (self.config.flow_type == "counter_current"):
                    return (b.solid_inlet_block[t].
                            get_enthalpy_flow_terms('Sol') ==
                            b.solid_emulsion._enthalpy_flow[t, 1, 'Sol'])

        elif self.config.energy_balance_type == EnergyBalanceType.none:
            # If energy balance is none fix gas and solid temperatures to inlet
            @self.Constraint(
                    self.flowsheet().config.time,
                    self.length_domain,
                    doc="Isothermal solid emulsion constraint")
            def isothermal_solid_emulsion(b, t, x):
                return (
                        b.solid_emulsion.properties[t, x].temperature ==
                        b.solid_inlet.temperature[t])

            @self.Constraint(
                    self.flowsheet().config.time,
                    self.length_domain,
                    doc="Isothermal gas emulsion constraint")
            def isothermal_gas_emulsion(b, t, x):
                return (
                        b.gas_emulsion.properties[t, x].temperature ==
                        b.gas_inlet.temperature[t])

            @self.Constraint(
                    self.flowsheet().config.time,
                    self.length_domain,
                    doc="Isothermal gas emulsion constraint")
            def isothermal_bubble(b, t, x):
                return (
                        b.bubble.properties[t, x].temperature ==
                        b.gas_inlet.temperature[t])

        # ---------------------------------------------------------------------
        # Outlet boundary Conditions
        # Gas outlet pressure
        @self.Constraint(self.flowsheet().config.time,
                         doc="Gas Outlet Pressure")
        def gas_pressure_out(b, t):
            return (1e2*b.gas_outlet.pressure[t] ==
                    1e2*b.gas_emulsion.properties[t, 1].pressure)

        # Gas outlet material balance
        @self.Constraint(
                self.flowsheet().config.time,
                gas_phase.property_package.phase_list,
                gas_phase.property_package.component_list,
                doc="Gas Outlet Material Balance")
        def gas_material_balance_out(b, t, p, j):
            return (b.gas_outlet_block[t].get_material_flow_terms(p, j) ==
                    b.bubble._flow_terms[t, 1, p, j] +
                    b.gas_emulsion._flow_terms[t, 1, p, j])

        # Solid outlet material balance
        @self.Constraint(
                self.flowsheet().config.time,
                solid_phase.property_package.phase_list,
                solid_phase.property_package.component_list,
                doc="Solid Outlet Material Balance")
        def solid_material_balance_out(b, t, p, j):
            if (self.config.flow_type == "co_current"):
                return (
                    b.solid_outlet_block[t].get_material_flow_terms(p, j) ==
                    b.solid_emulsion._flow_terms[t, 1, p, j])
            elif (self.config.flow_type == "counter_current"):
                return (
                    b.solid_outlet_block[t].get_material_flow_terms(p, j) ==
                    b.solid_emulsion._flow_terms[t, 0, p, j])

        # Solid outlet particle porosity
        @self.Constraint(
                self.flowsheet().config.time)
        def solid_particle_porosity_out(b, t):
            x = b.length_domain.last()
            return (
                b.solid_outlet_block[t].particle_porosity ==
                b.solid_emulsion.properties[t, x].particle_porosity)

        if self.config.energy_balance_type != EnergyBalanceType.none:
            # Gas outlet energy balance
            @self.Constraint(
                    self.flowsheet().config.time,
                    doc="Gas Outlet Energy Balance")
            def gas_energy_balance_out(b, t):
                return (b.gas_outlet_block[t].get_enthalpy_flow_terms('Vap') ==
                        b.bubble._enthalpy_flow[t, 1, 'Vap'] +
                        b.gas_emulsion._enthalpy_flow[t, 1, 'Vap'])

            # Solid outlet energy balance
            @self.Constraint(
                    self.flowsheet().config.time,
                    doc="Solid Outlet Energy Balance")
            def solid_energy_balance_out(b, t):
                if (self.config.flow_type == "co_current"):
                    return (
                        b.solid_outlet_block[t].get_enthalpy_flow_terms('Sol')
                        ==
                        b.solid_emulsion._enthalpy_flow[t, 1, 'Sol'])
                elif (self.config.flow_type == "counter_current"):
                    return (
                        b.solid_outlet_block[t].get_enthalpy_flow_terms('Sol')
                        ==
                        b.solid_emulsion._enthalpy_flow[t, 0, 'Sol'])

        elif self.config.energy_balance_type == EnergyBalanceType.none:
            # Gas outlet energy balance
            @self.Constraint(
                    self.flowsheet().config.time,
                    doc="Gas Outlet Energy Balance")
            def gas_energy_balance_out(b, t):
                return (b.gas_outlet_block[t].temperature ==
                        b.gas_emulsion.properties[t, 1].temperature)

            # Solid outlet energy balance
            @self.Constraint(
                    self.flowsheet().config.time,
                    doc="Solid Outlet Energy Balance")
            def solid_energy_balance_out(b, t):
                if (self.config.flow_type == "co_current"):
                    return (b.solid_outlet_block[t].temperature ==
                            b.solid_emulsion.properties[t, 1].temperature)
                elif (self.config.flow_type == "counter_current"):
                    return (b.solid_outlet_block[t].temperature ==
                            b.solid_emulsion.properties[t, 0].temperature)

    # =========================================================================
    # Model initialization routine

    def initialize(blk, gas_phase_state_args=None, solid_phase_state_args=None,
                   outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        """
        Initialisation routine for Bubbling Fluidized Bed unit

        Keyword Arguments:
            gas_phase_state_args : a dict of arguments to be passed to the
                        property package(s) to provide an initial state for
                        initialization (see documentation of the specific
                        property package) (default = None).
            solid_phase_state_args : a dict of arguments to be passed to the
                        property package(s) to provide an initial state for
                        initialization (see documentation of the specific
                        property package) (default = None).
            outlvl : sets output level of initialisation routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        # Set logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # local aliases used to shorten object names
        gas_phase = blk.config.gas_phase_config
        solid_phase = blk.config.solid_phase_config

        # Keep all unit model geometry constraints, derivative_var constraints,
        # and property block constraints active. Additionaly, in control
        # volumes - keep conservation linking constraints and
        # holdup calculation (for dynamic flowsheets) constraints active

        geometry_constraints_terms = ["orifice_area", "bed_area_eqn",
                                      "bubble_area", "gas_emulsion_area",
                                      "solid_emulsion_area", "bubble_length",
                                      "gas_emulsion_length",
                                      "solid_emulsion_length"]
        endswith_terms = ("_disc_eq", "linking_constraint",
                          "linking_constraints", "_holdup_calculation")
        startswith_terms = ("properties", "gas_inlet_block",
                            "gas_outlet_block", "solid_inlet_block",
                            "solid_outlet_block")

        for c in blk.component_objects(Constraint, descend_into=True):
            if not c.parent_block().local_name.startswith(startswith_terms) \
                and not c.local_name.endswith(endswith_terms) \
                    and c.local_name not in geometry_constraints_terms:
                c.deactivate()

        # Deactivate outlet blocks (activate at last initialization solve)
        blk.gas_outlet_block.deactivate()
        blk.solid_outlet_block.deactivate()

        # ---------------------------------------------------------------------
        # Fix Initial Values of State Variables and initialize unit model
        # and inlet property blocks
        init_log.info('Initialize Property Block Constraints')

        # Initialize inlet property blocks
        blk.gas_inlet_block.initialize(
            state_args=gas_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver)
        blk.solid_inlet_block.initialize(
            state_args=solid_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver)

        # Initialize bubble region block
        bubble_region_flags = (
                blk.bubble.properties.initialize(
                    state_args=gas_phase_state_args,
                    hold_state=True,
                    outlvl=outlvl,
                    optarg=optarg,
                    solver=solver))

        # Initialize gas_emulsion region block
        gas_emulsion_region_flags = (
                blk.gas_emulsion.properties.initialize(
                    state_args=gas_phase_state_args,
                    hold_state=True,
                    outlvl=outlvl,
                    optarg=optarg,
                    solver=solver))

        # Initialize solid_emulsion properties block
        solid_emulsion_region_flags = (
                blk.solid_emulsion.properties.initialize(
                    state_args=solid_phase_state_args,
                    hold_state=True,
                    outlvl=outlvl,
                    optarg=optarg,
                    solver=solver))

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Initialize geometric constraints, property block constraints
        # and reaction block constraints
        # Fix delta, delta_e and void_emul to inital values for square problem
        blk.delta.fix()
        blk.delta_e.fix()
        blk.voidage_emulsion.fix()
        blk.bubble_diameter.fix()  # Fix is needed because of its DerivativeVar

        init_log.info('Initialize Geometric Constraints')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
        if results.solver.termination_condition \
                == TerminationCondition.optimal:
            init_log.info_high(
                "Initialization Step 2 {}.".format(
                        idaeslog.condition(results))
                        )
        else:
            init_log.warning('{} Initialisation Step 2 Failed.'
                             .format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize hydrodynamics
        # vel_superficial_gas, delta are fixed during this stage
        for t in blk.flowsheet().config.time:
            for x in blk.length_domain:
                # Superficial velocity initialized at 3 * min fluidization vel.
                blk.velocity_superficial_gas[t, x].fix(
                    value(3 * blk.solid_inlet_block[t]._params.velocity_mf))

                # ve is fixed across all bed during initialization of
                # hydrodynamic sub-model. It will be unfixed during
                # the mass balance initialization when
                # ve = f(gas_emulsion flowrate)
                # is activated, and the ve=vmf at inlet boundary is activated
                blk.velocity_emulsion_gas[t, x].fix(value(
                        blk.solid_inlet_block[t]._params.velocity_mf))
                blk.bubble_diameter[t, x] = value(
                        1.38 * (constants.acceleration_gravity**(-0.2)) *
                        ((blk.velocity_superficial_gas[t, x] -
                          blk.velocity_emulsion_gas[t, x]) *
                            (1/blk.number_orifice))**0.4)
                blk.bubble_diameter_max[t, x] = value(
                        2.59 * (constants.acceleration_gravity**(-0.2)) *
                        ((blk.velocity_superficial_gas[t, x] -
                          blk.velocity_emulsion_gas[t, x]) *
                            ((constants.pi/4) * blk.bed_diameter**2))**0.4)
                blk.bubble_growth_coeff[t, x] = value(
                        2.56e-2 * sqrt(blk.bed_diameter /
                                       constants.acceleration_gravity) /
                        blk.solid_inlet_block[t]._params.velocity_mf)
                blk.velocity_bubble_rise[t, x] = value(
                        0.711 * sqrt(constants.acceleration_gravity *
                                     blk.bubble_diameter[t, x]))
                blk.velocity_bubble[t, x] = (
                        blk.velocity_superficial_gas[t, x].value +
                        blk.velocity_bubble_rise[t, x].value -
                        blk.solid_inlet_block[t]._params.velocity_mf.value)
                blk.delta[t, x].fix((blk.velocity_superficial_gas[t, x].value -
                                     blk.velocity_emulsion_gas[t, x].value) /
                                    blk.velocity_bubble[t, x].value)
                blk.delta_e[t, x] = (1 - blk.delta[t, x].value)
                blk.voidage_emulsion[t, x] = (
                        blk.solid_inlet_block[t]._params.voidage_mf.value)
                blk.voidage_average[t, x] = (
                        1 - (1 - blk.voidage_emulsion[t, x].value) *
                        (1 - blk.delta[t, x].value))
                blk._reform_var_1[t, x] = sqrt(
                        blk.bed_diameter.value *
                        blk.bubble_diameter[t, x].value)

        # Unfix variables to make problem square
        blk.delta.unfix()
        blk.delta_e.unfix()
        blk.voidage_emulsion.unfix()
        blk.bubble_diameter.unfix()

        # Activate relavant constraints
        blk.bubble_vol_frac_eqn.activate()  # delta
        blk.average_gas_density_eqn.activate()

        blk.emulsion_vol_frac.activate()
        blk.average_voidage.activate()
        blk.bubble_growth_coefficient.activate()
        blk.bubble_diameter_maximum.activate()
        blk._reformulation_eqn_1.activate()
        blk.bubble_diameter_eqn.activate()
        blk.bubble_velocity_rise.activate()
        blk.bubble_velocity.activate()
        blk.emulsion_voidage.activate()

        init_log.info('Initialize Hydrodynamics')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
        if results.solver.termination_condition \
                == TerminationCondition.optimal:
            init_log.info_high(
                "Initialization Step 3 {}.".format(
                        idaeslog.condition(results))
                        )
        else:
            init_log.warning('{} Initialisation Step 3 Failed.'
                             .format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize mass balance - no reaction

        # Initialize variables
        for t in blk.flowsheet().config.time:
            for x in blk.length_domain:
                calculate_variable_from_constraint(
                    blk._reform_var_3[t, x],
                    blk._reformulation_eqn_3[t, x])
                for j in (gas_phase.property_package.component_list):
                    calculate_variable_from_constraint(
                        blk._reform_var_2[t, x, j],
                        blk._reformulation_eqn_2[t, x, j])
                    calculate_variable_from_constraint(
                        blk.Kbe[t, x, j],
                        blk.bubble_cloud_mass_trans_coeff[t, x, j])

        # Unfix variables, and fix reaction rate variables (no rxns assumed)
        # Unfix material balance state variables but keep other states fixed
        blk.bubble.properties.release_state(bubble_region_flags)
        blk.gas_emulsion.properties.release_state(
                gas_emulsion_region_flags)
        blk.solid_emulsion.properties.release_state(
                solid_emulsion_region_flags)

        for t in blk.flowsheet().config.time:
            for x in blk.length_domain:
                blk.gas_emulsion.properties[t, x].pressure.fix()
                blk.bubble.properties[t, x].temperature.fix()
                blk.gas_emulsion.properties[t, x].temperature.fix()
                blk.solid_emulsion.properties[t, x].temperature.fix()

        # Fix reaction rate variables (no rxns assumed)
        # Homogeneous reactions (gas phase rxns)
        if gas_phase.reaction_package is not None:
            for t in blk.flowsheet().config.time:
                bubble_rxn_gen = blk.bubble.rate_reaction_generation
                gas_emulsion_rxn_gen = (
                    blk.gas_emulsion.rate_reaction_generation)
                for x in blk.length_domain:
                    for j in gas_phase.property_package.component_list:
                        # Bubble region
                        (bubble_rxn_gen[t, x, 'Vap', j].fix(0.0))
                        # Gas emulsion region
                        (gas_emulsion_rxn_gen[t, x, 'Vap', j].fix(0.0))

        # Heterogeneous rxns (solid phase rxns with gas phase interactions)
        if solid_phase.reaction_package is not None:
            # local alias
            solid_emulsion_rxn_gen = (
                blk.solid_emulsion.rate_reaction_generation)
            for t in blk.flowsheet().config.time:
                # Solid emulsion region
                for x in blk.length_domain:
                    for j in solid_phase.property_package.component_list:
                        (solid_emulsion_rxn_gen[t, x, 'Sol', j].fix(0.0))
                # Gas emulsion region
                for x in blk.length_domain:
                    for j in gas_phase.property_package.component_list:
                        blk.gas_emulsion_hetero_rxn[t, x, j].fix(0.0)

        # Unfix velocity_emulsion_gas
        blk.velocity_emulsion_gas.unfix()

        # Unfix velocity_superficial_gas
        blk.velocity_superficial_gas.unfix()

        # Activate relevant model level constraints
        blk.velocity_gas_superficial.activate()
        blk.solid_super_vel.activate()

        blk._reformulation_eqn_2.activate()
        blk._reformulation_eqn_3.activate()
        blk.bubble_cloud_mass_trans_coeff.activate()
        blk.bubble_cloud_bulk_mass_trans.activate()
        blk.bubble_mass_transfer.activate()
        blk.gas_emulsion_mass_transfer.activate()
        blk.bubble_gas_flowrate.activate()
        blk.emulsion_gas_flowrate.activate()

        # Activate relevant boundary constraints
        blk.gas_mole_flow_in.activate()
        blk.bubble_mole_frac_in.activate()
        blk.gas_emulsion_mole_frac_in.activate()
        blk.solid_emulsion_mass_flow_in.activate()
        blk.solid_emulsion_mass_frac_in.activate()
        blk.particle_porosity_in.activate()
        blk.emulsion_gas_velocity_in.activate()

        # Activate relevant control volume constraints
        blk.bubble.material_balances.activate()
        blk.gas_emulsion.material_balances.activate()
        blk.solid_emulsion.material_balances.activate()

        # Initialize, unfix and activate pressure drop related
        # variables and constraints
        if blk.config.has_pressure_change:
            blk.gas_emulsion_pressure_in.activate()
            blk.gas_emulsion_pressure_drop.activate()
            blk.gas_emulsion.pressure_balance.activate()

            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    blk.gas_emulsion.properties[t, x].pressure.unfix()
                    calculate_variable_from_constraint(
                        blk.gas_emulsion.deltaP[t, x],
                        blk.gas_emulsion_pressure_drop[t, x])

        elif blk.config.has_pressure_change is False:
            blk.isobaric_gas_emulsion.activate()

            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    blk.gas_emulsion.properties[t, x].pressure.unfix()

        init_log.info('Initialize Mass Balances')
        init_log.info_high('initialize mass balances with no reactions')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
        if results.solver.termination_condition \
                == TerminationCondition.optimal:
            init_log.info_high(
                "Initialization Step 4a {}.".format(
                        idaeslog.condition(results))
                        )
        else:
            init_log.warning('{} Initialisation Step 4a Failed.'
                             .format(blk.name))

        # Homogeneous reactions (gas phase rxns)
        if gas_phase.reaction_package is not None:
            # local aliases used to shorten object names
            bubble_rxn_gen = blk.bubble.rate_reaction_generation
            bubble_stoichiometry_eqn = (
                blk.bubble.rate_reaction_stoichiometry_constraint)
            gas_emulsion_rxn_gen = blk.gas_emulsion.rate_reaction_generation
            gas_emulsion_stoichiometry_eqn = (
                blk.gas_emulsion.rate_reaction_stoichiometry_constraint)

            # Initialize, unfix and activate relevant model and CV1D
            # variables and constraints
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    for r in (gas_phase.reaction_package.rate_reaction_idx):
                        calculate_variable_from_constraint(
                            blk.bubble.rate_reaction_extent[t, x, r],
                            blk.bubble_rxn_ext_constraint[t, x, r])
                        calculate_variable_from_constraint(
                            blk.gas_emulsion.rate_reaction_extent[t, x, r],
                            blk.gas_emulsion_rxn_ext_constraint[t, x, r])
                    for j in gas_phase.property_package.component_list:
                        blk.bubble.rxn_generation[t, x, 'Vap', j].unfix()
                        calculate_variable_from_constraint(
                            bubble_rxn_gen[t, x, 'Vap', j],
                            bubble_stoichiometry_eqn[t, x, 'Vap', j])
                        gas_emulsion_rxn_gen[t, x, 'Vap', j].unfix()
                        calculate_variable_from_constraint(
                            gas_emulsion_rxn_gen[t, x, 'Vap', j],
                            gas_emulsion_stoichiometry_eqn[t, x, 'Vap', j])

            bubble_stoichiometry_eqn.activate()
            blk.bubble_rxn_ext_constraint.activate()
            gas_emulsion_stoichiometry_eqn.activate()
            blk.gas_emulsion_rxn_ext_constraint.activate()

            # Initialize homogeneous reaction property packages
            blk.bubble.reactions.activate()
            blk.gas_emulsion.reactions.activate()

            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    bubble = blk.bubble.reactions[t, x]
                    for c in bubble.component_objects(
                            Constraint, descend_into=False):
                        c.activate()
                    gas_emulsion = blk.gas_emulsion.reactions[t, x]
                    for c in gas_emulsion.component_objects(
                            Constraint, descend_into=False):
                        c.activate()

            blk.bubble.reactions.initialize(outlvl=outlvl,
                                            optarg=optarg,
                                            solver=solver)

            blk.gas_emulsion.reactions.initialize(outlvl=outlvl,
                                                  optarg=optarg,
                                                  solver=solver)

        # Heterogeneous rxns (solid phase rxns with gas phase interactions)
        if solid_phase.reaction_package is not None:
            # local aliases used to shorten object names
            solid_emulsion_rxn_gen = (
                blk.solid_emulsion.rate_reaction_generation)
            solid_emulsion_stoichiometry_eqn = (
                blk.solid_emulsion.rate_reaction_stoichiometry_constraint)

            # Initialize, unfix and activate relevant model and
            # CV1D variables and constraints
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    for j in gas_phase.property_package.component_list:
                        blk.gas_emulsion_hetero_rxn[t, x, j].unfix()
                        calculate_variable_from_constraint(
                            blk.gas_emulsion_hetero_rxn[t, x, j],
                            blk.gas_emulsion_hetero_rxn_eqn[t, x, j])
                for x in blk.length_domain:
                    for r in solid_phase.reaction_package.rate_reaction_idx:
                        calculate_variable_from_constraint(
                            blk.solid_emulsion.rate_reaction_extent[t, x, r],
                            blk.solid_emulsion_rxn_ext_constraint[t, x, r])
                    for j in solid_phase.property_package.component_list:
                        solid_emulsion_rxn_gen[t, x, 'Sol', j].unfix()
                        if not ((blk.config.transformation_scheme != "FORWARD"
                                 and x == blk.length_domain.first()) or
                                (blk.config.transformation_scheme == "FORWARD"
                                 and x == blk.length_domain.last())):
                            calculate_variable_from_constraint(
                                solid_emulsion_rxn_gen[t, x, 'Sol', j],
                                solid_emulsion_stoichiometry_eqn[t, x,
                                                                 'Sol', j])

            blk.solid_emulsion_rxn_ext_constraint.activate()
            blk.gas_emulsion_hetero_rxn_eqn.activate()
            solid_emulsion_stoichiometry_eqn.activate()

            # Initialize heterogeneous reaction property package
            blk.solid_emulsion.reactions.activate()
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    obj = blk.solid_emulsion.reactions[t, x]
                    for c in obj.component_objects(
                            Constraint, descend_into=False):
                        c.activate()

            blk.solid_emulsion.reactions.initialize(outlvl=outlvl,
                                                    optarg=optarg,
                                                    solver=solver)

        if (gas_phase.reaction_package is not None or
                solid_phase.reaction_package is not None):
            init_log.info_high('initialize mass balances with reactions')
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                init_log.info_high(
                    "Initialization Step 4b {}.".format(
                            idaeslog.condition(results))
                            )
            else:
                init_log.warning('{} Initialisation Step 4b Failed.'
                                 .format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize energy balance
        if blk.config.energy_balance_type != EnergyBalanceType.none:
            # Initialize relevant heat transfer variables
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    calculate_variable_from_constraint(
                        blk._reform_var_4[t, x],
                        blk._reformulation_eqn_4[t, x])
                    calculate_variable_from_constraint(
                        blk._reform_var_5[t, x],
                        blk._reformulation_eqn_5[t, x])
                    calculate_variable_from_constraint(
                        blk.Hbe[t, x],
                        blk.bubble_cloud_heat_trans_coeff[t, x])
                    calculate_variable_from_constraint(
                            blk.htc_conv[t, x],
                            blk.convective_heat_trans_coeff[t, x])
                    calculate_variable_from_constraint(
                        blk.ht_conv[t, x],
                        blk.convective_heat_transfer[t, x])

            # Unfix temperature variables
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    blk.bubble.properties[t, x].temperature.unfix()
                    blk.gas_emulsion.properties[t, x].temperature.unfix()
                    blk.solid_emulsion.properties[t, x].temperature.unfix()

            # Activate relevant model level constraints
            blk._reformulation_eqn_4.activate()
            blk._reformulation_eqn_5.activate()
            blk.bubble_cloud_heat_trans_coeff.activate()
            blk.convective_heat_trans_coeff.activate()
            blk.convective_heat_transfer.activate()
            blk.bubble_cloud_bulk_heat_trans.activate()
            blk.bubble_heat_transfer.activate()
            blk.gas_emulsion_heat_transfer.activate()
            blk.solid_emulsion_heat_transfer.activate()

            # Activate energy balance equations
            blk.bubble.enthalpy_balances.activate()
            blk.gas_emulsion.enthalpy_balances.activate()
            blk.solid_emulsion.enthalpy_balances.activate()

            # Activate energy balance boundary conditions
            blk.gas_energy_balance_in.activate()
            blk.gas_emulsion_temperature_in.activate()
            blk.solid_energy_balance_in.activate()

            init_log.info('Initialize Energy Balances')
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                init_log.info_high(
                    "Initialization Step 5 {}.".format(
                            idaeslog.condition(results))
                            )
            else:
                init_log.warning('{} Initialisation Step 5 Failed.'
                                 .format(blk.name))

        # Initialize energy balance
        if blk.config.energy_balance_type == EnergyBalanceType.none:
            for t in blk.flowsheet().config.time:
                for x in blk.length_domain:
                    blk.bubble.properties[t, x].temperature.unfix()
                    blk.gas_emulsion.properties[t, x].temperature.unfix()
                    blk.solid_emulsion.properties[t, x].temperature.unfix()

            blk.isothermal_gas_emulsion.activate()
            blk.isothermal_bubble.activate()
            blk.isothermal_solid_emulsion.activate()

            init_log.info('Initialize Energy Balances')
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                init_log.info_high(
                    "Initialization Step 5 {}.".format(
                            idaeslog.condition(results))
                            )
            else:
                init_log.warning('{} Initialisation Step 5 Failed.'
                                 .format(blk.name))
        # ---------------------------------------------------------------------
        # Initialize outlet conditions
        # Initialize gas_outlet block
        blk.gas_outlet_block.activate()
        blk.solid_outlet_block.activate()

        # Activate outlet boundary conditions
        blk.gas_pressure_out.activate()
        blk.gas_material_balance_out.activate()
        blk.solid_material_balance_out.activate()
        blk.solid_particle_porosity_out.activate()

        blk.gas_energy_balance_out.activate()
        blk.solid_energy_balance_out.activate()

        blk.gas_outlet_block.initialize(
            state_args=gas_phase_state_args,
            hold_state=False,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver)
        blk.solid_outlet_block.initialize(
            state_args=solid_phase_state_args,
            hold_state=False,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver)

        init_log.info('Initialize Outlet Conditions')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
        if results.solver.termination_condition \
                == TerminationCondition.optimal:
            init_log.info_high(
                "Initialization Step 6 {}.".format(
                        idaeslog.condition(results))
                        )
        else:
            init_log.warning('{} Initialisation Step 6 Failed.'
                             .format(blk.name))

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

    def results_plot(blk):
        '''
        Plot method for common bubbling fluidized bed variables

        Variables plotted:
            Tge : temperature of gas in the emulsion region
            Tgb : temperature of gas in the bubble region
            Tse : temperature of solid in the emulsion region
            Ge : flowrate of gas in the emulsion region
            Gb : flowrate of gas in the bubble region
            cet : total concentration of gas in the emulsion region
            cbt : total concentration of gas in the bubble region
            y_b : mole fraction of gas components in the bubble region
            x_e : mass fraction of solid components in the emulsion region
        '''
        print()
        print('================================= Reactor plots ==============='
              '==================')
        # local aliases used to shorten object names
        gas_phase = blk.config.gas_phase_config
        solid_phase = blk.config.solid_phase_config

        Tge = []
        Tgb = []
        Tse = []
        Ge = []
        Gb = []
        cet = []
        cbt = []

        for t in blk.flowsheet().config.time:
            for x in blk.gas_emulsion.length_domain:
                Tge.append(value(
                        blk.gas_emulsion.properties[t, x].temperature))
                Tgb.append(value(
                        blk.bubble.properties[t, x].temperature))
                Tse.append(value(
                        blk.solid_emulsion.properties[t, x].temperature))
                Ge.append(value(blk.gas_emulsion.properties[t, x].flow_mol))
                Gb.append(value(blk.bubble.properties[t, x].flow_mol))
                cet.append(value(
                        blk.gas_emulsion.properties[t, x].dens_mol))
                cbt.append(value(
                        blk.bubble.properties[t, x].dens_mol))

        # Bed temperature profile
        plt.figure(1)
        plt.plot(blk.gas_emulsion.length_domain, Tge, label='Tge')
        plt.plot(blk.gas_emulsion.length_domain, Tgb, label='Tgb')

        plt.legend(loc=9, ncol=2)
        plt.grid()
        plt.xlabel("Bed height")
        plt.ylabel("Gas temperatures in bed regions (K)")

        plt.figure(2)
        plt.plot(blk.gas_emulsion.length_domain, Tse, label='Tse')

        plt.legend(loc=9, ncol=3)
        plt.grid()
        plt.xlabel("Bed height")
        plt.ylabel("Solid temperatures in bed regions (K)")

        plt.figure(3)
        plt.plot(blk.gas_emulsion.length_domain, Ge, label='Ge')
        plt.plot(blk.gas_emulsion.length_domain, Gb, label='Gb')

        plt.legend(loc=9, ncol=2)
        plt.grid()
        plt.xlabel("Bed height")
        plt.ylabel("Gas flow (mol/s)")

        plt.figure(4)
        plt.plot(blk.gas_emulsion.length_domain, cbt, label='cbt')
        plt.plot(blk.gas_emulsion.length_domain, cet, label='cet')

        plt.legend(loc=9, ncol=3)
        plt.grid()
        plt.xlabel("Bed height")
        plt.ylabel("gas conc. mol/s")

        # Gas phase mole composition
        for t in blk.flowsheet().config.time:
            for i in (gas_phase.property_package.component_list):
                y_b = []
                for x in blk.gas_emulsion.length_domain:
                    y_b.append(value(
                        blk.bubble.properties[t, x].mole_frac_comp[i]))
                plt.figure(5)
                plt.plot(blk.gas_emulsion.length_domain, y_b, label=i)
        plt.legend(loc=9, ncol=len(
                gas_phase.property_package.component_list))
        plt.grid()
        plt.xlabel("Bed height")
        plt.ylabel("Gas bubble mole frac. (-)")

        # Solid phase mass composition
        for t in blk.flowsheet().config.time:
            for i in (solid_phase.property_package.component_list):
                x_e = []
                for x in blk.solid_emulsion.length_domain:
                    x_e.append(value(
                        blk.solid_emulsion.properties[t, x].mass_frac_comp[i]))
                plt.figure(6)
                plt.plot(blk.solid_emulsion.length_domain, x_e, label=i)
        plt.legend(loc=9, ncol=len(
                solid_phase.property_package.component_list))
        plt.grid()
        plt.xlabel("Bed height")
        plt.ylabel("Solid emulsion mass frac. (-)")
