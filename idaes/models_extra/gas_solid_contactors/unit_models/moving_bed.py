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
IDAES Moving Bed Model.

The moving bed model is a 1D axially discretized model with a gas and
solid phase and a counter-current flow direction. The model captures
the gas-solid interaction between both phases through reaction, mass
and heat transfer.

Equations written in this model were derived from:
(1) A. Ostace, A. Lee, C.O. Okoli, A.P. Burgard, D.C. Miller, D. Bhattacharyya,
Mathematical modeling of a moving-bed reactor for chemical looping combustion
of methane, in: M.R. Eden, M. Ierapetritou, G.P. Towler (Eds.),13th Int. Symp.
Process Syst. Eng. (PSE 2018), Computer-Aided Chemical Engineering 2018,
pp. 325–330 , San Diego, CA.
(2) C.O. Okoli, A. Ostace, S. Nadgouda, A. Lee, A. Tong, A.P. Burgard,
D. Bhattacharyya, D.C. Miller, 2020.
A framework for the optimization of chemical looping combustion processes,
Journal of Powder Technology, Vol. 365, pp 149 - 162.

Assumptions:
Property package contains temperature and pressure variables.
Property package contains minimum fluidization velocity.

"""
from __future__ import division

# Import Python libraries
import matplotlib.pyplot as plt

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    Reals,
    value,
    TransformationFactory,
    Constraint,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.dae import ContinuousSet

# Import IDAES cores
from idaes.core import (
    ControlVolume1DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
    DistributedVars,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    BurntToast,
    InitializationError,
)
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants as constants
from idaes.core.util.math import smooth_abs
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

__author__ = "Chinedu Okoli", "Anca Ostace"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("MBR")
class MBRData(UnitModelBlockData):
    """Standard Moving Bed Unit Model Class."""

    # Create template for unit level config arguments
    CONFIG = UnitModelBlockData.CONFIG()

    # Unit level config arguments
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=10,
            domain=int,
            description="Number of finite elements length domain",
            doc="""Number of finite elements to use when discretizing length
domain (default=20)""",
        ),
    )
    CONFIG.declare(
        "length_domain_set",
        ConfigValue(
            default=[0.0, 1.0],
            domain=list,
            description="Number of finite elements length domain",
            doc="""length_domain_set - (optional) list of point to use to
initialize a new ContinuousSet if length_domain is not
provided (default = [0.0, 1.0])""",
        ),
    )
    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default="dae.finite_difference",
            description="Method to use for DAE transformation",
            doc="""Method to use to transform domain. Must be a method recognised
by the Pyomo TransformationFactory,
**default** - "dae.finite_difference".
**Valid values:** {
**"dae.finite_difference"** - Use a finite difference transformation method,
**"dae.collocation"** - use a collocation transformation method}""",
        ),
    )
    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
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
**"LAGRANGE-RADAU""** - use a collocation transformation method}""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=3,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)""",
        ),
    )
    CONFIG.declare(
        "flow_type",
        ConfigValue(
            default="counter_current",
            domain=In(["counter_current"]),
            description="Flow configuration of Moving Bed",
            doc="""Flow configuration of Moving Bed
- counter_current: gas side flows from 0 to 1
solid side flows from 1 to 0""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
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
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
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
        "has_pressure_change",
        ConfigValue(
            default=True,
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
        "pressure_drop_type",
        ConfigValue(
            default="simple_correlation",
            domain=In(["simple_correlation", "ergun_correlation"]),
            description="Construction flag for type of pressure drop",
            doc="""Indicates what type of pressure drop correlation should be used,
**default** - "simple_correlation".
**Valid values:** {
**"simple_correlation"** - Use a simplified pressure drop correlation,
**"ergun_correlation"** - Use the Ergun equation.}""",
        ),
    )

    # Create template for phase specific config arguments
    _PhaseTemplate = UnitModelBlockData.CONFIG()
    _PhaseTemplate.declare(
        "has_equilibrium_reactions",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Equilibrium reaction construction flag",
            doc="""Indicates whether terms for equilibrium controlled reactions
should be constructed,
**default** - True.
**Valid values:** {
**True** - include equilibrium reaction terms,
**False** - exclude equilibrium reaction terms.}""",
        ),
    )
    _PhaseTemplate.declare(
        "property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object""",
        ),
    )
    _PhaseTemplate.declare(
        "property_package_args",
        ConfigValue(
            default={},
            domain=dict,
            description="Arguments for constructing gas property package",
            doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)""",
        ),
    )
    _PhaseTemplate.declare(
        "reaction_package",
        ConfigValue(
            default=None,
            domain=is_reaction_parameter_block,
            description="Reaction package to use for control volume",
            doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
        ),
    )
    _PhaseTemplate.declare(
        "reaction_package_args",
        ConfigBlock(
            implicit=True,
            implicit_domain=ConfigBlock,
            description="Arguments to use for constructing reaction packages",
            doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
        ),
    )

    # Create individual config blocks for gas and solid sides
    CONFIG.declare("gas_phase_config", _PhaseTemplate(doc="gas phase config arguments"))
    CONFIG.declare(
        "solid_phase_config", _PhaseTemplate(doc="solid phase config arguments")
    )

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
        super().build()

        # Consistency check for transformation method and transformation scheme
        if (
            self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme is None
        ):
            self.config.transformation_scheme = "BACKWARD"
        elif (
            self.config.transformation_method == "dae.collocation"
            and self.config.transformation_scheme is None
        ):
            self.config.transformation_scheme = "LAGRANGE-RADAU"
        elif (
            self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme != "BACKWARD"
            and self.config.transformation_scheme != "FORWARD"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_scheme argument. "
                "Must be "
                "BACKWARD"
                " or "
                "FORWARD"
                " "
                "if transformation_method is"
                " "
                "dae.finite_difference"
                ".".format(self.name)
            )
        elif (
            self.config.transformation_method == "dae.collocation"
            and self.config.transformation_scheme != "LAGRANGE-RADAU"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_scheme argument."
                "Must be "
                "LAGRANGE-RADAU"
                " if "
                "transformation_method is"
                " "
                "dae.collocation"
                ".".format(self.name)
            )

        # Set flow directions for the control volume blocks
        # Gas flows from 0 to 1, solid flows from 1 to 0
        # An if statement is used here despite only one option to allow for
        # future extensions to other flow configurations
        if self.config.flow_type == "counter_current":
            set_direction_gas = FlowDirection.forward
            set_direction_solid = FlowDirection.backward
        else:
            raise BurntToast(
                "{} encountered unrecognized argument "
                "for flow type. Please contact the IDAES"
                " developers with this bug.".format(self.name)
            )

        # Set arguments for gas sides if homoogeneous reaction block
        if self.config.gas_phase_config.reaction_package is not None:
            has_rate_reaction_gas_phase = True
        else:
            has_rate_reaction_gas_phase = False

        # Set arguments for gas and solid sides if heterogeneous reaction block
        if self.config.solid_phase_config.reaction_package is not None:
            has_rate_reaction_solid_phase = True
            has_mass_transfer_gas_phase = True
        else:
            has_rate_reaction_solid_phase = False
            has_mass_transfer_gas_phase = False

        # Set heat transfer terms
        if self.config.energy_balance_type != EnergyBalanceType.none:
            has_heat_transfer = True
        else:
            has_heat_transfer = False

        # Set heat of reaction terms
        if (
            self.config.energy_balance_type != EnergyBalanceType.none
            and self.config.gas_phase_config.reaction_package is not None
        ):
            has_heat_of_reaction_gas_phase = True
        else:
            has_heat_of_reaction_gas_phase = False

        if (
            self.config.energy_balance_type != EnergyBalanceType.none
            and self.config.solid_phase_config.reaction_package is not None
        ):
            has_heat_of_reaction_solid_phase = True
        else:
            has_heat_of_reaction_solid_phase = False

        # local aliases used to shorten object names
        solid_phase = self.config.solid_phase_config

        # Get units meta data from property packages
        units_meta_solid = solid_phase.property_package.get_metadata().get_derived_units

        # Create a unit model length domain
        self.length_domain = ContinuousSet(
            bounds=(0.0, 1.0),
            initialize=self.config.length_domain_set,
            doc="Normalized length domain",
        )

        self.bed_height = Var(
            domain=Reals,
            initialize=1,
            doc="Bed length",
            units=units_meta_solid("length"),
        )

        # =========================================================================
        """ Build Control volume 1D for gas phase and
            populate gas control volume"""

        self.gas_phase = ControlVolume1DBlock(
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            area_definition=DistributedVars.variant,
            property_package=self.config.gas_phase_config.property_package,
            property_package_args=self.config.gas_phase_config.property_package_args,
            reaction_package=self.config.gas_phase_config.reaction_package,
            reaction_package_args=self.config.gas_phase_config.reaction_package_args,
        )

        self.gas_phase.add_geometry(
            length_domain=self.length_domain,
            length_domain_set=self.config.length_domain_set,
            length_var=self.bed_height,
            flow_direction=set_direction_gas,
        )

        self.gas_phase.add_state_blocks(
            information_flow=set_direction_gas, has_phase_equilibrium=False
        )

        if self.config.gas_phase_config.reaction_package is not None:
            self.gas_phase.add_reaction_blocks(
                has_equilibrium=self.config.gas_phase_config.has_equilibrium_reactions
            )

        self.gas_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
            has_mass_transfer=has_mass_transfer_gas_phase,
            has_rate_reactions=has_rate_reaction_gas_phase,
        )

        self.gas_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=has_heat_transfer,
            has_heat_of_reaction=has_heat_of_reaction_gas_phase,
        )

        self.gas_phase.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # =========================================================================
        """ Build Control volume 1D for solid phase and
            populate solid control volume"""

        # Set argument for heterogeneous reaction block
        self.solid_phase = ControlVolume1DBlock(
            transformation_method=self.config.transformation_method,
            transformation_scheme=self.config.transformation_scheme,
            finite_elements=self.config.finite_elements,
            collocation_points=self.config.collocation_points,
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            area_definition=DistributedVars.variant,
            property_package=self.config.solid_phase_config.property_package,
            property_package_args=self.config.solid_phase_config.property_package_args,
            reaction_package=self.config.solid_phase_config.reaction_package,
            reaction_package_args=self.config.solid_phase_config.reaction_package_args,
        )

        self.solid_phase.add_geometry(
            length_domain=self.length_domain,
            length_domain_set=self.config.length_domain_set,
            length_var=self.bed_height,
            flow_direction=set_direction_solid,
        )

        self.solid_phase.add_state_blocks(
            information_flow=set_direction_solid, has_phase_equilibrium=False
        )

        if self.config.solid_phase_config.reaction_package is not None:
            # TODO - a generalization of the heterogeneous reaction block
            # The heterogeneous reaction block does not use the
            # add_reaction_blocks in control volumes as control volumes are
            # currently setup to handle only homogeneous reaction properties.
            # Thus appending the heterogeneous reaction block to the
            # solid state block is currently hard coded here.

            tmp_dict = dict(**self.config.solid_phase_config.reaction_package_args)
            tmp_dict["gas_state_block"] = self.gas_phase.properties
            tmp_dict["solid_state_block"] = self.solid_phase.properties
            tmp_dict[
                "has_equilibrium"
            ] = self.config.solid_phase_config.has_equilibrium_reactions
            tmp_dict["parameters"] = self.config.solid_phase_config.reaction_package
            self.solid_phase.reactions = (
                self.config.solid_phase_config.reaction_package.reaction_block_class(
                    self.flowsheet().time,
                    self.length_domain,
                    doc="Reaction properties in control volume",
                    **tmp_dict,
                )
            )

        self.solid_phase.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=False,
            has_mass_transfer=False,
            has_rate_reactions=has_rate_reaction_solid_phase,
        )

        self.solid_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=has_heat_transfer,
            has_heat_of_reaction=has_heat_of_reaction_solid_phase,
        )

        self.solid_phase.add_momentum_balances(
            balance_type=MomentumBalanceType.none, has_pressure_change=False
        )

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
                " configuration argument. Please provide a valid value.".format(
                    self.name
                )
            )

        if self.config.transformation_method == "dae.finite_difference":
            self.discretizer = TransformationFactory(self.config.transformation_method)
            self.discretizer.apply_to(
                self,
                wrt=self.length_domain,
                nfe=self.config.finite_elements,
                scheme=self.config.transformation_scheme,
            )
        elif self.config.transformation_method == "dae.collocation":
            self.discretizer = TransformationFactory(self.config.transformation_method)
            self.discretizer.apply_to(
                self,
                wrt=self.length_domain,
                nfe=self.config.finite_elements,
                ncp=self.config.collocation_points,
                scheme=self.config.transformation_scheme,
            )

    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """
        # local aliases used to shorten object names
        gas_phase = self.config.gas_phase_config
        solid_phase = self.config.solid_phase_config

        # Get units meta data from property packages
        units_meta_gas = gas_phase.property_package.get_metadata().get_derived_units
        units_meta_solid = solid_phase.property_package.get_metadata().get_derived_units

        # Declare Mutable Parameters
        self.eps = Param(
            mutable=True,
            default=1e-8,
            doc="Smoothing Factor for Smooth IF Statements",
            units=pyunits.dimensionless,
        )

        # Unit Model variables
        self.bed_diameter = Var(
            initialize=1, doc="Reactor diameter", units=units_meta_solid("length")
        )
        self.bed_area = Var(
            domain=Reals,
            initialize=1,
            doc="Reactor cross-sectional area",
            units=units_meta_solid("area"),
        )

        # Phase specific variables
        self.velocity_superficial_gas = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=0.05,
            doc="Gas superficial velocity",
            units=units_meta_gas("velocity"),
        )
        self.velocity_superficial_solid = Var(
            self.flowsheet().time,
            domain=Reals,
            initialize=0.005,
            doc="Solid superficial velocity",
            units=units_meta_solid("velocity"),
        )

        # Dimensionless numbers, mass and heat transfer coefficients
        self.Re_particle = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=1.0,
            doc="Particle Reynolds number",
            units=pyunits.dimensionless,
        )

        self.Pr_particle = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=1.0,
            doc="Prandtl number of gas in bed",
            units=pyunits.dimensionless,
        )

        self.Nu_particle = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=1.0,
            doc="Particle Nusselt number",
            units=pyunits.dimensionless,
        )
        self.gas_solid_htc = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=1.0,
            doc="Gas-solid heat transfer coefficient",
            units=units_meta_gas("heat_transfer_coefficient"),
        )

        # Fixed variables (these are parameters that can be estimated)
        self.bed_voidage = Var(
            domain=Reals, initialize=0.4, doc="Bed voidage", units=pyunits.dimensionless
        )
        self.bed_voidage.fix()

        # =========================================================================
        # Add performance equations

        # ---------------------------------------------------------------------
        # Geometry contraints

        # Bed area
        @self.Constraint(doc="Bed area")
        def bed_area_eqn(b):
            return b.bed_area == (constants.pi * (0.5 * b.bed_diameter) ** 2)

        # Area of gas side, and solid side
        @self.Constraint(self.flowsheet().time, self.length_domain, doc="Gas side area")
        def gas_phase_area(b, t, x):
            return (
                b.gas_phase.area[t, x]
                == pyunits.convert(b.bed_area, to_units=units_meta_gas("area"))
                * b.bed_voidage
            )

        @self.Constraint(
            self.flowsheet().time, self.length_domain, doc="Solid side area"
        )
        def solid_phase_area(b, t, x):
            return b.solid_phase.area[t, x] == b.bed_area * (1 - b.bed_voidage)

        # ---------------------------------------------------------------------
        # Hydrodynamic contraints

        # Gas superficial velocity
        @self.Constraint(
            self.flowsheet().time, self.length_domain, doc="Gas superficial velocity"
        )
        def gas_super_vel(b, t, x):
            return (
                b.velocity_superficial_gas[t, x]
                * pyunits.convert(b.bed_area, to_units=units_meta_gas("area"))
                * b.gas_phase.properties[t, x].dens_mol
                == b.gas_phase.properties[t, x].flow_mol
            )

        # Solid superficial velocity
        @self.Constraint(
            self.flowsheet().time, self.length_domain, doc="Solid superficial velocity"
        )
        # This equation uses inlet values to compute the constant solid
        # superficial velocity, and then computes the solid particle density
        # through the rest of the bed.
        def solid_super_vel(b, t, x):
            return (
                b.velocity_superficial_solid[t]
                * b.bed_area
                * b.solid_phase.properties[t, x].dens_mass_particle
                == b.solid_phase.properties[t, x].flow_mass
            )

        # Gas side pressure drop calculation
        if (
            self.config.has_pressure_change
            and self.config.pressure_drop_type == "simple_correlation"
        ):
            # Simplified pressure drop
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Gas side pressure drop calculation -" "simplified pressure drop",
            )
            def gas_phase_config_pressure_drop(b, t, x):
                #  0.2/s is a unitted constant in the correlation
                return pyunits.convert(
                    b.gas_phase.deltaP[t, x],
                    to_units=units_meta_solid("pressure") / units_meta_solid("length"),
                ) == -(0.2 / pyunits.s) * (
                    b.velocity_superficial_gas[t, x]
                    * (
                        pyunits.convert(
                            b.solid_phase.properties[t, x].dens_mass_particle,
                            to_units=units_meta_solid("density_mass"),
                        )
                        - b.gas_phase.properties[t, x].dens_mass
                    )
                )

        elif (
            self.config.has_pressure_change
            and self.config.pressure_drop_type == "ergun_correlation"
        ):
            # Ergun equation
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Gas side pressure drop calculation -" "Ergun equation",
            )
            def gas_phase_config_pressure_drop(b, t, x):
                return -pyunits.convert(
                    b.gas_phase.deltaP[t, x],
                    to_units=units_meta_solid("pressure") / units_meta_solid("length"),
                ) == (
                    (150 * pyunits.dimensionless)
                    * (1 - b.bed_voidage) ** 2
                    * b.gas_phase.properties[t, x].visc_d
                    * (
                        b.velocity_superficial_gas[t, x]
                        + pyunits.convert(
                            b.velocity_superficial_solid[t],
                            to_units=units_meta_solid("velocity"),
                        )
                    )
                    / (
                        pyunits.convert(
                            b.solid_phase.properties[t, x]._params.particle_dia,
                            to_units=units_meta_solid("length"),
                        )
                        ** 2
                        * b.bed_voidage**3
                    )
                ) + (
                    (1.75 * pyunits.dimensionless)
                    * b.gas_phase.properties[t, x].dens_mass
                    * (1 - b.bed_voidage)
                    * (
                        b.velocity_superficial_gas[t, x]
                        + pyunits.convert(
                            b.velocity_superficial_solid[t],
                            to_units=units_meta_solid("velocity"),
                        )
                    )
                    ** 2
                    / (
                        pyunits.convert(
                            b.solid_phase.properties[t, x]._params.particle_dia,
                            to_units=units_meta_solid("length"),
                        )
                        * b.bed_voidage**3
                    )
                )

            # to_units=deltaP_units)
            # The above expression has no absolute values - assumes:
            # (velocity_superficial_gas + velocity_superficial_solid) > 0

        elif self.config.has_pressure_change is False:
            pass

        else:
            raise BurntToast(
                "{} encountered unrecognized argument for "
                "the pressure drop correlation. Please contact the IDAES"
                " developers with this bug.".format(self.name)
            )
        # ---------------------------------------------------------------------
        # Reaction contraints

        # Build homogeneous reaction constraints
        if gas_phase.reaction_package is not None:
            # Gas side rate reaction extent
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                gas_phase.reaction_package.rate_reaction_idx,
                doc="Gas side rate reaction extent",
            )
            def gas_phase_config_rxn_ext(b, t, x, r):
                return b.gas_phase.rate_reaction_extent[t, x, r] == (
                    b.gas_phase.reactions[t, x].reaction_rate[r]
                    * b.gas_phase.area[t, x]
                )

        # Build hetereogeneous reaction constraints
        if solid_phase.reaction_package is not None:
            # Solid side rate reaction extent
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                solid_phase.reaction_package.rate_reaction_idx,
                doc="Solid side rate reaction extent",
            )
            def solid_phase_config_rxn_ext(b, t, x, r):
                return b.solid_phase.rate_reaction_extent[t, x, r] == (
                    b.solid_phase.reactions[t, x].reaction_rate[r]
                    * b.solid_phase.area[t, x]
                )

            # Gas side heterogeneous rate reaction generation
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                gas_phase.property_package.phase_list,
                gas_phase.property_package.component_list,
                doc="Gas side heterogeneous" "rate reaction generation",
            )
            def gas_comp_hetero_rxn(b, t, x, p, j):
                return b.gas_phase.mass_transfer_term[t, x, p, j] == (
                    sum(
                        b.solid_phase.reactions[t, x].rate_reaction_stoichiometry[
                            r, p, j
                        ]
                        * b.solid_phase.reactions[t, x].reaction_rate[r]
                        for r in (solid_phase.reaction_package.rate_reaction_idx)
                    )
                    * b.solid_phase.area[t, x]
                )

        # ---------------------------------------------------------------------
        if self.config.energy_balance_type != EnergyBalanceType.none:
            # Solid phase - gas to solid heat transfer
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Solid phase - gas to solid heat transfer",
            )
            def solid_phase_heat_transfer(b, t, x):
                return pyunits.convert(
                    b.solid_phase.heat[t, x],
                    to_units=units_meta_gas("power") / units_meta_gas("length"),
                ) * pyunits.convert(
                    b.solid_phase.properties[t, x]._params.particle_dia,
                    to_units=units_meta_gas("length"),
                ) == 6 * b.gas_solid_htc[
                    t, x
                ] * (
                    b.gas_phase.properties[t, x].temperature
                    - pyunits.convert(
                        b.solid_phase.properties[t, x].temperature,
                        to_units=units_meta_gas("temperature"),
                    )
                ) * pyunits.convert(
                    b.solid_phase.area[t, x], to_units=units_meta_gas("area")
                )

            # Dimensionless numbers, mass and heat transfer coefficients
            # Particle Reynolds number
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Particle Reynolds number",
            )
            def reynolds_number_particle(b, t, x):
                return (
                    b.Re_particle[t, x] * b.gas_phase.properties[t, x].visc_d
                    == b.velocity_superficial_gas[t, x]
                    * pyunits.convert(
                        b.solid_phase.properties[t, x]._params.particle_dia,
                        to_units=units_meta_gas("length"),
                    )
                    * b.gas_phase.properties[t, x].dens_mass
                )

            # Prandtl number
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Prandtl number of gas in bed",
            )
            def prandtl_number(b, t, x):
                return (
                    b.Pr_particle[t, x] * b.gas_phase.properties[t, x].therm_cond
                    == b.gas_phase.properties[t, x].cp_mass
                    * b.gas_phase.properties[t, x].visc_d
                )

            # Particle Nusselt number
            @self.Constraint(
                self.flowsheet().time, self.length_domain, doc="Particle Nusselt number"
            )
            def nusselt_number_particle(b, t, x):
                return (
                    b.Nu_particle[t, x] ** 3
                    == (
                        2.0 + 1.1 * (smooth_abs(b.Re_particle[t, x], b.eps) ** 0.6) ** 3
                    )
                    * b.Pr_particle[t, x]
                )

            # Gas-solid heat transfer coefficient
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Gas-solid heat transfer coefficient",
            )
            def gas_solid_htc_eqn(b, t, x):
                return (
                    b.gas_solid_htc[t, x]
                    * pyunits.convert(
                        b.solid_phase.properties[t, x]._params.particle_dia,
                        to_units=units_meta_gas("length"),
                    )
                    == b.Nu_particle[t, x] * b.gas_phase.properties[t, x].therm_cond
                )

            # Gas phase - gas to solid heat transfer
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Gas phase - gas to solid heat transfer",
            )
            def gas_phase_heat_transfer(b, t, x):
                return b.gas_phase.heat[t, x] * pyunits.convert(
                    b.solid_phase.properties[t, x]._params.particle_dia,
                    to_units=units_meta_gas("length"),
                ) == -6 * b.gas_solid_htc[t, x] * (
                    b.gas_phase.properties[t, x].temperature
                    - pyunits.convert(
                        b.solid_phase.properties[t, x].temperature,
                        to_units=units_meta_gas("temperature"),
                    )
                ) * pyunits.convert(
                    b.solid_phase.area[t, x], to_units=units_meta_gas("area")
                )

        elif self.config.energy_balance_type == EnergyBalanceType.none:
            # If energy balance is none fix gas and solid temperatures to
            # solid inlet temperature
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Isothermal gas phase constraint",
            )
            def isothermal_gas_phase(b, t, x):
                if x == self.length_domain.first():
                    return Constraint.Skip
                else:
                    return (
                        b.gas_phase.properties[t, x].temperature
                        == b.solid_inlet.temperature[t]
                    )

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Isothermal solid phase constraint",
            )
            def isothermal_solid_phase(b, t, x):
                if x == self.length_domain.last():
                    return Constraint.Skip
                else:
                    return (
                        b.solid_phase.properties[t, x].temperature
                        == b.solid_inlet.temperature[t]
                    )

    # =========================================================================
    # Model initialization routine

    def initialize_build(
        blk,
        gas_phase_state_args=None,
        solid_phase_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for MB unit.

        Keyword Arguments:
            gas_phase_state_args : a dict of arguments to be passed to the
                        property package(s) to provide an initial state for
                        initialization (see documentation of the specific
                        property package) (default = None).
            solid_phase_state_args : a dict of arguments to be passed to the
                        property package(s) to provide an initial state for
                        initialization (see documentation of the specific
                        property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """

        # Set up logger for initialization and solve
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

        geometry_constraints_terms = [
            "bed_area_eqn",
            "solid_phase_area",
            "gas_phase_area",
            "gas_phase_length",
            "solid_phase_length",
        ]
        endswith_terms = (
            "_disc_eq",
            "linking_constraint",
            "linking_constraints",
            "_holdup_calculation",
        )
        startswith_terms = "properties"

        for c in blk.component_objects(Constraint, descend_into=True):
            if (
                not c.parent_block().local_name.startswith(startswith_terms)
                and not c.local_name.endswith(endswith_terms)
                and c.local_name not in geometry_constraints_terms
            ):
                c.deactivate()

        # ---------------------------------------------------------------------
        # Initialize thermophysical property constraints
        init_log.info("Initialize Thermophysical Properties")
        # Initialize gas_phase block
        gas_phase_flags = blk.gas_phase.properties.initialize(
            state_args=gas_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        # Initialize solid_phase properties block
        solid_phase_flags = blk.solid_phase.properties.initialize(
            state_args=solid_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Initialize hydrodynamics (gas velocity)
        for t in blk.flowsheet().time:
            for x in blk.length_domain:
                calculate_variable_from_constraint(
                    blk.velocity_superficial_gas[t, x], blk.gas_super_vel[t, x]
                )

        blk.gas_super_vel.activate()

        init_log.info("Initialize Hydrodynamics")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
        if check_optimal_termination(results):
            init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(results))
            )
        else:
            _log.warning("{} Initialization Step 2 Failed.".format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize mass balance - no reaction and no pressure drop

        # Unfix material balance state variables (including particle porosity)
        # but keep other states fixed
        blk.gas_phase.properties.release_state(gas_phase_flags)
        blk.solid_phase.properties.release_state(solid_phase_flags)
        for t in blk.flowsheet().time:
            for x in blk.length_domain:
                blk.gas_phase.properties[t, x].pressure.fix()
                blk.gas_phase.properties[t, x].temperature.fix()
                blk.solid_phase.properties[t, x].temperature.fix()

        blk.gas_phase.material_balances.activate()

        if gas_phase.reaction_package is not None:
            for t in blk.flowsheet().time:
                gas_rxn_gen = blk.gas_phase.rate_reaction_generation
                for x in blk.length_domain:
                    for p in gas_phase.property_package.phase_list:
                        for j in gas_phase.property_package.component_list:
                            (gas_rxn_gen[t, x, p, j].fix(0.0))

        blk.solid_phase.material_balances.activate()
        blk.solid_super_vel.activate()

        if solid_phase.reaction_package is not None:
            for t in blk.flowsheet().time:
                solid_rxn_gen = blk.solid_phase.rate_reaction_generation
                for x in blk.length_domain:
                    for p in solid_phase.property_package.phase_list:
                        for j in solid_phase.property_package.component_list:
                            (solid_rxn_gen[t, x, p, j].fix(0.0))

                # Gas side heterogeneous rate reaction generation
                for x in blk.length_domain:
                    for p in gas_phase.property_package.phase_list:
                        for j in gas_phase.property_package.component_list:
                            (blk.gas_phase.mass_transfer_term[t, x, p, j].fix(0.0))

        init_log.info("Initialize Mass Balances")
        init_log.info_high(
            "initialize mass balances - no reactions " "and no pressure drop"
        )
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)
        if check_optimal_termination(results):
            init_log.info_high(
                "Initialization Step 3a {}.".format(idaeslog.condition(results))
            )
        else:
            _log.warning("{} Initialization Step 3a Failed.".format(blk.name))

        # Initialize mass balance - with reaction and no pressure drop
        if gas_phase.reaction_package is not None:
            # local aliases used to shorten object names
            gas_rxn_gen = blk.gas_phase.rate_reaction_generation
            gas_phase_stoichiometry_eqn = (
                blk.gas_phase.rate_reaction_stoichiometry_constraint
            )

            # Initialize reaction property package
            blk.gas_phase.reactions.activate()
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    obj = blk.gas_phase.reactions[t, x]
                    for c in obj.component_objects(Constraint, descend_into=False):
                        c.activate()

            blk.gas_phase.reactions.initialize(
                outlvl=outlvl, optarg=optarg, solver=solver
            )

            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    for r in gas_phase.reaction_package.rate_reaction_idx:
                        calculate_variable_from_constraint(
                            blk.gas_phase.rate_reaction_extent[t, x, r],
                            blk.gas_phase_config_rxn_ext[t, x, r],
                        )
                    for p in gas_phase.property_package.phase_list:
                        for j in gas_phase.property_package.component_list:
                            (gas_rxn_gen[t, x, p, j].unfix())
                            if not (
                                (
                                    blk.gas_phase.config.transformation_scheme
                                    != "FORWARD"
                                    and x == blk.length_domain.first()
                                )
                                or (
                                    blk.gas_phase.config.transformation_scheme
                                    == "FORWARD"
                                    and x == blk.length_domain.last()
                                )
                            ):
                                calculate_variable_from_constraint(
                                    gas_rxn_gen[t, x, p, j],
                                    gas_phase_stoichiometry_eqn[t, x, p, j],
                                )

            gas_phase_stoichiometry_eqn.activate()
            blk.gas_phase_config_rxn_ext.activate()

        if solid_phase.reaction_package is not None:
            # local aliases used to shorten object names
            solid_rxn_gen = blk.solid_phase.rate_reaction_generation
            solid_phase_stoichiometry_eqn = (
                blk.solid_phase.rate_reaction_stoichiometry_constraint
            )
            gas_mass_transfer_term = blk.gas_phase.mass_transfer_term

            # Initialize reaction property package
            blk.solid_phase.reactions.activate()
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    obj = blk.solid_phase.reactions[t, x]
                    for c in obj.component_objects(Constraint, descend_into=False):
                        c.activate()

            blk.solid_phase.reactions.initialize(
                outlvl=outlvl, optarg=optarg, solver=solver
            )

            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    for p in gas_phase.property_package.phase_list:
                        for j in gas_phase.property_package.component_list:
                            (gas_mass_transfer_term[t, x, p, j].unfix())
                            calculate_variable_from_constraint(
                                gas_mass_transfer_term[t, x, p, j],
                                blk.gas_comp_hetero_rxn[t, x, p, j],
                            )
                for x in blk.length_domain:
                    for r in solid_phase.reaction_package.rate_reaction_idx:
                        calculate_variable_from_constraint(
                            blk.solid_phase.rate_reaction_extent[t, x, r],
                            blk.solid_phase_config_rxn_ext[t, x, r],
                        )
                    for p in solid_phase.property_package.phase_list:
                        for j in solid_phase.property_package.component_list:
                            (solid_rxn_gen[t, x, p, j].unfix())
                            if not (
                                (
                                    blk.solid_phase.config.transformation_scheme
                                    != "FORWARD"
                                    and x == blk.length_domain.first()
                                )
                                or (
                                    blk.solid_phase.config.transformation_scheme
                                    == "FORWARD"
                                    and x == blk.length_domain.last()
                                )
                            ):
                                calculate_variable_from_constraint(
                                    solid_rxn_gen[t, x, p, j],
                                    solid_phase_stoichiometry_eqn[t, x, p, j],
                                )

            blk.gas_comp_hetero_rxn.activate()
            blk.solid_phase.rate_reaction_stoichiometry_constraint.activate()
            blk.solid_phase_config_rxn_ext.activate()

        if (
            gas_phase.reaction_package is not None
            or solid_phase.reaction_package is not None
        ):
            init_log.info_high(
                "initialize mass balances - with reactions " "and no pressure drop"
            )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if check_optimal_termination(results):
                init_log.info_high(
                    "Initialization Step 3b {}.".format(idaeslog.condition(results))
                )
            else:
                _log.warning("{} Initialization Step 3b Failed.".format(blk.name))

        # Initialize mass balance - with pressure drop
        for t in blk.flowsheet().time:
            for x in blk.length_domain:
                # Unfix all pressure variables except at the inlet
                if blk.gas_phase.properties[t, x].config.defined_state is False:
                    blk.gas_phase.properties[t, x].pressure.unfix()

        blk.gas_phase.pressure_balance.activate()

        if blk.config.has_pressure_change:
            blk.gas_phase_config_pressure_drop.activate()

            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    calculate_variable_from_constraint(
                        blk.gas_phase.deltaP[t, x],
                        blk.gas_phase_config_pressure_drop[t, x],
                    )

            init_log.info_high(
                "initialize mass balances - with reactions " "and pressure drop"
            )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if check_optimal_termination(results):
                init_log.info_high(
                    "Initialization Step 3c {}.".format(idaeslog.condition(results))
                )
            else:
                _log.warning("{} Initialization Step 3c Failed.".format(blk.name))
        # ---------------------------------------------------------------------
        # Initialize energy balance
        if blk.config.energy_balance_type != EnergyBalanceType.none:
            # Initialize dimensionless numbers,
            # mass and heat transfer coefficients
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    calculate_variable_from_constraint(
                        blk.Re_particle[t, x], blk.reynolds_number_particle[t, x]
                    )
                    calculate_variable_from_constraint(
                        blk.Pr_particle[t, x], blk.prandtl_number[t, x]
                    )
                    calculate_variable_from_constraint(
                        blk.Nu_particle[t, x], blk.nusselt_number_particle[t, x]
                    )
                    calculate_variable_from_constraint(
                        blk.gas_solid_htc[t, x], blk.gas_solid_htc_eqn[t, x]
                    )
                    calculate_variable_from_constraint(
                        blk.gas_phase.heat[t, x], blk.gas_phase_heat_transfer[t, x]
                    )
                    calculate_variable_from_constraint(
                        blk.solid_phase.heat[t, x], blk.solid_phase_heat_transfer[t, x]
                    )

            # Unfix temperatures
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    # Unfix all gas temperature variables except at the inlet
                    if blk.gas_phase.properties[t, x].config.defined_state is False:
                        blk.gas_phase.properties[t, x].temperature.unfix()
                for x in blk.length_domain:
                    # Unfix all solid temperature variables except at the inlet
                    if blk.solid_phase.properties[t, x].config.defined_state is False:
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

            init_log.info("Initialize Energy Balances")
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if check_optimal_termination(results):
                init_log.info_high(
                    "Initialization Step 4 {}.".format(idaeslog.condition(results))
                )
            else:
                raise InitializationError(
                    f"{blk.name} failed to initialize successfully. Please "
                    f"check the output logs for more information."
                )

        # Initialize energy balance
        if blk.config.energy_balance_type == EnergyBalanceType.none:
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    # Unfix all gas temperature variables except at the inlet
                    if blk.gas_phase.properties[t, x].config.defined_state is False:
                        blk.gas_phase.properties[t, x].temperature.unfix()
                for x in blk.length_domain:
                    # Unfix all solid temperature variables except at the inlet
                    if blk.solid_phase.properties[t, x].config.defined_state is False:
                        blk.solid_phase.properties[t, x].temperature.unfix()

            blk.isothermal_gas_phase.activate()
            blk.isothermal_solid_phase.activate()

            init_log.info("Initialize Energy Balances")
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee)
            if check_optimal_termination(results):
                init_log.info_high(
                    "Initialization Step 4 {}.".format(idaeslog.condition(results))
                )
            else:
                raise InitializationError(
                    f"{blk.name} failed to initialize successfully. Please "
                    f"check the output logs for more information."
                )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # local aliases used to shorten object names
        gas_phase = self.config.gas_phase_config

        # Get units meta data from property packages
        units_meta_gas = gas_phase.property_package.get_metadata().get_derived_units

        # scale some variables
        if hasattr(self, "bed_diameter"):
            if iscale.get_scaling_factor(self.bed_diameter) is None:
                sf = 1 / value(self.bed_diameter)
                iscale.set_scaling_factor(self.bed_diameter, 1e-1 * sf)

        if hasattr(self, "bed_area"):
            if iscale.get_scaling_factor(self.bed_area) is None:
                sf = 1 / value(constants.pi * (0.5 * self.bed_diameter) ** 2)
                iscale.set_scaling_factor(self.bed_area, 1e-1 * sf)

        if hasattr(self.gas_phase, "area"):
            if iscale.get_scaling_factor(self.gas_phase.area) is None:
                sf = iscale.get_scaling_factor(self.bed_area)
                iscale.set_scaling_factor(self.gas_phase.area, sf)

        if hasattr(self.solid_phase, "area"):
            if iscale.get_scaling_factor(self.solid_phase.area) is None:
                sf = iscale.get_scaling_factor(self.bed_area)
                iscale.set_scaling_factor(self.solid_phase.area, sf)

        if self.config.energy_balance_type != EnergyBalanceType.none:
            if hasattr(self, "Re_particle"):
                for (t, x), v in self.Re_particle.items():
                    if iscale.get_scaling_factor(v) is None:
                        sf1 = iscale.get_scaling_factor(
                            self.gas_phase.properties[t, x].dens_mass
                        )
                        sf2 = 1 / value(
                            pyunits.convert(
                                self.solid_phase.properties[t, x]._params.particle_dia,
                                to_units=units_meta_gas("length"),
                            )
                        )
                        iscale.set_scaling_factor(v, 1e-1 * sf1 * sf2)

            if hasattr(self, "Pr_particle"):
                for (t, x), v in self.Pr_particle.items():
                    if iscale.get_scaling_factor(v) is None:
                        sf1 = iscale.get_scaling_factor(
                            self.gas_phase.properties[t, x].cp_mass
                        )
                        sf2 = iscale.get_scaling_factor(
                            self.gas_phase.properties[t, x].visc_d
                        )
                        iscale.set_scaling_factor(v, sf1 * sf2)

            if hasattr(self, "Nu_particle"):
                for (t, x), v in self.Nu_particle.items():
                    if iscale.get_scaling_factor(v) is None:
                        sf = 1 / (
                            value(
                                (
                                    (
                                        2.0
                                        + 1.1
                                        * (
                                            smooth_abs(self.Re_particle[t, x], self.eps)
                                            ** 0.6
                                        )
                                        ** 3
                                    )
                                    * self.Pr_particle[t, x]
                                )
                                ** (1 / 3)
                            )
                        )
                        iscale.set_scaling_factor(v, 1e-1 * sf)

            if hasattr(self, "gas_solid_htc"):
                if iscale.get_scaling_factor(self.gas_solid_htc) is None:
                    for (t, x), v in self.gas_solid_htc.items():
                        sf1 = iscale.get_scaling_factor(self.Nu_particle[t, x])
                        sf2 = iscale.get_scaling_factor(
                            self.gas_phase.properties[t, x].therm_cond
                        )
                        iscale.set_scaling_factor(v, sf1 * sf2)

            # scale heat variables in CV1D
            for (t, x), v in self.gas_phase.heat.items():
                sf1 = iscale.get_scaling_factor(
                    self.gas_phase.properties[t, x].enth_mol
                )
                sf2 = iscale.get_scaling_factor(
                    self.gas_phase.properties[t, x].flow_mol
                )
                iscale.set_scaling_factor(v, sf1 * sf2)

            for (t, x), v in self.solid_phase.heat.items():
                sf1 = iscale.get_scaling_factor(
                    self.solid_phase.properties[t, x].enth_mass
                )
                sf2 = iscale.get_scaling_factor(
                    self.solid_phase.properties[t, x].flow_mass
                )
                iscale.set_scaling_factor(v, sf1 * sf2)

        # Scale some constraints
        if hasattr(self, "bed_area_eqn"):
            for c in self.bed_area_eqn.values():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.bed_area), overwrite=False
                )

        if hasattr(self, "gas_phase_area"):
            for (t, x), c in self.gas_phase_area.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.bed_area), overwrite=False
                )

        if hasattr(self, "solid_phase_area"):
            for (t, x), c in self.solid_phase_area.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.bed_area), overwrite=False
                )

        if hasattr(self, "gas_super_vel"):
            for (t, x), c in self.gas_super_vel.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.gas_phase.properties[t, x].flow_mol),
                    overwrite=False,
                )

        if hasattr(self, "solid_super_vel"):
            for (t, x), c in self.solid_super_vel.items():
                iscale.constraint_scaling_transform(
                    c,
                    1e-1
                    * iscale.get_scaling_factor(
                        self.solid_phase.properties[t, x].flow_mass
                    ),
                    overwrite=False,
                )

        if hasattr(self, "gas_phase_config_pressure_drop"):
            for (t, x), c in self.gas_phase_config_pressure_drop.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.gas_phase.deltaP[t, x]),
                    overwrite=False,
                )

        if hasattr(self, "gas_phase_config_rxn_ext"):
            for (t, x, r), c in self.gas_phase_config_rxn_ext.items():
                sf1 = iscale.get_scaling_factor(
                    self.gas_phase.reactions[t, x].reaction_rate[r]
                )
                sf2 = iscale.get_scaling_factor(self.bed_area)
                iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

        if hasattr(self, "solid_phase_config_rxn_ext"):
            for (t, x, r), c in self.solid_phase_config_rxn_ext.items():
                sf1 = iscale.get_scaling_factor(
                    self.solid_phase.reactions[t, x].reaction_rate[r]
                )
                sf2 = iscale.get_scaling_factor(self.bed_area)
                iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

        if hasattr(self, "gas_comp_hetero_rxn"):
            for (t, x, p, j), c in self.gas_comp_hetero_rxn.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.gas_phase.mass_transfer_term[t, x, p, j]
                    ),
                    overwrite=False,
                )

        if self.config.energy_balance_type != EnergyBalanceType.none:
            if hasattr(self, "reynolds_number_particle"):
                for (t, x), c in self.reynolds_number_particle.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(self.Re_particle[t, x]),
                        overwrite=False,
                    )

            if hasattr(self, "prandtl_number"):
                for (t, x), c in self.prandtl_number.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(self.Pr_particle[t, x]),
                        overwrite=False,
                    )

            if hasattr(self, "nusselt_number_particle"):
                for (t, x), c in self.nusselt_number_particle.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(self.Nu_particle[t, x]),
                        overwrite=False,
                    )

            if hasattr(self, "gas_solid_htc_eqn"):
                for (t, x), c in self.gas_solid_htc_eqn.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(self.gas_solid_htc[t, x]),
                        overwrite=False,
                    )

            if hasattr(self, "gas_phase_heat_transfer"):
                for (t, x), c in self.gas_phase_heat_transfer.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(self.gas_phase.heat[t, x]),
                        overwrite=False,
                    )

            if hasattr(self, "solid_phase_heat_transfer"):
                for (t, x), c in self.solid_phase_heat_transfer.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(self.solid_phase.heat[t, x]),
                        overwrite=False,
                    )

        elif self.config.energy_balance_type == EnergyBalanceType.none:
            if hasattr(self, "isothermal_gas_phase"):
                for (t, x), c in self.isothermal_gas_phase.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(
                            self.gas_phase.properties[t, x].temperature
                        ),
                        overwrite=False,
                    )

            if hasattr(self, "isothermal_solid_phase"):
                for (t, x), c in self.isothermal_solid_phase.items():
                    iscale.constraint_scaling_transform(
                        c,
                        iscale.get_scaling_factor(
                            self.gas_phase.properties[t, x].temperature
                        ),
                        overwrite=False,
                    )

    def results_plot(blk):
        """
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
        """
        print()
        print(
            "================================= Reactor plots ==============="
            "=================="
        )
        # local aliases used to shorten object names
        gas_phase = blk.config.gas_phase_config
        solid_phase = blk.config.solid_phase_config

        Tg = []
        Ts = []
        P = []
        Ftotal = []
        Mtotal = []
        vg = []

        for t in blk.flowsheet().time:
            for x in blk.gas_phase.length_domain:
                vg.append(value(blk.velocity_superficial_gas[t, x]))

        fig_vg = plt.figure(1)
        plt.plot(blk.gas_phase.length_domain, vg, label="vg")
        plt.legend(loc=0, ncol=2)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Superficial gas velocity [m/s]")
        fig_vg.savefig("superficial_vel.png")

        # Pressure profile
        for t in blk.flowsheet().time:
            for x in blk.gas_phase.length_domain:
                P.append(blk.gas_phase.properties[t, x].pressure.value)

        fig_P = plt.figure(2)
        plt.plot(blk.gas_phase.length_domain, P, label="P")
        plt.legend(loc=0, ncol=2)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Total Pressure [Pa]")
        fig_P.savefig("Pressure.png")

        # Temperature profile
        for t in blk.flowsheet().time:
            for x in blk.gas_phase.length_domain:
                Tg.append(blk.gas_phase.properties[t, x].temperature.value)
            for x in blk.solid_phase.length_domain:
                Ts.append(blk.solid_phase.properties[t, x].temperature.value)
        fig_T = plt.figure(3)
        plt.plot(blk.gas_phase.length_domain, Tg, label="Tg")
        plt.plot(blk.solid_phase.length_domain, Ts, label="Ts")
        plt.legend(loc=0, ncol=2)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Temperature [K]")
        fig_T.savefig("Temperature.png")

        # Bulk gas phase total molar flow rate
        for t in blk.flowsheet().time:
            for x in blk.gas_phase.length_domain:
                Ftotal.append(blk.gas_phase.properties[t, x].flow_mol.value)
        fig_Ftotal = plt.figure(4)
        plt.plot(blk.gas_phase.length_domain, Ftotal)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Total molar gas flow rate [mol/s]")
        fig_Ftotal.savefig("Total_gas_flow.png")

        # Bulk solid phase total mass flow rate
        for t in blk.flowsheet().time:
            for x in blk.solid_phase.length_domain:
                Mtotal.append(blk.solid_phase.properties[t, x].flow_mass.value)
        fig_Mtotal = plt.figure(5)
        plt.plot(blk.solid_phase.length_domain, Mtotal)
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Solid total mass flow rate [kg/s]")
        fig_Mtotal.savefig("Total_solid_flow.png")

        # Gas phase mole fractions
        for t in blk.flowsheet().time:
            for j in gas_phase.property_package.component_list:
                y_frac = []
                for x in blk.gas_phase.length_domain:
                    y_frac.append(
                        value(blk.gas_phase.properties[t, x].mole_frac_comp[j])
                    )
                fig_y = plt.figure(6)
                plt.plot(blk.gas_phase.length_domain, y_frac, label=j)
        plt.legend(loc=0, ncol=len(gas_phase.property_package.component_list))
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("y_frac [-]")
        fig_y.savefig("Gas_mole_fractions.png")

        # Solid phase mass fractions
        for t in blk.flowsheet().time:
            for j in solid_phase.property_package.component_list:
                x_frac = []
                for x in blk.solid_phase.length_domain:
                    x_frac.append(
                        value(blk.solid_phase.properties[t, x].mass_frac_comp[j])
                    )
                fig_x = plt.figure(7)
                plt.plot(blk.solid_phase.length_domain, x_frac, label=j)
        plt.legend(loc=0, ncol=len(solid_phase.property_package.component_list))
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("x_frac [-]")
        fig_x.savefig("Solid_mass_fractions.png")

        # Gas phase concentrations
        for t in blk.flowsheet().time:
            for j in gas_phase.property_package.component_list:
                Cg = []
                for x in blk.gas_phase.length_domain:
                    Cg.append(blk.gas_phase.properties[t, x].dens_mol_comp[j].value)
                fig_Cg = plt.figure(8)
                plt.plot(blk.gas_phase.length_domain, Cg, label=j)
        plt.legend(loc=0, ncol=len(gas_phase.property_package.component_list))
        plt.grid()
        plt.xlabel("Bed height [-]")
        plt.ylabel("Concentration [mol/m3]")
        fig_Cg.savefig("Gas_concentration.png")

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
