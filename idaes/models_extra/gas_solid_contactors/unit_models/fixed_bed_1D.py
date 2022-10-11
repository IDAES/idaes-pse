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
IDAES 1D Fixed Bed model.

"""

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
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.util.subsystems import TemporarySubsystemManager
from pyomo.common.collections import ComponentSet
from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
    solve_strongly_connected_components,
)

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
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import get_index_set_except
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants as constants
from idaes.core.util.math import smooth_abs
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

__author__ = "Chinedu Okoli"

# Set up logger
_log = idaeslog.getLogger(__name__)

# Assumptions:
# Only gas and solid phases, which are explicitly named as such
# Perfect mixing in gas phase
# Axially varying gas and solid phase
# Reaction/adsorption on the solid is the rate limiting step


@declare_process_block_class("FixedBed1D")
class FixedBed1DData(UnitModelBlockData):
    """
    1D Fixed Bed Reactor Model Class
    """

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
            doc="""Method to use to transform domain. Must be a method recognized
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
            default="forward_flow",
            domain=In(["forward_flow", "reverse_flow"]),
            description="Flow configuration of Fixed Bed",
            doc="""Flow configuration of Fixed Bed
**default** - "forward_flow".
**Valid values:** {
**"forward_flow"** - gas flows from 0 to 1,
**"reverse_flow"** -  gas flows from 1 to 0.}""",
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
            default="ergun_correlation",
            domain=In(["ergun_correlation", "simple_correlation"]),
            description="Construction flag for type of pressure drop",
            doc="""Indicates what type of pressure drop correlation should be used,
**default** - "ergun_correlation".
**Valid values:** {
**"ergun_correlation"** - Use the Ergun equation,
**"simple_correlation"** - Use a simplified pressure drop correlation.}""",
        ),
    )

    # Create template for phase specific config arguments
    _PhaseTemplate = UnitModelBlockData.CONFIG()
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

        # local aliases used to shorten object names
        solid_phase = self.config.solid_phase_config
        gas_phase = self.config.gas_phase_config

        # Get units meta data from property packages
        units_meta_solid = solid_phase.property_package.get_metadata().get_derived_units

        # Set flow direction for the gas control volume
        # Gas flows from 0 to 1
        if self.config.flow_type == "forward_flow":
            set_direction_gas = FlowDirection.forward
        # Gas flows from 1 to 0
        if self.config.flow_type == "reverse_flow":
            set_direction_gas = FlowDirection.backward

        # Consistency check for flow direction, transformation method and
        # transformation scheme
        if (
            self.config.flow_type == "forward_flow"
            and self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme is None
        ):
            self.config.transformation_scheme = "BACKWARD"
        elif (
            self.config.flow_type == "reverse_flow"
            and self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme is None
        ):
            self.config.transformation_scheme = "FORWARD"
        elif (
            self.config.flow_type == "forward_flow"
            and self.config.transformation_method == "dae.collocation"
            and self.config.transformation_scheme is None
        ):
            self.config.transformation_scheme = "LAGRANGE-RADAU"
        elif (
            self.config.flow_type == "reverse_flow"
            and self.config.transformation_method == "dae.collocation"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_method argument."
                "Must be "
                "dae.finite_difference "
                "if "
                "flow_type is"
                " "
                "reverse_flow"
                ".".format(self.name)
            )
        elif (
            self.config.flow_type == "forward_flow"
            and self.config.transformation_scheme == "FORWARD"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_scheme argument. "
                "Must be "
                "BACKWARD "
                "if flow_type is"
                " "
                "forward_flow"
                ".".format(self.name)
            )
        elif (
            self.config.flow_type == "reverse_flow"
            and self.config.transformation_scheme == "BACKWARD"
        ):
            raise ConfigurationError(
                "{} invalid value for "
                "transformation_scheme argument."
                "Must be "
                "FORWARD "
                "if "
                "flow_type is"
                " "
                "reverse_flow"
                ".".format(self.name)
            )
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

        # Set arguments for gas sides if homoogeneous reaction block
        if gas_phase.reaction_package is not None:
            has_rate_reaction_gas_phase = True
        else:
            has_rate_reaction_gas_phase = False

        # Set arguments for gas and solid sides if heterogeneous reaction block
        if solid_phase.reaction_package is not None:
            has_mass_transfer_gas_phase = True
        else:
            has_mass_transfer_gas_phase = False

        # Set heat transfer terms
        if self.config.energy_balance_type != EnergyBalanceType.none:
            has_heat_transfer = True
        else:
            has_heat_transfer = False

        # Set heat of reaction terms
        if (
            self.config.energy_balance_type != EnergyBalanceType.none
            and gas_phase.reaction_package is not None
        ):
            has_heat_of_reaction_gas_phase = True
        else:
            has_heat_of_reaction_gas_phase = False

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
            dynamic=True,  # Fixed beds must be dynamic
            has_holdup=True,  # holdup must be True for fixed beds
            area_definition=DistributedVars.variant,
            property_package=gas_phase.property_package,
            property_package_args=gas_phase.property_package_args,
            reaction_package=gas_phase.reaction_package,
            reaction_package_args=gas_phase.reaction_package_args,
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

        if gas_phase.reaction_package is not None:
            self.gas_phase.add_reaction_blocks(
                has_equilibrium=gas_phase.has_equilibrium_reactions
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

        # Build Solid Phase StateBlock
        # As there is no solid flow, there is not need for a control volume as a
        # set of indexed state blocks will be sufficient.
        # Defined state is set to True so that the "sum(mass_frac)=1" eqn in
        # the solid state block is deactivated. This is done here as there is
        # currently no way to deactivate the constraint at the initial time
        # for batch systems (i.e. no inlet or outlet ports).
        # The "sum(mass_frac)=1 for all t neq 0" eqn is written in the
        # unit model instead
        self.solid_properties = solid_phase.property_package.state_block_class(
            self.flowsheet().time,
            self.length_domain,
            parameters=solid_phase.property_package,
            defined_state=True,
        )

        if solid_phase.reaction_package is not None:
            # TODO - a generalization of the heterogeneous reaction block
            # The heterogeneous reaction block does not use the
            # add_reaction_blocks in control volumes as control volumes are
            # currently setup to handle only homogeneous reaction properties.
            # Thus appending the heterogeneous reaction block to the
            # solid state block is currently hard coded here.

            tmp_dict = dict(**solid_phase.reaction_package_args)
            tmp_dict["gas_state_block"] = self.gas_phase.properties
            tmp_dict["solid_state_block"] = self.solid_properties
            tmp_dict["parameters"] = solid_phase.reaction_package
            self.solid_reactions = solid_phase.reaction_package.reaction_block_class(
                self.flowsheet().time,
                self.length_domain,
                doc="Reaction properties in control volume",
                **tmp_dict,
            )

        # =========================================================================
        """ Add ports"""
        # Add Ports for gas side
        self.add_inlet_port(name="gas_inlet", block=self.gas_phase)
        self.add_outlet_port(name="gas_outlet", block=self.gas_phase)

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

        # Solid phase variables
        self.solid_phase_area = Var(
            self.flowsheet().time,
            self.length_domain,
            domain=Reals,
            initialize=1,
            doc="Solid phase area",
            units=units_meta_solid("area"),
        )

        self.solid_material_holdup = Var(
            self.flowsheet().time,
            self.length_domain,
            solid_phase.property_package.component_list,
            domain=Reals,
            initialize=1.0,
            doc="Solid material holdup per unit length",
            units=units_meta_solid("mass") / units_meta_solid("length"),
        )

        self.solid_material_accumulation = DerivativeVar(
            self.solid_material_holdup,
            wrt=self.flowsheet().config.time,
            doc="Solid material accumulation",
            units=units_meta_solid("mass")
            / units_meta_solid("length")
            / units_meta_solid("time"),
        )

        if self.config.energy_balance_type != EnergyBalanceType.none:
            # Solid phase - gas to solid heat transfer
            self.solid_phase_heat = Var(
                self.flowsheet().time,
                self.length_domain,
                domain=Reals,
                initialize=0.0,
                doc="Heat transferred per unit length",
                units=units_meta_solid("power") / units_meta_solid("length"),
            )

            self.solid_energy_holdup = Var(
                self.flowsheet().time,
                self.length_domain,
                initialize=1,
                doc="Solid phase energy holdup",
                units=units_meta_solid("energy") / units_meta_solid("length"),
            )

            self.solid_energy_accumulation = DerivativeVar(
                self.solid_energy_holdup,
                initialize=0,
                wrt=self.flowsheet().config.time,
                doc="Solid energy accumulation",
                units=units_meta_solid("power") / units_meta_solid("length"),
            )

        # =========================================================================
        # Add performance equations

        # ---------------------------------------------------------------------
        # TODO - Deactivate/delete sum_component_eqns of initial conditions for
        #  all space in gas_phase as currently no way to tell CV1D to not build
        #  sum_eqns for initial conditions
        if self.gas_phase.config.dynamic:
            t0 = self.flowsheet().time.first()
            for x in self.length_domain:
                blk = self.gas_phase.properties[t0, x]
                if hasattr(blk, "sum_component_eqn"):
                    blk.sum_component_eqn.deactivate()

        # ---------------------------------------------------------------------
        # Geometry contraints

        # Bed area
        @self.Constraint(doc="Bed area")
        def bed_area_eqn(b):
            return b.bed_area == (constants.pi * (0.5 * b.bed_diameter) ** 2)

        # Area of gas side, and solid side
        @self.Constraint(self.flowsheet().time, self.length_domain, doc="Gas side area")
        def gas_phase_area_constraint(b, t, x):
            return (
                b.gas_phase.area[t, x]
                == pyunits.convert(b.bed_area, to_units=units_meta_gas("area"))
                * b.solid_properties[t, x]._params.voidage
            )

        @self.Constraint(
            self.flowsheet().time, self.length_domain, doc="Solid side area"
        )
        def solid_phase_area_constraint(b, t, x):
            return b.solid_phase_area[t, x] == b.bed_area * (
                1 - b.solid_properties[t, x]._params.voidage
            )

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
                            b.solid_properties[t, x].dens_mass_particle,
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
                    * (1 - b.solid_properties[t, x]._params.voidage) ** 2
                    * b.gas_phase.properties[t, x].visc_d
                    * b.velocity_superficial_gas[t, x]
                    / (
                        pyunits.convert(
                            b.solid_properties[t, x]._params.particle_dia,
                            to_units=units_meta_solid("length"),
                        )
                        ** 2
                        * b.solid_properties[t, x]._params.voidage ** 3
                    )
                ) + (
                    (1.75 * pyunits.dimensionless)
                    * b.gas_phase.properties[t, x].dens_mass
                    * (1 - b.solid_properties[t, x]._params.voidage)
                    * b.velocity_superficial_gas[t, x] ** 2
                    / (
                        pyunits.convert(
                            b.solid_properties[t, x]._params.particle_dia,
                            to_units=units_meta_solid("length"),
                        )
                        * b.solid_properties[t, x]._params.voidage ** 3
                    )
                )

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
                        b.solid_reactions[t, x].rate_reaction_stoichiometry[r, p, j]
                        * b.solid_reactions[t, x].reaction_rate[r]
                        for r in (solid_phase.reaction_package.rate_reaction_idx)
                    )
                    * b.solid_phase_area[t, x]
                )

        # ---------------------------------------------------------------------
        # Solid phase component balance

        # material holdup constraint
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            solid_phase.property_package.component_list,
            doc="Solid phase material holdup constraints",
        )
        def solid_material_holdup_calculation(b, t, x, j):
            return b.solid_material_holdup[t, x, j] == (
                b.solid_phase_area[t, x]
                * b.solid_properties[t, x].get_material_density_terms("Sol", j)
            )

        # Add component balances
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            solid_phase.property_package.component_list,
            doc="Material balances",
        )
        def solid_material_balances(b, t, x, j):
            if solid_phase.reaction_package is not None:
                return b.solid_material_accumulation[t, x, j] == (
                    b.solid_phase_area[t, x]
                    * b.solid_properties[t, x]._params.mw_comp[j]
                    * sum(
                        b.solid_reactions[t, x].reaction_rate[r]
                        * solid_phase.reaction_package.rate_reaction_stoichiometry[
                            r, "Sol", j
                        ]
                        for r in solid_phase.reaction_package.rate_reaction_idx
                    )
                )
            else:
                return b.solid_material_accumulation[t, x, j] == 0

        # Sum of mass fractions at all time equals 1
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Sum of mass fractions at all time",
        )
        def solid_sum_component_eqn(b, t, x):
            if t == b.flowsheet().config.time.first():
                return Constraint.Skip
            else:
                return 1 == sum(
                    b.solid_properties[t, x].mass_frac_comp[j]
                    for j in b.solid_properties[t, x]._params.component_list
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
                    b.solid_phase_heat[t, x],
                    to_units=units_meta_gas("power") / units_meta_gas("length"),
                ) * pyunits.convert(
                    b.solid_properties[t, x]._params.particle_dia,
                    to_units=units_meta_gas("length"),
                ) == 6 * b.gas_solid_htc[
                    t, x
                ] * (
                    b.gas_phase.properties[t, x].temperature
                    - pyunits.convert(
                        b.solid_properties[t, x].temperature,
                        to_units=units_meta_gas("temperature"),
                    )
                ) * pyunits.convert(
                    b.solid_phase_area[t, x], to_units=units_meta_gas("area")
                )

            # Solid phase energy balance
            # Accumulation equal to heat transfer
            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Solid phase energy holdup constraints",
            )
            def solid_energy_holdup_calculation(b, t, x):
                return b.solid_energy_holdup[t, x] == (
                    b.solid_phase_area[t, x]
                    * pyunits.convert(
                        b.solid_properties[t, x].get_energy_density_terms("Sol"),
                        to_units=units_meta_solid("energy")
                        / units_meta_solid("volume"),
                    )
                )

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Solid phase energy balances",
            )
            def solid_enthalpy_balances(b, t, x):
                if solid_phase.reaction_package is not None:
                    return b.solid_energy_accumulation[t, x] / units_meta_solid(
                        "time"
                    ) == b.solid_phase_heat[t, x] + -sum(
                        b.solid_reactions[t, x].reaction_rate[r]
                        * b.solid_phase_area[t, x]
                        * pyunits.convert(
                            b.solid_reactions[t, x].dh_rxn[r],
                            to_units=units_meta_solid("energy_mole"),
                        )
                        for r in solid_phase.reaction_package.rate_reaction_idx
                    )
                else:
                    return b.solid_energy_accumulation[t, x] == b.solid_phase_heat[t, x]

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
                        b.solid_properties[t, x]._params.particle_dia,
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
                        b.solid_properties[t, x]._params.particle_dia,
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
                    b.solid_properties[t, x]._params.particle_dia,
                    to_units=units_meta_gas("length"),
                ) == -6 * b.gas_solid_htc[t, x] * (
                    b.gas_phase.properties[t, x].temperature
                    - pyunits.convert(
                        b.solid_properties[t, x].temperature,
                        to_units=units_meta_gas("temperature"),
                    )
                ) * pyunits.convert(
                    b.solid_phase_area[t, x], to_units=units_meta_gas("area")
                )

        elif self.config.energy_balance_type == EnergyBalanceType.none:
            # if self.config.energy_balance_type == EnergyBalanceType.none:
            # If energy balance is none fix gas and solid temperatures to
            # initial solid temperature
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Isothermal gas phase constraint",
            )
            def isothermal_gas_phase(b, t, x):
                # Skip constraint at the initial and boundary points
                t0 = b.flowsheet().config.time.first()
                if b.config.flow_type == "forward_flow":
                    x_inlet = b.length_domain.first()
                else:
                    x_inlet = b.length_domain.last()
                if x == x_inlet or t == b.flowsheet().config.time.first():
                    return Constraint.Skip
                else:
                    return (
                        b.gas_phase.properties[t, x].temperature
                        == b.solid_properties[t0, x].temperature
                    )

            @self.Constraint(
                self.flowsheet().config.time,
                self.length_domain,
                doc="Isothermal solid phase constraint",
            )
            def isothermal_solid_phase(b, t, x):
                # Skip constraint at the initial point
                t0 = b.flowsheet().config.time.first()
                if t == b.flowsheet().config.time.first():
                    return Constraint.Skip
                else:
                    return (
                        b.solid_properties[t, x].temperature
                        == b.solid_properties[t0, x].temperature
                    )

        # Set default scaling values for some CV1D vars (to avoid warnings)
        if hasattr(self.gas_phase, "heat"):
            for (t, x), v in self.gas_phase.heat.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-6)
        if hasattr(self.gas_phase, "area"):
            for (t, x), v in self.gas_phase.area.items():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

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
        Initialization routine for 1DFixedBed unit.

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
        # and property block constraints active. Additionally, in control
        # volumes - keep conservation linking constraints and
        # holdup calculation (for dynamic flowsheets) constraints active

        geometry_constraints_terms = [
            "bed_area_eqn",
            "solid_phase_area_constraint",
            "gas_phase_area_constraint",
        ]
        endswith_terms = (
            "_disc_eq",
            "linking_constraint",
            "linking_constraints",
            "_holdup_calculation",
        )
        startswith_terms = ("properties", "solid_properties")

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
        solid_phase_flags = blk.solid_properties.initialize(
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
        blk.solid_properties.release_state(solid_phase_flags)
        for t in blk.flowsheet().time:
            for x in blk.length_domain:
                blk.gas_phase.properties[t, x].pressure.fix()
                blk.gas_phase.properties[t, x].temperature.fix()
                blk.solid_properties[t, x].temperature.fix()

        if blk.config.has_holdup is True:
            # Fix initial conditions of flowrate
            t0 = blk.flowsheet().time.first()
            if blk.config.flow_type == "forward_flow":
                x_inlet = blk.length_domain.first()
            else:
                x_inlet = blk.length_domain.last()
            for x in blk.length_domain:
                if x != x_inlet:
                    blk.gas_phase.properties[t0, x].flow_mol.fix()
                    blk.gas_phase.properties[t0, x].sum_component_eqn.deactivate()

        blk.gas_phase.material_balances.activate()

        if gas_phase.reaction_package is not None:
            for t in blk.flowsheet().time:
                gas_rxn_gen = blk.gas_phase.rate_reaction_generation
                for x in blk.length_domain:
                    for p in gas_phase.property_package.phase_list:
                        for j in gas_phase.property_package.component_list:
                            (gas_rxn_gen[t, x, p, j].fix(0.0))

        blk.solid_material_balances.activate()
        blk.solid_sum_component_eqn.activate()

        if solid_phase.reaction_package is not None:
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    for r in solid_phase.reaction_package.rate_reaction_idx:
                        blk.solid_reactions[t, x].reaction_rate[r].fix(0.0)

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
            results = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
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
            gas_mass_transfer_term = blk.gas_phase.mass_transfer_term

            # Initialize reaction property package
            blk.solid_reactions.activate()
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    blk.solid_reactions[t, x].reaction_rate.unfix()
                    obj = blk.solid_reactions[t, x]
                    for c in obj.component_objects(Constraint, descend_into=False):
                        c.activate()

            blk.solid_reactions.initialize(outlvl=outlvl, optarg=optarg, solver=solver)

            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    for p in gas_phase.property_package.phase_list:
                        for j in gas_phase.property_package.component_list:
                            (gas_mass_transfer_term[t, x, p, j].unfix())
                            calculate_variable_from_constraint(
                                gas_mass_transfer_term[t, x, p, j],
                                blk.gas_comp_hetero_rxn[t, x, p, j],
                            )

            blk.gas_comp_hetero_rxn.activate()

        if (
            gas_phase.reaction_package is not None
            or solid_phase.reaction_package is not None
        ):
            init_log.info_high(
                "initialize mass balances - with reactions " "and no pressure drop"
            )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
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
                results = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
            if check_optimal_termination(results):
                init_log.info_high(
                    "Initialization Step 3c {}.".format(idaeslog.condition(results))
                )
            else:
                _log.warning("{} Initialization Step 3c Failed.".format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize energy balance
        calc_var_kwds = {"eps": 5e-6}
        if blk.config.energy_balance_type != EnergyBalanceType.none:
            # Unfix temperatures
            if blk.config.flow_type == "forward_flow":
                x_inlet = blk.length_domain.first()
            else:
                x_inlet = blk.length_domain.last()
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    if t != blk.flowsheet().config.time.first() and x != x_inlet:
                        # Unfix gas temperature variables except at the inlet
                        blk.gas_phase.properties[t, x].temperature.unfix()

            for t in blk.flowsheet().time:
                if t != blk.flowsheet().config.time.first():
                    for x in blk.length_domain:
                        # Unfix solid temperature variables except at the inlet
                        blk.solid_properties[t, x].temperature.unfix()

            # Initialize dimensionless numbers,
            # mass and heat transfer coefficients
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    # Initialize gas temperature to solid temperature values
                    calculate_variable_from_constraint(
                        blk.Re_particle[t, x],
                        blk.reynolds_number_particle[t, x],
                        **calc_var_kwds,
                    )
                    calculate_variable_from_constraint(
                        blk.Pr_particle[t, x], blk.prandtl_number[t, x], **calc_var_kwds
                    )
                    calculate_variable_from_constraint(
                        blk.Nu_particle[t, x],
                        blk.nusselt_number_particle[t, x],
                        **calc_var_kwds,
                    )
                    calculate_variable_from_constraint(
                        blk.gas_solid_htc[t, x],
                        blk.gas_solid_htc_eqn[t, x],
                        **calc_var_kwds,
                    )
                    calculate_variable_from_constraint(
                        blk.gas_phase.heat[t, x],
                        blk.gas_phase_heat_transfer[t, x],
                        **calc_var_kwds,
                    )
                    calculate_variable_from_constraint(
                        blk.solid_phase_heat[t, x],
                        blk.solid_phase_heat_transfer[t, x],
                        **calc_var_kwds,
                    )

                    if t != blk.flowsheet().time.first():
                        calculate_variable_from_constraint(
                            blk.gas_phase.energy_accumulation[t, x, "Vap"],
                            blk.gas_phase.energy_accumulation_disc_eq[t, x, "Vap"],
                            **calc_var_kwds,
                        )
                        calculate_variable_from_constraint(
                            blk.solid_energy_accumulation[t, x],
                            blk.solid_energy_accumulation_disc_eq[t, x],
                            **calc_var_kwds,
                        )

            blk.reynolds_number_particle.activate()
            blk.prandtl_number.activate()
            blk.nusselt_number_particle.activate()
            blk.gas_solid_htc_eqn.activate()

            # Activate gas phase energy balance equations
            blk.gas_phase_heat_transfer.activate()
            blk.gas_phase.enthalpy_balances.activate()

            # Activate solid phase energy balance equations
            blk.solid_phase_heat_transfer.activate()
            blk.solid_enthalpy_balances.activate()

            init_log.info("Initialize Energy Balances")
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = opt.solve(blk, tee=slc.tee, symbolic_solver_labels=True)
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
            if blk.config.flow_type == "forward_flow":
                x_inlet = blk.length_domain.first()
            else:
                x_inlet = blk.length_domain.last()
            t0 = blk.flowsheet().config.time.first()
            for t in blk.flowsheet().time:
                for x in blk.length_domain:
                    if t != blk.flowsheet().config.time.first() and x != x_inlet:
                        # Unfix gas temperature variables except at the inlet
                        blk.gas_phase.properties[t, x].temperature.unfix()
                        blk.gas_phase.properties[t, x].temperature = value(
                            blk.solid_properties[t0, x].temperature
                        )

            for t in blk.flowsheet().time:
                t0 = blk.flowsheet().config.time.first()
                if t != blk.flowsheet().config.time.first():
                    for x in blk.length_domain:
                        # Unfix solid temperature variables except initial
                        blk.solid_properties[t, x].temperature = value(
                            blk.solid_properties[t0, x].temperature
                        )
                        blk.solid_properties[t, x].temperature.unfix()

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

    def block_triangularization_initialize(
        blk,
        gas_phase_state_args=None,
        solid_phase_state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        calc_var_kwds=None,
    ):
        """
        Block triangularization (BT) initialization routine for 1D FixedBed unit.

        The BT initialization routine initializes the unit model in the following steps:

        Step 1: Convert the system of equations from a partial differential set of
                equations (PDAE) to an algebraic set of equations (AE) by deactivating
                the discretization equations and sum_component_equations (if present),
                and fixing the state variables to a guess (typically inlet conditions).
        Step 2: Decompose the AE into strongly connected components
                via block triangularization and solve individual blocks
        Step 3: revert the fixed variables and deactivated constraints to their
                original states.

        The initialization procedure is wrapped in the TemporarySubsystemManager utility
        to ensure that the deactivation of constraints and fixing of variables are only
        temporary and the final model structure is unchanged after initialization.

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
            calc_var_kwds: Dictionary of keyword arguments for
                calculate_variable_from_constraint

        Returns:
            None
        """

        # Set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")

        # =========================================================================
        # Check that the original model has zero degrees of freedom
        if degrees_of_freedom(blk) != 0:
            msg = (
                "Original model has nonzero degrees of freedom. This was "
                "unexpected. "
                + "\n"
                + "Degrees of freedom: "
                + str(degrees_of_freedom(blk))
            )
            init_log.error(msg)
            raise ValueError("Nonzero degrees of freedom.")

        # =========================================================================
        # Create list of constraints to deactivate
        to_deactivate = []
        # Deactivate gas/solid phase sum equations and discretization equations
        # (mass/energy flow, mass/energy accumulation, pressure)
        endswith_terms = ("disc_eq", "sum_component_eqn")
        for c in blk.component_objects(Constraint, descend_into=True):
            if c.local_name.endswith(endswith_terms):
                c.attribute_errors_generate_exceptions = False
                to_deactivate.append(c)

        # Additional equations to deactivate if energy balance type is none
        if hasattr(blk, "isothermal_gas_phase"):
            isothermal_gas_phase_eqn = blk.isothermal_gas_phase[...]
            isothermal_gas_phase_eqn.attribute_errors_generate_exceptions = False
            to_deactivate.extend(isothermal_gas_phase_eqn)
        if hasattr(blk, "isothermal_solid_phase"):
            isothermal_solid_phase_eqn = blk.isothermal_solid_phase[...]
            isothermal_solid_phase_eqn.attribute_errors_generate_exceptions = False
            to_deactivate.extend(isothermal_solid_phase_eqn)

        # =========================================================================
        # Variables to fix during initialization

        # List of state vars to fix (typically fixed to boundary/initial conditions)
        gas_flags = fix_state_vars(blk.gas_phase.properties, gas_phase_state_args)
        solid_flags = fix_state_vars(blk.solid_properties, solid_phase_state_args)

        # List of additional variables to fix
        to_fix = []
        # Flow derivative variables are fixed to ensure square AE
        set_value = 0
        time_set = blk.flowsheet().time
        for t in blk.flowsheet().time:
            for var in ComponentSet(blk.component_objects(Var)):
                if var.name.endswith("_flow_dx"):
                    n = var.index_set().dimen
                    if n == 1:
                        var.set_value(set_value)
                        to_fix.append(var[t])
                    elif n >= 2:
                        info = get_index_set_except(var, time_set)
                        non_time_set = info["set_except"]
                        index_getter = info["index_getter"]
                        for non_time_index in non_time_set:
                            index = index_getter(non_time_index, t)
                            var[index].set_value(set_value)
                            to_fix.append(var[index])

        # =========================================================================
        with TemporarySubsystemManager(to_fix=to_fix, to_deactivate=to_deactivate):

            # Check if the system is structurally singular
            igraph = IncidenceGraphInterface(blk)
            N = len(igraph.constraints)
            M = len(igraph.variables)
            matching = igraph.maximum_matching()
            if not (N == len(matching) and M == len(matching)):
                msg = (
                    "Model is structurally singular as maximum matching "
                    "does not include all costraints and variables "
                    "\n"
                    "Number of constraints: "
                    + str(N)
                    + "\n"
                    + "Number of variables: "
                    + str(M)
                    + "\n"
                    + "Maximum matching: "
                    + str(len(matching))
                )
                init_log.warning(msg)
                raise ValueError("Structural singularity")

            # Decompose AE system into strongly connected components and solve
            if calc_var_kwds is None:
                calc_var_kwds = {"eps": 1e-8}
            try:
                solve_strongly_connected_components(
                    blk, calc_var_kwds=calc_var_kwds, solver=solver
                )
            except RuntimeError:
                msg = (
                    "Initial value for variable results in a derivative value "
                    "that is very close to zero. Please provide a different initial "
                    "guess and/or adjust calc_var_kwds tolerance setting."
                )
                init_log.warning(msg)
                raise

        # Revert the state vars to their original state
        revert_state_vars(blk.gas_phase.properties, gas_flags)
        revert_state_vars(blk.solid_properties, solid_flags)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # local aliases used to shorten object names
        gas_phase = self.config.gas_phase_config
        solid_phase = self.config.solid_phase_config

        # Get units meta data from property packages
        units_meta_gas = gas_phase.property_package.get_metadata().get_derived_units

        # scale some variables
        if hasattr(self, "bed_diameter"):
            sf = 1 / value(self.bed_diameter)
            iscale.set_scaling_factor(self.bed_diameter, 1e-1 * sf)

        if hasattr(self, "bed_area"):
            sf = 1 / value(constants.pi * (0.5 * self.bed_diameter) ** 2)
            iscale.set_scaling_factor(self.bed_area, 1e-1 * sf)

        if hasattr(self, "solid_phase_area"):
            sf = iscale.get_scaling_factor(self.bed_area)
            iscale.set_scaling_factor(self.solid_phase_area, sf)

        if hasattr(self.gas_phase, "area"):
            sf = iscale.get_scaling_factor(self.bed_area)
            iscale.set_scaling_factor(self.gas_phase.area, sf)

        if self.config.energy_balance_type != EnergyBalanceType.none:
            if hasattr(self, "Re_particle"):
                for (t, x), v in self.Re_particle.items():
                    if iscale.get_scaling_factor(v) is None:
                        sf1 = iscale.get_scaling_factor(
                            self.gas_phase.properties[t, x].dens_mass
                        )
                        sf2 = 1 / value(
                            pyunits.convert(
                                self.solid_properties[t, x]._params.particle_dia,
                                to_units=units_meta_gas("length"),
                            )
                        )
                        iscale.set_scaling_factor(v, 1e-1 * sf1 * sf2)

        if hasattr(self, "solid_material_accumulation"):
            for (t, x, j), v in self.solid_material_accumulation.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.solid_properties[t, x].get_material_density_terms(
                            "Sol", j
                        ),
                        default=1,
                        warning=True,
                    )
                    iscale.set_scaling_factor(v, sf)

        if self.config.energy_balance_type != EnergyBalanceType.none:
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
                        iscale.set_scaling_factor(v, 1e-2 * sf)

            if hasattr(self, "gas_solid_htc"):
                for (t, x), v in self.gas_solid_htc.items():
                    if iscale.get_scaling_factor(v) is None:
                        sf1 = iscale.get_scaling_factor(self.Nu_particle[t, x])
                        sf2 = iscale.get_scaling_factor(
                            self.gas_phase.properties[t, x].therm_cond
                        )
                        iscale.set_scaling_factor(v, sf1 * sf2)

            if hasattr(self, "solid_energy_accumulation"):
                for (t, x), v in self.solid_energy_accumulation.items():
                    if iscale.get_scaling_factor(v) is None:
                        sf = iscale.get_scaling_factor(
                            self.solid_properties[t, x].get_energy_density_terms("Sol"),
                            default=1,
                            warning=True,
                        )
                        iscale.set_scaling_factor(v, sf)

            if hasattr(self, "solid_energy_holdup"):
                for (t, x), v in self.solid_energy_holdup.items():
                    sf1 = iscale.get_scaling_factor(self.solid_phase_area[t, x])
                    sf2 = iscale.get_scaling_factor(
                        self.solid_properties[t, x].get_energy_density_terms("Sol")
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

            for (t, x), v in self.solid_phase_heat.items():
                sf1 = iscale.get_scaling_factor(self.solid_properties[t, x].enth_mass)
                sf2 = iscale.get_scaling_factor(
                    self.solid_properties[t, x].dens_mass_particle
                )
                sf3 = iscale.get_scaling_factor(self.solid_phase_area[t, x])
                iscale.set_scaling_factor(v, sf1 * sf2 * sf3)

        # Scale some constraints
        if hasattr(self, "bed_area_eqn"):
            for c in self.bed_area_eqn.values():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.bed_area), overwrite=False
                )

        if hasattr(self, "gas_phase_area_constraint"):
            for (t, x), c in self.gas_phase_area_constraint.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.bed_area), overwrite=False
                )

        if hasattr(self, "solid_phase_area_constraint"):
            for (t, x), c in self.solid_phase_area_constraint.items():
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

        if hasattr(self, "gas_comp_hetero_rxn"):
            for (t, x, p, j), c in self.gas_comp_hetero_rxn.items():
                iscale.constraint_scaling_transform(
                    c,
                    1e-3
                    * iscale.get_scaling_factor(
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
                        1e-2 * iscale.get_scaling_factor(self.gas_phase.heat[t, x]),
                        overwrite=False,
                    )

            if hasattr(self, "solid_phase_heat_transfer"):
                for (t, x), c in self.solid_phase_heat_transfer.items():
                    iscale.constraint_scaling_transform(
                        c,
                        1e-2 * iscale.get_scaling_factor(self.solid_phase_heat[t, x]),
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
                            self.solid_properties[t, x].temperature
                        ),
                        overwrite=False,
                    )

        # Scale some constraints
        if hasattr(self, "solid_material_holdup_calculation"):
            for (t, x, j), c in self.solid_material_holdup_calculation.items():
                sf1 = iscale.get_scaling_factor(self.solid_phase_area[t, x])
                sf2 = iscale.get_scaling_factor(
                    self.solid_properties[t, x].dens_mass_particle
                )
                iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

        if hasattr(self, "solid_material_accumulation_disc_eq"):
            for (t, x, j), c in self.solid_material_accumulation_disc_eq.items():
                iscale.constraint_scaling_transform(
                    c,
                    1e-6
                    * iscale.get_scaling_factor(
                        self.solid_material_accumulation[t, x, j]
                    ),
                    overwrite=False,
                )

        if hasattr(self, "solid_material_balances"):
            component_list = solid_phase.property_package.component_list
            # Get a single representative value in component list
            for j in component_list:
                break
            if solid_phase.reaction_package is not None:
                rate_reaction_idx = solid_phase.reaction_package.rate_reaction_idx
                # Get a single representative value in rate_reaction index list
                for r in rate_reaction_idx:
                    break
                for (t, x, j), c in self.solid_material_balances.items():
                    sf1 = iscale.get_scaling_factor(self.solid_phase_area[t, x])
                    sf2 = iscale.get_scaling_factor(
                        self.solid_properties[t, x]._params.mw_comp[j],
                        default=1e1,
                    )
                    sf3 = iscale.get_scaling_factor(
                        self.solid_reactions[t, x].reaction_rate[r]
                    )
                    iscale.constraint_scaling_transform(
                        c, sf1 * sf2 * sf3, overwrite=False
                    )
            else:
                for (t, x, j), c in self.solid_material_balances.items():
                    sf1 = iscale.get_scaling_factor(self.solid_phase_area[t, x])
                    sf2 = iscale.get_scaling_factor(
                        self.solid_properties[t, x]._params.mw_comp[j],
                        default=1e1,
                    )
                    iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

        if hasattr(self, "solid_sum_component_eqn"):
            component_list = solid_phase.property_package.component_list
            # Get a single representative value in component list
            for j in component_list:
                break
            for (t, x), c in self.solid_sum_component_eqn.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.solid_properties[t, x].mass_frac_comp[j]
                    ),
                    overwrite=False,
                )

        if self.config.energy_balance_type != EnergyBalanceType.none:

            if hasattr(self, "solid_energy_holdup_calculation"):
                for (t, x), c in self.solid_energy_holdup_calculation.items():
                    sf1 = iscale.get_scaling_factor(self.solid_phase_area[t, x])
                    sf2 = iscale.get_scaling_factor(
                        self.solid_properties[t, x].get_energy_density_terms("Sol")
                    )
                    iscale.constraint_scaling_transform(c, sf1 * sf2, overwrite=False)

            if hasattr(self, "solid_enthalpy_balances"):
                for (t, x), c in self.solid_enthalpy_balances.items():
                    sf1 = iscale.get_scaling_factor(self.solid_phase_area[t, x])
                    sf2 = iscale.get_scaling_factor(
                        self.solid_properties[t, x].enth_mass
                    )
                    iscale.constraint_scaling_transform(
                        c, 1e-2 * sf1 * sf2, overwrite=False
                    )
            if hasattr(self.gas_phase, "enthalpy_balances"):
                phase_list = self.gas_phase.properties.phase_list
                for (t, x), c in self.gas_phase.enthalpy_balances.items():
                    sf = iscale.min_scaling_factor(
                        [self.gas_phase._enthalpy_flow[t, x, p] for p in phase_list]
                    )
                    iscale.constraint_scaling_transform(c, 1e-3 * sf, overwrite=True)

            if hasattr(self, "solid_energy_accumulation_disc_eq"):
                for (t, x), c in self.solid_energy_accumulation_disc_eq.items():
                    iscale.constraint_scaling_transform(
                        c,
                        1e-3
                        * iscale.get_scaling_factor(
                            self.solid_energy_accumulation[t, x]
                        ),
                        overwrite=False,
                    )

            if hasattr(self.gas_phase, "energy_accumulation_disc_eq"):
                for (t, x, p), c in self.gas_phase.energy_accumulation_disc_eq.items():
                    iscale.constraint_scaling_transform(
                        c,
                        1e-3
                        * iscale.get_scaling_factor(
                            self.gas_phase.energy_accumulation[t, x, p]
                        ),
                        overwrite=True,
                    )

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {"Gas Inlet": self.gas_inlet, "Gas Outlet": self.gas_outlet},
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Bed Height"] = self.bed_height
        var_dict["Bed Area"] = self.bed_area
        var_dict["Gas Inlet Velocity"] = self.velocity_superficial_gas[time_point, 0]
        var_dict["Gas Outlet Velocity"] = self.velocity_superficial_gas[time_point, 1]

        return {"vars": var_dict}
