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
Base class for control volumes
"""

# Import Python libraries
import copy
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    Param,
    Reals,
    TransformationFactory,
    units as pyunits,
    Var,
    Reference,
    value,
)
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    ControlVolumeBlockData,
    FlowDirection,
    MaterialFlowBasis,
    MaterialBalanceType,
)
from idaes.core.util.exceptions import (
    BalanceTypeNotSupportedError,
    ConfigurationError,
    PropertyNotSupportedError,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.config import is_transformation_method, is_transformation_scheme
from idaes.core.util import scaling as iscale

import idaes.logger as idaeslog

__author__ = "Andrew Lee, Jaffer Ghouse"


_log = idaeslog.getLogger(__name__)

# TODO : Custom terms in material balances, other types of material balances
# Diffusion terms need to be added


# Enumerate options for area
class DistributedVars(Enum):
    variant = 0
    uniform = 1


@declare_process_block_class(
    "ControlVolume1DBlock",
    doc="""
    ControlVolume1DBlock is a specialized Pyomo block for IDAES control volume
    blocks discretized in one spatial direction, and contains instances of
    ControlVolume1DBlockData.

    ControlVolume1DBlock should be used for any control volume with a defined
    volume and distinct inlets and outlets where there is a single spatial
    domain parallel to the material flow direction. This encompases unit
    operations such as plug flow reactors and pipes.""",
)
class ControlVolume1DBlockData(ControlVolumeBlockData):
    """
    1-Dimensional ControlVolume Class

    This class forms the core of all 1-D IDAES models. It provides
    methods to build property and reaction blocks, and add mass, energy and
    momentum balances. The form of the terms used in these constraints is
    specified in the chosen property package.
    """

    CONFIG = ControlVolumeBlockData.CONFIG()
    CONFIG.declare(
        "area_definition",
        ConfigValue(
            default=DistributedVars.uniform,
            domain=In(DistributedVars),
            description="Argument for defining form of area variable",
            doc="""Argument defining whether area variable should be spatially
        variant or not. **default** - DistributedVars.uniform.
        **Valid values:** {
        DistributedVars.uniform - area does not vary across spatial domian,
        DistributedVars.variant - area can vary over the domain and is indexed
        by time and space.}""",
        ),
    )
    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default=None,
            domain=is_transformation_method,
            description="DAE transformation method",
            doc="""Method to use to transform domain. Must be a method recognised
by the Pyomo TransformationFactory.""",
        ),
    )
    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default=None,
            domain=is_transformation_scheme,
            description="DAE transformation scheme",
            doc="""Scheme to use when transforming domain. See Pyomo
documentation for supported schemes.""",
        ),
    )
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=None,
            domain=int,
            description="Number of finite elements",
            doc="""Number of finite elements to use in transformation (equivalent
to Pyomo nfe argument).""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=None,
            domain=int,
            description="Number of collocation points",
            doc="""Number of collocation points to use (equivalent to Pyomo ncp
argument).""",
        ),
    )

    def build(self):
        """
        Build method for ControlVolume1DBlock blocks.

        Returns:
            None
        """
        # Call build method from base class
        super(ControlVolume1DBlockData, self).build()

        self._validate_config_args()

    def _validate_config_args(self):
        # Validate DAE config arguments
        if self.config.transformation_method is None:
            raise ConfigurationError(
                "{} was not provided a value for the transformation_method"
                " configuration argument. Please provide a valid value.".format(
                    self.name
                )
            )

        if self.config.transformation_scheme is None:
            raise ConfigurationError(
                "{} was not provided a value for the transformation_scheme"
                " configuration argument. Please provide a valid value.".format(
                    self.name
                )
            )
        elif (
            self.config.transformation_method == "dae.finite_difference"
            and self.config.transformation_scheme not in ["BACKWARD", "FORWARD"]
        ) or (
            self.config.transformation_method == "dae.collocation"
            and self.config.transformation_scheme
            not in ["LAGRANGE-LEGENDRE", "LAGRANGE-RADAU"]
        ):
            raise ConfigurationError(
                "{} transformation_scheme configuration argument is not "
                "consistent with transformation_method argument. See Pyomo"
                " documentation for argument options.".format(self.name)
            )

    def add_geometry(
        self,
        length_domain=None,
        length_domain_set=None,
        length_var=None,
        flow_direction=FlowDirection.forward,
    ):
        """
        Method to create spatial domain and volume Var in ControlVolume.

        Args:
            length_domain - (optional) a ContinuousSet to use as the length
                            domain for the ControlVolume. If not provided, a
                            new ContinuousSet will be created (default=None).
                            ContinuousSet should be normalized to run between
                            0 and 1.
            length_domain_set - (optional) list of point to use to initialize
                            a new ContinuousSet if length_domain is not
                            provided (default = [0.0, 1.0]).
            length_var - (optional) external variable to use for the length of
                         the spatial domain,. If a variable is provided, a
                         reference will be made to this in place of the length
                         Var.
            flow_direction - argument indicating direction of material flow
                            relative to length domain. Valid values:
                                - FlowDirection.forward (default), flow goes
                                  from 0 to 1.
                                - FlowDirection.backward, flow goes from 1 to 0

        Returns:
            None
        """
        units = self.config.property_package.get_metadata().get_derived_units

        if length_domain is not None:
            # Validate domain and make a reference
            if isinstance(length_domain, ContinuousSet):
                add_object_reference(self, "length_domain", length_domain)
            else:
                raise ConfigurationError(
                    "{} length_domain argument must be a Pyomo "
                    "ContinuousSet object".format(self.name)
                )
        else:
            # Create new length domain
            if length_domain_set is None:
                length_domain_set = [0.0, 1.0]

            self.length_domain = ContinuousSet(
                bounds=(0.0, 1.0),
                initialize=length_domain_set,
                doc="Normalized length domain",
            )

        # Validated and create flow direction attribute
        if flow_direction in (flwd for flwd in FlowDirection):
            self._flow_direction = flow_direction
        else:
            raise ConfigurationError(
                "{} invalid value for flow_direction "
                "argument. Must be a FlowDirection Enum.".format(self.name)
            )
        if flow_direction is FlowDirection.forward:
            self._flow_direction_term = -1
        else:
            self._flow_direction_term = 1

        # Add geomerty variables and constraints
        if self.config.area_definition == DistributedVars.variant:
            self.area = Var(
                self.flowsheet().time,
                self.length_domain,
                initialize=1.0,
                doc="Cross-sectional area of Control Volume",
                units=units("area"),
            )
        else:
            self.area = Var(
                initialize=1.0,
                doc="Cross-sectional area of Control Volume",
                units=units("area"),
            )

        if length_var is not None:
            # Validate length_Var and add a reference
            if not isinstance(length_var, (Var, Param, Expression)):
                raise ConfigurationError(
                    f"{self.name} length_var must be a Pyomo Var, Param or "
                    "Expression."
                )
            elif length_var.is_indexed():
                raise ConfigurationError(
                    f"{self.name} length_var must be a scalar (unindexed) " "component."
                )
            add_object_reference(self, "length", length_var)
        else:
            self.length = Var(
                initialize=1.0, doc="Length of Control Volume", units=units("length")
            )

    def add_state_blocks(
        self, information_flow=FlowDirection.forward, has_phase_equilibrium=None
    ):
        """
        This method constructs the state blocks for the
        control volume.

        Args:
            information_flow: a FlowDirection Enum indicating whether
                               information flows from inlet-to-outlet or
                               outlet-to-inlet
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
            package_arguments: dict-like object of arguments to be passed to
                                state blocks as construction arguments
        Returns:
            None
        """
        if has_phase_equilibrium is None:
            raise ConfigurationError(
                "{} add_state_blocks method was not provided with a "
                "has_phase_equilibrium argument.".format(self.name)
            )
        elif has_phase_equilibrium not in [True, False]:
            raise ConfigurationError(
                "{} add_state_blocks method was provided with an invalid "
                "has_phase_equilibrium argument. Must be True or False".format(
                    self.name
                )
            )

        # d0 is config for defined state d1 is config for not defined state
        d0 = dict(**self.config.property_package_args)
        d0.update(has_phase_equilibrium=has_phase_equilibrium, defined_state=True)
        d1 = copy.copy(d0)
        d1["defined_state"] = False

        def idx_map(i):  # i = (t, x)
            if (
                information_flow == FlowDirection.forward
                and i[1] == self.length_domain.first()
            ):
                return 0
            elif (
                information_flow == FlowDirection.backward
                and i[1] == self.length_domain.last()
            ):
                return 0
            else:
                return 1

        self.properties = self.config.property_package.build_state_block(
            self.flowsheet().time,
            self.length_domain,
            doc="Material properties",
            initialize={0: d0, 1: d1},  # TODO: What if the domain has different bounds?
            idx_map=idx_map,
        )

    def add_reaction_blocks(self, has_equilibrium=None):
        """
        This method constructs the reaction block for the control volume.

        Args:
            has_equilibrium: indicates whether equilibrium calculations
                              will be required in reaction block
            package_arguments: dict-like object of arguments to be passed to
                                reaction block as construction arguments

        Returns:
            None
        """
        if has_equilibrium is None:
            raise ConfigurationError(
                "{} add_reaction_blocks method was not provided with a "
                "has_equilibrium argument.".format(self.name)
            )
        elif has_equilibrium not in [True, False]:
            raise ConfigurationError(
                "{} add_reaction_blocks method was provided with an "
                "invalid has_equilibrium argument. Must be True or False".format(
                    self.name
                )
            )

        # TODO : Should not have ReactionBlock at inlet
        tmp_dict = dict(**self.config.reaction_package_args)
        tmp_dict["state_block"] = self.properties
        tmp_dict["has_equilibrium"] = has_equilibrium

        self.reactions = self.config.reaction_package.build_reaction_block(
            self.flowsheet().time,
            self.length_domain,
            doc="Reaction properties in control volume",
            **tmp_dict,
        )  # TODO: Do we need something similar to above to skip equilibrium at bounds?

    def _add_material_balance_common(
        self,
        balance_type,
        has_rate_reactions,
        has_equilibrium_reactions,
        has_phase_equilibrium,
        has_mass_transfer,
        custom_molar_term,
        custom_mass_term,
    ):
        # Get dynamic and holdup flags from config block
        dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup

        component_list = self.properties.component_list
        phase_list = self.properties.phase_list
        pc_set = self.properties.phase_component_set

        # Check that reaction block exists if required
        if has_rate_reactions or has_equilibrium_reactions:
            try:
                rblock = self.reactions
            except AttributeError:
                raise ConfigurationError(
                    "{} does not contain a Reaction Block, but material "
                    "balances have been set to contain reaction terms. "
                    "Please construct a reaction block before adding "
                    "balance equations.".format(self.name)
                )

        if has_equilibrium_reactions:
            # Check that reaction block is set to calculate equilibrium
            for t in self.flowsheet().time:
                for x in self.length_domain:
                    if self.reactions[t, x].config.has_equilibrium is False:
                        raise ConfigurationError(
                            "{} material balance was set to include "
                            "equilibrium reactions, however the associated "
                            "ReactionBlock was not set to include equilibrium "
                            "constraints (has_equilibrium_reactions=False). "
                            "Please correct your configuration arguments.".format(
                                self.name
                            )
                        )

        if has_phase_equilibrium:
            # Check that state blocks are set to calculate equilibrium
            for t in self.flowsheet().time:
                for x in self.length_domain:
                    if not self.properties[t, x].config.has_phase_equilibrium:
                        raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the associated "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments.".format(self.name)
                        )

        # Get units from property package
        units = self.config.property_package.get_metadata().get_derived_units

        if units("length") is not None:
            if (
                self.properties[
                    self.flowsheet().time.first(), self.length_domain.first()
                ].get_material_flow_basis()
                == MaterialFlowBasis.molar
            ):
                holdup_l_units = units("amount") / units("length")
                flow_units = units("flow_mole")
                flow_l_units = units("flow_mole") / units("length")
            elif (
                self.properties[
                    self.flowsheet().time.first(), self.length_domain.first()
                ].get_material_flow_basis()
                == MaterialFlowBasis.mass
            ):
                holdup_l_units = units("mass") / units("length")
                flow_units = units("flow_mass")
                flow_l_units = units("flow_mass") / units("length")
            else:
                holdup_l_units = None
                flow_units = None
                flow_l_units = None
        else:
            holdup_l_units = None
            flow_units = None
            flow_l_units = None

        # Get units for accumulation term if required
        if self.config.dynamic:
            f_time_units = self.flowsheet().time_units
            if (f_time_units is None) ^ (units("time") is None):
                raise ConfigurationError(
                    "{} incompatible time unit specification between "
                    "flowsheet and property package. Either both must use "
                    "units, or neither.".format(self.name)
                )

            if f_time_units is None:
                acc_units = None
            elif (
                self.properties[
                    self.flowsheet().time.first(), self.length_domain.first()
                ].get_material_flow_basis()
                == MaterialFlowBasis.other
            ):
                acc_units = None
            else:
                acc_units = holdup_l_units / f_time_units

        # Check if reaction package exists, and get units
        if hasattr(self.config, "reaction_package"):
            if self.config.reaction_package is not None:
                if (
                    self.reactions[
                        self.flowsheet().time.first(), self.length_domain.first()
                    ].get_reaction_rate_basis()
                    == MaterialFlowBasis.molar
                ):
                    rxn_flow_l_units = units("flow_mole") / units("length")
                elif (
                    self.reactions[
                        self.flowsheet().time.first(), self.length_domain.first()
                    ].get_reaction_rate_basis()
                    == MaterialFlowBasis.mass
                ):
                    rxn_flow_l_units = units("flow_mass") / units("length")
                else:  # reaction basis not defined
                    rxn_flow_l_units = None
            else:  # reaction package is NoneType object
                rxn_flow_l_units = None
        else:  # reaction package not defined
            rxn_flow_l_units = None

        # Material holdup and accumulation
        if has_holdup:
            self.material_holdup = Var(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                domain=Reals,
                initialize=1.0,
                doc="Material holdup per unit length",
                units=holdup_l_units,
            )
        if dynamic:
            self.material_accumulation = DerivativeVar(
                self.material_holdup,
                wrt=self.flowsheet().time,
                doc="Material accumulation per unit length",
                units=acc_units,
            )

        # Create material balance terms as required
        # Flow terms and derivatives
        self._flow_terms = Var(
            self.flowsheet().time,
            self.length_domain,
            pc_set,
            initialize=1.0,
            doc="Flow terms for material balance equations",
            units=flow_units,
        )

        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            pc_set,
            doc="Material flow linking constraints",
        )
        def material_flow_linking_constraints(b, t, x, p, j):
            return b._flow_terms[t, x, p, j] == b.properties[
                t, x
            ].get_material_flow_terms(p, j)

        self.material_flow_dx = DerivativeVar(
            self._flow_terms,
            wrt=self.length_domain,
            doc="Parital derivative of material flow " "wrt to normalized length",
            units=flow_units,
        )

        # Kinetic reaction generation
        if has_rate_reactions:
            if not hasattr(self.config.reaction_package, "rate_reaction_idx"):
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of rate "
                    "reactions (rate_reaction_idx), thus does not support "
                    "rate-based reactions.".format(self.name)
                )
            self.rate_reaction_generation = Var(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                domain=Reals,
                initialize=0.0,
                doc="Amount of component generated in "
                "by kinetic reactions per unit length",
                units=rxn_flow_l_units,
            )  # use reaction package flow basis

        # Equilibrium reaction generation
        if has_equilibrium_reactions:
            if not hasattr(self.config.reaction_package, "equilibrium_reaction_idx"):
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of "
                    "equilibrium reactions (equilibrium_reaction_idx), thus "
                    "does not support equilibrium-based reactions.".format(self.name)
                )
            self.equilibrium_reaction_generation = Var(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                domain=Reals,
                initialize=0.0,
                doc="Amount of component generated by equilibrium "
                "reactions per unit length",
                units=rxn_flow_l_units,
            )  # use reaction package flow basis

        # Inherent reaction generation
        if self.properties.include_inherent_reactions:
            if not hasattr(self.config.property_package, "inherent_reaction_idx"):
                raise PropertyNotSupportedError(
                    "{} Property package does not contain a list of "
                    "inherent reactions (inherent_reaction_idx), but "
                    "include_inherent_reactions is True.".format(self.name)
                )
            self.inherent_reaction_generation = Var(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                domain=Reals,
                initialize=0.0,
                doc="Amount of component generated in control volume "
                "by inherent reactions",
                units=flow_l_units,
            )  # use property package flow basis

        # Phase equilibrium generation
        if has_phase_equilibrium and balance_type == MaterialBalanceType.componentPhase:
            if not hasattr(self.config.property_package, "phase_equilibrium_idx"):
                raise PropertyNotSupportedError(
                    "{} Property package does not contain a list of phase "
                    "equilibrium reactions (phase_equilibrium_idx), thus does "
                    "not support phase equilibrium.".format(self.name)
                )
            self.phase_equilibrium_generation = Var(
                self.flowsheet().time,
                self.length_domain,
                self.config.property_package.phase_equilibrium_idx,
                domain=Reals,
                initialize=0.0,
                doc="Amount of generation in unit by phase "
                "equilibria per unit length",
                units=flow_l_units,
            )  # use property package flow basis

        # Material transfer term
        if has_mass_transfer:
            self.mass_transfer_term = Var(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                domain=Reals,
                initialize=0.0,
                doc="Component material transfer into unit per unit " "length",
                units=flow_l_units,
            )

        # Create rules to substitute material balance terms
        # Accumulation term
        def accumulation_term(b, t, x, p, j):
            return (
                pyunits.convert(
                    b.material_accumulation[t, x, p, j], to_units=flow_l_units
                )
                if dynamic
                else 0
            )

        def kinetic_term(b, t, x, p, j):
            return b.rate_reaction_generation[t, x, p, j] if has_rate_reactions else 0

        def equilibrium_term(b, t, x, p, j):
            return (
                b.equilibrium_reaction_generation[t, x, p, j]
                if has_equilibrium_reactions
                else 0
            )

        def inherent_term(b, t, x, p, j):
            return (
                b.inherent_reaction_generation[t, x, p, j]
                if b.properties.include_inherent_reactions
                else 0
            )

        def phase_equilibrium_term(b, t, x, p, j):
            if (
                has_phase_equilibrium
                and balance_type == MaterialBalanceType.componentPhase
            ):
                sd = {}
                for r in b.config.property_package.phase_equilibrium_idx:
                    if b.config.property_package.phase_equilibrium_list[r][0] == j:
                        if (
                            b.config.property_package.phase_equilibrium_list[r][1][0]
                            == p
                        ):
                            sd[r] = 1
                        elif (
                            b.config.property_package.phase_equilibrium_list[r][1][1]
                            == p
                        ):
                            sd[r] = -1
                        else:
                            sd[r] = 0
                    else:
                        sd[r] = 0

                return sum(
                    b.phase_equilibrium_generation[t, x, r] * sd[r]
                    for r in b.config.property_package.phase_equilibrium_idx
                )
            else:
                return 0

        def transfer_term(b, t, x, p, j):
            return b.mass_transfer_term[t, x, p, j] if has_mass_transfer else 0

        # TODO: Need to set material_holdup = 0 for non-present component-phase
        # pairs. Not ideal, but needed to close DoF. Is there a better way?

        # Material Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                doc="Material holdup calculations",
            )
            def material_holdup_calculation(b, t, x, p, j):
                if (p, j) in pc_set:
                    return b.material_holdup[t, x, p, j] == (
                        b._area_func(t, x)
                        * self.phase_fraction[t, x, p]
                        * b.properties[t, x].get_material_density_terms(p, j)
                    )

        if has_rate_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.rate_reaction_extent = Var(
                self.flowsheet().time,
                self.length_domain,
                self.config.reaction_package.rate_reaction_idx,
                domain=Reals,
                initialize=0.0,
                doc="Extent of kinetic reactions at point x",
                units=rxn_flow_l_units,
            )  # use reaction package flow basis

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                doc="Kinetic reaction stoichiometry constraint",
            )
            def rate_reaction_stoichiometry_constraint(b, t, x, p, j):
                if (
                    b.config.transformation_scheme != "FORWARD"
                    and x == b.length_domain.first()
                ) or (
                    b.config.transformation_scheme == "FORWARD"
                    and x == b.length_domain.last()
                ):
                    return Constraint.Skip
                elif (p, j) in pc_set:
                    rparam = rblock[t, x].config.parameters
                    return b.rate_reaction_generation[t, x, p, j] == (
                        sum(
                            rparam.rate_reaction_stoichiometry[r, p, j]
                            * b.rate_reaction_extent[t, x, r]
                            for r in b.config.reaction_package.rate_reaction_idx
                        )
                    )
                else:
                    return Constraint.Skip

        if has_equilibrium_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.equilibrium_reaction_extent = Var(
                self.flowsheet().time,
                self.length_domain,
                self.config.reaction_package.equilibrium_reaction_idx,
                domain=Reals,
                initialize=0.0,
                doc="Extent of equilibrium reactions at point x",
                units=rxn_flow_l_units,
            )  # use reaction package flow basis

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                doc="Equilibrium reaction stoichiometry",
            )
            def equilibrium_reaction_stoichiometry_constraint(b, t, x, p, j):
                if (
                    b.config.transformation_scheme != "FORWARD"
                    and x == b.length_domain.first()
                ) or (
                    b.config.transformation_scheme == "FORWARD"
                    and x == b.length_domain.last()
                ):
                    return Constraint.Skip
                elif (p, j) in pc_set:
                    return b.equilibrium_reaction_generation[t, x, p, j] == (
                        sum(
                            rblock[
                                t, x
                            ].config.parameters.equilibrium_reaction_stoichiometry[
                                r, p, j
                            ]
                            * b.equilibrium_reaction_extent[t, x, r]
                            for r in b.config.reaction_package.equilibrium_reaction_idx
                        )
                    )
                else:
                    return Constraint.Skip

        if self.properties.include_inherent_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.inherent_reaction_extent = Var(
                self.flowsheet().time,
                self.length_domain,
                self.config.property_package.inherent_reaction_idx,
                domain=Reals,
                initialize=0.0,
                doc="Extent of inherent reactions at point x",
                units=flow_l_units,
            )  # use property package flow basis

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                doc="Inherent reaction stoichiometry",
            )
            def inherent_reaction_stoichiometry_constraint(b, t, x, p, j):
                if (
                    b.config.transformation_scheme != "FORWARD"
                    and x == b.length_domain.first()
                ) or (
                    b.config.transformation_scheme == "FORWARD"
                    and x == b.length_domain.last()
                ):
                    return Constraint.Skip
                elif (p, j) in pc_set:
                    return b.inherent_reaction_generation[t, x, p, j] == (
                        sum(
                            b.properties[
                                t, x
                            ].config.parameters.inherent_reaction_stoichiometry[r, p, j]
                            * b.inherent_reaction_extent[t, x, r]
                            for r in b.config.property_package.inherent_reaction_idx
                        )
                    )
                else:
                    return Constraint.Skip

        # Add custom terms and material balances
        if balance_type == MaterialBalanceType.componentPhase:

            def user_term_mol(b, t, x, p, j):
                if custom_molar_term is not None:
                    flow_basis = b.properties[t, x].get_material_flow_basis()
                    if flow_basis == MaterialFlowBasis.molar:
                        return custom_molar_term(t, x, p, j)
                    elif flow_basis == MaterialFlowBasis.mass:
                        try:
                            return (
                                custom_molar_term(t, x, p, j)
                                * b.properties[t, x].mw_comp[j]
                            )
                        except AttributeError:
                            raise PropertyNotSupportedError(
                                "{} property package does not support "
                                "molecular weight (mw), which is required for "
                                "using custom terms in material balances.".format(
                                    self.name
                                )
                            )
                    else:
                        raise ConfigurationError(
                            "{} contained a custom_molar_term argument, but "
                            "the property package used an undefined basis "
                            "(MaterialFlowBasis.other). Custom terms can "
                            "only be used when the property package declares "
                            "a molar or mass flow basis.".format(self.name)
                        )
                else:
                    return 0

            def user_term_mass(b, t, x, p, j):
                if custom_mass_term is not None:
                    flow_basis = b.properties[t, x].get_material_flow_basis()
                    if flow_basis == MaterialFlowBasis.mass:
                        return custom_mass_term(t, x, p, j)
                    elif flow_basis == MaterialFlowBasis.molar:
                        try:
                            return (
                                custom_mass_term(t, x, p, j)
                                / b.properties[t, x].mw_comp[j]
                            )
                        except AttributeError:
                            raise PropertyNotSupportedError(
                                "{} property package does not support "
                                "molecular weight (mw), which is required for "
                                "using custom terms in material balances.".format(
                                    self.name
                                )
                            )
                    else:
                        raise ConfigurationError(
                            "{} contained a custom_mass_term argument, but "
                            "the property package used an undefined basis "
                            "(MaterialFlowBasis.other). Custom terms can "
                            "only be used when the property package declares "
                            "a molar or mass flow basis.".format(self.name)
                        )
                else:
                    return 0

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                pc_set,
                doc="Material balances",
            )
            def material_balances(b, t, x, p, j):
                if (
                    b.config.transformation_scheme != "FORWARD"
                    and x == b.length_domain.first()
                ) or (
                    b.config.transformation_scheme == "FORWARD"
                    and x == b.length_domain.last()
                ):
                    return Constraint.Skip
                else:
                    if (p, j) in pc_set:
                        return b.length * accumulation_term(b, t, x, p, j) == (
                            b._flow_direction_term * b.material_flow_dx[t, x, p, j]
                            + b.length
                            * kinetic_term(b, t, x, p, j)
                            * b._rxn_rate_conv(t, x, j, has_rate_reactions)
                            + b.length * equilibrium_term(b, t, x, p, j)
                            + b.length * inherent_term(b, t, x, p, j)
                            + b.length * phase_equilibrium_term(b, t, x, p, j)
                            + b.length * transfer_term(b, t, x, p, j)
                            +
                            # b.area*diffusion_term(b, t, x, p, j)/b.length +
                            b.length * user_term_mol(b, t, x, p, j)
                            + b.length * user_term_mass(b, t, x, p, j)
                        )
                    else:
                        return Constraint.Skip

        elif balance_type == MaterialBalanceType.componentTotal:

            def user_term_mol(b, t, x, j):
                if custom_molar_term is not None:
                    flow_basis = b.properties[t, x].get_material_flow_basis()
                    if flow_basis == MaterialFlowBasis.molar:
                        return custom_molar_term(t, x, j)
                    elif flow_basis == MaterialFlowBasis.mass:
                        try:
                            return (
                                custom_molar_term(t, x, j)
                                * b.properties[t, x].mw_comp[j]
                            )
                        except AttributeError:
                            raise PropertyNotSupportedError(
                                "{} property package does not support "
                                "molecular weight (mw), which is required for "
                                "using custom terms in material balances.".format(
                                    self.name
                                )
                            )
                    else:
                        raise ConfigurationError(
                            "{} contained a custom_molar_term argument, but "
                            "the property package used an undefined basis "
                            "(MaterialFlowBasis.other). Custom terms can "
                            "only be used when the property package declares "
                            "a molar or mass flow basis.".format(self.name)
                        )
                else:
                    return 0

            def user_term_mass(b, t, x, j):
                if custom_mass_term is not None:
                    flow_basis = b.properties[t, x].get_material_flow_basis()
                    if flow_basis == MaterialFlowBasis.mass:
                        return custom_mass_term(t, x, j)
                    elif flow_basis == MaterialFlowBasis.molar:
                        try:
                            return (
                                custom_mass_term(t, x, j)
                                / b.properties[t, x].mw_comp[j]
                            )
                        except AttributeError:
                            raise PropertyNotSupportedError(
                                "{} property package does not support "
                                "molecular weight (mw), which is required for "
                                "using custom terms in material balances.".format(
                                    self.name
                                )
                            )
                    else:
                        raise ConfigurationError(
                            "{} contained a custom_mass_term argument, but "
                            "the property package used an undefined basis "
                            "(MaterialFlowBasis.other). Custom terms can "
                            "only be used when the property package declares "
                            "a molar or mass flow basis.".format(self.name)
                        )
                else:
                    return 0

            # Add component balances
            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                component_list,
                doc="Material balances",
            )
            def material_balances(b, t, x, j):
                if (
                    b.config.transformation_scheme != "FORWARD"
                    and x == b.length_domain.first()
                ) or (
                    b.config.transformation_scheme == "FORWARD"
                    and x == b.length_domain.last()
                ):
                    return Constraint.Skip
                else:
                    cplist = []
                    for p in phase_list:
                        if (p, j) in pc_set:
                            cplist.append(p)
                    return b.length * sum(
                        accumulation_term(b, t, x, p, j) for p in cplist
                    ) == b._flow_direction_term * sum(
                        b.material_flow_dx[t, x, p, j] for p in cplist
                    ) + b.length * sum(
                        kinetic_term(b, t, x, p, j) for p in cplist
                    ) * b._rxn_rate_conv(
                        t, x, j, has_rate_reactions
                    ) + b.length * sum(
                        equilibrium_term(b, t, x, p, j) for p in cplist
                    ) + b.length * sum(
                        inherent_term(b, t, x, p, j) for p in cplist
                    ) + b.length * sum(
                        transfer_term(b, t, x, p, j) for p in cplist
                    ) + b.length * user_term_mol(
                        b, t, x, j
                    ) + b.length * user_term_mass(
                        b, t, x, j
                    )

    def add_phase_component_balances(
        self,
        has_rate_reactions=False,
        has_equilibrium_reactions=False,
        has_phase_equilibrium=False,
        has_mass_transfer=False,
        custom_molar_term=None,
        custom_mass_term=None,
    ):
        """
        This method constructs a set of 1D material balances indexed by time,
        length, phase and component.

        Args:
            has_rate_reactions: whether default generation terms for rate
                    reactions should be included in material balances
            has_equilibrium_reactions: whether generation terms should for
                    chemical equilibrium reactions should be included in
                    material balances
            has_phase_equilibrium: whether generation terms should for phase
                    equilibrium behaviour should be included in material
                    balances
            has_mass_transfer: whether generic mass transfer terms should be
                    included in material balances
            custom_molar_term: a Pyomo Expression representing custom terms to
                    be included in material balances on a molar basis.
                    Expression must be indexed by time, length domain, phase
                    list and component list
            custom_mass_term: a Pyomo Expression representing custom terms to
                    be included in material balances on a mass basis.
                    Expression must be indexed by time, length domain, phase
                    list and component list

        Returns:
            Constraint object representing material balances
        """
        self._add_material_balance_common(
            balance_type=MaterialBalanceType.componentPhase,
            has_rate_reactions=has_rate_reactions,
            has_equilibrium_reactions=has_equilibrium_reactions,
            has_phase_equilibrium=has_phase_equilibrium,
            has_mass_transfer=has_mass_transfer,
            custom_molar_term=custom_molar_term,
            custom_mass_term=custom_mass_term,
        )

        return self.material_balances

    def add_total_component_balances(
        self,
        has_rate_reactions=False,
        has_equilibrium_reactions=False,
        has_phase_equilibrium=False,
        has_mass_transfer=False,
        custom_molar_term=None,
        custom_mass_term=None,
    ):
        """
        This method constructs a set of 1D material balances indexed by time
        length and component.

        Args:
            has_rate_reactions: whether default generation terms for rate
                    reactions should be included in material balances
            has_equilibrium_reactions: whether generation terms should for
                    chemical equilibrium reactions should be included in
                    material balances
            has_phase_equilibrium: whether generation terms should for phase
                    equilibrium behaviour should be included in material
                    balances
            has_mass_transfer: whether generic mass transfer terms should be
                    included in material balances
            custom_molar_term: a Pyomo Expression representing custom terms to
                    be included in material balances on a molar basis.
                    Expression must be indexed by time, length domain and
                    component list
            custom_mass_term: a Pyomo Expression representing custom terms to
                    be included in material balances on a mass basis.
                    Expression must be indexed by time, length domain and
                    component list

        Returns:
            Constraint object representing material balances
        """
        self._add_material_balance_common(
            balance_type=MaterialBalanceType.componentTotal,
            has_rate_reactions=has_rate_reactions,
            has_equilibrium_reactions=has_equilibrium_reactions,
            has_phase_equilibrium=has_phase_equilibrium,
            has_mass_transfer=has_mass_transfer,
            custom_molar_term=custom_molar_term,
            custom_mass_term=custom_mass_term,
        )

        return self.material_balances

    def add_total_element_balances(
        self,
        has_rate_reactions=False,
        has_equilibrium_reactions=False,
        has_phase_equilibrium=False,
        has_mass_transfer=False,
        custom_elemental_term=None,
    ):
        """
        This method constructs a set of 1D element balances indexed by time and
        length.

        Args:
            has_rate_reactions - whether default generation terms for rate
                    reactions should be included in material balances
            has_equilibrium_reactions - whether generation terms should for
                    chemical equilibrium reactions should be included in
                    material balances
            has_phase_equilibrium - whether generation terms should for phase
                    equilibrium behaviour should be included in material
                    balances
            has_mass_transfer - whether generic mass transfer terms should be
                    included in material balances
            custom_elemental_term - a Pyomo Expression representing custom
                    terms to be included in material balances on a molar
                    elemental basis. Expression must be indexed by time, length
                    and element list

        Returns:
            Constraint object representing material balances
        """
        # Get dynamic and holdup flags from config block
        dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup

        component_list = self.properties.component_list
        pc_set = self.properties.phase_component_set

        # Check that property package supports element balances
        if not hasattr(self.config.property_package, "element_list"):
            raise PropertyNotSupportedError(
                "{} property package provided does not contain a list of "
                "elements (element_list), and thus does not support "
                "elemental material balances. Please choose another type "
                "of material balance or a property package which supports "
                "elemental balances."
            )

        # Check validity of arguments to write the total elemental balance
        if has_rate_reactions:
            raise ConfigurationError(
                "{} add_total_element_balances method as provided with "
                "argument has_rate_reactions = True. Total element "
                "balances do not support rate based reactions, "
                "please correct your configuration arguments".format(self.name)
            )

        if has_equilibrium_reactions:
            raise ConfigurationError(
                "{} add_total_element_balances method as provided with "
                "argument has_equilibrium_reactions = True. Total element "
                "balances do not support equilibrium based reactions, "
                "please correct your configuration arguments".format(self.name)
            )

        if has_phase_equilibrium:
            # Check that state blocks are set to calculate equilibrium
            raise ConfigurationError(
                "{} add_total_element_balances method as provided with "
                "argument has_phase_equilibrium = True. Total element "
                "balances do not support equilibrium based reactions, "
                "please correct your configuration arguments".format(self.name)
            )

        # Get units from property package
        units = self.config.property_package.get_metadata().get_derived_units

        if units("amount") is not None:
            amount_l_units = units("amount") / units("length")
            flow_units = units("flow_mole")
            flow_l_units = units("flow_mole") / units("length")
        else:
            amount_l_units = None
            flow_units = None
            flow_l_units = None

        # Get units for accumulation term if required
        if self.config.dynamic:
            f_time_units = self.flowsheet().time_units
            if (f_time_units is None) ^ (units("time") is None):
                raise ConfigurationError(
                    "{} incompatible time unit specification between "
                    "flowsheet and property package. Either both must use "
                    "units, or neither.".format(self.name)
                )

            if f_time_units is None:
                acc_units = None
            else:
                acc_units = amount_l_units / f_time_units

        # Identify linearly dependent elements
        # It is possible for there to be linearly dependent element balances
        # e.g. if a single species is the only source of two different elements
        linearly_dependent = []

        # Get a representative time point
        rtime = self.flowsheet().time.first()
        rdomain = self.length_domain.first()

        # For each component in the material, search for elements which are
        # unique to it
        for i in component_list:
            unique_elements = []
            for e in self.config.property_package.element_list:
                if self.properties[rtime, rdomain].params.element_comp[i][e] != 0:
                    # Assume unique until shown otherwise
                    unique = True

                    for j in component_list:
                        if j == i:
                            continue

                        # If element appears in any other component, not unique
                        if (
                            self.properties[rtime, rdomain].params.element_comp[j][e]
                            != 0
                        ):
                            unique = False

                    if unique:
                        unique_elements.append(e)

            # If more than 1 unique element, they are linearly dependent
            if len(unique_elements) > 1:
                # Add all but the first to the list of linearly dependent
                linearly_dependent.extend(unique_elements[1:])

        # Set indexing set for element balances
        if len(linearly_dependent) == 0:
            # No linearly depednet equations, so use full element list
            e_index = self.config.property_package.element_list
        else:
            # Otherwise, use only non-dependent elements, and log a message
            _log.info_low(
                "{} detected linearly dependent element balance "
                "equations. Element balances will NOT be written "
                "for the following elements: {}".format(self.name, linearly_dependent)
            )
            e_index = self.config.property_package.element_list - linearly_dependent

        # Add Material Balance terms
        if has_holdup:
            self.element_holdup = Var(
                self.flowsheet().time,
                self.length_domain,
                self.config.property_package.element_list,
                domain=Reals,
                initialize=1.0,
                doc="Elemental holdup per unit length",
                units=amount_l_units,
            )

        if dynamic:
            self.element_accumulation = DerivativeVar(
                self.element_holdup,
                wrt=self.flowsheet().time,
                doc="Elemental accumulation per unit length",
                units=acc_units,
            )

        self.elemental_flow_term = Var(
            self.flowsheet().time,
            self.length_domain,
            self.config.property_package.element_list,
            initialize=1.0,
            doc="Elemental flow terms",
            units=flow_units,
        )

        # Method to convert mass flow basis to mole flow basis
        def conv_factor(b, t, x, j):
            flow_basis = b.properties[t, x].get_material_flow_basis()
            if flow_basis == MaterialFlowBasis.molar:
                return 1
            elif flow_basis == MaterialFlowBasis.mass:
                return 1 / b.properties[t, x].mw_comp[j]
            else:
                raise BalanceTypeNotSupportedError(
                    "{} property package MaterialFlowBasis == 'other'. Cannot "
                    "automatically generate elemental balances.".format(self.name)
                )

        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            self.config.property_package.element_list,
            doc="Elemental flow constraints",
        )
        def elemental_flow_constraint(b, t, x, e):
            return b.elemental_flow_term[t, x, e] == (
                sum(
                    conv_factor(b, t, x, j)
                    * b.properties[t, x].get_material_flow_terms(p, j)
                    * b.properties[t, x].config.parameters.element_comp[j][e]
                    for p, j in pc_set
                )
            )

        self.elemental_flow_dx = DerivativeVar(
            self.elemental_flow_term,
            wrt=self.length_domain,
            doc="Partial derivative of " "elemental flow wrt normalized " "length",
            units=flow_units,
        )

        # Create material balance terms as needed
        if has_mass_transfer:
            self.elemental_mass_transfer_term = Var(
                self.flowsheet().time,
                self.length_domain,
                e_index,
                domain=Reals,
                initialize=0.0,
                doc="Element material transfer into unit per unit " "length",
                units=flow_l_units,
            )

        # Create rules to substitute material balance terms
        # Accumulation term
        def accumulation_term(b, t, x, e):
            return (
                pyunits.convert(b.element_accumulation[t, x, e], to_units=flow_l_units)
                if dynamic
                else 0
            )

        # Mass transfer term
        def transfer_term(b, t, x, e):
            return b.elemental_mass_transfer_term[t, x, e] if has_mass_transfer else 0

        # Custom term
        def user_term(t, x, e):
            if custom_elemental_term is not None:
                return custom_elemental_term(t, x, e)
            else:
                return 0

        # Element balances
        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            e_index,
            doc="Elemental material balances",
        )
        def element_balances(b, t, x, e):
            if (
                b.config.transformation_scheme != "FORWARD"
                and x == b.length_domain.first()
            ) or (
                b.config.transformation_scheme == "FORWARD"
                and x == b.length_domain.last()
            ):
                return Constraint.Skip
            else:
                return b.length * accumulation_term(b, t, x, e) == (
                    b._flow_direction_term * b.elemental_flow_dx[t, x, e]
                    + b.length * transfer_term(b, t, x, e)
                    + b.length * user_term(t, x, e)
                )  # +
                # TODO : Add diffusion terms
                # b.area*diffusion_term(b, t, x, e)/b.length)

        # Elemental Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                self.config.property_package.element_list,
                doc="Elemental holdup calculation",
            )
            def elemental_holdup_calculation(b, t, x, e):
                return b.element_holdup[t, x, e] == (
                    b._area_func(t, x)
                    * sum(
                        conv_factor(b, t, x, j)
                        * b.phase_fraction[t, x, p]
                        * b.properties[t, x].get_material_density_terms(p, j)
                        * b.properties[t, x].config.parameters.element_comp[j][e]
                        for p, j in pc_set
                    )
                )

        return self.element_balances

    def add_total_material_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
            "{} OD control volumes do not support "
            "add_total_material_balances (yet).".format(self.name)
        )

    def add_total_enthalpy_balances(
        self,
        has_heat_of_reaction=False,
        has_heat_transfer=False,
        has_work_transfer=False,
        has_enthalpy_transfer=False,
        custom_term=None,
    ):
        """
        This method constructs a set of 1D enthalpy balances indexed by time
        and phase.

        Args:
            has_heat_of_reaction - whether terms for heat of reaction should
                    be included in enthalpy balance
            has_heat_transfer - whether terms for heat transfer should be
                    included in enthalpy balances
            has_work_transfer - whether terms for work transfer should be
                    included in enthalpy balances
            has_enthalpy_transfer - whether terms for enthalpy transfer due to
                    mass transfer should be included in enthalpy balance. This
                    should generally be the same as the has_mass_transfer
                    argument in the material balance methods
            custom_term - a Python method which returns Pyomo expressions representing
                    custom terms to be included in enthalpy balances.
                    Method should accept time, length and phase list as arguments.

        Returns:
            Constraint object representing enthalpy balances
        """
        # Get dynamic and holdup flags from config block
        dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup

        phase_list = self.properties.phase_list

        # Test for components that must exist prior to calling this method
        if has_heat_of_reaction:
            if not (
                hasattr(self, "rate_reaction_extent")
                or hasattr(self, "equilibrium_reaction_extent")
            ):
                raise ConfigurationError(
                    "{} extent of reaction terms must exist in order to "
                    "have heat of reaction terms. Please ensure that "
                    "add_material_balance (or equivalent) is called before"
                    " adding energy balances.".format(self.name)
                )

        # Get units from property package
        units = self.config.property_package.get_metadata().get_derived_units

        if units("energy") is not None:
            energy_l_units = units("energy") / units("length")
            power_l_units = units("power") / units("length")
        else:
            energy_l_units = None
            power_l_units = None

        # Get units for accumulation term if required
        if self.config.dynamic:
            f_time_units = self.flowsheet().time_units
            if (f_time_units is None) ^ (units("time") is None):
                raise ConfigurationError(
                    "{} incompatible time unit specification between "
                    "flowsheet and property package. Either both must use "
                    "units, or neither.".format(self.name)
                )

            if f_time_units is None:
                acc_units = None
            else:
                acc_units = energy_l_units / f_time_units

        # Create variables
        self._enthalpy_flow = Var(
            self.flowsheet().time,
            self.length_domain,
            phase_list,
            initialize=1.0,
            doc="Enthalpy flow terms",
            units=units("power"),
        )

        @self.Constraint(
            self.flowsheet().time,
            self.length_domain,
            phase_list,
            doc="Enthalpy flow linking constraints",
        )
        def enthalpy_flow_linking_constraint(b, t, x, p):
            return b._enthalpy_flow[t, x, p] == b.properties[
                t, x
            ].get_enthalpy_flow_terms(p)

        self.enthalpy_flow_dx = DerivativeVar(
            self._enthalpy_flow,
            wrt=self.length_domain,
            doc="Partial derivative of " "enthalpy flow wrt normalized " "length",
            units=units("power"),
        )

        if has_holdup:
            self.energy_holdup = Var(
                self.flowsheet().time,
                self.length_domain,
                phase_list,
                domain=Reals,
                initialize=1.0,
                doc="Enthalpy holdup per unit length",
                units=energy_l_units,
            )

        if dynamic is True:
            self.energy_accumulation = DerivativeVar(
                self.energy_holdup,
                wrt=self.flowsheet().time,
                doc="Energy accumulation per unit length",
                units=acc_units,
            )

        # Create energy balance terms as needed
        # Heat transfer term
        if has_heat_transfer:
            self.heat = Var(
                self.flowsheet().time,
                self.length_domain,
                domain=Reals,
                initialize=0.0,
                doc="Heat transferred per unit length",
                units=power_l_units,
            )

        # Work transfer
        if has_work_transfer:
            self.work = Var(
                self.flowsheet().time,
                self.length_domain,
                domain=Reals,
                initialize=0.0,
                doc="Work transfered per unit length",
                units=power_l_units,
            )

        # Enthalpy transfer
        if has_enthalpy_transfer:
            self.enthalpy_transfer = Var(
                self.flowsheet().time,
                self.length_domain,
                domain=Reals,
                initialize=0.0,
                doc="Enthalpy transferred due to mass transfer per unit length",
                units=power_l_units,
            )

        # Heat of Reaction
        if has_heat_of_reaction:

            @self.Expression(
                self.flowsheet().time,
                self.length_domain,
                doc="Heat of reaction term at point x",
            )
            def heat_of_reaction(b, t, x):
                if hasattr(self, "rate_reaction_extent"):
                    rate_heat = -sum(
                        b.rate_reaction_extent[t, x, r] * b.reactions[t, x].dh_rxn[r]
                        for r in self.config.reaction_package.rate_reaction_idx
                    )
                else:
                    rate_heat = 0

                if hasattr(self, "equilibrium_reaction_extent"):
                    equil_heat = -sum(
                        b.equilibrium_reaction_extent[t, x, e]
                        * b.reactions[t, x].dh_rxn[e]
                        for e in self.config.reaction_package.equilibrium_reaction_idx
                    )
                else:
                    equil_heat = 0

                return rate_heat + equil_heat

        # Create rules to substitute energy balance terms
        # Accumulation term
        def accumulation_term(b, t, x, p):
            return (
                pyunits.convert(b.energy_accumulation[t, x, p], to_units=power_l_units)
                if dynamic
                else 0
            )

        def heat_term(b, t, x):
            return b.heat[t, x] if has_heat_transfer else 0

        def work_term(b, t, x):
            return b.work[t, x] if has_work_transfer else 0

        def enthalpy_transfer_term(b, t, x):
            return b.enthalpy_transfer[t, x] if has_enthalpy_transfer else 0

        def rxn_heat_term(b, t, x):
            return b.heat_of_reaction[t, x] if has_heat_of_reaction else 0

        # Custom term
        def user_term(t, x):
            if custom_term is not None:
                return custom_term(t, x)
            else:
                return 0

        # Energy balance equation
        @self.Constraint(
            self.flowsheet().time, self.length_domain, doc="Energy balances"
        )
        def enthalpy_balances(b, t, x):
            if (
                b.config.transformation_scheme != "FORWARD"
                and x == b.length_domain.first()
            ) or (
                b.config.transformation_scheme == "FORWARD"
                and x == b.length_domain.last()
            ):
                return Constraint.Skip
            else:
                return b.length * sum(
                    accumulation_term(b, t, x, p) for p in phase_list
                ) == (
                    b._flow_direction_term
                    * sum(b.enthalpy_flow_dx[t, x, p] for p in phase_list)
                    + b.length * heat_term(b, t, x)
                    + b.length * work_term(b, t, x)
                    + b.length * enthalpy_transfer_term(b, t, x)
                    + b.length * rxn_heat_term(b, t, x)
                    + b.length * user_term(t, x)
                )
                # TODO : Add conduction/dispersion term

        # Energy Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                phase_list,
                doc="Enthalpy holdup constraint",
            )
            def energy_holdup_calculation(b, t, x, p):
                return b.energy_holdup[t, x, p] == (
                    b._area_func(t, x)
                    * self.phase_fraction[t, x, p]
                    * b.properties[t, x].get_energy_density_terms(p)
                )

        return self.enthalpy_balances

    def add_phase_enthalpy_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
            "{} OD control volumes do not support "
            "add_phase_enthalpy_balances.".format(self.name)
        )

    def add_phase_energy_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
            "{} OD control volumes do not support "
            "add_phase_energy_balances.".format(self.name)
        )

    def add_total_energy_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
            "{} OD control volumes do not support "
            "add_total_energy_balances.".format(self.name)
        )

    def add_total_pressure_balances(self, has_pressure_change=False, custom_term=None):
        """
        This method constructs a set of 1D pressure balances indexed by time.

        Args:
            has_pressure_change - whether terms for pressure change should be
                    included in enthalpy balances
            custom_term - a Pyomo Expression representing custom terms to
                    be included in pressure balances.
                    Expression must be indexed by time and length domain
            custom_term - a Python method which returns Pyomo expressions representing
                    custom terms to be included in enthalpy balances.
                    Method should accept time and length as arguments.

        Returns:
            Constraint object representing pressure balances
        """

        # Get units from property package
        units = self.config.property_package.get_metadata().get_derived_units

        if units("pressure") is not None:
            pressure_l_units = units("pressure") / units("length")
        else:
            pressure_l_units = None

        # Create dP/dx terms
        self.pressure = Reference(self.properties[:, :].pressure)

        self.pressure_dx = DerivativeVar(
            self.pressure,
            wrt=self.length_domain,
            doc="Partial derivative of pressure wrt " "normalized length domain",
            units=units("pressure"),
        )

        # Add Momentum Balance Variables as necessary
        if has_pressure_change:
            self.deltaP = Var(
                self.flowsheet().time,
                self.length_domain,
                domain=Reals,
                initialize=0.0,
                doc="Pressure difference per unit length " "of domain",
                units=pressure_l_units,
            )

        # Create rules to substitute energy balance terms
        # Pressure change term
        def deltaP_term(b, t, x):
            return b.deltaP[t, x] if has_pressure_change else 0

        # Custom term
        def user_term(t, x):
            if custom_term is not None:
                return custom_term(t, x)
            else:
                return 0

        # Momentum balance equation
        @self.Constraint(
            self.flowsheet().time, self.length_domain, doc="Momentum balance"
        )
        def pressure_balance(b, t, x):
            if (
                b.config.transformation_scheme != "FORWARD"
                and x == b.length_domain.first()
            ) or (
                b.config.transformation_scheme == "FORWARD"
                and x == b.length_domain.last()
            ):
                return Constraint.Skip
            else:
                return 0 == (
                    b._flow_direction_term * b.pressure_dx[t, x]
                    + b.length * deltaP_term(b, t, x)
                    + b.length * user_term(t, x)
                )

        return self.pressure_balance

    def add_phase_pressure_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
            "{} OD control volumes do not support "
            "add_phase_pressure_balances.".format(self.name)
        )

    def add_phase_momentum_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
            "{} OD control volumes do not support "
            "add_phase_momentum_balances.".format(self.name)
        )

    def add_total_momentum_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
            "{} OD control volumes do not support "
            "add_total_momentum_balances.".format(self.name)
        )

    def apply_transformation(self):
        """
        Method to apply DAE transformation to the Control Volume length domain.
        Transformation applied will be based on the Control Volume
        configuration arguments.
        """
        if self.length_domain.parent_block() != self:
            raise ConfigurationError(
                "{} tried to apply a DAE transformation to an external "
                "domain. To avoid complications, the apply_transformation "
                "method only supports transformation of local domains.".format(
                    self.name
                )
            )

        if self.config.finite_elements is None:
            raise ConfigurationError(
                "{} was not provided a value for the finite_elements"
                " configuration argument. Please provide a valid value.".format(
                    self.name
                )
            )

        if (
            self.config.collocation_points is None
            and self.config.transformation_method == "dae.collocation"
        ):
            raise ConfigurationError(
                "{} was not provided a value for the collocation_points"
                " configuration argument. Please provide a valid value.".format(
                    self.name
                )
            )

        if self.config.transformation_method == "dae.finite_difference":
            self.discretizer = TransformationFactory(self.config.transformation_method)
            self.discretizer.apply_to(
                self,
                nfe=self.config.finite_elements,
                wrt=self.length_domain,
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
        else:
            raise ConfigurationError(
                "{} unrecognised transformation_method, "
                "must match one of the Transformations "
                "supported by Pyomo's "
                "TransformationFactory.".format(self.name)
            )

    def model_check(blk):
        """
        This method executes the model_check methods on the associated state
        blocks (if they exist). This method is generally called by a unit model
        as part of the unit's model_check method.

        Args:
            None

        Returns:
            None
        """
        # Try property block model check
        for t in blk.flowsheet().time:
            for x in blk.length_domain:
                try:
                    blk.properties[t, x].model_check()
                except AttributeError:
                    _log.warning(
                        "{} ControlVolume StateBlock has no "
                        "model checks. To correct this, add a model_check"
                        " method to the associated StateBlock class.".format(blk.name)
                    )

                try:
                    blk.reactions[t, x].model_check()
                except AttributeError:
                    _log.warning(
                        "{} ControlVolume outlet reaction block has no "
                        "model check. To correct this, add a "
                        "model_check method to the associated "
                        "ReactionBlock class.".format(blk.name)
                    )

    def initialize(
        blk,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        optarg=None,
        solver=None,
        hold_state=True,
    ):
        """
        Initialization routine for 1D control volume.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None)
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - True. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization else the release state is
            triggered.
        """
        if optarg is None:
            optarg = {}

        # Get inlet state if not provided
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="control_volume")

        # Get source block
        if blk._flow_direction == FlowDirection.forward:
            source_idx = blk.length_domain.first()
        else:
            source_idx = blk.length_domain.last()
        source = blk.properties[blk.flowsheet().time.first(), source_idx]

        # Fix source state and get state_args if not provided
        source_flags = {}
        if state_args is None:
            # No state args, create whilst fixing vars
            state_args = {}
            # Should be checking flow direction
            state_dict = source.define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    source_flags[k] = {}
                    for m in state_dict[k].keys():
                        source_flags[k][m] = state_dict[k][m].fixed
                        if state_dict[k][m].value is not None:
                            state_dict[k][m].fix()
                            state_args[k][m] = state_dict[k][m].value
                        else:
                            raise Exception(
                                "State variables have not been "
                                "fixed nor have been given "
                                "initial values."
                            )
                else:
                    source_flags[k] = state_dict[k].fixed
                    if state_dict[k].value is not None:
                        state_dict[k].fix()
                        state_args[k] = state_dict[k].value
                    else:
                        raise Exception(
                            "State variables have not been "
                            "fixed nor have been given "
                            "initial values."
                        )
        else:
            # State  args provided
            state_dict = source.define_port_members()

            for k in state_dict.keys():
                source_flags[k] = {}
                if state_dict[k].is_indexed():
                    for m in state_dict[k].keys():
                        source_flags[k][m] = state_dict[k][m].fixed
                        if not state_dict[k][m].fixed:
                            state_dict[k][m].fix(state_args[k][m])
                else:
                    source_flags[k] = state_dict[k].fixed
                    if state_dict[k].value is not None:
                        state_dict[k].fix()
                        if not state_dict[k].fixed:
                            state_dict[k].fix(state_args[k])

        # Initialize state blocks
        flags = blk.properties.initialize(
            state_args=state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        try:
            # TODO: setting state_vars_fixed may not work for heterogeneous
            # systems where a second control volume is involved, as we cannot
            # assume those state vars are also fixed. For now, heterogeneous
            # reactions should ignore the state_vars_fixed argument and always
            # check their state_vars.
            blk.reactions.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_vars_fixed=True,
            )
        except AttributeError:
            pass

        init_log.info("Initialization Complete")

        # Unfix state variables except for source block
        blk.properties.release_state(flags)

        if hold_state is True:
            return source_flags
        else:
            blk.release_state(source_flags)

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state = True.
            outlvl : sets output level of logging

        Returns:
            None
        """
        # Get source block
        if blk._flow_direction == FlowDirection.forward:
            source_idx = blk.length_domain.first()
        else:
            source_idx = blk.length_domain.last()
        source = blk.properties[blk.flowsheet().time.first(), source_idx]

        # Set fixed attribute on state vars based on flags
        state_dict = source.define_port_members()

        for k in state_dict.keys():
            if state_dict[k].is_indexed():
                for m in state_dict[k].keys():
                    state_dict[k][m].fixed = flags[k][m]
            else:
                state_dict[k].fixed = flags[k]

    def _add_phase_fractions(self):
        """
        This method constructs the phase_fraction variables for the control
        volume, and the associated constraint on the sum of phase_fractions
        == 1. For systems with only one phase, phase_fraction is created as a
        Pyomo Expression with a value of 1.

        Args:
            None

        Returns:
            None
        """
        phase_list = self.properties.phase_list

        if len(phase_list) > 1:
            self.phase_fraction = Var(
                self.flowsheet().time,
                self.length_domain,
                phase_list,
                initialize=1 / len(phase_list),
                doc="Volume fraction of holdup by phase",
            )

            @self.Constraint(
                self.flowsheet().time,
                self.length_domain,
                doc="Sum of phase fractions == 1",
            )
            def sum_of_phase_fractions(b, t, x):
                return 1 == sum(b.phase_fraction[t, x, p] for p in phase_list)

        else:

            @self.Expression(
                self.flowsheet().time,
                self.length_domain,
                phase_list,
                doc="Volume fraction of holdup by phase",
            )
            def phase_fraction(self, t, x, p):
                return 1

    def _area_func(b, t, x):
        if b.config.area_definition == DistributedVars.uniform:
            return b.area
        return b.area[t, x]

    def report(self, time_point=0, dof=False, ostream=None, prefix=""):
        """
        No report method defined for ControlVolume1D class. This is due to the
        difficulty of presenting spatially discretized data in a readable form
        without plotting.
        """
        raise NotImplementedError(
            """
                Due ot the difficultly in presenting spatially distributed data
                in a clean format, ControlVolume1D does not currently support
                the report method."""
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        phase_list = self.properties.phase_list
        phase_component_set = self.properties.phase_component_set

        # Default scale factors
        # If the parent component of an indexed component has a scale factor,
        # but some of the data objects don't, propagate the indexed component
        # scale factor to the missing scaling factors.
        iscale.propagate_indexed_component_scaling_factors(self)

        # Set scaling for geometry variables
        if hasattr(self, "area"):
            for v in self.area.values():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(self.area, default=1, warning=True)
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "length"):
            if iscale.get_scaling_factor(self.length) is None:
                iscale.set_scaling_factor(self.length, 1)

        if hasattr(self, "phase_fraction"):
            for v in self.phase_fraction.values():
                if iscale.get_scaling_factor(v) is None:
                    # phase fraction typically between 0.1 and 1
                    iscale.set_scaling_factor(v, 10)

        # Set scaling factors for common material balance variables
        if hasattr(self, "material_holdup"):
            for (t, x, p, j), v in self.material_holdup.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(self._area_func(t, x))
                    sf *= iscale.get_scaling_factor(self.phase_fraction[t, x, p])
                    sf *= iscale.get_scaling_factor(
                        self.properties[t, x].get_material_density_terms(p, j),
                        default=1,
                        warning=True,
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "material_accumulation"):
            for (t, x, p, j), v in self.material_accumulation.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.properties[t, x].get_material_flow_terms(p, j),
                        default=1,
                        warning=True,
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "_flow_terms"):
            for (t, x, p, j), v in self._flow_terms.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.properties[t, x].get_material_flow_terms(p, j),
                        default=1,
                        warning=True,
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "material_flow_dx"):
            for (t, x, p, j), v in self.material_flow_dx.items():
                if iscale.get_scaling_factor(v) is None:
                    # As domain is normalized, derivative should have same
                    # scale as flow
                    sf = iscale.get_scaling_factor(self._flow_terms[t, x, p, j])
                    iscale.set_scaling_factor(v, sf)

        # Control Volume has no way of knowing how best to scale
        # reaction extents - this is something only the unit model can provide
        # This also applies to the phase_equilibrium_generation term.

        if hasattr(self, "rate_reaction_generation"):
            for (t, x, p, j), v in self.rate_reaction_generation.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.min_scaling_factor(self.rate_reaction_extent[t, x, :])
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "equilibrium_reaction_generation"):
            for (t, x, p, j), v in self.equilibrium_reaction_generation.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.min_scaling_factor(
                        self.equilibrium_reaction_extent[t, x, ...]
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "inherent_reaction_generation"):
            for (t, x, p, j), v in self.inherent_reaction_generation.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.min_scaling_factor(
                        self.inherent_reaction_extent[t, x, ...]
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "mass_transfer_term"):
            for (t, x, p, j), v in self.mass_transfer_term.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.properties[t, x].get_material_flow_terms(p, j),
                        default=1,
                        warning=True,
                    )
                    iscale.set_scaling_factor(v, sf)

        # Set scaling factors for element balance variables
        if hasattr(self, "elemental_flow_term"):
            for (t, x, e), v in self.elemental_flow_term.items():
                flow_basis = self.properties[t, x].get_material_flow_basis()

                sf = iscale.min_scaling_factor(
                    [
                        self.properties[t, x].get_material_density_terms(p, j)
                        for (p, j) in phase_component_set
                    ],
                    default=1,
                    warning=True,
                )
                if flow_basis == MaterialFlowBasis.molar:
                    sf *= 1
                elif flow_basis == MaterialFlowBasis.mass:
                    # MW scaling factor is the inverse of its value
                    sf *= value(self.properties[t, x].mw_comp[j])

                iscale.set_scaling_factor(v, sf)

        if hasattr(self, "elemental_flow_dx"):
            for (t, x, e), v in self.elemental_flow_dx.items():
                if iscale.get_scaling_factor(v) is None:
                    # As domain is normalized, scale should be equal to flow
                    sf = iscale.get_scaling_factor(self.elemental_flow_term[t, x, e])
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "element_holdup"):
            for (t, x, e), v in self.element_holdup.items():
                flow_basis = self.properties[t, x].get_material_flow_basis()
                sf_list = []
                for p, j in phase_component_set:
                    if flow_basis == MaterialFlowBasis.molar:
                        sf = 1
                    elif flow_basis == MaterialFlowBasis.mass:
                        # MW scaling factor is the inverse of its value
                        sf = value(self.properties[t, x].mw_comp[j])
                    sf *= iscale.get_scaling_factor(self.phase_fraction[t, x, p])
                    sf *= iscale.get_scaling_factor(
                        self.properties[t, x].get_material_density_terms(p, j),
                        default=1,
                        warning=True,
                    )
                    sf *= value(self.properties[t, x].params.element_comp[j][e]) ** -1
                    sf_list.append(sf)
                sf_h = min(sf_list) * iscale.get_scaling_factor(self._area_func(t, x))
                iscale.set_scaling_factor(v, sf_h)

        if hasattr(self, "element_accumulation"):
            for (t, x, e), v in self.element_accumulation.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.min_scaling_factor(
                        self.elemental_flow_term[t, x, ...], default=1, warning=True
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "elemental_mass_transfer_term"):
            for (t, x, e), v in self.elemental_mass_transfer_term.items():
                # minimum scaling factor for elemental_flow terms
                sf_list = []
                flow_basis = self.properties[t, x].get_material_flow_basis()
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.min_scaling_factor(
                        self.elemental_flow_term[t, x, ...], default=1, warning=True
                    )
                    iscale.set_scaling_factor(v, sf)

        # Set scaling factors for enthalpy balance variables
        if hasattr(self, "_enthalpy_flow"):
            for (t, x, p), v in self._enthalpy_flow.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.properties[t, x].get_enthalpy_flow_terms(p),
                        default=1,
                        warning=True,
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "enthalpy_flow_dx"):
            for (t, x, p), v in self.enthalpy_flow_dx.items():
                if iscale.get_scaling_factor(v) is None:
                    # Normalized domain, so scale should be the same as flow
                    sf = iscale.get_scaling_factor(self._enthalpy_flow[t, x, p])
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "energy_holdup"):
            for (t, x, p), v in self.energy_holdup.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(self._area_func(t, x))
                    sf = iscale.get_scaling_factor(self.phase_fraction[t, x, p])
                    sf *= iscale.get_scaling_factor(
                        self.properties[t, x].get_energy_density_terms(p),
                        default=1,
                        warning=True,
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "energy_accumulation"):
            for (t, x, p), v in self.energy_accumulation.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.properties[t, x].get_enthalpy_flow_terms(p),
                        default=1,
                        warning=True,
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "heat"):
            for v in self.heat.values():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.heat, default=1e-6, warning=True
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "work"):
            for v in self.work.values():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.work, default=1e-6, warning=True
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "enthalpy_transfer"):
            for (t, x), v in self.enthalpy_transfer.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.min_scaling_factor(
                        [
                            self.properties[t, x].get_enthalpy_flow_terms(p)
                            for p in phase_list
                        ]
                    )
                    iscale.set_scaling_factor(v, sf)

        # Set scaling for momentum balance variables
        if hasattr(self, "pressure_dx"):
            for (t, x), v in self.pressure_dx.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = iscale.get_scaling_factor(
                        self.properties[t, x].pressure, default=1, warning=True
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "deltaP"):
            for (t, x), v in self.deltaP.items():
                if iscale.get_scaling_factor(v) is None:
                    sf = 10 * iscale.get_scaling_factor(
                        self.properties[t, x].pressure, default=1, warning=True
                    )
                    iscale.set_scaling_factor(v, sf)

        # Transform constraints in order of appearance
        if hasattr(self, "material_flow_linking_constraints"):
            for (t, x, p, j), c in self.material_flow_linking_constraints.items():
                sf = iscale.get_scaling_factor(self._flow_terms[t, x, p, j])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "material_holdup_calculation"):
            for (t, x, p, j), c in self.material_holdup_calculation.items():
                sf = iscale.get_scaling_factor(self.material_holdup[t, x, p, j])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "rate_reaction_stoichiometry_constraint"):
            for (t, x, p, j), c in self.rate_reaction_stoichiometry_constraint.items():
                sf = iscale.get_scaling_factor(
                    self.rate_reaction_generation[t, x, p, j]
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "equilibrium_reaction_stoichiometry_constraint"):
            for (
                t,
                x,
                p,
                j,
            ), c in self.equilibrium_reaction_stoichiometry_constraint.items():
                sf = iscale.get_scaling_factor(
                    self.equilibrium_reaction_generation[t, x, p, j]
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "material_balances"):
            mb_type = self._constructed_material_balance_type
            if mb_type == MaterialBalanceType.componentPhase:
                for (t, x, p, j), c in self.material_balances.items():
                    sf = iscale.get_scaling_factor(self._flow_terms[t, x, p, j])
                    iscale.constraint_scaling_transform(c, sf, overwrite=False)
            elif mb_type == MaterialBalanceType.componentTotal:
                for (t, x, j), c in self.material_balances.items():
                    sf = iscale.min_scaling_factor(
                        [
                            self._flow_terms[t, x, p, j]
                            for p in phase_list
                            if (p, j) in phase_component_set
                        ]
                    )
                    iscale.constraint_scaling_transform(c, sf, overwrite=False)
            else:
                _log.warning(f"Unknown material balance type {mb_type}")

        if hasattr(self, "elemental_flow_constraint"):
            for (t, x, e), c in self.elemental_flow_constraint.items():
                sf = iscale.get_scaling_factor(self.elemental_flow_term[t, x, e])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "element_balances"):
            for (t, x, e), c in self.element_balances.items():
                sf = iscale.get_scaling_factor(self.elemental_flow_dx[t, x, e])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "elemental_holdup_calculation"):
            for (t, x, e), c in self.elemental_holdup_calculation.items():
                sf = iscale.get_scaling_factor(self.element_holdup[t, x, e])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "enthalpy_flow_linking_constraint"):
            for (t, x, p), c in self.enthalpy_flow_linking_constraint.items():
                sf = iscale.get_scaling_factor(self._enthalpy_flow[t, x, p])
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "enthalpy_balances"):
            for (t, x), c in self.enthalpy_balances.items():
                sf = iscale.min_scaling_factor(
                    [self._enthalpy_flow[t, x, p] for p in phase_list]
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "energy_holdup_calculation"):
            for (t, x, p), c in self.energy_holdup_calculation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.energy_holdup[t, x, p]),
                    overwrite=False,
                )

        if hasattr(self, "pressure_balance"):
            for (t, x), c in self.pressure_balance.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.properties[t, x].pressure, default=1e-5
                    ),
                    overwrite=False,
                )

        if hasattr(self, "sum_of_phase_fractions"):
            for (t, x), c in self.sum_of_phase_fractions.items():
                sf = iscale.min_scaling_factor(
                    [self.phase_fraction[t, x, p] for p in phase_list]
                )
                iscale.constraint_scaling_transform(c, sf, overwrite=False)

        if hasattr(self, "material_flow_dx_disc_eq"):
            for (t, x, p, j), c in self.material_flow_dx_disc_eq.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.material_flow_dx[t, x, p, j]),
                    overwrite=False,
                )

        if hasattr(self, "material_accumulation_disc_eq"):
            for (t, x, p, j), c in self.material_accumulation_disc_eq.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.material_accumulation[t, x, p, j]),
                    overwrite=False,
                )

        # Scaling for discretization equations
        if hasattr(self, "enthalpy_flow_dx_disc_eq"):
            for (t, x, p), c in self.enthalpy_flow_dx_disc_eq.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.enthalpy_flow_dx[t, x, p]),
                    overwrite=False,
                )

        if hasattr(self, "energy_accumulation_disc_eq"):
            for (t, x, p), c in self.energy_accumulation_disc_eq.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.energy_accumulation[t, x, p]),
                    overwrite=False,
                )

        if hasattr(self, "pressure_dx_disc_eq"):
            for (t, x), c in self.pressure_dx_disc_eq.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.pressure_dx[t, x]),
                    overwrite=False,
                )

        if hasattr(self, "elemental_flow_dx_disc_eq"):
            for (t, x, e), c in self.elemental_flow_dx_disc_eq.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.elemental_flow_dx[t, x, e]),
                    overwrite=False,
                )

        if hasattr(self, "element_accumulation_disc_eq"):
            for (t, x, e), c in self.element_accumulation_disc_eq.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.element_accumulation[t, x, e]),
                    overwrite=False,
                )
