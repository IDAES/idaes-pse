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
Framework for generic reaction packages
"""
from math import log

# Import Pyomo libraries
from pyomo.environ import Block, Constraint, Expression, Set, units as pyunits, Var
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    ReactionBlockBase,
    MaterialFlowBasis,
)
from idaes.models.properties.modular_properties.base.utility import ConcentrationForm
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyPackageError,
)
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


class GenericReactionPackageError(PropertyPackageError):
    # Error message for when a property is called for but no option provided
    def __init__(self, block, prop):
        self.prop = prop
        self.block = block

    def __str__(self):
        return (
            f"Generic Reaction Package instance {self.block} called for "
            f"{self.prop}, but was not provided with a method "
            f"for this property. Please add a method for this property "
            f"in the reaction parameter configuration."
        )


rxn_config = ConfigBlock()
rxn_config.declare(
    "stoichiometry",
    ConfigValue(
        domain=dict,
        description="Stoichiometry of reaction",
        doc="Dict describing stoichiometry of reaction",
    ),
)
rxn_config.declare(
    "heat_of_reaction",
    ConfigValue(
        description="Method for calculating specific heat of reaction",
        doc="Valid Python class containing instructions on how to calculate "
        "the heat of reaction for this reaction.",
    ),
)
rxn_config.declare(
    "concentration_form",
    ConfigValue(
        default=None,
        domain=In(ConcentrationForm),
        description="Form to use for concentration terms in reaction equation",
        doc="ConcentrationForm Enum indicating what form to use for concentration "
        "terms when constructing reaction equation.",
    ),
)
rxn_config.declare(
    "parameter_data",
    ConfigValue(
        default={},
        domain=dict,
        description="Dict containing initialization data for parameters",
    ),
)

rate_rxn_config = rxn_config()
rate_rxn_config.declare(
    "rate_constant",
    ConfigValue(
        description="Expression form describing rate constant",
        doc="Valid Python class containing instructions on how to construct "
        "the rate constant for this reaction.",
    ),
)
rate_rxn_config.declare(
    "rate_form",
    ConfigValue(
        description="Expression form describing rate of reaction",
        doc="Valid Python class containing instructions on how to construct "
        "the rate expression for this reaction.",
    ),
)

equil_rxn_config = rxn_config()
equil_rxn_config.declare(
    "equilibrium_constant",
    ConfigValue(
        description="Expression form describing equilibrium constant",
        doc="Valid Python class containing instructions on how to construct "
        "the equilibrium constant for this reaction.",
    ),
)
equil_rxn_config.declare(
    "equilibrium_form",
    ConfigValue(
        description="Expression form describing reaction equilibrium",
        doc="Valid Python class containing instructions on how to construct "
        "the equilibrium constraint for this reaction.",
    ),
)


@declare_process_block_class("GenericReactionParameterBlock")
class GenericReactionParameterData(ReactionParameterBlock):
    """
    General Reaction Parameter Block Class
    """

    CONFIG = ReactionParameterBlock.CONFIG()

    CONFIG.declare(
        "reaction_basis",
        ConfigValue(
            default=MaterialFlowBasis.molar,
            domain=In(MaterialFlowBasis),
            doc="Basis of reactions",
            description="Argument indicating basis of reaction terms. Should be "
            "an instance of a MaterialFlowBasis Enum",
        ),
    )

    CONFIG.declare(
        "rate_reactions", ConfigBlock(implicit=True, implicit_domain=rate_rxn_config)
    )

    CONFIG.declare(
        "equilibrium_reactions",
        ConfigBlock(implicit=True, implicit_domain=equil_rxn_config),
    )

    # Base units of measurement
    CONFIG.declare(
        "base_units",
        ConfigValue(
            default={},
            domain=dict,
            description="Base units for property package",
            doc="Dict containing definition of base units of measurement to use "
            "with property package.",
        ),
    )

    # User-defined default scaling factors
    CONFIG.declare(
        "default_scaling_factors",
        ConfigValue(
            domain=dict,
            description="User-defined default scaling factors",
            doc="Dict of user-defined properties and associated default "
            "scaling factors",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        # Call super.build() to initialize Block
        # In this case we are replicating the super.build to get around a
        # chicken-and-egg problem
        # The super.build tries to validate units, but they have not been set
        # and cannot be set until the config block is created by super.build
        super(ReactionParameterBlock, self).build()
        self.default_scaling_factor = {}

        # Set base units of measurement
        self.get_metadata().add_default_units(self.config.base_units)

        # TODO: Need way to tie reaction package to a specfic property package
        self._validate_property_parameter_units()
        self._validate_property_parameter_properties()

        # Call configure method to set construction arguments
        self.configure()

        # Build core components
        self._reaction_block_class = GenericReactionBlock

        # Alias associated property package to keep line length down
        ppack = self.config.property_package

        if not hasattr(ppack, "_electrolyte") or not ppack._electrolyte:
            pc_set = ppack._phase_component_set
        elif ppack.config.state_components.name == "true":
            pc_set = ppack.true_phase_component_set
        elif ppack.config.state_components.name == "apparent":
            pc_set = ppack.apparent_phase_component_set
        else:
            raise BurntToast()

        # Construct rate reaction attributes if required
        if len(self.config.rate_reactions) > 0:
            # Construct rate reaction index
            self.rate_reaction_idx = Set(initialize=self.config.rate_reactions.keys())

            # Construct rate reaction stoichiometry dict
            self.rate_reaction_stoichiometry = {}
            for r, rxn in self.config.rate_reactions.items():
                for p, j in pc_set:
                    self.rate_reaction_stoichiometry[(r, p, j)] = 0

                if rxn.stoichiometry is None:
                    raise ConfigurationError(
                        "{} rate reaction {} was not provided with a "
                        "stoichiometry configuration argument.".format(self.name, r)
                    )
                else:
                    for k, v in rxn.stoichiometry.items():
                        if k[0] not in ppack.phase_list:
                            raise ConfigurationError(
                                "{} stoichiometry for rate reaction {} "
                                "included unrecognised phase {}.".format(
                                    self.name, r, k[0]
                                )
                            )
                        if k[1] not in ppack.component_list:
                            raise ConfigurationError(
                                "{} stoichiometry for rate reaction {} "
                                "included unrecognised component {}.".format(
                                    self.name, r, k[1]
                                )
                            )
                        self.rate_reaction_stoichiometry[(r, k[0], k[1])] = v

                # Check that a method was provided for the rate form
                if rxn.rate_form is None:
                    _log.debug(
                        "{} rate reaction {} was not provided with a "
                        "rate_form configuration argument. This is suitable "
                        "for processes using stoichiometric reactors, but not "
                        "for those using unit operations which rely on "
                        "reaction rate.".format(self.name, r)
                    )

        # Construct equilibrium reaction attributes if required
        if len(self.config.equilibrium_reactions) > 0:
            # Construct equilibrium reaction index
            self.equilibrium_reaction_idx = Set(
                initialize=self.config.equilibrium_reactions.keys()
            )

            # Construct equilibrium reaction stoichiometry dict
            self.equilibrium_reaction_stoichiometry = {}
            for r, rxn in self.config.equilibrium_reactions.items():
                for p, j in pc_set:
                    self.equilibrium_reaction_stoichiometry[(r, p, j)] = 0

                if rxn.stoichiometry is None:
                    raise ConfigurationError(
                        "{} equilibrium reaction {} was not provided with a "
                        "stoichiometry configuration argument.".format(self.name, r)
                    )
                else:
                    for k, v in rxn.stoichiometry.items():
                        if k[0] not in ppack.phase_list:
                            raise ConfigurationError(
                                "{} stoichiometry for equilibrium reaction {} "
                                "included unrecognised phase {}.".format(
                                    self.name, r, k[0]
                                )
                            )
                        if k[1] not in ppack.component_list:
                            raise ConfigurationError(
                                "{} stoichiometry for equilibrium reaction {} "
                                "included unrecognised component {}.".format(
                                    self.name, r, k[1]
                                )
                            )
                        self.equilibrium_reaction_stoichiometry[(r, k[0], k[1])] = v

                # Check that a method was provided for the equilibrium form
                if rxn.equilibrium_form is None:
                    raise ConfigurationError(
                        "{} equilibrium reaction {} was not provided with a "
                        "equilibrium_form configuration argument.".format(self.name, r)
                    )

        # Add a master reaction index which includes both types of reactions
        if (
            len(self.config.rate_reactions) > 0
            and len(self.config.equilibrium_reactions) > 0
        ):
            self.reaction_idx = Set(
                initialize=(self.rate_reaction_idx | self.equilibrium_reaction_idx)
            )
        elif len(self.config.rate_reactions) > 0:
            self.reaction_idx = Set(initialize=self.rate_reaction_idx)
        elif len(self.config.equilibrium_reactions) > 0:
            self.reaction_idx = Set(initialize=self.equilibrium_reaction_idx)
        else:
            raise BurntToast(
                "{} Generic property package failed to construct "
                "master reaction Set. This should not happen. "
                "Please contact the IDAES developers with this "
                "bug".format(self.name)
            )

        # Construct blocks to contain parameters for each reaction
        for r in self.reaction_idx:
            self.add_component("reaction_" + str(r), Block())

        # Build parameters
        if len(self.config.rate_reactions) > 0:
            for r in self.rate_reaction_idx:
                rblock = getattr(self, "reaction_" + r)
                r_config = self.config.rate_reactions[r]

                order_init = {}
                for p, j in pc_set:
                    if "reaction_order" in r_config.parameter_data:
                        try:
                            order_init[p, j] = r_config.parameter_data[
                                "reaction_order"
                            ][p, j]
                        except KeyError:
                            order_init[p, j] = 0
                    else:
                        # Assume elementary reaction and use stoichiometry
                        try:
                            if r_config.stoichiometry[p, j] < 0:
                                # These are reactants, but order is -ve stoic
                                order_init[p, j] = -r_config.stoichiometry[p, j]
                            else:
                                # Anything else is a product, not be included
                                order_init[p, j] = 0
                        except KeyError:
                            order_init[p, j] = 0

                rblock.reaction_order = Var(
                    pc_set, initialize=order_init, doc="Reaction order", units=None
                )

                for val in self.config.rate_reactions[r].values():
                    try:
                        val.build_parameters(rblock, self.config.rate_reactions[r])
                    except AttributeError:
                        pass

        if len(self.config.equilibrium_reactions) > 0:
            for r in self.equilibrium_reaction_idx:
                rblock = getattr(self, "reaction_" + r)
                r_config = self.config.equilibrium_reactions[r]

                order_init = {}
                for p, j in pc_set:
                    if "reaction_order" in r_config.parameter_data:
                        try:
                            order_init[p, j] = r_config.parameter_data[
                                "reaction_order"
                            ][p, j]
                        except KeyError:
                            order_init[p, j] = 0
                    else:
                        # Assume elementary reaction and use stoichiometry
                        try:
                            # Here we use the stoic. coeff. directly
                            # However, solids should be excluded as they
                            # normally do not appear in the equilibrium
                            # relationship
                            pobj = ppack.get_phase(p)
                            if not pobj.is_solid_phase():
                                order_init[p, j] = r_config.stoichiometry[p, j]
                            else:
                                order_init[p, j] = 0
                        except KeyError:
                            order_init[p, j] = 0

                rblock.reaction_order = Var(
                    pc_set, initialize=order_init, doc="Reaction order", units=None
                )

                for val in self.config.equilibrium_reactions[r].values():
                    try:
                        val.build_parameters(
                            rblock, self.config.equilibrium_reactions[r]
                        )
                    except AttributeError:
                        pass
                    except KeyError as err:
                        # This likely arises from mismatched true and apparent
                        # species sets. Reaction packages must use the same
                        # basis as the associated thermo properties
                        # Raise an exception to inform the user
                        raise PropertyPackageError(
                            "{} KeyError encountered whilst constructing "
                            "reaction parameters. This may be due to "
                            "mismatched state_components between the "
                            "Reaction Package and the associated Physical "
                            "Property Package - Reaction Packages must use the"
                            "same basis (true or apparent species) as the "
                            "Physical Property Package.".format(self.name),
                            err,
                        )

        # As a safety check, make sure all Vars in reaction blocks are fixed
        for v in self.component_objects(Var, descend_into=True):
            for i in v:
                if v[i].value is None:
                    raise ConfigurationError(
                        "{} parameter {} was not assigned"
                        " a value. Please check your configuration "
                        "arguments.".format(self.name, v.local_name)
                    )
                v[i].fix()

        # Set default scaling factors
        if self.config.default_scaling_factors is not None:
            self.default_scaling_factor.update(self.config.default_scaling_factors)
        # Finally, call populate_default_scaling_factors method to fill blanks
        iscale.populate_default_scaling_factors(self)

    def configure(self):
        """
        Placeholder method to allow users to specify config arguments via a
        class. The user class should inherit from this one and implement a
        configure() method which sets the values of the desired config
        arguments.

        Args:
            None

        Returns:
            None
        """
        pass

    def parameters(self):
        """
        Placeholder method to allow users to specify parameters via a
        class. The user class should inherit from this one and implement a
        parameters() method which creates the required components.

        Args:
            None

        Returns:
            None
        """
        pass

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {
                "dh_rxn": {"method": "_dh_rxn"},
                "k_eq": {"method": "_k_eq"},
                "log_k_eq": {"method": "_log_k_eq"},
                "k_rxn": {"method": "_k_rxn"},
                "reaction_rate": {"method": "_reaction_rate"},
            }
        )


class _GenericReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(blk, outlvl=idaeslog.NOTSET, **kwargs):
        """
        Initialization routine for reaction package.

        Keyword Arguments:
            outlvl : sets output level of initialization routine

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        init_log.info("Initialization Complete.")


@declare_process_block_class("GenericReactionBlock", block_class=_GenericReactionBlock)
class GenericReactionBlockData(ReactionBlockDataBase):
    def build(self):
        # TODO: Need a different error here
        super(GenericReactionBlockData, self).build()

        if self.config.has_equilibrium:
            if len(self.params.config.equilibrium_reactions) == 0:
                raise PropertyPackageError(
                    "{} Generic Reaction Block was set to include equilibrium "
                    "reactions, however no equilibrium reactions were "
                    "defined. Either set has_equilibrium to be False, or "
                    "include equilibrium reactions in the package definition.".format(
                        self.name
                    )
                )
            self._equilibrium_constraint()

    def calculate_scaling_factors(self):
        # Get default scale factors and do calculations from base classes
        super().calculate_scaling_factors()

        # Lock attribute creation to avoid hasattr triggering builds
        self._lock_attribute_creation = True

        # Due to the exponential relationship between most reaction properties
        # and temeprature, it is very hard to calculate good scaling factors
        # from order-of-magnitude guesses. Thus ,reaction scaling will always
        # require a lot of user input. Here we will calculate scaling factors
        # for those properties that have non-exponential relationships to T.
        if hasattr(self, "dh_rxn"):
            for r, v in self.dh_rxn.items():
                if iscale.get_scaling_factor(v) is None:
                    rblock = getattr(self.params, "reaction_" + r)
                    try:
                        carg = self.params.config.rate_reactions[r]
                    except (AttributeError, KeyError):
                        carg = self.params.config.equilibrium_reactions[r]
                    sf = carg["heat_of_reaction"].calculate_scaling_factors(
                        self, rblock
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "k_eq"):
            for r, v in self.k_eq.items():
                if iscale.get_scaling_factor(v) is None:
                    rblock = getattr(self.params, "reaction_" + r)
                    carg = self.params.config.equilibrium_reactions[r]
                    sf = carg["equilibrium_constant"].calculate_scaling_factors(
                        self, rblock
                    )
                    iscale.set_scaling_factor(v, sf)

        if hasattr(self, "log_k_eq"):
            for r, v in self.log_k_eq.items():
                if iscale.get_scaling_factor(v) is None:
                    rblock = getattr(self.params, "reaction_" + r)
                    carg = self.params.config.equilibrium_reactions[r]
                    sf = carg["equilibrium_constant"].calculate_scaling_factors(
                        self, rblock
                    )
                    iscale.set_scaling_factor(v, log(sf))

        if hasattr(self, "equilibrium_constraint"):
            for r, v in self.equilibrium_constraint.items():
                carg = self.params.config.equilibrium_reactions[r]
                if iscale.get_scaling_factor(v) is None:
                    if carg["equilibrium_form"].__name__.startswith("log_"):
                        # Log form constraint
                        sf = iscale.get_scaling_factor(
                            self.log_k_eq[r], default=1, warning=True
                        )
                    else:
                        sf = iscale.get_scaling_factor(
                            self.k_eq[r], default=1, warning=True
                        )

                sf_const = carg["equilibrium_form"].calculate_scaling_factors(self, sf)

                iscale.constraint_scaling_transform(v, sf_const, overwrite=False)

        # Unlock attribute creation when done
        self._lock_attribute_creation = False

    def _dh_rxn(self):
        def dh_rule(b, r):
            rblock = getattr(b.params, "reaction_" + r)
            try:
                carg = b.params.config.rate_reactions[r]
            except (AttributeError, KeyError):
                carg = b.params.config.equilibrium_reactions[r]
            return carg["heat_of_reaction"].return_expression(
                b, rblock, r, b.state_ref.temperature
            )

        self.dh_rxn = Expression(
            self.params.reaction_idx, doc="Specific heat of reaction", rule=dh_rule
        )

    def _k_rxn(self):
        def krxn_rule(b, r):
            rblock = getattr(b.params, "reaction_" + r)

            carg = b.params.config.rate_reactions[r]

            return carg["rate_constant"].return_expression(
                b, rblock, r, b.state_ref.temperature
            )

        self.k_rxn = Expression(
            self.params.rate_reaction_idx, doc="Reaction rate constant", rule=krxn_rule
        )

    def _reaction_rate(self):
        def rate_rule(b, r):
            rblock = getattr(b.params, "reaction_" + r)

            carg = b.params.config.rate_reactions[r]

            if carg["rate_form"] is None:
                raise ConfigurationError(
                    f"{b.name} Generic Reaction {r} was not provided with a "
                    f"rate_form configuration argument."
                )

            return carg["rate_form"].return_expression(
                b, rblock, r, b.state_ref.temperature
            )

        self.reaction_rate = Expression(
            self.params.rate_reaction_idx, doc="Reaction rate", rule=rate_rule
        )

    def _k_eq(self):
        def keq_rule(b, r):
            rblock = getattr(b.params, "reaction_" + r)

            carg = b.params.config.equilibrium_reactions[r]

            return carg["equilibrium_constant"].return_expression(
                b, rblock, r, b.state_ref.temperature
            )

        self.k_eq = Expression(
            self.params.equilibrium_reaction_idx,
            doc="Equilibrium constant",
            rule=keq_rule,
        )

    def _log_k_eq(self):
        def log_keq_rule(b, r):
            rblock = getattr(b.params, "reaction_" + r)

            carg = b.params.config.equilibrium_reactions[r]

            return carg["equilibrium_constant"].return_log_expression(
                b, rblock, r, b.state_ref.temperature
            )

        self.log_k_eq = Var(
            self.params.equilibrium_reaction_idx,
            initialize=1,
            doc="Log of equilibrium constant",
            units=pyunits.dimensionless,
        )

        self.log_k_eq_constraint = Constraint(
            self.params.equilibrium_reaction_idx,
            rule=log_keq_rule,
            doc="Constraint for log of K_eq",
        )

    def _equilibrium_constraint(self):
        def equil_rule(b, r):
            rblock = getattr(b.params, "reaction_" + r)

            carg = b.params.config.equilibrium_reactions[r]

            return carg["equilibrium_form"].return_expression(
                b, rblock, r, b.state_ref.temperature
            )

        self.equilibrium_constraint = Constraint(
            self.params.equilibrium_reaction_idx,
            doc="Equilibrium constraint",
            rule=equil_rule,
        )

    def get_reaction_rate_basis(b):
        return b.params.config.reaction_basis
