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
Tool to interrogate IDAES flowsheets and list the reaction properties
required to simulate it.
"""
import sys

# Import Pyomo libraries
from pyomo.environ import Set, Var, units as pyunits

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    ReactionParameterBlock,
    ReactionBlockDataBase,
    ReactionBlockBase,
    UnitModelBlockData,
)
import idaes.logger as idaeslog

# Some more information about this module
__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ReactionInterrogatorBlock")
class ReactionInterrogatorData(ReactionParameterBlock):
    """
    Interrogator Parameter Block Class

    This class contains the methods and attributes for recording and displaying
    the reaction properties requried by the flowsheet.
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super(ReactionInterrogatorData, self).build()

        self._reaction_block_class = InterrogatorReactionBlock

        # List of valid phases in property package
        p_list = []
        for p in self.config.property_package.phase_list:
            p_list.append(p)
        self.phase_list = Set(initialize=p_list)

        # Component list - a list of component identifiers
        c_list = []
        for j in self.config.property_package.component_list:
            c_list.append(j)
        self.component_list = Set(initialize=c_list)

        # Set up dict to record property calls
        self.required_properties = {}

        # Dummy reaction definition
        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1"])

        # Reaction Stoichiometry - fake the dict with all 1's
        self.rate_reaction_stoichiometry = {}
        for p in p_list:
            for j in c_list:
                self.rate_reaction_stoichiometry[("R1", p, j)] = 1

    def list_required_properties(self):
        """
        Method to list all reaction properties required by the flowsheet.

        Args:
            None

        Returns:
            A list of properties required
        """
        return list(self.required_properties)

    def list_models_requiring_property(self, prop):
        """
        Method to list all models in the flowsheet requiring the given
        property.

        Args:
            prop : the property of interest

        Returns:
            A list of unit model names which require prop
        """
        try:
            return self.required_properties[prop]
        except KeyError:
            raise KeyError(
                "Property {} does not appear in required_properties. "
                "Please check the spelling of the property that you are "
                "interested in.".format(prop)
            )

    def list_properties_required_by_model(self, model):
        """
        Method to list all reaction properties required by a given unit model.

        Args:
            model : the unit model of interest. Can be given as either a model
                    component or the unit name as a string

        Returns:
            A list of reaction properties required by model
        """
        prop_list = []
        if not isinstance(model, str):
            model = model.name

        for k, v in self.required_properties.items():
            if model in v:
                prop_list.append(k)

        if len(prop_list) < 1:
            raise ValueError(
                "Model {} does not appear in the flowsheet. Please check "
                "the spelling of the model provided."
            )
        else:
            return prop_list

    def print_required_properties(self, ostream=None):
        """
        Method to print a summary of the reaction properties required by the
        flowsheet.

        Args:
            ostream : output stream to print to. If not provided will print to
                      sys.stdout

        Returns:
            None
        """
        if ostream is None:
            ostream = sys.stdout

        # Write header
        max_str_length = 74
        tab = " " * 4
        ostream.write("\n" + "=" * max_str_length + "\n")
        ostream.write("Reaction Property Interrogator Summary" + "\n")
        ostream.write(
            "\n"
            + "The Flowsheet requires the following reaction properties "
            + "(times required):"
            + "\n"
            + "\n"
        )
        for k, v in self.required_properties.items():
            lead_str = tab + k
            trail_str = str(len(v))
            mid_str = " " * (max_str_length - len(lead_str) - len(trail_str))
            ostream.write(lead_str + mid_str + trail_str + "\n")
        ostream.write(
            "\n"
            + "Note: User constraints may require additional properties "
            + "which are not"
            + "\n"
            + "reported here."
            + "\n"
        )

    def print_models_requiring_property(self, prop, ostream=None):
        """
        Method to print a summary of the models in the flowsheet requiring a
        given property.

        Args:
            prop : the property of interest.
            ostream : output stream to print to. If not provided will print to
                      sys.stdout

        Returns:
            None
        """
        if ostream is None:
            ostream = sys.stdout

        tab = " " * 4

        ostream.write("\n")
        ostream.write(
            f"The following models in the Flowsheet " f"require {prop}:" + "\n"
        )

        for m in self.required_properties[prop]:
            ostream.write(tab + m + "\n")

    def print_properties_required_by_model(self, model, ostream=None):
        """
        Method to print a summary of the reaction properties required by
        a given unit model.

        Args:
            model : the unit model of interest.
            ostream : output stream to print to. If not provided will print to
                      sys.stdout

        Returns:
            None
        """
        if not isinstance(model, str):
            model = model.name

        if ostream is None:
            ostream = sys.stdout

        tab = " " * 4

        ostream.write("\n")
        ostream.write(
            f"The following reaction properties are required by "
            f"model {model}:" + "\n"
        )

        for m in self.list_properties_required_by_model(model):
            ostream.write(tab + m + "\n")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


class _InterrogatorReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """

    def initialize(blk, *args, **kwargs):
        """
        Dummy initialization routine, This will raise an TypeError if a user
        tries to initialize a model using the Interrogator Reaction Package
        and tell them that the model cannot be solved.
        """
        raise TypeError(
            "Models constructed using the Reaction Interrogator package "
            "cannot be used to solve a flowsheet. Please rebuild your "
            "flowsheet using a valid reaction package."
        )


@declare_process_block_class(
    "InterrogatorReactionBlock", block_class=_InterrogatorReactionBlock
)
class InterrogatorReactionBlockData(ReactionBlockDataBase):
    """
    A dummy reaction block for interrogating flowsheets and recording reaction
    properties called for during construction.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(InterrogatorReactionBlockData, self).build()

        # Add dummy vars for returning expressions
        self._dummy_var = Var(initialize=1)
        self._dummy_var_phase = Var(self.params.phase_list, initialize=1)
        self._dummy_var_comp = Var(self.params.component_list, initialize=1)
        self._dummy_var_phase_comp = Var(
            self.params.phase_list, self.params.component_list, initialize=1
        )
        self._dummy_reaction_idx = Var(self.params.rate_reaction_idx, initialize=1)

    # Set default values for required attributes so construction doesn't fail
    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar

    def __getattr__(self, prop):
        """
        Overload getattr to log each call for an unknown attribute, assuming
        these are all properties.
        Then, return a dummy variable with the correct indexing set.
        """
        # Log call
        self._log_call(prop)

        # Return dummy var
        if prop in ["reaction_rate", "dh_rxn"]:
            return self._dummy_reaction_idx
        elif prop.endswith("_phase_comp"):
            return self._dummy_var_phase_comp
        elif prop.endswith("_phase"):
            return self._dummy_var_phase
        elif prop.endswith("_comp"):
            return self._dummy_var_comp
        else:
            return self._dummy_var

    def _log_call(self, prop):
        """
        Method to log calls for properties in required_properties dict
        """
        # Get the required_properties dict from parameter block
        prop_dict = self.params.required_properties

        # Get name of parent unit to record in required_properties
        name = self._get_parent_unit_name()

        try:
            # If name is not listed for current property, add to list
            if name not in prop_dict[prop]:
                prop_dict[prop].append(name)
        except KeyError:
            # If a KeyError occurs, it means property has not been logged
            # before, so add new entry to dict
            prop_dict[prop] = [name]

    def _get_parent_unit_name(self):
        """
        Method to find the parent unit of the current StateBlock (if one
        exists) and return this so it can be logged as in required_properties.

        If current StateBlock has no parent unit, it is a stand-alone
        StateBlock, so log the name of this instead.
        """
        # Start with current block (i.e. a StateBlock)
        parent = self

        # Search up the parent tree until we find a UnitModel or top of tree
        while True:
            if isinstance(parent, UnitModelBlockData):
                # If parent is a UnitModel, we have found our target
                # Return parent name
                return parent.name
            else:
                if parent.parent_block() is None:
                    # Check if the parent object has no parent, i.e. is top of
                    # tree. If so, we are dealling with a stand-alone
                    # StateBlock.
                    # Return name of parent_component to strip indices
                    return self.parent_component().name
                else:
                    # Otherwise continue searching up tree
                    parent = parent.parent_block()
