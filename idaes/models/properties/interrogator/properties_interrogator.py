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
Tool to interrogate IDAES flowsheets and list the physical properties
required to simulate it.
"""
import sys
from inspect import isclass

# Import Pyomo libraries
from pyomo.environ import Set, Var, units as pyunits
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    UnitModelBlockData,
    Phase,
    LiquidPhase,
    VaporPhase,
    Component,
)
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

# Some more information about this module
__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PropertyInterrogatorBlock")
class PropertyInterrogatorData(PhysicalParameterBlock):
    """
    Interrogator Parameter Block Class

    This class contains the methods and attributes for recording and displaying
    the properties requried by the flowsheet.
    """

    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare(
        "phase_list",
        ConfigValue(
            domain=dict,
            description="User defined phase list. Dict with form {name: Type}",
        ),
    )
    CONFIG.declare(
        "component_list",
        ConfigValue(
            domain=dict,
            description="User defined component list. Dict with form {name: Type}",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        super(PropertyInterrogatorData, self).build()

        self._state_block_class = InterrogatorStateBlock

        # Phase objects
        if self.config.phase_list is None:
            self.Liq = LiquidPhase()
            self.Vap = VaporPhase()
        else:
            for p, t in self.config.phase_list.items():
                if t is None:
                    t = Phase
                elif not isclass(t) or not issubclass(t, Phase):
                    raise ConfigurationError(
                        f"{self.name} invalid phase type {t} (for phase {p})."
                        f" Type must be a subclass of Phase."
                    )
                self.add_component(p, t())

        # Component objects
        if self.config.component_list is None:
            self.A = Component()
            self.B = Component()
        else:
            for j, t in self.config.component_list.items():
                if t is None:
                    t = Component
                elif not isclass(t) or not issubclass(t, Component):
                    raise ConfigurationError(
                        f"{self.name} invalid component type {t} (for "
                        f"component {j}). Type must be a subclass of "
                        f"Component."
                    )
                self.add_component(j, t())

        # Set up dict to record property calls
        self.required_properties = {}

        # Dummy phase equilibrium definition so we can handle flash cases
        self.phase_equilibrium_idx = Set(initialize=[1])

        self.phase_equilibrium_list = {1: ["A", ("Vap", "Liq")]}

    def list_required_properties(self):
        """
        Method to list all thermophysical properties required by the flowsheet.

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
        Method to list all thermophysical properties required by a given unit
        model.

        Args:
            model : the unit model of interest. Can be given as either a model
                    component or the unit name as a string

        Returns:
            A list of thermophysical properties required by model
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
        Method to print a summary of the thermophysical properties required by
        the flowsheet.

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
        ostream.write("Property Interrogator Summary" + "\n")
        ostream.write(
            "\n"
            + "The Flowsheet requires the following properties "
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
        Method to print a summary of the thermophysical properties required by
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
            f"The following properties are required by model " f"{model}:" + "\n"
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


class _InterrogatorStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(blk, *args, **kwargs):
        """
        Dummy initialization routine, This will raise an TypeError if a user
        tries to initialize a model using the Interrogator Property Package
        and tell them that the model cannot be solved.
        """
        raise TypeError(
            "Models constructed using the Property Interrogator package "
            "cannot be used to solve a flowsheet. Please rebuild your "
            "flowsheet using a valid property package."
        )


@declare_process_block_class(
    "InterrogatorStateBlock", block_class=_InterrogatorStateBlock
)
class InterrogatorStateBlockData(StateBlockData):
    """
    A dummy state block for interrogating flowsheets and recording physical
    properties called for during construction.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(InterrogatorStateBlockData, self).build()

        # Add dummy vars for building Ports and returning expressions
        self._dummy_var = Var(initialize=1)
        self._dummy_var_phase = Var(self.params.phase_list, initialize=1)
        self._dummy_var_comp = Var(self.params.component_list, initialize=1)
        self._dummy_var_phase_comp = Var(
            self.params.phase_list, self.params.component_list, initialize=1
        )
        # T and P are often involved in unit conversion checks, so need ot have units
        self._dummy_var_T = Var(initialize=1, units=pyunits.K)
        self._dummy_var_P = Var(initialize=1, units=pyunits.Pa)

    # Define standard methods and log calls before returning dummy variable
    def get_material_flow_terms(self, p, j):
        self._log_call("material flow terms")
        return self._dummy_var

    def get_enthalpy_flow_terms(self, p):
        self._log_call("enthalpy flow terms")
        return self._dummy_var

    def get_material_density_terms(self, p, j):
        self._log_call("material density terms")
        return self._dummy_var

    def get_energy_density_terms(self, p):
        self._log_call("energy density terms")
        return self._dummy_var

    # Set default values for required attributes so construction doesn't fail
    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(b):
        return {"_dummy_var": b._dummy_var}

    def define_display_vars(b):
        raise TypeError(
            "Models constructed using the Property Interrogator package "
            "should not be used for report methods."
        )

    def get_material_flow_basis(b):
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
        if prop.endswith("_phase_comp"):
            return self._dummy_var_phase_comp
        elif prop.endswith("_phase"):
            return self._dummy_var_phase
        elif prop.endswith("_comp"):
            return self._dummy_var_comp
        elif prop == "temperature":
            # Need this for unit conversion checks in some models
            return self._dummy_var_T
        elif prop == "pressure":
            # Need this for unit conversion checks in some models
            return self._dummy_var_P
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
