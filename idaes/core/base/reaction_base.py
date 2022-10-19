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
This module contains classes for reaction blocks and reaction parameter blocks.
"""

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, Bool
from pyomo.environ import Var, Constraint, Expression

# Import IDAES cores
from idaes.core.base.process_block import ProcessBlock
from idaes.core import ProcessBlockData, MaterialFlowBasis
from idaes.core.base import property_meta
from idaes.core.util.exceptions import (
    BurntToast,
    PropertyNotSupportedError,
    PropertyPackageError,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
    is_state_block,
)
from idaes.core.util.misc import add_object_reference

# WHY on Python 3.6, using the alternate syntax "import idaes.core.util.scaling as iscale"
# fails with "AttributeError: module 'idaes' has no attribute 'core'"
# this is likely due to a bug/limitation in how the Python import mechanism resolves circular imports
# for more information, see https://stackoverflow.com/questions/24807434
# and the official Python bug report: http://bugs.python.org/issue30024
from idaes.core.util import scaling as iscale
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

# Some more information about this module
__author__ = "Andrew Lee, John Eslick"

__all__ = [
    "ReactionBlock",  # pylint: disable=undefined-all-variable
    "ReactionParameterBlock",  # pylint: disable=undefined-all-variable
]


class _lock_attribute_creation_context(object):
    """Context manager to lock creation of new attributes on a state block"""

    def __init__(self, block):
        self.block = block

    def __enter__(self):
        self.block._lock_attribute_creation = True

    def __exit__(self, exc_type, exc_value, traceback):
        self.block._lock_attribute_creation = False


class ReactionParameterBlock(ProcessBlockData, property_meta.HasPropertyClassMetadata):
    """
    This is the base class for reaction parameter blocks. These are blocks
    that contain a set of parameters associated with a specific reaction
    package, and are linked to by all instances of that reaction package.
    """

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "property_package",
        ConfigValue(
            description="Reference to associated PropertyPackageParameter " "object",
            domain=is_physical_parameter_block,
        ),
    )
    CONFIG.declare(
        "default_arguments",
        ConfigBlock(
            description="Default arguments to use with Property Package", implicit=True
        ),
    )

    def __init__(self, *args, **kwargs):
        self.__reaction_block_class = None
        super().__init__(*args, **kwargs)

    def build(self):
        """
        General build method for ReactionParameterBlocks. Inheriting models
        should call super().build.

        Args:
            None

        Returns:
            None
        """
        super(ReactionParameterBlock, self).build()

        if not hasattr(self, "_reaction_block_class"):
            self._reaction_block_class = None

        # TODO: Need way to tie reaction package to a specfic property package
        self._validate_property_parameter_units()
        self._validate_property_parameter_properties()

        # This is a dict to store default property scaling factors. They are
        # defined in the parameter block to provide a universal default for
        # quantities in a particular kind of state block.  For example, you can
        # set flow scaling once instead of for every state block. Some of these
        # may be left for the user to set and some may be defined in a property
        # module where reasonable defaults can be defined a priori. See
        # set_default_scaling, get_default_scaling, and unset_default_scaling
        self.default_scaling_factor = {}

    def set_default_scaling(self, attrbute, value, index=None):
        """Set a default scaling factor for a property.

        Args:
            attribute: property attribute name
            value: default scaling factor
            index: for indexed properties, if this is not provied the scaling
                factor default applies to all indexed elements where specific
                indexes are no specifcally specified.

        Returns:
            None
        """
        self.default_scaling_factor[(attrbute, index)] = value

    def unset_default_scaling(self, attrbute, index=None):
        """Remove a previously set default value

        Args:
            attribute: property attribute name
            index: optional index for indexed properties

        Returns:
            None
        """
        try:
            del self.default_scaling_factor[(attrbute, index)]
        except KeyError:
            pass

    def get_default_scaling(self, attrbute, index=None):
        """Returns a default scale factor for a property

        Args:
            attribute: property attribute name
            index: optional index for indexed properties

        Returns:
            None
        """
        try:
            # If a specific component data index exists
            return self.default_scaling_factor[(attrbute, index)]
        except KeyError:
            try:
                # indexed, but no specifc index?
                return self.default_scaling_factor[(attrbute, None)]
            except KeyError:
                # Can't find a default scale factor for what you asked for
                return None

    @property
    def reaction_block_class(self):
        if self._reaction_block_class is not None:
            return self._reaction_block_class
        else:
            raise AttributeError(
                "{} has not assigned a ReactionBlock class to be associated "
                "with this reaction package. Please contact the developer of "
                "the reaction package.".format(self.name)
            )

    def build_reaction_block(self, *args, **kwargs):
        """
        Methods to construct a ReactionBlock assoicated with this
        ReactionParameterBlock. This will automatically set the parameters
        construction argument for the ReactionBlock.

        Returns:
            ReactionBlock

        """
        default = kwargs.pop("default", {})
        initialize = kwargs.pop("initialize", {})

        if initialize == {}:
            default["parameters"] = self
        else:
            for i in initialize.keys():
                initialize[i]["parameters"] = self

        return self.reaction_block_class(  # pylint: disable=not-callable
            *args, **kwargs, **default, initialize=initialize
        )

    def _validate_property_parameter_units(self):
        """
        Checks that the property parameter block associated with the
        reaction block uses the same set of default units.
        """
        r_units = self.get_metadata().default_units
        prop_units = self.config.property_package.get_metadata().default_units
        if not r_units.unitset_is_consistent(prop_units):
            raise PropertyPackageError(
                "{} the property package associated with this "
                "reaction package does not use the same set of "
                "units of measurement. Please choose a "
                "property package which uses the same units.".format(self.name)
            )

    def _validate_property_parameter_properties(self):
        """
        Checks that the property parameter block associated with the
        reaction block supports the necessary properties with correct units.
        """
        req_props = self.get_metadata().required_properties
        supp_props = self.config.property_package.get_metadata().properties

        for p in req_props:
            if p not in supp_props:
                raise PropertyPackageError(
                    "{} the property package associated with this "
                    "reaction package does not support the necessary "
                    "property, {}. Please choose a property package "
                    "which supports all required properties.".format(self.name, p)
                )
            elif supp_props[p]["method"] is False:
                raise PropertyPackageError(
                    "{} the property package associated with this "
                    "reaction package does not support the necessary "
                    "property, {}. Please choose a property package "
                    "which supports all required properties.".format(self.name, p)
                )

            # Check property units
            if req_props[p]["units"] != supp_props[p]["units"]:
                raise PropertyPackageError(
                    "{} the units associated with property {} in this "
                    "reaction package ({}) do not match with the units "
                    "used in the assoicated property package ({}). Please "
                    "choose a property package which used the same "
                    "units for all properties.".format(
                        self.name, p, req_props[p]["units"], supp_props[p]["units"]
                    )
                )


class ReactionBlockBase(ProcessBlock):
    """
    This is the base class for reaction block objects. These are used when
    constructing the SimpleBlock or IndexedBlock which will contain the
    PropertyData objects, and contains methods that can be applied to
    multiple ReactionBlockData objects simultaneously.
    """

    def initialize(self, *args):
        """
        This is a default initialization routine for ReactionBlocks to ensure
        that a routine is present. All ReactionBlockData classes should
        overload this method with one suited to the particular reaction package

        Args:
            None

        Returns:
            None
        """
        raise NotImplementedError(
            "{} reaction package has not implemented the"
            " initialize method. Please contact "
            "the reaction package developer".format(self.name)
        )

    def report(self, index=(0), true_state=False, dof=False, ostream=None, prefix=""):
        raise NotImplementedError(
            """The current Reaction Package has not implemented a report
                method. Please contact the package developer about this."""
        )


class ReactionBlockDataBase(ProcessBlockData):
    """
    This is the base class for reaction block data objects. These are
    blocks that contain the Pyomo components associated with calculating a
    set of reacion properties for a given material.
    """

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "parameters",
        ConfigValue(
            domain=is_reaction_parameter_block,
            description="""A reference to an instance of the Reaction Parameter
Block associated with this property package.""",
        ),
    )
    CONFIG.declare(
        "state_block",
        ConfigValue(
            domain=is_state_block,
            description="""A reference to an instance of a StateBlock with
which this reaction block should be associated.""",
        ),
    )
    CONFIG.declare(
        "has_equilibrium",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Equilibrium constraint flag",
            doc="""Flag indicating whether equilibrium constraints
should be constructed in this reaction block,
**default** - True.
**Valid values:** {
**True** - ReactionBlock should enforce equilibrium constraints,
**False** - ReactionBlock should not enforce equilibrium constraints.}""",
        ),
    )

    def __init__(self, *args, **kwargs):
        self._lock_attribute_creation = False
        super().__init__(*args, **kwargs)

    @property
    def component_list(self):
        return self.state_ref.component_list

    @property
    def phase_list(self):
        return self.state_ref.phase_list

    @property
    def phase_component_set(self):
        return self.state_ref.phase_component_set

    def lock_attribute_creation_context(self):
        """Returns a context manager that does not allow attributes to be created
        while in the context and allows attributes to be created normally outside
        the context.
        """
        return _lock_attribute_creation_context(self)

    def is_property_constructed(self, attr):
        """Returns True if the attribute ``attr`` already exists, or false if it
        would be added in ``__getattr__``, or does not exist.

        Args:
            attr (str): Attribute name to check

        Return:
            True if the attribute is already constructed, False otherwise
        """
        with self.lock_attribute_creation_context():
            return hasattr(self, attr)

    def build(self):
        """
        General build method for PropertyBlockDatas. Inheriting models should
        call super().build.

        Args:
            None

        Returns:
            None
        """
        super(ReactionBlockDataBase, self).build()
        add_object_reference(self, "_params", self.config.parameters)

        self._validate_state_block()

    @property
    def params(self):
        return self._params

    def _validate_state_block(self):
        """
        Method to validate that the associated state block matches with the
        PropertyParameterBlock assoicated with the ReactionParameterBlock.
        """
        # Add a reference to the corresponding state block data for later use
        add_object_reference(self, "state_ref", self.config.state_block[self.index()])

        # Validate that property package of state matches that of reaction pack
        if (
            self.config.parameters.config.property_package
            != self.state_ref.config.parameters
        ):
            raise PropertyPackageError(
                "{} the StateBlock associated with this "
                "ReactionBlock does not match with the "
                "PropertyParamterBlock associated with the "
                "ReactionParameterBlock. The modelling framework "
                "does not support mixed associations of property "
                "and reaction packages.".format(self.name)
            )

    def get_reaction_rate_basis(self):
        """
        Method which returns an Enum indicating the basis of the reaction rate
        term.
        """
        return MaterialFlowBasis.other

    def __getattr__(self, attr):
        """
        This method is used to avoid generating unnecessary property
        calculations in reaction blocks. __getattr__ is called whenever a
        property is called for, and if a propery does not exist, it looks for
        a method to create the required property, and any associated
        components.

        Create a property calculation if needed. Return an attrbute error if
        attr == 'domain' or starts with a _ . The error for _ prevents a
        recursion error if trying to get a function to create a property and
        that function doesn't exist.  Pyomo also ocasionally looks for things
        that start with _ and may not exist.  Pyomo also looks for the domain
        attribute, and it may not exist.
        This works by creating a property calculation by calling the "_"+attr
        function.

        A list of __getattr__ calls is maintained in self.__getattrcalls to
        check for recursive loops which maybe useful for debugging. This list
        is cleared after __getattr__ completes successfully.

        Args:
            attr: an attribute to create and return. Should be a property
                  component.
        """
        if self._lock_attribute_creation:
            raise AttributeError(
                f"{attr} does not exist, and attribute creation is locked."
            )

        def clear_call_list(self, attr):
            """Local method for cleaning up call list when a call is handled.

            Args:
                attr: attribute currently being handled
            """
            if self.__getattrcalls[-1] == attr:
                if len(self.__getattrcalls) <= 1:
                    del self.__getattrcalls
                else:
                    del self.__getattrcalls[-1]
            else:
                raise PropertyPackageError(
                    "{} Trying to remove call {} from __getattr__"
                    " call list, however this is not the most "
                    "recent call in the list ({}). This indicates"
                    " a bug in the __getattr__ calls. Please "
                    "contact the IDAES developers with this bug.".format(
                        self.name, attr, self.__getattrcalls[-1]
                    )
                )

        # Check that attr is not something we shouldn't touch
        if attr == "domain" or attr.startswith("_"):
            # Don't interfere with anything by getting attributes that are
            # none of my business
            raise PropertyPackageError(
                "{} {} does not exist, but is a protected "
                "attribute. Check the naming of your "
                "components to avoid any reserved names".format(self.name, attr)
            )

        if attr == "config":
            try:
                self._get_config_args()
                return self.config
            except:
                raise BurntToast(
                    "{} getattr method was triggered by a call "
                    "to the config block, but _get_config_args "
                    "failed. This should never happen."
                )

        # Check for recursive calls
        try:
            # Check if __getattrcalls is initialized
            self.__getattrcalls
        except AttributeError:
            # Initialize it
            self.__getattrcalls = [attr]
        else:
            # Check to see if attr already appears in call list
            if attr in self.__getattrcalls:
                # If it does, indicates a recursive loop.
                if attr == self.__getattrcalls[-1]:
                    # attr method is calling itself
                    self.__getattrcalls.append(attr)
                    raise PropertyPackageError(
                        "{} _{} made a recursive call to "
                        "itself, indicating a potential "
                        "recursive loop. This is generally "
                        "caused by the {} method failing to "
                        "create the {} component.".format(self.name, attr, attr, attr)
                    )
                else:
                    self.__getattrcalls.append(attr)
                    raise PropertyPackageError(
                        "{} a potential recursive loop has been "
                        "detected whilst trying to construct {}. "
                        "A method was called, but resulted in a "
                        "subsequent call to itself, indicating a "
                        "recursive loop. This may be caused by a "
                        "method trying to access a component out "
                        "of order for some reason (e.g. it is "
                        "declared later in the same method). See "
                        "the __getattrcalls object for a list of "
                        "components called in the __getattr__ "
                        "sequence.".format(self.name, attr)
                    )
            # If not, add call to list
            self.__getattrcalls.append(attr)

        # Get property information from get_supported_properties
        try:
            m = self.config.parameters.get_metadata().properties

            if m is None:
                raise PropertyPackageError(
                    "{} reaction package get_supported_properties"
                    " method returned None when trying to create "
                    "{}. Please contact the developer of the "
                    "property package".format(self.name, attr)
                )
        except KeyError:
            # If attr not in get_supported_properties, assume package does not
            # support property
            clear_call_list(self, attr)
            raise PropertyNotSupportedError(
                "{} {} is not supported by reaction package (property is "
                "not listed in get_supported_properties).".format(self.name, attr)
            )

        # Get method name from get_supported_properties
        try:
            if m[attr]["method"] is None:
                # If method is none, property should be constructed
                # by property package, so raise PropertyPackageError
                clear_call_list(self, attr)
                raise PropertyPackageError(
                    "{} {} should be constructed automatically "
                    "by reaction package, but is not present. "
                    "This can be caused by methods being called "
                    "out of order.".format(self.name, attr)
                )
            elif m[attr]["method"] is False:
                # If method is False, package does not support property
                # Raise NotImplementedError
                clear_call_list(self, attr)
                raise PropertyNotSupportedError(
                    "{} {} is not supported by reaction package "
                    "(property method is listed as False in "
                    "get_supported_properties).".format(self.name, attr)
                )
            elif isinstance(m[attr]["method"], str):
                # Try to get method name in from PropertyBlock object
                try:
                    f = getattr(self, m[attr]["method"])
                except AttributeError:
                    # If fails, method does not exist
                    clear_call_list(self, attr)
                    raise PropertyPackageError(
                        "{} {} get_supported_properties method "
                        "returned a name that does not correspond"
                        " to any method in the reaction package. "
                        "Please contact the developer of the "
                        "reaction package.".format(self.name, attr)
                    )
            else:
                # Otherwise method name is invalid
                clear_call_list(self, attr)
                raise PropertyPackageError(
                    "{} {} get_supported_properties method "
                    "returned invalid value for method name. "
                    "Please contact the developer of the "
                    "reaction package.".format(self.name, attr)
                )
        except KeyError:
            # No method key - raise Exception
            # Need to use an AttributeError so Pyomo.DAE will handle this
            clear_call_list(self, attr)
            raise PropertyNotSupportedError(
                "{} get_supported_properties method "
                "does not contain a method for {}. "
                "Please select a package which supports "
                "the necessary properties for your process.".format(self.name, attr)
            )

        # Call attribute if it is callable
        # If this fails, it should return a meaningful error.
        if callable(f):
            try:
                f()
            except Exception:
                # Clear call list and reraise error
                clear_call_list(self, attr)
                raise
        else:
            # If f is not callable, inform the user and clear call list
            clear_call_list(self, attr)
            raise PropertyPackageError(
                "{} tried calling attribute {} in order to create "
                "component {}. However the method is not callable.".format(
                    self.name, f, attr
                )
            )

        # Clear call list, and return
        comp = getattr(self, attr)
        clear_call_list(self, attr)
        return comp

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # Get scaling factor defaults, if no scaling factor set
        for v in self.component_data_objects(
            (Constraint, Var, Expression), descend_into=False
        ):
            if iscale.get_scaling_factor(v) is None:  # don't replace if set
                name = v.getname().split("[")[0]
                index = v.index()
                sf = self.config.parameters.get_default_scaling(name, index)
                if sf is not None:
                    iscale.set_scaling_factor(v, sf)
