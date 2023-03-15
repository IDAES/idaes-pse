#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
# TODO: Missing doc strings
# pylint: disable=missing-module-docstring

# Build on demand needs to play with some private attributes
# pylint: disable=protected-access

from idaes.core.util.exceptions import (
    PropertyPackageError,
    PropertyNotSupportedError,
    BurntToast,
)


__author__ = "John Eslick, Andrew Lee"


def build_on_demand(self, attr):
    """
    This method is used to avoid generating unnecessary property
    calculations in state blocks. __getattr__ is called whenever a
    property is called for, and if a property does not exist, it looks for
    a method to create the required property, and any associated
    components.

    Create a property calculation if needed. Return an attribute error if
    attr == 'domain' or starts with a _ . The error for _ prevents a
    recursion error if trying to get a function to create a property and
    that function doesn't exist.  Pyomo also occasionally looks for things
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
                    f"{self.name} _{attr} made a recursive call to "
                    f"itself, indicating a potential "
                    f"recursive loop. This is generally "
                    f"caused by the {attr} method failing to "
                    f"create the {attr} component."
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

    # Get property information from properties metadata
    try:
        m = self.config.parameters.get_metadata().properties

        if m is None:
            raise PropertyPackageError(
                "{} property package get_metadata()"
                " method returned None when trying to create "
                "{}. Please contact the developer of the "
                "property package".format(self.name, attr)
            )
        n, i = m.get_name_and_index(attr)
        meta = getattr(m, n)[i]
    except (ValueError, AttributeError):
        # If attr not in metadata, assume package does not
        # support property
        clear_call_list(self, attr)
        raise PropertyNotSupportedError(
            "{} {} is not supported by property package (property is "
            "not listed in package metadata properties).".format(self.name, attr)
        )

    # Get method name from resulting properties
    try:
        if not meta.supported:
            # If method is False, package does not support property
            # Raise NotImplementedError
            clear_call_list(self, attr)
            raise PropertyNotSupportedError(
                f"{self.name} {attr} is not supported by property package."
            )
        elif meta.method is None:
            # If method is none, property should be constructed
            # by property package, so raise PropertyPackageError
            clear_call_list(self, attr)
            raise PropertyPackageError(
                "{} {} should be constructed automatically "
                "by property package, but is not present. "
                "This can be caused by methods being called "
                "out of order.".format(self.name, attr)
            )
        elif isinstance(meta.method, str):
            # Try to get method name in from PropertyBlock object
            try:
                f = getattr(self, meta.method)
            except AttributeError:
                # If fails, method does not exist
                clear_call_list(self, attr)
                raise PropertyPackageError(
                    "{} {} package property metadata method "
                    "returned a name that does not correspond"
                    " to any method in the property package. "
                    "Please contact the developer of the "
                    "property package.".format(self.name, attr)
                )
        else:
            # Otherwise method name is invalid
            clear_call_list(self, attr)
            raise PropertyPackageError(
                "{} {} package property metadata method "
                "returned invalid value for method name. "
                "Please contact the developer of the "
                "property package.".format(self.name, attr)
            )
    except KeyError:
        # No method key - raise Exception
        # Need to use an AttributeError so Pyomo.DAE will handle this
        clear_call_list(self, attr)
        raise PropertyNotSupportedError(
            "{} package property metadata method "
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
