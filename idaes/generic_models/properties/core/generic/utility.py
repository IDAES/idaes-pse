##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Common methods used by generic framework

Author: A Lee
"""

import types

from pyomo.environ import units as pyunits

from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


class GenericPropertyPackageError(PropertyPackageError):
    # Error message for when a property is called for but no option provided
    def __init__(self, block, prop):
        self.prop = prop
        self.block = block

    def __str__(self):
        return f"Generic Property Package instance {self.block} called for " \
               f"{self.prop}, but was not provided with a method " \
               f"for this property. Please add a method for this property " \
               f"in the property parameter configuration."


def get_method(self, config_arg, comp=None):
    """
    Method to inspect configuration argument and return the user-defined
    construction method associated with it.

    This method checks whether the value provided is a method or a
    module. If the value is a module, it looks in the module for a method
    with the same name as the config argument and returns this. If the
    value is a method, the method is returned. If the value is neither a
    module or a method, an ConfigurationError is raised.

    Args:
        config_arg : the configuration argument to look up

    Returns:
        A callable method or a ConfigurationError
    """
    if comp is None:
        source_block = self.params.config
    else:
        source_block = self.params.get_component(comp).config

    try:
        c_arg = getattr(source_block, config_arg)
    except AttributeError:
        raise AttributeError("{} Generic Property Package called for invalid "
                             "configuration option {}. Please contact the "
                             "developer of the property package."
                             .format(self.name, config_arg))

    if c_arg is None:
        raise GenericPropertyPackageError(self, config_arg)

    if isinstance(c_arg, types.ModuleType):
        c_arg = getattr(c_arg, config_arg)

    try:
        mthd = c_arg.return_expression
    except AttributeError:
        mthd = c_arg

    if callable(mthd):
        return mthd
    else:
        raise ConfigurationError(
                "{} Generic Property Package received invalid value "
                "for argument {}. Value must be a method, a class with a "
                "method named expression or a module containing one of the "
                "previous.".format(self.name, config_arg))


def get_component_object(self, comp):
    """
    Utility method to get a component object from the property parameter block.
    This code is used frequently throughout the generic property pacakge
    libraries.

    Args:
        comp: name of the component object to be returned.

    Returns:
        Component: Component object with name comp.

    """
    return self.params.get_component(comp)


def get_bounds_from_config(b, state, base_units):
    """
    Method to take a 3- or 4-tuple state definition config argument and return
    tuples for the bounds and default value of the Var object.

    Expects the form (lower, default, upper, units) where units is optional

    Args:
        b - StateBlock on which the state vars are to be constructed
        state - name of state var as a string (to be matched with config dict)
        base_units - base units of state var to be used if conversion required

    Returns:
        bounds - 2-tuple of state var bounds in base units
        default_val - default value of state var in base units
    """
    try:
        var_config = b.params.config.state_bounds[state]
    except (KeyError, TypeError):
        # State definition missing
        return (None, None), None

    if len(var_config) == 4:
        # Units provided, need to convert values
        bounds = (pyunits.convert_value(var_config[0],
                                        from_units=var_config[3],
                                        to_units=base_units),
                  pyunits.convert_value(var_config[2],
                                        from_units=var_config[3],
                                        to_units=base_units))
        default_val = pyunits.convert_value(var_config[1],
                                            from_units=var_config[3],
                                            to_units=base_units)
    else:
        bounds = (var_config[0], var_config[2])
        default_val = var_config[1]

    return bounds, default_val
