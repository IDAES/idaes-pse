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
This module contains miscellaneous utility functions for use in IDAES models.
"""
from enum import Enum

import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


# Author: Andrew Lee
def add_object_reference(self, local_name, remote_object):
    """
    Method to create a reference in the local model to a remote Pyomo object.
    This method should only be used where Pyomo Reference objects are not
    suitable (such as for referencing scalar Pyomo objects where the None
    index is undesirable).

    Args:
        local_name : name to use for local reference (str)
        remote_object : object to make a reference to

    Returns:
        None
    """
    try:
        object.__setattr__(self, local_name, remote_object)
    except AttributeError:
        raise AttributeError(
            "{} failed to construct reference to {} - remote "
            "object does not exist.".format(self.name, remote_object)
        )


# Author: Jaffer Ghouse
def extract_data(data_dict):
    """
    General method that returns a rule to extract data from a python
    dictionary. This method allows the param block to have a database for
    a parameter but extract a subset of this data to initialize a Pyomo
    param object.
    """

    def _rule_initialize(m, *args):
        if len(args) > 1:
            return data_dict[args]
        else:
            return data_dict[args[0]]

    return _rule_initialize


def set_param_from_config(b, param, config=None, index=None):
    """
    Utility method to set parameter value from a config block. This allows for
    converting units if required. This method directly sets the value of the
    parameter.

    This method supports three forms for defining the parameter value:

    1. a 2-tuple of the form (value, units) where units are the units that
    the value are defined in
    2. a float where the float is assumed to be the value of the parameter
    value in the base units of the property package
    3. a Python Class which has a get_parameter_value method which will
    return a 2-tuple of (value, units) based on a lookup of the parameter name

    Args:
        b - block on which parameter and config block are defined
        param - name of parameter as str. Used to find param and config arg
        units - units of param object (used if conversion required)
        config - (optional) config block to get parameter data from. If
                unset, assumes b.config.
        index - (optional) used for pure component properties where a single
                property may have multiple parameters associated with it.

    Returns:
        None
    """
    if config is None:
        try:
            config = b.config
        except AttributeError:
            raise AttributeError(
                "{} - set_param_from_config method was not provided with a "
                "config argument, but no default Config block exists. Please "
                "specify the Config block to use via the config argument.".format(
                    b.name
                )
            )

    # Check that config is an instance of a Config Block
    if not isinstance(config, ConfigBlock):
        raise TypeError(
            "{} - set_param_from_config - config argument provided is not an "
            "instance of a Config Block.".format(b.name)
        )

    if index is None:
        try:
            param_obj = getattr(b, param)
        except AttributeError:
            raise AttributeError(
                "{} - set_param_from_config method was provided with param "
                "argument {}, but no attribute of that name exists.".format(
                    b.name, param
                )
            )

        try:
            p_data = config.parameter_data[param]
        except (KeyError, AttributeError):
            raise KeyError(
                "{} - set_param_from_config method was provided with param "
                "argument {}, but the config block does not contain a "
                "value for this parameter.".format(b.name, param)
            )
    else:
        try:
            param_obj = getattr(b, param + "_" + index)
        except AttributeError:
            raise AttributeError(
                "{} - set_param_from_config method was provided with param and"
                " index arguments {} {}, but no attribute with that "
                "combination ({}_{}) exists.".format(b.name, param, index, param, index)
            )

        try:
            p_data = config.parameter_data[param][index]
        except (KeyError, AttributeError):
            raise KeyError(
                "{} - set_param_from_config method was provided with param and"
                " index arguments {} {}, but the config block does not contain"
                " a value for this parameter and index.".format(b.name, param, index)
            )

    units = param_obj.get_units()

    # Check to see if p_data is callable, and if so, try to call the
    # get_parameter_value method to get 2-tuple
    if hasattr(p_data, "get_parameter_value"):
        p_data = p_data.get_parameter_value(b.local_name, param)

    if isinstance(p_data, tuple):
        # 11 Dec 2020 - There is currently a bug in Pyomo where trying to
        # convert the units of a unitless quantity results in a TypeError.
        # To avoid this, we check here for cases where both the parameter and
        # user provided value are unitless and bypass unit conversion.
        if (units is None or units is pyo.units.dimensionless) and (
            p_data[1] is None or p_data[1] is pyo.units.dimensionless
        ):
            param_obj.value = p_data[0]
        else:
            param_obj.value = pyo.units.convert_value(
                p_data[0], from_units=p_data[1], to_units=units
            )
    else:
        _log.debug(
            "{} no units provided for parameter {} - assuming default "
            "units".format(b.name, param)
        )
        param_obj.value = p_data


class StrEnum(str, Enum):
    """
    Multiple inheritance string-Enum for representing Enums with string values
    """

    def __str__(self):
        return str(self.value)
