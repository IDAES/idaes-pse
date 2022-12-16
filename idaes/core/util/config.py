# -*- coding: UTF-8 -*-
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
This module contains utility functions useful for validating arguments to
IDAES modeling classes. These functions are primarily designed to be used as
the `domain` argument in ConfigBlocks.
"""

__author__ = "Andrew Lee"

from pyomo.common.config import Bool

from pyomo.environ import Set
from pyomo.dae import ContinuousSet
from pyomo.network import Port
from idaes.core import useDefault
from idaes.core.util.exceptions import ConfigurationError

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def is_physical_parameter_block(val):
    """Domain validator for property package attributes

    Args:
        val : value to be checked

    Returns:
        ConfigurationError if val is not an instance of PhysicalParameterBlock
        or useDefault
    """
    from idaes.core.base.property_base import PhysicalParameterBlock

    if isinstance(val, PhysicalParameterBlock) or val == useDefault:
        return val
    else:
        _log.error(
            "Property package argument {} should == useDefault or "
            "be an instance of PhysicalParameterBlock".format(val)
        )
        raise ConfigurationError(
            """Property package argument should be an instance
                of a PhysicalParameterBlock or useDefault"""
        )


def is_reaction_parameter_block(val):
    """Domain validator for reaction package attributes

    Args:
        val : value to be checked

    Returns:
        ConfigurationError if val is not an instance of ReactionParameterBlock
    """
    from idaes.core.base.reaction_base import ReactionParameterBlock

    if isinstance(val, ReactionParameterBlock):
        return val
    else:
        raise ConfigurationError(
            """Reaction package argument should be an instance
                of a ReactionParameterBlock"""
        )


def is_state_block(val):
    """Domain validator for state block as an argument

    Args:
        val : value to be checked

    Returns:
        ConfigurationError if val is not an instance of StateBlock
        or None
    """
    from idaes.core.base.property_base import StateBlock

    if isinstance(val, StateBlock) or val is None:
        return val
    else:
        raise ConfigurationError(
            """State block should be an instance of a StateBlock or
                None"""
        )


def is_port(arg):
    """Domain validator for ports

    Args:
        arg : argument to be checked as a Port

    Returns:
        Port object or Exception
    """
    if not isinstance(arg, Port):
        raise ConfigurationError(
            "Invalid argument type. Expected an instance " "of a Pyomo Port object"
        )
    return arg


def is_time_domain(arg):
    """Domain validator for time domains

    Args:
        arg : argument to be checked as a time domain (i.e. Set or
        ContinuousSet)

    Returns:
        Set, ContinuousSet or Exception
    """
    if not isinstance(arg, (Set, ContinuousSet)):
        raise ConfigurationError(
            "Invalid argument type. Expected an instance "
            "of a Pyomo Set or ContinuousSet object"
        )
    return arg


def is_transformation_method(arg):
    """Domain validator for transformation methods

    Args:
        arg : argument to be checked for membership in recognized strings

    Returns:
        Recognised string or Exception
    """
    if arg in ["dae.finite_difference", "dae.collocation"]:
        return arg
    else:
        raise ConfigurationError(
            "Invalid value provided for transformation_method. "
            "Please check the value and spelling of the argument provided."
        )


def is_transformation_scheme(arg):
    """Domain validator for transformation scheme

    Args:
        arg : argument to be checked for membership in recognized strings

    Returns:
        Recognised string or Exception
    """
    if arg in ["BACKWARD", "FORWARD", "LAGRANGE-RADAU", "LAGRANGE-LEGENDRE"]:
        return arg
    else:
        raise ConfigurationError(
            "Invalid value provided for transformation_scheme. "
            "Please check the value and spelling of the argument provided."
        )


def DefaultBool(arg):
    """
    Domain validator for bool-like arguments with a 'use default' option.
    Relies on Pyomo's Bool() validator for identifying bool-like arguments.

    Args:
        arg : argument to be validated.

    Returns:
        A bool or useDefault

    Raises:
        ValueError if arg is not bool-like or useDefault
    """
    if arg == useDefault:
        return arg
    else:
        return Bool(arg)
