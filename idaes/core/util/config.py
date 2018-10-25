# -*- coding: UTF-8 -*-
##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
# 
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
This module contains utility functions useful for validating arguments to
IDAES modeling classes. These functions are primarily designed to be used as
the `domain` argument in ConfigBlocks.
"""

__author__ = "Andrew Lee"

from pyomo.network import Port
from idaes.core.property_base import PropertyParameterBase


def is_parameter_block(val):
    '''Domain validator for property package attributes

    Args:
        val : value to be checked

    Returns:
        TypeError if val is not an instance of PropertyParameterBase,
        'use_parameter_block' or None
    '''
    if (isinstance(val, PropertyParameterBase) or
            val in (None, 'use_parent_value')):
        return val
    else:
        raise TypeError("""Property package argument should be an instance
                        of a Property Parameter Block, or have a value of
                        'use_parent_value' or None""")


def list_of_floats(arg):
    '''Domain validator for lists of floats

    Args:
        arg : argument to be cast to list of floats and validated

    Returns:
        List of strings
    '''
    try:
        # Assume arg is iterable
        lst = [float(i) for i in arg]
    except TypeError:
        # arg is not iterable
        lst = [float(arg)]
    return lst


def list_of_strings(arg):
    '''Domain validator for lists of strings

    Args:
        arg : argument to be cast to list of strings and validated

    Returns:
        List of strings
    '''
    if isinstance(arg, dict):
        raise ValueError("list_of_strings cannot cast dict")

    try:
        # Assume arg is iterable
        if isinstance(arg, str):
            lst = [arg]
        else:
            lst = [str(i) for i in arg]
    except TypeError:
        # arg is not iterable
        lst = [str(arg)]
    return lst


def is_port(arg):
    '''Domain validator for ports

    Args:
        arg : argument to be checked as a Port

    Returns:
        Port object or Exception
    '''
    if not isinstance(arg, Port):
        raise TypeError('One side of the Stream is not a Port. Check '
                        'that both source and destination are Port '
                        'objects.')
    return arg
