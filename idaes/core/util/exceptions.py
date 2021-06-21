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
This module contains custom IDAES exceptions.
"""

__author__ = "Andrew Lee"


class BalanceTypeNotSupportedError(NotImplementedError):
    """
    IDAES exception to be used when a control volumedoes not support a given
    type of balance equation.
    """
    pass  # Tried to put bagel in normal toaster


class ConfigurationError(ValueError):
    """
    IDAES exception to be used when configuration arguments are incorrect
    or inconsistent.
    """
    pass  # Too many buttons, burnt toast


class DynamicError(ValueError):
    """
    IDAES exception for cases where settings associated with dynamic models
    are incorrect.
    """
    pass  # Incorrect browness setting


class BurntToast(Exception):
    """
    General exception for when something breaks badly in the core.
    """
    pass  # Toaster on fire


class PropertyNotSupportedError(AttributeError):
    """
    IDAES exception for cases when a models calls for a property which is
    not supported by the chosen property package.

    Needs to inherit from AttributeError for Pyomo interactions.
    """
    pass  # Could not find bread


class PropertyPackageError(AttributeError):
    """
    IDAES exception for generic errors arising from property packages.

    Needs to inherit from AttributeError for Pyomo interactions.
    """
    pass  # Bread stuck
