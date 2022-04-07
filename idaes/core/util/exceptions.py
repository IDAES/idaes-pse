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


class IdaesError(Exception):
    """
    General exception for all IDAES errors. Intended to be used for multiple
    inheritance in derived exceptions to allow catching of all IDAES-related
    exceptions.
    """

    pass  # Problem with toaster


class BalanceTypeNotSupportedError(NotImplementedError, IdaesError):
    """
    IDAES exception to be used when a control volume does not support a given
    type of balance equation.
    """

    pass  # Tried to put bagel in normal toaster


class ConfigurationError(ValueError, IdaesError):
    """
    IDAES exception to be used when configuration arguments are incorrect
    or inconsistent.
    """

    pass  # Too many buttons, burnt toast


class DynamicError(ValueError, IdaesError):
    """
    IDAES exception for cases where settings associated with dynamic models
    are incorrect.
    """

    pass  # Incorrect browness setting


class BurntToast(IdaesError):
    """
    General exception for when something breaks badly in the core.
    """

    pass  # Toaster on fire


class PropertyNotSupportedError(AttributeError, IdaesError):
    """
    IDAES exception for cases when a models calls for a property which is
    not supported by the chosen property package.

    Needs to inherit from AttributeError for Pyomo interactions.
    """

    pass  # Could not find bread


class PropertyPackageError(AttributeError, IdaesError):
    """
    IDAES exception for generic errors arising from property packages.

    Needs to inherit from AttributeError for Pyomo interactions.
    """

    pass  # Bread stuck


class InitializationError(ArithmeticError, IdaesError):
    """
    IDAES exception to be used when initialization routines fail. All
    initialization routines should raise this Exception if the final step
    fails to converge, and can raise this exception earlier if the routine
    enters a state from which recovery is impossible.
    """

    pass


class UserModelError(ValueError, IdaesError):
    """
    IDAES exception for when a user model returns unphysical values that
    prevent further execution of code.
    """

    pass
