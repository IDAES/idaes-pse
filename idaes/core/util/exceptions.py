# -*- coding: UTF-8 -*-
###############################################################################
# ** Copyright Notice **
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
###############################################################################
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
