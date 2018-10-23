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
This module contains custom IDAES exceptions.
"""

__author__ = "Andrew Lee"


def ConfigurationError(ValueError):
    """
    IDAES exception to be used when configuration argumnet are incorrect
    or inconsistent.
    """
    pass  # Too many buttons, burnt toast


def DynamicError(ValueError):
    """
    IDAES exception for cases where settings associated with dynamic models
    are incorrect.
    """
    pass  # Incorrect browness setting


def PropertyNotSupportedError(NotImplementedError):
    """
    IDAES exception for cases when a models calls for a property which is
    not supported by the chosen property package.
    """
    pass  # Could not find bread


def PropertyPackageError(Exception):
    """
    IDAES exception for generic errors arising from property packages.
    """
    pass  # Bread stuck
