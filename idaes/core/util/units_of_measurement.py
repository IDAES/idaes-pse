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
This module contains utility functions for converting units of measurement and
reporting model outputs.
"""
from pyomo.environ import as_quantity, units
import idaes


def report_quantity(c):
    q = as_quantity(c)

    return convert_quantity_to_reporting_units(q)


def convert_quantity_to_reporting_units(q):
    """
    Converts a pint quantity to the units defined in the IDAES config block.

    Args:
        q: pint quantity to be converted

    Returns:
        A new pint quantity in the units defined by the IDAES config block.
    """
    # First, check if quantity is dimensionless
    if q.dimensionless:
        # No need to do anything here
        return q

    # Get definition of desired units from config block
    def_units = idaes.cfg.reporting_units

    # Get dimensionality of q
    dim = q.dimensionality

    # Iterate through unit definition to try to find matching dimensionality
    for u in def_units.values():
        # Get pint unit object from string
        u_obj = getattr(units.pint_registry, u)
        if str(u_obj.dimensionality) == str(dim):
            # Found matching dimensionality
            q2 = q.to(u_obj)
            break
    else:
        # No matching dimensionality found, fall back to default system of units
        q2 = q.to_base_units()

    return q2
