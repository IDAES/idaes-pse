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
# pylint: disable=missing-class-docstring

from pyomo.environ import units as pyunits
from idaes.core.base.property_base import PhysicalParameterBlock


class IndexMePlease1(PhysicalParameterBlock):
    @classmethod
    def define_metadata(cls, m):
        m.add_default_units(
            {
                "temperature": pyunits.K,
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
            }
        )
        m.add_properties(
            {
                "pressure": {"units": "Pa", "method": "foo"},
                "temperature": {"method": "bar"},
            }
        )
