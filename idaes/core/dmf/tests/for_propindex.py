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
from idaes.core.base.property_base import PhysicalParameterBlock
from pyomo.environ import units as pyunits


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
