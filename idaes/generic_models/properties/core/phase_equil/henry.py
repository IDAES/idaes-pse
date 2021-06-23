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
Methods for calculating fugacity of Henry's Law components

For now, only support mole fraction basis
"""
from pyomo.environ import Var


class ConstantH():

    @staticmethod
    def build_parameters(cobj, p):
        b = cobj.parent_block()
        units = b.get_metadata().derived_units

        cobj.add_component(
            "henry_ref_"+p,
            Var(initialize=cobj.config.parameter_data["henry_ref"][p],
                doc="Henry coefficient (mole fraction basis) at reference "
                "state for phase "+p,
                units=units["pressure"]))

    @staticmethod
    def return_expression(b, p, j, T=None):
        cobj = b.params.get_component(j)
        H = getattr(cobj, "henry_ref_"+p)

        return H

    @staticmethod
    def dT_expression(b, p, j, T=None):
        return 0
