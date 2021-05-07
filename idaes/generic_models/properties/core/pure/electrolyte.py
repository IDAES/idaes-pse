##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
General pure component electrolyte methods
"""
from pyomo.environ import Var, units as pyunits

from idaes.core.util.misc import set_param_from_config


# -----------------------------------------------------------------------------
# Method for constant relative permittivity
# RElative permittivity was referred to as the dielectric constant in the past
class relative_permittivity_constant():

    @staticmethod
    def build_parameters(cobj):
        cobj.relative_permittivity_liq_comp = Var(
            doc="Relative permittivity",
            units=pyunits.dimensionless)
        set_param_from_config(cobj, param="relative_permittivity_liq_comp")

    @staticmethod
    def return_expression(b, cobj, T):
        return cobj.relative_permittivity_liq_comp
