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
