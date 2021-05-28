###############################################################################
# Copyright
# =========
#
# Institute for the Design of Advanced Energy Systems Process Systems Engineering
# Framework (IDAES PSE Framework) was produced under the DOE Institute for the
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
# perform publicly and display publicly, and to permit other to do so. Copyright
# (C) 2018-2019 IDAES - All Rights Reserved
#
###############################################################################
"""
Methods for calculating fugacity of Henry's Law components
"""
from pyomo.environ import Var


# -----------------------------------------------------------------------------
# Molarity basis
class Hcp_constant():

    @staticmethod
    def build_parameters(cobj, p):
        cobj.add_component(
            "henry_ref_"+p,
            Var(initialize=cobj.config.parameter_data["henry_ref"][p],
                doc="Henry coeffiicient (molarity basis) at reference state "
                "for phase "+p))

    @staticmethod
    def return_expression(b, p, j, T=None):
        if T is None:
            T = b.temperature

        cobj = b.params.get_component(j)
        H = getattr(cobj, "henry_ref_"+p)

        return (b.mole_frac_phase_comp[p, j]*b.dens_mol_phase[p]*H)


# -----------------------------------------------------------------------------
# Mole fraction basis
class Hxp_constant():

    @staticmethod
    def build_parameters(cobj, p):
        cobj.add_component(
            "henry_ref_"+p,
            Var(initialize=cobj.config.parameter_data["henry_ref"][p],
                doc="Henry coeffiicient (mole fraction basis) at reference "
                "state for phase "+p))

    @staticmethod
    def return_expression(b, p, j, T=None):
        if T is None:
            T = b.temperature

        cobj = b.params.get_component(j)
        H = getattr(cobj, "henry_ref_"+p)

        return (b.mole_frac_phase_comp[p, j]*H)
