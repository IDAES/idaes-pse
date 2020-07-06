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
Mock-up EoS module for testing generic property packages
"""
from pyomo.environ import Var


def common(b, pobj):
    # Create dummy var to be returned by expression calls
    # This Var is used to create expressions where required.
    if not hasattr(b, "dummy_var"):
        b.dummy_var = Var(initialize=42)

    # Counter for how many times this method is called
    # This is used to ensure that the method has been called by checking that
    # the counter has advanced
    if hasattr(b, "eos_common"):
        b.eos_common += 1
    else:
        b.eos_common = 1


def build_parameters(b):
    if not hasattr(b, "dummy_param"):
        b.dummy_param = Var(initialize=42)


def dens_mass_phase(b, p):
    return b.dummy_var


def dens_mol_phase(b, p):
    return b.dummy_var


def enth_mol_phase(b, p):
    return b.dummy_var


def enth_mol_phase_comp(b, p, j):
    return b.dummy_var


def entr_mol_phase(b, p):
    return b.dummy_var


def entr_mol_phase_comp(b, p, j):
    return b.dummy_var


def fug_phase_comp(b, p, j):
    return b.dummy_var


def fug_coeff_phase_comp(b, p, j):
    return b.dummy_var


def gibbs_mol_phase(b, p):
    return b.dummy_var


def gibbs_mol_phase_comp(b, p, j):
    return b.dummy_var
