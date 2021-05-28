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
Mock-up EoS module for testing generic property packages
"""
from pyomo.environ import Var

from idaes.generic_models.properties.core.eos.eos_base import EoSBase


class DummyEoS(EoSBase):
    # Add attribute indicating support for electrolyte systems
    electrolyte_support = True

    def common(b, pobj):
        # Create dummy var to be returned by expression calls
        # This Var is used to create expressions where required.
        if not hasattr(b, "dummy_var"):
            b.dummy_var = Var(initialize=42)

        # Counter for how many times this method is called
        # This is used to ensure that the method has been called by checking
        # that the counter has advanced
        if hasattr(b, "eos_common"):
            b.eos_common += 1
        else:
            b.eos_common = 1

    def calculate_scaling_factors(b, pobj):
        pass

    def build_parameters(b):
        if not hasattr(b, "dummy_param"):
            b.dummy_param = Var(initialize=42)

    def dens_mass_phase(b, p):
        return 42

    def dens_mol_phase(b, p):
        return 55e3

    def energy_internal_mol_phase(b, p):
        return 2e2*b.temperature

    def energy_internal_mol_phase_comp(b, p, j):
        return 2e2*b.temperature

    def enth_mol_phase(b, p):
        return 1e2*b.temperature

    def enth_mol_phase_comp(b, p, j):
        return 1e2*b.temperature

    def entr_mol_phase(b, p):
        return 42

    def entr_mol_phase_comp(b, p, j):
        return 42

    def fug_phase_comp(b, p, j):
        return 42

    def fug_coeff_phase_comp(b, p, j):
        return 42

    def gibbs_mol_phase(b, p):
        return 42

    def gibbs_mol_phase_comp(b, p, j):
        return 42
