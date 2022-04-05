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
Mock-up EoS module for testing generic property packages
"""
from pyomo.environ import Var, sqrt, units as pyunits

from idaes.models.properties.modular_properties.eos.eos_base import EoSBase


class DummyEoS(EoSBase):
    # Add attribute indicating support for electrolyte systems
    electrolyte_support = True

    @staticmethod
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

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def build_parameters(b):
        if not hasattr(b, "dummy_param"):
            b.dummy_param = Var(initialize=42)

    @staticmethod
    def act_phase_comp(b, p, j):
        return 42

    @staticmethod
    def act_coeff_phase_comp(b, p, j):
        return 1

    @staticmethod
    def compress_fact_phase(b, p):
        return 42

    @staticmethod
    def cp_mol_phase(b, p):
        return 42

    @staticmethod
    def cp_mol_phase_comp(b, p, j):
        return 42

    @staticmethod
    def cv_mol_phase(b, p):
        return 42

    @staticmethod
    def cv_mol_phase_comp(b, p, j):
        return 42

    @staticmethod
    def dens_mass_phase(b, p):
        return 42

    @staticmethod
    def dens_mol_phase(b, p):
        return 55e3 * pyunits.mol / pyunits.m**3

    @staticmethod
    def energy_internal_mol_phase(b, p):
        return 2e2 * b.temperature

    @staticmethod
    def energy_internal_mol_phase_comp(b, p, j):
        return 2e2 * b.temperature

    @staticmethod
    def enth_mol_phase(b, p):
        return 1e2 * b.temperature

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        return 1e2 * b.temperature

    @staticmethod
    def entr_mol_phase(b, p):
        return 42

    @staticmethod
    def entr_mol_phase_comp(b, p, j):
        return 42

    @staticmethod
    def fug_phase_comp(b, p, j):
        return 42

    @staticmethod
    def fug_coeff_phase_comp(b, p, j):
        return 42

    @staticmethod
    def gibbs_mol_phase(b, p):
        return 42

    @staticmethod
    def gibbs_mol_phase_comp(b, p, j):
        return 42

    @staticmethod
    def isothermal_speed_sound_phase(b, p):
        return 250

    @staticmethod
    def isentropic_speed_sound_phase(b, p):
        return sqrt(b.heat_capacity_ratio_phase[p]) * b.isothermal_speed_sound_phase[p]

    @staticmethod
    def vol_mol_phase(b, p):
        return 42

    @staticmethod
    def vol_mol_phase_comp(b, p, j):
        return 42

    @staticmethod
    def log_act_phase_comp(b, p, j):
        return 1

    @staticmethod
    def log_act_phase_solvents(b, p):
        return 1
