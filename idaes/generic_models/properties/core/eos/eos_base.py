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
Base class for EoS modules.

Raises NotImplementedErrors for all expected methods in case developer misses
some. EoS developers should overload all these methods.
"""
from pyomo.environ import units as pyunits
from idaes.core.util.constants import Constants as const


class EoSBase():

    @staticmethod
    def gas_constant(b):
        # Utility method to convert gas constant to base units
        base_units = b.params.get_metadata().default_units

        r_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["temperature"]**-1 *
                   base_units["amount"]**-1 *
                   base_units["time"]**-2)

        return pyunits.convert(const.gas_constant, to_units=r_units)

    @staticmethod
    def common(b, pobj):
        raise NotImplementedError(_msg(b, "common"))

    @staticmethod
    def build_parameters(b):
        raise NotImplementedError(_msg(b, "build_parameters"))

    @staticmethod
    def dens_mass_phase(b, p):
        raise NotImplementedError(_msg(b, "dens_mass_phase"))

    @staticmethod
    def dens_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "dens_mol_phase"))

    @staticmethod
    def enth_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "enth_mol_phase"))

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "enth_mol_phase_comp"))

    @staticmethod
    def entr_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "entr_mol_phase"))

    @staticmethod
    def entr_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "entr_mol_phase_comp"))

    @staticmethod
    def fug_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "fug_phase_comp"))

    @staticmethod
    def fug_phase_comp_eq(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_eq"))

    @staticmethod
    def log_fug_phase_comp_eq(b, p, j, pp):
        raise NotImplementedError(_msg(b, "log_fug_phase_comp_eq"))

    @staticmethod
    def fug_coeff_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "fug_coeff_phase_comp"))

    @staticmethod
    def fug_coeff_phase_comp_eq(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_coeff_phase_comp_eq"))

    @staticmethod
    def fug_phase_comp_Tbub(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_Tbub"))

    @staticmethod
    def fug_phase_comp_Tdew(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_Tdew"))

    @staticmethod
    def fug_phase_comp_Pbub(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_Pbub"))

    @staticmethod
    def fug_phase_comp_Pdew(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_Pdew"))

    @staticmethod
    def gibbs_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "gibbs_mol_phase"))

    @staticmethod
    def gibbs_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "gibbs_mol_phase_comp"))


def _msg(b, attr):
    return ("{} Equation of State module has not implemented a method for {}. "
            "Please contact the EoS developer or use a different module."
            .format(b.name, attr))
