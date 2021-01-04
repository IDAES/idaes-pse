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
Methods for ideal equations of state.

Currently only supports liquid and vapor phases
"""
from pyomo.environ import log

from idaes.core.util.exceptions import PropertyNotSupportedError
from idaes.generic_models.properties.core.generic.utility import (
    get_method, get_component_object as cobj)
from .eos_base import EoSBase


# TODO: Add support for ideal solids
class Ideal(EoSBase):

    @staticmethod
    def common(b, pobj):
        # No common components required for ideal property calculations
        pass

    @staticmethod
    def build_parameters(b):
        # No EoS specific parameters required
        pass

    @staticmethod
    def compress_fact_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase() or pobj.is_liquid_phase():
            return 1
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def dens_mass_phase(b, p):
        return b.dens_mol_phase[p]*b.mw_phase[p]

    @staticmethod
    def dens_mol_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return b.pressure/(Ideal.gas_constant(b)*b.temperature)
        elif pobj.is_liquid_phase():
            return sum(b.mole_frac_phase_comp[p, j] *
                       get_method(b, "dens_mol_liq_comp", j)(
                           b, cobj(b, j), b.temperature)
                       for j in b.components_in_phase(p))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def cp_mol_phase(b, p):
        return sum(b.mole_frac_phase_comp[p, j]*b.cp_mol_phase_comp[p, j]
                   for j in b.components_in_phase(p))

    @staticmethod
    def cp_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return get_method(b, "cp_mol_ig_comp", j)(
                b, cobj(b, j), b.temperature)
        elif pobj.is_liquid_phase():
            return get_method(b, "cp_mol_liq_comp", j)(
                b, cobj(b, j), b.temperature)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def enth_mol_phase(b, p):
        return sum(b.mole_frac_phase_comp[p, j]*b.enth_mol_phase_comp[p, j]
                   for j in b.components_in_phase(p))

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return get_method(b, "enth_mol_ig_comp", j)(
                b, cobj(b, j), b.temperature)
        elif pobj.is_liquid_phase():
            return get_method(b, "enth_mol_liq_comp", j)(
                b, cobj(b, j), b.temperature)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def entr_mol_phase(b, p):
        return sum(b.mole_frac_phase_comp[p, j]*b.entr_mol_phase_comp[p, j]
                   for j in b.components_in_phase(p))

    @staticmethod
    def entr_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return (get_method(b, "entr_mol_ig_comp", j)(
                b, cobj(b, j), b.temperature) -
                Ideal.gas_constant(b)*log(
                    b.mole_frac_phase_comp[p, j]*b.pressure /
                    b.params.pressure_ref))
        elif pobj.is_liquid_phase():
            # Assume no pressure/volume dependecy of entropy for ideal liquids
            return (get_method(b, "entr_mol_liq_comp", j)(
                b, cobj(b, j), b.temperature))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def fug_phase_comp(b, p, j):
        return _fug_phase_comp(b, p, j, b.temperature)

    @staticmethod
    def fug_phase_comp_eq(b, p, j, pp):
        return _fug_phase_comp(b, p, j, b._teq[pp])

    @staticmethod
    def log_fug_phase_comp_eq(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return log(b.mole_frac_phase_comp[p, j]) + log(b.pressure)
        elif pobj.is_liquid_phase():
            return (log(b.mole_frac_phase_comp[p, j]) +
                    log(get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature)))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def fug_coeff_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return 1

    @staticmethod
    def fug_coeff_phase_comp_eq(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return 1

    @staticmethod
    def log_fug_coeff_phase_comp_Tbub(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return log(1)

    @staticmethod
    def log_fug_coeff_phase_comp_Tdew(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return log(1)

    @staticmethod
    def log_fug_coeff_phase_comp_Pbub(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return log(1)

    @staticmethod
    def log_fug_coeff_phase_comp_Pdew(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return log(1)

    @staticmethod
    def gibbs_mol_phase(b, p):
        return sum(b.mole_frac_phase_comp[p, j]*b.gibbs_mol_phase_comp[p, j]
                   for j in b.components_in_phase(p))

    @staticmethod
    def gibbs_mol_phase_comp(b, p, j):
        return (b.enth_mol_phase_comp[p, j] -
                b.entr_mol_phase_comp[p, j] *
                b.temperature)


def _invalid_phase_msg(name, phase):
    return ("{} received unrecognised phase name {}. Ideal property "
            "libray only supports Vap and Liq phases."
            .format(name, phase))


def _fug_phase_comp(b, p, j, T):
    pobj = b.params.get_phase(p)
    if pobj.is_vapor_phase():
        return b.mole_frac_phase_comp[p, j] * b.pressure
    elif pobj.is_liquid_phase():
        return (b.mole_frac_phase_comp[p, j] *
                get_method(b, "pressure_sat_comp", j)(
                    b, cobj(b, j), T))
    else:
        raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
