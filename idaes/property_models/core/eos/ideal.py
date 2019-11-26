##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
"""
from idaes.core.util.exceptions import PropertyNotSupportedError
from idaes.property_models.core.generic.generic_property import \
        GenericPropertyPackageError


def common(b):
    # No common components required for ideal property calculations
    pass


def dens_mass_phase(b, p):
    return b.dens_mol_phase[p]*b.mw_phase[p]


def dens_mol_phase(b, p):
    if p == "Vap":
        return b.pressure/(b._params.gas_const*b.temperature)
    elif p == "Liq":
        if b._params.config.dens_mol_comp_liq.dens_mol_liq is None:
            raise GenericPropertyPackageError(b, "dens_mol_comp_liq")
        return sum(b.mole_frac_phase_comp[p, j] *
                   b._params.config.dens_mol_comp_liq.dens_mol_liq(
                           b, b.temperature, j)
                   for j in b._params.component_list)
    else:
        raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))


def enth_mol_phase(b, p):
    return sum(b.mole_frac_phase_comp[p, j]*b.enth_mol_phase_comp[p, j]
               for j in b._params.component_list)


def enth_mol_phase_comp(b, p, j):
    if p == "Vap":
        if b._params.config.enth_mol_comp_ig.enth_mol_ig is None:
            raise GenericPropertyPackageError(b, "enth_mol_comp_vap")
        return b._params.config.enth_mol_comp_ig.enth_mol_ig(
                    b, j, b.temperature) + \
               b._params.dh_vap_ref[j]
    elif p == "Liq":
        if b._params.config.enth_mol_comp_liq.enth_mol_liq is None:
            raise GenericPropertyPackageError(b, "enth_mol_comp_liq")
        return b._params.config.enth_mol_comp_liq.enth_mol_liq(
                b, j, b.temperature)
    else:
        raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))


def entr_mol_phase(b, p):
    return sum(b.mole_frac_phase_comp[p, j]*b.entr_mol_phase_comp[p, j]
               for j in b._params.component_list)


def entr_mol_phase_comp(b, p, j):
    if p == "Vap":
        if b._params.config.entr_mol_comp_ig.entr_mol_ig is None:
            raise GenericPropertyPackageError(b, "entr_mol_comp_ig")
        return b._params.config.entr_mol_comp_ig.entr_mol_ig(
                b, j, b.temperature) + \
               b._params.ds_vap_ref[j]
    elif p == "Liq":
        if b._params.config.entr_mol_comp_liq.entr_mol_liq is None:
            raise GenericPropertyPackageError(b, "entr_mol_comp_liq")
        return b._params.config.entr_mol_comp_liq.entr_mol_liq(
                b, j, b.temperature)
    else:
        raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))


def entr_mol_phase_comp_ref(b, p, j):
    if p == "Vap":
        if b._params.config.entr_mol_comp_ig.entr_mol_ig is None:
            raise GenericPropertyPackageError(b, "entr_mol_comp_ig")
        return b._params.config.entr_mol_comp_ig.entr_mol_ig(
                b, j, b._params.temperature_ref) +\
               b._params.ds_vap_ref[j]
    elif p == "Liq":
        if b._params.config.entr_mol_comp_liq.entr_mol_liq is None:
            raise GenericPropertyPackageError(b, "entr_mol_comp_liq")
        return b._params.config.entr_mol_comp_liq.entr_mol_liq(
                b, j, b._params.temperature_ref)
    else:
        raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))


def fugacity(b, p, j):
    if p == "Vap":
        return b.mole_frac_phase_comp[p, j]*b.pressure
    elif p == "Liq":
        if b._params.config.pressure_sat_comp.pressure_sat is None:
            raise GenericPropertyPackageError(b, "pressure_sat_comp")
        return b.mole_frac_phase_comp[p, j] * \
               b._params.config.pressure_sat_comp.pressure_sat(
                       b, b.temperature, j)
    else:
        raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))


def fug_coeff(b, p, j):
    return 1


def gibbs_mol_phase(b, p):
    return sum(b.mole_frac_phase_comp[p, j]*b.gibbs_mol_phase_comp[p, j]
               for j in b._params.component_list)


def gibbs_mol_phase_comp(b, p, j):
    return (b.enth_mol_phase_comp[p, j] -
            (b.entr_mol_phase_comp[p, j] -
             b._params.config.equation_of_state[p]
             .entr_mol_phase_comp_ref(b, p, j)) *
            b.temperature)


def _invalid_phase_msg(name, phase):
    return ("{} recieved unrecognised phase name {}. Ideal property "
            "libray only supports Vap and Liq phases."
            .format(name, phase))
