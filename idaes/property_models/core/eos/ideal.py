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


def common(b):
    # No common componets required for ideal property calculations
    pass


def dens_mass(b, p):
    return b.dens_mol_phase[p]*b.mw_phase[p]


def dens_mol(b, p):
    if p == "Vap":
        return b.pressure/(b._params.gas_const*b.temperature)
    elif p == "Liq":
        return sum(b.mole_frac_phase[p, j] *
                   b._params.config.dens_mol_liq(b, b.temperature, j)
                   for j in b._params.component_list)
    else:
        raise PropertyNotSupportedError(
                "{} recieved unrecognised phase name {}. Ideal property "
                "libray only supports Vap and Liq phases."
                .format(b.name, p))


def enth_mol_comp(b, p, j):
    if p == "Vap":
        return b._params.config.enth_mol_vap(b, j, b.temperature) + \
                b._params.dh_vap_ref[j]
    elif p == "Liq":
        return b._params.config.enth_mol_liq(b, j, b.temperature)
    else:
        raise PropertyNotSupportedError(
                "{} recieved unrecognised phase name {}. Ideal property "
                "libray only supports Vap and Liq phases."
                .format(b.name, p))


def entr_mol_comp(b, p, j):
    if p == "Vap":
        return b._params.config.entr_mol_vap(b, j, b.temperature) + \
                b._params.ds_vap_ref[j]
    elif p == "Liq":
        return b._params.config.entr_mol_liq(b, j, b.temperature)
    else:
        raise PropertyNotSupportedError(
                "{} recieved unrecognised phase name {}. Ideal property "
                "libray only supports Vap and Liq phases."
                .format(b.name, p))


def entr_mol_comp_ref(b, p, j):
    if p == "Vap":
        return b._params.config.entr_mol_vap(b,
                                             j,
                                             b._params.temperature_ref) + \
               b._params.ds_vap_ref[j]
    elif p == "Liq":
        return b._params.config.entr_mol_liq(b, j, b._params.temperature_ref)
    else:
        raise PropertyNotSupportedError(
                "{} recieved unrecognised phase name {}. Ideal property "
                "libray only supports Vap and Liq phases."
                .format(b.name, p))


def fugacity(b, p, j):
    if p == "Vap":
        return b.mole_frac_phase[p, j]*b.pressure
    elif p == "Liq":
        return b.mole_frac_phase[p, j] * \
                b._params.config.pressure_sat(b, b.temperature, j)
    else:
        raise PropertyNotSupportedError(
                "{} recieved unrecognised phase name {}. Ideal property "
                "libray only supports Vap and Liq phases."
                .format(b.name, p))


def fug_coeff(b, p, j):
    return 1
