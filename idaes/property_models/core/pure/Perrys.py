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
Methods for calculating pure component properties from:

Perry's Chemical Engineers' Handbook, 7th Edition
Perry, Green, Maloney, 1997, McGraw-Hill

All parameter indicies based on conventions used by the source
"""
from pyomo.environ import log


# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
def Perry_cp_liq(b, j, T):
    # Specific enthalpy
    return (b._params.cp_liq[j, "5"]*T**4 +
            b._params.cp_liq[j, "4"]*T**3 +
            b._params.cp_liq[j, "3"]*T**2 +
            b._params.cp_liq[j, "2"]*T +
            b._params.cp_liq[j, "1"])


def Perry_enth_liq(b, j, T):
    # Specific enthalpy
    return ((b._params.cp_liq[j, "5"]/5) *
            (T**5-b._params.temperature_ref**5) +
            (b._params.cp_liq[j, "4"]/4) *
            (T**4-b._params.temperature_ref**4) +
            (b._params.cp_liq[j, "3"]/3) *
            (T**3-b._params.temperature_ref**3) +
            (b._params.cp_liq[j, "2"]/2) *
            (T**2-b._params.temperature_ref**2) +
            b._params.cp_liq[j, "1"] *
            (T-b._params.temperature_ref))


def Perry_entr_liq(b, j, T):
    # Specific entropy
    return ((b._params.cp_liq[j, '5']/4)*T**4 +
            (b._params.cp_liq[j, '4']/3)*T**3 +
            (b._params.cp_liq[j, '3']/2)*T**2 +
            b._params.cp_liq[j, '2']*T +
            b._params.cp_liq[j, '1']*log(T))


# -----------------------------------------------------------------------------
# Densities
def Perry_dens_liq(b, T, j):
    # pg. 2-98
    return (b._params.dens_mol_liq_coeff[j, '1'] /
            b._params.dens_mol_liq_coeff[j, '2']**(
                    1 + (1-T/b._params.dens_mol_liq_coeff[j, '3']) **
                    b._params.dens_mol_liq_coeff[j, '4']))
