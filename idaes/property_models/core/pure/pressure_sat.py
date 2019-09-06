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
from pyomo.environ import exp


def antoine(b, T, j):
    return 10**(b._params.antoine_coeff[j, 'A'] -
                b._params.antoine_coeff[j, 'B'] /
                (T + b._params.antoine_coeff[j, 'C']))


def RPP_Psat1(b, T, j):
    x = 1 - T/b._params.temperature_crit[j]

    return (exp((1-x)**-1 * (b._params.pressure_sat_coeff[j, 'A']*x +
                             b._params.pressure_sat_coeff[j, 'B']*x**1.5 +
                             b._params.pressure_sat_coeff[j, 'C']*x**3 +
                             b._params.pressure_sat_coeff[j, 'D']*x**6)) *
            b._params.pressure_crit[j])
