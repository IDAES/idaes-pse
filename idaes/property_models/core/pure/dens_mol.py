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


def Perry_dens_liq(b, T, j):
    # Correlation from "Perry's Chemical Engineers Handbook by Robert H. Perry"
    # 7th Edition, pg. 2-98
    return (b._params.dens_mol_liq_coeff[j, '1'] /
            b._params.dens_mol_liq_coeff[j, '2']**(
                    1 + (1-T/b._params.dens_mol_liq_coeff[j, '3']) **
                    b._params.dens_mol_liq_coeff[j, '4']))
