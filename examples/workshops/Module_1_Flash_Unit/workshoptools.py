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

from pyomo.environ import value, TerminationCondition, SolverStatus
from pyomo.core.base.var import IndexedVar

def solve_successful(status):
    if (status.solver.termination_condition == TerminationCondition.optimal or
        status.solver.status == SolverStatus.ok):
        return True
    return False

def print_ports_summary(ports):
    # find the list of variables to print
    vars_idxs = []
    varname_set = set()
    for prt in ports:
        for vname in prt.vars:
            if vname in varname_set:
                continue

            varname_set.add(vname)
            var = prt.vars[vname]

            for k in var.keys():
                vars_idxs.append((vname, k, '{}[{}]'.format(vname, k)))
    vname_len = max([len(vidx[2]) for vidx in vars_idxs])
    vname_len = max(vname_len, len('Variable'))

    # print the header
    print('{0:>{1}s}'.format('Variable', vname_len), end="\t")

    for prt in ports:
        print(prt.name, end='\t')
    print()

    print(''.join(['-']*vname_len), end='\t')
    for prt in ports:
        print(''.join(['-']*len(prt.name)), end='\t')
    print()

    for vname, idx, vname_idx in vars_idxs:
        print('{0:>{1}s}'.format(vname_idx, vname_len), end="\t")
        for prt in ports:
            print('{0:{1}f}'.format(value(prt.vars[vname][idx]), len(prt.name)), end='\t')
        print()
