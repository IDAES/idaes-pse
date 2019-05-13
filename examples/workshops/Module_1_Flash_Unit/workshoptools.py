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

    
