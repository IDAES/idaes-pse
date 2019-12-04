
import sys

def solveMyCPLEXSolnPool(opt,tee=True,logfile=None):
    """

    Args:
        opt: param tee:  (Default value = True)
        logfile: Default value = None)
        tee:  (Default value = True)

    Returns:

    """
    for pyomo_key,val in opt.options.items():
        # In this loop below, we build up the attribute path to the CPLEX option
        CPLEXOptionAttr = opt._solver_model.parameters
        for pyomo_key_part in pyomo_key.split('_'):
            CPLEXOptionAttr = getattr(CPLEXOptionAttr,pyomo_key_part)
        CPLEXOptionAttr.set(val)
    if(tee):
        def _process_stream(arg):
            """

            Args:
                arg: 

            Returns:

            """
            sys.stdout.write(arg)
            return arg
        if(logfile is None):
            # NOTE: This case is necessary to correctly print output if logfile is None
            #       Else, CPLEX disables all output
            opt._solver_model.set_results_stream(sys.stdout)
        else:
            opt._solver_model.set_results_stream(logfile,_process_stream)
    else:
        opt._solver_model.set_results_stream(logfile)
    opt._solver_model.set_problem_type(opt._solver_model.problem_type.MILP)
    opt._solver_model.populate_solution_pool()

def getMyCPLEXSolnPoolSize(opt):
    """

    Args:
        opt: 

    Returns:

    """
    return opt._solver_model.solution.pool.get_num()

def getMyCPLEXSolnPoolVal(opt,v,iSoln):
    """

    Args:
        opt: param v:
        iSoln: 
        v: 

    Returns:

    """
    var_map = opt._pyomo_var_to_solver_var_map
    cplex_var_to_load = var_map[v]
    val = opt._solver_model.solution.pool.get_values(iSoln,cplex_var_to_load)
    return val

def loadMyCPLEXSolnPoolVals(opt,iSoln):
    """

    Args:
        opt: param iSoln:
        iSoln: 

    Returns:

    """
    var_map = opt._pyomo_var_to_solver_var_map
    ref_vars = opt._referenced_variables
    vars_to_load = var_map.keys()
    cplex_vars_to_load = [var_map[pyomo_var] for pyomo_var in vars_to_load]
    vals = opt._solver_model.solution.pool.get_values(iSoln,cplex_vars_to_load)
    for i, pyomo_var in enumerate(vars_to_load):
        if ref_vars[pyomo_var] > 0:
            pyomo_var.stale = False
            pyomo_var.value = vals[i]

def printMyCPLEXSolnStatus(opt,stream=sys.stdout):
    """

    Args:
        opt: param stream:  (Default value = sys.stdout)
        stream:  (Default value = sys.stdout)

    Returns:

    """
    cplex = opt._solver_model
    stream.write('Status={}\n'.format(cplex.solution.get_status_string()))
    stream.write('Objective={:.6f}\n'.format(cplex.solution.get_objective_value()))
    stream.write('BestObjective={:.6f}\n'.format(cplex.solution.MIP.get_best_objective()))
    stream.write('RelGap={:.2f}%\n'.format(100.0*cplex.solution.MIP.get_mip_relative_gap()))
