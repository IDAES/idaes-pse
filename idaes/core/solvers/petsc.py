#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################

import sys
import enum
import idaes
import pyomo.environ as pyo
from pyomo.opt.base.solvers import UnknownSolver
import pyomo.dae as pyodae
from pyomo.common import Executable
from pyomo.common.collections import ComponentSet
from pyomo.dae.flatten import flatten_dae_components
from pyomo.util.subsystems import (
    TemporarySubsystemManager,
    create_subsystem_block,
)
from pyomo.common.errors import ApplicationError
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
from idaes.core.util import get_solver


class DaeVarTypes(enum.Enum):
    """Enum DAE variable type for suffix
    """
    ALGEBRAIC = 0
    DIFFERENTIAL = 1
    DERIVATIVE = 2
    TIME = 3


class PetscSolverType(enum.Enum):
    """PETSc solver types, TAO (optimization) is not yet implimented
    """
    SNES = 0
    TS = 1
    TAO = 2


def get_petsc_solver(options=None, wsl=None, solver_type=None):
    """Get a Pyomo PETSc solver object. The IDAES solver distribution does not
    contain a PETSc executable for Windows, so the recomended method of using
    PETSc on Windows is to use the WSL to run the Linux executable.  This
    function provides a wrapper for the SovlerFactory that allows the same
    function to be used to get the PETSc solver whether using the WSL or not.
    For more information on how to set the PETSc solver up on Windows see the
    IDAES documentation.

    Args:
        options (dict): Solver options, default=None
        wsl (bool): If True force WSL version, if False force not WSL version,
            if None, try non-WSL version then try WSL version
        solver_type (PetscSolverType): If a type is provided the default options
            dictionary for the specified type will used and updated with any
            options provided.

    Returns:
        PETSc Pyomo solver object
    """
    if solver_type is not None:
        opt = petsc_default_options(solver_type=solver_type)
        if options is not None:
            opt.update(options)
        options = opt
    if not wsl or wsl is None:
        if Executable("petsc"): # checking this first avoids lot of error log
            return pyo.SolverFactory("petsc", type="asl", options=options)
    if sys.platform.startswith('win32') and (wsl or wsl is None):
        # On Windows, assume running the solver using WSL, user will need to
        # add batch file to `idaes bin-directory`
        if Executable("petsc_wsl.bat"):
            return pyo.SolverFactory(
                "petsc_wsl",
                executable="petsc_wsl.bat",
                type='asl',
                options=options)
    return None # didn't find PETSc at all


def petsc_available(wsl=None):
    """Check if the IDAES AMPL solver wrapper for PETSc is available.

    Args:
        wsl (bool): If True force WSL version, if False force not WSL version,
            if None, try non-WSL version then try WSL version

    Returns:
        (bool): True is PETSc is available
    """
    solver = get_petsc_solver(wsl=wsl)
    if solver is not None:
        return solver.available()
    return False


def petsc_default_options(solver_type=PetscSolverType.SNES):
    """Get a default options dictionary for a PETSc solver.  This provides a
    convenient way to get a resonable set of solver options for a nonlinear
    system or DAE.

    Args:
        solver_type (PetscSolverType): In {PetscSolverType.SNES,
            PetscSolverType.DAE}

    Returns:
        (dict): Solver options dictionary
    """

    if solver_type == PetscSolverType.SNES:
        return {
            "--snes_monitor":"",
            "--pc_type":"lu",             #direct solve MUMPS default LU fact
            "--ksp_type":"preonly",       #no ksp used direct solve preconditioner
        }
    elif solver_type == PetscSolverType.TS:
        return {
            "--dae_solve":"",             #tell solver to expect dae problem
            "--ts_monitor":"",            #show progess of TS solver
            "--ts_max_snes_failures":40,  #max nonlin solve fails before give up
            "--ts_max_reject":20,         #max steps to reject
            "--ts_type":"alpha",          #ts_solver
            "--snes_monitor":"",          #show progress on nonlinear solves
            "--pc_type":"lu",             #direct solve MUMPS default LU fact
            "--ksp_type":"preonly",       #no ksp used direct solve preconditioner
            #"--show_jac":"",
            #"--show_initial":"",
            "--snes_type":"newtonls",     # newton line search for nonliner solver
            "--ts_adapt_type":"basic",
            "--ts_init_time":0,           # initial time
            "--ts_max_time":1,            # final time
            #"--ts_save_trajectory":1,
            #"--ts_trajectory_type":"visualization",
            #"--ts_exact_final_time":"stepover",
            "--ts_exact_final_time":"matchstep",
            #"--ts_exact_final_time":"interpolate",
            #"--ts_view":""
        }
    elif solver_type == PetscSolverType.TAO:
        raise NotImplementedError("Tao optimization solvers not yet implimented")
    raise ValueError(f"{solver_type} is not a supported solver type")


def _copy_time(time_vars, t_from, t_to):
    for v in time_vars:
        if not v[t_to].fixed:
            v[t_to].value = v[t_from].value

def _generate_time_discretization(m, time):
    for var in m.component_objects(pyo.Var):
        if isinstance(var, pyodae.DerivativeVar):
            if time in ComponentSet(var.get_continuousset_list()):
                parent = var.parent_block()
                name = var.local_name + "_disc_eq"
                disc_eq = getattr(parent, name)
                yield disc_eq

def _set_dae_suffixes_from_variables(m, time, variables):
    """
    This is an alternative to the above method that allows us to build
    a "local" suffix, which only contains variables at a single point.

    This method is unfortunately very slow and somewhat fragile.
    There is a better way to get corresponding differential and derivative
    variables, partitioned by time, but it relies on the
    flatten_component_along_sets function added in Pyomo PR #2141

    """
    m.dae_suffix = pyo.Suffix(direction=pyo.Suffix.EXPORT, datatype=pyo.Suffix.INT)
    m.dae_link = pyo.Suffix(direction=pyo.Suffix.EXPORT, datatype=pyo.Suffix.INT)
    dae_var_link_index = 1
    differential_vars = []
    for var in variables:
        # Getting an index from a vardata is unfortunately very slow
        t = var.index()
        # this recovers the original component from the reference
        parent = var.parent_component()
        if (isinstance(parent, pyodae.DerivativeVar)
                and time in ComponentSet(parent.get_continuousset_list())):
            # This works because the original component is also indexed
            # only by time.
            difvar = parent.get_state_var()[t]
            differential_vars.append(difvar)
            m.dae_suffix[difvar] = 1
            m.dae_suffix[var] = 2
            m.dae_link[difvar] = dae_var_link_index
            m.dae_link[var] = dae_var_link_index
            dae_var_link_index += 1
    return differential_vars

def petsc_dae_by_time_element(
    m,
    time,
    initial_constraints=None,
    initial_variables=None,
    snes_options=None,
    ts_options=None,
    wsl=None,
):
    """Solve a DAE problem step by step using the PETSc DAE solver.  This
    integrates from one time point to the next.

    Args:
        m (Block): Pyomo model to solve
        time (ContinuousSet): Time set
        initial_constraints (list): Constraints to solve with in the initial
            condition solve step.  This can include initial condition constraints
            as well as other non-time-indexed constraints
        initial_variables (list): This is a list of variables to fix after the
            initial condition solve step.  If these variables were origially
            unfixed, they will be unfixed at the end of the solve.
        snes_options (dict): PETSc nonlinear equation solver options
        ts_options (dict): PETSc time-stepping solver options.
        wsl (bool): if True use WSL to run PETSc, if False don't use WSL to run
            PETSc, if None automatic. The WSL is only for Windows.

    Returns:
        None (for now, should return some status)
    """
    solve_log = idaeslog.getSolveLogger("petsc-dae")
    regular_vars, time_vars = flatten_dae_components(m, time, pyo.Var)
    regular_cons, time_cons = flatten_dae_components(m, time, pyo.Constraint)
    time_disc = list(_generate_time_discretization(m, time))

    #_set_dae_suffixes(m, time)
    solver_snes = get_petsc_solver(
        solver_type=PetscSolverType.SNES,
        options=snes_options,
        wsl=wsl)
    solver_dae = get_petsc_solver(
        solver_type=PetscSolverType.TS,
        options=ts_options,
        wsl=wsl)

    # Context manager so we don't forget to reactivate
    with TemporarySubsystemManager(to_deactivate=time_disc):
        # Solver time steps
        var_unfix = []
        for t in time:
            constraints = [con[t] for con in time_cons if t in con]
            variables = [var[t] for var in time_vars]

            # Create a temporary block with references to original
            # constraints and variables so we can integrate this "subsystem"
            # without altering the rest of the model.
            if t == time.first():
                t_block = create_subsystem_block(
                    constraints + initial_constraints,
                    variables
                )
                with idaeslog.solver_log(solve_log, idaeslog.INFO) as slc:
                    res = solver_snes.solve(t_block, tee=slc.tee)
                for var in initial_variables:
                    if not var.fixed:
                        var.fix()
                        var_unfix.append(var)
            else:
                t_block = create_subsystem_block(constraints, variables)
                differential_vars = _set_dae_suffixes_from_variables(
                    t_block, time, variables)
                _copy_time(time_vars, tprev, t)
                with idaeslog.solver_log(solve_log, idaeslog.INFO) as slc:
                    res = solver_dae.solve(
                        t_block,
                        tee=slc.tee,
                        #keepfiles=True,
                        #symbolic_solver_labels=True,
                        export_nonlinear_variables=differential_vars,
                        options={"--ts_init_time":tprev, "--ts_max_time":t}
                    )
            tprev = t
            for var in var_unfix:
                var.unfix()
