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
import copy
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
from pyomo.solvers.plugins.solvers.ASL import ASL
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


@pyo.SolverFactory.register('petsc', doc='ASL PETSc interface')
class Petsc(ASL):
    def __init__(self, **kwds):
        self._wsl = kwds.pop("wsl", None)
        super().__init__(**kwds)
        self.options.solver = "petsc"

    def _default_executable(self):
        executable = False
        if not self._wsl or self._wsl is None:
            executable = Executable("petsc")
        if sys.platform.startswith('win32') and (not executable):
            # On Windows, if wsl was requested or a normal petsc solver
            # executable was not found, look for a batch file to run it with WSL
            executable = Executable("petsc_wsl.bat")
        if not executable:
            raise RuntimeError("No PETSc executable found.")
        return executable.path()

@pyo.SolverFactory.register('petsc_snes', doc='ASL PETSc SNES interface')
class PetscSNES(Petsc):
    def __init__(self, **kwds):
        if "options" in kwds and kwds["options"] is not None:
            kwds["options"] = copy.deepcopy(kwds["options"])
        else:
            kwds["options"] = {}
        kwds["options"]["--snes_monitor"] = ""
        if "petsc_snes" in idaes.cfg:
            if "options" in idaes.cfg["petsc_snes"]:
                default_options = dict(idaes.cfg["petsc_snes"]["options"])
                default_options.update(kwds["options"])
                kwds["options"] = default_options
        super().__init__(**kwds)
        self.options.solver = "petsc_snes"

@pyo.SolverFactory.register('petsc_ts', doc='ASL PETSc TS interface')
class PetscTS(Petsc):
    def __init__(self, **kwds):
        if "options" in kwds and kwds["options"] is not None:
            kwds["options"] = copy.deepcopy(kwds["options"])
        else:
            kwds["options"] = {}
        kwds["options"]["--dae_solve"] = ""
        kwds["options"]["--ts_monitor"] = ""
        if "petsc_ts" in idaes.cfg:
            if "options" in idaes.cfg["petsc_ts"]:
                default_options = dict(idaes.cfg["petsc_ts"]["options"])
                default_options.update(kwds["options"])
                kwds["options"] = default_options
        super().__init__(**kwds)
        self.options.solver = "petsc_ts"

@pyo.SolverFactory.register('petsc_tao', doc='ASL PETSc TAO interface')
class PetscTAO(Petsc):
    def __init__(self, **kwds):
        raise NotImplementedError(
            "The PETSc TAO interface has not yet been implimented")



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
    return pyo.SolverFactory("petsc", wsl=wsl, options=options)



def petsc_available(wsl=None):
    """Check if the IDAES AMPL solver wrapper for PETSc is available.

    Args:
        wsl (bool): If True force WSL version, if False force not WSL version,
            if None, try non-WSL version then try WSL version

    Returns:
        (bool): True is PETSc is available
    """
    solver = pyo.SolverFactory("petsc", wsl=wsl)
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
            "--ts_type":"beuler",         #ts_solver
            #"--ts_dt":0.1,
            #"--ts_adapt_monitor":"",
            #"--snes_monitor":"",          #show progress on nonlinear solves
            #"--pc_type":"lu",             #direct solve MUMPS default LU fact
            #"--ksp_type":"preonly",       #no ksp used direct solve preconditioner
            #"--show_jac":"",
            #"--show_initial":"",
            #"--snes_type":"newtonls",     # newton line search for nonliner solver
            "--ts_adapt_type":"basic",
            #"--ts_init_time":0,           # initial time
            #"--ts_max_time":1,            # final time
            #"--ts_save_trajectory":1,
            #"--ts_trajectory_type":"visualization",
            #"--ts_exact_final_time":"stepover",
            "--ts_exact_final_time":"matchstep",
            #"--ts_exact_final_time":"interpolate",
            #"--ts_view":"",
        }
    elif solver_type == PetscSolverType.TAO:
        raise NotImplementedError("Tao optimization solvers not yet implimented")
    raise ValueError(f"{solver_type} is not a supported solver type")


def _copy_time(time_vars, t_from, t_to):
    """PRIVATE FUNCTION:

    This is used on the flattened (only indexed by time) variable
    representations to copy variable values that are unfixed at the "to" time
    from the value at the "from" time. The PETSc DAE solver uses the initial
    variable values as the initial condition, so this is used to copy the
    previous time in as the initial condition for the next step.

    Args:
        time_vars (list): list of variables or refernces to variables indexed
            only by time
        t_from (float): time point to copy from
        t_to (float): time point to copy to, only unfix vars will be overwritten

    Returns:
        None
    """
    for v in time_vars:
        if not v[t_to].fixed:
            v[t_to].value = v[t_from].value

def _generate_time_discretization(m, time):
    """This is a generator for time discretization equations. Since we aren't
    solving the whole time period simulataneously, we'll want to deactivate
    these constraints.

    Args:
        m (Block): model or block to search for constraints
        time (ContinuousSet):

    Yields:
        time discretization constraints
    """
    for var in m.component_objects(pyo.Var):
        if isinstance(var, pyodae.DerivativeVar):
            if time in ComponentSet(var.get_continuousset_list()):
                parent = var.parent_block()
                name = var.local_name + "_disc_eq"
                disc_eq = getattr(parent, name)
                yield disc_eq

def _set_dae_suffixes_from_variables(m, variables):
    """Write suffixes used by the solver to identify different variable types
    and assocciate derivative and differential variables.

    Note:  This is a method that allows us to build a "local" suffix, which only
    contains variables at a single point. This method is unfortunately very
    slow and somewhat fragile. There is a better way to get corresponding
    differential and derivative variables, partitioned by time, but it relies on
    the flatten_component_along_sets function added in Pyomo PR #2141

    Args:
        m: model search for varaibles and write suffixes to
        variables (list): List of time indexed variables at a sepcific time point

    Returns:
        None
    """
    m.dae_suffix = pyo.Suffix(direction=pyo.Suffix.EXPORT, datatype=pyo.Suffix.INT)
    m.dae_link = pyo.Suffix(direction=pyo.Suffix.EXPORT, datatype=pyo.Suffix.INT)
    dae_var_link_index = 1
    differential_vars = []
    n = 0
    for var in variables:
        # Getting an index from a vardata is unfortunately very slow
        t = var.index()
        # this recovers the original component from the reference
        parent = var.parent_component()
        if isinstance(parent, pyodae.DerivativeVar):
            n += 1
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
    timevar=None,
    initial_constraints=None,
    initial_variables=None,
    skip_initial=False,
    snes_options=None,
    ts_options=None,
    wsl=None,
    keepfiles=False,
    symbolic_solver_labels=False,
):
    """Solve a DAE problem step by step using the PETSc DAE solver.  This
    integrates from one time point to the next.

    Args:
        m (Block): Pyomo model to solve
        time (ContinuousSet): Time set
        timevar (Var): Optional sepcification of a time variable
        initial_constraints (list): Constraints to solve with in the initial
            condition solve step.  Since the time indexed constraints are picked
            up automaticly, this generally inlcudes non-time-inded constraints.
        initial_variables (list): This is a list of variables to fix after the
            initial condition solve step.  If these variables were origially
            unfixed, they will be unfixed at the end of the solve. This usually
            includes non-time-indexed varaibles that are calculated allong with
            the initial condition calculations.
        skip_initial (bool): Don't do the initial condition calculation step, and
            assume that the initial condition values have already been calculated.
            This can be usful, for example, if you read initial conditions from a
            speratly solved steady state problem, or otherwise have a know initial
            condition.
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
    tdisc = list(_generate_time_discretization(m, time))

    #_set_dae_suffixes(m, time)
    solver_snes = pyo.SolverFactory("petsc_snes", options=snes_options, wsl=wsl)
    solver_dae = pyo.SolverFactory("petsc_ts", options=ts_options, wsl=wsl)

    # First calculate the inital conditions and non-time-tindexed constraints
    t = time.first()
    if not skip_initial:
        with TemporarySubsystemManager(to_deactivate=tdisc, to_fix=initial_variables):
            constraints = [con[t] for con in time_cons if t in con]
            variables = [var[t] for var in time_vars]
            t_block = create_subsystem_block(
                constraints + initial_constraints,
                variables
            )
            with idaeslog.solver_log(solve_log, idaeslog.INFO) as slc:
                res = solver_snes.solve(t_block, tee=slc.tee)
    tprev = t
    with TemporarySubsystemManager(to_deactivate=tdisc):
        # Solver time steps
        for t in time:
            constraints = [con[t] for con in time_cons if t in con]
            variables = [var[t] for var in time_vars]
            if t == time.first():
                continue
            # Create a temporary block with references to original constraints
            # and variables so we can integrate this "subsystem" without altering
            # the rest of the model.
            t_block = create_subsystem_block(constraints, variables)
            differential_vars = _set_dae_suffixes_from_variables(t_block, variables)
            if timevar is not None:
                t_block.dae_suffix[timevar[t]] = 3
            t_block.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
            if hasattr(m, "scaling_factor"):
                for c in t_block.component_data_objects((pyo.Var, pyo.Constraint)):
                    if c in m.scaling_factor:
                        t_block.scaling_factor[c] = m.scaling_factor[c]
            # Take initial conditions for this step from the result of previous
            _copy_time(time_vars, tprev, t)
            with idaeslog.solver_log(solve_log, idaeslog.INFO) as slc:
                res = solver_dae.solve(
                    t_block,
                    tee=slc.tee,
                    keepfiles=keepfiles,
                    symbolic_solver_labels=symbolic_solver_labels,
                    export_nonlinear_variables=differential_vars,
                    options={"--ts_init_time":tprev, "--ts_max_time":t}
                )
            tprev = t
