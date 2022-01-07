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

import os
import sys
import shutil
import enum
import copy
import json
import idaes
import pyomo.environ as pyo
from pyomo.common.collections import ComponentSet
import pyomo.dae as pyodae
from pyomo.common import Executable
from pyomo.common.collections import ComponentSet
from pyomo.dae.flatten import flatten_dae_components
from pyomo.util.subsystems import (
    TemporarySubsystemManager,
    create_subsystem_block,
)
from pyomo.solvers.plugins.solvers.ASL import ASL
from pyomo.opt.solver import SystemCallSolver
from pyomo.common.tempfiles import TempfileManager
from pyomo.common.errors import ApplicationError
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
from idaes.core.util import get_solver
import idaes.config as icfg

PetscBinaryIOTrajectory = None
PetscBinaryIO = None

def _import_petsc_binary_io():
    global PetscBinaryIOTrajectory
    global PetscBinaryIO
    try:
        import PetscBinaryIOTrajectory
        import PetscBinaryIO
    except ImportError:
        petsc_dir = os.path.join(icfg.bin_directory, "petscpy")
        if not os.path.isdir(petsc_dir):
            return
        sys.path.append(petsc_dir)
        try:
            import PetscBinaryIOTrajectory
            import PetscBinaryIO
        except ImportError:
            pass

_import_petsc_binary_io()

class DaeVarTypes(enum.IntEnum):
    """Enum DAE variable type for suffix"""

    ALGEBRAIC = 0
    DIFFERENTIAL = 1
    DERIVATIVE = 2
    TIME = 3


@pyo.SolverFactory.register("petsc", doc="ASL PETSc interface")
class Petsc(ASL):
    """ASL solver plugin for the PETSc solver.  This adds the option to use an
    alternative executable batch file to run the solver through the WSL."""

    def __init__(self, **kwds):
        self._wsl = kwds.pop("wsl", None)
        super().__init__(**kwds)
        self.options.solver = "petsc"

    def _default_executable(self):
        """In addition to looking for the petsc executable, optionally check for
        a WSL batch file on Windows. Users could potentially also compile a
        cygwin exectable on Windows, so WSL isn't the only option, but it is the
        easiest for Windows."""
        executable = False
        if not self._wsl or self._wsl is None:
            executable = Executable("petsc")
        if sys.platform.startswith("win32") and (not executable):
            # On Windows, if wsl was requested or a normal petsc solver
            # executable was not found, look for a batch file to run it with WSL
            executable = Executable("petsc_wsl.bat")
        if not executable:
            raise RuntimeError("No PETSc executable found.")
        return executable.path()


@pyo.SolverFactory.register("petsc_snes", doc="ASL PETSc SNES interface")
class PetscSNES(Petsc):
    """PETSc solver plugin that sets options for SNES solver.  This turns on
    SNES monitoring, and checks the config for default options"""

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


@pyo.SolverFactory.register("petsc_ts", doc="ASL PETSc TS interface")
class PetscTS(Petsc):
    """PETSc solver plugin that sets options for SNES solver.  This turns on
    TS monitoring, sets the DAE flag, and checks the config for default options."""

    def __init__(self, **kwds):
        self._ts_vars_stub = kwds.pop("vars_stub", None)
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

    def _postsolve(self):
        stub = os.path.splitext(self._soln_file)[0]
        # There is a type file created by the solver to give the varaible types
        # this is needed to read the trajectory data, and we want to delete it
        # with other tmp files
        typ_file = stub + ".typ"
        TempfileManager.add_tempfile(typ_file)
        # If the vars_stub option was specified, copy the col and typ files to
        # the working directory. These files are need to get the names and
        # types of variables, and are need to make sense of trajectory data.
        if self._ts_vars_stub is not None:
            try:
                shutil.copyfile(f"{stub}.col", f"{self._ts_vars_stub}.col")
            except:
                pass
            try:
                shutil.copyfile(f"{stub}.typ", f"{self._ts_vars_stub}.typ")
            except:
                pass
        return ASL._postsolve(self)


@pyo.SolverFactory.register("petsc_tao", doc="ASL PETSc TAO interface")
class PetscTAO(Petsc):
    """This is a place holder for optimization solvers"""

    def __init__(self, **kwds):
        raise NotImplementedError(
            "The PETSc TAO interface has not yet been implimented"
        )


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
        try:
            return solver.available()
        except RuntimeError:
            return False
    return False


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
        t_to (float): time point to copy to, only unfixed vars will be overwritten

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
        m: model to search for variables and write suffixes to
        variables (list): List of time indexed variables at a specific time point

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
            m.dae_suffix[difvar] = int(DaeVarTypes.DIFFERENTIAL)
            m.dae_suffix[var] = int(DaeVarTypes.DERIVATIVE)
            m.dae_link[difvar] = dae_var_link_index
            m.dae_link[var] = dae_var_link_index
            dae_var_link_index += 1
    return differential_vars


def _sub_problem_scaling_suffix(m, t_block):
    """Copy scaling factors from the full model to the submodel.  This assumes
    the scaling suffixes could be in two places.  First check the parent block
    of the compoent (typical place for idaes models) then check the top-level
    model.  The top level model will take precedence.
    """
    if not hasattr(t_block, "scaling_factor"):
        # if the subsystem block doesn't already have a scaling suffix, make one
        t_block.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    # first check the parent block for scaling factors
    for c in t_block.component_data_objects((pyo.Var, pyo.Constraint)):
        if hasattr(c.parent_block(), "scaling_factor"):
            if c in c.parent_block().scaling_factor:
                t_block.scaling_factor[c] = c.parent_block().scaling_factor[c]
    # now check the top level model
    if hasattr(m, "scaling_factor"):
        for c in t_block.component_data_objects((pyo.Var, pyo.Constraint)):
            if c in m.scaling_factor:
                t_block.scaling_factor[c] = m.scaling_factor[c]


def petsc_dae_by_time_element(
    m,
    time,
    timevar=None,
    initial_constraints=None,
    initial_variables=None,
    detect_initial=True,
    skip_initial=False,
    snes_options=None,
    ts_options=None,
    wsl=None,
    keepfiles=False,
    symbolic_solver_labels=True,
    vars_stub=None,
):
    """Solve a DAE problem step by step using the PETSc DAE solver.  This
    integrates from one time point to the next.

    Args:
        m (Block): Pyomo model to solve
        time (ContinuousSet): Time set
        timevar (Var): Optional specification of a time variable, which can be
            used to write constraints that are an explicit function of time.
        initial_constraints (list): Constraints to solve in the initial
            condition solve step.  Since the time-indexed constraints are picked
            up automatically, this generally includes non-time-indexed constraints.
        initial_variables (list): This is a list of variables to fix after the
            initial condition solve step.  If these variables were originally
            unfixed, they will be unfixed at the end of the solve. This usually
            includes non-time-indexed variables that are calculated along with
            the initial conditions.
        detect_initial (bool): If True, add non-time-indexed variables and
            constraints to initial_variables and initial_constraints.
        skip_initial (bool): Don't do the initial condition calculation step, and
            assume that the initial condition values have already been calculated.
            This can be useful, for example, if you read initial conditions from a
            separately solved steady state problem, or otherwise know the initial
            conditions.
        snes_options (dict): PETSc nonlinear equation solver options
        ts_options (dict): PETSc time-stepping solver options
        wsl (bool): if True use WSL to run PETSc, if False don't use WSL to run
            PETSc, if None automatic. The WSL is only for Windows.
        keepfiles (bool): pass to keepfiles arg for solvers
        symbolic_solver_labels (bool): pass to symoblic_solver_labels argument for
            solvers. If you want to read trajectory data from the time-stepping
            solver, this should be True.
        vars_stub (str or None): Copy the `*.col` and `*.typ` files to the working
            directory using this stub if not None.  These are needed to
            interpret the trajectory data.

    Returns:
        List of solver results objects from each solve. If there are initial
        condition constraints and they are not skipped, the fist object will
        be from the initial condition solve.  Then there should be one for each
        time elemenet for each TS solve.
    """
    solve_log = idaeslog.getSolveLogger("petsc-dae")
    regular_vars, time_vars = flatten_dae_components(m, time, pyo.Var)
    regular_cons, time_cons = flatten_dae_components(m, time, pyo.Constraint)
    tdisc = list(_generate_time_discretization(m, time))

    # _set_dae_suffixes(m, time)
    solver_snes = pyo.SolverFactory("petsc_snes", options=snes_options, wsl=wsl)
    solver_dae = pyo.SolverFactory(
        "petsc_ts", options=ts_options, wsl=wsl, vars_stub=vars_stub
    )

    if initial_variables is None:
        initial_variables = []
    if initial_constraints is None:
        initial_constraints = []

    if detect_initial:
        rvset = ComponentSet(regular_vars)
        rcset = ComponentSet(regular_cons)
        icset = ComponentSet(initial_constraints)
        ivset = ComponentSet(initial_variables)
        initial_variables = list(ivset | rvset)
        initial_constraints = list(icset | rcset)

    # First calculate the inital conditions and non-time-tindexed constraints
    res_list = []
    t = time.first()
    if not skip_initial:
        with TemporarySubsystemManager(to_deactivate=tdisc):
            constraints = [
                con[t] for con in time_cons if t in con
            ] + initial_constraints
            variables = [var[t] for var in time_vars] + initial_variables
            if len(constraints) > 0:
                # if the initial condition is specified and there are no
                # initial constraints, don't try to solve.
                t_block = create_subsystem_block(
                    constraints,
                    variables,
                )
                # set up the scaling factor suffix
                _sub_problem_scaling_suffix(m, t_block)
                with idaeslog.solver_log(solve_log, idaeslog.INFO) as slc:
                    res = solver_snes.solve(t_block, tee=slc.tee)
    res_list.append(res)
    tprev = t
    with TemporarySubsystemManager(to_deactivate=tdisc, to_fix=initial_variables):
        # Solver time steps
        for t in time:
            if t == time.first():
                continue
            constraints = [con[t] for con in time_cons if t in con]
            variables = [var[t] for var in time_vars]
            # Create a temporary block with references to original constraints
            # and variables so we can integrate this "subsystem" without altering
            # the rest of the model.
            t_block = create_subsystem_block(constraints, variables)
            differential_vars = _set_dae_suffixes_from_variables(t_block, variables)
            if timevar is not None:
                t_block.dae_suffix[timevar[t]] = int(DaeVarTypes.TIME)
            # Set up the scaling factor suffix
            _sub_problem_scaling_suffix(m, t_block)
            # Take initial conditions for this step from the result of previous
            _copy_time(time_vars, tprev, t)
            with idaeslog.solver_log(solve_log, idaeslog.INFO) as slc:
                res = solver_dae.solve(
                    t_block,
                    tee=slc.tee,
                    keepfiles=keepfiles,
                    symbolic_solver_labels=symbolic_solver_labels,
                    export_nonlinear_variables=differential_vars,
                    options={"--ts_init_time": tprev, "--ts_max_time": t},
                )
            tprev = t
            res_list.append(res)
    return res_list


class PetscTrajectory(object):
    def __init__(self, stub, pth=None, vis_dir="Visualization-data", delete_on_read=False):
        """Class to read PETSc TS solver trajectory data.

        Args:
            stub (str): file name stub for variable info
            pth (str): path where variable info and trajectory data are stored
                if None, use current working directory
            vis_dir (str): subdirectory where visualization data is stored
            delete_on_read (bool): if true delete trajectory data after reading
        """
        if PetscBinaryIOTrajectory is None:
            raise RuntimeError("PetscBinaryIOTrajectory could not be imported")
        if pth is not None:
            stub = os.path.join(pth, stub)
            vis_dir = os.path.join(pth, vis_dir)
        self.stub = stub
        self.vis_dir = vis_dir
        self.path = pth
        self._read()
        if delete_on_read:
            self.delete_files()

    def _read(self):
        with open(f'{self.stub}.col') as f:
            names = list(map(str.strip, f.readlines()))
        with open(f'{self.stub}.typ') as f:
            typ = list(map(int,f.readlines()))
        self.vars = [name for i, name in enumerate(names) if typ[i] in [0,1]]
        (t, v, names) = PetscBinaryIOTrajectory.ReadTrajectory("Visualization-data")
        self.time = t
        self.vecs_by_time = v
        self.vecs = dict.fromkeys(self.vars, None)
        for k in self.vecs.keys():
            self.vecs[k] = [0]*len(self.time)
        self.vecs["_time"] = list(self.time)
        for i, v in enumerate(self.vars):
            for j, vt in enumerate(self.vecs_by_time):
                self.vecs[v][j] = vt[i]

    def _var_name_list(self, vars):
        if isinstance(vars, str):
            vars = [vars]
        if hasattr(vars, "ctype"):
            if issubclass(vars.ctype, pyo.Var):
                vars = [vars]
        vars = list(map(str, vars))
        for name in vars:
            if name not in self.vars:
                raise RuntimeError(f"Variable {name} not found.")
        return vars

    def get_vec(self, var):
        """Return the time vector for var.

        Args:
            var (str or Var): Variable to get vector for.

        Retruns:
            (list): time vector

        """
        var = str(var)
        return self.vecs[var]

    def get_dt(self):
        """Get a list of time steps

        Args:
            None

        Returns:
            (list)
        """
        dt = [None]*(len(self.time) - 1)
        for i in range(len(self.time)-1):
            dt[i] = self.time[i + 1] - self.time[i]
        return dt

    def interpolate_vecs(self, times):
        """Create a new vector dictionary interpolated at times.

        Args:
            vars (list): list of variables or variable names to plot.  The
                varibles include the time index of the final time.

        Returns:
            (dict)
        """
        vecs = dict.fromkeys(self.vars, None)
        vecs["_time"] = copy.copy(times)
        for k in vecs.keys():
            vecs[k] = [0]*len(times)
        j = 0
        for i, t in enumerate(times):
            while t > self.time[j] and j < len(self.time) - 1:
                j += 1
            if j == 0:
                j = 1
            dt = self.time[j] - self.time[j - 1]
            w1 = 1 - (t - self.time[j - 1])/dt
            w2 = (t - self.time[j - 1])/dt
            for k, v in self.vecs.items():
                vecs[k][i] = w1*v[j - 1] + w2*v[j]
        return vecs

    def plot_trajectory_original(self, vars):
        """Plot function privided by PETSc.

        Args:
            vars (list): list of variables or variable names to plot.  The
                varibles include the time index of the final time.

        Returns:
            None
        """
        vars = self._var_name_list(vars)
        PetscBinaryIOTrajectory.PlotTrajectories(
            self.time, self.vecs_by_time, self.vars, self._var_name_list(vars)
        )

    def delete_files(self):
        """Delete the trajectory data and varible information files.

        Args:
            None

        Returns:
            None
        """
        shutil.rmtree(self.vis_dir)
        os.remove(f'{self.stub}.col')
        os.remove(f'{self.stub}.typ')

    def to_json(self, pth):
        """Dump the trajectory data to a json file in the form of a dictionary
        with varaible name keys and '_time' with time vectors values.

        Args:
            pth (str): path for json file to write

        Returns:
            None
        """
        with open(pth, "w") as fp:
            json.dump(self.vecs, fp)
