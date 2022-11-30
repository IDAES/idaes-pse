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
import gzip
import numpy as np

import idaes
import pyomo.environ as pyo
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.core.expr.visitor import identify_variables
import pyomo.dae as pyodae
from pyomo.common import Executable
from pyomo.dae.flatten import flatten_dae_components
from pyomo.util.subsystems import (
    TemporarySubsystemManager,
    create_subsystem_block,
)
from pyomo.solvers.plugins.solvers.ASL import ASL
from pyomo.common.tempfiles import TempfileManager
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import idaes.logger as idaeslog
import idaes.config as icfg


def petsc_binary_io():
    if petsc_binary_io.PetscBinaryIOTrajectory is not None:
        return petsc_binary_io.PetscBinaryIOTrajectory

    # First see if the python IO helpers are directly importable
    try:
        import PetscBinaryIOTrajectory
        import PetscBinaryIO

        petsc_binary_io.PetscBinaryIOTrajectory = PetscBinaryIOTrajectory
        return PetscBinaryIOTrajectory
    except ImportError:
        pass
    # Next, look for a 'petscpy' directory alongside the 'petsc'
    # executable: first look at the petsc we found on the path, then
    # look in the IDAES bin dir.  Casting the Executable path to a
    # string will map None to '' in the case where there is no petsc
    # executable found.
    for petsc_exe_dir in (
        os.path.dirname(str(Executable("petsc").path())),
        icfg.bin_directory,
    ):
        if not petsc_exe_dir:
            continue
        petscpy_dir = os.path.join(petsc_exe_dir, "petscpy")
        try:
            sys.path.insert(0, petscpy_dir)
            import PetscBinaryIOTrajectory
            import PetscBinaryIO

            # Import the petsc_conf so that it is cached in sys.modules
            # (because we will be removing petscpy from sys.path)
            import petsc_conf

            petsc_binary_io.PetscBinaryIOTrajectory = PetscBinaryIOTrajectory
            return PetscBinaryIOTrajectory
        except ImportError:
            pass
        finally:
            sys.path.remove(petscpy_dir)
    return None


petsc_binary_io.PetscBinaryIOTrajectory = None


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
        super().__init__(**kwds)
        self.options.solver = "petsc"

    def _default_executable(self):
        """In addition to looking for the petsc executable, optionally check for
        a WSL batch file on Windows. Users could potentially also compile a
        cygwin exectable on Windows, so WSL isn't the only option, but it is the
        easiest for Windows."""
        executable = Executable("petsc")
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
    """PETSc solver plugin that sets options for TS solver.  This turns on
    TS monitoring, sets the DAE flag, and checks the config for default options.
    """

    def __init__(self, **kwds):
        self._ts_vars_stub = kwds.pop("vars_stub", "tmp_vars_stub")
        if "options" in kwds and kwds["options"] is not None:
            kwds["options"] = copy.deepcopy(kwds["options"])
        else:
            kwds["options"] = {}
        # Force some options.
        kwds["options"]["--dae_solve"] = ""  # is DAE solver
        kwds["options"]["--ts_monitor"] = ""  # show TS solver progress
        # We're assuming trajectory will be written in the visualization
        # style, so just set that here.
        kwds["options"]["--ts_trajectory_type"] = "visualization"
        # Get the options set in the IDAES config
        if "petsc_ts" in idaes.cfg:
            if "options" in idaes.cfg["petsc_ts"]:
                default_options = dict(idaes.cfg["petsc_ts"]["options"])
                default_options.update(kwds["options"])
                kwds["options"] = default_options
        super().__init__(**kwds)
        self.options.solver = "petsc_ts"

    def _postsolve(self):
        stub = os.path.splitext(self._soln_file)[0]
        # There is a type file created by the solver to give the variable types
        # this is needed to read the trajectory data, and we want to include it
        # with other tmp files
        typ_file = stub + ".typ"
        TempfileManager.add_tempfile(typ_file)
        # If trajectory is saved, copy the col and typ files to
        # the working directory. These files are needed to get the names and
        # types of variables and to make sense of trajectory data.
        if self.options.get("--ts_save_trajectory", 0):
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
            "The PETSc TAO interface has not yet been implemented"
        )


def petsc_available():
    """Check if the IDAES AMPL solver wrapper for PETSc is available.

    Args:
        None

    Returns (bool):
        True if PETSc is available
    """
    solver = pyo.SolverFactory("petsc")
    if solver is not None:
        try:
            return solver.available()
        except RuntimeError:
            return False
    return False


def _copy_time(time_vars, t_from, t_to):
    """PRIVATE FUNCTION:

    This is used on the flattened (indexed only by time) variable
    representations to copy variable values that are unfixed at the "to" time
    from the value at the "from" time. The PETSc DAE solver uses the initial
    variable values as the initial condition, so this is used to copy the
    previous time in as the initial condition for the next step.

    Args:
        time_vars (list): list of variables or references to variables indexed
            only by time
        t_from (float): time point to copy from
        t_to (float): time point to copy to, only unfixed vars will be
            overwritten

    Returns:
        None
    """
    for v in time_vars:
        if not v[t_to].fixed:
            v[t_to].value = v[t_from].value


def find_discretization_equations(m, time):
    """This returns a list of time discretization equations. Since we aren't
    solving the whole time period simultaneously, we'll want to deactivate
    these constraints.

    Args:
        m (Block): model or block to search for constraints
        time (ContinuousSet):

    Returns:
        list of time discretization constraints
    """
    disc_eqns = []
    for var in m.component_objects(pyo.Var):
        if isinstance(var, pyodae.DerivativeVar):
            cont_set_set = ComponentSet(var.get_continuousset_list())
            if time in ComponentSet(cont_set_set):
                if len(cont_set_set) > 1:
                    raise NotImplementedError(
                        f"IDAES presently does not support PETSc for second order or higher derivatives like {var.name} "
                        "that are differentiated at least once with respect to time. Please reformulate your model so "
                        "it does not contain such a derivative (such as by introducing intermediate variables)."
                    )
                parent = var.parent_block()
                name = var.local_name + "_disc_eq"
                disc_eq = getattr(parent, name)
                disc_eqns.append(disc_eq)
    return disc_eqns


def _set_dae_suffixes_from_variables(m, variables, deriv_diff_map):
    """Write suffixes used by the solver to identify different variable types
    and associated derivative and differential variables.

    Args:
        m: model to search for variables and write suffixes to
        variables (list): List of time indexed variables at a specific time
            point
        deriv_diff_map (ComponentMap): Maps DerivativeVar data objects to
            differential variable data objects

    Returns:
        None
    """
    # The dae_suffix provides the solver information about variables types
    # algebraic, differential, derivative, or time, see DaeVarTypes
    m.dae_suffix = pyo.Suffix(
        direction=pyo.Suffix.EXPORT,
        datatype=pyo.Suffix.INT,
    )
    # The dae_link suffix provides the solver a link between the differential
    # and derivative variable, by assigning the pair a unique integer index.
    m.dae_link = pyo.Suffix(
        direction=pyo.Suffix.EXPORT,
        datatype=pyo.Suffix.INT,
    )
    dae_var_link_index = 1
    differential_vars = []
    i = 0
    for var in variables:
        if var in deriv_diff_map:
            deriv = var
            diffvar = deriv_diff_map[deriv]
            m.dae_suffix[diffvar] = int(DaeVarTypes.DIFFERENTIAL)
            m.dae_suffix[deriv] = int(DaeVarTypes.DERIVATIVE)
            m.dae_link[diffvar] = dae_var_link_index
            m.dae_link[deriv] = dae_var_link_index
            i += 1
            dae_var_link_index += 1
            if not diffvar.fixed:
                differential_vars.append(diffvar)
            else:
                raise RuntimeError(
                    f"Problem cannot contain a fixed differential variable and "
                    f"unfixed derivative. Consider either fixing the "
                    f"corresponding derivative or adding a constraint for the "
                    f"differential variable {diffvar} possibly using an "
                    f"explicit time variable."
                )
    return differential_vars


def _get_derivative_differential_data_map(m, time):
    """Get a map from data objects of derivative variables to the
    corresponding data objects of differential variables.

    Args:
        m: Model to search for DerivativeVars
        time: Set with respect to which DerivativeVars must be differentiated

    Returns:
        (ComponentMap): Map from derivative data objects to differential
            data objects
    """
    # Get corresponding derivative and differential data objects,
    # with no attention paid to fixed or active status.
    deriv_diff_list = []
    for var in m.component_objects(pyo.Var):
        if isinstance(var, pyodae.DerivativeVar) and time in ComponentSet(
            var.get_continuousset_list()
        ):
            deriv = var
            diffvar = deriv.get_state_var()
            for idx in var:
                if deriv[idx].fixed and pyo.value(abs(deriv[idx])) > 1e-10:
                    raise RuntimeError(
                        f"{deriv[idx]} is fixed to a nonzero value "
                        f"{pyo.value(deriv[idx])}. This is "
                        f"most likely a modeling error. Instead of fixing the "
                        f"derivative consider adding a constraint like "
                        f"dxdt = constant"
                    )
                deriv_diff_list.append((deriv[idx], diffvar[idx]))

    # Get unfixed variables in active constraints
    active_con_vars = ComponentSet()
    for con in m.component_data_objects(pyo.Constraint, active=True):
        for var in identify_variables(con.expr, include_fixed=False):
            active_con_vars.add(var)

    # Filter out derivatives that are fixed or not in an active constraint
    filtered_deriv_diff_list = []
    for deriv, diff in deriv_diff_list:
        if deriv in active_con_vars:
            filtered_deriv_diff_list.append((deriv, diff))
    return ComponentMap(filtered_deriv_diff_list)


def _sub_problem_scaling_suffix(m, t_block):
    """Copy scaling factors from the full model to the submodel.  This assumes
    the scaling suffixes could be in two places.  First check the parent block
    of the component (typical place for idaes models) then check the top-level
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


class PetscDAEResults(object):
    """This class stores the results of ``petsc_dae_by_time_element()`` it has
    two attributes ``results`` and ``trajectory``.  Results is a list of Pyomo
    solver results objects.  Since this PETSc DAE solver utility solves the
    initial conditions first then may integrate over the time domain in one or
    more steps, the results of multiple solves are returned.  The trajectory is
    a ``PetscTrajectory`` object which gives the full trajectory of the solution
    for all time steps taken by the PETSc solver. This is generally finer than
    the Pyomo.DAE discretization.
    """

    def __init__(self, results=None, trajectory=None):
        self.results = results
        self.trajectory = trajectory


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
    keepfiles=False,
    symbolic_solver_labels=True,
    between=None,
    interpolate=True,
    calculate_derivatives=True,
    previous_trajectory=None,
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
            up automatically, this generally includes non-time-indexed
            constraints.
        initial_variables (list): This is a list of variables to fix after the
            initial condition solve step.  If these variables were originally
            unfixed, they will be unfixed at the end of the solve. This usually
            includes non-time-indexed variables that are calculated along with
            the initial conditions.
        detect_initial (bool): If True, add non-time-indexed variables and
            constraints to initial_variables and initial_constraints.
        skip_initial (bool): Don't do the initial condition calculation step,
            and assume that the initial condition values have already been
            calculated. This can be useful, for example, if you read initial
            conditions from a separately solved steady state problem, or
            otherwise know the initial conditions.
        snes_options (dict): PETSc nonlinear equation solver options
        ts_options (dict): PETSc time-stepping solver options
        keepfiles (bool): pass to keepfiles arg for solvers
        symbolic_solver_labels (bool): pass to symbolic_solver_labels argument
            for solvers. If you want to read trajectory data from the
            time-stepping solver, this should be True.
        between (list or tuple): List of time points to integrate between. If
            None use all time points in the model. Generally the list should
            start with the first time point and end with the last time point;
            however, this is not a requirement and there are situations where you
            may not want to include them. If you are not including the first
            time point, you should take extra care with the initial conditions.
            If the initial conditions are already correct consider using the
            ``skip_initial`` option if the first time point is not included.
        interpolate (bool): if True and trajectory is read, interpolate model
            values from the trajectory
        calculate_derivatives: (bool) if True, calculate the derivative values
            based on the values of the differential variables in the discretized
            Pyomo model.
        previous_trajectory: (PetscTrajectory) Trajectory from previous integration
            of this model. New results will be appended to this trajectory object.

    Returns (PetscDAEResults):
        See PetscDAEResults documentation for more informations.
    """
    if interpolate:
        if ts_options is None:
            ts_options = {}
        ts_options["--ts_save_trajectory"] = 1
    if timevar:
        for t in time:
            timevar[t].set_value(t)
    if between is None:
        between = time
    else:
        bad_times = between - time
        if bad_times:
            raise RuntimeError(
                "Elements of the 'between' argument must be in the time set"
            )
        between = pyo.Set(initialize=sorted(between))
        between.construct()

    solve_log = idaeslog.getSolveLogger("petsc-dae")
    regular_vars, time_vars = flatten_dae_components(m, time, pyo.Var)
    regular_cons, time_cons = flatten_dae_components(m, time, pyo.Constraint)
    tdisc = find_discretization_equations(m, time)

    solver_dae = pyo.SolverFactory("petsc_ts", options=ts_options)
    save_trajectory = solver_dae.options.get("--ts_save_trajectory", 0)

    # First calculate the inital conditions and non-time-indexed constraints
    res_list = []
    t0 = between.first()
    # list of variables to add to initial condition problem
    if initial_variables is None:
        initial_variables = []
    if detect_initial:
        rvset = ComponentSet(regular_vars)
        ivset = ComponentSet(initial_variables)
        initial_variables = list(ivset | rvset)

    if not skip_initial:
        # Nonlinear equation solver for initial conditions
        solver_snes = pyo.SolverFactory("petsc_snes", options=snes_options)
        # list of constraints to add to the initial condition problem
        if initial_constraints is None:
            initial_constraints = []

        if detect_initial:
            # If detect_initial, solve the non-time-indexed variables and
            # constraints with the initial conditions
            rcset = ComponentSet(regular_cons)
            icset = ComponentSet(initial_constraints)
            initial_constraints = list(icset | rcset)

        with TemporarySubsystemManager(to_deactivate=tdisc):
            constraints = [
                con[t0] for con in time_cons if t0 in con
            ] + initial_constraints
            variables = [var[t0] for var in time_vars] + initial_variables
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

    tprev = t0
    count = 1
    fix_derivs = []
    tj = previous_trajectory
    if tj is not None:
        variables_prev = [var[t0] for var in time_vars]

    with TemporarySubsystemManager(
        to_deactivate=tdisc,
        to_fix=initial_variables + fix_derivs,
    ):
        # Solver time steps
        deriv_diff_map = _get_derivative_differential_data_map(m, time)
        for t in between:
            if t == between.first():
                # t == between.first() was handled above
                continue
            constraints = [con[t] for con in time_cons if t in con]
            variables = [var[t] for var in time_vars]
            # Create a temporary block with references to original constraints
            # and variables so we can integrate this "subsystem" without
            # altering the rest of the model.
            t_block = create_subsystem_block(constraints, variables)
            differential_vars = _set_dae_suffixes_from_variables(
                t_block,
                variables,
                deriv_diff_map,
            )
            # We need to check if there are derivatives in the problem before
            # sending this to the solver.  We'll assume that if you are using
            # this and don't have any differential equations, you are making a
            # mistake.
            if len(differential_vars) < 1:
                raise RuntimeError(
                    f"No differential equations found at t = {t}, not a DAE"
                )
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
            if save_trajectory:
                tj_prev = tj
                tj = PetscTrajectory(
                    stub="tmp_vars_stub",
                    delete_on_read=True,
                    unscale=True,
                    model=t_block,
                )
                # add fixed vars to the trajectory. this does two things 1)
                # helps users looking for fixed var trajectory and 2) lets
                # us concatenate trajectories with section that are mixed fixed
                # and unfixed
                for i, v in enumerate(variables):
                    if isinstance(v.parent_component(), pyodae.DerivativeVar):
                        continue  # skip derivative vars
                    try:
                        vec = tj.get_vec(v)
                    except KeyError:
                        tj._set_vec(v, [pyo.value(v)] * len(tj.time))
                if tj_prev is not None:
                    # due to the way variables is generated we know variables
                    # have corresponding positions in the list
                    no_repeat = set()
                    for i, v in enumerate(variables):
                        vp = variables_prev[i]
                        if id(v) in no_repeat:
                            continue  # variables can be repeated in list
                        if isinstance(v.parent_component(), pyodae.DerivativeVar):
                            continue  # skip derivative vars
                        no_repeat.add(id(v))
                        # We'll add fixed vars in case they aren't fixed in
                        # another section. Fixed vars don't go to the solver
                        # so they don't show up in the trajectory data
                        vec = tj.get_vec(v)
                        vec_prev = tj_prev.get_vec(vp)
                        tj._set_vec(v, vec_prev + vec)
                    tj._set_time_vec(tj_prev.time + tj.time)
                variables_prev = variables
            tprev = t
            res_list.append(res)
        # If the interpolation option is True and the trajectory is available
        # interpolate the values any skipped time points from the trajectory
        if interpolate and tj is not None:
            t0 = between.first()
            tlast = between.last()
            itime = [t for t in time if not (t < t0 or t > tlast)]
            no_repeat = set()
            if timevar is not None:
                for t in itime:
                    timevar[t] = t
                no_repeat.add(id(timevar[tlast]))
            for var in time_vars:
                if id(var[tlast]) in no_repeat:
                    continue
                no_repeat.add(id(var[tlast]))
                if isinstance(var[t0].parent_component(), pyodae.DerivativeVar):
                    continue  # skip derivative vars
                vec = tj.interpolate_vec(itime, var[tlast])
                for i, t in enumerate(itime):
                    if t in between:
                        # Time is already set
                        continue
                    if not var[t].fixed:
                        # May not have trajectory from fixed variables and they
                        # shouldn't change anyway, so only set not fixed vars
                        var[t].value = vec[i]
        if calculate_derivatives:
            # the petsc solver interface does not currently return time
            # derivatives, and if it did, they would be estimated based on a
            # smaller time step.  This option uses Pyomo.DAE's discretization
            # equations to calculate the time derivative values
            calculate_time_derivatives(m, time)
        # return the solver results and trajectory if available
    return PetscDAEResults(results=res_list, trajectory=tj)


def calculate_time_derivatives(m, time):
    """Calculate the derivative values from the discretization equations.

    Args:
        m (Block): Pyomo model block
        time (ContinuousSet): Time set

    Returns:
        None
    """
    for var in m.component_objects(pyo.Var):
        if isinstance(var, pyodae.DerivativeVar):
            if time in ComponentSet(var.get_continuousset_list()):
                parent = var.parent_block()
                name = var.local_name + "_disc_eq"
                disc_eq = getattr(parent, name)
                for i, v in var.items():
                    try:
                        if disc_eq[i].active:
                            v.value = 0  # Make sure there is a value
                            calculate_variable_from_constraint(v, disc_eq[i])
                    except KeyError:
                        pass  # discretization equation may not exist at first time


class PetscTrajectory(object):
    def __init__(
        self,
        stub=None,
        vecs=None,
        json=None,
        pth=None,
        vis_dir="Visualization-data",
        delete_on_read=False,
        unscale=None,
        model=None,
        no_read=False,
    ):
        """Class to read PETSc TS solver trajectory data.  This can either read
        PETSc output by providing the ``stub`` argument, a trajectory dict by
        providing ``vecs`` or a json file by providing ``json``.

        Args:
            stub (str): file name stub for variable info
            pth (str): path where variable info and trajectory data are stored
                if None, use current working directory
            vis_dir (str): subdirectory where visualization data is stored
            delete_on_read (bool): if true delete trajectory data after reading
            unscale (bool or Block): if True and model is specified, use model
                to obtain scale factors and unscale the trajectory. If a Block
                is specified use the specified block to unscale the model. If
                False or None do not unscale.
            model (Block): if specified use for unscaling
            no_read (bool): if True make an uninitialized trajectory object
        """
        if no_read:
            return
        if petsc_binary_io() is None and stub is not None:
            raise RuntimeError("PetscBinaryIOTrajectory could not be imported")
        # if unscale is True, use model as scale factor source
        if model is not None and unscale is True:
            unscale = model
        self.model = model
        self.id_map = {}
        if pth is not None:
            stub = os.path.join(pth, stub)
            vis_dir = os.path.join(pth, vis_dir)
        if stub is not None:
            self.stub = stub
            self.vis_dir = vis_dir
            self.path = pth
            self.unscale = unscale
            self._read()
            if delete_on_read:
                self.delete_files()
            if unscale is not None:
                self._unscale(unscale)
        elif vecs is not None:
            self.vecs = vecs
            self.time = vecs["_time"]
        elif json is not None:
            self.from_json(json)
        else:
            raise RuntimeError("To read trajectory, provide stub, vecs, or json")

    def _read(self):
        with open(f"{self.stub}.col") as f:
            names = list(map(str.strip, f.readlines()))
        with open(f"{self.stub}.typ") as f:
            typ = list(map(int, f.readlines()))
        vars = [name for i, name in enumerate(names) if typ[i] in [0, 1]]
        (t, v, names) = petsc_binary_io().ReadTrajectory("Visualization-data")
        self.time = t
        self.vecs_by_time = v
        self.vecs = dict.fromkeys(vars, None)
        for k in self.vecs.keys():
            self.vecs[k] = [0] * len(self.time)
        self.vecs["_time"] = list(self.time)
        for i, v in enumerate(vars):
            for j, vt in enumerate(self.vecs_by_time):
                self.vecs[v][j] = vt[i]

    def _set_vec(self, var, vec):
        try:
            var = self.id_map[id(var)]
        except KeyError:
            var_str = str(var)
            self.id_map[id(var)] = var_str
            var = var_str
        self.vecs[var] = vec

    def _set_time_vec(self, vec):
        self.vecs["_time"] = vec
        self.time = self.vecs["_time"]

    def get_vec(self, var):
        """Return the vector of variable values at each time point for var.

        Args:
            var (str or Var): Variable to get vector for.
            time (Set): Time index set

        Retruns (list):
            vector of variable values at each time point

        """
        try:
            var = self.id_map[id(var)]
        except KeyError:
            var_str = str(var)
            self.id_map[id(var)] = var_str
            var = var_str
        return self.vecs[var]

    def get_dt(self):
        """Get a list of time steps

        Args:
            None

        Returns:
            (list)
        """
        dt = [None] * (len(self.time) - 1)
        for i in range(len(self.time) - 1):
            dt[i] = self.time[i + 1] - self.time[i]
        return dt

    def interpolate(self, times):
        """Create a new vector dictionary interpolated at times. This method
        will also extrapolate values outside the original time range, so care
        should be taken not to specify times too far outside the range.

        Args:
            times (list): list of times to interpolate. These must be in
                increasing order.

        Returns (dict):
            Dictionary of vectors for values at interpolated points
        """
        tj = PetscTrajectory(no_read=True)
        tj.id_map = copy.copy(self.id_map)
        tj.vecs = dict.fromkeys(self.vecs.keys(), None)
        tj.vecs["_time"] = copy.copy(times)
        tj.time = tj.vecs["_time"]
        for var in self.vecs:
            if var == "_time":
                continue
            tj.vecs[var] = np.interp(tj.time, self.time, self.vecs[var])
        return tj

    def interpolate_vec(self, times, var):
        """Create a new vector dictionary interpolated at times. This method
        will also extrapolate values outside the original time range, so care
        should be taken not to specify times too far outside the range.

        Args:
            times (list): list of times to interpolate. These must be in
                increasing order.

        Returns (dict):
            Dictionary of vectors for values at interpolated points
        """
        return np.interp(times, self.time, self.get_vec(var))

    def _unscale(self, m):
        """If variable scale factors are used, the solver will see scaled
        variables, and the scaled trajectory will be written. This function
        uses variable scaling factors from the given model to unscale the
        trajectory.

        Args:
            m (Block): model or block to read scale factors from.

        Returns:
            None
        """
        # Variables might show up more than once because of References
        already_scaled = set()
        for var in m.component_data_objects():
            try:
                vec = self.get_vec(var)
            except KeyError:
                continue
            vid = id(var)
            if vid in already_scaled:
                continue
            s = None
            if hasattr(var.parent_block(), "scaling_factor"):
                s = var.parent_block().scaling_factor.get(var, s)
            if hasattr(m, "scaling_factor"):
                s = m.scaling_factor.get(var, s)
            if s is not None:
                for i, x in enumerate(vec):
                    vec[i] = x / s
            already_scaled.add(vid)

    def delete_files(self):
        """Delete the trajectory data and variable information files.

        Args:
            None

        Returns:
            None
        """
        shutil.rmtree(self.vis_dir)
        os.remove(f"{self.stub}.col")
        os.remove(f"{self.stub}.typ")

    def to_json(self, pth):
        """Dump the trajectory data to a json file in the form of a dictionary
        with variable name keys and '_time' with vectors of values at each time.

        Args:
            pth (str): path for json file to write

        Returns:
            None
        """
        if pth.endswith(".gz"):
            with gzip.open(pth, "w") as fp:
                fp.write(json.dumps(self.vecs).encode("utf-8"))
        else:
            with open(pth, "w") as fp:
                json.dump(self.vecs, fp)

    def from_json(self, pth):
        """Read the trajectory data from a json file in the form of a dictionary.

        Args:
            pth (str): path for json file to write

        Returns:
            None
        """
        if pth.endswith(".gz"):
            with gzip.open(pth, "r") as fp:
                self.vecs = json.loads(fp.read())
        else:
            with open(pth, "r") as fp:
                self.vecs = json.load(fp)
        self.time = self.vecs["_time"]
