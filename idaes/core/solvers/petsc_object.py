import os
import sys
import shutil
import enum
import copy
import json
import gzip
import numpy as np

from pyomo.core.base import BlockData
from pyomo.environ import Constraint, Set, SolverFactory, Suffix, value, Var
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.core.expr.visitor import identify_variables
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common import Executable
from pyomo.common.config import ConfigDict, ConfigValue, document_kwargs_from_configdict, ListOf
from pyomo.dae.flatten import flatten_dae_components, slice_component_along_sets
from pyomo.util.subsystems import (
    TemporarySubsystemManager,
    create_subsystem_block,
)
from pyomo.solvers.plugins.solvers.ASL import ASL
from pyomo.common.tempfiles import TempfileManager
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.deprecation import deprecation_warning

import idaes
import idaes.logger as idaeslog
import idaes.config as icfg

from idaes.core.scaling import get_scaling_factor, set_scaling_factor
from idaes.core.solvers.petsc import (
    calculate_time_derivatives,
    DaeVarTypes,
    find_discretization_equations,
    _get_derivative_differential_data_map,
    PetscDAEResults,
    PetscTrajectory,
    _set_dae_suffixes_from_variables,
)

from idaes.core.util.model_serializer import StoreSpec, from_json, to_json

# Importing a few things here so that they are cached
# pylint: disable=unused-import
# pylint: disable=import-outside-toplevel
# pylint: disable=protected-access

DAE_DISC_SUFFIX = "_disc_eq"

CONFIG = ConfigDict()

CONFIG.declare(
    "time_var",
    ConfigValue(
        default=None,
        domain=Var,
        description="Optional specification of a time variable, which can be"
            "used to write constraints that are an explicit function of time.",
    ),
)
# CONFIG.declare(
#     "initial_constraints",
#     ConfigValue(
#         default=None,
#         domain=ListOf(Constraint),
#         description="Constraints to solve in the initial "
#             "condition solve step.  Since the time-indexed constraints are picked "
#             "up automatically, this generally includes non-time-indexed "
#             "constraints.",
#     )
# )
# CONFIG.declare(
#     "initial_variables",
#     ConfigValue(
#         default=None,
#         domain=ListOf(Var),
#         description="initial_variables (list): This is a list of variables to fix after the "
#             "initial condition solve step.  If these variables were originally "
#             "unfixed, they will be unfixed at the end of the solve. This usually "
#             "includes non-time-indexed variables that are calculated along with "
#             "the initial conditions."
#     )
# )
# CONFIG.declare(
#     "detect_initial",
#     ConfigValue(
#         default=True,
#         domain=bool,
#         description="If True, add non-time-indexed variables and "
#             "constraints to initial_variables and initial_constraints."
#     )
# )
# CONFIG.declare(
#     "skip_initial",
#     ConfigValue(
#         default=False,
#         domain=bool,
#         description="Don't do the initial condition calculation step, "
#             "and assume that the initial condition values have already been "
#             "calculated. This can be useful, for example, if you read initial "
#             "conditions from a separately solved steady state problem, or "
#             "otherwise know the initial conditions."
#     )
# )
CONFIG.declare(
    "initial_solver",
    ConfigValue(
        default="petsc_snes",
        domain=str,
        description="The nonlinear equations solver "
            "to use for the initial conditions (e.g. petsc_snes, ipopt, ...).",
    )
)
CONFIG.declare(
    "initial_solver_options",
    ConfigValue(
        default=None,
        domain=dict,
        description="Solver options to use with initial_solver."
    )
)
CONFIG.declare(
    "initial_solver_writer_config",
    ConfigValue(
        default=None,
        domain=dict,
        description="Configurations to use with the .nl writer "
            "in the initial condition problem."
    )
)
CONFIG.declare(
    "ts_options",
    ConfigValue(
        default=None,
        domain=dict,
        description="Options to use with the PETSc integrator TS."
    )
)
# CONFIG.declare(
#     "between",
#     ConfigValue(
#         default=None,
#         domain=ListOf(float),
#         description="List of time points to integrate between. If "
#             "None use all time points in the model."
#     )
# )
CONFIG.declare(
    "interpolate_results",
    ConfigValue(
        default=True,
        domain=bool,
        description="If True and trajectory is read, interpolate model "
            "values from the trajectory"
    )
)
CONFIG.declare(
    "calculate_derivatives",
    ConfigValue(
        # TODO do we want to change this to True?
        default=False,
        domain=bool,
        description="If True, calculate the derivative values "
            "based on the values of the differential variables in the discretized "
            "Pyomo model."
    )
)
# CONFIG.declare(
#     "previous_trajectory",
#     ConfigValue(
#         default=None,
#         domain=PetscTrajectory,
#         description="Trajectory from previous integration "
#             "of this model. New results will be appended to this trajectory object."
#     )
# )
CONFIG.declare(
    "representative_time",
    ConfigValue(
        default=None,
        domain=float,
        description="Time when all equations necessary to solve DAE are active. Often equations need "
            "to be deactivated for the initial condition problem (for example, mole fractions summing to one) "
            "because state variables (the individual mole fractions) are fixed. representative_time is a time "
            "after the initial condition problem is solved and these equations are reactivated. Note that the "
            "equations not active at this point are excluded from the DAE at future points. If no "
            "representative_time is specified, it is assumed to be the second element of between. "
            "Must be an element of between."
    )
)

@document_kwargs_from_configdict(CONFIG)
class PETScIntegrator(object):

    def __init__(
            self,
            model: BlockData,
            time_set: ContinuousSet,
            detect_initial: bool=True,
            initial_variables: list=None,
            initial_constraints: list=None,
            **kwargs
    ):

        self._model = model
        self._time_set = time_set

        if "scheme" not in time_set.get_discretization_info():
            raise RuntimeError("The ContinuousSet time has not been discretized")

        scheme = time_set.get_discretization_info()["scheme"]
        if scheme != "BACKWARD Difference" and scheme != "LAGRANGE-RADAU":
            raise RuntimeError(
                "The PETSc interface supports only the backward difference and Lagrange-Radau collocation "
                f"schemes. Found instead the {scheme} scheme."
            )

        self.config = CONFIG(kwargs)

        if self.config.representative_time is None:
            self.config.representative_time = time_set.at(2) # 2 because Pyomo sets start at 1

        if self.config.interpolate_results:
            if self.config.ts_options is None:
                self.config.ts_options = {}
            if "--ts_save_trajectory" in self.config.ts_options:
                # I think ts requires "1" to be passed instead of True to save the trajectory.
                if self.config.ts_options["--ts_save_trajectory"] is True:
                    self.config.ts_options["--ts_save_trajectory"] = 1

                if self.config.ts_options["--ts_save_trajectory"] != 1:
                    raise ValueError(
                        "In order to interpolate model values from the PETSc trajectory, "
                        "the trajectory must be saved. Either set interpolate = False in the "
                        "config or set '--ts_save_trajectory' = 1 in the ts_options dictionary."
                    ) 

                self.config.ts_options["--ts_save_trajectory"] = 1

        # Debugging option to keep the input/output files
        # from calling PETSc. Create it as a private attribute
        # because users should not need to use it.
        self._keepfiles = False
        # This used to be a kwarg for petsc_dae_by_time_element.
        # Leave it as a semiprivate variable in case there is a
        # use case that I am not aware of.
        self._symbolic_solver_labels = True


        self.refresh_flattened_model()

        if initial_variables is None:
            initial_variables = []

        if initial_constraints is None:
            initial_constraints = []

        if detect_initial:
            self.detect_initial_variables_and_constraints(
                initial_variables=initial_variables,
                initial_constraints=initial_constraints
            )
        else:
            self._initial_variables = initial_variables
            self._initial_constraints = initial_constraints

        tdisc = find_discretization_equations(self._model, self._time_set)

        # Workaround for Pyomo bug in which the context manager activates
        # all ConstraintData children of IndexedConstraint object upon
        # exit of the context manager. See Pyomo issue #3734
        disc_condata = []
        for con in tdisc:
            for condata in con.values():
                disc_condata.append(condata)
        self._disc_condata = disc_condata

        self._derivative_differential_vardata_map = _get_derivative_differential_data_map(model, time_set)

    @property
    def atemporal_variables(self):
        return self._atemporal_variables.copy()

    @property
    def time_variables(self):
        return self._time_variables.copy()
    
    @property
    def atemporal_constraints(self):
        return self._atemporal_constraints.copy()
    
    @property
    def time_constraints(self):
        return self._time_constraints.copy()

    @property
    def initial_variables(self):
        return self._initial_variables.copy()
    
    @initial_variables.setter
    def initial_variables(self, initial_variables):
        self._initial_variables = initial_variables

    @property
    def initial_constraints(self):
        return self._initial_constraints.copy()
    
    @initial_constraints.setter
    def initial_constraints(self, initial_constraints):
        self._initial_constraints = initial_constraints

    def refresh_flattened_model(self):
        """
        Flatten the DAE model. Reflattening the model is 
        necessary if components have been activated or deactivated.
        """
        atemporal_variables, time_variables = flatten_dae_components(
            self._model, self._time_set, Var, active=True, indices=(self.config.representative_time,)
        )
        self._atemporal_variables = atemporal_variables
        self._time_variables = time_variables

        atemporal_constraints, time_constraints = flatten_dae_components(
            self._model, self._time_set, Constraint, active=True, indices=(self.config.representative_time,)
        )
        self._atemporal_constraints = atemporal_constraints
        self._time_constraints = time_constraints

    def detect_initial_variables_and_constraints(
            self,
            initial_variables: list=None,
            initial_constraints: list=None
        ):
        """
        Finds variables and constraints that are not indexed by time and
        adds them to user-specified variables and constraints that are included in the 
        initial condition problem but fixed and deactivated, respectively, in the
        time stepping problem.
        """
        rvset = ComponentSet(self._atemporal_variables)
        ivset = ComponentSet(initial_variables)
        self.initial_variables = list(ivset | rvset)

        acset = ComponentSet(self._atemporal_constraints)
        icset = ComponentSet(initial_constraints)
        self.initial_constraints = list(icset | acset)        

    def get_initial_condition_problem(
        self,
        time_point: float = None,
    ):
        if time_point is None:
            time_point = self._time_set.at(1)  # One because Pyomo sets use one-based indexing

        constraints = self.initial_constraints # Returns a copy
        for con in self.time_constraints:
            if time_point in con and con[time_point] not in self._disc_condata:
                constraints.append(con[time_point])

        variables = [var[time_point] for var in self.time_variables] + self.initial_variables
        if len(constraints) > 0:
            t_block = create_subsystem_block(
                constraints,
                variables,
            )
            self._sub_problem_scaling_suffix(self._model, t_block)
        else:
            t_block = None

        return t_block
    
    def get_timestep_problem(
        self,
        time_point: float = None,    
    ):

        variables = [var[time_point] for var in self.time_variables]
        constraints = []
        for con in self.time_constraints:
            if time_point in con and con[time_point] not in self._disc_condata:
                constraints.append(con[time_point])

        variable_fixedness = to_json(
            self._model,
            wts=StoreSpec.isfixed(),
            return_dict=True
        )

        t_block = create_subsystem_block(constraints, variables)

        # Set up the scaling factor suffix
        _sub_problem_scaling_suffix(self._model, t_block)
        
        differential_vars = _set_dae_suffixes_from_variables(
            t_block,
            variables,
            self._derivative_differential_vardata_map,
        )
        # We need to check if there are derivatives in the problem before
        # sending this to the solver.  We'll assume that if you are using
        # this and don't have any differential equations, you are making a
        # mistake.
        if len(differential_vars) < 1:
            raise RuntimeError(
                f"No differential equations found at t = {time_point}, not a DAE"
            )
        if self.config.time_var is not None:
            t_block.dae_suffix[self.config.time_var[time_point]] = int(DaeVarTypes.TIME)
        
        return t_block, differential_vars, variable_fixedness

    def _copy_time(self, t_from, t_to):
        # TODO this can probably be made public
        # but we should take stock of the number of near-identical
        # functions propagated through the codebase
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
        for v in self.time_variables:
            if not v[t_to].fixed:
                v[t_to].value = v[t_from].value

    def dae_by_time_element(self, between=None, skip_initial=False, previous_trajectory=None):
        if between is None:
            between = self._time_set
        else:
            bad_times = between - self._time_set
            if bad_times:
                raise ValueError(
                    "Elements of the 'between' argument must be in the time set"
                )
            between = Set(initialize=sorted(between))
            between.construct()

        # First calculate the initial conditions and non-time-indexed constraints
        res_list = []
        t0 = between.first()

        solve_log = idaeslog.getSolveLogger("petsc-dae")

        solver_dae = SolverFactory("petsc_ts", options=self.config.ts_options)
        save_trajectory = solver_dae.options.get("--ts_save_trajectory", 0)

        if not skip_initial:
            ic_block = self.get_initial_condition_problem(time_point=t0)
            if ic_block is not None:
                # TODO Why aren't we using get_solver here?
                initial_solver_obj = SolverFactory(
                    self.config.initial_solver,
                    options=self.config.initial_solver_options,
                    writer_config=self.config.initial_solver_writer_config,
                )

                with idaeslog.solver_log(solve_log, idaeslog.INFO) as slc:
                    res = initial_solver_obj.solve(
                        ic_block,
                        tee=slc.tee,
                        symbolic_solver_labels=self._symbolic_solver_labels,
                    )
                res_list.append(res)

        tprev = t0
        trajectory = previous_trajectory
        if trajectory is not None:
            variables_prev = [var[t0] for var in self.time_variables]
        else:
            variables_prev = None

        for t in between:
            if t == between.first():
                # t == between.first() was handled above
                continue
            timestep_block, differential_vars, variable_fixedness = self.get_timestep_problem(t)
            self._copy_time(tprev, t)

            with idaeslog.solver_log(solve_log, idaeslog.INFO) as slc:
                res = solver_dae.solve(
                    timestep_block,
                    tee=slc.tee,
                    keepfiles=self._keepfiles,
                    symbolic_solver_labels=self._symbolic_solver_labels,
                    export_nonlinear_variables=differential_vars,
                    options={"--ts_init_time": tprev, "--ts_max_time": t},
                )
                res_list.append(res)

            if save_trajectory:
                trajectory, variables_prev = self._save_trajectory(
                    timestep_block=timestep_block,
                    time_point=t,
                    previous_trajectory=trajectory,
                    previous_variables=variables_prev
                )

            # Revert variable fixedness after solving the timestep
            from_json(self._model, sd=variable_fixedness, wts=StoreSpec.isfixed())
            tprev = t
            res_list.append(res)
        if self.config.interpolate_results:
            self._interpolate_results(between=between, trajectory=trajectory)

        if self.config.calculate_derivatives:
            calculate_time_derivatives(self._model, self._time_set, between=between)

        return PetscDAEResults(results=res_list, trajectory=trajectory)

    def _save_trajectory(self, timestep_block, time_point, previous_trajectory, previous_variables):
        tj = PetscTrajectory(
            stub="tmp_vars_stub",
            delete_on_read=True,
            unscale=True,
            model=timestep_block,
        )
        # add fixed vars to the trajectory. this does two things 1)
        # helps users looking for fixed var trajectory and 2) lets
        # us concatenate trajectories with section that are mixed fixed
        # and unfixed
        variables = [var[time_point] for var in self.time_variables]
        for i, v in enumerate(variables):
            if isinstance(v.parent_component(), DerivativeVar):
                continue  # skip derivative vars
            try:
                vec = tj.get_vec(v)
            except KeyError:
                tj._set_vec(v, [value(v)] * len(tj.time))
        if previous_trajectory is not None:
            # due to the way variables is generated we know variables
            # have corresponding positions in the list
            no_repeat = set()
            for i, v in enumerate(variables):
                vp = previous_variables[i]
                if id(v) in no_repeat:
                    continue  # variables can be repeated in list
                if isinstance(v.parent_component(), DerivativeVar):
                    continue  # skip derivative vars
                no_repeat.add(id(v))
                # We'll add fixed vars in case they aren't fixed in
                # another section. Fixed vars don't go to the solver
                # so they don't show up in the trajectory data
                vec = tj.get_vec(v)
                vec_prev = previous_trajectory.get_vec(vp)
                tj._set_vec(v, vec_prev + vec)
            tj._set_time_vec(previous_trajectory.time + tj.time)
        return tj, variables
    
    def _interpolate_results(self, between, trajectory):
        t0 = between.first()
        tlast = between.last()
        itime = [t for t in self._time_set if not (t < t0 or t > tlast)]
        no_repeat = set()
        if self.config.time_var is not None:
            for t in itime:
                self.config.time_var[t] = t
            no_repeat.add(id(self.config.time_var[tlast]))
        for var in self.time_variables:
            if id(var[tlast]) in no_repeat:
                continue
            no_repeat.add(id(var[tlast]))
            if isinstance(var[t0].parent_component(), DerivativeVar):
                continue  # skip derivative vars
            vec = trajectory.interpolate_vec(itime, var[tlast])
            for i, t in enumerate(itime):
                if t in between:
                    # Time is already set
                    continue
                if not var[t].fixed:
                    # May not have trajectory from fixed variables and they
                    # shouldn't change anyway, so only set not fixed vars
                    var[t].value = vec[i]
