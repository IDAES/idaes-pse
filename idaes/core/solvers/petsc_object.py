#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains the PETScIntegrator object to manage calls to the PETSc
TS integrator. It caches a flattened version of the model in order to speed up
repeated calls to the integrator and includes the get_initial_condition_problem
and get_timestep_problem functions to help troubleshoot issues with dynamic
models.

This module also contains the get_derivative_differential_vardata_map function,
which builds a map between derivative variables and their corresponding
differential variables and discretization equations.

Authors: Douglas Allan, John Eslick
"""

from pyomo.core.base import BlockData
from pyomo.environ import Constraint, Set, SolverFactory, Suffix, value, Var
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.config import (
    ConfigDict,
    ConfigValue,
    document_kwargs_from_configdict,
    IsInstance,
    ListOf,
)
from pyomo.dae.flatten import flatten_dae_components
from pyomo.util.subsystems import create_subsystem_block

import idaes.logger as idaeslog

from idaes.core.scaling.util import propagate_scaling_factors_to_temporary_block
from idaes.core.solvers.petsc import (
    calculate_time_derivatives,
    DaeVarTypes,
    PetscDAEResults,
    PetscTrajectory,
)

from idaes.core.util.model_serializer import StoreSpec, from_json, to_json

DAE_DISC_SUFFIX = "_disc_eq"
DAE_CONT_SUFFIX = "_cont_eq"


def get_derivative_differential_vardata_map(
    blk: BlockData,
    cont_set: ContinuousSet,
    raise_higher_derivative_exception: bool = False,
):
    """
    Creates or updates the _derivative_differential_vardata_map attribute
    to contain a map from data objects of derivative variables to the
    corresponding data objects of differential variables, discretization
    equations and/or continuity equations (for Lagrange Legendre collocation).
    "Differential variables" are DerivativeVar objects who have a corresponding
    active discretization or continuity equation.

    This function should be re-run if the model structurally changes (i.e., a
    variable is fixed or unfixed, or a constraint is activated/deactivated).

    Args:
        blk (BlockData): Block which will be searched for DerivativeVars
        cont_set (ContinuousSet): This function searches for DerivativeVars
            which have been differentiated by this set.
        raise_higher_derivative_exception(Exception): Raise a ValueError if
            this function encounters a derivative with respect to cont_set
            that is of order greater than one. Higher derivatives
            that do not involve cont_set are ignored.

    Returns:
        deriv_diff_map (ComponentMap): Map from derivative variables to a
            dictionary containing several keys: "diff_var", which is mapped
            to the corresponding differential variable, "disc_eq", which
            is mapped to the corresponding discretization equation (if it
            exists), and "cont_eq", which is mapped to a continuity equation
            (if it exists).

    """
    cont_set_name = cont_set.name
    deriv_diff_map = ComponentMap()
    for var in blk.component_objects(Var):
        if isinstance(var, DerivativeVar) and cont_set in var.get_continuousset_list():
            if (
                raise_higher_derivative_exception
                and len(var.get_continuousset_list()) > 1
            ):
                raise ValueError(
                    f"Unexpectedly encountered {var.name}, which is differentiated multiple times."
                )

            deriv = var
            diffvar = deriv.get_state_var()
            block = deriv.parent_block()
            disc_eq = block.find_component(var.local_name + DAE_DISC_SUFFIX)
            cont_eq = block.find_component(
                diffvar.local_name + "_" + cont_set_name + DAE_CONT_SUFFIX
            )

            for idx in var:
                out_map = {}
                if (
                    disc_eq is not None
                    and disc_eq.active
                    and idx in disc_eq
                    and disc_eq[idx].active
                ):
                    out_map["disc_eq"] = disc_eq[idx]
                if (
                    cont_eq is not None
                    and cont_eq.active
                    and idx in cont_eq
                    and cont_eq[idx].active
                ):
                    out_map["cont_eq"] = cont_eq[idx]
                if out_map:
                    out_map["diff_var"] = diffvar[idx]
                    deriv_diff_map[deriv[idx]] = out_map

    # The old function had a step to filter out derivative variables that
    # are not present in active constraints. Checking the activity of the
    # corresponding discretization equation should be sufficient to
    # determine whether a variable is active. We have an edge case with
    # Lagrange-Legendre collocation, because the presence of an active
    # continuity equation does not guarantee the derivative variable appears
    # in any equation. For now, assume that, if a continuity equation exists,
    # the user will provide additional validation when appropriate.
    return deriv_diff_map


CONFIG = ConfigDict()

CONFIG.declare(
    "time_var",
    ConfigValue(
        default=None,
        domain=IsInstance(Var),
        description="Optional specification of a time variable, which can be"
        "used to write constraints that are an explicit function of time.",
    ),
)
CONFIG.declare(
    "initial_constraints",
    ConfigValue(
        default=None,
        domain=ListOf(Constraint),
        description="Constraints to solve in the initial "
        "condition solve step.  Since the time-indexed constraints are picked "
        "up automatically, this generally includes non-time-indexed "
        "constraints.",
    ),
)
CONFIG.declare(
    "initial_variables",
    ConfigValue(
        default=None,
        domain=ListOf(Var),
        description="initial_variables (list): This is a list of variables to fix after the "
        "initial condition solve step.  If these variables were originally "
        "unfixed, they will be unfixed at the end of the solve. This usually "
        "includes non-time-indexed variables that are calculated along with "
        "the initial conditions.",
    ),
)
CONFIG.declare(
    "detect_initial",
    ConfigValue(
        default=True,
        domain=bool,
        description="If True, add non-time-indexed variables and "
        "constraints to initial_variables and initial_constraints.",
    ),
)
CONFIG.declare(
    "initial_solver",
    ConfigValue(
        default="petsc_snes",
        domain=str,
        description="The nonlinear equations solver "
        "to use for the initial conditions (e.g. petsc_snes, ipopt, ...).",
    ),
)
CONFIG.declare(
    "initial_solver_options",
    ConfigValue(
        default=None,
        domain=dict,
        description="Solver options to use with initial_solver.",
    ),
)
CONFIG.declare(
    "initial_solver_writer_config",
    ConfigValue(
        default=None,
        domain=dict,
        description="Configurations to use with the .nl writer "
        "in the initial condition problem.",
    ),
)
CONFIG.declare(
    "ts_options",
    ConfigValue(
        default=None,
        domain=dict,
        description="Options to use with the PETSc integrator TS.",
    ),
)
CONFIG.declare(
    "interpolate_results",
    ConfigValue(
        default=True,
        domain=bool,
        description="If True and trajectory is read, interpolate model "
        "values from the trajectory",
    ),
)
CONFIG.declare(
    "calculate_derivatives",
    ConfigValue(
        # TODO do we want to change the default to True?
        default=False,
        domain=bool,
        description="If True, calculate the derivative values "
        "based on the values of the differential variables in the discretized "
        "Pyomo model.",
    ),
)
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
        "Must be an element of between.",
    ),
)


@document_kwargs_from_configdict(CONFIG)
class PETScIntegrator(object):
    """
    Class to interface with the PETSc TS integrator.

    Args:
        model (BlockData): Pyomo model or block which will be integrated.
        time_set (ContinuousSet): The set over which the model will be integrated.
            This is typically time, but a length domain could also be used.
        detect_initial (bool): Whether to guess which variables and constraints
            should go on the initial condition problem or whether only those set
            by the user should be added. See the detect_initial method for more
            information. Default is True.
        initial_variables (list): This is a list of variables to fix after the
            initial condition solve step.  If these variables were originally
            unfixed, they will be unfixed at the end of the solve. This usually
            includes non-time-indexed variables that are calculated along with
            the initial conditions.
        initial_constraints (list): Constraints to solve in the initial
            condition solve step.  Since the time-indexed constraints are picked
            up automatically, this generally includes non-time-indexed
            constraints.

    """

    def __init__(
        self,
        model: BlockData,
        time_set: ContinuousSet,
        detect_initial: bool = True,
        initial_variables: list = None,
        initial_constraints: list = None,
        **kwargs,
    ):

        if not isinstance(time_set, ContinuousSet):
            raise ValueError("Argument time_set is not a Pyomo ContinuousSet")

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
            self.config.representative_time = time_set.at(
                2
            )  # 2 because Pyomo sets start at 1

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

        self.refresh_model()

        if initial_variables is None:
            initial_variables = []

        if initial_constraints is None:
            initial_constraints = []

        if detect_initial:
            self.detect_initial_variables_and_constraints(
                initial_variables=initial_variables,
                initial_constraints=initial_constraints,
            )
        else:
            self._initial_variables = initial_variables
            self._initial_constraints = initial_constraints

        disc_condata = []
        for data_dict in self._derivative_differential_vardata_map.values():
            if "disc_eq" in data_dict:
                disc_condata.append(data_dict["disc_eq"])
        self._disc_condata = disc_condata

    @property
    def atemporal_variables(self):
        """
        List of variables not indexed by time.
        """
        return self._atemporal_variables.copy()

    @property
    def time_variables(self):
        """
        List of variables indexed by time.
        """
        return self._time_variables.copy()

    @property
    def atemporal_constraints(self):
        """
        List of constraints not indexed by time.
        """
        return self._atemporal_constraints.copy()

    @property
    def time_constraints(self):
        """
        List of constraints indexed by time.
        """
        return self._time_constraints.copy()

    @property
    def initial_variables(self):
        """
        List of additional variables to add to
        the initial condition problem.
        """
        return self._initial_variables.copy()

    @initial_variables.setter
    def initial_variables(self, initial_variables):
        self._initial_variables = list(initial_variables)

    @property
    def initial_constraints(self):
        """
        List of additional constraints to add to
        initial condition problem.
        """
        return self._initial_constraints.copy()

    @initial_constraints.setter
    def initial_constraints(self, initial_constraints):
        self._initial_constraints = list(initial_constraints)

    def _validate_no_fixed_nonzero_derivatives(self):
        for deriv_vardata in self._derivative_differential_vardata_map.keys():
            if deriv_vardata.fixed and value(deriv_vardata) != 0:
                raise RuntimeError(
                    f"{deriv_vardata} is fixed to a nonzero value {value(deriv_vardata)}. "
                    "This is most likely a modeling error. Instead of fixing the "
                    "derivative consider adding a constraint like dxdt = constant."
                )

    def _set_dae_suffixes_from_variables(self, blk, variables):
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
        blk.dae_suffix = Suffix(
            direction=Suffix.EXPORT,
            datatype=Suffix.INT,
        )
        # The dae_link suffix provides the solver a link between the differential
        # and derivative variable, by assigning the pair a unique integer index.
        blk.dae_link = Suffix(
            direction=Suffix.EXPORT,
            datatype=Suffix.INT,
        )
        dae_var_link_index = 1
        differential_vars = []
        i = 0
        for var in variables:
            if var in self._derivative_differential_vardata_map:
                deriv = var
                diffvar = self._derivative_differential_vardata_map[deriv]["diff_var"]
                blk.dae_suffix[diffvar] = int(DaeVarTypes.DIFFERENTIAL)
                blk.dae_suffix[deriv] = int(DaeVarTypes.DERIVATIVE)
                blk.dae_link[diffvar] = dae_var_link_index
                blk.dae_link[deriv] = dae_var_link_index
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

    def refresh_model(self):
        """
        Flatten the DAE model and recreates the derivative-differential variable map.
        Refreshing the model is necessary if variables have been fixed or unfixed or
        components have been activated or deactivated.

        Args:
            None
        Returns:
            None
        """
        atemporal_variables, time_variables = flatten_dae_components(
            self._model,
            self._time_set,
            Var,
            active=True,
            indices=(self.config.representative_time,),
        )
        self._atemporal_variables = atemporal_variables
        self._time_variables = time_variables

        atemporal_constraints, time_constraints = flatten_dae_components(
            self._model,
            self._time_set,
            Constraint,
            active=True,
            indices=(self.config.representative_time,),
        )
        self._atemporal_constraints = atemporal_constraints
        self._time_constraints = time_constraints

        try:
            self._derivative_differential_vardata_map = (
                get_derivative_differential_vardata_map(
                    self._model, self._time_set, True
                )
            )
        except ValueError as err:
            if "differentiated multiple times" in str(err):
                raise NotImplementedError(
                    "IDAES presently does not support PETSc for second order or higher derivatives "
                    "that are differentiated at least once with respect to time. Please reformulate your model so "
                    "it does not contain such a derivative (such as by introducing intermediate variables)."
                ) from err
            else:
                raise err

    def detect_initial_variables_and_constraints(
        self, initial_variables: list = None, initial_constraints: list = None
    ):
        """
        Finds variables and constraints that are not indexed by time and
        adds them to user-specified variables and constraints that are included in the
        initial condition problem but fixed and deactivated, respectively, in the
        time stepping problem. It stores these values on the initial_variables
        and initial_constraints attributes.

        Args:
            initial_variables (list): List of user-specified variables to include in
                initial condition problem.
            initial_constraints (list): List of user-specified constraints to include
                in the initial condition problem.
        Returns:
            None
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
        """
        From the cached, flattened model, construct an initial condition
        problem at a given time point. Solving this problem allows us
        to ensure that atemporal constraints are satisfied (such as those
        involving unit geometry, which does not change over time) as well
        as ensuring all the algebraic equations are satisfied at the
        initial condition. This includes the variables specified on
        the initial_variables and initial_constraint attributes, as well
        as all time-indexed variables at time_point.

        This function is useful for diagnosing failures during solution
        of the initial condition problem.

        Args:
            time_point (float): Time point at which to create initial
                condition problem

        Returns:
            Block: Subsystem block containing References to the variables and
                constraints used in the initial condition problem. This block
                also contains the scaling_factor Suffix containing variable
                and constraint scaling factors from the original model.

        """
        if time_point is None:
            time_point = self._time_set.at(
                1
            )  # One because Pyomo sets use one-based indexing

        constraints = self.initial_constraints  # Returns a copy
        for con in self.time_constraints:
            if time_point in con and con[time_point] not in self._disc_condata:
                constraints.append(con[time_point])

        variables = [
            var[time_point] for var in self.time_variables
        ] + self.initial_variables
        if len(constraints) > 0:
            t_block = create_subsystem_block(
                constraints,
                variables,
            )
            propagate_scaling_factors_to_temporary_block(t_block)
        else:
            t_block = None

        return t_block

    def get_timestep_problem(
        self,
        time_point: float = None,
    ):
        """
        From the cached, flattened model, construct a problem to solve
        for (time) derivative variables and algebraic variables at a
        given time point. Any atemporal variables that are included
        in the constraints are fixed.  If the user fixes the
        differential variables (which are not fixed by default), the
        block should constitute a square system which can be solved
        for the (time) derivative variables and algebraic variables.

        This function is useful for diagnosing model issues when
        the TS integrator fails to solve or says the model is nonsquare.

        Args:
            time_point (float): Time at which to create the timestep
                problem.

        Returns:
            t_block (BlockData): Pyomo block with the variables and
                constraints necessary to solve the initial condition
                problem as References. Variables that occur in these
                constraints that do not occur at this time point are
                included in the input_vars list and fixed.
            differential_vars (List): List of differential variables
                at this time. If these are fixed by the user, t_block
                should contain a square system.
            variable_fixedness (dict): Dictionary serialization of the
                model containing the fixedness status of all variables
                and constraints. This can be used with the from_json
                function from idaes.core.util in order to restore the
                fixedness of the input_vars. The differential_vars
                are not fixed by this function.

        """
        self._validate_no_fixed_nonzero_derivatives()

        # Use a ComponentSet to avoid variables from being counted
        # twice if a Reference to that variable exists
        variables = ComponentSet()
        for var in self.time_variables:
            variables.add(var[time_point])

        constraints = []
        for con in self.time_constraints:
            if time_point in con and con[time_point] not in self._disc_condata:
                constraints.append(con[time_point])

        variable_fixedness = to_json(
            self._model, wts=StoreSpec.isfixed(), return_dict=True
        )
        t_block = create_subsystem_block(constraints, list(variables))
        # Fix input variables (which shouldn't be indexed by time)
        for vardata in t_block.input_vars.values():
            vardata.fix()

        # Set up the scaling factor suffix
        propagate_scaling_factors_to_temporary_block(t_block)

        differential_vars = self._set_dae_suffixes_from_variables(
            t_block,
            variables,
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
        """
        This is used on the cached, flattened (indexed only by time) variable
        representations to copy variable values that are unfixed at the "to" time
        from the value at the "from" time. The PETSc DAE solver uses the initial
        variable values as the initial condition, so this is used to copy the
        previous time in as the initial condition for the next step.

        Args:
            t_from (float): time point to copy from
            t_to (float): time point to copy to, only unfixed vars will be
                overwritten

        Returns:
            None
        """
        for v in self.time_variables:
            if not v[t_to].fixed:
                v[t_to].value = v[t_from].value

    def dae_by_time_element(
        self, between=None, skip_initial=False, previous_trajectory=None
    ):
        """Solve a DAE problem step by step using the PETSc DAE solver.  This
        integrates from one time point to the next.

        Args:
            between (list or tuple): List of time points to integrate between. If
                None use all time points in the model.
            skip_initial (bool): Don't do the initial condition calculation step,
                and assume that the initial condition values have already been
                calculated. This can be useful, for example, if you read initial
                conditions from a separately solved steady state problem, or
                otherwise know the initial conditions.
            previous_trajectory: (PetscTrajectory) Trajectory from previous integration
                of this model. New results will be appended to this trajectory object.

        Returns (PetscDAEResults):
            See PetscDAEResults documentation for more information.
        """
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
            timestep_block, differential_vars, variable_fixedness = (
                self.get_timestep_problem(t)
            )
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
                    previous_variables=variables_prev,
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

    def _save_trajectory(
        self, timestep_block, time_point, previous_trajectory, previous_variables
    ):
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
                # pylint: disable-next = protected-access
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
                tj._set_vec(v, vec_prev + vec)  # pylint: disable = protected-access
            # pylint: disable-next = protected-access
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
