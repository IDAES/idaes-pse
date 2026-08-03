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

import numpy as np
from scipy.linalg import expm
import scipy.sparse.linalg as spla

from pyomo.environ import Constraint, value, Var
from pyomo.common.collections import ComponentSet, ComponentMap, OrderedSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.dae.flatten import flatten_dae_components
from pyomo.util.subsystems import (
    TemporarySubsystemManager,
    create_subsystem_block,
)

from idaes.core.solvers.petsc import _sub_problem_scaling_suffix
from idaes.core.solvers.petsc_object import PETScIntegrator
from idaes.core.util import from_json, StoreSpec, to_json
from idaes.core.util.exceptions import BurntToast
import idaes.core.util.scaling as iscale
from idaes.core.scaling.util import (
    propagate_scaling_factors_to_temporary_block,
    get_jacobian,
)
import idaes.apps.caprese.nmpc_var as nmpc_var
from idaes.apps.caprese.common.config import (
    VariableCategory,
    ConstraintCategory,
)
from idaes.apps.caprese.categorize import (
    _get_state_vardata,
    _get_disc_eq,
    _is_derivative_wrt,
)

DAE_DISC_SUFFIX = "_disc_eq"

VC = VariableCategory
CC = ConstraintCategory

CATEGORY_TYPE_MAP = {
    VariableCategory.DIFFERENTIAL: nmpc_var.DiffVar,
    VariableCategory.ALGEBRAIC: nmpc_var.AlgVar,
    VariableCategory.DERIVATIVE: nmpc_var.DerivVar,
    VariableCategory.INPUT: nmpc_var.InputVar,
    VariableCategory.FIXED: nmpc_var.FixedVar,
    VariableCategory.MEASUREMENT: nmpc_var.MeasuredVar,
}


def _unfix_vars(var_list):
    flags = ComponentMap()
    for var in var_list:
        flags[var] = var.fixed
        var.unfix()
    return flags


def _restore_fixedness(flags):
    for var, fix in flags.items():
        if fix:
            var.fix()


def _identify_derivatives_if_differential(condata, wrt, include_fixed=False):
    # Originally from caprese.categorize, but we want to allow more than one
    # DerivativeVar per equation
    parent = condata.parent_component()
    if parent.local_name.endswith(DAE_DISC_SUFFIX):
        return False, None
    derivs = None
    for var in identify_variables(condata.expr, include_fixed=include_fixed):
        if _is_derivative_wrt(var, wrt):
            if derivs is None:
                derivs = [var]
            else:
                derivs.append(var)
    is_diff = False if derivs is None else True
    return is_diff, derivs


def _deduplicate_components(comp_list, idx, ctype=None):
    """
    Iterates through a list of Pyomo components and removes any duplicate
    components (e.g., References, object references, etc.). In order to
    determine whether the elements of comp_list refer to the same component,
    they need to be evaluated at some index to get a child component data
    object. If they refer to the same data object, they are considered to
    be the same component

    Args:
        comp_list (list): List of components that need to be deduplicated
        idx: Index at which to evaluate each component to get the
            associated data objects
        ctype (Class): Expected type of Pyomo component. If ctype is not
            set to None, then raise an exception if a component is not
            an instance of ctype.
    Returns:
        filtered_list (list): List of components with duplicates removed
        duplicate_list (list): List of duplicate components in comp_list
    """
    visited = set()
    filtered_list = []
    duplicate_list = []
    for comp in comp_list:
        if ctype is not None and not isinstance(comp, ctype):
            raise ValueError(
                f"Encounted component {comp.name} with unexpected type "
                f"{type(comp)}, which is not a subclass of the "
                f"expected type {ctype}."
            )
        _id = id(comp[idx])
        if _id not in visited:
            visited.add(_id)
            filtered_list.append(comp)
        else:
            duplicate_list.append(comp)


def _validate_user_specified_variables(
    vardata_set, input_vardata_set, disturbance_vardata_set, output_vardata_set
):
    """
    Iterates through a set of variables to ensure the user-specified
    input, disturbance, and output variables are present in that set.
    If some variables are not present, then all the absent variables
    are collected to raise a ValueError.

    Args:
        var_set: Variable set which we are checking for inclusion of
            the other variable set.
        input_vardata_set: Set of user-specified input variables
        disturbance_vardata_set: Set of user-specified disturbance variables
        output_vardata_set: Set of user-specified output variables
    Returns:
        None
    """
    missing_inputs = []
    for input in input_vardata_set:
        if not input in vardata_set:
            missing_inputs.append(input.name)

    missing_disturbances = []
    for dist in disturbance_vardata_set:
        if not dist in vardata_set:
            missing_disturbances.append(dist.name)

    missing_outputs = []
    for out in output_vardata_set:
        if not out in vardata_set:
            missing_outputs.append(out.name)

    if missing_inputs or missing_disturbances or missing_outputs:
        err_message = f"User specified variables not present in any active constraint at time {t1}: \n\n"
        if missing_inputs:
            err_message += (
                "Missing input variables:\n" + "\n".join(missing_inputs) + "\n\n"
            )
        if missing_disturbances:
            err_message += (
                "Missing disturbance variables:\n"
                + "\n".join(missing_disturbances)
                + "\n\n"
            )
        if missing_outputs:
            err_message += (
                "Missing output variables:\n" + "\n".join(missing_outputs) + "\n\n"
            )

        raise ValueError(err_message)


# This function was created with the assistence of Google Gemini 3.5 Flash
def _validate_vardata_collections(vardata_lists, vardata_sets):
    """
    In order for linearization to work correctly, we want to
    ensure that the variable categories satisfy some important
    properties.

    1. The lists of variables must have no duplicate elements---
        i.e., the lengths of the sets must equal the lengths
        of the lists
    2. The sets of differential variables, derivative variables,
        input variables, and disturbance variables must be disjoint
    3. All input, disturbance, and output variables must appear
        in the overall set of total variables.

    Args:
        vardata_lists (dict): Dictionary mapping the categories
            of variable to a list of those variables. It must
            have the keys: "total", "deriv", "diff", "alg",
            "input", "output", and "dist"
        vardata_sets (dict): Dictionary mapping the categories
            of variables to a ComponentSet of those variables.
            It must contain the same keys as vardata_lists.
    Returns:
        output_data_lists (dict): Dictionary partitioning the
            output variables into lists of derivative, differential,
            algebraic, input, and disturbance variables.
    """
    # 1. Ensure no duplicate elements exist within individual lists
    for key, lst in vardata_lists.items():
        if len(lst) != len(vardata_sets[key]):
            # Find specific duplicated variables
            seen = ComponentSet()
            duplicates = []
            for var in lst:
                if var in seen:
                    duplicates.append(var)
                else:
                    seen.add(var)
            raise ValueError(
                f"Duplicate variables found in the '{key}' variable list: {duplicates}"
            )

    # 2. Ensure differential, derivative, input, and disturbance sets are disjoint
    disjoint_categories = ["diff", "deriv", "input", "dist"]
    for i, cat_a in enumerate(disjoint_categories):
        for cat_b in disjoint_categories[i + 1 :]:
            if not vardata_sets[cat_a].isdisjoint(vardata_sets[cat_b]):
                intersection = [
                    v for v in vardata_lists[cat_a] if v in vardata_sets[cat_b]
                ]
                raise ValueError(
                    f"The variable sets '{cat_a}' and '{cat_b}' are not disjoint. "
                    f"Shared variables: {intersection}"
                )

    # 3. Ensure input, disturbance, and output variables are in the total set
    total_set = vardata_sets["total"]
    for cat in ["input", "dist", "output"]:
        for var in vardata_lists[cat]:
            if var not in total_set:
                raise ValueError(
                    f"Variable {var} in '{cat}' is not present in the total variables set."
                )

    # Partition output variables into lists preserving their original order
    output_data_lists = {
        "deriv": ComponentSet(),
        "diff": ComponentSet(),
        "alg": ComponentSet(),
        "input": ComponentSet(),
        "dist": ComponentSet(),
    }

    for var in vardata_lists["output"]:
        if var in vardata_sets["deriv"]:
            output_data_lists["deriv"].add(var)
        elif var in vardata_sets["diff"]:
            output_data_lists["diff"].add(var)
        elif var in vardata_sets["alg"]:
            output_data_lists["alg"].add(var)
        elif var in vardata_sets["input"]:
            output_data_lists["input"].add(var)
        elif var in vardata_sets["dist"]:
            output_data_lists["dist"].add(var)
        else:
            # Since "alg" is mathematically defined as "total - (diff + deriv + input + dist)",
            # any output variable belonging to "total" must fall into one of these 5 sets.
            raise BurntToast(
                f"Output variable {var} is in the total variables set but could not be "
                f"partitioned into derivative, differential, algebraic, input, or disturbance subsets."
                "Encountering this error indicates either a developer oversight or that the "
                "implementation of ComponentSet has been changed. Please open an issue on the "
                "IDAES Github with steps to reproduce this error."
            )

    return output_data_lists


def _validate_steady_state(
    t_block, deriv_var_set, steady_state_derivative_tolerance, constraint_tolerance
):
    """
    Validate whether the dynamic system is at steady state at
    the point where linearization occurs. If a constraint has
    a large (scaled) residual or a derivative variable has a
    (scaled) value sufficiently far away from zero.

    Raises a ValueError if it detects the system is not at
    steady state.

    Args:
        t_block: Temporary block containing timestep problem
        deriv_var_set: Set of time DerivativeVars
        steady_state_derivative_tolerance: Tolerance for how far
            the (scaled) time derivatives can be away from zero
            before an exception is raised.
        steady_state_constraint_tolerance: Tolerance for how large
            the (scaled) constraint residual can be before an
            exception is raised.

    """

    def _get_scaling_factor(compdata):
        if compdata in t_block.scaling_factor:
            return t_block.scaling_factor[compdata]
        else:
            return 1

    nonzero_deriv_variables = []
    for vardata in deriv_var_set:
        sf = _get_scaling_factor(vardata)
        if sf * abs(value(vardata)) > steady_state_derivative_tolerance:
            nonzero_deriv_variables.append(
                vardata.name + f": {value(vardata):.2e} (scaling factor={sf:.2e})\n"
            )

    nonzero_constraint_residuals = []
    for condata in t_block.component_data_objects(ctype=Constraint):
        sf = _get_scaling_factor(condata)
        if sf * abs(value(condata)) > constraint_tolerance:
            nonzero_constraint_residuals.append(
                condata.name + f": {value(condata):.2e} (scaling factor={sf:.2e})\n"
            )

    err_msg = ""
    if nonzero_deriv_variables:
        err_msg += (
            "Nonzero derivative variables: \n\n"
            + "\n".join(nonzero_deriv_variables)
            + "\n\n"
        )
    if nonzero_constraint_residuals:
        err_msg += (
            "Nonzero constraint residuals: \n\n"
            + "\n".join(nonzero_constraint_residuals)
            + "\n\n"
        )

    if err_msg:
        raise ValueError("The system was not at steady state: \n\n" + err_msg)


def linearize_system(
    model,
    time,
    representative_time=None,
    input_vars=None,
    disturbance_vars=None,
    output_vars=None,
    scaled=False,
    steady_state_derivative_tolerance=1e-6,
    constraint_tolerance=1e-6,
):
    """
    Function to obtain a linear time invariant state space model by
    linearizing a set of differential algebraic equations. It takes
    the general system
    0 = f(x, xdot, u, d)
    y = h(x, u, d)
    and turns it into the linear system
    xdot = A @ x + B @ u + B_d @ d
    y = C @ X + D @ u + D_d @ d

    Note that the sets of input, disturbance, and output variables
    should be indexed only by time. This indexing can be achieved
    by creating a Reference to a slice of a variable that is indexed
    by multiple sets or is only indirectly indexed by time, e.g., the
    variable itself is a scalar but its parent block is indexed by time.

    Args:
        model (Block): Pyomo model that will be linearized
        time (ContinuousSet): Time domain for set of DAEs
        representative_time (float): Time index at which the correct
            equations are active and variables are fixed to create
            a well-defined system of DAEs. If the differential
            variables, input variables, and disturbance variables are
            fixed, the system should be square and be able to be solved
            for the algebraic and derivative variables.
        scaled (bool): Whether to return a linear system with variables
            adjusted by their scaling factors, or to unscale the
            system matrices before returning.
        input_variables (list): List of input (u) variables.
        disturbance_variables (list): List of disturbance (d) variables.
        output_variables (list): List of output (y) variables.
        steady_state_derivative_tolerance (float): Tolerance used to
            test whether the time derivative variables are close
            enough to zero to consider that the system is at steady state.
        constraint_tolerance (float): Tolerance used to test whether the
        DAE constraints are close enough to being satisfied to consider
            that the system 0 = f(x, xdot, u, d) is satisfied.
    Returns:
        Dictionary containing the following fields:
            "diff_vars": The list of differential VarData in the order
                they appear in the A and C matrices
            "alg_vars": The list of algebraic VarData
            "input_vars": The list of input VarData in the order they
                appear in the B and D matrices
            "disturbance_vars": The list of disturbance VarData in the
                order they appear in the B_d and D_d matrices
            "output_vars": The list of output VarData in the order
                they appear in the C, D, and D_d matrices
            "A" (numpy.array): The A matrix
            "B" (numpy.array): The B matrix
            "B_d" (numpy.array): The B_d matrix
            "C" (numpy.array): The C matrix
            "D" (numpy.array): The D matrix
            "D_d" (numpy.array): The D_d matrix
    """
    if representative_time is not None:
        t1 = representative_time
    else:
        # Use the first non-initial time point as a "representative
        # index." Don't use get_finite_elements so this will be valid
        # for general ordered sets.
        t1 = time.at(2)

    if input_vars is None:
        input_vars = []
    if disturbance_vars is None:
        disturbance_vars = []
    if output_vars is None:
        output_vars = []

    vardata_lists = {}
    vardata_sets = {}

    vardata_lists["input"] = [inp[t1] for inp in input_vars]
    vardata_lists["dist"] = [dist[t1] for dist in disturbance_vars]
    vardata_lists["output"] = [out[t1] for out in output_vars]

    fixedness = to_json(model, return_dict=True, wts=StoreSpec.isfixed())
    vartypes = ["input", "dist", "output"]

    try:
        # Need inputs, disturbances, and outputs to be unfixed to ensure that
        # they show up in the jacobian
        for vartype in vartypes:
            for vardata in vardata_lists[vartype]:
                vardata.unfix()
        integrator = PETScIntegrator(model, time, representative_time=t1)
        # t_block is a Pyomo block containing the equations necessary to solve
        # for the derivative and algebraic variables so long as the input,
        # disturbance, and differential variables are fixed.
        # In other words, it contains the equations to formulate
        # the DAE system 0 = f(x, xdot, u, d)
        t_block, _, _ = integrator.get_timestep_problem(t1)
        jac, nlp = get_jacobian(t_block, equality_constraints_only=True)
    finally:
        from_json(model, fixedness, wts=StoreSpec.isfixed())

    vardata_lists["total"] = nlp.get_pyomo_variables()

    vardata_lists["deriv"] = []
    vardata_lists["diff"] = []
    for deriv, out in integrator.derivative_differential_vardata_map.items():
        vardata_lists["deriv"].append(deriv)
        vardata_lists["diff"].append(out["diff_var"])

    vartypes += ["deriv", "diff", "total"]
    vardata_sets = {}
    # Cast lists to a ComponentSet to remove any duplicate elements
    for vartype in vartypes:
        vardata_sets[vartype] = ComponentSet(vardata_lists[vartype])

    vartypes.append("alg")
    vardata_sets["alg"] = vardata_sets["total"] - (
        vardata_sets["diff"]
        | vardata_sets["deriv"]
        | vardata_sets["input"]
        | vardata_sets["dist"]
    )
    # Iterate through vardata_lists["total"] instead
    # of directly converting vardata_sets["alg"] to
    # ensure a consistent variabel order is maintained
    vardata_lists["alg"] = []
    for vardata in vardata_lists["total"]:
        if vardata in vardata_sets["alg"]:
            vardata_lists["alg"].append(vardata)

    output_lists = _validate_vardata_collections(vardata_lists, vardata_sets)

    dof = (
        nlp.n_primals()
        - nlp.n_constraints()
        - len(vardata_lists["input"])
        - len(vardata_lists["dist"])
        - len(vardata_lists["diff"])
    )
    if dof != 0:
        raise RuntimeError(
            f"Found {dof} degrees of freedom at {t1}. "
            f"The timestep problem at {t1} with the "
            "input, disturbance, and differential fixed "
            "and the output variables unfixed should be "
            "square."
        )

    _validate_steady_state(
        t_block,
        vardata_sets["deriv"],
        steady_state_derivative_tolerance,
        constraint_tolerance,
    )

    raw_jac_indices = {}
    for key, idx_list in vardata_lists.items():
        raw_jac_indices[key] = nlp.get_primal_indices(idx_list)

    out = {
        "scaled_jac": jac,
        "nlp": nlp,
        "diff_vars": vardata_lists["diff"],
        "alg_vars": vardata_lists["alg"],
    }

    jac_deriv_alg = jac[:, raw_jac_indices["deriv"] + raw_jac_indices["alg"]]
    jac_rest = jac[
        :, raw_jac_indices["diff"] + raw_jac_indices["input"] + raw_jac_indices["dist"]
    ]

    sys_raw = spla.spsolve(-jac_deriv_alg, jac_rest.tocsc()).todense()

    # We seek an LTI system of the form:
    # xdot = A @ x + B @ u + B_d @ d
    # y = C @ x + D @ u + D_d @ d
    # We don't actually need estimates of the algebraic states in the LTI system,
    # however some outputs (y's) may be algebraic states.
    # Most of these matrices come from indexing and rearranging sys_raw, with a
    # few elements of C being appropriate unit vectors

    nx = len(vardata_lists["diff"])
    nz = len(vardata_lists["alg"])
    nu = len(vardata_lists["input"])
    nd = len(vardata_lists["dist"])
    ny = len(vardata_lists["output"])

    AB_raw = sys_raw[:nx, :]

    out["A"] = AB_raw[:, :nx]
    out["B"] = AB_raw[:, nx : nx + nu]
    out["Bd"] = AB_raw[:, nx + nu : nx + nu + nd]

    out["C"] = np.zeros((ny, nx))
    out["D"] = np.zeros((ny, nu))
    out["Dd"] = np.zeros((ny, nd))

    for Crow, output, out_idx in zip(
        range(ny), vardata_lists["output"], raw_jac_indices["output"]
    ):
        if output in output_lists["diff"]:
            # Because we've rearranged the rows/columns when creating
            # A, we need to find the *new* index corresponding to the
            # differential variable
            di = raw_jac_indices["diff"].index(out_idx)
            out["C"][Crow, di] = 1
        elif output in output_lists["alg"]:
            ai = len(raw_jac_indices["deriv"]) + raw_jac_indices["alg"].index(out_idx)
            out["C"][Crow, :] = sys_raw[ai, raw_jac_indices["diff"]]
            out["D"][Crow, :] = sys_raw[ai, raw_jac_indices["input"]]
            out["Dd"][Crow, :] = sys_raw[ai, raw_jac_indices["dist"]]
        elif output in output_lists["deriv"]:
            # Grab the correct rows of A, B, and Bd
            drvi = raw_jac_indices["deriv"].index(out_idx)
            out["C"][Crow, :] = out["A"][drvi, :]
            out["D"][Crow, :] = out["B"][drvi, :]
            out["Dd"][Crow, :] = out["Bd"][drvi, :]
        elif output in output_lists["input"]:
            # The C row is full of zeros
            ini = len(raw_jac_indices["diff"]) + vardata_lists["input"].index(out_idx)
            out["D"][Crow, ini] = 1
            # The Dd row is full of zeros
        elif output in output_lists["dist"]:
            # The C and D rows are full of zeros
            ind = vardata_lists["dist"].index(out_idx)
            out["Dd"][Crow, ind] = 1
        else:
            raise BurntToast(
                "This branch should be inaccessible. Please open an issue on the IDAES "
                "Github with steps to reproduce this exception so that the IDAES developers "
                "can address this problem."
            )
    if not scaled:
        all_sfs = nlp.get_primals_scaling()
        sfx = all_sfs[raw_jac_indices["diff"]]
        sfxdot = all_sfs[raw_jac_indices["deriv"]]
        sfy = all_sfs[raw_jac_indices["output"]]
        sfu = all_sfs[raw_jac_indices["input"]]
        sfd = all_sfs[raw_jac_indices["dist"]]

        # This step could be accomplished more efficiently using
        # array broadcasting, but using matrices avoids bugs
        # involving matrices broadcast along the wrong axis
        out["A"] = np.diag(1 / sfxdot) @ out["A"] @ np.diag(sfx)
        out["B"] = np.diag(1 / sfxdot) @ out["B"] @ np.diag(sfu)
        out["Bd"] = np.diag(1 / sfxdot) @ out["Bd"] @ np.diag(sfd)

        out["C"] = np.diag(1 / sfy) @ out["C"] @ np.diag(sfx)
        out["D"] = np.diag(1 / sfy) @ out["D"] @ np.diag(sfu)
        out["Dd"] = np.diag(1 / sfy) @ out["Dd"] @ np.diag(sfd)

    return out


def c2d(sys: dict, dt: float):
    r"""
    Convert a linear time invariant system from continuous time to discrete time.
    Start with a system of the form

    .. math:: \dot{x} = A x + B u + B_d d

    and create a system

    ..math:: x(t + \Delta t) = \tilde{A} x +  \tilde{B} u + \tilde{B}_d d

    Note that the output matrices :math:`C`, :math:`D`, and :math:`D_d`
    do not need to be converted.

    Args:
        sys (dict): A dictionary containing the fields:
            "A" (numpy array)
            "B" (numpy array)
            "Bd" (numpy array)
    Returns:
        out (dict): A dictionary containing the following fields:
            "A" (numpy array)
            "B" (numpy array)
            "Bd" (numpy array)
    """

    A = sys["A"]
    if "B" in sys:
        B = sys["B"]
    else:
        B = np.zeros((A.shape[0], 0))

    if "Bd" in sys:
        Bd = sys["Bd"]
    else:
        Bd = np.zeros((A.shape[0], 0))

    nx = A.shape[0]
    nu = B.shape[1]
    nd = Bd.shape[1]

    A_aug = np.block([[A, B, Bd], [np.zeros((nu + nd, nx + nu + nd))]])

    A_tilde_aug = expm(A_aug * dt)

    out = {}
    out["A"] = A_tilde_aug[:nx, :nx]
    out["B"] = A_tilde_aug[:nx, nx : nx + nu]
    out["Bd"] = A_tilde_aug[:nx, nx + nu : nx + nu + nd]

    return out
