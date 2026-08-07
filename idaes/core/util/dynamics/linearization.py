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
from numpy import any, block, diag, isnan, zeros
from scipy.linalg import expm
from scipy.sparse.linalg import spsolve

from pyomo.environ import Constraint, value
from pyomo.common.collections import ComponentSet, ComponentMap

from idaes.core.solvers.petsc_object import PETScIntegrator
from idaes.core.util import from_json, StoreSpec, to_json
from idaes.core.util.exceptions import BurntToast
from idaes.core.scaling.util import get_jacobian


# This function was created with the assistance of Google Gemini 3.5 Flash
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
        output_sets (dict): Dictionary partitioning the output
            variables into ComponentSets of derivative, differential,
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
                    duplicates.append(var.name)
                else:
                    seen.add(var)
            raise ValueError(
                f"Duplicate variables found in the '{key}' "
                "variable list: {duplicates}"
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
                    f"The variable sets '{cat_a}' and '{cat_b}' are not "
                    "disjoint. Shared variables: {intersection}"
                )

    # 3. Ensure input, disturbance, and output variables are in the total set
    total_set = vardata_sets["total"]
    for cat in ["input", "dist", "output"]:
        for var in vardata_lists[cat]:
            if var not in total_set:
                raise ValueError(
                    f"Variable {var} in '{cat}' is not present in the total "
                    "variables set."
                )

    # Partition output variables into lists preserving their original order
    output_sets = {
        "deriv": ComponentSet(),
        "diff": ComponentSet(),
        "alg": ComponentSet(),
        "input": ComponentSet(),
        "dist": ComponentSet(),
    }

    for var in vardata_lists["output"]:
        if var in vardata_sets["deriv"]:
            output_sets["deriv"].add(var)
        elif var in vardata_sets["diff"]:
            output_sets["diff"].add(var)
        elif var in vardata_sets["alg"]:
            output_sets["alg"].add(var)
        elif var in vardata_sets["input"]:
            output_sets["input"].add(var)
        elif var in vardata_sets["dist"]:
            output_sets["dist"].add(var)
        else:
            # Since "alg" is mathematically defined as "total - (diff + deriv + input + dist)",
            # any output variable belonging to "total" must fall into one of these 5 sets.
            raise BurntToast(
                f"Output variable {var} is in the total variables set but "
                f"could not be partitioned into derivative, differential, "
                "algebraic, input, or disturbance subsets. Encountering this "
                "error indicates either a developer oversight or that the "
                "implementation of ComponentSet has been changed. Please open "
                "an issue on the IDAES Github with steps to reproduce this "
                "error."
            )

    return output_sets


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
    for condata in t_block.component_data_objects(ctype=Constraint, active=True):
        sf = _get_scaling_factor(condata)
        if condata.lb != condata.ub:
            raise BurntToast(
                f"Encountered the inequality constraint {condata.name} "
                "while trying to validate that the system is at a steady "
                "state. This should never happen. Please open an issue on "
                "the IDAES Github with steps to reproduce this error."
            )
        if sf * abs(value(condata) - condata.lb) > constraint_tolerance:
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

    .. math:: 0 = f(x, \\dot{x}, u, d)
    .. math:: y = h(x, u, d)

    and turns it into the linear system

    .. math:: \\dot{x} = A x + B u + B_d d
    .. math:: y = C x + D u + D_d d

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
        input_variables (list): List of input (:math:`u`) variables.
        disturbance_variables (list): List of disturbance (:math:`d`)
        variables.
        output_variables (list): List of output (:math:`y`) variables.
        steady_state_derivative_tolerance (float): Tolerance used to
            test whether the time derivative variables are close
            enough to zero to consider that the system is at steady
            state.
        constraint_tolerance (float): Tolerance used to test whether the
            DAE constraints are close enough to being satisfied to consider
            that the system :math:`0 = f(x, xdot, u, d)` is satisfied.

    Returns:
        dict: A dictionary containing the keys:
            * "scaled_jac" (`scipy.sparse.csr_matrix`): Jacobian of scaled
              system evaluated at the point of linearization
            * "nlp" (`PyomoNLP`): NLP object used to evaluate scaled_jac
            * "diff_vars (`list`): The list of differential `VarData` in the
              order they appear in the :math:`A` and :math:`C` matrices
            * "alg_vars" (`list`): The list of algebraic `VarData`
            * "input_vars" (`list`): The list of input `VarData` in the order
              they appear in the :math:`B` and :math:`D` matrices
            * "disturbance_vars (`list`): The list of disturbance `VarData` in
              the order they appear in the :math:`B_d` and :math:`D_d` matrices
            * "output_vars" (`list`): The list of output `VarData` in the order
              they appear in the :math:`C`, :math:`D`, and :math:`D_d` matrices
            * "A" (`numpy.ndarray`): The :math:`A` matrix
            * "B" (`numpy.ndarray`): The :math:`B` matrix
            * "Bd" (`numpy.ndarray`): The :math:`B_d` matrix
            * "C" (`numpy.ndarray`): The :math:`C` matrix
            * "D" (`numpy.ndarray`): The :math:`D` matrix
            * "Dd" (`numpy.ndarray`): The :math:`D_d` matrix
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
        t_block, diff_vars, _ = integrator.get_timestep_problem(t1)
        jac, nlp = get_jacobian(t_block, equality_constraints_only=True)
    finally:
        from_json(model, fixedness, wts=StoreSpec.isfixed())

    vardata_lists["total"] = nlp.get_pyomo_variables()

    vardata_lists["deriv"] = []
    vardata_lists["diff"] = diff_vars

    # Invert the map to get derivative vars in the same order as diff_vars
    diff_deriv_map = ComponentMap()
    for deriv, out in integrator.derivative_differential_vardata_map.items():
        diff_deriv_map[out["diff_var"]] = deriv

    for diff_var in vardata_lists["diff"]:
        vardata_lists["deriv"].append(diff_deriv_map[diff_var])

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
    # ensure a consistent variable order is maintained
    vardata_lists["alg"] = []
    for vardata in vardata_lists["total"]:
        if vardata in vardata_sets["alg"]:
            vardata_lists["alg"].append(vardata)

    output_sets = _validate_vardata_collections(vardata_lists, vardata_sets)

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
        "input_vars": vardata_lists["input"],
        "disturbance_vars": vardata_lists["dist"],
        "output_vars": vardata_lists["output"],
    }

    jac_deriv_alg = jac[:, raw_jac_indices["deriv"] + raw_jac_indices["alg"]]
    jac_rest = jac[
        :, raw_jac_indices["diff"] + raw_jac_indices["input"] + raw_jac_indices["dist"]
    ]

    try:
        sys_raw = spsolve(-jac_deriv_alg.tocsc(), jac_rest.tocsc())
    except RuntimeError as err:
        if "Factor is exactly singular" in repr(err):
            raise RuntimeError(
                "Derivative-algebraic variable jacobian is singular. "
                "This typically indicates a DAE with elevated index "
                "or some other modeling issue."
            ) from err
        else:
            raise

    if hasattr(sys_raw, "todense"):
        sys_raw = sys_raw.todense()
    else:
        # In the edge case where sys_raw is a 1x1 matrix, it's cast to
        # an ndarray and calling todense() would result in an error
        if len(sys_raw.shape) != 2:
            if len(sys_raw) != 1:
                raise RuntimeError(
                    "Unexpected behavior from results returned by spsolve. Please "
                    "open an issue on the IDAES Github with steps to reproduce "
                    "this error."
                )
            else:
                sys_raw = sys_raw.reshape((1, 1))
    if any(isnan(sys_raw)):
        # Should we also check for nans in the jacobian pre-solve?
        raise RuntimeError(
            "nan values encountered when solving for reduced jacobian. "
            "This typically means that the derivative-algebraic variable "
            "jacobian is singular and a modeling problem exists, like the "
            "DAE having an elevated index."
        )

    # We seek an LTI system of the form:
    # xdot = A @ x + B @ u + B_d @ d
    # y = C @ x + D @ u + D_d @ d
    # We don't actually need estimates of the algebraic states in the LTI system,
    # however some outputs (y's) may be algebraic states.
    # Most of these matrices come from indexing and rearranging sys_raw, with a
    # few elements of C being appropriate unit vectors

    nx = len(vardata_lists["diff"])
    # nz = len(vardata_lists["alg"])
    nu = len(vardata_lists["input"])
    nd = len(vardata_lists["dist"])
    ny = len(vardata_lists["output"])

    AB_raw = sys_raw[:nx, :]

    out["A"] = AB_raw[:, :nx]
    out["B"] = AB_raw[:, nx : nx + nu]
    out["Bd"] = AB_raw[:, nx + nu : nx + nu + nd]

    out["C"] = zeros((ny, nx))
    out["D"] = zeros((ny, nu))
    out["Dd"] = zeros((ny, nd))

    for Crow, output, out_idx in zip(
        range(ny), vardata_lists["output"], raw_jac_indices["output"]
    ):
        if output in output_sets["diff"]:
            # Because we've rearranged the rows/columns when creating
            # A, we need to find the *new* index corresponding to the
            # differential variable
            di = raw_jac_indices["diff"].index(out_idx)
            out["C"][Crow, di] = 1
        elif output in output_sets["alg"]:
            ai = len(raw_jac_indices["deriv"]) + raw_jac_indices["alg"].index(out_idx)
            out["C"][Crow, :] = sys_raw[ai, :nx]
            out["D"][Crow, :] = sys_raw[ai, nx : nx + nu]
            out["Dd"][Crow, :] = sys_raw[ai, nx + nu : nx + nu + nd]
        elif output in output_sets["deriv"]:
            # Grab the correct rows of A, B, and Bd
            drvi = raw_jac_indices["deriv"].index(out_idx)
            out["C"][Crow, :] = out["A"][drvi, :]
            out["D"][Crow, :] = out["B"][drvi, :]
            out["Dd"][Crow, :] = out["Bd"][drvi, :]
        elif output in output_sets["input"]:
            # The C row is full of zeros
            ini = raw_jac_indices["input"].index(out_idx)
            out["D"][Crow, ini] = 1
            # The Dd row is full of zeros
        elif output in output_sets["dist"]:
            # The C and D rows are full of zeros
            ind = raw_jac_indices["dist"].index(out_idx)
            out["Dd"][Crow, ind] = 1
        else:
            raise BurntToast(
                "This branch should be inaccessible. Please open an issue on "
                " the IDAES Github with steps to reproduce this exception so "
                "that the IDAES developers can address this problem."
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
        out["A"] = diag(1 / sfxdot) @ out["A"] @ diag(sfx)
        out["B"] = diag(1 / sfxdot) @ out["B"] @ diag(sfu)
        out["Bd"] = diag(1 / sfxdot) @ out["Bd"] @ diag(sfd)

        out["C"] = diag(1 / sfy) @ out["C"] @ diag(sfx)
        out["D"] = diag(1 / sfy) @ out["D"] @ diag(sfu)
        out["Dd"] = diag(1 / sfy) @ out["Dd"] @ diag(sfd)

    return out


def c2d(sys: dict, dt: float):
    r"""Convert a linear time invariant system from continuous time to discrete
    time. Start with a system of the form:

    .. math:: \dot{x} = A x + B u + B_d d

    and create a system

    .. math:: x(t + \Delta t) = \tilde{A} x +  \tilde{B} u + \tilde{B}_d d

    Note that the output matrices :math:`C`, :math:`D`, and :math:`D_d`
    do not need to be converted.

    Args:
        sys (dict): A dictionary containing the keys "A", "B", and "Bd"
            "A" (numpy.ndarray)
            "B" (numpy.ndarray)
            "Bd" (numpy.ndarray)
        dt (float): The timestep  :math:`\Delta t`

    Returns:
        dict: A dictionary containing the keys:
            * A (numpy.ndarray): Discrete time :math:`\tilde{A}` matrix.
            * B (numpy.ndarray): Discrete time :math:`\tilde{B}` matrix.
            * Bd (numpy.ndarray): Discrete time :math:`\tilde{B}_d` matrix.
    """

    A = sys["A"]
    if "B" in sys:
        B = sys["B"]
    else:
        B = zeros((A.shape[0], 0))

    if "Bd" in sys:
        Bd = sys["Bd"]
    else:
        Bd = zeros((A.shape[0], 0))

    nx = A.shape[0]
    if B.shape[0] != nx:
        raise ValueError(f"B matrix has {B.shape[0]} rows, but {nx} were expected.")
    nu = B.shape[1]
    if Bd.shape[0] != nx:
        raise ValueError(f"Bd matrix has {Bd.shape[0]} rows, but {nx} were expected.")
    nd = Bd.shape[1]

    A_aug = block([[A, B, Bd], [zeros((nu + nd, nx + nu + nd))]])

    A_tilde_aug = expm(A_aug * dt)

    out = {}
    out["A"] = A_tilde_aug[:nx, :nx]
    out["B"] = A_tilde_aug[:nx, nx : nx + nu]
    out["Bd"] = A_tilde_aug[:nx, nx + nu : nx + nu + nd]

    return out
