#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tools for profiling scaling alternatives.
"""
import sys

from pyomo.environ import check_optimal_termination, Constraint
from pyomo.common.tempfiles import TempfileManager

from idaes.core.util.scaling import jacobian_cond
from idaes.core.scaling import AutoScaler, CustomScalerBase
from idaes.core.solvers import get_solver


class ScalingProfiler:
    """
    Class for running a set of constraint scaling methods on a model and reporting the
    effect on model condition number and solver behavior.

    Users should call the profile_scaling_methods method to generate a dict of results or
    the report_scaling_profiles method for a stream-based output.

    Users are expected to provide callback functions to 1) construct an initialized model
    that will be used for profiling, 2) apply user-defined variable scaling (used for the
    imperfect information case) and 3) perturb the model from the initialized state to
    test how well the model solves (optional).

    Users may also provide a dict of scaling methods they wish to apply using the scaling_methods
    argument. If this is not provided, the tool will default to applying all the scaling methods
    defined by the AutoScaler and CustomScalerBase classes.

    **NOTE** methods from the AutoScaler class are applied to Pyomo Blocks, whilst those from
    CustomScalerBase are applied to individual ConstraintDatas. The profiling tool assumes that
    methods will be applied to ConstraintDatas unless the `block_based` keyword argument is set to True
    for the scaling method.

    Args:
        build_model: callback to use to construct initialized model for testing
        user_scaling: callback to use to apply user-defined scaling to initialized model
        perturb_states: (optional) callback to use to perturb model state for re-solve tests
        scaling_methods: (optional) dict of constraint scaling methods to profile. {"Name": (method, kwargs)}
        solver: (optional) Pyomo solver object to use for re-solve tests
    """

    def __init__(
        self,
        build_model,
        user_scaling,
        perturb_state=None,
        scaling_methods: dict = None,
        solver=None,
    ):
        """
        Sets up a framework for applying different scaling methods to a model and compiling a
        report of their effects on the Jacobian condition number and how easily the scaled
        model can be solved for a perturbed state.
        """
        self._build_model = build_model
        self._user_scaling = user_scaling
        self._perturb_state = perturb_state
        self._scaling_methods = scaling_methods
        self._solver = solver

        if self._solver is None:
            self._solver = get_solver("ipopt_v2", writer_config={"scale_model": True})

        if self._scaling_methods is None:
            ascaler = AutoScaler()
            cscaler = CustomScalerBase()

            self._scaling_methods = {
                "Vars Only": (None, {}),
                "Harmonic": (
                    cscaler.scale_constraint_by_nominal_value,
                    {"scheme": "harmonic_mean"},
                ),
                "Inverse Sum": (
                    cscaler.scale_constraint_by_nominal_value,
                    {"scheme": "inverse_sum"},
                ),
                "Inverse Root Sum Squares": (
                    cscaler.scale_constraint_by_nominal_value,
                    {"scheme": "inverse_root_sum_squared"},
                ),
                "Inverse Maximum": (
                    cscaler.scale_constraint_by_nominal_value,
                    {"scheme": "inverse_maximum"},
                ),
                "Inverse Minimum": (
                    cscaler.scale_constraint_by_nominal_value,
                    {"scheme": "inverse_minimum"},
                ),
                "Nominal L1 Norm": (
                    cscaler.scale_constraint_by_nominal_derivative_norm,
                    {"norm": 1},
                ),
                "Nominal L2 Norm": (
                    cscaler.scale_constraint_by_nominal_derivative_norm,
                    {"norm": 2},
                ),
                "Actual L1 Norm": (
                    ascaler.scale_constraints_by_jacobian_norm,
                    {"norm": 1, "block_based": True},
                ),
                "Actual L2 Norm": (
                    ascaler.scale_constraints_by_jacobian_norm,
                    {"norm": 2, "block_based": True},
                ),
            }

    def profile_scaling_methods(self):
        """
        Generate results for all provided scaling methods.

        For each scaling method, calculate the Jacobian condition number and re-solve
        the model with both user-provided variable scaling and perfect variable scaling
        (scaling by inverse magnitude). A base case with no scaling applied is also
        run for reference.

        Args:
            None

        Returns:
            dict with results for all scaling methods
        """
        # Generate data for unscaled model
        m = self._build_model()
        unscaled = jacobian_cond(m, scaled=False)
        stats = self._solved_perturbed_state(m)

        results = {
            "Unscaled": {"Manual": {"condition_number": unscaled, **stats}},
        }

        # Run other cases
        for case, (meth, margs) in self._scaling_methods.items():
            results[case] = self.run_case(meth, **margs)

        return results

    def report_scaling_profiles(self, stream=None):
        """
        Run scaling profile workflow nad report results to a stream.

        Args:
            stream: StringIO object to write result to (default=stdout)

        Returns:
            None
        """
        results = self.profile_scaling_methods()

        self.write_profile_report(results, stream)

    def write_profile_report(self, results: dict, stream=None):
        """
        Write a report on the comparison of scaling methods to a stream based on
        existing results dict.

        Args:
            results: dict containing results from a scaler profiling run
            stream: StringIO object to write result to (default=stdout)

        Returns:
            None
        """
        # If stream is None, default to stdout
        if stream is None:
            stream = sys.stdout

        # Get length of longest method name for padding
        max_str = max([len(i) for i in results.keys()])
        if max_str < len("Scaling Method"):
            max_str = len("Scaling Method")
        # Length of each stats field is 22 characters, plus 4 for column dividers
        # Max line length is thus longest string name + 2 columns of 22+4 characters
        max_line = max_str + 26 * 2

        # Write header rows
        stream.write(f"\n{'='*max_line}\n")
        stream.write("Scaling Profile Report\n")
        stream.write(f"{'-' * max_line}\n")
        # Pad User Scaling columns to full column width (22)
        stream.write(
            f"{'Scaling Method': <{max_str}} || {'User Scaling': <{22}} || Perfect Scaling\n"
        )

        # Iterate over keys in results and write summary for each scaling method
        for k, v in results.items():
            if v["Manual"]["solved"]:
                msolved = "Solved"
            else:
                msolved = "Failed"
            # Pad iterations to 3 characters - hopefully we don;t see more than 999 iterations
            miters = f"{v['Manual']['iterations']: <{3}}"
            stream.write(
                f"{k: <{max_str}} || {v['Manual']['condition_number']:.3E} | {msolved} {miters} "
            )

            if "Auto" in v.keys():
                if v["Auto"]["solved"]:
                    asolved = "Solved"
                else:
                    asolved = "Failed"
                # Pad iterations again
                aiters = f"{v['Auto']['iterations']: <{3}}"
                stream.write(
                    f"|| {v['Auto']['condition_number']:.3E} | {asolved} {aiters}\n"
                )
            else:
                # Add training column divider but no auto column
                stream.write("||\n")

        # Write footer row
        stream.write(f"{'=' * max_line}\n")

    def run_case(self, scaling_method, **kwargs):
        """
        Run case for a given scaling method with both perfect and imperfect scaling information.

        Args:
            scaling_method: constraint scaling method to be tested
            kwargs: keyword argument to be passed to scaling method

        Returns:
            dict summarising results of scaling case
        """
        block_based = kwargs.pop("block_based", False)

        # Imperfect information
        manual = self._run_scenario(
            scaling_method, perfect=False, block_based=block_based, **kwargs
        )

        # Perfect information
        perfect = self._run_scenario(
            scaling_method, perfect=True, block_based=block_based, **kwargs
        )

        return {"Manual": manual, "Auto": perfect}

    def _scale_vars(self, model, perfect=False):
        """
        Apply variable scaling to model.
        """
        if perfect:
            scaler = AutoScaler()
            scaler.scale_variables_by_magnitude(model)
            return model

        self._user_scaling(model)

        return model

    def _apply_scaling(self, model, scaling_method, block_based, **kwargs):
        """
        Collect stats for a given scaling method.
        """
        if scaling_method is not None:
            if block_based:
                scaling_method(model, **kwargs)
            else:
                for c in model.component_data_objects(
                    ctype=Constraint, descend_into=True
                ):
                    scaling_method(c, **kwargs)

    def _run_scenario(self, scaling_method, block_based, perfect, **kwargs):
        """
        Run a single scenario
        """
        m = self._build_model()
        self._scale_vars(m, perfect=perfect)
        self._apply_scaling(m, scaling_method, block_based=block_based, **kwargs)

        cond = jacobian_cond(m, scaled=True)
        stats = self._solved_perturbed_state(m)

        return {"condition_number": cond, **stats}

    def _solved_perturbed_state(self, model):
        """
        Run re-solve tests if perturb_state callback provided.
        """
        if self._perturb_state is None:
            return {}

        self._perturb_state(model)

        TempfileManager.push()
        tempfile = TempfileManager.create_tempfile(suffix="ipopt_out", text=True)
        opts = {"output_file": tempfile}

        try:
            status_obj = self._solver.solve(model, options=opts, tee=True)
            solved = True
            if not check_optimal_termination(status_obj):
                solved = False

            iters, iters_in_restoration, iters_w_regularization = (
                self._parse_ipopt_output(tempfile)
            )

            return {
                "solved": solved,
                "termination_message": status_obj.solver.termination_message,
                "iterations": iters,
                "iters_in_restoration": iters_in_restoration,
                "iters_w_regularization": iters_w_regularization,
            }
        except RuntimeError as err:
            # Likely a critical solver failure
            return {
                "solved": False,
                "termination_message": str(err),
                "iterations": -1,
                "iters_in_restoration": -1,
                "iters_w_regularization": -1,
            }

    def _parse_ipopt_output(self, ipopt_file):
        """
        Parse an IPOPT output file and return:

        * number of iterations
        * time in IPOPT

        Returns
        -------
           Returns a tuple with (solve status object, bool (solve successful or
           not), number of iters, solve time)
        """
        # ToDO: Check for final iteration with regularization or restoration

        iters = 0
        iters_in_restoration = 0
        iters_w_regularization = 0
        # parse the output file to get the iteration count, solver times, etc.
        with open(ipopt_file, "r") as f:
            parseline = False
            for line in f:
                if line.startswith("iter"):
                    # This marks the start of the iteration logging, set parseline True
                    parseline = True
                elif line.startswith("Number of Iterations....:"):
                    # Marks end of iteration logging, set parseline False
                    parseline = False
                    tokens = line.split()
                    iters = int(tokens[3])
                elif parseline:
                    # Line contains details of an iteration, look for restoration or regularization
                    tokens = line.split()
                    try:
                        if not tokens[6] == "-":
                            # Iteration with regularization
                            iters_w_regularization += 1
                        if tokens[0].endswith("r"):
                            # Iteration in restoration
                            iters_in_restoration += 1
                    except IndexError:
                        # Blank line at end of iteration list, so assume we hit this
                        pass

        return iters, iters_in_restoration, iters_w_regularization
