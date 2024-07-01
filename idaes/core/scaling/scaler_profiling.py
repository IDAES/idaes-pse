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
Tools for profiling scaling tools.
"""

from pyomo.environ import check_optimal_termination, Constraint
from pyomo.common.tempfiles import TempfileManager

from idaes.core.util.scaling import jacobian_cond
from idaes.core.scaling import AutoScaler, CustomScalerBase
from idaes.core.solvers import get_solver


# TODO: Needs tests
class ScalingProfiler:
    def __init__(
        self,
        build_model,
        scale_variables,
        perturb_state=None,
        scaling_methods=None,
        solver=None,
    ):
        self._build_model = build_model
        self._scale_variables = scale_variables
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
                "Root Sum Squares": (
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
                    ascaler.constraints_by_jacobian_norm,
                    {"norm": 1, "auto": True},
                ),
                "Actual L2 Norm": (
                    ascaler.constraints_by_jacobian_norm,
                    {"norm": 2, "auto": True},
                ),
            }

    def _scale_vars(self, model, perfect=False):
        if perfect:
            scaler = AutoScaler()
            scaler.variables_by_magnitude(model)
            return model

        self._scale_variables(model)

        return model

    def _collect_data(self, scaling_method, auto, perfect, **kwargs):
        m = self._build_model()
        self._scale_vars(m, perfect=perfect)
        if scaling_method is not None:
            if auto:
                scaling_method(m, **kwargs)
            else:
                for c in m.component_data_objects(ctype=Constraint, descend_into=True):
                    scaling_method(c, **kwargs)

        cond = jacobian_cond(m, scaled=True)

        stats = self._solved_perturbed_state(m)

        return {"condition_number": cond, **stats}

    def run_case(self, scaling_method, **kwargs):
        auto = kwargs.pop("auto", False)

        # Imperfect information
        manual = self._collect_data(scaling_method, perfect=False, auto=auto, **kwargs)

        # Perfect information
        perfect = self._collect_data(scaling_method, perfect=True, auto=auto, **kwargs)

        return {"Manual": manual, "Auto": perfect}

    def _solved_perturbed_state(self, model):
        if self._perturb_state is None:
            return {}

        self._perturb_state(model)

        TempfileManager.push()
        tempfile = TempfileManager.create_tempfile(suffix="ipopt_out", text=True)
        opts = {"output_file": tempfile}

        status_obj = self._solver.solve(model, options=opts, tee=True)
        solved = True
        if not check_optimal_termination(status_obj):
            solved = False

        iters, iters_in_restoration, iters_w_regularization, time = (
            self._parse_ipopt_output(tempfile)
        )

        print(status_obj)

        return {
            "solved": solved,
            "termination_message": status_obj.solver.termination_message,
            "iterations": iters,
            "iters_in_restoration": iters_in_restoration,
            "iters_w_regularization": iters_w_regularization,
        }

    def profile_scaling_methods(self):
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
        time = 0
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
                elif line.startswith(
                    "Total CPU secs in IPOPT (w/o function evaluations)   ="
                ):
                    tokens = line.split()
                    time += float(tokens[9])
                elif line.startswith(
                    "Total CPU secs in NLP function evaluations           ="
                ):
                    tokens = line.split()
                    time += float(tokens[8])

        return iters, iters_in_restoration, iters_w_regularization, time
