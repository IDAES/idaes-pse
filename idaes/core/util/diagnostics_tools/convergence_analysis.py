# -*- coding: utf-8 -*-
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
This module contains a tool for analyzing convergence characteristics of a model using IPOPT.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"

import json
import sys

from math import isclose

from pyomo.environ import (
    check_optimal_termination,
    SolverFactory,
)
from pyomo.core.base.block import BlockData
from pyomo.common.config import (
    ConfigDict,
    ConfigValue,
)
from pyomo.common.tempfiles import TempfileManager

from idaes.core.util.model_statistics import (
    greybox_block_set,
)
from idaes.core.util.parameter_sweep import (
    SequentialSweepRunner,
    ParameterSweepBase,
    is_psweepspec,
)
from idaes.core.util.diagnostics_tools.diagnostics_toolbox import (
    DiagnosticsToolbox,
)


def psweep_runner_validator(val):
    """Domain validator for Parameter Sweep runners

    Args:
        val : value to be checked

    Returns:
        TypeError if val is not a valid callback
    """
    if issubclass(val, ParameterSweepBase):
        return val

    raise ValueError("Workflow runner must be a subclass of ParameterSweepBase.")


CACONFIG = ConfigDict()
CACONFIG.declare(
    "input_specification",
    ConfigValue(
        domain=is_psweepspec,
        doc="ParameterSweepSpecification object defining inputs to be sampled",
    ),
)
CACONFIG.declare(
    "workflow_runner",
    ConfigValue(
        default=SequentialSweepRunner,
        domain=psweep_runner_validator,
        doc="Parameter sweep workflow runner",
    ),
)
CACONFIG.declare(
    "solver_options",
    ConfigValue(
        domain=None,
        description="Options to pass to IPOPT.",
    ),
)
CACONFIG.declare(
    "halt_on_error",
    ConfigValue(
        default=False,
        domain=bool,
        doc="Whether to halt execution of parameter sweep on encountering a solver error (default=False).",
    ),
)


class IpoptConvergenceAnalysis:
    """
    Tool to perform a parameter sweep of model checking for numerical issues and
    convergence characteristics. Users may specify an IDAES ParameterSweep class to
    perform the sweep (default is SequentialSweepRunner).
    """

    CONFIG = CACONFIG()

    def __init__(self, model, solver_obj=None, **kwargs):
        # TODO: In future may want to generalise this to accept indexed blocks
        # However, for now some of the tools do not support indexed blocks
        if not isinstance(model, BlockData):
            raise TypeError(
                "model argument must be an instance of a Pyomo BlockData object "
                "(either a scalar Block or an element of an indexed Block)."
            )
        if len(greybox_block_set(model)) != 0:
            raise NotImplementedError(
                "Model contains Greybox models, which are not supported by Diagnostics toolbox at the moment"
            )
        self.config = self.CONFIG(kwargs)

        self._model = model

        if solver_obj is None:
            solver_obj = SolverFactory("ipopt")
        self._solver_obj = solver_obj

        if self.config.solver_options is not None:
            solver_obj.options = self.config.solver_options

        self._psweep = self.config.workflow_runner(
            input_specification=self.config.input_specification,
            build_model=self._build_model,
            rebuild_model=True,
            run_model=self._run_model,
            build_outputs=self._build_outputs,
            halt_on_error=self.config.halt_on_error,
            handle_solver_error=self._recourse,
            solver=solver_obj,
        )

    @property
    def results(self):
        """
        Returns the results of the IpoptConvergenceAnalysis run
        """
        return self._psweep.results

    @property
    def samples(self):
        """
        Returns the set of input samples for convergence analysis (pandas DataFrame)
        """
        return self._psweep.get_input_samples()

    def run_convergence_analysis(self):
        """
        Execute convergence analysis sweep by calling execute_parameter_sweep
        method in specified runner.

        Returns:
            dict of results from parameter sweep
        """
        return self._psweep.execute_parameter_sweep()

    def run_convergence_analysis_from_dict(self, input_dict: dict):
        """
        Execute convergence analysis sweep using specification defined in dict.

        Args:
            input_dict: dict to load specification from

        Returns:
            dict of results from parameter sweep
        """
        self.from_dict(input_dict)
        return self.run_convergence_analysis()

    def run_convergence_analysis_from_file(self, filename: str):
        """
        Execute convergence analysis sweep using specification defined in json file.

        Args:
            filename: name of file to load specification from as string

        Returns:
            dict of results from parameter sweep
        """
        self.from_json_file(filename)
        return self.run_convergence_analysis()

    def compare_convergence_to_baseline(
        self, filename: str, rel_tol: float = 0.1, abs_tol: float = 1
    ):
        """
        Run convergence analysis and compare results to those defined in baseline file.

        Args:
            filename: name of baseline file to load specification from as string
            rel_tol: relative tolerance to use for comparing number of iterations
            abs_tol: absolute tolerance to use for comparing number of iterations

        Returns:
            dict containing lists of differences between convergence analysis run and baseline
        """
        with open(filename, "r") as f:
            # Load file manually so we have saved results
            jdict = json.load(f)
        f.close()

        # Run convergence analysis from dict
        self.run_convergence_analysis_from_dict(jdict)

        # Compare results
        return self._compare_results_to_dict(jdict, rel_tol=rel_tol, abs_tol=abs_tol)

    def assert_baseline_comparison(
        self, filename: str, rel_tol: float = 0.1, abs_tol: float = 1
    ):
        """
        Run convergence analysis and assert no differences in results to those defined
        in baseline file.

        Args:
            filename: name of baseline file to load specification from as string
            rel_tol: relative tolerance to use for comparing number of iterations
            abs_tol: absolute tolerance to use for comparing number of iterations

        Raises:
            AssertionError if results of convergence analysis do not match baseline
        """
        diffs = self.compare_convergence_to_baseline(
            filename, rel_tol=rel_tol, abs_tol=abs_tol
        )

        if any(len(v) != 0 for v in diffs.values()):
            raise AssertionError("Convergence analysis does not match baseline")

    def report_convergence_summary(self, stream=None):
        """
        Reports a brief summary of the model convergence run.

        Args:
            stream: Optional output stream to print results to.

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        successes = 0
        failures = 0
        runs_w_restoration = 0
        runs_w_regulariztion = 0
        runs_w_num_iss = 0

        for v in self.results.values():
            # Check for differences
            if v["success"]:
                successes += 1
            else:
                failures += 1

            if v["results"]["iters_in_restoration"] > 0:
                runs_w_restoration += 1
            if v["results"]["iters_w_regularization"] > 0:
                runs_w_regulariztion += 1
            if v["results"]["numerical_issues"] > 0:
                runs_w_num_iss += 1

        stream.write(
            f"Successes: {successes}, Failures {failures} ({100*successes/(successes+failures)}%)\n"
        )
        stream.write(f"Runs with Restoration: {runs_w_restoration}\n")
        stream.write(f"Runs with Regularization: {runs_w_regulariztion}\n")
        stream.write(f"Runs with Numerical Issues: {runs_w_num_iss}\n")

    def to_dict(self):
        """
        Serialize specification and current results to dict form

        Returns:
            dict
        """
        return self._psweep.to_dict()

    def from_dict(self, input_dict):
        """
        Load specification and results from dict.

        Args:
            input_dict: dict to load from

        Returns:
            None
        """
        return self._psweep.from_dict(input_dict)

    def to_json_file(self, filename):
        """
        Write specification and results to json file.

        Args:
            filename: name of file to write to as string

        Returns:
            None
        """
        return self._psweep.to_json_file(filename)

    def from_json_file(self, filename):
        """
        Load specification and results from json file.

        Args:
            filename: name of file to load from as string

        Returns:
            None
        """
        return self._psweep.from_json_file(filename)

    def _build_model(self):
        # Create new instance of model by cloning
        return self._model.clone()

    def _run_model(self, model, solver):
        # Run model using IPOPT and collect stats
        (
            status,
            iters,
            iters_in_restoration,
            iters_w_regularization,
            time,
        ) = self._run_ipopt_with_stats(model, solver)

        run_stats = [
            iters,
            iters_in_restoration,
            iters_w_regularization,
            time,
        ]

        success = check_optimal_termination(status)

        return success, run_stats

    @staticmethod
    def _build_outputs(model, run_stats):
        # Run model diagnostics numerical checks
        dt = DiagnosticsToolbox(model=model)

        warnings = False
        try:
            dt.assert_no_numerical_warnings()
        except AssertionError:
            warnings = True

        # Compile Results
        return {
            "iters": run_stats[0],
            "iters_in_restoration": run_stats[1],
            "iters_w_regularization": run_stats[2],
            "time": run_stats[3],
            "numerical_issues": warnings,
        }

    @staticmethod
    def _recourse(model):
        # Return a default dict indicating no results
        return {
            "iters": -1,
            "iters_in_restoration": -1,
            "iters_w_regularization": -1,
            "time": -1,
            "numerical_issues": -1,
        }

    @staticmethod
    def _parse_ipopt_output(ipopt_file):
        # PArse IPOPT logs and return key metrics
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

    def _run_ipopt_with_stats(self, model, solver, max_iter=500, max_cpu_time=120):
        # Solve model using provided solver (assumed to be IPOPT) and parse logs
        # ToDo: Check that the "solver" is, in fact, IPOPT

        TempfileManager.push()
        tempfile = TempfileManager.create_tempfile(suffix="ipopt_out", text=True)
        opts = {
            "output_file": tempfile,
            "max_iter": max_iter,
            "max_cpu_time": max_cpu_time,
        }

        status_obj = solver.solve(model, options=opts, tee=True)

        (
            iters,
            iters_in_restoration,
            iters_w_regularization,
            time,
        ) = self._parse_ipopt_output(tempfile)

        TempfileManager.pop(remove=True)
        return status_obj, iters, iters_in_restoration, iters_w_regularization, time

    def _compare_results_to_dict(
        self,
        compare_dict: dict,
        rel_tol: float = 0.1,
        abs_tol: float = 1,
    ):
        # Compare results
        success_diff = []
        iters_diff = []
        restore_diff = []
        reg_diff = []
        num_iss_diff = []

        for k, v in self.results.items():
            # Get sample result from compare dict
            try:
                comp = compare_dict["results"][k]
            except KeyError:
                # Reading from json often converts ints to strings
                # Check to see if the index as a string works
                comp = compare_dict["results"][str(k)]

            # Check for differences
            if v["success"] != comp["success"]:
                success_diff.append(k)
            if not isclose(
                v["results"]["iters"],
                comp["results"]["iters"],
                rel_tol=rel_tol,
                abs_tol=abs_tol,
            ):
                iters_diff.append(k)
            if not isclose(
                v["results"]["iters_in_restoration"],
                comp["results"]["iters_in_restoration"],
                rel_tol=rel_tol,
                abs_tol=abs_tol,
            ):
                restore_diff.append(k)
            if not isclose(
                v["results"]["iters_w_regularization"],
                comp["results"]["iters_w_regularization"],
                rel_tol=rel_tol,
                abs_tol=abs_tol,
            ):
                reg_diff.append(k)
            if v["results"]["numerical_issues"] != comp["results"]["numerical_issues"]:
                num_iss_diff.append(k)

        return {
            "success": success_diff,
            "iters": iters_diff,
            "iters_in_restoration": restore_diff,
            "iters_w_regularization": reg_diff,
            "numerical_issues": num_iss_diff,
        }
