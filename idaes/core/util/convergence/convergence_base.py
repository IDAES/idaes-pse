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
"""
This module provides the base classes and methods for running convergence
evaluations on IDAES models. The convergence evaluation runs a given model over
a set of sample points to ensure reliable convergence over the parameter space.

The module requires the user to provide:
  - a set of inputs along with their lower bound, upper bound, mean,
and standard deviation.
  - an initialized Pyomo model
  - a Pyomo solver with appropriate options

The module executes convergence evaluation in two steps. In the first step, a
json file is created that containsa set of points sampled from the provided
inputs. This step only needs to be done once - up front. The second step, which
should be executed any time there is a major code change that could impact the
model, takes that set of sampled points and solves the model at each of the
points, collecting convergence statistics (success/failure, iterations, and
solution time).

This can be used as a tool to evaluate model convergence reliability over the
defined input space, or to verify that convergence performance is not
decreasing with framework and/or model changes.

In order to write a convergence evaluation for your model, you must inherit a
class from ConvergenceEvaluation, and implement three methods:

- get_specification: This method should create and return a
    ConvergenceEvaluationSpecification object. There are methods on
    ConvergenceEvaluationSpecification to add inputs. These inputs contain a
    string that identifies a Pyomo Param or Var object, the lower and upper
    bounds, and the mean and standard deviation to be used for sampling. When
    samples are generated, they are drawn from a normal distribution, and then
    truncated by the lower or upper bounds.
- get_initialized_model: This method should create and return a Pyomo model
    object that is already initialized and ready to be solved. This model will
    be modified according to the sampled inputs, and then it will be solved.
- get_solver: This method should return an instance of the Pyomo solver that
    will be used for the analysis.

There are methods to create the sample points file (on
ConvergenceEvaluationSpecification), to run a convergence evaluation
(run_convergence_evaluation), and print the results in table form
(print_convergence_statistics).

However, this package can also be executed using the command-line interface.
See the documentation in convergence.py for more information.
"""
# stdlib
from collections import OrderedDict
import getpass
import importlib as il
import json
import logging
import numpy as np
import sys
from io import StringIO

# pyomo
from pyomo.common.tempfiles import TempfileManager
from pyomo.common.tee import capture_output
from pyomo.core import Param, Var
from pyomo.opt import TerminationCondition
from pyomo.common.log import LoggingIntercept

# idaes
import idaes.core.util.convergence.mpi_utils as mpiu
from idaes.core.dmf import resource
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


convergence_classes = {}


def register_convergence_class(name):
    def _register_convergence_class(cls):
        if name in convergence_classes:
            raise KeyError(f"Convergence class {name} already registered.")
        convergence_classes[name] = ".".join([cls.__module__, cls.__name__])
        return cls

    return _register_convergence_class


class ConvergenceEvaluationSpecification(object):
    def __init__(self):
        self._inputs = OrderedDict()

    def add_sampled_input(
        self, name, pyomo_path, lower, upper, mean=None, std=None, distribution="normal"
    ):
        """
        Add an input that should be sampled when forming the set of
        specific points that need to be run for the convergence evaluation

        The input will be sampled assuming a normal distribution
        (with given mean and standard devation)
        truncated to the values given by lower and upper bounds

        Parameters
        ----------
        name : str
           The name of the input.
        pyomo_path : str
           A string representation of the path to the variable or parameter to
           be sampled. This string will be executed to retrieve the Pyomo
           component.
        lower : float
           A lower bound on the input variable or parameter.
        upper : float
           An upper bound on the input variable or parameter
        mean : float
           The mean value to use when generating normal distribution samples
        std : float
           The standard deviation to use when generating normal distribution samples
        distribution : str
            The Distribution type {"normal", "uniform"}

        Returns
        -------
           N/A
        """
        # ToDo: put some error checking here ... Maybe we should have the model
        # ToDo: already? Can use to check if the pyomo_path is valid? check if
        # ToDo: bounds will be violated?
        spec = OrderedDict()
        spec["pyomo_path"] = pyomo_path
        spec["lower"] = lower
        spec["upper"] = upper
        spec["mean"] = mean
        spec["std"] = std
        spec["distribution"] = distribution
        self._inputs[name] = spec

    @property
    def inputs(self):
        return self._inputs


class ConvergenceEvaluation(object):
    def __init__(self):
        self._sampling_specifications = list()

    def get_specification(self):
        """
        User should override this method to return an instance of the
        ConvergenceEvaluationSpecification for this particular model and test
        set.

        The basic flow for this method is:
           - Create a ConvergenceEvaluationSpecification
           - Call add_sampled_input for every input that should be varied.
           - return the ConvergenceEvaluationSpecification

        Returns
        -------
           ConvergenceEvaluationSpecification
        """
        raise NotImplementedError(
            "Not implemented in the base class. This"
            " should be overridden in the derived class"
        )

    def get_initialized_model(self):
        """
        User should override this method to return an initialized model that is
        ready to solve. The convergence evaluation methods will change the
        values of parameters or variables according to the sampling
        specifications.

        Returns
        -------
           Pyomo model : return a Pyomo model object that is initialized and
                        ready to solve. This is the model object that will be
                        used in the evaluation.
        """
        raise NotImplementedError(
            "Not implemented in the base class. This"
            " should be overridden in the derived class"
        )

    def get_solver(self):
        """
        User should create and return the solver that will be used for the
        convergence evaluation (including any necessary options)

        Returns
        -------
           Pyomo solver

        """
        # ToDo: We may want this to be standardized across all the models
        # within IDAES (should not need different solver options for different
        # unit models)
        raise NotImplementedError(
            "Not implemented in the base class. This"
            " should be overridden in the derived class"
        )


def _class_import(class_path):
    # Note, this method assumes that everything in front
    # of the last dot is a module, followed by one class
    # I don't think this will work for nested classes
    tokens = class_path.split(".")
    modpath = ".".join(tokens[0:-1])
    mod = il.import_module(modpath)
    ret_class = getattr(mod, tokens[-1])
    return ret_class


def _run_ipopt_with_stats(model, solver, max_iter=500, max_cpu_time=120):
    """
    Run the solver (must be ipopt) and return the convergence statistics

    Parameters
    ----------
    model : Pyomo model
       The pyomo model to be solved

    solver : Pyomo solver
       The pyomo solver to use - it must be ipopt, but with whichever options
       are preferred

    max_iter : int
       The maximum number of iterations to allow for ipopt

    max_cpu_time : int
       The maximum cpu time to allow for ipopt (in seconds)

    Returns
    -------
       Returns a tuple with (solve status object, bool (solve successful or
       not), number of iters, solve time)
    """
    # ToDo: Check that the "solver" is, in fact, IPOPT

    TempfileManager.push()
    tempfile = TempfileManager.create_tempfile(suffix="ipopt_out", text=True)
    opts = {"output_file": tempfile, "max_iter": max_iter, "max_cpu_time": max_cpu_time}

    status_obj = solver.solve(model, options=opts, tee=True)
    solved = True
    if status_obj.solver.termination_condition != TerminationCondition.optimal:
        solved = False

    iters = 0
    time = 0
    # parse the output file to get the iteration count, solver times, etc.
    with open(tempfile, "r") as f:
        for line in f:
            if line.startswith("Number of Iterations....:"):
                tokens = line.split()
                iters = int(tokens[3])
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

    TempfileManager.pop(remove=True)
    return status_obj, solved, iters, time


def _progress_bar(fraction, msg, length=20):
    length = length - 2
    n_complete = int(length * fraction)
    n_remaining = length - n_complete
    characters = ["*"] * n_complete  # ['*' for i in range(n_complete)]
    characters.extend(["-"] * n_remaining)  # ['-' for i in range(n_remaining)])
    sys.stdout.write(
        "%5.1f%s [%s] %s\n" % (fraction * 100.0, "%", "".join(characters), msg)
    )
    sys.stdout.flush()


def _set_model_parameters_from_sample(model, inputs, sample_point):
    """
    This method takes the parameter values in sample_points and sets them on
    the model to prepare it for the new solve

    Parameters
    ----------
    model : Pyomo model
       The pyomo model that will be modified with values from the new
       sample_points

    inputs : dict
       The inputs dictionary from the sample-file

    sample_point : dict
       The dictionary of var and/or parameter values for the sample point

    Returns
    -------
       N/A
    """
    for k, v in sample_point.items():
        if k == "_name":
            # because parallel task manager does not handle dictionaries,
            # we add an _name entry to the sample point dictionary
            # ToDo: remove this when parallel task manager is fixed
            continue
        # k stores the "name" of the input
        # need to get the pyomo path of the input from the inputs structure
        pyomo_path = inputs[k]["pyomo_path"]

        comp = model.find_component(pyomo_path)
        try:
            ctype = comp.ctype
        except AttributeError:
            ctype = None

        if ctype is Param:
            if comp.is_constant():
                raise ValueError(
                    "Convergence testing found an input of type Param that"
                    "was not mutable. Please make sure all sampled inputs"
                    " are either mutable params or fixed vars."
                )
            comp.set_value(v)
        elif ctype is Var:
            if not comp.is_fixed():
                raise ValueError(
                    "Convergence testing found an input of type Var that"
                    "was not fixed. Please make sure all sampled inputs"
                    " are either mutable params or fixed vars."
                )
            comp.set_value(float(v))
        else:
            raise ValueError(
                "Failed to find a valid input component (must be"
                " a fixed Var or a mutable Param). Instead,"
                " pyomo_path: {} returned: {}".format(pyomo_path, comp)
            )


def write_sample_file(
    eval_spec, filename, convergence_evaluation_class_str, n_points, seed=None
):
    """
    Samples the space of the inputs defined in the eval_spec, and creates a
    json file with all the points to be used in executing a convergence
    evaluation

    Parameters
    ----------
    filename : str
       The filename for the json file that will be created containing all the
       points to be run
    eval_spec : ConvergenceEvaluationSpecification
       The convergence evaluation specification object that we would like to
       sample
    convergence_evaluation_class_str : str
       Python string that identifies the convergence evaluation class for this
       specific evaluation. This is usually in the form of module.class_name.
    n_points : int
       The total number of points that should be created
    seed : int or None
       The seed to be used when generating samples. If set to None, then the
       seed is not set
    Returns
    -------
       N/A
    """
    if seed is not None:
        np.random.seed(seed)

    # build the samples
    samples = OrderedDict()
    for i in range(n_points):
        sample = samples["Sample-{}".format(i + 1)] = OrderedDict()
        for k, v in eval_spec.inputs.items():
            if v["distribution"] == "normal":
                s = np.random.normal(loc=v["mean"], scale=v["std"])
                s = v["lower"] if s < v["lower"] else s
                s = v["upper"] if s > v["upper"] else s
            elif v["distribution"] == "uniform":
                s = np.random.uniform(low=v["lower"], high=v["upper"])
            sample[k] = s

    # create the dictionary storing all the necessary information
    jsondict = OrderedDict()
    jsondict["inputs"] = OrderedDict(eval_spec.inputs)
    jsondict["n_points"] = len(samples)
    jsondict["seed"] = seed
    jsondict["samples"] = samples
    jsondict["convergence_evaluation_class_str"] = convergence_evaluation_class_str

    with open(filename, "w") as fd:
        json.dump(jsondict, fd, indent=3)


def run_convergence_evaluation_from_sample_file(sample_file):
    # load the sample file
    try:
        with open(sample_file, "r") as fd:
            jsondict = json.load(fd, object_pairs_hook=OrderedDict)
    except Exception as e:
        _log.exception(f"Error reading json file {sample_file}")
        raise ValueError(f"Problem reading json file: {sample_file}")

    # create the convergence evaluation object
    convergence_evaluation_class_str = jsondict["convergence_evaluation_class_str"]
    try:
        convergence_evaluation_class = _class_import(convergence_evaluation_class_str)
        conv_eval = convergence_evaluation_class()
    except Exception as e:
        _log.exception(f"Error creating class: {convergence_evaluation_class_str}")
        raise ValueError(
            f"Invalid value specified for convergence_evaluation_class_str:"
            "{convergence_evaluation_class_str} in sample file: {sample_file}"
        )
    return run_convergence_evaluation(jsondict, conv_eval)


def run_single_sample_from_sample_file(sample_file, name):
    # load the sample file
    try:
        with open(sample_file, "r") as fd:
            jsondict = json.load(fd, object_pairs_hook=OrderedDict)
    except Exception as e:
        _log.exception(f"Error reading json file {sample_file}")
        raise ValueError(f"Problem reading json file: {sample_file}")

    # create the convergence evaluation object
    convergence_evaluation_class_str = jsondict["convergence_evaluation_class_str"]
    try:
        convergence_evaluation_class = _class_import(convergence_evaluation_class_str)
        conv_eval = convergence_evaluation_class()
    except Exception as e:
        _log.exception(f"Error creating class: {convergence_evaluation_class_str}")
        raise ValueError(
            f"Invalid value specified for convergence_evaluation_class_str:"
            "{convergence_evaluation_class_str} in sample file: {sample_file}"
        )
    return run_single_sample(jsondict, conv_eval, name)


def run_single_sample(sample_file_dict, conv_eval, name):
    inputs = sample_file_dict["inputs"]
    samples = sample_file_dict["samples"]
    model = conv_eval.get_initialized_model()
    _set_model_parameters_from_sample(model, inputs, sample_file_dict["samples"][name])
    solver = conv_eval.get_solver()
    return _run_ipopt_with_stats(model, solver)


def run_convergence_evaluation(sample_file_dict, conv_eval):
    """
    Run convergence evaluation and generate the statistics based on information
    in the sample_file.

    Parameters
    ----------
    sample_file_dict : dict
        Dictionary created by ConvergenceEvaluationSpecification that contains
        the input and sample point information

    conv_eval : ConvergenceEvaluation
        The ConvergenceEvaluation object that should be used

    Returns
    -------
       N/A
    """
    inputs = sample_file_dict["inputs"]
    samples = sample_file_dict["samples"]

    # current parallel task manager code does not work with dictionaries, so
    # convert samples to a list
    # ToDo: fix and test parallel task manager with dictionaries and change
    # this
    samples_list = list()
    for k, v in samples.items():
        v["_name"] = k
        samples_list.append(v)
    n_samples = len(samples_list)

    task_mgr = mpiu.ParallelTaskManager(n_samples)
    local_samples_list = task_mgr.global_to_local_data(samples_list)

    results = list()
    for (si, ss) in enumerate(local_samples_list):
        sample_name = ss["_name"]
        # print progress on the rank-0 process
        if task_mgr.is_root():
            _progress_bar(
                float(si) / float(len(local_samples_list)),
                "Root Process: {}".format(sample_name),
            )

        # capture the output
        # ToDo: make this an option and turn off for single sample execution
        output_buffer = StringIO()
        with LoggingIntercept(output_buffer, "idaes", logging.ERROR):
            with capture_output():  # as str_out:
                model = conv_eval.get_initialized_model()
                _set_model_parameters_from_sample(model, inputs, ss)
                solver = conv_eval.get_solver()
                (status_obj, solved, iters, time) = _run_ipopt_with_stats(model, solver)

        if not solved:
            _log.error(f"Sample: {sample_name} failed to converge.")

        results_dict = OrderedDict()
        results_dict["name"] = sample_name
        results_dict["sample_point"] = ss
        results_dict["solved"] = solved
        results_dict["iters"] = iters
        results_dict["time"] = time
        results.append(results_dict)

    global_results = task_mgr.gather_global_data(results)
    return inputs, samples, global_results


def save_convergence_statistics(
    inputs, results, dmf=None, display=True, json_path=None, report_path=None
):
    """ """
    s = Stats(inputs, results)
    if display:
        s.report()
    if report_path:
        with open(report_path, "w") as f:
            s.report(f)
    if json_path is not None:
        with open(json_path, "w") as f:
            s.to_json(f)
    if dmf is not None:
        s.to_dmf(dmf)
    return s


class Stats(object):
    def __init__(self, inputs=None, results=None, from_dict=None, from_json=None):
        """A convergence stats and results object.  This class stores the
        convergence test evaluation results and generates reports and storage
        formats.

        Args:
            inputs: sample inputs
            results: convergence evaluation results
            from_dict: load a stats object from a given dict
            from_json: load a stats object from a json file path
        """
        # Reload from a dict or json file.  This can be used to compare results
        # or generate new reports from.
        if from_json:
            with open(from_json, "r") as f:
                from_dict = json.load(f)
        if from_dict:
            self.from_dict(from_dict)
            return
        assert inputs is not None
        assert results is not None

        self.inputs = inputs
        self.results = results
        self.notable_cases, self.failed_cases = [], []
        self.iters_successful = list()
        self.time_successful = list()

        # loop through and gather some data
        for r in results:
            if r["solved"] is True:
                self.iters_successful.append(r["iters"])
                self.time_successful.append(r["time"])
            else:
                self.failed_cases.append(r)
        # data for summary table
        self.iters_min = 0
        self.iters_mean = 0
        self.iters_std = 0
        self.iters_max = 0
        if len(self.iters_successful) > 0:
            self.iters_min = int(np.min(self.iters_successful))
            self.iters_mean = int(np.mean(self.iters_successful))
            self.iters_std = int(np.std(self.iters_successful))
            self.iters_max = int(np.max(self.iters_successful))
        self.time_min = 0
        self.time_mean = 0
        self.time_std = 0
        self.time_max = 0
        if len(self.time_successful) > 0:
            self.time_min = float(np.min(self.time_successful))
            self.time_mean = float(np.mean(self.time_successful))
            self.time_std = float(np.std(self.time_successful))
            self.time_max = float(np.max(self.time_successful))

        for r in results:
            flag = ""
            if r["solved"] is not True:
                flag = "F"
            else:
                if (
                    r["time"] > self.time_mean + 2.0 * self.time_std
                    and r["time"] > self.time_mean + 5.0
                ):
                    # add a more absolute check when s.time_std is small
                    flag += "T"
                if (
                    r["iters"] > self.iters_mean + 2.0 * self.iters_std
                    and r["iters"] > self.iters_mean + 5
                ):
                    # add a more absolute check when s.iters_std is small
                    flag += "I"
            r["flag"] = flag
            if flag != "":
                self.notable_cases.append(r)

    def from_dict(self, d):
        for k, v in d.items():
            setattr(self, k, v)

    def to_dict(self):
        keys = [
            a
            for a in dir(self)
            if not a.startswith("_") and not callable(getattr(self, a))
        ]
        d = {}
        for k in keys:
            d[k] = getattr(self, k)
        return d

    def to_json(self, fp):
        json.dump(self.to_dict(), fp, indent=4)

    def to_dmf(self, dmf):
        # PYLINT-TODO-FIX fix error due to undefined variable "stats"
        rsrc = resource.Resource(
            value={
                "name": "convergence_results",
                "desc": "statistics returned from run_convergence_evaluation",
                "creator": {"name": getpass.getuser()},
                # pylint: disable=undefined-variable
                "data": stats.to_dict(),
            },
            type_=resource.ResourceTypes.data,
        )
        # pylint: enable=undefined-variable
        dmf.add(rsrc)

    def report(self, fp=sys.stdout):
        s = self
        res = self.results
        n = len(res)
        fp.write(f"\n{'='*24}{'Scenario Statistics':^24s}{'='*24}\n\n")
        fp.write(
            f"{'Parameter':>20s}{'Min':>10s}{'Mean':>10s}"
            f"{'Stdev':>10s}{'Max':>10s}\n"
        )
        fp.write(f"{'-'*60}\n")
        for k, v in self.inputs.items():
            values = [res[i]["sample_point"][k] for i in range(n)]
            fp.write(
                f"{k:>20s}{float(np.min(values)):10.3g}"
                f"{float(np.mean(values)):10.3g}"
                f"{float(np.std(values)):10.3g}"
                f"{float(np.max(values)):10.3}\n"
            )
        fp.write(f"{'-'*60}\n\n")
        fp.write(f"\n{'='*24}{'Summary':^24s}{'='*24}\n\n")
        nsuc = n - len(s.failed_cases)
        fp.write(f"Number of Successful Cases (solved=True): {nsuc}/{n}\n\n")
        fp.write(
            f"{'':20s}{'min':>10s}{'-1std':>10s}{'mean':>10s}"
            f"{'+1std':>10s}{'max':>10s}\n"
        )
        fp.write(f"{'-'*70}\n")
        fp.write(
            f"{'Iterations':>20s}{s.iters_min:10.3g}"
            f"{(s.iters_mean - s.iters_std):10.3g}"
            f"{s.iters_mean:10.3g}{(s.iters_mean + s.iters_std):10.3g}"
            f"{s.iters_max:10.3g}\n"
        )
        fp.write(
            f"{'Solver Time (s)':>20s}{s.time_min:10.3g}"
            f"{(s.time_mean - s.time_std):10.3g}"
            f"{s.time_mean:10.3g}{(s.time_mean + s.time_std):10.3g}"
            f"{s.time_max:10.3g}\n"
        )
        fp.write(f"{'-'*70}\n\n")
        # print the detailed table
        fp.write(f"\n{'='*24}{'Table of Results':^24s}{'='*24}\n\n")
        fp.write(
            f"{'Flag':>5s}{'Name':>20s}{'Solved':>10s}" f"{'Iters':>10s}{'Time':>10s}\n"
        )
        fp.write(f"{'-'*55}\n")
        for r in res:
            fp.write(
                f"{r['flag']:>5s}{r['name']:>20s}{str(r['solved']):>10s}"
                f"{r['iters']:>10d}{r['time']:>10.2f}\n"
            )
        fp.write(f"{'-'*55}\n\n")
        fp.write(f"\n{'='*24}{'Notable Cases':^24s}{'='*24}\n\n")
        if len(s.notable_cases) == 0:
            fp.write("... None\n")
        for c in s.notable_cases:
            msg = ""
            if c["flag"] == "F":
                msg = "failed solve"
            elif c["flag"] == "T":
                msg = "long runtime"
            elif c["flag"] == "I":
                msg = "high iteration count"
            elif c["flag"] == "TI":
                msg = "high iteration count; long runtime"

            fp.write(f"{c['name']} : {msg}\n")

            for k, v in c["sample_point"].items():
                if k == "_name":
                    # we need this because the parallel task manager does not
                    # handle dicts  and we add _name to identify the sample name
                    # ToDo: fix this when the parallel task manager is extended to
                    # handle dicts
                    continue
                if self.inputs[k]["distribution"] == "normal":
                    mean = self.inputs[k]["mean"]
                    std = self.inputs[k]["std"]
                    alpha = (v - mean) / std
                    fp.write(
                        f"  {k:20s}: {v:10g} ({alpha:5.2f} standard "
                        "deviations above/below the mean)\n"
                    )
