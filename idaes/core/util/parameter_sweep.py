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

import json
from collections import OrderedDict

from pandas import DataFrame
from pandas.testing import assert_frame_equal

from pyomo.core import Param, Var
from pyomo.environ import check_optimal_termination, SolverFactory
from pyomo.common.config import ConfigDict, ConfigValue, document_kwargs_from_configdict

import idaes.logger as idaeslog
from idaes.core.surrogate.pysmo.sampling import SamplingMethods, UniformSampling
from idaes.core.util.exceptions import ConfigurationError

# Set up logger
_log = idaeslog.getLogger(__name__)


class ParameterSweepSpecification(object):
    def __init__(self):
        self._inputs = OrderedDict()
        self._sampling_method = None
        self._samples = None
        self._sample_size = None

    def add_sampled_input(
        self,
        name: str,
        pyomo_path: str,
        lower: float,
        upper: float,
    ):
        """
        Add an input that should be sampled when forming the set of
        specific points that need to be run for the convergence evaluation

        The input will be sampled assuming a normal distribution
        (with given mean and standard deviation)
        truncated to the values given by lower and upper bounds

        Args:
            name : the name of the input.
            pyomo_path A string representation of the path to the variable or parameter to
                be sampled. This string will be executed to retrieve the Pyomo
            component.
            lower : lower bound on the input variable or parameter.
            upper : upper bound on the input variable or parameter

        Returns:
            None
        """
        # ToDo: put some error checking here ... Maybe we should have the model
        # ToDo: already? Can use to check if the pyomo_path is valid? check if
        # ToDo: bounds will be violated?
        spec = OrderedDict()
        spec["pyomo_path"] = pyomo_path
        spec["lower"] = lower
        spec["upper"] = upper
        self._inputs[name] = spec

    def set_sampling_method(self, sampling_method):
        try:
            if issubclass(sampling_method, SamplingMethods):
                self._sampling_method = sampling_method
            else:
                raise TypeError
        except TypeError:
            raise TypeError(
                f"Sampling method must be an instance of a Pysmo SamplingMethod "
                f"(received {sampling_method})"
            )

    def set_sample_size(self, sample_size):
        self._sample_size = sample_size

    def _generate_pysmo_data_input(self):
        lb_list = []
        ub_list = []

        for var_def in self.inputs.values():
            lb_list.append(var_def["lower"])
            ub_list.append(var_def["upper"])

        return [lb_list, ub_list]

    def generate_samples(self):
        """
        Generate set of samples of specified size from inputs.

        Args:
            sample_size - [int or list of ints] size of sample to be generated.
                For UniformSampling this must be a list of ints indicating number of
                graduations for each input. For all other sampling methods this must be
                an int indicating the number of samples to generate.

        Returns:
            Pandas DataFrame of samples
        """
        if self.sampling_method is None:
            raise ValueError(
                "Please choose a sampling method to use for sample generation. "
                "Sampling method must be chosen from those available in Pysmo."
            )

        if len(self._inputs) == 0:
            raise ValueError("Please identify at least on input variable to sample.")

        if self._sample_size is None:
            raise ValueError("Please set a sample size.")
        elif self.sampling_method is UniformSampling:
            if not isinstance(self._sample_size, list):
                raise TypeError(
                    "For UniformSampling, sample_size must be list of integers."
                )
            sample_size = self._sample_size
        else:
            try:
                sample_size = int(self._sample_size)
            except ValueError:
                raise ValueError("sample_size must be an integer.")
            if sample_size <= 0:
                raise ValueError("sample_size must be an integer greater than 1.")

        bounds_list = self._generate_pysmo_data_input()

        space_init = self.sampling_method(
            bounds_list,
            sample_size,
            sampling_type="creation",
        )

        self._samples = DataFrame(
            space_init.sample_points(), columns=[_ for _ in self.inputs.keys()]
        )

        return self.samples

    @property
    def inputs(self):
        return self._inputs

    @property
    def sampling_method(self):
        return self._sampling_method

    @property
    def sample_size(self):
        return self._sample_size

    @property
    def samples(self):
        return self._samples

    def to_dict(self):
        outdict = {}

        outdict["inputs"] = self._inputs

        # Pysmo sampling methods are not serializable, so get name instead
        outdict["sampling_method"] = self._sampling_method.__name__
        outdict["sample_size"] = self._sample_size
        outdict["samples"] = self._samples.to_dict(orient="tight")

        return outdict

    def from_dict(self, input_dict):
        self._inputs = OrderedDict()
        for k, v in input_dict["inputs"].items():
            self._inputs[k] = OrderedDict(v)

        # Hack to get Pysom sampling method from string name
        mod = __import__(
            "idaes.core.surrogate.pysmo.sampling",
            fromlist=[input_dict["sampling_method"]],
        )
        self._sampling_method = getattr(mod, input_dict["sampling_method"])

        self._sample_size = input_dict["sample_size"]
        self._samples = DataFrame().from_dict(
            input_dict["samples"],
            orient="tight",
        )

    def to_json_file(self, filename):
        with open(filename, "w") as fd:
            json.dump(self.to_dict(), fd, indent=3)

    def from_json_file(self, filename):
        with open(filename, "r") as f:
            self.from_dict(json.load(f))
        f.close()


def is_psweepspec(val):
    if isinstance(val, ConvergenceEvaluationSpecification):
        return val
    _log.error(
        f"Input configuration {val} must be an instance of ConvergenceEvaluationSpecification."
    )
    raise ValueError(
        "Input configuration must be an instance of ConvergenceEvaluationSpecification."
    )


CONFIG = ConfigDict()
CONFIG.declare(
    "build_model",
    ConfigValue(doc="Callback method to construct initialized model for execution."),
)
CONFIG.declare(
    "execute_model",
    ConfigValue(doc="Callback method to use when running model."),
)
CONFIG.declare(
    "collect_results",
    ConfigValue(
        doc="Callback method to use to collect results from model after execution."
    ),
)
CONFIG.declare(
    "recourse",
    ConfigValue(
        doc="Callback method to use if an error occurs during normal execution."
    ),
)
CONFIG.declare(
    "workflow_runner",
    ConfigValue(
        domain=None,
        doc="Callback setting up workflow manager for parameter sweep",
    ),
)
CONFIG.declare(
    "input_specification",
    ConfigValue(
        domain=is_psweepspec,
        doc="ConvergenceEvaluationSpecification object defining inputs to be sampled",
    ),
)
CONFIG.declare(
    "solver",
    ConfigValue(
        default=None,
        domain=str,
        doc="Name of solver to use (recognised by Pyomo SolverFactory)",
    ),
)
CONFIG.declare(
    "solver_options",
    ConfigValue(
        default=None,
        domain=dict,
        doc="Dict of options to pass to solver",
    ),
)


@document_kwargs_from_configdict(CONFIG)
class ParameterSweep:
    def __init__(self, **kwargs):
        self.config = CONFIG(kwargs)
        self._input_spec = None
        self._results = None

    @property
    def results(self):
        return self._results

    def get_specification(self):
        """
        Returns input specification to use for parameter sweep from the
        configuration block.

        For legacy purposes, user can override this method to return an
        instance of the ConvergenceEvaluationSpecification for this particular
        model and test set.

        The basic flow for this method is:
           - Create a ConvergenceEvaluationSpecification.
           - Call add_sampled_input for every input that should be varied.
           - Call generate_samples().
           - return the ConvergenceEvaluationSpecification.

        Returns:
           ConvergenceEvaluationSpecification
        """
        if self._input_spec is None:
            if self.config.input_specification is not None:
                self._input_spec = self.config.input_specification
            else:
                raise ConfigurationError(
                    "Please specify an input specification to use for sampling."
                )
        return self._input_spec

    def get_initialized_model(self):
        """
        Returns the initialized Pyomo model to be executed. By default, this calls
        the build_model callback method defined in the configuration block. The
        convergence evaluation methods will change the values of parameters or
        variables according to the sampling specifications.

        For legacy purposes, user can also overload this method to return an
        initialized model that is ready to solve.

        Returns:
           Pyomo model - return a Pyomo model object that is initialized and
                ready to solve. This is the model object that will be
                used in the evaluation.
        """
        if self.config.build_model is not None:
            return self.config.build_model()
        raise ConfigurationError(
            "Please specify a method to construct the model of interest."
        )

    def get_input_samples(self):
        if self._input_spec is None:
            self.get_specification()

        samples = self._input_spec.samples

        if samples is None:
            samples = self._input_spec.generate_samples()

        return samples

    def _set_input_values(self, model, sample_id):
        samples = self.get_input_samples()
        inputs = self.get_specification().inputs

        for k, i in inputs.items():
            v = samples[k][sample_id]

            # k stores the "name" of the input
            # need to get the pyomo path of the input from the inputs structure
            pyomo_path = i["pyomo_path"]

            comp = model.find_component(pyomo_path)
            try:
                ctype = comp.ctype
            except AttributeError:
                ctype = None

            # TODO: Validate bounds

            if ctype is Param:
                if comp.is_constant():
                    raise ValueError(
                        f"Convergence testing found an input of type Param that "
                        f"was not mutable ({comp.name}). Please make sure all "
                        f"sampled inputs are either mutable params or fixed vars."
                    )
                comp.set_value(v)
            elif ctype is Var:
                if not comp.is_fixed():
                    raise ValueError(
                        f"Convergence testing found an input of type Var that "
                        f"was not fixed ({comp.name}). Please make sure all "
                        f"sampled inputs are either mutable params or fixed vars."
                    )
                comp.set_value(float(v))
            else:
                raise ValueError(
                    f"Failed to find a valid input component (must be "
                    f"a fixed Var or a mutable Param). Instead, "
                    f"pyomo_path: {pyomo_path} returned: {comp}."
                )

    def _get_solver(self):
        solver = SolverFactory(self.config.solver)
        if self.config.solver_options is not None:
            solver.options = self.config.solver_options
        return solver

    def execute_single_sample(self, sample_id):
        # TODO: Add recourse for critical failure
        # Try/except to catch any failures that occur
        try:
            model = self.get_initialized_model()

            # Load sample values
            self._set_input_values(model, sample_id)

            # Solve model and capture output
            solver = self._get_solver()
            status = self.config.execute_model(model, solver)

            solved = check_optimal_termination(status)
            if not solved:
                _log.error(f"Sample: {sample_id} failed to converge.")

            # Compile Results
            results = self.config.collect_results(model)
        except:
            # Catch any exception here
            solved = False
            results = self.config.recourse(model)

        return results, solved

    def to_dict(self):
        # TODO : Need to serialize build_model method somehow?
        outdict = {}
        outdict["specification"] = self.get_specification().to_dict()
        outdict["results"] = self._results

        return outdict

    def from_dict(self, input_dict):
        self._input_spec = ConvergenceEvaluationSpecification()
        self._input_spec.from_dict(input_dict["specification"])

        self._results = OrderedDict()
        # Need to iterate to convert string indices to int
        # Converting to json turns the indices to strings
        for k, v in input_dict["results"].items():
            self._results[int(k)] = v

    def to_json_file(self, filename):
        with open(filename, "w") as fd:
            json.dump(self.to_dict(), fd, indent=3)

    def from_json_file(self, filename):
        with open(filename, "r") as f:
            self.from_dict(json.load(f))
        f.close()

    def _verify_samples_from_dict(self, compare_dict):
        comp_samples = DataFrame().from_dict(
            compare_dict["specification"]["samples"],
            orient="tight",
        )
        try:
            assert_frame_equal(self.get_input_samples(), comp_samples)
        except AssertionError:
            raise ValueError(
                "Samples in comparison evaluation do not match current evaluation"
            )
