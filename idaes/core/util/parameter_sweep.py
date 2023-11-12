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

from pyomo.core import Param, Var
from pyomo.environ import check_optimal_termination, SolverFactory
from pyomo.common.config import ConfigDict, ConfigValue, document_kwargs_from_configdict

import idaes.logger as idaeslog
from idaes.core.surrogate.pysmo.sampling import SamplingMethods, UniformSampling
from idaes.core.util.exceptions import ConfigurationError

__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


class ParameterSweepSpecification(object):
    """Defines a set of input variables/parameters and values to be used in
    a parameter sweep study.

    Users can define inputs to be sampled, the number of samples and the method to
    use to generate the samples (chosen from Pysmo's sampling methods) and then
    call the generate_samples method to generate the desired set of samples.
    """

    # TODO: Consider supporting sampling from data sets in the future
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
        """Defines the method to use to generate samples. This must be a Pysmo
        SamplingMethod

        Args:
            sampling_method: SamplingMethod class to use.

        Returns:
            None
        """
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
        """Sets the number of samples to be generated.

        The format depends on the sampling method chosen. For UniformSampling, this must
        be a list of ints of length equal to the number of inputs, otherwise it must be
        an int.

        Args:
            sample_size: int or list-of-ints

        Returns:
            None
        """
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
        """Write specification to a dict

        Returns:
            Input specification as a dict
        """
        outdict = {}

        outdict["inputs"] = self._inputs

        # Pysmo sampling methods are not serializable, so get name instead
        outdict["sampling_method"] = self._sampling_method.__name__
        outdict["sample_size"] = self._sample_size
        outdict["samples"] = self._samples.to_dict(orient="tight")

        return outdict

    def from_dict(self, input_dict: dict):
        """
        Load specification from dict

        Args:
            input_dict: dict to load into specification

        Returns:
            None
        """
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

    def to_json_file(self, filename: str):
        """
        Serialize specification to file in json format.

        Args:
            filename: name of file to write to as string

        Returns:
            None
        """
        with open(filename, "w") as fd:
            json.dump(self.to_dict(), fd, indent=3)

    def from_json_file(self, filename: str):
        """
        Load specification from json file.

        Args:
            filename: name of file to load as string

        Returns:
            None
        """
        with open(filename, "r") as f:
            self.from_dict(json.load(f))
        f.close()


def is_psweepspec(val):
    if isinstance(val, ParameterSweepSpecification):
        return val
    _log.error(
        f"Input configuration {val} must be an instance of ParameterSweepSpecification."
    )
    raise ValueError(
        "Input configuration must be an instance of ParameterSweepSpecification."
    )


CONFIG = ConfigDict()
CONFIG.declare(
    "build_model",
    ConfigValue(doc="Callback method to construct initialized model for execution."),
)
CONFIG.declare(
    "run_model",
    ConfigValue(
        doc="Callback method to use when running model. If None, model is run with solver.solve()"
    ),
)
CONFIG.declare(
    "collect_results",
    ConfigValue(
        doc="Callback method to use to collect results from model after execution."
    ),
)
CONFIG.declare(
    "failure_recourse",
    ConfigValue(
        doc="Callback method to use in case of exception when solving a sample."
    ),
)
CONFIG.declare(
    "input_specification",
    ConfigValue(
        domain=is_psweepspec,
        doc="ParameterSweepSpecification object defining inputs to be sampled",
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
class ParameterSweepBase:
    """
    IDAES base class for defining parameter sweep runners.

    Thus base class defines a basic API for setting up parameter sweep studies
    and allows the end-user to provide instructions for how to set up and run
    the model to be studied using a set of callbacks.
    """

    def __init__(self, **kwargs):
        self.config = CONFIG(kwargs)
        self._results = OrderedDict()

    @property
    def results(self):
        return self._results

    def execute_parameter_sweep(self):
        """
        Placeholder method for parameter sweep runners.

        Developers should overload this with a method to execute the parameter sweep
        using their preferred workflow manager.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "Derived classes should overload this method with a " "workflow manager."
        )

    def execute_single_sample(self, sample_id: int):
        """
        Executes a single run of the model using input values from sample_id.

        Args:
            sample_id: int indicating row in specification.samples to load into model
                for run.

        Returns:
            results: results generated by collect_results callback
            solved: bool indicating whether solver reported optimal convergence
        """
        model = self.get_initialized_model()

        # Load sample values
        self.set_input_values(model, sample_id)

        # Solve model and capture output
        solver = self._get_solver()

        # Try/except to catch any critical failures that occur
        try:
            status, run_stats = self.run_model(model, solver)

            solved = check_optimal_termination(status)
            if not solved:
                _log.error(f"Sample: {sample_id} failed to converge.")

            # Compile Results
            results = self.collect_results(model, status, run_stats)
        except:
            # Catch any Exception for recourse
            results, solved = self.execute_recourse(model)

        return results, solved

    def get_initialized_model(self):
        """
        Get instance of model to be run by calling build_model callback.

        Returns:
            Pyomo model to be run
        """
        if self.config.build_model is None:
            raise ConfigurationError(
                "Please specify a method to construct the model of interest."
            )

        model = self.config.build_model()

        # TODO: Verify model is actually a model?

        return model

    def get_input_specification(self):
        """
        Get input specification from config block

        Returns:
            ParameterSweepConfiguration object

        Raises:
            ConfigurationError if no input specification has been defined
        """
        if self.config.input_specification is None:
            raise ConfigurationError(
                "Please specify an input specification to use for sampling."
            )

        return self.config.input_specification

    def get_input_samples(self):
        """
        Get dataframe of samples from input specification.

        Returns:
            pandas dataframe of samples

        Raises:
            ConfigurationError if no input specification has been defined
        """
        spec = self.get_input_specification()

        samples = spec.samples

        if samples is None:
            samples = spec.generate_samples()

        return samples

    def set_input_values(self, model, sample_id: int):
        """
        Set values of input variables/parameters in instance of model using values
        from sample_id.

        Args:
            model: instance of model to be executed
            sample_id: int representing a row in the specification.samples dataframe

        Returns:
            None

        Raises:
            ValueError if an input cannot be found or if it is not a fixed Var or
            mutable Param
        """
        samples = self.get_input_samples()
        inputs = self.config.input_specification.inputs

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

    def run_model(self, model, solver):
        """
        Executes run of model by calling the run_model callback

        Args:
            model: instance of model to be run
            solver: Pyomo solver object to use to solve model

        Returns:
            status: a Pyomo solver status object
            run_stats: additional output collected by run_model callback
        """
        if self.config.run_model is None:
            return solver.solve(model), None

        return self.config.run_model(model, solver)

    def collect_results(self, model, status, run_stats):
        """
        Collects desired results from instance of model by calling collect_results callback.

        Args:
            model: instance of model to collect results from
            status: Pyomo solver status object returned from solving model
            run_stats: additional output from run_model callback to be collected in results

        Returns:
            Output of collect_results callback
        """
        if self.config.collect_results is None:
            raise ConfigurationError(
                "Please provide a method to collect results from sample run."
            )

        return self.config.collect_results(model, status, run_stats)

    def execute_recourse(self, model):
        """
        Call failure_recourse callback. This method is used in case the solver encounters
        a critical error whilst solving a sample.

        Args:
            model: instance of model for performing recourse

        Returns:
            Output of failure_recourse callback
        """
        if self.config.failure_recourse is None:
            # No recourse specified, so solved=False and results=None
            return None, False
        else:
            return self.config.failure_recourse(model)

    def _get_solver(self):
        if self.config.solver is None:
            raise ConfigurationError("Please specify a solver to use.")
        solver = SolverFactory(self.config.solver)
        if self.config.solver_options is not None:
            solver.options = self.config.solver_options
        return solver

    def to_dict(self):
        """
        Serialize specification and current results to dict form

        Returns:
            dict
        """
        # TODO : Need to serialize build_model method somehow?
        outdict = {}
        outdict["specification"] = self.get_input_specification().to_dict()
        outdict["results"] = OrderedDict(self.results)

        return outdict

    def from_dict(self, input_dict: dict):
        """
        Load specification and results from dict.

        Args:
            input_dict: dict to load from

        Returns:
            None
        """
        if self.config.input_specification is not None:
            # Log a warning about overwriting
            _log.debug("Overwriting existing input specification")

        self.config.input_specification = ParameterSweepSpecification()
        self.config.input_specification.from_dict(input_dict["specification"])

        self._results = OrderedDict()
        # Need to iterate to convert string indices to int
        # Converting to json turns the indices to strings
        for k, v in input_dict["results"].items():
            self._results[int(k)] = v

    def to_json_file(self, filename: str):
        """
        Write specification and results to json file.

        Args:
            filename: name of file to write to as string

        Returns:
            None
        """
        with open(filename, "w") as fd:
            json.dump(self.to_dict(), fd, indent=3)

    def from_json_file(self, filename: str):
        """
        Load specification and results from json file.

        Args:
            filename: name of file to load from as string

        Returns:
            None
        """
        with open(filename, "r") as f:
            self.from_dict(json.load(f))
        f.close()


class SequentialSweepRunner(ParameterSweepBase):
    """
    Sequential runner for parameter sweeps.

    This class executes a parameter sweep by running all samples sequentially
    in a for loop.
    """

    def execute_parameter_sweep(self):
        """
        Execute sequential parameter sweep.

        Returns:
            OrderedDict of results indexed by sample ID.
        """
        self._results = OrderedDict()
        samples = self.get_input_samples()

        for s in samples.index:
            sresults, solved = self.execute_single_sample(s)
            self._results[s] = {"solved": solved, "results": sresults}

        return self.results
