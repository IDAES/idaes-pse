#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
IDAES Parameter Sweep API and sequential workflow runner.
"""

import sys
import json

from pandas import DataFrame

from pyomo.core import Param, Var
from pyomo.environ import check_optimal_termination
from pyomo.common.config import ConfigDict, ConfigValue, document_kwargs_from_configdict

import idaes.logger as idaeslog
from idaes.core.surrogate.pysmo.sampling import SamplingMethods, UniformSampling
from idaes.core.util.exceptions import ConfigurationError

__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)

# TODO: Timeouts
# TODO: Re-initialize option/callback


class ParameterSweepSpecification(object):
    """Defines a set of input variables/parameters and values to be used in
    a parameter sweep study.

    Users can define inputs to be sampled, the number of samples and the method to
    use to generate the samples (chosen from Pysmo's sampling methods) and then
    call the generate_samples method to generate the desired set of samples.
    """

    # TODO: Consider supporting sampling from data sets in the future
    def __init__(self):
        self._inputs = {}
        self._sampling_method = None
        self._samples = None
        self._sample_size = None

    def add_sampled_input(
        self,
        pyomo_path: str,
        lower: float,
        upper: float,
        name: str = None,
    ):
        """
        Add an input that should be sampled when forming the set of
        specific points that need to be run for the convergence evaluation

        The input will be sampled assuming a normal distribution
        (with given mean and standard deviation)
        truncated to the values given by lower and upper bounds

        Args:
            pyomo_path A string representation of the path to the variable or parameter to
                be sampled. This string will be executed to retrieve the Pyomo
            component.
            lower : lower bound on the input variable or parameter.
            upper : upper bound on the input variable or parameter
            name : (Optional) display name for the input. If not provided,
                pyomo_pth will be used.

        Returns:
            None
        """
        # ToDo: put some error checking here ... Maybe we should have the model
        # ToDo: already? Can use to check if the pyomo_path is valid? check if
        # ToDo: bounds will be violated?
        spec = {}
        spec["pyomo_path"] = pyomo_path
        spec["lower"] = lower
        spec["upper"] = upper

        if name is None:
            # If no name given, use Pyomo path
            name = pyomo_path

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
            raise ValueError("Please identify at least one input variable to sample.")

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

        # pylint: disable=not-callable
        space_init = self.sampling_method(
            bounds_list,
            sample_size,
            sampling_type="creation",
        )

        self._samples = DataFrame(
            space_init.sample_points(), columns=[_ for _ in self.inputs]
        )

        return self.samples

    @property
    def inputs(self):
        """
        Returns an dict containing the declared inputs.
        """
        return self._inputs

    @property
    def sampling_method(self):
        """
        Returns the declared sampling method.
        """
        return self._sampling_method

    @property
    def sample_size(self):
        """
        Returns the declared sample size.
        """
        return self._sample_size

    @property
    def samples(self):
        """
        Returns the generated set of samples (pandas DataFrame).
        """
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
        self._inputs = {}
        for k, v in input_dict["inputs"].items():
            self._inputs[k] = v

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
    """
    Config validator for ParameterSweepSpecifications

    Args:
        val: value to validate

    Returns:
        ParameterSweepSpecification

    Raises:
        ValueError if val is not a ParameterSweepSpecification
    """
    if isinstance(val, ParameterSweepSpecification):
        return val
    _log.error(
        f"Input configuration {val} must be an instance of ParameterSweepSpecification."
    )
    raise ValueError(
        "Input configuration must be an instance of ParameterSweepSpecification."
    )


def _is_solver(val):
    """
    Config validator for Solver object.

    Use duck-typing and just look for a solve method

    Args:
        val: value to validate

    Returns:
        val if it has a solve method

    Raises:
        ValueError if val does not have a solve method
    """
    try:
        if callable(val.solve):
            return val
    except AttributeError:
        raise ValueError("solver must be an instance of a Pyomo Solver object.")


CONFIG = ConfigDict()
CONFIG.declare(
    "rebuild_model",
    ConfigValue(
        default=True,
        domain=bool,
        doc="Whether to rebuild model for each sample, or reuse the same instance for "
        "all runs (default=True, rebuild for all runs).",
    ),
)
CONFIG.declare(
    "build_model",
    ConfigValue(doc="Callback method to construct initialized model for execution."),
)
CONFIG.declare(
    "build_model_arguments",
    ConfigValue(
        domain=dict,
        doc="Arguments to pass to build_model callback.",
    ),
)
CONFIG.declare(
    "run_model",
    ConfigValue(
        doc="Callback method to use when running model. If None, model is run with solver.solve()"
    ),
)
CONFIG.declare(
    "run_model_arguments",
    ConfigValue(
        domain=dict,
        doc="Arguments to pass to run_model callback.",
    ),
)
CONFIG.declare(
    "build_outputs",
    ConfigValue(
        doc="Callback method to use to collect results from model after execution."
    ),
)
CONFIG.declare(
    "build_outputs_arguments",
    ConfigValue(
        domain=dict,
        doc="Arguments to pass to build_outputs callback.",
    ),
)
CONFIG.declare(
    "handle_solver_error",
    ConfigValue(
        doc="Callback method to use in case of exception when solving a sample."
    ),
)
CONFIG.declare(
    "handle_solver_error_arguments",
    ConfigValue(
        domain=dict,
        doc="Arguments to pass to handle_solver_error callback.",
    ),
)
CONFIG.declare(
    "halt_on_error",
    ConfigValue(
        default=False,
        domain=bool,
        doc="Whether to halt execution of parameter sweep on encountering a solver error (default=False).",
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
        domain=_is_solver,
        doc="Pyomo solver object to use when solving model",
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
        self._results = {}
        self._model = None  # used to store model instance if rebuild_model is False

    @property
    def results(self):
        """
        Returns dict containing the results from the parameter sweep.
        """
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
            "Derived classes should overload this method with a workflow manager."
        )

    def execute_single_sample(self, sample_id: int):
        """
        Executes a single run of the model using input values from sample_id.

        Args:
            sample_id: int indicating row in specification.samples to load into model
                for run.

        Returns:
            results: results generated by build_outputs callback
            success: bool indicating whether execution was successful
            error: str if error occurs during solve else None
        """
        model = self.get_initialized_model()

        # Load sample values
        self.set_input_values(model, sample_id)

        # Try/except to catch any critical failures that occur
        error = None
        try:
            success, run_stats = self.run_model(model, self.config.solver)

            if not success:
                _log.warning(f"Sample {sample_id} did not report success.")
        except Exception as e:  # pylint: disable=broad-except
            if self.config.halt_on_error:
                raise

            success = False
            error = str(e)  # Cast to string for storage

        # Compile Results
        if error is None:
            results = self.build_outputs(model, run_stats)
        else:
            # Catch any Exception for recourse
            results = self.handle_error(model)

        _log.info(f"Sample {sample_id} finished.")

        return results, success, error

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

        if not self.config.rebuild_model:
            # If reusing model, see if instance has been constructed yet
            if self._model is not None:
                # If yes, return and done
                return self._model

        # Otherwise, build instance of model
        args = self.config.build_model_arguments
        if args is None:
            args = {}
        model = self.config.build_model(**args)
        # TODO: Verify model is actually a model?

        if not self.config.rebuild_model:
            # If reusing model, store instance for reuse
            self._model = model

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

    def get_sample_values(self, sample_id: int):
        """
        Get inputs for a specific sample run indicated by sample_id.

        Args:
            sample_id: int representing a row in the specification.samples dataframe

        Returns:
            Pandas Series of input values for chosen sample
        """
        samples = self.get_input_samples()

        return samples.iloc[sample_id]

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
                try:
                    if not comp.is_fixed():
                        raise ValueError(
                            f"Convergence testing found an input of type Var that "
                            f"was not fixed ({comp.name}). Please make sure all "
                            f"sampled inputs are either mutable params or fixed vars."
                        )
                    comp.set_value(float(v))
                except AttributeError:
                    # Component might be indexed, try iterating and setting value
                    for i, c in comp.items():
                        if not c.is_fixed():
                            raise ValueError(
                                f"Convergence testing found an input of type IndexedVar that "
                                f"was not fixed ({comp.name}, index {i}). Please make sure all "
                                f"sampled inputs are either mutable params or fixed vars."
                            )
                        c.set_value(float(v))
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
            success: bool indicating whether execution was successful or not
            run_stats: output collected by run_model callback (default is Pyomo SolverResults object)
        """
        if self.config.run_model is None:
            res = solver.solve(model)

            success = check_optimal_termination(res)
            return success, res

        args = self.config.run_model_arguments
        if args is None:
            args = {}

        return self.config.run_model(model, solver, **args)

    def build_outputs(self, model, run_stats):
        """
        Collects desired results from instance of model by calling build_outputs callback.

        Args:
            model: instance of model to collect results from
            run_stats: output from run_model callback to be collected in results
                (default is Pyomo results object)

        Returns:
            Output of build_outputs callback. If no build_outputs callback is provided,
            returns run_stats instead.
        """
        if self.config.build_outputs is None:
            return run_stats
        else:
            args = self.config.build_outputs_arguments
            if args is None:
                args = {}

            return self.config.build_outputs(model, run_stats, **args)

    def handle_error(self, model):
        """
        Call handle_solver_error callback. This method is used in case the solver encounters
        a critical error whilst solving a sample.

        Args:
            model: instance of model for performing recourse

        Returns:
            Output of handle_solver_error callback
        """
        if self.config.handle_solver_error is None:
            # No recourse specified, so results=None
            return None
        else:
            args = self.config.handle_solver_error_arguments
            if args is None:
                args = {}
            return self.config.handle_solver_error(model, **args)

    def to_dict(self):
        """
        Serialize specification and current results to dict form

        Returns:
            dict
        """
        # TODO : Need to serialize build_model method somehow?
        outdict = {}
        outdict["specification"] = self.get_input_specification().to_dict()
        outdict["results"] = self.results

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

        self._results = {}
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

    @staticmethod
    def progress_bar(fraction, msg: str, length: int = 20):
        """
        Prints a progress bar to stdout

        Args:
            fraction: fraction of total samples executed
            msg: string to append to progress bar
            length: length of progress bar (default=20)

        Returns:
            None
        """
        n_complete = int(length * fraction)
        n_remaining = length - n_complete
        characters = f"{'*' * n_complete}{'-' * n_remaining}"
        sys.stdout.write(f"{fraction * 100:.1f}% {characters} {msg}\n")
        sys.stdout.flush()


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
            dict of results indexed by sample ID.
        """
        self._results = {}
        samples = self.get_input_samples()

        count = 1
        for s in samples.index:
            sresults, success, error = self.execute_single_sample(s)
            self._results[s] = {"success": success, "results": sresults, "error": error}

            self.progress_bar(float(count) / float(len(samples)), "Complete")
            count += 1

        return self.results
