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
from enum import Enum
import subprocess
from io import StringIO
import sys
import os
import numpy as np
import pandas as pd
import json

from pyomo.environ import Constraint, value, sin, cos, log, exp, Set, Reals
from pyomo.common.config import ConfigValue, In, Path, ListOf, Bool
from pyomo.common.tee import TeeStream
from pyomo.common.fileutils import Executable
from pyomo.common.tempfiles import TempfileManager

from idaes.core.surrogate.base.surrogate_base import SurrogateTrainer, SurrogateBase
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)

# TODO: Adaptive sampling

alamo = Executable("alamo")

# Define mapping of Pyomo function names for expression evaluation
GLOBAL_FUNCS = {"sin": sin, "cos": cos, "log": log, "exp": exp}


# The values associated with these must match those expected in the .alm file
class Modelers(Enum):
    BIC = 1
    MallowsCp = 2
    AICc = 3
    HQC = 4
    MSE = 5
    SSEP = 6
    RIC = 7
    MADp = 8


class Screener(Enum):
    none = 0
    lasso = 1
    SIS = 2


supported_options = [
    "xfactor",
    "xscaling",
    "scalez",
    "monomialpower",
    "multi2power",
    "multi3power",
    "ratiopower",
    "expfcns",
    "linfcns",
    "logfcns",
    "sinfcns",
    "cosfcns",
    "constant",
    "grbfcns",
    "rbfparam",
    "modeler",
    "builder",
    "backstepper",
    "convpen",
    "screener",
    "ncvf",
    "sismult",
    "maxtime",
    "maxiter",
    "datalimitterms",
    "maxterms",
    "minterms",
    "numlimitbasis",
    "exclude",
    "ignore",
    "xisint",
    "zisint",
    "tolrelmetric",
    "tolabsmetric",
    "tolmeanerror",
    "tolsse",
    "mipoptca",
    "mipoptcr",
    "linearerror",
    "GAMS",
    "GAMSSOLVER",
    "solvemip",
    "print_to_screen",
]

"""
Excluded Options:
NSAMPLE, NVALSAMPLE, INITIALIZER, PRINT TO FILE, FUNFORM, NTRANS

Excluded Simulator options:
MAXSIM, MINPOINTS, MAXPOINTS, SAMPLER, SIMULATOR, PRESET, TOLMAXERROR, SIMIN,
SIMOUT
"""


# Headers from ALAMO trace file that should be common for all outputs
common_trace = [
    "filename",
    "NINPUTS",
    "NOUTPUTS",
    "INITIALPOINTS",
    "SET",
    "INITIALIZER",
    "SAMPLER",
    "MODELER",
    "BUILDER",
    "GREEDYBUILD",
    "BACKSTEPPER",
    "GREEDYBACK",
    "REGULARIZER",
    "SOLVEMIP",
]


class AlamoTrainer(SurrogateTrainer):
    """
    Standard SurrogateTrainer for ALAMO.

    This defines a set of configuration options for ALAMO along with
    methods to read and write the ALAMO input and output files and to call
    the ALAMO executable.

    Generally, options default to None to indicate that no entry will be
    written in the ALAMO input file. In this case, the default ALAMO settings
    will be used.
    """

    # The following ALAMO options are not (yet) supported
    # Returning model prediction, as we can do that in IDAES
    # Iniital sampling of data, as SurrogateModelTrainer should do that
    # Similarly, adaptive sampling is better handled in Python
    # Single validation set, due to limitations of current API
    # Custom basis functions are not yet implemented

    CONFIG = SurrogateTrainer.CONFIG()

    CONFIG.declare(
        "xfactor",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="List of scaling factors for input variables.",
        ),
    )
    CONFIG.declare(
        "xscaling",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Option to scale input variables.",
            doc="Option to scale input variables. If True and xfactors are not "
            "provided, ALAMO sets XFACTORS equal to the range of each input "
            "variable.",
        ),
    )
    CONFIG.declare(
        "scalez",
        ConfigValue(
            default=None, domain=Bool, description="Option to scale output variables."
        ),
    )

    # Basis function options
    CONFIG.declare(
        "monomialpower",
        ConfigValue(
            default=None,
            domain=ListOf(
                int, In(Reals - {0, 1})
            ),  # allow any float except for 0 and 1
            description="Vector of monomial powers considered in basis "
            "functions - cannot include 0 or 1.",
        ),
    )
    CONFIG.declare(
        "multi2power",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="Vector of powers to be considered for pairwise "
            "combinations in basis functions.",
        ),
    )
    CONFIG.declare(
        "multi3power",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="Vector of three variable combinations of powers to be "
            "considered as basis functions.",
        ),
    )
    CONFIG.declare(
        "ratiopower",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="Vector of ratio combinations of powers to be considered "
            "in the basis functions.",
        ),
    )
    CONFIG.declare(
        "constant",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Include constant basis function if True. Default = True",
        ),
    )
    CONFIG.declare(
        "linfcns",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Include linear basis functions if True. Default = True",
        ),
    )
    CONFIG.declare(
        "expfcns",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Include exponential basis functions if True.",
        ),
    )
    CONFIG.declare(
        "logfcns",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Include logarithmic basis functions if True.",
        ),
    )
    CONFIG.declare(
        "sinfcns",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Include sine basis functions if True.",
        ),
    )
    CONFIG.declare(
        "cosfcns",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Include cosine basis functions if True.",
        ),
    )
    CONFIG.declare(
        "grbfcns",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Include Gaussian radial basis functions if True.",
        ),
    )
    CONFIG.declare(
        "rbfparam",
        ConfigValue(
            default=None,
            domain=float,
            description="Multiplicative constant used in the Gaussian radial basis"
            " functions.",
        ),
    )
    CONFIG.declare(
        "custom_basis_functions",
        ConfigValue(
            default=None,
            domain=ListOf(str),
            description="List of custom basis functions to include in surrogate "
            "fitting.",
            doc="List of custom basis functions to include in surrogate model "
            "fitting. These should be in a form that can be rendered as a string "
            "that meets ALAMO's requirements.",
        ),
    )

    # Other fitting options
    CONFIG.declare(
        "modeler",
        ConfigValue(
            default=None,
            domain=In(Modelers),
            description="Fitness metric to be used for model building. Must be an "
            "instance of Modelers Enum.",
        ),
    )
    CONFIG.declare(
        "builder",
        ConfigValue(
            default=None,
            domain=Bool,
            description="If True, a greedy heuristic builds up a model "
            "by adding one variable at a time.",
            doc="If True, a greedy heuristic builds up a model by adding one "
            "variable at a time. This model is used as a starting point for "
            "solving an integer programming formulation according to the choice "
            "of modeler. If an optimizer is not available, the heuristic model "
            "will be the final model to be returned.",
        ),
    )
    CONFIG.declare(
        "backstepper",
        ConfigValue(
            default=None,
            domain=Bool,
            description="If set to 1, a greedy heuristic builds down a model by "
            "starting from the least squares model and removing one variable at "
            "a time.",
        ),
    )
    CONFIG.declare(
        "convpen",
        ConfigValue(
            default=None,
            domain=float,
            description="Convex penalty term to use if Modeler == SSEP or MADp.",
            doc="When MODELER is set to 6 or 8, a penalty consisting of the sum "
            "of square errors (SSEP) or the maximum absolute error (MADp) and a "
            "term penalizing model size is used for model building. The size of "
            "the model is weighted by convpen. If convpen=0, this metric reduces "
            "to the classical sum of square errors (SSEP) or the maximum absolute "
            "deviation (MADp).",
        ),
    )
    CONFIG.declare(
        "screener",
        ConfigValue(
            default=None,
            domain=In(Screener),
            description="Regularization method used to reduce the number of "
            "potential basis functions before optimization. Must be instance of "
            "Screener Enum.",
        ),
    )
    CONFIG.declare(
        "ncvf",
        ConfigValue(
            default=None,
            domain=int,
            description="Number of folds to be used for cross validation by the "
            "lasso screener. ALAMO will use a two-fold validation if fewer than "
            "10 data points are available. NCVF must be a nonnegative integer.",
        ),
    )
    CONFIG.declare(
        "sismult",
        ConfigValue(
            default=None,
            domain=int,
            description="This parameter must be non-negative and is used to "
            "determine the number of basis functions retained by the SIS "
            "screener. The number of basis functions retained equals the floor "
            "of SSISmult n ln(n), where n is the number of measurements "
            "available at the current ALAMO iteration.",
        ),
    )

    CONFIG.declare(
        "maxiter",
        ConfigValue(
            default=None,
            domain=int,
            description="Maximum number of ALAMO iterations. 1 = no adaptive "
            "sampling, 0 = no limit.",
            doc="Maximum number of ALAMO iterations. Each iteration begins with "
            "a model-building step. An adaptive sampling step follows if maxiter "
            "does not equal 1. If maxiter is set to a number less than or equal "
            "to 0, ALAMO will enforce no limit on the number of iterations.",
        ),
    )
    CONFIG.declare(
        "maxtime",
        ConfigValue(
            default=1000,
            domain=float,
            description="Maximum total execution time allowed in seconds. "
            "Default = 1000.",
            doc="Maximum total execution time allowed in seconds. This time "
            "includes all steps of the algorithm, including time to read problem, "
            "preprocess data, solve optimization subproblems, and print results.",
        ),
    )
    CONFIG.declare(
        "datalimitterms",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Limit model terms to number of measurements.",
            doc="If True, ALAMO will limit the number of terms in the model to be "
            "no more than the number of data measurements; otherwise, no limit "
            "based on the number of data measurements will be placed. The user "
            "may provide an additional limit on the number of terms in the model "
            "through the maxterms and minterms options.",
        ),
    )
    CONFIG.declare(
        "maxterms",
        ConfigValue(
            default=None,
            domain=ListOf(int),
            description="List of maximum number of model terms to per output.",
            doc="Row vector of maximum terms allowed in the modeling of output "
            "variables. One per output variable, space separated. A −1 signals "
            "that no limit is imposed.",
        ),
    )
    CONFIG.declare(
        "minterms",
        ConfigValue(
            default=None,
            domain=ListOf(int),
            description="List of minimum number of model terms to per output.",
            doc="Row vector of minimum terms required in the modeling of output "
            "variables. One per output variable, space separated. A 0 signals "
            "that no limit is imposed.",
        ),
    )
    CONFIG.declare(
        "numlimitbasis",
        ConfigValue(
            default=True,  # default this to true to avoid numerical issues
            domain=Bool,
            description="Eliminate infeasible basis functions. Default = True",
            doc="If True, ALAMO will eliminate basis functions that are not "
            "numerically acceptable (e.g., log(x) will be eliminated if x may be "
            "negative); otherwise, no limit based on the number of data "
            "measurements will be placed. The user may provide additional limits "
            "on the the type and number of selected basis functions through the "
            "options exclude and groupcon.",
        ),
    )
    CONFIG.declare(
        "exclude",
        ConfigValue(
            default=None,
            domain=ListOf(int),
            description="List of inputs to exclude during building,",
            doc="Row vector of 0/1 flags that specify which input variables, if "
            "any, ALAMO should exclude during the model building process. All "
            "input variables must be present in the data but ALAMO will not "
            "include basis functions that involve input variables for which "
            "exclude equals 1. This feature does not apply to custom basis "
            "functions or RBFs.",
        ),
    )
    CONFIG.declare(
        "ignore",
        ConfigValue(
            default=None,
            domain=ListOf(int),
            description="List of outputs to ignore during building.",
            doc="Row vector of 0/1 flags that specify which output variables, "
            "if any, ALAMO should ignore. All output variables must be present in "
            "the data but ALAMO does not model output variables for which ignore "
            "equals 1.",
        ),
    )
    CONFIG.declare(
        "xisint",
        ConfigValue(
            default=None,
            domain=ListOf(int),
            description="List of inputs that should be treated as integers.",
            doc="Row vector of 0/1 flags that specify which input variables, if "
            "any, ALAMO should treat as integers. For integer inputs, ALAMO’s "
            "sampling will be restricted to integer values.",
        ),
    )
    CONFIG.declare(
        "zisint",
        ConfigValue(
            default=None,
            domain=ListOf(int),
            description="List of outputs that should be treated as integers.",
            doc="Row vector of 0/1 flags that specify which output variables, if "
            "any, ALAMO should treat as integers. For integer variables, ALAMO’s "
            "model will include the rounding of a function to the nearest integer "
            "(equivalent to the nint function in Fortran.)",
        ),
    )

    CONFIG.declare(
        "tolrelmetric",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="Relative tolerance for outputs.",
            doc="Relative convergence tolerance for the chosen fitness metric for "
            "the modeling of output variables. One per output variable, space "
            "separated. Incremental model building will stop if two consecutive "
            "iterations do not improve the chosen metric by at least this amount.",
        ),
    )
    CONFIG.declare(
        "tolabsmetric",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="Absolute tolerance for outputs.",
            doc="Absolute convergence tolerance for the chosen fitness metric for "
            "the modeling of output variables. One per output variable, space "
            "separated. Incremental model building will stop if two consecutive "
            "iterations do not improve the chosen metric by at least this amount.",
        ),
    )
    CONFIG.declare(
        "tolmeanerror",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="Convergence tolerance for mean errors in outputs.",
            doc="Row vector of convergence tolerances for mean errors in the "
            "modeling of output variables. One per output variable, space "
            "separated. Incremental model building will stop if tolmeanerror, "
            "tolrelmetric, or tolabsmetric is satisfied.",
        ),
    )
    CONFIG.declare(
        "tolsse",
        ConfigValue(
            default=None,
            domain=float,
            description="Absolute tolerance on SSE",
            doc="Absolute tolerance on sum of square errors (SSE). ALAMO will "
            "terminate if it finds a solution whose SSE is within tolsse from "
            "the SSE of the full least squares problem.",
        ),
    )

    CONFIG.declare(
        "mipoptca",
        ConfigValue(
            default=None, domain=float, description="Absolute tolerance for MIP."
        ),
    )
    CONFIG.declare(
        "mipoptcr",
        ConfigValue(
            default=None, domain=float, description="Relative tolerance for MIP."
        ),
    )
    CONFIG.declare(
        "linearerror",
        ConfigValue(
            default=None,
            domain=Bool,
            description="If True, a linear objective is used when solving "
            "mixed-integer optimization problems; otherwise, a squared error will "
            "be employed.",
        ),
    )
    CONFIG.declare(
        "GAMS",
        ConfigValue(
            default=None,
            domain=str,
            description="Complete path of GAMS executable (or name if GAMS is in "
            "the user path).",
        ),
    )
    CONFIG.declare(
        "GAMSSOLVER",
        ConfigValue(
            default=None,
            domain=str,
            description="Name of preferred GAMS solver for solving ALAMO’s "
            "mixed-integer quadratic subproblems. Special facilities have been "
            "implemented in ALAMO and BARON that make BARON the preferred "
            "selection for this option. However, any mixed-integer quadratic "
            "programming solver available under GAMS can be used.",
        ),
    )
    CONFIG.declare(
        "solvemip",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Whether to use an optimizer to solve MIP.",
        ),
    )
    CONFIG.declare(
        "print_to_screen",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Send ALAMO output to stdout. Output is returned by the "
            "call_alamo method.",
        ),
    )

    # I/O file options
    CONFIG.declare(
        "alamo_path",
        ConfigValue(
            default=None, domain=Path, doc="Path to ALAMO executable (if not in path)."
        ),
    )
    CONFIG.declare(
        "filename",
        ConfigValue(
            default=None,
            domain=str,
            description="File name to use for ALAMO files - must be full path of a"
            " .alm file. Other files will be defined from this pattern. If this "
            "option is not None, then working files will not be deleted.",
        ),
    )
    CONFIG.declare(
        "working_directory",
        ConfigValue(
            default=None,
            domain=str,
            description="Full path to working directory for ALAMO to use. "
            "If this option is not None, then working files will not be deleted.",
        ),
    )
    CONFIG.declare(
        "overwrite_files",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Flag indicating whether existing files can be " "overwritten.",
        ),
    )

    # TODO: We need to do some processing of the labels since ALAMO is
    # restrictive about the labels
    # TODO: We need to think more carefully about "input_bounds".
    # Alamo uses bounds during training, but we also want to consider
    # "valid" bounds for the surrogate
    def __init__(self, **settings):
        super().__init__(**settings)

        self._temp_context = None
        self._almfile = None
        self._trcfile = None

        if self.config.alamo_path is not None:
            alamo.executable = self.config.alamo_path

        self._results = None

    def train_surrogate(self):
        """
        General workflow method for training an ALAMO surrogate.

        Takes the existing data set and executes the ALAMO workflow to create
        an AlamoSurrogate based on the current configuration arguments.

        Args:
            None

        Returns:
            tuple : (success, AlamoSurrogate, message) where success indicates
            whether ALAMO was usccessfully executed, an instance of an
            AlamoSurrogate representing the trained surrogate, and message is
            the final status line from the ALAMO output log.
        """
        # Get paths for temp files
        self._get_files()

        return_code = None
        alamo_log = None
        alamo_object = None

        # Check for any issues in arguments that might cause problems
        self._verify_inputs()

        try:
            # Write .alm file
            self._write_alm_file()

            # Call ALAMO executable
            return_code, alamo_log = self._call_alamo()

            # Read back results
            trace_dict = self._read_trace_file(self._trcfile)

            # Populate results and SurrogateModel object
            self._populate_results(trace_dict)
            alamo_object = self._build_surrogate_object()

        finally:
            # Clean up temporary files if required
            self._remove_temp_files()

        success = False
        if return_code == 0 and alamo_object is not None:
            # Non-zero return code implies an error
            # specifics returned in the msg
            success = True
        alamo_msg = alamo_log.split("\n")[-3]

        return success, alamo_object, alamo_msg

    # TODO: let's generalize this under the metrics?
    def get_alamo_results(self):
        return self._results

    def _get_files(self):
        """
        Method to get/set paths for .alm and .trc files based on filename
        configuration argument.

        If filename is None, temporary files will be created.

        Args:
            None

        Returns:
            None
        """
        if self._temp_context is None:
            self._temp_context = TempfileManager.new_context()

        if self.config.filename is None:
            # Get a temporary file from the manager
            almfile = self._temp_context.create_tempfile(suffix=".alm")
        else:
            almfile = self.config.filename

            if not self.config.overwrite_files:
                # It is OK if the trace file exists, as ALAMO will append to it
                # The trace file reader handles this case
                if os.path.isfile(almfile):
                    raise FileExistsError(
                        f"A file with the name {almfile} already exists. "
                        f"Either choose a new file name or set "
                        f"overwrite_files = True."
                    )

            self._temp_context.add_tempfile(almfile, exists=False)

        trcfile = os.path.splitext(almfile)[0] + ".trc"
        self._temp_context.add_tempfile(trcfile, exists=False)

        if self.config.working_directory is None:
            wrkdir = self._temp_context.create_tempdir()
        else:
            wrkdir = self.config.working_directory

        # Set attributes to track file names
        self._almfile = almfile
        self._trcfile = trcfile
        self._wrkdir = wrkdir

    def _verify_inputs(self):
        """
        Check for issues in input arguments that will cause issues later.
        """
        # ALAMO does not support spaces in variable names
        # TODO: Should this be a check for all IDAES surrogate models?
        for i in self._input_labels:
            if " " in i:
                raise ValueError(
                    f"ALAMO does not support the presence of spaces in variable names. "
                    f"Found space in input names: {i}."
                )
        for o in self._output_labels:
            if " " in o:
                raise ValueError(
                    f"ALAMO does not support the presence of spaces in variable names. "
                    f"Found space in output names: {o}."
                )

    def _write_alm_to_stream(
        self, stream, trace_fname=None, training_data=None, validation_data=None
    ):
        """
        Method to write an ALAMO input file (.alm) to a stream.
        Users may provide specific data sets for training and validation.
        If no data sets are provided, the data sets contained in the
        AlamoModelTrainer are used.

        Args:
            stream: stream that data should be writen to
            trace_fname: name for trace file (.trc) to be included in .alm file
            training_data: Pandas dataframe to use for training surrogate
            validation_data: Pandas dataframe to use for validating surrogate

        Returns:
            None
        """
        if training_data is None:
            training_data = self._training_dataframe
        if validation_data is None:
            validation_data = self._validation_dataframe

        # Check bounds on inputs to avoid potential ALAMO failures
        input_max = list()
        input_min = list()
        for k, b in self._input_bounds.items():
            if b is None or b[0] is None or b[1] is None:
                raise ConfigurationError(
                    f"ALAMO configuration error: invalid bounds on input {k} " f"({b})."
                )
            elif b[0] == b[1]:
                raise ConfigurationError(
                    f"ALAMO configuration error: upper and lower bounds on "
                    f"input {k} are equal."
                )
            elif b[1] < b[0]:
                raise ConfigurationError(
                    f"ALAMO configuration error: upper bound is less than "
                    f"lower bound for input {k}."
                )
            input_min.append(b[0])
            input_max.append(b[1])

        # Get number of data points to build alm file
        n_rdata, n_inputs = training_data.shape

        if validation_data is not None:
            n_vdata, n_inputs = validation_data.shape
        else:
            n_vdata = 0

        stream.write("# IDAES Alamopy input file\n")
        stream.write(f"NINPUTS {len(self._input_labels)}\n")
        stream.write(f"NOUTPUTS {len(self._output_labels)}\n")
        stream.write(f"XLABELS {' '.join(map(str, self._input_labels))}\n")
        stream.write(f"ZLABELS {' '.join(map(str, self._output_labels))}\n")
        stream.write(f"XMIN {' '.join(map(str, input_min))}\n")
        stream.write(f"XMAX {' '.join(map(str, input_max))}\n")
        stream.write(f"NDATA {n_rdata}\n")
        if validation_data is not None:
            stream.write(f"NVALDATA {n_vdata}\n")
        stream.write("\n")

        # Other options for config
        # Can be bool, list of floats, None, Enum, float
        # Special cases
        if self.config.monomialpower is not None:
            stream.write(f"MONO {len(self.config.monomialpower)}\n")
        if self.config.multi2power is not None:
            stream.write(f"MULTI2 {len(self.config.multi2power)}\n")
        if self.config.multi3power is not None:
            stream.write(f"MULTI3 {len(self.config.multi3power)}\n")
        if self.config.ratiopower is not None:
            stream.write(f"RATIOS {len(self.config.ratiopower)}\n")
        if self.config.custom_basis_functions is not None:
            stream.write(f"NCUSTOMBAS {len(self.config.custom_basis_functions)}\n")

        for o in supported_options:
            if self.config[o] is None:
                # Write nothing to alm file
                continue
            elif isinstance(self.config[o], Enum):
                # Need to write value of Enum
                stream.write(f"{o} {self.config[o].value}\n")
            elif isinstance(self.config[o], bool):
                # Cast bool to int
                stream.write(f"{o} {int(self.config[o])}\n")
            elif isinstance(self.config[o], (str, float, int)):
                # Write value to file
                stream.write(f"{o} {self.config[o]}\n")
            else:
                # Assume the argument is a list
                stream.write(f"{o} {' '.join(map(str, self.config[o]))}\n")

        stream.write("\nTRACE 1\n")
        if trace_fname is not None:
            stream.write(f"TRACEFNAME {trace_fname}\n")

        def _df_to_data_fragment(df, **kwargs):
            text = df.to_string(
                header=False,
                index=False,
                float_format=lambda x: str(x).format(":g"),
                **kwargs,
            )
            return text

        stream.write("\nBEGIN_DATA\n")
        # Columns will be writen in order in input and output lists
        training_data_str = _df_to_data_fragment(
            training_data,
            columns=self._input_labels + self._output_labels,
        )
        stream.write(training_data_str)
        stream.write("\nEND_DATA\n")

        if validation_data is not None:
            # Add validation data defintion
            stream.write("\nBEGIN_VALDATA\n")
            val_data_str = _df_to_data_fragment(
                validation_data,
                columns=self._input_labels + self._output_labels,
            )
            stream.write(val_data_str)
            stream.write("\nEND_VALDATA\n")

        if self.config.custom_basis_functions is not None:
            stream.write("\nBEGIN_CUSTOMBAS\n")
            for i in self.config.custom_basis_functions:
                stream.write(f"{str(i)}\n")
            stream.write("END_CUSTOMBAS\n")

    def _write_alm_file(self, training_data=None, validation_data=None):
        """
        Method to write an ALAMO input file (.alm) using the current settings.
        Users may provide specific data sets for training and validation.
        If no data sets are provided, the data sets contained in the
        AlamoModelTrainer are used.

        Args:
            training_data: Pandas dataframe to use for training surrogate
            validation_data: Pandas dataframe to use for validating surrogate

        Returns:
            None
        """
        f = open(self._almfile, "w")
        self._write_alm_to_stream(f, self._trcfile, training_data, validation_data)
        f.close()

    def _call_alamo(self):
        """
        Method to call ALAMO executable from Python, passing the current .alm
        file as an argument.

        Args:
            None

        Returns:
            ALAMO: return code
            log: string of the text output from ALAMO
        """
        if alamo.executable is None:
            raise FileNotFoundError(
                "Could not find ALAMO executable. Please ensure that ALAMO "
                "is installed and in the system path, or provide a path to "
                "the executable."
            )

        ostreams = [StringIO(), sys.stdout]

        if self._temp_context is None:
            self._get_files()

        # Set working directory
        cwd = os.getcwd()
        os.chdir(self._wrkdir)

        try:
            # Add lst file to temp file manager
            lstfname = os.path.splitext(os.path.basename(self._almfile))[0] + ".lst"
            lstpath = os.path.join(cwd, lstfname)
            self._temp_context.add_tempfile(lstpath, exists=False)

            with TeeStream(*ostreams) as t:
                results = subprocess.run(
                    [alamo.executable, str(self._almfile)],
                    stdout=t.STDOUT,
                    stderr=t.STDERR,
                    universal_newlines=True,
                )

                t.STDOUT.flush()
                t.STDERR.flush()

            return_code = results.returncode
            alamo_log = ostreams[0].getvalue()

        except OSError:
            _log.error(
                f"Could not execute the command: alamo {str(self._almfile)}. ",
                f"Error message: {sys.exc_info()[1]}.",
            )
            raise

        finally:
            # Revert cwd to where it started
            os.chdir(cwd)

        if "ALAMO terminated with termination code " in alamo_log:
            self._remove_temp_files()
            _log.warn(
                "ALAMO executable returned non-zero return code. Check "
                "the ALAMO output for more information."
            )

        return return_code, alamo_log

    def _read_trace_file(self, trcfile, has_validation_data=False):
        """
        Method to read the results of an ALAMO run from a trace (.trc) file.
        The name location of the trace file is tored on the AlamoModelTrainer
        object and is generally set automatically.

        Args:
            trcfile : str
               Path to the trcfile to read
            output_labels : list of str
               List of strings of the output_labels (in order)
            has_validation_data : bool
                Bool indicating whether valdiation data was included in ALAMO run

        Returns:
            trace_dict: contents of trace file as a dict
        """
        try:
            with open(trcfile, "r") as f:
                lines = f.readlines()
            f.close()
        except FileNotFoundError:
            # Trace file does not exist
            raise FileNotFoundError(
                "Error occured when trying to read the ALAMO trace file - this probably "
                "indicates that a trace file was not created by the ALAMO executable. "
                "Please check the ALAMO output logs."
            )

        output_labels = self.output_labels()

        trace_read = {}
        # Get headers from first line in trace file
        headers = lines[0].split(", ")
        for i in range(len(headers)):
            header = headers[i].strip("#\n")
            if header in common_trace:
                trace_read[header] = None
            else:
                trace_read[header] = {}

        # Get trace output from final line(s) of file
        # ALAMO will append new lines to existing trace files
        # For multiple outputs, each output has its own line in trace file
        if self._validation_dataframe is not None or has_validation_data:
            omult = 2
        else:
            omult = 1

        for j in range(len(output_labels)):
            trace = lines[(-len(output_labels) + j) * omult].split(", ")

            for i in range(len(headers)):
                header = headers[i].strip("#\n")
                trace_val = trace[i].strip("\n")

                # Replace Fortran powers (^) with Python powers (**)
                trace_val = trace_val.replace("^", "**")
                # Replace = with ==
                trace_val = trace_val.replace("=", "==")

                if header in common_trace:
                    # These should be common across all output
                    if trace_read[header] is None:
                        # No value yet, so set value
                        trace_read[header] = trace_val
                    else:
                        # Check that current value matches the existng value
                        if trace_read[header] != trace_val:
                            raise RuntimeError(
                                f"Mismatch in values when reading ALAMO trace "
                                f"file. Values for {header}: "
                                f"{trace_read[header]}, {header}"
                            )
                else:
                    trace_read[header][output_labels[j]] = trace_val

                # Do some final sanity checks
                if header == "OUTPUT":
                    # OUTPUT should be equal to the current index of outputs
                    if trace_val != str(j + 1):
                        raise RuntimeError(
                            f"Mismatch when reading ALAMO trace file. "
                            f"Expected OUTPUT = {j+1}, found {trace_val}."
                        )
                elif header == "SET":
                    # SET should always be 0 - higher numbers are for
                    # validation data sets
                    if trace_val != "0":
                        raise RuntimeError(
                            f"Mismatch when reading ALAMO trace file. "
                            f"Expected SET = 0, found {trace_val}. "
                            f"This likely indicates the presence of "
                            f"unexpected validation data sets."
                        )
                elif header == "Model":
                    # Var label on LHS should match output label
                    vlabel = trace_val.split("==")[0].strip()
                    if vlabel != output_labels[j]:
                        raise RuntimeError(
                            f"Mismatch when reading ALAMO trace file. "
                            f"Label of output variable in expression "
                            f"({vlabel}) does not match expected label "
                            f"({output_labels[j]})."
                        )

        return trace_read

    def _populate_results(self, trace_dict):
        """
        Method to populate the results object with data from a trace file.

        Args:
            trace_dict: trace file data in form of a dict

        Returns:
            None
        """
        self._results = trace_dict

    def _build_surrogate_object(self):
        """
        Method to construct an AlmaoObject from the current results
        object.

        Args:
            None

        Returns:
            AlamoSurrogate
        """
        return AlamoSurrogate(
            surrogate_expressions=self._results["Model"],
            input_labels=self._input_labels,
            output_labels=self._output_labels,
            input_bounds=self._input_bounds,
        )

    def _remove_temp_files(self):
        """
        Method to remove temporary files created during the ALAMO workflow,
        i.e. the .alm and .trc files.

        Args:
            None

        Returns:
            None
        """
        if self._temp_context is None:
            # Assume we have already cleaned up
            return

        remove = True
        if self.config.filename is not None:
            remove = False

        self._temp_context.release(remove=remove)
        # Release tempfile context
        self._temp_context = None


class AlamoSurrogate(SurrogateBase):
    """
    Standard SurrogateObject for surrogates trained using ALAMO.

    Contains methods to both populate a Pyomo Block with constraints
    representing the surrogate and to evaluate the surrogate a set of user
    provided points.
    """

    def __init__(
        self, surrogate_expressions, input_labels, output_labels, input_bounds=None
    ):
        super().__init__(input_labels, output_labels, input_bounds)
        self._surrogate_expressions = surrogate_expressions
        self._fcn = None

    def evaluate_surrogate(self, inputs):
        """
        Method to method to evaluate the ALAMO surrogate model at a set of user
        provided values.

        Args:
           dataframe: pandas DataFrame
              The dataframe of input values to be used in the evaluation. The dataframe
              needs to contain a column corresponding to each of the input labels. Additional
              columns are fine, but are not used.

        Returns:
            output: pandas Dataframe
              Returns a dataframe of the the output values evaluated at the provided inputs.
              The index of the output dataframe should match the index of the provided inputs.
        """
        # Create a set of lambda functions for evaluating the surrogate.
        if self._fcn is None:
            fcn = dict()
            for o in self._output_labels:
                fcn[o] = eval(
                    f"lambda {', '.join(self._input_labels)}: "
                    f"{self._surrogate_expressions[o].split('==')[1]}",
                    GLOBAL_FUNCS,
                )
            self._fcn = fcn

        # Use numpy to do the calculations as it is faster
        inputdata = inputs[self._input_labels].to_numpy()
        outputs = np.zeros(shape=(inputs.shape[0], len(self._output_labels)))

        for i in range(inputdata.shape[0]):
            for o in range(len(self._output_labels)):
                o_name = self._output_labels[o]
                outputs[i, o] = value(self._fcn[o_name](*inputdata[i, :]))

        return pd.DataFrame(
            data=outputs, index=inputs.index, columns=self._output_labels
        )

    def populate_block(self, block, additional_options=None):
        """
        Method to populate a Pyomo Block with surrogate model constraints.

        Args:
            block: Pyomo Block component to be populated with constraints.
            additional_options: None
               No additional options are required for this surrogate object
        Returns:
            None
        """

        # TODO: do we need to add the index_set stuff back in?
        output_set = Set(initialize=self._output_labels, ordered=True)

        def alamo_rule(b, o):
            lvars = block.input_vars_as_dict()
            lvars.update(block.output_vars_as_dict())
            return eval(self._surrogate_expressions[o], GLOBAL_FUNCS, lvars)

        block.alamo_constraint = Constraint(output_set, rule=alamo_rule)

    def save(self, strm):
        """
        Save an instance of this surrogate to the strm so the model can be used later.

        Args:
           strm: IO.TextIO
              This is the python stream like a file object or StringIO that will be used
              to serialize the surrogate object. This method writes a string
              of json data to the stream.
        """
        json.dump(
            {
                "surrogate": self._surrogate_expressions,
                "input_labels": self._input_labels,
                "output_labels": self._output_labels,
                "input_bounds": self._input_bounds,
            },
            strm,
        )

    @classmethod
    def load(cls, strm):
        """
        Create an instance of a surrogate from a stream.

        Args:
           strm: stream
              This is the python stream containing the data required to load the surrogate.
              This is often, but does not need to be a string of json data.

        Returns: an instance of the derived class or None if it failed to load
        """
        d = json.load(strm)

        surrogate_expressions = d["surrogate"]
        input_labels = d["input_labels"]
        output_labels = d["output_labels"]

        # Need to convert list of bounds to tuples
        input_bounds = {}
        for k, v in d["input_bounds"].items():
            input_bounds[k] = tuple(v)

        return AlamoSurrogate(
            surrogate_expressions=surrogate_expressions,
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=input_bounds,
        )
