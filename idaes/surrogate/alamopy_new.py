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
import os

from idaes.surrogate.my_surrogate_base import Surrogate, SurrogateModelObject
from idaes.core.util.config import list_of_ints, list_of_floats
from pyomo.environ import Constraint
from pyomo.common.config import ConfigValue, In


# TODO: Test calling ALAMO
# TODO: Default file and path names
# TODO: Pre-existing files
# TODO: Custom basis functions
# TODO: Generate expression from string representation
# TODO: Log ALAMO output


DEFAULTPATH = os.path.dirname(__file__)
DEFAULTFNAME = "alamopy.alm"


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
    "solvemip"]

"""
Excluded Options:
NSAMPLE, NVALSAMPLE, INITIALIZER, PRINT TO FILE, FUNFORM, NTRANS
PRINT TO SCREEN - try to handle through logger

Excluded Simulator options:
MAXSIM, MINPOINTS, MAXPOINTS, SAMPLER, SIMULATOR, PRESET, TOLMAXERROR, SIMIN,
SIMOUT

Other options not yet included:
NCUSTOMBAS
"""


# Headers from ALAMO trace file that should be common for all outputs
common_trace = [
    "filename", "NINPUTS", "NOUTPUTS", "INITIALPOINTS", "SET", "INITIALIZER",
    "SAMPLER", "MODELER", "BUILDER", "GREEDYBUILD", "BACKSTEPPER",
    "GREEDYBACK", "REGULARIZER", "SOLVEMIP"]


class Alamopy(Surrogate):
    """
    Standard SurrogateModelTrainer for ALAMO.

    This mainly defines a set of configuration options for ALAMO along with
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

    CONFIG = Surrogate.CONFIG()

    CONFIG.declare('xfactor', ConfigValue(
        default=None,
        domain=list_of_floats,
        description="List of scaling factors for input variables."))
    CONFIG.declare('xscaling', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Option to scale input variables.",
        doc="Option to scale input variables. If True and xfactors are not "
        "provided, ALAMO sets XFACTORS equal to the range of each input "
        "variable."))
    CONFIG.declare('scalez', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Option to scale output variables."))

    # Basis function options
    CONFIG.declare('monomialpower', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="Vector of monomial powers considered in basis "
        "functions."))
    CONFIG.declare('multi2power', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="Vector of powers to be considered for pairwise "
        "combinations in basis functions."))
    CONFIG.declare('multi3power', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="Vector of three variable combinations of powers to be "
        "considered as basis functions."))
    CONFIG.declare('ratiopower', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="Vector of ratio combinations of powers to be considered "
        "in the basis functions."))
    CONFIG.declare('constant', ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Include constant basis function if True. Default = True"))
    CONFIG.declare('linfcns', ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Include linear basis functions if True. Default = True"))
    CONFIG.declare('expfcns', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Include exponential basis functions if True."))
    CONFIG.declare('logfcns', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Include logarithmic basis functions if True."))
    CONFIG.declare('sinfcns', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Include sine basis functions if True."))
    CONFIG.declare('cosfcns', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Include cosine basis functions if True."))
    CONFIG.declare('grbfcns', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Include Gaussian radial basis functions if True."))
    CONFIG.declare('rbfparam', ConfigValue(
        default=None,
        domain=float,
        description="Multiplicative constant used in the Gaussian radial basis"
        " functions."))

    # Other fitting options
    CONFIG.declare('modeler', ConfigValue(
        default=None,
        domain=In(Modelers),
        description="Fitness metric to be used for model building. Must be an "
        "instance of Modelers Enum."))
    CONFIG.declare('builder', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="If True, a greedy heuristic builds up a model "
        "by adding one variable at a time.",
        doc="If True, a greedy heuristic builds up a model by adding one "
        "variable at a time. This model is used as a starting point for "
        "solving an integer programming formulation according to the choice "
        "of modeler. If an optimizer is not available, the heuristic model "
        "will be the final model to be returned."))
    CONFIG.declare('backstepper', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="If set to 1, a greedy heuristic builds down a model by "
        "starting from the least squares model and removing one variable at "
        "a time."))
    CONFIG.declare('convpen', ConfigValue(
        default=None,
        domain=float,
        description="Convex penalty term to use if Modeler == SSEP or MADp.",
        doc="When MODELER is set to 6 or 8, a penalty consisting of the sum "
        "of square errors (SSEP) or the maximum absolute error (MADp) and a "
        "term penalizing model size is used for model building. The size of "
        "the model is weighted by convpen. If convpen=0, this metric reduces "
        "to the classical sum of square errors (SSEP) or the maximum absolute "
        "deviation (MADp)."))
    CONFIG.declare('screener', ConfigValue(
        default=None,
        domain=In(Screener),
        description="Regularization method used to reduce the number of "
        "potential basis functions before optimization. Must be instance of "
        "Screener Enum."))
    CONFIG.declare('ncvf', ConfigValue(
        default=None,
        domain=int,
        description="Number of folds to be used for cross validation by the "
        "lasso screener. ALAMO will use a two-fold validation if fewer than "
        "10 data points are available. NCVF must be a nonnegative integer."))
    CONFIG.declare('sismult', ConfigValue(
        default=None,
        domain=int,
        description="This parameter must be non-negative and is used to "
        "determine the number of basis functions retained by the SIS "
        "screener. The number of basis functions retained equals the floor "
        "of SSISmult n ln(n), where n is the number of measurements "
        "available at the current ALAMO iteration."))

    CONFIG.declare('maxiter', ConfigValue(
        default=None,
        domain=int,
        description="Maximum number of ALAMO iterations. 1 = no adaptive "
        "sampling, 0 = no limit.",
        doc="Maximum number of ALAMO iterations. Each iteration begins with "
        "a model-building step. An adaptive sampling step follows if maxiter "
        "does not equal 1. If maxiter is set to a number less than or equal "
        "to 0, ALAMO will enforce no limit on the number of iterations."))
    CONFIG.declare('maxtime', ConfigValue(
        default=1000,
        domain=float,
        description="Maximum total execution time allowed in seconds. "
        "Default = 1000.",
        doc="Maximum total execution time allowed in seconds. This time "
        "includes all steps of the algorithm, including time to read problem, "
        "preprocess data, solve optimization subproblems, and print results."))
    CONFIG.declare('datalimitterms', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Limit model terms to number of measurements.",
        doc="If True, ALAMO will limit the number of terms in the model to be "
        "no more than the number of data measurements; otherwise, no limit "
        "based on the number of data measurements will be placed. The user "
        "may provide an additional limit on the number of terms in the model "
        "through the maxterms and minterms options."))
    CONFIG.declare('maxterms', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="List of maximum number of model terms to per output.",
        doc="Row vector of maximum terms allowed in the modeling of output "
        "variables. One per output variable, space separated. A −1 signals "
        "that no limit is imposed."))
    CONFIG.declare('minterms', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="List of minimum number of model terms to per output.",
        doc="Row vector of minimum terms required in the modeling of output "
        "variables. One per output variable, space separated. A 0 signals "
        "that no limit is imposed."))
    CONFIG.declare('numlimitbasis', ConfigValue(
        default=True,  # default this to true to avoid numerical issues
        domain=In([True, False]),
        description="Eliminate infeasible basis functions. Default = True",
        doc="If True, ALAMO will eliminate basis functions that are not "
        "numerically acceptable (e.g., log(x) will be eliminated if x may be "
        "negative); otherwise, no limit based on the number of data "
        "measurements will be placed. The user may provide additional limits "
        "on the the type and number of selected basis functions through the "
        "options exclude and groupcon."))
    CONFIG.declare('exclude', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="List of inputs to exclude during building,",
        doc="Row vector of 0/1 flags that specify which input variables, if "
        "any, ALAMO should exclude during the model building process. All "
        "input variables must be present in the data but ALAMO will not "
        "include basis functions that involve input variables for which "
        "exclude equals 1. This feature does not apply to custom basis "
        "functions or RBFs."))
    CONFIG.declare('ignore', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="List of outputs to ignore during building.",
        doc="Row vector of 0/1 flags that specify which output variables, "
        "if any, ALAMO should ignore. All output variables must be present in "
        "the data but ALAMO does not model output variables for which ignore "
        "equals 1."))
    CONFIG.declare('xisint', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="List of inputs that should be treated as integers.",
        doc="Row vector of 0/1 flags that specify which input variables, if "
        "any, ALAMO should treat as integers. For integer inputs, ALAMO’s "
        "sampling will be restricted to integer values."))
    CONFIG.declare('zisint', ConfigValue(
        default=None,
        domain=list_of_ints,
        description="List of outputs that should be treated as integers.",
        doc="Row vector of 0/1 flags that specify which output variables, if "
        "any, ALAMO should treat as integers. For integer variables, ALAMO’s "
        "model will include the rounding of a function to the nearest integer "
        "(equivalent to the nint function in Fortran.)"))

    CONFIG.declare('tolrelmetric', ConfigValue(
        default=None,
        domain=list_of_floats,
        description="Relative tolerance for outputs.",
        doc="Relative convergence tolerance for the chosen fitness metric for "
        "the modeling of output variables. One per output variable, space "
        "separated. Incremental model building will stop if two consecutive "
        "iterations do not improve the chosen metric by at least this amount."
        ))
    CONFIG.declare('tolabsmetric', ConfigValue(
        default=None,
        domain=list_of_floats,
        description="Absolute tolerance for outputs.",
        doc="Absolute convergence tolerance for the chosen fitness metric for "
        "the modeling of output variables. One per output variable, space "
        "separated. Incremental model building will stop if two consecutive "
        "iterations do not improve the chosen metric by at least this amount."
        ))
    CONFIG.declare('tolmeanerror', ConfigValue(
        default=None,
        domain=list_of_floats,
        description="Convergence tolerance for mean errors in outputs.",
        doc="Row vector of convergence tolerances for mean errors in the "
        "modeling of output variables. One per output variable, space "
        "separated. Incremental model building will stop if tolmeanerror, "
        "tolrelmetric, or tolabsmetric is satisfied."))
    CONFIG.declare('tolsse', ConfigValue(
        default=None,
        domain=float,
        description="Absolute tolerance on SSE",
        doc="Absolute tolerance on sum of square errors (SSE). ALAMO will "
        "terminate if it finds a solution whose SSE is within tolsse from "
        "the SSE of the full least squares problem."))

    CONFIG.declare('mipoptca', ConfigValue(
        default=None,
        domain=float,
        description="Absolute tolerance for MIP."))
    CONFIG.declare('mipoptcr', ConfigValue(
        default=None,
        domain=float,
        description="Relative tolerance for MIP."))
    CONFIG.declare('linearerror', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="If True, a linear objective is used when solving "
        "mixed-integer optimization problems; otherwise, a squared error will "
        "be employed."))
    CONFIG.declare('GAMS', ConfigValue(
        default=None,
        domain=str,
        description="Complete path of GAMS executable (or name if GAMS is in "
        "the user path)."))
    CONFIG.declare('GAMSSOLVER', ConfigValue(
        default=None,
        domain=str,
        description="Name of preferred GAMS solver for solving ALAMO’s "
        "mixed-integer quadratic subproblems. Special facilities have been "
        "implemented in ALAMO and BARON that make BARON the preferred "
        "selection for this option. However, any mixed-integer quadratic "
        "programming solver available under GAMS can be used."))
    CONFIG.declare('solvemip', ConfigValue(
        default=None,
        domain=In([True, False]),
        description="Whether to use an optimizer to solve MIP."))
    # CONFIG.declare('log_output', ConfigValue(
    #     default=False,
    #     domain=In([True, False]),
    #     description="Whether to log ALAMO output."))

    # I/O file options
    CONFIG.declare("temp_path", ConfigValue(
        default=DEFAULTPATH,
        domain=str,
        description="Path for directory to store temp files."))
    CONFIG.declare("filename", ConfigValue(
        default=DEFAULTFNAME,
        domain=str,
        description="File name to use for temp files."))
    CONFIG.declare("keep_files", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Flag indicating whether to retain temp files."))
    CONFIG.declare("overwrite_files", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Flag indicating whether existing files can be "
        "overwritten."))

    def __init__(self, **settings):
        super().__init__(**settings)

    def build_model(self):
        super().build_model()

        # Get paths for temp files
        almname, trcname = self.get_paths()

        # Write .alm file
        self.write_alm_file(almname, trcname)

        # Call ALAMO executable
        self.call_alamo(almname)

        # Read back results
        trace_dict = self.read_trace_file(trcname)

        # Populate results and SurrogateModel object
        self.populate_results(trace_dict)
        self.build_surrogate_model_object()

        # Clean up temporary files if required
        if not self.config.keep_files:
            self.remove_temp_files(almname, trcname)

    def get_paths(self):
        path = self.config.temp_path
        if path is None:
            path = DEFAULTPATH

        path = os.chdir(path)

        almname = self.config.filename
        if almname is None:
            almname = DEFAULTFNAME

        trcname = almname.split(".")[0] + ".trc"

        return almname, trcname

    def write_alm_to_stream(
            self, stream=None, trace_fname=None, x_reg=None,
            z_reg=None, x_val=None, z_val=None):
        """Write the input file for the ALAMO executable to a stream."""
        if stream is None:
            raise TypeError("No stream provided to ALAMO writer.")

        if x_reg is None:
            x_reg = self._rdata_in
        if z_reg is None:
            z_reg = self._rdata_out
        if x_val is None:
            x_val = self._vdata_in
        if z_val is None:
            z_val = self._vdata_out

        # Get number of data points to build alm file
        if x_reg is not None:
            n_inputs, n_rdata = x_reg.shape
        else:
            n_rdata = 0
        if x_val is not None:
            n_inputs, n_vdata = x_val.shape
        else:
            n_vdata = 0

        stream.write("# IDAES Alamopy input file\n")
        stream.write(f"NINPUTS {self._n_inputs}\n")
        stream.write(f"NOUTPUTS {self._n_outputs}\n")
        stream.write(f"XLABELS {' '.join(map(str, self._input_labels))}\n")
        stream.write(f"ZLABELS {' '.join(map(str, self._output_labels))}\n")
        stream.write(f"XMIN {' '.join(map(str, self._input_min))}\n")
        stream.write(f"XMAX {' '.join(map(str, self._input_max))}\n")
        stream.write(f"NDATA {n_rdata}\n")
        stream.write(f"NVALDATA {n_vdata}\n\n")

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
            stream.write(f"TRACEFNAME {trace_fname}")

        # TODO : Set FUNFORM if required

        stream.write("\nBEGIN_DATA\n")
        nin = self._n_inputs
        nout = self._n_outputs
        for row in range(n_rdata):
            stream.write(
                f"{' '.join(map(str, (x_reg[x][row] for x in range(nin))))}"
                f" {' '.join(map(str, (z_reg[z][row] for z in range(nout))))}\n")
        stream.write("END_DATA\n")

        if x_val is not None:
            # Add validation data defintion
            stream.write("\nBEGIN_VALDATA\n")
            for row in range(n_vdata):
                stream.write(
                    f"{' '.join(map(str, (x_val[x][row] for x in range(nin))))}"
                    f" {' '.join(map(str, (z_val[z][row] for z in range(nout))))}\n")
            stream.write("END_VALDATA\n")

        # Custom basis functions if required

    def write_alm_file(self, alm_fname=None, trace_fname=None,
                       x_reg=None, z_reg=None, x_val=None, z_val=None):
        f = open(alm_fname, "x")
        self.write_alm_to_stream(f, trace_fname, x_reg, z_reg, x_val, z_val)
        f.close()

    def call_alamo(self, almname):
        os.system("alamo " + almname)

    def read_trace_file(self, trace_file):
        with open(trace_file, "r") as f:
            lines = f.readlines()
        f.close()

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
        for j in range(self._n_outputs):
            trace = lines[-self._n_outputs+j].split(", ")

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
                                f"{trace_read[header]}, {header}")
                else:
                    trace_read[header][self._output_labels[j]] = trace_val

                # Do some final sanity checks
                if header == "OUTPUT":
                    # OUTPUT should be equal to the current index of outputs
                    if trace_val != str(j+1):
                        raise RuntimeError(
                            f"Mismatch when reading ALAMO trace file. "
                            f"Expected OUTPUT = {j+1}, found {trace_val}.")
                elif header == "Model":
                    # Var label on LHS should match output label
                    vlabel = trace_val.split("==")[0].strip()
                    if vlabel != self._output_labels[j]:
                        raise RuntimeError(
                            f"Mismatch when reading ALAMO trace file. "
                            f"Label of output variable in expression "
                            f"({vlabel}) does not match expected label "
                            f"({self._output_labels[j]}).")

        return trace_read

    def populate_results(self, trace_dict):
        self._results = trace_dict

    def build_surrogate_model_object(self):
        self._model = AlamoModelObject(
            surrogate=self._results["Model"],
            input_labels=self._input_labels,
            output_labels=self._output_labels)

    def remove_temp_files(self, almname, trcname):
        os.remove(almname)
        os.remove(trcname)


class AlamoModelObject(SurrogateModelObject):

    def evaluate_surroagte(self, inputs):
        values = {}
        for o in self._output_labels:
            expr = self._surrogate[o].split("==")[1]
            values[o] = eval(expr, locals=inputs)
        return values

    def populate_block(self, block, variables=None):
        if variables is None:
            variables = self.construct_variables(block)

        def alamo_rule(b, o):
            return eval(self._surrogate[o], locals=variables)

        block.add_component(
            "alamo_constraint",
            Constraint(self._output_labels,
                       rule=alamo_rule))
