# -*- coding: utf-8 -*-
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
""" Block-like object meant for models where almost everything is indexed
by time.
"""

from collections import OrderedDict

import idaes.logger as idaeslog
from idaes.apps.caprese.util import initialize_by_element_in_range
from idaes.apps.caprese.common.config import (
    InputOption,
)
from idaes.apps.caprese.common.config import VariableCategory as VC
from idaes.apps.caprese.categorize import (
    categorize_dae_variables,
    CATEGORY_TYPE_MAP,
)
from idaes.apps.caprese.nmpc_var import (
    NmpcVar,
    _NmpcVector,
    DiffVar,
    DerivVar,
    AlgVar,
    InputVar,
    FixedVar,
    MeasuredVar,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from pyomo.environ import (
    Block,
    Reference,
    TransformationFactory,
    Var,
    Set,
    ComponentUID,
    Suffix,
)
from pyomo.core.base.initializer import Initializer
from pyomo.core.base.block import _BlockData, SubclassOf
from pyomo.core.base.indexed_component import UnindexedComponent_set
from pyomo.common.collections import ComponentMap
from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.core.base.range import remainder
from pyomo.dae.set_utils import deactivate_model_at
from pyomo.dae.flatten import flatten_dae_components


class _DynamicBlockData(_BlockData):
    """This class adds methods and data structures that are useful
    for working with dynamic models. These include methods for
    initialization and references to time-indexed variables.
    """

    # TODO: This class should probably give the option to clone
    # the user's model.

    logger = idaeslog.getLogger("nmpc")

    CONFIG = ConfigDict()

    CONFIG.declare(
        "tee",
        ConfigValue(
            default=True,
            domain=bool,
            doc="tee option for embedded solver calls",
        ),
    )

    CONFIG.declare(
        "outlvl",
        ConfigValue(
            default=idaeslog.INFO,
            doc="Output level for IDAES logger",
        ),
    )

    def _construct(self):
        """Generates time-indexed references and categorizes them."""
        model = self.mod
        time = self.time
        try:
            inputs = self._inputs
        except AttributeError:
            inputs = self._inputs = None
        try:
            measurements = self._measurements
        except AttributeError:
            measurements = self._measurements = None

        # TODO: Give the user the option to provide their own
        # category_dict (they know the structure of their model
        # better than I do...)
        if self._category_dict is not None:
            category_dict = self.category_dict = self._category_dict
            if VC.INPUT not in category_dict and inputs is not None:
                self.category_dict[VC.INPUT] = inputs
            if VC.MEASUREMENT not in category_dict and measurements is not None:
                self.category_dict[VC.MEASUREMENT] = measurements
            self.dae_vars = []
            for categ, varlist in category_dict.items():
                if categ is not VC.MEASUREMENT:
                    # Assume that measurements are duplicates
                    self.dae_vars.extend(varlist)

        else:
            scalar_vars, dae_vars = flatten_dae_components(
                model,
                time,
                ctype=Var,
            )
            self.scalar_vars = scalar_vars
            self.dae_vars = dae_vars
            category_dict = categorize_dae_variables(
                dae_vars,
                time,
                inputs,
                measurements=measurements,
            )
            self.category_dict = category_dict

        keys = list(category_dict)
        for categ in keys:
            if not category_dict[categ]:
                # Empty categories cause problems for us down the road
                # due empty (unknown dimension) slices.
                category_dict.pop(categ)

        self._add_category_blocks()
        self._add_category_references()

        self.categories = set(category_dict)

        # Now that we don't rely on knowing what the categories
        # will be, we do not need these attributes. They are, however,
        # used in the tests, so for now, we add them if possible.
        if VC.DIFFERENTIAL in category_dict:
            self.differential_vars = category_dict[VC.DIFFERENTIAL]
        if VC.ALGEBRAIC in category_dict:
            self.algebraic_vars = category_dict[VC.ALGEBRAIC]
        if VC.DERIVATIVE in category_dict:
            self.derivative_vars = category_dict[VC.DERIVATIVE]
        if VC.INPUT in category_dict:
            self.input_vars = category_dict[VC.INPUT]
        if VC.FIXED in category_dict:
            self.fixed_vars = category_dict[VC.FIXED]
        if VC.MEASUREMENT in category_dict:
            self.measurement_vars = category_dict.pop(VC.MEASUREMENT)

        # The categories in category_dict now form a partition of the
        # time-indexed variables. This is necessary to have a well-defined
        # vardata map, which maps each vardata to a unique component indexed
        # only by time.

        # Maps each vardata (of a time-indexed var) to the NmpcVar
        # that contains it.
        self.vardata_map = ComponentMap(
            (var[t], var)
            for var in self.component_objects(SubclassOf(NmpcVar))
            # for varlist in category_dict.values()
            # for var in varlist
            for t in time
            if var.ctype is not MeasuredVar
        )
        # NOTE: looking up var[t] instead of iterating over values()
        # appears to be ~ 5x faster

        # These should be overridden by a call to `set_sample_time`
        # The defaults assume that the entire model is one sample.
        self.sample_points = [time.first(), time.last()]
        self.sample_point_indices = [1, len(time)]

    _var_name = "var"
    _block_suffix = "_BLOCK"
    _set_suffix = "_SET"

    @classmethod
    def get_category_block_name(cls, categ):
        """Gets block name from name of enum entry"""
        return categ.name + cls._block_suffix

    @classmethod
    def get_category_set_name(cls, categ):
        """Gets set name from name of enum entry"""
        return categ.name + cls._set_suffix

    def _add_category_blocks(self):
        """Adds an indexed block for each category of variable and
        attach a reference to each variable to one of the BlockDatas.
        """
        category_dict = self.category_dict
        var_name = self._var_name
        for categ, varlist in category_dict.items():
            ctype = CATEGORY_TYPE_MAP.get(categ, NmpcVar)
            # These names are e.g. 'DIFFERENTIAL_BLOCK', 'DIFFERENTIAL_SET'
            # They serve as a way to access all the "differential variables"
            block_name = self.get_category_block_name(categ)
            set_name = self.get_category_set_name(categ)
            set_range = range(len(varlist))

            # Construct a set that indexes, eg, the "differential variables"
            category_set = Set(initialize=set_range)
            self.add_component(set_name, category_set)

            # Construct an IndexedBlock, each data object of which
            # will contain a single reference-to-timeslice of that
            # category, and with the corresponding custom ctype
            category_block = Block(category_set)
            self.add_component(block_name, category_block)

            # Don't want these blocks sent to any solver.
            category_block.deactivate()

            for i, var in enumerate(varlist):
                # Add reference-to-timeslices to new blocks:
                #
                # Create a new reference if var is not a reference or
                # it has the wrong ctype...
                # We may want to reconsider overriding the user's ctype...
                # We may also want to make sure these references are not
                # attached to anything. Otherwise we will fail here.
                ref = (
                    var
                    if var.is_reference() and var.ctype is ctype
                    else Reference(var, ctype=ctype)
                )
                category_block[i].add_component(var_name, ref)
                # These vars were created by the categorizer
                # and have custom ctypes.

    _vectors_name = "vectors"

    def _add_category_references(self):
        """Create a "time-indexed vector" for each category of variables."""
        category_dict = self.category_dict

        # Add a deactivated block to store all my `_NmpcVector`s
        # These be will vars, named by category, indexed by the index
        # into the list of that category and by time. E.g.
        # self.vectors.differential
        self.add_component(self._vectors_name, Block())
        self.vectors.deactivate()

        for categ in category_dict:
            ctype = CATEGORY_TYPE_MAP.get(categ, NmpcVar)
            # Get the block that holds this category of var,
            # and the name of the attribute that holds the
            # custom-ctype var (this attribute is the same
            # for all blocks).
            block_name = self.get_category_block_name(categ)
            var_name = self._var_name

            # Get a slice of the block, e.g. self.DIFFERENTIAL_BLOCK[:]
            block = getattr(self, block_name)
            _slice = block[:]
            # _slice = self.__getattribute__(block_name)[:]
            # Why does this work when self.__getattr__(block_name) does not?
            # __getattribute__ appears to work just fine...

            # Get a slice of the block and var, e.g.
            # self.DIFFERENTIAL_BLOCK[:].var[:]
            _slice = getattr(_slice, var_name)[:]

            # Add a reference to this slice to the `vectors` block.
            # This will be, e.g. `self.vectors.differential` and
            # can be accessed with its two indices, e.g.
            # `self.vectors.differential[i,t0]`
            # to get the "ith coordinate" of the vector of differential
            # variables at time t0.
            ref = Reference(_slice, ctype=_NmpcVector)
            self.vectors.add_component(
                categ.name.lower(),  # Lowercase of the enum name
                ref,
            )

    # Time is added in DynamicBlock.construct but this is nice if the user wants
    # to add time in a rule without a long messy line of super().__setattr__.
    def add_time(self):
        # Do this because I can't add a reference to a set
        super(_BlockData, self).__setattr__("time", self.time)

    def set_sample_time(self, sample_time, tolerance=1e-8):
        """Validates and sets sample time"""
        self.validate_sample_time(sample_time, tolerance)
        self.sample_time = sample_time

    def validate_sample_time(self, sample_time, tolerance=1e-8):
        """Makes sure sample points, or integer multiple of sample time-offsets
        from time.first(), lie on finite element boundaries, and that the
        horizon of each model is an integer multiple of sample time. Assembles
        a list of sample points and a dictionary mapping sample points to the
        number of finite elements in the preceding sampling period, and adds
        them as attributes to _NMPC_NAMESPACE.

        Args:
            sample_time: Sample time to check
            tolerance: Tolerance within which time points must be integer
                       multiples of sample time

        """
        time = self.time
        horizon_length = time.last() - time.first()
        n_t = len(time)

        # TODO: This should probably be a DAE utility
        min_spacing = horizon_length
        for t in time:
            if t == time.first():
                continue
            prev = time.prev(t)
            if t - prev < min_spacing:
                min_spacing = t - prev
        # Sanity check:
        assert min_spacing > 0
        # Required so only one point can satisfy equality to tolerance
        if tolerance >= min_spacing / 2:
            raise ValueError(
                "ContinuousSet tolerance is larger than half the minimum "
                "spacing. An element of this set will not necessarily be "
                "unique within this tolerance."
            )

        off_by = abs(remainder(horizon_length, sample_time))
        if off_by > tolerance:
            raise ValueError(
                "Sampling time must be an integer divider of "
                "horizon length within tolerance %f" % tolerance
            )
        n_samples = round(horizon_length / sample_time)
        self.samples_per_horizon = n_samples

        finite_elements = time.get_finite_elements()
        fe_set = set(finite_elements)
        finite_element_indices = [i for i in range(1, n_t + 1) if time.at(i) in fe_set]
        sample_points = [time.first()]
        sample_indices = [1]  # Indices of sample points with in time set
        sample_no = 1
        fe_per = 0
        fe_per_sample_dict = {}
        for i, t in zip(finite_element_indices, finite_elements):
            if t == time.first():
                continue
            fe_per += 1
            time_since = t - time.first()
            sp = sample_no * sample_time
            diff = abs(sp - time_since)
            if diff < tolerance:
                sample_points.append(t)
                sample_indices.append(i)
                sample_no += 1
                fe_per_sample_dict[sample_no] = fe_per
                fe_per = 0
            if time_since > sp:
                raise ValueError(
                    "Could not find a time point for the %ith "
                    "sample point" % sample_no
                )
        assert len(sample_points) == n_samples + 1
        self.fe_per_sample = fe_per_sample_dict
        self.sample_points = sample_points
        self.sample_point_indices = sample_indices

    def initialize_sample_to_setpoint(
        self,
        sample_idx,
        ctype=(DiffVar, AlgVar, InputVar, DerivVar),
    ):
        """Set values to setpoint values for variables of the
        specified variable ctypes in the specified sample.
        """
        time = self.time
        sample_point_indices = self.sample_point_indices
        i_0 = sample_point_indices[sample_idx - 1]
        i_s = sample_point_indices[sample_idx]
        for var in self.component_objects(ctype):
            # `type(var)` is a subclass of `NmpcVar`, so I can
            # access the `setpoint` attribute.
            #
            # Would like:
            # var[t1:ts].set_value(var.setpoint)
            for i in range(i_0 + 1, i_s + 1):
                # Want to exclude first time point of sample,
                # but include last time point of sample.
                t = time.at(i)
                var[t].set_value(var.setpoint)

    def initialize_sample_to_initial(
        self,
        sample_idx,
        ctype=(DiffVar, AlgVar, DerivVar),
    ):
        """Set values to initial values for variables of the
        specified variable ctypes in the specified sample.
        """
        time = self.time
        sample_point_indices = self.sample_point_indices
        i_0 = sample_point_indices[sample_idx - 1]
        i_s = sample_point_indices[sample_idx]
        t0 = time.at(i_0)
        for var in self.component_objects(ctype):
            # Would be nice if I could use a slice with
            # start/stop indices to make this more concise.
            for i in range(i_0 + 1, i_s + 1):
                t = time.at(i)
                var[t].set_value(var[t0].value)

    def initialize_to_setpoint(
        self,
        ctype=(DiffVar, AlgVar, InputVar, DerivVar),
    ):
        """Sets values to setpoint values for specified variable
        ctypes for all time points.
        """
        # There should be negligible overhead to initializing
        # in many small loops as opposed to one big loop here.
        for i in range(len(self.sample_points)):
            self.initialize_sample_to_setpoint(i, ctype=ctype)

    def initialize_to_initial_conditions(
        self,
        ctype=(DiffVar, AlgVar, DerivVar),
    ):
        """Sets values to initial values for specified variable
        ctypes for all time points.
        """
        # There should be negligible overhead to initializing
        # in many small loops as opposed to one big loop here.
        for i in range(len(self.sample_points)):
            self.initialize_sample_to_initial(i, ctype=ctype)

    def initialize_by_solving_elements(self, solver, **kwargs):
        """Solve the square problem with fixed inputs in each
        of the time finite elements individually. This can be
        thought of as a time integration.
        """
        strip_var_bounds = kwargs.pop("strip_var_bounds", True)
        input_option = kwargs.pop("input_option", InputOption.CURRENT)
        config = self.CONFIG(kwargs)
        square_solve_context = SquareSolveContext(
            self,
            strip_var_bounds=strip_var_bounds,
            input_option=input_option,
        )
        model = self.mod
        time = self.time
        # There is a significant amount of overhead when calling
        # initialize_by_element_in_range multiple times, so
        # this method does not call `initialize_samples_by_element`
        # in a loop.
        if VC.DIFFERENTIAL in self.categories:
            time_linking_vars = self.differential_vars
        else:
            time_linking_vars = None
        with square_solve_context as sqs:
            initialize_by_element_in_range(
                model,
                time,
                time.first(),
                time.last(),
                dae_vars=self.dae_vars,
                # time_linking_vars=list(self.differential_vars[:]),
                time_linking_vars=time_linking_vars,
                outlvl=config.outlvl,
                solver=solver,
            )

    def initialize_samples_by_element(self, samples, solver, **kwargs):
        """Solve the square problem with fixed inputs for the specified samples"""
        # TODO: ConfigBlock for this class
        strip_var_bounds = kwargs.pop("strip_var_bounds", True)
        input_option = kwargs.pop("input_option", InputOption.CURRENT)
        config = self.CONFIG(kwargs)

        if type(samples) not in {list, tuple}:
            samples = (samples,)

        # Create a context manager that will temporarily strip bounds
        # and fix inputs, preparing the model for a "square solve."
        square_solve_context = SquareSolveContext(
            self,
            samples=samples,
            strip_var_bounds=strip_var_bounds,
            input_option=input_option,
        )
        sample_points = self.sample_points
        model = self.mod
        time = self.time
        with square_solve_context as sqs:
            for s in samples:
                t0 = sample_points[s - 1]
                t1 = sample_points[s]
                # Really I would like an `ElementInitializer` context manager
                # class that deactivates the model once, then allows me to
                # activate the elements I want to solve one at a time.
                # This would allow me to not repeat so much work when
                # initializing multiple samples.
                initialize_by_element_in_range(
                    model,
                    time,
                    t0,
                    t1,
                    dae_vars=self.dae_vars,
                    time_linking_vars=list(self.differential_vars[:]),
                    outlvl=config.outlvl,
                    solver=solver,
                )

    def set_variance(self, variance_list):
        """Set variance for corresponding NmpcVars to the values provided

        Arguments:
            variance_list: List of vardata, value tuples. The vardatas
                           correspond to time-indexed references, and values
                           are the variances.
        """
        t0 = self.time.first()
        variance_map = ComponentMap(variance_list)
        for var, val in variance_list:
            nmpc_var = self.vardata_map[var]
            nmpc_var.variance = val
        # MeasurementVars will not have their variance set since they are
        # not mapped to in vardata_map
        for var in self.component_objects(MeasuredVar):
            # component_objects is fine here because we don't need
            # to access measurements in any particular order.
            if var[t0] in variance_map:
                var.variance = variance_map[var[t0]]

    def generate_inputs_at_time(self, t):
        if VC.INPUT in self.categories:
            for val in self.vectors.input[:, t].value:
                yield val
        else:
            raise RuntimeError(
                "Trying to generate inputs but no input " "category has been specified."
            )

    def generate_measurements_at_time(self, t):
        if VC.MEASUREMENT in self.categories:
            for val in self.vectors.measurement[:, t].value:
                yield val
        else:
            raise RuntimeError(
                "Trying to generate measurements but no measurement "
                "category has been specified."
            )

    def inject_inputs(self, inputs):
        # To simulate computational delay, this function would
        # need an argument for the start time of inputs.
        if VC.INPUT in self.categories:
            for var, val in zip(self.INPUT_BLOCK[:].var, inputs):
                # Would like:
                # self.vectors.input[:,:].fix(inputs)
                # This is an example of setting a matrix from a vector.
                # Could even aspire towards:
                # self.vectors.input[:,t0:t1].fix(inputs)
                var[:].fix(val)
        else:
            raise RuntimeError(
                "Trying to set input values but no input "
                "category has been specified."
            )

    def load_measurements(self, measured):
        t0 = self.time.first()
        if VC.MEASUREMENT in self.categories:
            # Want: self.measured_vars[:,t0].fix(measured)
            for var, val in zip(self.MEASUREMENT_BLOCK[:].var, measured):
                var[t0].fix(val)
        else:
            raise RuntimeError(
                "Trying to set measurement values but no measurement "
                "category has been specified."
            )

    def advance_by_time(
        self,
        t_shift,
        ctype=(DiffVar, DerivVar, AlgVar, InputVar, FixedVar),
        # Fixed variables are included as I expect disturbances
        # should shift in time as well.
        tolerance=1e-8,
    ):
        """Set values for the variables of the specified ctypes
        to their values `t_shift` in the future.
        """
        time = self.time
        # The outer loop is over time so we don't have to call
        # `find_nearest_index` for every variable.
        # I am assuming that `find_nearest_index` is slower than
        # accessing `component_objects`
        for t in time:
            ts = t + t_shift
            idx = time.find_nearest_index(ts, tolerance)
            if idx is None:
                # t + sample_time is outside the model's "horizon"
                continue
            ts = time.at(idx)
            for var in self.component_objects(ctype):
                var[t].set_value(var[ts].value)

    def advance_one_sample(
        self,
        ctype=(DiffVar, DerivVar, AlgVar, InputVar, FixedVar),
        tolerance=1e-8,
    ):
        """Set values for the variables of the specified ctypes
        to their values one sample time in the future.
        """
        sample_time = self.sample_time
        self.advance_by_time(
            sample_time,
            ctype=ctype,
            tolerance=tolerance,
        )

    def generate_time_in_sample(
        self,
        ts,
        t0=None,
        include_t0=False,
        tolerance=1e-8,
    ):
        """Generate time points between the provided time point
        and one sample time in the past.
        """
        # TODO: Need to address the question of whether I want users
        # passing around time points or the integer index of samples.
        time = self.time
        idx_s = time.find_nearest_index(ts, tolerance=tolerance)
        ts = time.at(idx_s)
        if t0 is None:
            t0 = ts - self.sample_time
        idx_0 = time.find_nearest_index(t0, tolerance=tolerance)
        idx_start = idx_0 if include_t0 else idx_0 + 1
        for i in range(idx_start, idx_s + 1):
            # Don't want to include first point in sample
            yield time.at(i)

    def get_data_from_sample(
        self,
        ts,
        variables=(
            VC.DIFFERENTIAL,
            VC.INPUT,
        ),
        tolerance=1e-8,
        include_t0=False,
    ):
        """Creates an `OrderedDict` that maps the time-indexed reference
        of each variable provided to a list of its values over the sample
        preceding the specified time point.
        """
        time = self.time
        sample_time = self.sample_time
        category_dict = self.category_dict
        vardata_map = self.vardata_map

        data = OrderedDict()
        queue = list(variables)
        for var in queue:
            if type(var) is VC:
                category = var
                varlist = category_dict[category]
                queue.extend(var[ts] for var in varlist)
                continue
            _slice = vardata_map[var]
            cuid = ComponentUID(_slice.referent)
            if include_t0:
                i0 = time.find_nearest_index(ts - sample_time, tolerance=tolerance)
                t0 = time.at(i0)
                data[cuid] = [_slice[t0].value]
            else:
                data[cuid] = []
            data[cuid].extend(
                _slice[t].value
                for t in self.generate_time_in_sample(ts, tolerance=tolerance)
            )

        return data

    def add_ipopt_suffixes(self):
        """Adds suffixes for communicating dual variables with IPOPT"""
        # Maybe there should be some helper class to do solver-specific
        # stuff like this...
        self.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
        self.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)

        self.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
        self.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

        self.dual = Suffix(direction=Suffix.IMPORT_EXPORT)

    def update_ipopt_multipliers(self):
        self.ipopt_zL_in.update(self.ipopt_zL_out)
        self.ipopt_zU_in.update(self.ipopt_zU_out)

    def advance_ipopt_multipliers(
        self,
        t_shift,
        ctype=(
            DiffVar,
            AlgVar,
            InputVar,
        ),
        tolerance=1e-8,
    ):
        """Set the values of bound multipliers to the corresponding
        values a time `t_shift` in the future.
        """
        zL = self.ipopt_zL_in
        zU = self.ipopt_zU_in
        time = self.time
        # The outer loop is over time so we don't have to call
        # `find_nearest_index` for every variable.
        # I am assuming that `find_nearest_index` is slower than
        # accessing `component_objects`
        for t in time:
            ts = t + t_shift
            idx = time.find_nearest_index(ts, tolerance)
            if idx is None:
                # t + sample_time is outside the model's "horizon"
                continue
            ts = time.at(idx)
            for var in self.component_objects(ctype):
                if var[t] in zL and var[ts] in zL:
                    zL[var[t]] = zL[var[ts]]
                if var[t] in zU and var[ts] in zU:
                    zU[var[t]] = zU[var[ts]]

    def advance_ipopt_multipliers_one_sample(
        self,
        ctype=(
            DiffVar,
            AlgVar,
            InputVar,
        ),
        tolerance=1e-8,
    ):
        """Set the values of bound multipliers to the corresponding
        values one sample time in the future.
        """
        sample_time = self.sample_time
        self.advance_ipopt_multipliers(
            sample_time,
            ctype=ctype,
            tolerance=tolerance,
        )


class DynamicBlock(Block):
    """This is the DynamicBlock class that the user will instantiate."""

    _ComponentDataClass = _DynamicBlockData

    def __new__(cls, *args, **kwds):
        # Decide what class to allocate
        if cls != DynamicBlock:
            target_cls = cls
        elif not args or (args[0] is UnindexedComponent_set and len(args) == 1):
            target_cls = SimpleDynamicBlock
        else:
            target_cls = IndexedDynamicBlock
        return super(DynamicBlock, cls).__new__(target_cls)

    def __init__(self, *args, **kwds):
        # This will get called regardless of what class we are instantiating
        self._init_model = Initializer(kwds.pop("model", None))
        self._init_time = Initializer(
            kwds.pop("time", None), treat_sequences_as_mappings=False
        )
        self._init_inputs = Initializer(
            kwds.pop("inputs", None), treat_sequences_as_mappings=False
        )
        self._init_measurements = Initializer(
            kwds.pop("measurements", None), treat_sequences_as_mappings=False
        )
        self._init_category_dict = Initializer(
            kwds.pop("category_dict", None), treat_sequences_as_mappings=False
        )
        Block.__init__(self, *args, **kwds)

    def _getitem_when_not_present(self, idx):
        block = super(DynamicBlock, self)._getitem_when_not_present(idx)
        parent = self.parent_block()

        if self._init_model is not None:
            block.mod = self._init_model(parent, idx)

        if self._init_time is not None:
            super(_BlockData, block).__setattr__("time", self._init_time(parent, idx))

        if self._init_inputs is not None:
            block._inputs = self._init_inputs(parent, idx)

        if self._init_measurements is not None:
            block._measurements = self._init_measurements(parent, idx)

        if self._init_category_dict is not None:
            block._category_dict = self._init_category_dict(parent, idx)
        else:
            block._category_dict = None

        block._construct()


class SimpleDynamicBlock(_DynamicBlockData, DynamicBlock):
    def __init__(self, *args, **kwds):
        _DynamicBlockData.__init__(self, component=self)
        DynamicBlock.__init__(self, *args, **kwds)

    # Pick up the display() from Block and not BlockData
    display = DynamicBlock.display


class IndexedDynamicBlock(DynamicBlock):
    def __init__(self, *args, **kwargs):
        DynamicBlock.__init__(self, *args, **kwargs)


class SquareSolveContext(object):
    """
    Utility class to prepare DynamicBlock for a square solve.
    """

    def __init__(
        self,
        dynamic_block,
        samples=None,
        strip_var_bounds=True,
        input_option=InputOption.CURRENT,
    ):
        """
        Parameters
        ----------
            dynamic_block: A _DynamicBlockData object
            samples: A list of integers corresponding to the samples that
                     will be solved

        """
        self.block = dynamic_block
        self.samples = samples
        self.strip_var_bounds = strip_var_bounds
        self.input_option = input_option

        # Get indices for the time points where we need to fix inputs.
        if samples is None:
            # Assume we need to fix inputs for all non-initial time
            self.time_indices = range(2, len(dynamic_block.time) + 1)
        else:
            # Get the indices of all time points in the samples specified
            sample_points = dynamic_block.sample_points
            time = dynamic_block.time
            time_indices = []
            already_visited = set()
            for s in samples:
                # s is an integer in range(1, len(sample_points)+1)
                # the ith sample is the interval (ts_{i-1}, ts_i]
                t0 = sample_points[s - 1]
                ts = sample_points[s]
                idx_0 = time.find_nearest_index(t0)
                idx_s = time.find_nearest_index(ts)
                for i in range(idx_0 + 1, idx_s + 1):  # Pyomo sets are 1-indexed
                    if i not in already_visited:
                        # Want to make sure each index gets added at most
                        # once in the case of repeated samples...
                        time_indices.append(i)
                        already_visited.add(i)
            self.time_indices = time_indices

    def __enter__(self):
        # Strip bounds:
        if self.strip_var_bounds:
            self.strip_bounds = TransformationFactory("contrib.strip_var_bounds")
            self.strip_bounds.apply_to(self.block.mod, reversible=True)

        # Fix inputs:
        time = self.block.time
        t0 = time.first()
        time_indices = self.time_indices
        input_vars = self.block.input_vars
        input_option = self.input_option
        if input_option is InputOption.CURRENT:
            input_vals = [None for _ in input_vars]
        elif input_option is InputOption.INITIAL:
            input_vals = [v[t0] for v in input_vars]
        elif input_option is InputOption.SETPOINT:
            input_vals = [v.setpoint for v in input_vars]
        else:
            raise NotImplementedError("Unrecognized input option")

        for var, val in zip(input_vars, input_vals):
            for i in time_indices:
                t = time.at(i)
                if val is None:
                    var[t].fix()
                else:
                    var[t].fix(val)

        # Model should now be square and bound-free
        return self

    def __exit__(self, ex_type, ex_val, ex_tb):
        # Unfix inputs:
        time = self.block.time
        time_indices = self.time_indices
        input_vars = self.block.input_vars
        for var in input_vars:
            for i in time_indices:
                t = time.at(i)
                var[t].unfix()

        if self.strip_var_bounds:
            self.strip_bounds.revert(self.block.mod)
