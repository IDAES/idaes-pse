import time as time_module
from collections import OrderedDict

from six import iteritems

import idaes.logger as idaeslog
from idaes.apps.caprese.base_class import DynamicBase
from idaes.apps.caprese.util import (
        initialize_by_element_in_range,
        NMPCVarLocator,
        cuid_from_timeslice,
        )
from idaes.apps.caprese.common.config import (
        VariableCategory,
        ControlPenaltyType,
        ControlInitOption,
        )
from idaes.apps.caprese.categorize import (
        categorize_dae_variables,
        CATEGORY_TYPE_MAP,
        )
from idaes.apps.caprese.nmpc_var import NmpcVar, _NmpcVector
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (
        find_comp_in_block_at_time,
        )
from pyomo.environ import (
        value,
        Objective,
        TerminationCondition,
        Constraint,
        Block,
        Reference,
        TransformationFactory,
        Var,
        Set,
        )
from pyomo.core.base.util import Initializer, ConstantInitializer
from pyomo.core.base.block import _BlockData, declare_custom_block
from pyomo.core.base.indexed_component import UnindexedComponent_set
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.common.modeling import unique_component_name
from pyomo.common.timing import ConstructionTimer
from pyomo.core.base.range import remainder
from pyomo.dae.set_utils import deactivate_model_at
from pyomo.dae.flatten import flatten_dae_components

class _DynamicBlockData(_BlockData):
    # TODO: This class should probably give the option to clone
    # the user's model.

    def _construct(self):
        # Is it "bad practice" to have a construct method on a data object?
        # ^ Yes, because this method needs to be called via DynamicBlock,
        # not immediately when _DynamicBlockData.construct is invoked...
        model = self.mod
        time = self.time
        inputs = self._inputs
        measurements = self._measurements

        # TODO: Give the user the option to provide their own
        # category_dict (they know the structure of their model
        # better than I do...)
        scalar_vars, dae_vars = flatten_dae_components(
                model,
                time,
                ctype=Var,
                )
        self.dae_vars = dae_vars
        category_dict = categorize_dae_variables(dae_vars, time, inputs)
        self.category_dict = category_dict

        self._add_category_blocks()
        self._add_category_references()

        # Why are these attributes not getting set...?
        self.measured_vars = category_dict.pop(VariableCategory.MEASUREMENT)
        # The categories in category_dict now form a partition of the
        # time-indexed variables. This is necessary to have a well-defined
        # vardata map, which maps each vardata to a unique component indexed
        # only by time.

        # Maps each vardata (of a time-indexed var) to the NmpcVar
        # that contains it.
        self.vardata_map = ComponentMap((var[t], var) 
                for varlist in category_dict.values()
                for var in varlist
                for t in time
                )
        # NOTE: looking up var[t] instead of iterating over values() 
        # appears to be ~ 5x faster

    _var_name = 'var'
    _block_suffix = '_BLOCK'
    _set_suffix = '_SET'

    @classmethod
    def get_category_block_name(cls, categ):
        categ_name = str(categ).split('.')[1]
        return categ_name + cls._block_suffix

    @classmethod
    def get_category_set_name(cls, categ):
        categ_name = str(categ).split('.')[1]
        return categ_name + cls._set_suffix

    def _add_category_blocks(self):
        category_dict = self.category_dict
        var_name = self._var_name
        for categ, varlist in category_dict.items():
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
                category_block[i].add_component(var_name, var)
                # These vars were created by the categorizer
                # and have custom ctypes.

    _vectors_name = 'vectors'

    def _add_category_references(self):
        category_dict = self.category_dict

        # Add a deactivated block to store all my `_NmpcVector`s
        # These will vars, named by category, indexed by the index 
        # into the list of that category and by time. E.g.
        # self.vectors.differential
        self.add_component(self._vectors_name, Block())
        self.vectors.deactivate()

        for categ in category_dict:
            # TODO: use ctype as the key in category_dict
            ctype = CATEGORY_TYPE_MAP[categ]
            # Get the block that holds this category of var,
            # and the name of the attribute that holds the
            # custom-ctype var (this attribute is the same
            # for all blocks).
            block_name = self.get_category_block_name(categ)
            var_name = self._var_name

            # Get a slice of the block, e.g.
            # self.DIFFERENTIAL_BLOCK[:]
            _slice = getattr(self, block_name)[:]
            #_slice = self.__getattribute__(block_name)[:]
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
            self.vectors.add_component(
                    ctype._attr,
                    # ^ I store the name I want this attribute to have,
                    # e.g. 'differential', on the custom ctype.
                    Reference(_slice, ctype=_NmpcVector),
                    )

    # Should be unnecessary here as time is added in DynamicBlock.construct
    # But this is nice if the user wants to add time in a rule and doesn't
    # want a long messy line of super().__setattr__.
    def add_time(self):
        # Do this because I can't add a reference to a set
        super(_BlockData, self).__setattr__('time', self.time)

    def find_components(self, model, comps, time):
        # TODO: Block has a very similar method. This method
        # is different as it is "time-aware."
        # I should probably rewrite this to use CUIDs, as my
        # underlying method is doing something almost equivalent.
        #
        # _slices = [c.referent for c in comps]
        # cuids = [ComponentUID(s, context=src_model) for s in _slices]
        # for cuid in cuids:
        #     for comp in cuid.list_components(self.mod):
        #         break
        #     comp = self.vardata_map[comp]
        #     tgt_comps.append(comp)
        # ^ This seems sufficiently short to do in NmpcManager.
        # Once I do this I can remove this method.
        # If I could get a slice from a CUID, I could do:
        #
        # _slices = [c.referent for c in comps]
        # cuids = [ComponentUID(s, context=src_model) for s in _slices]
        # for cuid in cuids:
        #     _slice = cuid.find_component(self.mod)
        #     tgt_comps.append(Reference(_slice))
        # ^ This won't return components with ctype NmpcVar, but I can
        # always pass in a ctype argument.
        t0_src = time.first()
        t0_tgt = self.time.first()
        tgt_comps = []
        for comp in comps:
            src_vardata = comp[t0_src]
            tgt_vardata = find_comp_in_block_at_time(
                    self.mod,
                    model,
                    src_vardata,
                    self.time,
                    t0_tgt,
                    )
            tgt_comps.append(self.vardata_map[tgt_vardata])
        return tgt_comps

    def set_sample_time(self, sample_time, tolerance=1e-8):
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

        """
        time = self.time
        horizon_length = time.last() - time.first()

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
        if tolerance >= min_spacing/2:
            raise ValueError(
                'ContinuousSet tolerance is larger than half the minimum '
                'spacing. An element of this set will not necessarily be '
                'unique within this tolerance.')

        off_by = abs(remainder(horizon_length, sample_time))
        if off_by > tolerance:
            raise ValueError(
                'Sampling time must be an integer divider of '
                'horizon length within tolerance %f' % tolerance)
        n_samples = round(horizon_length/sample_time)
        self.samples_per_horizon = n_samples

        finite_elements = time.get_finite_elements()
        sample_points = [time.first()]
        sample_no = 1
        fe_per = 0
        fe_per_sample_dict = {}
        for t in finite_elements:
            if t == time.first():
                continue
            fe_per += 1
            time_since = t - time.first()
            sp = sample_no*sample_time
            diff = abs(sp-time_since)
            if diff < tolerance:
                sample_points.append(t)
                sample_no += 1
                fe_per_sample_dict[sample_no] = fe_per
                fe_per = 0
            if time_since > sp:
                raise ValueError(
                        'Could not find a time point for the %ith '
                        'sample point' % sample_no)
        assert len(sample_points) == n_samples + 1
        self.fe_per_sample = fe_per_sample_dict
        self.sample_points = sample_points

    def initialize(self, option):
        time = self.time

        if option == ControlInitOption.FROM_PREVIOUS:
            pass
        elif option == ControlInitOption.BY_TIME_ELEMENT:
            pass
        elif option == ControlInitOption.FROM_INITIAL_CONDITIONS:
            pass
        elif option == ControlInitOption.SETPOINT:
            pass

    def initialize_from_setpoint(self,
            categories=[
                VariableCategory.DIFFERENTIAL,
                VariableCategory.ALGEBRAIC,
                VariableCategory.DERIVATIVE,
                VariableCategory.INPUT,
                ],
            ):
        time = self.time
        t0 = time.first()
        category_dict = self.category_dict
        for categ, varlist in categories.items():
            # Would like:
            # vector[:,t1:].set_value(vector[:].setpoint)
            for var in varlist:
                for t in time:
                    if t == t0:
                        continue
                    var[t].set_value(var.setpoint)

    def initialize_from_initial_conditions(self,
            categories=[
                VariableCategory.DERIVATIVE,
                VariableCategory.DIFFERENTIAL,
                VariableCategory.ALGEBRAIC,
                ],
            ):
        """ 
        Set values of differential, algebraic, and derivative variables to
        their values at the initial conditions.
        An implicit assumption here is that the initial conditions are
        consistent.

        Args:
            model : Flowsheet model whose variables are initialized
            categories : List of VariableCategory enum items to 
                         initialize. Default contains DERIVATIVE, DIFFERENTIAL,
                         and ALGEBRAIC.

        """
        time = self.time
        t0 = time.first()
        cat_dict = self.category_dict
        for categ in categories:
            for v in cat_dict[categ]:
                v[:].set_value(v[t0].value)
        # for var in component_objects(categories):
        #     var[:].set_value(var[t0].value)
        # This should work without a custom var class.

    def initialize_by_solving_elements(self, solver, fix_inputs=False, strip_bounds=False):
        time = self.time
        model = self.model

        if strip_bounds:
            strip_bounds = TransformationFactory('contrib.strip_var_bounds')
            strip_bounds.apply_to(model, reversible=True)

        if fix_inputs:
            input_vars = self.input_vars
            for var in input_vars[:]:
                sp = var.setpoint
                # Would like:
                # var[t1:].fix(sp)
                for t in time:
                    if t != time.first():
                        var[t].fix(sp)
                    else:
                        var[t].fix()
        
        if hasattr(self, 'tracking_objective'):
            self.tracking_objective.deactivate()
        if hasattr(self, 'pwc_constraint'):
            self.pwc_constraint.deactivate()

        initialize_by_element_in_range(
                model,
                time,
                time.first(),
                time.last(),
                dae_vars=self.dae_vars,
                time_linking_vars=list(self.diff_vars[:]),
                outlvl=idaeslog.DEBUG,
                solver=solver,
                )

        if hasattr(self, 'tracking_objective'):
            self.tracking_objective.activate()
        if hasattr(self, 'pwc_constraint'):
            self.pwc_constraint.activate()

        if fix_inputs:
            for var in self.input_vars[:]:
                for t in time:
                    if t == time.first():
                        continue
                    var[t].unfix()

        if strip_bounds:
            strip_bounds.revert(model)
        # TODO: Can check for violated bounds if I'm really worried about it.
        # Should have a general method to deal with violated bounds that I
        # can use for noise as well.

    def initialize_sample_by_element(self, solver, ts, fix_inputs=False):
        time = self.time
        model = self.model
        strip_bounds = TransformationFactory('contrib.strip_var_bounds')
        strip_bounds.apply_to(model, reversible=True)

        t0 = ts - self.sample_time
        idx_0 = time.find_nearest_index(t0)
        t0 = time[idx_0]
        idx_s = time.find_nearest_index(ts)

        if fix_inputs:
            input_vars = self.input_vars
            for var in input_vars[:]:
                sp = var.setpoint
                # Would like:
                # var[t1:].fix(sp)
                for i in range(idx_0+1, idx_s+1):
                    t = time[i]
                    var[t].fix(sp)
        
        if hasattr(self, 'tracking_objective'):
            self.tracking_objective.deactivate()
        if hasattr(self, 'pwc_constraint'):
            self.pwc_constraint.deactivate()

        initialize_by_element_in_range(
                model,
                time,
                t0,
                ts,
                dae_vars=self.dae_vars,
                time_linking_vars=list(self.diff_vars[:]),
                outlvl=idaeslog.DEBUG,
                solver=solver,
                )

        if hasattr(self, 'tracking_objective'):
            self.tracking_objective.activate()
        if hasattr(self, 'pwc_constraint'):
            self.pwc_constraint.activate()

        if fix_inputs:
            for var in self.input_vars[:]:
                for i in range(idx_0+1, idx_s+1):
                    t = time[i]
                    var[t].unfix()

        strip_bounds.revert(model)

    def generate_inputs_at_time(self, t):
        for val in self.input_vars[:][t].value:
            yield val

    def generate_measurements_at_time(self, t):
        for var in self.measured_vars:
            yield var[t].value

    def inject_inputs(self, inputs):
        # To simulate computational delay, this function would 
        # need an argument for the start time of inputs.
        for var, val in zip(self.input_vars[:], inputs):
            # Would like:
            # self.input_vars[:,:].fix(inputs)
            # This is an example of setting a matrix from a vector.
            # Could even aspire towards:
            # self.input_vars[:,t0:t1].fix(inputs[t1])
            var[:].fix(val)

    def load_measurements(self, measured):
        t0 = self.time.first()
        # Want: self.measured_vars[:,t0].fix(measured)
        for var, val in zip(self.measured_vars, measured):
            var[t0].fix(val)

    def cycle_by(self,
            t_cycle,
            categories=(
                VariableCategory.DIFFERENTIAL,
                VariableCategory.DERIVATIVE,
                VariableCategory.ALGEBRAIC,
                VariableCategory.INPUT,
                VariableCategory.FIXED,
                ),
            tolerance=1e-8,
            ):
        time = self.time
        category_dict = self.category_dict
        for t in time:
            ts = t + t_cycle
            idx = time.find_nearest_index(ts, tolerance)
            if idx is None:
                # t + sample_time is outside the model's "horizon"
                continue
            ts = time[idx]
            for categ in categories:
                # Would like:
                # for var in component_objects(categories):
                for var in category_dict[categ]:
                    var[t].set_value(var[ts].value)

    def cycle(self,
            categories=(
                VariableCategory.DIFFERENTIAL,
                VariableCategory.DERIVATIVE,
                VariableCategory.ALGEBRAIC,
                VariableCategory.INPUT,
                VariableCategory.FIXED,
                ),
            tolerance=1e-8,
            ):
        #horizon = self.time.last() - self.time.last()
        sample_time = self.sample_time
        self.cycle_by(
                sample_time,
                categories=categories,
                tolerance=tolerance,
                )

    def generate_time_in_sample(self,
            ts,
            t0=None,
            tolerance=1e-8,
            ):
        time = self.time
        idx_s = time.find_nearest_index(ts, tolerance=tolerance)
        ts = time[idx_s]
        if t0 is None:
            t0 = ts - self.sample_time
        idx_0 = time.find_nearest_index(t0, tolerance=tolerance)
        for i in range(idx_0+1, idx_s+1):
            # Don't want to include first point in sample
            yield time[i]

    def get_data(self,
            ts,
            variables=(
                VariableCategory.DIFFERENTIAL,
                VariableCategory.INPUT,
                ),
            tolerance=1e-8,
            include_t0=False,
            ):
        time = self.time
        sample_time = self.sample_time
        category_dict = self.category_dict
        vardata_map = self.vardata_map

        data = OrderedDict()
        queue = list(variables)
        for var in queue:
            if var in VariableCategory:
                category = var
                varlist = category_dict[category]
                queue.extend(var[ts] for var in varlist)
                continue
            _slice = vardata_map[var]
            cuid = cuid_from_timeslice(_slice, time)
            if include_t0:
                i0 = time.find_nearest_index(ts-sample_time, tolerance=tolerance)
                t0 = time[i0]
                data[cuid] = [_slice[t0].value]
            else:
                data[cuid] = []
            data[cuid].extend(_slice[t].value for t in
                    self.generate_time_in_sample(ts, tolerance=tolerance))

        return data

class DynamicBlock(Block):
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
        self._init_model = ConstantInitializer(kwds.pop('model', None))
        self._init_time = ConstantInitializer(kwds.pop('time', None))
        self._init_inputs = ConstantInitializer(kwds.pop('inputs', None))
        self._init_measurements = ConstantInitializer(kwds.pop('measurements', None))
        Block.__init__(self, *args, **kwds)

    def construct(self, data=None):
        """
        Construct the DynamicBlockDatas
        """
        if self._constructed:
            return
        # Do not set the constructed flag - Block.construct() will do that

        timer = ConstructionTimer(self)

        # What is data, here?
        super(DynamicBlock, self).construct(data)

        # TODO: Should I raise an error if model, etc. was not provided?
        # Probably not, because maybe the user wants to provide them
        # manually as part of an initialization rule, rather than as
        # arguments to the component.
        self.to_dense_data()
        # When constructing an IndexedDynamicBlock, self._data is empty
        # here. Do I have to explicitly tell it to construct dense?

        # FIXME: Need to figure out a better way to initialize...
        if self._init_model.val is not None:
            if self.is_indexed():
                for index, data in iteritems(self):
                    # If indexed and model arg was provided,
                    # expect the arg to be a dict.
                    data.mod = self._init_model.val[index]
            else:
                # Else just expect it to be a model.
                self.mod = self._init_model.val
        if self._init_time.val is not None:
            if self.is_indexed():
                for index, data in iteritems(self):
                    time = self._init_time.val[index]
                    super(_BlockData, data).__setattr__('time', time)
            else:
                time = self._init_time.val
                super(_BlockData, self).__setattr__('time', time)
        if self._init_inputs.val is not None:
            if self.is_indexed():
                for index, data in iteritems(self):
                    data._inputs = self._init_inputs.val[index]
            else:
                self._inputs = self._init_inputs.val
        if self._init_measurements.val is not None:
            if self.is_indexed():
                for index, data in iteritems(self):
                    data._measurements = self._init_measurements.val[index]
            else:
                self._measurements = self._init_measurements.val

        # Will calling super().construct after adding attributes
        # break the logic in Block.construct?
        # It may not even be possible to add attributes before
        # constructing...
        #super(DynamicBlock, self).construct(data)
        for index, data in iteritems(self):
            data._construct()


class SimpleDynamicBlock(_DynamicBlockData, DynamicBlock):
    def __init__(self, *args, **kwds):
        _DynamicBlockData.__init__(self, component=self)
        DynamicBlock.__init__(self, *args, **kwds)

    # Pick up the display() from Block and not BlockData
    display = DynamicBlock.display


class IndexedDynamicBlock(DynamicBlock):

    def __init__(self, *args, **kwargs):
        DynamicBlock.__init__(self, *args, **kwargs)
