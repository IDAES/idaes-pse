# -*- coding: utf-8 -*-
##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Base class for dynamic simulation objects.
"""
from pyutilib.misc.config import ConfigDict, ConfigValue
from pyomo.core.base.block import Block
from pyomo.core.base.var import Var
from pyomo.dae import DerivativeVar
from pyomo.dae.flatten import flatten_dae_components
from pyomo.common.collections import ComponentSet, ComponentMap
from idaes.apps.caprese.common import config as dyn_config
from idaes.apps.caprese.common.config import (VariableCategory)
from idaes.apps.caprese.util import NMPCVarGroup, NMPCVarLocator
import idaes.core.util.dyn_utils as dyn_utils
import  idaes.logger as idaeslog

__author__ = "Robert Parker"

class DynamicBase(object):
    """
    Base class for dynamic simulations objects.
    """

    CONFIG = ConfigDict()

    CONFIG.declare(
            'outlvl',
            ConfigValue(
                default=idaeslog.INFO,
                doc='Severity threshold for IDAES logger',
                )
            )
    # Would a time_resolution_option be useful for MHE?
    CONFIG.declare(
            'continuous_set_tolerance',
            ConfigValue(
                default=1e-8,
                domain=float,
                doc=('Tolerance used for determining whether a float is a '
                    'member of a ContinuousSet'),
                )
            )
    CONFIG.declare(
            'solver',
            ConfigValue(
                default='ipopt',
                domain=dyn_config.validate_solver,
                doc='Pyomo solver object to be used to solve generated NLPs',
                )
            )
    CONFIG.declare(
            'tolerance',
            ConfigValue(
                default=1e-8,
                domain=float,
                doc='Tolerance for checking constraint violation',
                )
            )

    @classmethod
    def get_namespace_name(cls):
        return '_DYNAMIC_NAMESPACE'

    @classmethod
    def add_namespace_to(cls, model, time):
        """
        """
        name = DynamicBase.get_namespace_name()
        derived_name = cls.namespace_name
        if hasattr(model, name):
            # Return if namespace has already been added. Don't throw an error
            # as this is expected if the user, say wants to use the same model
            # for NMPC and MHE.
            return
        if time.model() != model.model():
            raise ValueError(
                'time must belong to same top-level model as model')
        model.add_component(name, Block())
        namespace = getattr(model, name)
        derived_namespace = getattr(model, derived_name)

        def get_time():
            return time
        namespace.get_time = get_time
        derived_namespace.get_time = namespace.get_time

        # Validate discretization scheme and get ncp:
        namespace.ncp = dyn_config.get_ncp(time)

        namespace.variables_categorized = False

    @classmethod
    def remove_namespace_from(cls, model):
        """
        """
        # TODO: add remove_namespace_from to derived classes
        name = DynamicBase.get_namespace_name()
        if not hasattr(model, name):
            raise RuntimeError(
                'Trying to delete block %s that does not exist on model'
                % name)
        model.del_component(name, Block())

    @classmethod
    def get_logger_name(cls):
        return 'dynamic'

    def _populate_namespace(self, model):
        """
        Given a model with categorized variables, a category_dict, and a
        var_locator, each referenced through the DynamicBase namespace,
        adds references to each of these objects through model's namespace
        corresponding to this particular instance of DynamicBase
        (i.e. NMPCSim or MHESim).
        """
        if self.get_namespace_name() == DynamicBase.get_namespace_name():
            return
        derived_namespace = getattr(model, self.get_namespace_name())
        base_namespace = getattr(model, DynamicBase.get_namespace_name())

        derived_namespace.dae_vars = base_namespace.dae_vars
        derived_namespace.diff_vars = base_namespace.diff_vars
        derived_namespace.n_diff_vars = base_namespace.n_diff_vars
        derived_namespace.deriv_vars = base_namespace.deriv_vars
        derived_namespace.n_deriv_vars = base_namespace.n_deriv_vars
        derived_namespace.input_vars = base_namespace.input_vars
        derived_namespace.n_input_vars = base_namespace.n_input_vars
        derived_namespace.alg_vars = base_namespace.alg_vars
        derived_namespace.n_alg_vars = base_namespace.n_alg_vars
        derived_namespace.fixed_vars = base_namespace.fixed_vars
        derived_namespace.n_fixed_vars = base_namespace.n_fixed_vars
        derived_namespace.scalar_vars = base_namespace.scalar_vars
        derived_namespace.n_scalar_vars = base_namespace.n_scalar_vars

        derived_namespace.ic_vars = base_namespace.ic_vars
        derived_namespace.n_ic_vars = base_namespace.n_ic_vars

        derived_namespace.variables_categorized = \
                base_namespace.variables_categorized

        derived_namespace.get_time = base_namespace.get_time
        derived_namespace.ncp = base_namespace.ncp

        derived_namespace.category_dict = base_namespace.category_dict
        derived_namespace.var_locator = base_namespace.var_locator


    def __init__(self, plant, plant_time, controller, controller_time,
            inputs_at_t0, **kwargs):

        self.config = self.CONFIG(kwargs)
        dyn_config.validate_list_of_vardata(inputs_at_t0)

        # How, if at all, should I configure logger here?

        # add namespace
        self.add_namespace_to(plant, plant_time)
        self.add_namespace_to(controller, controller_time)
        self.plant = plant
        self.plant_time = plant_time
        self.controller = controller
        self.controller_time = controller_time

        self.sample_time = None

        # categorize, create category dicts, create locator
        # ^should this be done separately for NMPC/MHE?
        self.categorize_variables(self.plant, inputs_at_t0)
        namespace = getattr(self.plant, DynamicBase.get_namespace_name())

        namespace.category_dict = {
                VariableCategory.DIFFERENTIAL: namespace.diff_vars,
                VariableCategory.DERIVATIVE: namespace.deriv_vars,
                VariableCategory.ALGEBRAIC: namespace.alg_vars,
                VariableCategory.INPUT: namespace.input_vars,
                VariableCategory.FIXED: namespace.fixed_vars,
                VariableCategory.SCALAR: namespace.scalar_vars,
                }
        self.build_variable_locator(self.plant,
                namespace.category_dict,
                ic_vars=namespace.ic_vars)
#                measurement_vars=namespace.measurement_vars)

        # find input and measurement vars in controller model
        init_controller_inputs = [
                dyn_utils.find_comp_in_block(self.controller, self.plant, comp)
                for comp in inputs_at_t0]

        self.categorize_variables(self.controller, init_controller_inputs)
        namespace = getattr(self.controller,
                DynamicBase.get_namespace_name())

        namespace.category_dict = {
                VariableCategory.DIFFERENTIAL: namespace.diff_vars,
                VariableCategory.DERIVATIVE: namespace.deriv_vars,
                VariableCategory.ALGEBRAIC: namespace.alg_vars,
                VariableCategory.INPUT: namespace.input_vars,
                VariableCategory.FIXED: namespace.fixed_vars,
                VariableCategory.SCALAR: namespace.scalar_vars,
                }
        self.build_variable_locator(self.controller,
                namespace.category_dict,
                ic_vars=namespace.ic_vars)

        # Populate derived-class namespaces with attributes just constructed
        # on the base class namespace
        self._populate_namespace(self.plant)
        self._populate_namespace(self.controller)


    def set_sample_time(self, sample_time):
        """
        Validates sample time and adds as attribute to model.
        This method exists because providing sample time should not be
        required in constructor, and could change during simulation.
        """
        self.validate_sample_time(sample_time)
        self.sample_time = sample_time

    def validate_sample_time(self, sample_time):
        raise NotImplementedError(
            'Derived class must implement method for validating sample time')

    @staticmethod
    def categorize_variables(model, initial_inputs):
        """Creates lists of time-only-slices of the different types of variables
        in a model, given knowledge of which are inputs. These lists are added 
        as attributes to the model's namespace.

        Possible variable categories are:

            - INPUT --- Those specified by the user to be inputs
            - DERIVATIVE --- Those declared as Pyomo DerivativeVars, whose 
                             "state variable" is not fixed, except possibly as an
                             initial condition
            - DIFFERENTIAL --- Those referenced as the "state variable" by an
                               unfixed (except possibly as an initial condition)
                               DerivativeVar
            - FIXED --- Those that are fixed at non-initial time points. These
                        are typically disturbances, design variables, or 
                        uncertain parameters.
            - ALGEBRAIC --- Unfixed, time-indexed variables that are neither
                            inputs nor referenced by an unfixed derivative.
            - SCALAR --- Variables unindexed by time. These could be variables
                         that refer to a specific point in time (initial or
                         final conditions), averages over time, or truly 
                         time-independent variables like diameter.

        Args:
            model : Model whose variables will be flattened and categorized
            initial_inputs : List of VarData objects that are input variables
                             at the initial time point

        """
        namespace = getattr(model,
                DynamicBase.get_namespace_name())
        time = namespace.get_time()
        t0 = time.first()
        t1 = time.get_finite_elements()[1]
        deriv_vars = []
        diff_vars = []
        input_vars = []
        alg_vars = []
        fixed_vars = []

        ic_vars = []

        # Create list of time-only-slices of time indexed variables
        # (And list of VarData objects for scalar variables)
        scalar_vars, dae_vars = flatten_dae_components(model, time, Var)

        dae_map = ComponentMap([(v[t0], v) for v in dae_vars])
        t0_vardata = list(dae_map.keys())
        namespace.dae_vars = list(dae_map.values())
        namespace.scalar_vars = \
            NMPCVarGroup(
                list(ComponentMap([(v, v) for v in scalar_vars]).values()),
                index_set=None, is_scalar=True)
        namespace.n_scalar_vars = \
                namespace.scalar_vars.n_vars
        input_set = ComponentSet(initial_inputs)
        updated_input_set = ComponentSet(initial_inputs)

        # Iterate over initial vardata, popping from dae map when an input,
        # derivative, or differential var is found.
        for var0 in t0_vardata:
            if var0 in updated_input_set:
                input_set.remove(var0)
                time_slice = dae_map.pop(var0)
                input_vars.append(time_slice)

            parent = var0.parent_component()
            if not isinstance(parent, DerivativeVar):
                continue
            if not time in ComponentSet(parent.get_continuousset_list()):
                continue
            index0 = var0.index()
            var1 = dae_map[var0][t1]
            index1 = var1.index()
            state = parent.get_state_var()

            if state[index1].fixed:
                # Assume state var is fixed everywhere, so derivative
                # 'isn't really' a derivative.
                # Should be safe to remove state from dae_map here
                state_slice = dae_map.pop(state[index0])
                fixed_vars.append(state_slice)
                continue
            if state[index0] in input_set:
                # If differential variable is an input, then this DerivativeVar
                # is 'not really a derivative'
                continue

            deriv_slice = dae_map.pop(var0)

            if var1.fixed:
                # Assume derivative has been fixed everywhere.
                # Add to list of fixed variables, and don't remove its state variable.
                fixed_vars.append(deriv_slice)
            elif var0.fixed:
                # In this case the derivative has been used as an initial condition. 
                # Still want to include it in the list of derivatives.
                ic_vars.append(deriv_slice)
                state_slice = dae_map.pop(state[index0])
                if state[index0].fixed:
                    ic_vars.append(state_slice)
                deriv_vars.append(deriv_slice)
                diff_vars.append(state_slice)
            else:
                # Neither is fixed. This should be the most common case.
                state_slice = dae_map.pop(state[index0])
                if state[index0].fixed:
                    ic_vars.append(state_slice)
                deriv_vars.append(deriv_slice)
                diff_vars.append(state_slice)

        if not updated_input_set:
            raise RuntimeError('Not all inputs could be found')
        assert len(deriv_vars) == len(diff_vars)

        for var0, time_slice in dae_map.items():
            var1 = time_slice[t1]
            # If the variable is still in the list of time-indexed vars,
            # it must either be fixed (not a var) or be an algebraic var
            if var1.fixed:
                fixed_vars.append(time_slice)
            else:
                if var0.fixed:
                    ic_vars.append(time_slice)
                alg_vars.append(time_slice)

        namespace.deriv_vars = NMPCVarGroup(deriv_vars, time)
        namespace.diff_vars = NMPCVarGroup(diff_vars, time)
        namespace.n_diff_vars = len(diff_vars)
        namespace.n_deriv_vars = len(deriv_vars)
        assert (namespace.n_diff_vars == 
                namespace.n_deriv_vars)
                
        # ic_vars will not be stored as a NMPCVarGroup - don't want to store
        # all the info twice
        namespace.ic_vars = ic_vars
        namespace.n_ic_vars = len(ic_vars)
        #assert model.n_dv == len(ic_vars)
        # Would like this to be true, but accurately detecting differential
        # variables that are not implicitly fixed (by fixing some input)
        # is difficult
        # Also, a categorization can have no input vars and still be
        # valid for MHE

        namespace.input_vars = NMPCVarGroup(input_vars, time)
        namespace.n_input_vars = len(input_vars)

        namespace.alg_vars = NMPCVarGroup(alg_vars, time)
        namespace.n_alg_vars = len(alg_vars)

        namespace.fixed_vars = NMPCVarGroup(fixed_vars, time)
        namespace.n_fixed_vars = len(fixed_vars)

        namespace.variables_categorized = True

    @staticmethod
    def build_variable_locator(model, category_dict, ic_vars=[],
            measurement_vars=[]):
        """Constructs a ComponentMap mapping each VarData object
        to a NMPCVarLocator object. This dictionary is added as an attribute to
        the model's namespace.

        Args:
            model : Flowsheet model containing the variables provided
            category_dict : Dictionary mapping VariableCategory enum items
                            to NMPCVarGroup instances
            ic_vars : List of variables (time-only slices) that are fixed
                      only at the initial time point

        """
        namespace = getattr(model,
                DynamicBase.get_namespace_name())
        time = namespace.get_time()
        ic_list = ic_vars

        locator = ComponentMap()
        for categ, vargroup in category_dict.items():
            varlist = vargroup.varlist
            if categ == VariableCategory.SCALAR:
                for i, var in enumerate(varlist):
                    locator[var] = NMPCVarLocator(categ, vargroup, i)
            else:
                for i, var in enumerate(varlist):
                    for t in time:
                        locator[var[t]] = NMPCVarLocator(categ, vargroup, i)

        # Since these variables already have NMPCVarLocator objects,
        # just set the desired attribute.
        for i, _slice in enumerate(ic_list):
            for t in time:
                locator[_slice[t]].is_ic = True

        for i, _slice in enumerate(measurement_vars):
            for t in time:
                locator[_slice[t]].is_measurement = True

        namespace.var_locator = locator


