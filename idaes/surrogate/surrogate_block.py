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
from pyomo.core.base.block import declare_custom_block, _BlockData
from pyomo.environ import Var, Set
from pyomo.core.base.var import ScalarVar, IndexedVar

@declare_custom_block(name='SurrogateBlock')
class SurrogateBlockData(_BlockData):
    def __init__(self, component):
        super(SurrogateBlockData,self).__init__(component)

    def build_model(self, surrogate_object, input_vars=None, output_vars=None,
                    use_surrogate_bounds=True, **kwargs):
        self._setup_inputs_outputs(
            n_inputs=surrogate_object.n_inputs,
            n_outputs=surrogate_object.n_outputs,
            input_vars=input_vars, input_labels=surrogate_object.input_labels,
            output_vars=output_vars, output_labels=surrogate_object.output_labels)

        # set the input bounds if desired
        input_bounds = surrogate_object.input_bounds
        if use_surrogate_bounds is True and input_bounds is not None:
            input_vars_as_dict = self._input_vars_as_dict()
            for k,bnd in input_bounds.items():
                v = input_vars_as_dict[k]
                lb = bnd[0]
                if v.lb is not None:
                    lb = max(lb, v.lb)
                ub = bnd[1]
                if v.ub is not None:
                    ub = min(ub, v.ub)
                v.setlb(lb)
                v.setub(ub)

        # call populate block to fill-in the constraints
        surrogate_object.populate_block(self, **kwargs)

    def _setup_inputs_outputs(self, n_inputs, n_outputs,
                              input_vars=None, input_labels=None,
                              output_vars=None, output_labels=None):
        if n_inputs < 1 or n_outputs < 1:
            raise ValueError('Attempting to create a Surrogate block with no inputs '
                             'and/or no outputs. A SurrogateBlock must have at least '
                             'one input and output')

        # copy or create the labels
        if input_labels is None:
            self._input_labels = list(range(n_inputs))
        else:
            self._input_labels = list(input_labels)
            if len(self._input_labels) != n_inputs:
                raise ValueError('Specifying input_labels to a SurrogateBlock, but the length'
                                 ' does not match n_inputs')

        # create or extract the variables
        if input_vars is None:
            # we need to create our own inputs
            self._inputs_set = Set(initialize=self._input_labels, ordered=True)
            self._inputs = Var(self._inputs_set, initialize=0)
            self._input_vars = list(self._inputs.values())
        else:
            # we were provided vars
            self._input_vars = _extract_var_data(input_vars)
            if len(self._input_vars) != n_inputs:
                raise ValueError('Specifying input_vars to a SurrogateBlock, but the'
                                 ' length of the input_vars (after extracting all'
                                 ' the VarData objects) does not match n_inputs')

        if output_labels is None:
            self._output_labels = list(range(n_outputs))
        else:
            self._output_labels = list(output_labels)
            if len(self._output_labels) != n_outputs:
                raise ValueError('Specifying output_labels to a SurrogateBlock, but the length'
                                 ' does not match n_outputs')

        # create or extract the output variables
        if output_vars is None:
            # we need to create our own outputs
            self._outputs_set = Set(initialize=self._output_labels, ordered=True)
            self._outputs = Var(self._outputs_set, initialize=0)
            self._output_vars = list(self._outputs.values())
        else:
            # we were provided vars
            self._output_vars = _extract_var_data(output_vars)
            if len(self._output_vars) != n_outputs:
                raise ValueError('Specifying output_vars to a SurrogateBlock, but the'
                                 ' length of the output_vars (after extracting all'
                                 ' the VarData objects) does not match n_outputs')

    def _input_vars_as_list(self):
        return self._input_vars

    def _output_vars_as_list(self):
        return self._output_vars

    def _input_vars_as_dict(self):
        return {self._input_labels[i]:v for i,v in enumerate(self._input_vars)}

    def _output_vars_as_dict(self):
        return {self._output_labels[i]:v for i,v in enumerate(self._output_vars)}


def _extract_var_data(vars):
    if vars is None:
        return None
    elif isinstance(vars, ScalarVar):
        return [vars]
    elif isinstance(vars, IndexedVar):
        if vars.indexed_set().is_ordered():
            return list(vars.values())
        raise ValueError('Expected IndexedVar: {} to be indexed over an ordered set.'.format(vars))
    elif isinstance(vars, list):
        varlist = list()
        for v in vars:
            if v.is_indexed():
                varlist.extend(v.values())
            else:
                varlist.append(v)
        return varlist
    else:
        raise ValueError("Unknown variable type {}".format(vars))
