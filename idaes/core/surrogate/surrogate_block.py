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
import collections


@declare_custom_block(name="SurrogateBlock")
class SurrogateBlockData(_BlockData):
    def __init__(self, component):
        """
        This custom block is used for importing surrogates into IDAES
        
        Example usage without specifying inputs and outputs:
        > # Load the surrogate object
        > alamo_surrogate = AlamoSurrogate.load_from_file('alamo_reformer.json')
        >
        > # create the Pyomo/IDAES model
        > m = ConcreteModel()
        >
        > # create the Pyomo/IDAES block that corresponds to the surrogate
        > m.reformer = SurrogateBlock()
        > m.reformer.build_model(alamo_surrogate)
        >
        > # maximize the concentration of H2
        > m.obj = Objective(expr=m.reformer.outputs['H2'], sense=maximize)
        >
        > # constrain the outlet concentration of N2
        > m.reformer.outputs['N2'].setub(0.34)
        >
        > status = SolverFactory('ipopt').solve(m, tee=True)
        > m.reformer.display()

        Note that you can also provide variables from your IDAES model as inputs
        and outputs. For example:
        > m.reformer.build_model(alamo_surrogate, \
        >     input_vars=[m.fs.reformer_duty, m.fs.product.inlet.flow_mol],
        >     output_vars=[m.fs.total_cost])

        See the specific methods for more documentation.
        """
        super(SurrogateBlockData, self).__init__(component)

    def build_model(
        self,
        surrogate_object,
        input_vars=None,
        output_vars=None,
        use_surrogate_bounds=True,
        **kwargs
    ):
        """
        Build an EO model of the surrogate on the block. This method will build the
        necessary equations, and will construct inputs and outputs if necessary (i.e.,
        they are not provided)

        Args:
           surrogate_object: class derived from SurrogateBase (e.g., AlamoSurrogate)
              This is the surrogate object that is used to build the equations.
           input_vars : None or list of scalar or indexed variables
              If None, then the input variables are automatically created as block.inputs with keys
              coming from the input_labels of the surrogate_object. If a list is provided, that list
              can contain Var objects (indexed, scalar, or VarData) where the order must match the order
              of the inputs on the surrogate model. If input_vars are provided (not None), then
              block.inputs will not be created on the model.
           output_vars : None or list of scalar or indexed variables
              If None, then the output variables are automatically created as block.outputs with keys
              coming from the output_labels of the surrogate_object. If a list is provided, that list
              can contain Var objects (indexed, scalar, or VarData) where the order must match the order
              of the outputs on the surrogate model.  If output_vars are provided (not None), then
              block.outputs will not be created on the model.
           use_surrogate_bounds : bool
              If True, then the bounds on the input variables are updated using the bounds
              in the surrogate object. If False, then the bounds on the input variables are not overwritten.

        Any additional keyword arguments are passed to the populate_block method of the underlying surrogate object.
        Please see the documentation for the specific surrogate object for details.
        """
        self._setup_inputs_outputs(
            n_inputs=surrogate_object.n_inputs(),
            n_outputs=surrogate_object.n_outputs(),
            input_vars=input_vars,
            input_labels=surrogate_object.input_labels(),
            output_vars=output_vars,
            output_labels=surrogate_object.output_labels(),
        )

        # set the input bounds if desired
        input_bounds = surrogate_object.input_bounds()
        if use_surrogate_bounds is True and input_bounds is not None:
            input_vars_as_dict = self.input_vars_as_dict()
            for k, bnd in input_bounds.items():
                v = input_vars_as_dict[k]
                lb = bnd[0]
                if v.lb is not None:
                    lb = max(lb, v.lb)
                ub = bnd[1]
                if v.ub is not None:
                    ub = min(ub, v.ub)
                print("Setting bound of {} to {}.".format(v, (lb, ub)))
                v.setlb(lb)
                v.setub(ub)

        # call populate block to fill-in the constraints
        surrogate_object.populate_block(self, additional_options=kwargs)

        # test that kwargs is empty
        # derived classes should call .pop when they use a keyword argument
        if len(kwargs) > 0:
            raise ValueError(
                "Error in keyword arguments passed to build_model."
                " The following arguments were not used: {}".format(
                    [k for k in kwargs.keys()]
                )
            )

    def _setup_inputs_outputs(
        self,
        n_inputs,
        n_outputs,
        input_vars=None,
        input_labels=None,
        output_vars=None,
        output_labels=None,
    ):
        if n_inputs < 1 or n_outputs < 1:
            raise ValueError(
                "Attempting to create a Surrogate block with no inputs "
                "and/or no outputs. A SurrogateBlock must have at least "
                "one input and output"
            )

        # copy or create the labels
        if input_labels is None:
            self._input_labels = list(range(n_inputs))
        else:
            self._input_labels = list(input_labels)
            if len(self._input_labels) != n_inputs:
                raise ValueError(
                    "Specifying input_labels to a SurrogateBlock, but the "
                    "length does not match n_inputs"
                )

        # create or extract the variables
        if input_vars is None:
            # we need to create our own inputs
            self.inputs_set = Set(initialize=self._input_labels, ordered=True)
            self.inputs = Var(self.inputs_set, initialize=0)
            self._input_vars = list(self.inputs.values())
        else:
            # we were provided vars
            self._input_vars = _extract_var_data(input_vars)
            if len(self._input_vars) != n_inputs:
                raise ValueError(
                    "Specifying input_vars to a SurrogateBlock, but the"
                    " length of the input_vars (after extracting all"
                    " the VarData objects) does not match n_inputs"
                )

        if output_labels is None:
            self._output_labels = list(range(n_outputs))
        else:
            self._output_labels = list(output_labels)
            if len(self._output_labels) != n_outputs:
                raise ValueError(
                    "Specifying output_labels to a SurrogateBlock, but the "
                    "length does not match n_outputs"
                )

        # create or extract the output variables
        if output_vars is None:
            # we need to create our own outputs
            self.outputs_set = Set(initialize=self._output_labels, ordered=True)
            self.outputs = Var(self.outputs_set, initialize=0)
            self._output_vars = list(self.outputs.values())
        else:
            # we were provided vars
            self._output_vars = _extract_var_data(output_vars)
            if len(self._output_vars) != n_outputs:
                raise ValueError(
                    "Specifying output_vars to a SurrogateBlock, but the"
                    " length of the output_vars (after extracting all"
                    " the VarData objects) does not match n_outputs"
                )

    def _input_vars_as_list(self):
        return self._input_vars

    def _output_vars_as_list(self):
        return self._output_vars

    def input_vars_as_dict(self):
        """
        Returns a dictionary of the input variables (VarData objects) with
        the labels as the keys
        """
        return {self._input_labels[i]: v for i, v in enumerate(self._input_vars)}

    def output_vars_as_dict(self):
        """
        Returns a dictionary of the output variables (VarData objects) with
        the labels as the keys
        """
        return {self._output_labels[i]: v for i, v in enumerate(self._output_vars)}


def _extract_var_data_gen(vars):
    if vars is None:
        return
    elif getattr(vars, "ctype", None) is Var:
        if vars.is_indexed():
            if not vars.index_set().isordered():
                raise ValueError(
                    "Expected IndexedVar: {} to be indexed over an ordered set.".format(
                        vars
                    )
                )
            yield from vars.values()
        else:
            yield vars
    elif isinstance(vars, collections.abc.Sequence):
        for v in vars:
            yield from _extract_var_data_gen(v)
    else:
        raise ValueError("Unknown variable type of {} for {}".format(type(vars), vars))


def _extract_var_data(vars):
    if vars is None:
        return None
    return list(_extract_var_data_gen(vars))
