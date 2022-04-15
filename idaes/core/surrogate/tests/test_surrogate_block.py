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
"""
Tests for SurrogateBlock
"""
import pytest

from pyomo.environ import ConcreteModel, Set, Var

from idaes.core.surrogate.surrogate_block import SurrogateBlock, _extract_var_data


@pytest.mark.unit
def test_extract_var_data_None():
    assert _extract_var_data(None) is None


@pytest.mark.unit
def test_extract_var_data_scalar():
    m = ConcreteModel()
    m.v = Var()
    assert _extract_var_data(m.v) == [m.v]


@pytest.mark.unit
def test_extract_var_data_indexed():
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3, 4])
    m.v = Var(m.s)
    assert _extract_var_data(m.v) == [m.v[1], m.v[2], m.v[3], m.v[4]]


@pytest.mark.unit
def test_extract_var_data_indexed_not_ordered():
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3, 4], ordered=False)
    m.v = Var(m.s)
    with pytest.raises(
        ValueError, match="Expected IndexedVar: v to be indexed over " "an ordered set."
    ):
        _extract_var_data(m.v)


@pytest.mark.unit
def test_extract_var_data_multiple():
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3, 4])
    m.vs = Var()
    m.vi = Var(m.s)
    assert _extract_var_data([m.vs, m.vi, m.vs]) == [
        m.vs,
        m.vi[1],
        m.vi[2],
        m.vi[3],
        m.vi[4],
        m.vs,
    ]


@pytest.mark.unit
def test_extract_var_data_not_var():
    m = ConcreteModel()
    with pytest.raises(
        ValueError,
        match="Unknown variable type of <class "
        "'pyomo.core.base.PyomoModel.ConcreteModel'> for unknown",
    ):
        _extract_var_data(m)


@pytest.mark.unit
def test_setup_inputs_outputs_zero_inputs():
    m = ConcreteModel()
    m.sb = SurrogateBlock()
    with pytest.raises(
        ValueError,
        match="Attempting to create a Surrogate block with no inputs "
        "and/or no outputs. A SurrogateBlock must have at least "
        "one input and output",
    ):
        m.sb._setup_inputs_outputs(0, 1)


@pytest.mark.unit
def test_setup_inputs_outputs_zero_outputs():
    m = ConcreteModel()
    m.sb = SurrogateBlock()
    with pytest.raises(
        ValueError,
        match="Attempting to create a Surrogate block with no inputs "
        "and/or no outputs. A SurrogateBlock must have at least "
        "one input and output",
    ):
        m.sb._setup_inputs_outputs(1, 0)


@pytest.mark.unit
def test_setup_inputs_outputs_no_vars_or_labels():
    m = ConcreteModel()
    m.sb = SurrogateBlock()

    m.sb._setup_inputs_outputs(4, 5)

    assert m.sb._input_labels == [0, 1, 2, 3]
    assert isinstance(m.sb.inputs_set, Set)
    assert isinstance(m.sb.inputs, Var)
    for i in m.sb.inputs_set:
        assert i in m.sb._input_labels
        assert m.sb.inputs[i].value == 0
    assert m.sb._input_vars == [
        m.sb.inputs[0],
        m.sb.inputs[1],
        m.sb.inputs[2],
        m.sb.inputs[3],
    ]

    assert m.sb._output_labels == [0, 1, 2, 3, 4]
    assert isinstance(m.sb.outputs_set, Set)
    assert isinstance(m.sb.outputs, Var)
    for i in m.sb.outputs_set:
        assert i in m.sb._output_labels
        assert m.sb.outputs[i].value == 0
    assert m.sb._output_vars == [
        m.sb.outputs[0],
        m.sb.outputs[1],
        m.sb.outputs[2],
        m.sb.outputs[3],
        m.sb.outputs[4],
    ]


@pytest.mark.unit
def test_setup_inputs_outputs_input_label_mismatch():
    m = ConcreteModel()
    m.sb = SurrogateBlock()

    with pytest.raises(
        ValueError,
        match="Specifying input_labels to a SurrogateBlock, but the "
        "length does not match n_inputs",
    ):
        m.sb._setup_inputs_outputs(4, 5, input_labels=["x1"])


@pytest.mark.unit
def test_setup_inputs_outputs_output_label_mismatch():
    m = ConcreteModel()
    m.sb = SurrogateBlock()

    with pytest.raises(
        ValueError,
        match="Specifying output_labels to a SurrogateBlock, but the "
        "length does not match n_outputs",
    ):
        m.sb._setup_inputs_outputs(4, 5, output_labels=["x1"])


@pytest.mark.unit
def test_setup_inputs_outputs_w_labels_no_vars():
    m = ConcreteModel()
    m.sb = SurrogateBlock()

    m.sb._setup_inputs_outputs(
        4,
        5,
        input_labels=["x1", "x2", "x3", "x4"],
        output_labels=["z1", "z2", "z3", "z4", "z5"],
    )

    assert m.sb._input_labels == ["x1", "x2", "x3", "x4"]
    assert isinstance(m.sb.inputs_set, Set)
    assert isinstance(m.sb.inputs, Var)
    for i in m.sb.inputs_set:
        assert i in m.sb._input_labels
        assert m.sb.inputs[i].value == 0
    assert m.sb._input_vars == [
        m.sb.inputs["x1"],
        m.sb.inputs["x2"],
        m.sb.inputs["x3"],
        m.sb.inputs["x4"],
    ]

    assert m.sb._output_labels == ["z1", "z2", "z3", "z4", "z5"]
    assert isinstance(m.sb.outputs_set, Set)
    assert isinstance(m.sb.outputs, Var)
    for i in m.sb.outputs_set:
        assert i in m.sb._output_labels
        assert m.sb.outputs[i].value == 0
    assert m.sb._output_vars == [
        m.sb.outputs["z1"],
        m.sb.outputs["z2"],
        m.sb.outputs["z3"],
        m.sb.outputs["z4"],
        m.sb.outputs["z5"],
    ]


@pytest.fixture(scope="module")
def surr_block():
    m = ConcreteModel()
    m.input_set = Set(initialize=["x1", "x2", "x3", "x4"])
    m.output_set = ["z1", "z2", "z3", "z4", "z5"]
    m.inputs = Var(m.input_set)
    m.outputs = Var(m.output_set)

    m.sb = SurrogateBlock()

    m.sb._setup_inputs_outputs(4, 5, input_vars=m.inputs, output_vars=m.outputs)

    return m


@pytest.mark.unit
def test_setup_inputs_outputs_w_vars_no_labels(surr_block):
    assert surr_block.sb._input_labels == [0, 1, 2, 3]
    assert not hasattr(surr_block.sb, "inputs_set")
    assert not hasattr(surr_block.sb, "inputs")
    assert surr_block.sb._input_vars == [
        surr_block.inputs["x1"],
        surr_block.inputs["x2"],
        surr_block.inputs["x3"],
        surr_block.inputs["x4"],
    ]

    assert surr_block.sb._output_labels == [0, 1, 2, 3, 4]
    assert not hasattr(surr_block.sb, "outputs_set")
    assert not hasattr(surr_block.sb, "outputs")
    assert surr_block.sb._output_vars == [
        surr_block.outputs["z1"],
        surr_block.outputs["z2"],
        surr_block.outputs["z3"],
        surr_block.outputs["z4"],
        surr_block.outputs["z5"],
    ]


@pytest.mark.unit
def test_setup_inputs_outputs_input_vars_mismatch():
    m = ConcreteModel()
    m.input_set = Set(initialize=["x1", "x2", "x3"])
    m.output_set = ["z1", "z2", "z3", "z4", "z5"]
    m.inputs = Var(m.input_set)
    m.outputs = Var(m.output_set)

    m.sb = SurrogateBlock()

    with pytest.raises(
        ValueError,
        match="Specifying input_vars to a SurrogateBlock, but "
        "the length of the input_vars \(after extracting all"
        " the VarData objects\) does not match n_inputs",
    ):
        m.sb._setup_inputs_outputs(4, 5, input_vars=m.inputs, output_vars=m.outputs)


@pytest.mark.unit
def test_setup_inputs_outputs_input_vars_empty_list():
    m = ConcreteModel()
    m.output_set = ["z1", "z2", "z3", "z4", "z5"]
    m.outputs = Var(m.output_set)

    m.sb = SurrogateBlock()

    with pytest.raises(
        ValueError,
        match="Specifying input_vars to a SurrogateBlock, but "
        "the length of the input_vars \(after extracting all"
        " the VarData objects\) does not match n_inputs",
    ):
        m.sb._setup_inputs_outputs(4, 5, input_vars=[], output_vars=m.outputs)


@pytest.mark.unit
def test_setup_inputs_outputs_output_vars_mismatch():
    m = ConcreteModel()
    m.input_set = Set(initialize=["x1", "x2", "x3", "x4"])
    m.output_set = ["z1", "z2", "z3", "z4"]
    m.inputs = Var(m.input_set)
    m.outputs = Var(m.output_set)

    m.sb = SurrogateBlock()

    with pytest.raises(
        ValueError,
        match="Specifying output_vars to a SurrogateBlock, but "
        "the length of the output_vars \(after extracting all"
        " the VarData objects\) does not match n_outputs",
    ):
        m.sb._setup_inputs_outputs(4, 5, input_vars=m.inputs, output_vars=m.outputs)


@pytest.mark.unit
def test_setup_inputs_outputs_output_vars_empty_list():
    m = ConcreteModel()
    m.input_set = Set(initialize=["x1", "x2", "x3", "x4"])
    m.inputs = Var(m.input_set)

    m.sb = SurrogateBlock()

    with pytest.raises(
        ValueError,
        match="Specifying output_vars to a SurrogateBlock, but "
        "the length of the output_vars \(after extracting all"
        " the VarData objects\) does not match n_outputs",
    ):
        m.sb._setup_inputs_outputs(4, 5, input_vars=m.inputs, output_vars=[])


@pytest.mark.unit
def test_input_vars_as_list(surr_block):
    assert surr_block.sb._input_vars_as_list() == [
        surr_block.inputs["x1"],
        surr_block.inputs["x2"],
        surr_block.inputs["x3"],
        surr_block.inputs["x4"],
    ]


@pytest.mark.unit
def test_output_vars_as_list(surr_block):
    assert surr_block.sb._output_vars_as_list() == [
        surr_block.outputs["z1"],
        surr_block.outputs["z2"],
        surr_block.outputs["z3"],
        surr_block.outputs["z4"],
        surr_block.outputs["z5"],
    ]


@pytest.mark.unit
def test_input_vars_as_dict(surr_block):
    assert surr_block.sb.input_vars_as_dict() == {
        0: surr_block.inputs["x1"],
        1: surr_block.inputs["x2"],
        2: surr_block.inputs["x3"],
        3: surr_block.inputs["x4"],
    }


@pytest.mark.unit
def test_output_vars_as_dict(surr_block):
    assert surr_block.sb.output_vars_as_dict() == {
        0: surr_block.outputs["z1"],
        1: surr_block.outputs["z2"],
        2: surr_block.outputs["z3"],
        3: surr_block.outputs["z4"],
        4: surr_block.outputs["z5"],
    }


class DummySurrogate:
    # Dummy a SurrogateObject to use to test build method
    def populate_block(self, blk, **kwargs):
        blk._populate_called = True

    def n_inputs(self):
        return 4

    def n_outputs(self):
        return 5

    def input_labels(self):
        return ["x1", "x2", "x3", "x4"]

    def output_labels(self):
        return ["z1", "z2", "z3", "z4", "z5"]

    def input_bounds(self):
        return {"x1": (1, 10), "x2": (2, 20), "x3": (3, 30), "x4": (4, 40)}


@pytest.mark.unit
def test_build_model():
    m = ConcreteModel()

    m.dummy = DummySurrogate()

    m.sb = SurrogateBlock()

    m.sb.build_model(m.dummy, use_surrogate_bounds=False)

    # Check that the _populate_called attribute was added to the SurrogateBlock
    assert m.sb._populate_called

    # Check variables and bounds
    assert m.sb._input_labels == ["x1", "x2", "x3", "x4"]
    assert isinstance(m.sb.inputs_set, Set)
    assert isinstance(m.sb.inputs, Var)
    for i in m.sb.inputs_set:
        assert i in m.sb._input_labels
        assert m.sb.inputs[i].value == 0
        assert m.sb.inputs[i]._lb is None
        assert m.sb.inputs[i]._ub is None
    assert m.sb._input_vars == [
        m.sb.inputs["x1"],
        m.sb.inputs["x2"],
        m.sb.inputs["x3"],
        m.sb.inputs["x4"],
    ]

    assert m.sb._output_labels == ["z1", "z2", "z3", "z4", "z5"]
    assert isinstance(m.sb.outputs_set, Set)
    assert isinstance(m.sb.outputs, Var)
    for i in m.sb.outputs_set:
        assert i in m.sb._output_labels
        assert m.sb.outputs[i].value == 0
    assert m.sb._output_vars == [
        m.sb.outputs["z1"],
        m.sb.outputs["z2"],
        m.sb.outputs["z3"],
        m.sb.outputs["z4"],
        m.sb.outputs["z5"],
    ]


@pytest.mark.unit
def test_build_model_w_bounds():
    m = ConcreteModel()

    m.dummy = DummySurrogate()

    m.sb = SurrogateBlock()

    m.sb.build_model(m.dummy, use_surrogate_bounds=True)

    # Check that the _populate_called attribute was added to the SurrogateBlock
    assert m.sb._populate_called

    # Check variables and bounds
    assert m.sb._input_labels == ["x1", "x2", "x3", "x4"]
    assert isinstance(m.sb.inputs_set, Set)
    assert isinstance(m.sb.inputs, Var)
    for i in m.sb.inputs_set:
        assert i in m.sb._input_labels
        assert m.sb.inputs[i].value == 0
        assert m.sb.inputs[i]._lb == m.dummy.input_bounds()[i][0]
        assert m.sb.inputs[i]._ub == m.dummy.input_bounds()[i][1]
    assert m.sb._input_vars == [
        m.sb.inputs["x1"],
        m.sb.inputs["x2"],
        m.sb.inputs["x3"],
        m.sb.inputs["x4"],
    ]

    assert m.sb._output_labels == ["z1", "z2", "z3", "z4", "z5"]
    assert isinstance(m.sb.outputs_set, Set)
    assert isinstance(m.sb.outputs, Var)
    for i in m.sb.outputs_set:
        assert i in m.sb._output_labels
        assert m.sb.outputs[i].value == 0
    assert m.sb._output_vars == [
        m.sb.outputs["z1"],
        m.sb.outputs["z2"],
        m.sb.outputs["z3"],
        m.sb.outputs["z4"],
        m.sb.outputs["z5"],
    ]


@pytest.mark.unit
def test_build_model_unused_kwarg():
    m = ConcreteModel()

    m.dummy = DummySurrogate()

    m.sb = SurrogateBlock()

    with pytest.raises(
        ValueError,
        match="Error in keyword arguments passed to "
        "build_model. The following arguments were not used: "
        "\['foo'\]",
    ):
        m.sb.build_model(m.dummy, use_surrogate_bounds=False, foo=True)
