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

from pandas import isna
import pytest

from pyomo.environ import (
    ConcreteModel,
    Expression,
    TransformationFactory,
    Var,
    value
)
from pyomo.network import Arc
from idaes.core import (
    FlowsheetBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class
)
from idaes.core.util.tables import (
    arcs_to_stream_dict,
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
    generate_table,
    tag_state_quantities,
    stream_states_dict,
)
import idaes.generic_models.properties.examples.saponification_thermo as thermo_props
import idaes.generic_models.properties.examples.saponification_reactions as rxn_props
from idaes.generic_models.unit_models import CSTR


@pytest.fixture()
def m():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.thermo_params = thermo_props.SaponificationParameterBlock()
    m.fs.reaction_params = rxn_props.SaponificationReactionParameterBlock(
        default={"property_package": m.fs.thermo_params})

    m.fs.tank1 = CSTR(default={"property_package": m.fs.thermo_params,
                               "reaction_package": m.fs.reaction_params})
    m.fs.tank2 = CSTR(default={"property_package": m.fs.thermo_params,
                               "reaction_package": m.fs.reaction_params})
    m.fs.tank_array = CSTR(
        range(3),
        default={"property_package": m.fs.thermo_params,
                 "reaction_package": m.fs.reaction_params})

    m.fs.stream = Arc(source=m.fs.tank1.outlet,
                      destination=m.fs.tank2.inlet)

    def stream_array_rule(b, i):
        return (b.tank_array[i].outlet, b.tank_array[i+1].inlet)
    m.fs.stream_array = Arc(range(2), rule=stream_array_rule)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


@pytest.mark.unit
def test_create_stream_table_dataframe_from_StateBlock(m):
    d = arcs_to_stream_dict(m, descend_into=True, prepend="model")
    assert "model.stream" in d
    assert "model.stream_array" in d
    assert d["model.stream"] == m.fs.stream
    assert d["model.stream_array"] == m.fs.stream_array

@pytest.mark.unit
def test_stream_states_dict(m):
    d = stream_states_dict(arcs_to_stream_dict(m, descend_into=True))
    assert "stream" in d
    assert "stream_array[0]" in d
    assert "stream_array[1]" in d
    assert d["stream_array[0]"] == m.fs.tank_array[1].control_volume.properties_in[0]
    assert d["stream_array[1]"] == m.fs.tank_array[2].control_volume.properties_in[0]

@pytest.mark.unit
def test_create_stream_table_dataframe_from_StateBlock_2(m):
    df = create_stream_table_dataframe(arcs_to_stream_dict(m, descend_into=True))

    assert df.loc["Pressure"]["stream"] == 101325
    assert df.loc["Temperature"]["stream"] == 298.15
    assert df.loc["Volumetric Flowrate"]["stream"] == 1.0
    assert df.loc["Molar Concentration H2O"]["stream"] == 100.0
    assert df.loc["Molar Concentration NaOH"]["stream"] == 100.0
    assert df.loc["Molar Concentration EthylAcetate"]["stream"] == 100.0
    assert df.loc["Molar Concentration SodiumAcetate"]["stream"] == 100.0
    assert df.loc["Molar Concentration Ethanol"]["stream"] == 100.0

    assert df.loc["Pressure"]["stream_array[0]"] == 101325
    assert df.loc["Temperature"]["stream_array[0]"] == 298.15
    assert df.loc["Volumetric Flowrate"]["stream_array[0]"] == 1.0
    assert df.loc["Molar Concentration H2O"]["stream_array[0]"] == 100.0
    assert df.loc["Molar Concentration NaOH"]["stream_array[0]"] == 100.0
    assert df.loc["Molar Concentration EthylAcetate"]["stream_array[0]"] == 100.0
    assert df.loc["Molar Concentration SodiumAcetate"]["stream_array[0]"] == 100.0
    assert df.loc["Molar Concentration Ethanol"]["stream_array[0]"] == 100.0

    assert df.loc["Pressure"]["stream_array[1]"] == 101325
    assert df.loc["Temperature"]["stream_array[1]"] == 298.15
    assert df.loc["Volumetric Flowrate"]["stream_array[1]"] == 1.0
    assert df.loc["Molar Concentration H2O"]["stream_array[1]"] == 100.0
    assert df.loc["Molar Concentration NaOH"]["stream_array[1]"] == 100.0
    assert df.loc["Molar Concentration EthylAcetate"]["stream_array[1]"] == 100.0
    assert df.loc["Molar Concentration SodiumAcetate"]["stream_array[1]"] == 100.0
    assert df.loc["Molar Concentration Ethanol"]["stream_array[1]"] == 100.0


@pytest.mark.unit
def test_create_stream_table_dataframe_from_StateBlock_true_state(m):
    df = create_stream_table_dataframe({
            "state": m.fs.tank1.control_volume.properties_out},
            true_state=True)

    assert df.loc["pressure"]["state"] == 101325
    assert df.loc["temperature"]["state"] == 298.15
    assert df.loc["flow_vol"]["state"] == 1.0
    assert df.loc["conc_mol_comp H2O"]["state"] == 100.0
    assert df.loc["conc_mol_comp NaOH"]["state"] == 100.0
    assert df.loc["conc_mol_comp EthylAcetate"]["state"] == 100.0
    assert df.loc["conc_mol_comp SodiumAcetate"]["state"] == 100.0
    assert df.loc["conc_mol_comp Ethanol"]["state"] == 100.0


@pytest.mark.unit
def test_create_stream_table_dataframe_from_StateBlock_orient(m):
    df = create_stream_table_dataframe({
            "state": m.fs.tank1.control_volume.properties_out},
            orient='index')

    assert df.loc["state"]["Pressure"] == 101325
    assert df.loc["state"]["Temperature"] == 298.15
    assert df.loc["state"]["Volumetric Flowrate"] == 1.0
    assert df.loc["state"]["Molar Concentration H2O"] == 100.0
    assert df.loc["state"]["Molar Concentration NaOH"] == 100.0
    assert df.loc["state"]["Molar Concentration EthylAcetate"] == 100.0
    assert df.loc["state"]["Molar Concentration SodiumAcetate"] == 100.0
    assert df.loc["state"]["Molar Concentration Ethanol"] == 100.0


@pytest.mark.unit
def test_create_stream_table_dataframe_from_StateBlock_time():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False, "time_set": [3]})
    m.fs.thermo_params = thermo_props.SaponificationParameterBlock()
    m.fs.reaction_params = rxn_props.SaponificationReactionParameterBlock(
        default={"property_package": m.fs.thermo_params})

    m.fs.tank1 = CSTR(default={"property_package": m.fs.thermo_params,
                               "reaction_package": m.fs.reaction_params})
    m.fs.tank2 = CSTR(default={"property_package": m.fs.thermo_params,
                               "reaction_package": m.fs.reaction_params})

    m.fs.stream = Arc(source=m.fs.tank1.outlet,
                      destination=m.fs.tank2.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    df = create_stream_table_dataframe({
            "state": m.fs.tank1.control_volume.properties_out},
            time_point=3)

    assert df.loc["Pressure"]["state"] == 101325
    assert df.loc["Temperature"]["state"] == 298.15
    assert df.loc["Volumetric Flowrate"]["state"] == 1.0
    assert df.loc["Molar Concentration H2O"]["state"] == 100.0
    assert df.loc["Molar Concentration NaOH"]["state"] == 100.0
    assert df.loc["Molar Concentration EthylAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration SodiumAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration Ethanol"]["state"] == 100.0


@pytest.mark.unit
def test_create_stream_table_dataframe_from_Port(m):
    df = create_stream_table_dataframe({
            "state": m.fs.tank1.outlet})

    assert df.loc["Pressure"]["state"] == 101325
    assert df.loc["Temperature"]["state"] == 298.15
    assert df.loc["Volumetric Flowrate"]["state"] == 1.0
    assert df.loc["Molar Concentration H2O"]["state"] == 100.0
    assert df.loc["Molar Concentration NaOH"]["state"] == 100.0
    assert df.loc["Molar Concentration EthylAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration SodiumAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration Ethanol"]["state"] == 100.0


@pytest.mark.unit
def test_create_stream_table_dataframe_from_Arc(m):
    df = create_stream_table_dataframe({
            "state": m.fs.stream})

    assert df.loc["Pressure"]["state"] == 101325
    assert df.loc["Temperature"]["state"] == 298.15
    assert df.loc["Volumetric Flowrate"]["state"] == 1.0
    assert df.loc["Molar Concentration H2O"]["state"] == 100.0
    assert df.loc["Molar Concentration NaOH"]["state"] == 100.0
    assert df.loc["Molar Concentration EthylAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration SodiumAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration Ethanol"]["state"] == 100.0


@pytest.mark.unit
def test_create_stream_table_dataframe_wrong_type(m):
    with pytest.raises(TypeError):
        create_stream_table_dataframe({"state": m.fs.tank1})


@pytest.mark.unit
def test_create_stream_table_dataframe_ordering(m):
    state_dict = {"state1": m.fs.stream,
                  "state3": m.fs.tank1.control_volume.properties_out,
                  "state2": m.fs.tank1.outlet}
    df = create_stream_table_dataframe(state_dict)

    columns = list(df)
    assert columns[0] == "state1"
    assert columns[1] == "state3"
    assert columns[2] == "state2"


@pytest.mark.unit
def test_stream_table_dataframe_to_string(m):
    df = create_stream_table_dataframe({
            "state": m.fs.tank1.control_volume.properties_out})

    stream_table_dataframe_to_string(df)


# -----------------------------------------------------------------------------
###
# Create a dummy StateBlock class
##
@declare_process_block_class("TestStateBlock", block_class=StateBlock)
class StateTestBlockData(StateBlockData):
    def build(self):
        self.x = Var(initialize=0)
        self.pressure = Var()
        self.enth_mol = Var()
        self.temperature = Expression(expr=self.enth_mol*2.0)
        self.div0 = Expression(expr=4.0/self.x)
        self.flow_mol = Var(["CO2", "H2O"])

@declare_process_block_class("TestStateBlock2", block_class=StateBlock)
class StateTestBlockData(StateBlockData):
    def build(self):
        self.pressure = Var()
        self.temperature = Var()
        self.flow_mol = Var()
        self.flow_vol = Var()

@pytest.fixture
def gtmodel():
    time_set = set([0,1])
    m = ConcreteModel()
    m.state_a = TestStateBlock(time_set)
    m.state_b = TestStateBlock(time_set, [1,2,3])
    m.state_c = TestStateBlock2(time_set)

    m.state_a[0].pressure = 11000
    m.state_a[0].enth_mol = 1100
    m.state_a[0].flow_mol["CO2"] = 110
    m.state_a[0].flow_mol["H2O"] = 111

    m.state_b[0,1].pressure = 10000
    m.state_b[0,1].enth_mol = 1000
    m.state_b[0,1].flow_mol["CO2"] = 100
    m.state_b[0,1].flow_mol["H2O"] = 101

    m.state_b[0,2].pressure = 20000
    m.state_b[0,2].enth_mol = 2000
    m.state_b[0,2].flow_mol["CO2"] = 200
    m.state_b[0,2].flow_mol["H2O"] = 201

    m.state_b[0,3].pressure = 30000
    m.state_b[0,3].enth_mol = 3000
    m.state_b[0,3].flow_mol["CO2"] = 300
    m.state_b[0,3].flow_mol["H2O"] = 301

    m.state_c[0].pressure = 1000
    m.state_c[0].flow_mol = 10
    m.state_c[0].flow_vol = 30
    m.state_c[0].temperature = 300

    return m

@pytest.mark.unit
def test_generate_table(gtmodel):
    m = gtmodel

    sd = {"a": m.state_a[0], "b1": m.state_b[(0,1)]}
    # This tests what happens if one of the requested attributes gives a division
    # by zero error and if one of the attributes doesn't exist.  With flow it
    # tests indexed attributes.
    st = generate_table(
        sd,
        attributes=["pressure", "temperature", "div0",
                    ("flow_mol", "CO2"), ("flow_mol", "H2O"),
                    "not_there", ("not_there_array", "hi")],
        heading=["P", "T", "ERR", "F_CO2", "F_H2O", "Miss", "Miss[hi]"])

    assert st.loc["a"]["P"] == 11000
    assert st.loc["a"]["F_CO2"] == 110
    assert st.loc["a"]["F_H2O"] == 111
    assert st.loc["a"]["T"] == 1100*2
    assert isna(st.loc["a"]["ERR"])
    assert isna(st.loc["a"]["Miss"])

    assert st.loc["b1"]["P"] == 10000
    assert st.loc["b1"]["F_CO2"] == 100
    assert st.loc["b1"]["F_H2O"] == 101
    assert st.loc["b1"]["T"] == 1000*2
    assert isna(st.loc["b1"]["ERR"])
    assert isna(st.loc["b1"]["Miss"])

@pytest.mark.unit
def test_generate_table_errors(gtmodel):
    m = gtmodel

    sd = {"a": m.state_a[0], "b1": m.state_b[0,1]}
    heading=["F"]

    with pytest.raises(AssertionError):
        st = generate_table(sd, attributes=[("flow_mol",)], heading=heading)

    with pytest.raises(KeyError):
        st = generate_table(sd, attributes=[("flow_mol", "coffee")], heading=heading)

    with pytest.raises(TypeError):
        st = generate_table(sd, attributes=["flow_mol"], heading=heading)

@pytest.mark.unit
def test_mixed_table(gtmodel):
    m = gtmodel

    sd = {
        "a": m.state_a[0],
        "b[1]": m.state_b[(0, 1)],
        "b[2]": m.state_b[(0, 2)],
        "b[3]": m.state_b[(0, 3)],
        "c": m.state_c[0],
    }


    st = generate_table(
        sd,
        exception=False,
        attributes=(
            "pressure",
            "temperature",
            "flow_mol",
            "flow_vol",
            "enth_mol",
            ("flow_mol", "CO2"),
            ("flow_mol", "H2O"),
        ),
        heading=(
            "P",
            "T",
            "F",
            "Fvol",
            "h",
            "F[CO2]",
            "F[H2O]",
        )
    )

    assert st.loc["a"]["P"] == 11000
    assert st.loc["a"]["F[CO2]"] == 110
    assert st.loc["a"]["F[H2O]"] == 111
    assert isna(st.loc["a"]["Fvol"])
    assert isna(st.loc["a"]["F"])
    assert st.loc["a"]["T"] == 1100*2

    assert st.loc["c"]["P"] == 1000
    assert isna(st.loc["c"]["F[CO2]"])
    assert isna(st.loc["c"]["F[H2O]"])
    assert st.loc["c"]["Fvol"] == 30
    assert st.loc["c"]["F"] == 10
    assert st.loc["c"]["T"] == 300


@pytest.mark.unit
def test_tag_states(gtmodel):
    m = gtmodel

    sd = {
        "a": m.state_a[0],
        "b[1]": m.state_b[(0, 1)],
        "b[2]": m.state_b[(0, 2)],
        "b[3]": m.state_b[(0, 3)],
        "c": m.state_c[0],
    }


    tags = tag_state_quantities(
        sd,
        attributes=(
            "pressure",
            "temperature",
            "flow_mol",
            "flow_vol",
            "enth_mol",
            ("flow_mol", "CO2"),
            ("flow_mol", "H2O"),
        ),
        labels=(
            "_pressure",
            "_temperature",
            "_flow_mol",
            "_flow_vol",
            "_enth_mol_differ",
            "_flow_mol[CO2]",
            "_flow_mol[H2O]",
        )
    )

    assert value(tags["a_pressure"]) == 11000
    assert value(tags["a_enth_mol_differ"]) == 1100
    assert value(tags["a_temperature"]) == 1100*2
    assert value(tags["a_flow_mol[CO2]"]) == 110
    assert value(tags["a_flow_mol[H2O]"]) == 111

    assert value(tags["b[1]_pressure"]) == 10000
    assert value(tags["b[1]_enth_mol_differ"]) == 1000
    assert value(tags["b[1]_flow_mol[CO2]"]) == 100
    assert value(tags["b[1]_flow_mol[H2O]"]) == 101

    assert value(tags["c_pressure"]) == 1000
    assert value(tags["c_temperature"]) == 300
    assert value(tags["c_flow_mol"]) == 10
    assert value(tags["c_flow_vol"]) == 30
    assert "c_flow_mol[H2O]" not in tags

    # check that I can change things
    tags["a_enth_mol_differ"].value = 1200
    assert value(m.state_a[0].enth_mol) == 1200
    assert value(m.state_a[0].temperature) == 1200*2
