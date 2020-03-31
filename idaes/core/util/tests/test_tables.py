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

from pandas import isna
import pytest

from pyomo.environ import (ConcreteModel,
                           Expression,
                           TransformationFactory,
                           Var)
from pyomo.network import Arc
from idaes.core import (FlowsheetBlock,
                        StateBlock,
                        StateBlockData,
                        declare_process_block_class)
from idaes.core.util.tables import (arcs_to_stream_dict,
                                    create_stream_table_dataframe,
                                    stream_table_dataframe_to_string,
                                    generate_table)

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

    m.fs.stream = Arc(source=m.fs.tank1.outlet,
                      destination=m.fs.tank2.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def test_create_stream_table_dataframe_from_StateBlock(m):
    d = arcs_to_stream_dict(m, descend_into=True)
    assert "stream" in d
    assert d["stream"] == m.fs.stream


def test_create_stream_table_dataframe_from_StateBlock_2(m):
    df = create_stream_table_dataframe({
            "state": m.fs.tank1.control_volume.properties_out})

    assert df.loc["Pressure"]["state"] == 101325
    assert df.loc["Temperature"]["state"] == 298.15
    assert df.loc["Volumetric Flowrate"]["state"] == 1.0
    assert df.loc["Molar Concentration H2O"]["state"] == 100.0
    assert df.loc["Molar Concentration NaOH"]["state"] == 100.0
    assert df.loc["Molar Concentration EthylAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration SodiumAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration Ethanol"]["state"] == 100.0


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


def test_create_stream_table_dataframe_wrong_type(m):
    with pytest.raises(TypeError):
        create_stream_table_dataframe({"state": m.fs.tank1})


def test_create_stream_table_dataframe_ordering(m):
    state_dict = {"state1": m.fs.stream,
                  "state3": m.fs.tank1.control_volume.properties_out,
                  "state2": m.fs.tank1.outlet}
    df = create_stream_table_dataframe(state_dict)

    columns = list(df)
    assert columns[0] == "state1"
    assert columns[1] == "state3"
    assert columns[2] == "state2"


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
        self.flow_mol = Var()


def test_generate_table():
    m = ConcreteModel()
    m.state_a = TestStateBlock()
    m.state_b = TestStateBlock([1,2,3])

    m.state_a.pressure = 11000
    m.state_a.enth_mol = 1100
    m.state_a.flow_mol = 110

    m.state_b[1].pressure = 10000
    m.state_b[1].enth_mol = 1000
    m.state_b[1].flow_mol = 100

    m.state_b[2].pressure = 20000
    m.state_b[2].enth_mol = 2000
    m.state_b[2].flow_mol = 200

    m.state_b[3].pressure = 30000
    m.state_b[3].enth_mol = 3000
    m.state_b[3].flow_mol = 300

    sd = {"a": m.state_a, "b1": m.state_b[1]}
    st = generate_table(
        sd,
        attributes=["pressure", "temperature", "div0",
                    "flow_mol", "not_there"],
        heading=["P", "T", "ERR", "F", "Miss"])

    assert(st.loc["a"]["P"] == 11000)
    assert(abs(st.loc["a"]["T"] - 1100*2) < 0.001)
    assert(isna(st.loc["a"]["ERR"]))
    assert(isna(st.loc["a"]["Miss"]))

    assert(st.loc["b1"]["P"] == 10000)
    assert(abs(st.loc["b1"]["T"] - 1000*2) < 0.001)
    assert(isna(st.loc["b1"]["ERR"]))
    assert(isna(st.loc["b1"]["Miss"]))
