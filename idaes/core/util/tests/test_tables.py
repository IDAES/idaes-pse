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
    value,
    units as pyunits,
)
from pyomo.network import Arc, Port
from idaes.core import (
    FlowsheetBlock,
    StateBlock,
    StateBlockData,
    declare_process_block_class,
)
from idaes.core.util.tables import (
    arcs_to_stream_dict,
    create_stream_table_dataframe,
    create_stream_table_ui,
    stream_table_dataframe_to_string,
    generate_table,
    stream_states_dict,
)
import idaes.models.properties.examples.saponification_thermo as thermo_props
import idaes.models.properties.examples.saponification_reactions as rxn_props
from idaes.models.unit_models import CSTR, Flash
from idaes.models.unit_models.heat_exchanger_1D import HeatExchanger1D as HX1D
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.models_extra.column_models import TrayColumn
from idaes.models_extra.column_models.condenser import CondenserType, TemperatureSpec
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)


@pytest.fixture()
def m():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.thermo_params = thermo_props.SaponificationParameterBlock()
    m.fs.reaction_params = rxn_props.SaponificationReactionParameterBlock(
        property_package=m.fs.thermo_params
    )

    m.fs.tank1 = CSTR(
        property_package=m.fs.thermo_params, reaction_package=m.fs.reaction_params
    )
    m.fs.tank2 = CSTR(
        property_package=m.fs.thermo_params, reaction_package=m.fs.reaction_params
    )
    m.fs.tank_array = CSTR(
        range(3),
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
    )

    m.fs.stream = Arc(source=m.fs.tank1.outlet, destination=m.fs.tank2.inlet)

    def stream_array_rule(b, i):
        return (b.tank_array[i].outlet, b.tank_array[i + 1].inlet)

    m.fs.stream_array = Arc(range(2), rule=stream_array_rule)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


@pytest.fixture()
def m_with_variable_types():
    """Flash unit model. Use '.fs' attribute to get the flowsheet."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # Flash properties
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
    )
    # Flash unit
    m.fs.flash = Flash(property_package=m.fs.properties)
    # Adding fixed and unfixed variables
    m.fs.flash.inlet.pressure.fix(3.14)
    m.fs.flash.inlet.pressure.unfix()
    m.fs.flash.inlet.temperature.fix(368)

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
    df = create_stream_table_dataframe(
        {"state": m.fs.tank1.control_volume.properties_out}, true_state=True
    )

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
    df = create_stream_table_dataframe(
        {"state": m.fs.tank1.control_volume.properties_out}, orient="index"
    )

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
    m.fs = FlowsheetBlock(dynamic=False, time_set=[3])
    m.fs.thermo_params = thermo_props.SaponificationParameterBlock()
    m.fs.reaction_params = rxn_props.SaponificationReactionParameterBlock(
        property_package=m.fs.thermo_params
    )

    m.fs.tank1 = CSTR(
        property_package=m.fs.thermo_params, reaction_package=m.fs.reaction_params
    )
    m.fs.tank2 = CSTR(
        property_package=m.fs.thermo_params, reaction_package=m.fs.reaction_params
    )

    m.fs.stream = Arc(source=m.fs.tank1.outlet, destination=m.fs.tank2.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    df = create_stream_table_dataframe(
        {"state": m.fs.tank1.control_volume.properties_out}, time_point=3
    )

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
    df = create_stream_table_dataframe({"state": m.fs.tank1.outlet})

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
    df = create_stream_table_dataframe({"state": m.fs.stream})

    assert df.loc["Pressure"]["state"] == 101325
    assert df.loc["Temperature"]["state"] == 298.15
    assert df.loc["Volumetric Flowrate"]["state"] == 1.0
    assert df.loc["Molar Concentration H2O"]["state"] == 100.0
    assert df.loc["Molar Concentration NaOH"]["state"] == 100.0
    assert df.loc["Molar Concentration EthylAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration SodiumAcetate"]["state"] == 100.0
    assert df.loc["Molar Concentration Ethanol"]["state"] == 100.0


@pytest.mark.unit
def test_create_stream_table_ui(m_with_variable_types):
    m = m_with_variable_types

    state_name = "state"
    state_dict = {state_name: m.fs.flash.inlet}
    df = create_stream_table_ui(state_dict)

    assert df.loc["pressure"][state_name] == (3.14, "unfixed")
    assert df.loc["temperature"][state_name] == (368, "fixed")


@pytest.mark.unit
def test_create_stream_table_dataframe_wrong_type(m):
    with pytest.raises(TypeError):
        create_stream_table_dataframe({"state": m.fs.tank1})


@pytest.mark.unit
def test_create_stream_table_dataframe_ordering(m):
    state_dict = {
        "state1": m.fs.stream,
        "state3": m.fs.tank1.control_volume.properties_out,
        "state2": m.fs.tank1.outlet,
    }
    df = create_stream_table_dataframe(state_dict)

    columns = list(df)
    assert columns[0] == "Units"
    assert columns[1] == "state1"
    assert columns[2] == "state3"
    assert columns[3] == "state2"


@pytest.mark.unit
def test_stream_table_dataframe_to_string(m):
    df = create_stream_table_dataframe(
        {"state": m.fs.tank1.control_volume.properties_out}
    )

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
        self.temperature = Expression(expr=self.enth_mol * 2.0)
        self.div0 = Expression(expr=4.0 / self.x)
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
    time_set = set([0, 1])
    m = ConcreteModel()
    m.state_a = TestStateBlock(time_set)
    m.state_b = TestStateBlock(time_set, [1, 2, 3])
    m.state_c = TestStateBlock2(time_set)

    m.state_a[0].pressure = 11000
    m.state_a[0].enth_mol = 1100
    m.state_a[0].flow_mol["CO2"] = 110
    m.state_a[0].flow_mol["H2O"] = 111

    m.state_b[0, 1].pressure = 10000
    m.state_b[0, 1].enth_mol = 1000
    m.state_b[0, 1].flow_mol["CO2"] = 100
    m.state_b[0, 1].flow_mol["H2O"] = 101

    m.state_b[0, 2].pressure = 20000
    m.state_b[0, 2].enth_mol = 2000
    m.state_b[0, 2].flow_mol["CO2"] = 200
    m.state_b[0, 2].flow_mol["H2O"] = 201

    m.state_b[0, 3].pressure = 30000
    m.state_b[0, 3].enth_mol = 3000
    m.state_b[0, 3].flow_mol["CO2"] = 300
    m.state_b[0, 3].flow_mol["H2O"] = 301

    m.state_c[0].pressure = 1000
    m.state_c[0].flow_mol = 10
    m.state_c[0].flow_vol = 30
    m.state_c[0].temperature = 300

    return m


@pytest.mark.unit
def test_generate_table(gtmodel):
    m = gtmodel

    sd = {"a": m.state_a[0], "b1": m.state_b[(0, 1)]}
    # This tests what happens if one of the requested attributes gives a division
    # by zero error and if one of the attributes doesn't exist.  With flow it
    # tests indexed attributes.
    st = generate_table(
        sd,
        attributes=[
            "pressure",
            "temperature",
            "div0",
            ("flow_mol", "CO2"),
            ("flow_mol", "H2O"),
            "not_there",
            ("not_there_array", "hi"),
        ],
        heading=["P", "T", "ERR", "F_CO2", "F_H2O", "Miss", "Miss[hi]"],
    )

    assert st.loc["a"]["P"] == 11000
    assert st.loc["a"]["F_CO2"] == 110
    assert st.loc["a"]["F_H2O"] == 111
    assert st.loc["a"]["T"] == 1100 * 2
    assert isna(st.loc["a"]["ERR"])
    assert isna(st.loc["a"]["Miss"])

    assert st.loc["b1"]["P"] == 10000
    assert st.loc["b1"]["F_CO2"] == 100
    assert st.loc["b1"]["F_H2O"] == 101
    assert st.loc["b1"]["T"] == 1000 * 2
    assert isna(st.loc["b1"]["ERR"])
    assert isna(st.loc["b1"]["Miss"])


@pytest.mark.unit
def test_generate_table_errors(gtmodel):
    m = gtmodel

    sd = {"a": m.state_a[0], "b1": m.state_b[0, 1]}
    heading = ["F"]

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
        ),
    )

    assert st.loc["a"]["P"] == 11000
    assert st.loc["a"]["F[CO2]"] == 110
    assert st.loc["a"]["F[H2O]"] == 111
    assert isna(st.loc["a"]["Fvol"])
    assert isna(st.loc["a"]["F"])
    assert st.loc["a"]["T"] == 1100 * 2

    assert st.loc["c"]["P"] == 1000
    assert isna(st.loc["c"]["F[CO2]"])
    assert isna(st.loc["c"]["F[H2O]"])
    assert st.loc["c"]["Fvol"] == 30
    assert st.loc["c"]["F"] == 10
    assert st.loc["c"]["T"] == 300


@pytest.fixture()
def HX1D_array_model():
    # An example of maximum perversity. Dynamics, 1D control volumes, and
    # an indexed unit
    unit_set = range(3)
    time_set = [0, 5]
    time_units = pyunits.s
    fs_cfg = {"dynamic": True, "time_set": time_set, "time_units": time_units}

    m = ConcreteModel()
    m.fs = FlowsheetBlock(**fs_cfg)

    m.fs.properties = thermo_props.SaponificationParameterBlock()

    m.fs.unit_array = HX1D(
        unit_set,
        hot_side_name="shell_side",
        cold_side_name="tube_side",
        shell_side={"property_package": m.fs.properties},
        tube_side={"property_package": m.fs.properties},
    )

    def tube_stream_array_rule(b, i):
        return {
            "source": b.unit_array[i].tube_side_outlet,
            "destination": b.unit_array[i + 1].tube_side_inlet,
        }

    m.fs.tube_stream_array = Arc(range(2), rule=tube_stream_array_rule)

    def shell_stream_array_rule(b, i):
        return {
            "source": b.unit_array[i + 1].shell_side_outlet,
            "destination": b.unit_array[i].shell_side_inlet,
        }

    m.fs.shell_stream_array = Arc(range(1, -1, -1), rule=shell_stream_array_rule)

    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


@pytest.mark.unit
def test_create_stream_table_dataframe_from_Port_HX1D(HX1D_array_model):
    m = HX1D_array_model
    df = create_stream_table_dataframe({"state": m.fs.unit_array[0].tube_side_inlet})

    assert df.loc["Pressure"]["state"] == pytest.approx(101325)
    assert df.loc["Temperature"]["state"] == pytest.approx(298.15)
    assert df.loc["Volumetric Flowrate"]["state"] == pytest.approx(1.0)
    assert df.loc["Molar Concentration H2O"]["state"] == pytest.approx(100.0)
    assert df.loc["Molar Concentration NaOH"]["state"] == pytest.approx(100.0)
    assert df.loc["Molar Concentration EthylAcetate"]["state"] == pytest.approx(100.0)
    assert df.loc["Molar Concentration SodiumAcetate"]["state"] == pytest.approx(100.0)
    assert df.loc["Molar Concentration Ethanol"]["state"] == pytest.approx(100.0)

    df = create_stream_table_dataframe({"state": m.fs.unit_array[0].tube_side_outlet})

    assert df.loc["Pressure"]["state"] == pytest.approx(101325)
    assert df.loc["Temperature"]["state"] == pytest.approx(298.15)
    assert df.loc["Volumetric Flowrate"]["state"] == pytest.approx(1.0)
    assert df.loc["Molar Concentration H2O"]["state"] == pytest.approx(100.0)
    assert df.loc["Molar Concentration NaOH"]["state"] == pytest.approx(100.0)
    assert df.loc["Molar Concentration EthylAcetate"]["state"] == pytest.approx(100.0)
    assert df.loc["Molar Concentration SodiumAcetate"]["state"] == pytest.approx(100.0)
    assert df.loc["Molar Concentration Ethanol"]["state"] == pytest.approx(100.0)


@pytest.mark.unit
def test_create_stream_table_dataframe_from_Arc_HX1D(HX1D_array_model):
    m = HX1D_array_model
    df = create_stream_table_dataframe({"state": m.fs.tube_stream_array}, time_point=5)

    for i in range(2):
        stg = f"state[{i}]"
    assert df.loc["Pressure"][stg] == pytest.approx(101325)
    assert df.loc["Temperature"][stg] == pytest.approx(298.15)
    assert df.loc["Volumetric Flowrate"][stg] == pytest.approx(1.0)
    assert df.loc["Molar Concentration H2O"][stg] == pytest.approx(100.0)
    assert df.loc["Molar Concentration NaOH"][stg] == pytest.approx(100.0)
    assert df.loc["Molar Concentration EthylAcetate"][stg] == pytest.approx(100.0)
    assert df.loc["Molar Concentration SodiumAcetate"][stg] == pytest.approx(100.0)
    assert df.loc["Molar Concentration Ethanol"][stg] == pytest.approx(100.0)


@pytest.fixture()
def flash_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()  # default={"valid_phase":
    # ('Liq', 'Vap')})
    m.fs.flash = Flash(property_package=m.fs.properties)

    return m


@pytest.mark.unit
def test_state_block_retrieval_fail(flash_model):
    # The flash unit does not have real state blocks associated with the
    # outlet ports. There is a mixture of references, expressions, and multiple
    # state blocks. Therefore we don't want any state block getting through.
    m = flash_model
    with pytest.raises(
        RuntimeError,
        match="No block could be retrieved from Port fs.flash.liq_outlet "
        "because components are derived from multiple blocks.",
    ):
        df = create_stream_table_dataframe({"state": m.fs.flash.liq_outlet})


@pytest.mark.unit
def test_state_block_retrieval_empty_port():
    m = ConcreteModel()
    m.p = Port()
    with pytest.raises(
        ValueError,
        match="No block could be retrieved from Port p because it contains "
        "no components.",
    ):
        df = create_stream_table_dataframe({"state": m.p})


@pytest.fixture()
def distillation_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
    )
    m.fs.unit = TrayColumn(
        number_of_trays=3,
        feed_tray_location=2,
        condenser_type=CondenserType.totalCondenser,
        condenser_temperature_spec=TemperatureSpec.atBubblePoint,
        property_package=m.fs.properties,
    )
    return m


@pytest.mark.unit
def test_extended_port_retrieval(distillation_model):
    m = distillation_model
    df = create_stream_table_dataframe({"state": m.fs.unit.feed})
    stg = "state"
    assert df.loc["pressure"][stg] == pytest.approx(101325)
    assert df.loc["temperature"][stg] == pytest.approx(298.15)
    assert df.loc["mole_frac_comp benzene"][stg] == pytest.approx(0.5)
    assert df.loc["mole_frac_comp toluene"][stg] == pytest.approx(0.5)
    assert df.loc["flow_mol"][stg] == pytest.approx(1)
