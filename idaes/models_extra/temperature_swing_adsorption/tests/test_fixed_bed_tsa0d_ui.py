import pytest
import logging
import pandas as pd
from pyomo.environ import check_optimal_termination
from pyomo.environ import value
from idaes.models_extra.temperature_swing_adsorption.util import tsa_summary
from idaes.models_extra.temperature_swing_adsorption.fixed_bed_tsa0d_ui import (
    export_to_ui,
    build,
    initialize,
    solve,
)
from idaes.core.util.tables import stream_table_dataframe_to_string

_log = logging.getLogger(__name__)


@pytest.mark.component
def test_export():
    ui = export_to_ui()
    assert ui is not None


@pytest.mark.component
def test_solve_base():
    ui = export_to_ui()
    ui.build(build_options=ui.fs_exp.build_options)
    r = ui.solve()
    assert check_optimal_termination(r)


@pytest.mark.component
def test_default_build_options():
    """test default build options, option values from jupyter notebook example"""
    default_build_options = {
        "adsorbent": "zeolite_13x",
        "number_of_beds": 1,
        "transformation_method": "dae.collocation",
        "transformation_scheme": "lagrangeRadau",
        "finite_elements": 20,
        "collocation_points": 6,
    }
    ui = export_to_ui()
    ui.build(build_options=ui.fs_exp.build_options)

    for key, val in ui.fs_exp.build_options.items():
        assert default_build_options[key] == val.value
    assert True


@pytest.mark.integration
def test_ui_output():
    """test ui data"""
    # get model options from UI export to ui
    ui_model_options = export_to_ui()
    build_options = ui_model_options.fs_exp.build_options

    # build and solve model
    m = build(build_options=build_options)
    initialize(m.fs)
    solve(m.fs)

    # base on UI output build data frame
    var_dict = m.fs.tsa.get_var_dict()
    ui_df = tsa_summary(m.fs.tsa)

    # base on jupyter notebook example data, build data frame
    jupyter_tsa_model_output_data = {
        "Adsorption temperature [K]": 310.00,
        "Desorption temperature [K]": 430.00,
        "Heating temperature [K]": 440.00,
        "Cooling temperature [K]": 300.00,
        "Column diameter [m]": 0.030000,
        "Column length [m]": 1.2000,
        "Column volume [m3]": 0.00084823,
        "CO2 mole fraction at feed [%]": 12.000,
        "Feed flow rate [mol/s]": 0.0096000,
        "Feed velocity [m/s]": 0.50008,
        "Minimum fluidization velocity [m/s]": 1.5207,
        "Time of heating step [h]": 0.37030,
        "Time of cooling step [h]": 0.20826,
        "Time of pressurization step [h]": 0.0051098,
        "Time of adsorption step [h]": 0.25221,
        "Cycle time [h]": 0.83588,
        "Purity [-]": 0.90219,
        "Recovery [-]": 0.89873,
        "Productivity [kg CO2/ton/h]": 84.085,
        "Specific energy [MJ/kg CO2]": 3.6532,
        "Heat duty per bed [MW]": 5.1244e-05,
        "Heat duty total [MW]": 0.00016646,
        "Pressure drop [Pa]": 5263.6,
        "Number of beds": 3.2484,
        "CO2 captured in one cycle per bed [kg/cycle]": 0.042210,
        "Cycles per year": 10480.0,
        "Total CO2 captured per year [tonne/year]": 1.4369,
        "Amount of flue gas processed per year [Gmol/year]": 0.00030275,
        "Amount of flue gas processed per year (target) [Gmol/year]": 0.00030275,
        "Amount of CO2 to atmosphere [mol/s]": 0.00011667,
        "Concentration of CO2 emitted to atmosphere [ppm]": 13803.0,
    }

    # build jupyter notebook data frame
    jupyter_df = pd.DataFrame.from_dict(
        jupyter_tsa_model_output_data, orient="index", columns=["Value"]
    )

    # compare data frame string
    jupyter_df_string = stream_table_dataframe_to_string(jupyter_df)
    ui_df_string = stream_table_dataframe_to_string(ui_df)
    assert jupyter_df_string == ui_df_string
