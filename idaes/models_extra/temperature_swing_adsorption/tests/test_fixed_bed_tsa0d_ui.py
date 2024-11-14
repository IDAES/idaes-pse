import pytest
import logging
from pyomo.environ import check_optimal_termination
from idaes.models_extra.temperature_swing_adsorption.fixed_bed_tsa0d_ui import (
    export_to_ui,
    build,
    initialize,
    solve,
)

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
