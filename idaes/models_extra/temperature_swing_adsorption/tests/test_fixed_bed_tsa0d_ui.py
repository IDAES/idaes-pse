import pytest
from pyomo.environ import check_optimal_termination
from ..fixed_bed_tsa0d_ui import export_to_ui, build


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
