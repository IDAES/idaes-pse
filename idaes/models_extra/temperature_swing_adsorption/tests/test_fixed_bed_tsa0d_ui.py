import pytest
from ..fixed_bed_tsa0d_ui import export_to_ui

@pytest.mark.component
def test_export():
    ui = export_to_ui()
    assert ui is not None
