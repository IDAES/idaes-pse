import pytest

from idaes.dmf.ui import icon_mapping


@pytest.mark.parametrize(
    "test_input,expected",
    [("cstr", "reactor_c.svg"),
    ("equilibrium_reactor", "reactor_e.svg"),
    ("gibbs_reactor", "reactor_g.svg"),
    ("plug_flow_reactor", "reactor_pfr.svg"),
    ("stoichiometric_reactor", "reactor_s.svg"),
    ("flash", "flash.svg"),
    ("mixer", "mixer.svg"),
    ("feed", "feed.svg"),
    ("feed_flash", "feed.svg"),
    ("product", "product.svg"),
    ("separator", "splitter.svg"), 
    ("heater", "heater_2.svg"),
    ("pressure_changer", "compressor.svg"),
    ("heat_exchanger", "heat_exchanger_1.svg"),
    ("heat_exchanger_1D", "heat_exchanger_1.svg"),
    ("statejunction", "NONE"),
    ("translator", "NONE"),
    ("packed_column", "packed_column_1.svg"),
    ("tray_column", "tray_column_1.svg")]
)
def test_icon_mapping(test_input, expected):
    assert icon_mapping.icon_mapping[test_input] == expected
