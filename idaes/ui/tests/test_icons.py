###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
#
###############################################################################
import pytest

from idaes.ui.icons import UnitModelIcon


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
@pytest.mark.unit
def test_icon_mapping(test_input, expected):
    assert UnitModelIcon(test_input).icon == expected

@pytest.mark.parametrize(
    "model_name,expected",
    [(
        "default",
        "{'groups': {'in': {'position': {'name': 'left', 'args': {'x': 2, 'y': 0, 'dx': 1, 'dy': 1}}, "
        "'attrs': {'rect': {'stroke': '#000000', 'stroke-width': 0, 'width': 0, 'height': 0}}, 'markup': "
        "'<g><rect/></g>'}, 'out': {'position': {'name': 'left', 'args': {'x': 48, 'y': 50, 'dx': 1, 'dy': 1}}, "
        "'attrs': {'rect': {'stroke': '#000000', 'stroke-width': 0, 'width': 0, 'height': 0}}, 'markup': "
        "'<g><rect/></g>'}}, 'items':[{'group': 'in', 'id': 'in'}, {'group': 'out', 'id': 'out'}]}"
    ),
    (
        "cstr",
        "{'groups': {'in': {'position': {'name': 'left', 'args': {'x': 15, 'y': 0, 'dx': 1, 'dy': 1}}, 'attrs': "
        "{'rect': {'stroke': '#000000', 'stroke-width': 0, 'width': 0, 'height': 0}}, 'markup': '<g><rect/></g>'}, "
        "'out': {'position': {'name': 'left', 'args': {'x': 48, 'y': 45, 'dx': 1, 'dy': 1}}, 'attrs': {'rect': "
        "{'stroke': '#000000', 'stroke-width': 0, 'width': 0, 'height': 0}}, 'markup': '<g><rect/></g>'}}, 'items': "
        "[{'group': 'in', 'id': 'in'}, {'group': 'out', 'id': 'out'}]}"
    )]
)
@pytest.mark.unit
def test_link_position_mapping(model_name, expected):
    expected_dict = eval(expected)
    positions = UnitModelIcon(model_name).link_positions
    assert positions == expected_dict
