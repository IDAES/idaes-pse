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
import pytest

from idaes.core.ui.icons.icons import UnitModelIcon


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ("default", "default.svg"),
        ("cstr", "reactor_c.svg"),
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
        ("statejunction", "default.svg"),
        ("translator", "default.svg"),
        ("packed_column", "packed_column_1.svg"),
        ("tray_column", "tray_column_1.svg"),
    ],
)
@pytest.mark.unit
def test_icon_mapping(test_input, expected):
    assert UnitModelIcon(test_input).icon == expected
    with pytest.raises(ValueError):
        UnitModelIcon("unregistered_model", "not_default")


@pytest.mark.parametrize(
    "model_name,expected",
    [
        (
            "default",
            "{'groups': {'in': {'position': {'name': 'left', 'args': {'x': 10, 'y': 35}}, "
            "'attrs': {'rect': {'stroke': '#000000', 'stroke-width': 0, 'width': 0, 'height': 0}}, 'markup': "
            "'<g><rect/></g>'}, 'out': {'position': {'name': 'left', 'args': {'x': 41, 'y': 35}}, "
            "'attrs': {'rect': {'stroke': '#000000', 'stroke-width': 0, 'width': 0, 'height': 0}}, 'markup': "
            "'<g><rect/></g>'}}, 'items':[]}",
        ),
        (
            "cstr",
            "{'groups': {'in': {'position': {'name': 'left', 'args': {'x': 15, 'y': 0, 'dx': 1, 'dy': 1}}, 'attrs': "
            "{'rect': {'stroke': '#000000', 'stroke-width': 0, 'width': 0, 'height': 0}}, 'markup': '<g><rect/></g>'}, "
            "'out': {'position': {'name': 'left', 'args': {'x': 48, 'y': 45, 'dx': 1, 'dy': 1}}, 'attrs': {'rect': "
            "{'stroke': '#000000', 'stroke-width': 0, 'width': 0, 'height': 0}}, 'markup': '<g><rect/></g>'}}, 'items': []}",
        ),
    ],
)
@pytest.mark.unit
def test_link_position_mapping(model_name, expected):
    expected_dict = eval(expected)
    positions = UnitModelIcon(model_name).link_positions
    assert positions == expected_dict
