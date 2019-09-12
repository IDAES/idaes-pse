import pytest

from idaes.dmf.ui import link_position_mapping


@pytest.mark.parametrize(
    "test_input,expected",
    [
        (
            "cstr",
            {
                "inlet_anchors": {
                    "name": "topLeft",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "bottomRight",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "equilibrium_reactor",
            {
                "inlet_anchors": {
                    "name": "topLeft",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "bottomRight",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "gibbs_reactor",
            {
                "inlet_anchors": {
                    "name": "topLeft",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "bottomRight",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "pfr",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "stoichiometric_reactor",
            {
                "inlet_anchors": {
                    "name": "topLeft",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "bottomRight",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "flash",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "top_outlet_anchor": {
                    "name": "top",
                    "args": {"rotate": "false", "padding": 0},
                },
                "bottom_outlet_anchor": {
                    "name": "bottom",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "mixer",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "feed",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "feed_flash",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "product",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "separator",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "heater",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "pressure_changer",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "heat_exchanger",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "heat_exchanger1d",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "outlet_anchors": {
                    "name": "right",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "packed_column",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "top_outlet_anchor": {
                    "name": "top",
                    "args": {"rotate": "false", "padding": 0},
                },
                "bottom_outlet_anchor": {
                    "name": "bottom",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
        (
            "tray_column",
            {
                "inlet_anchors": {
                    "name": "left",
                    "args": {"rotate": "false", "padding": 0},
                },
                "top_outlet_anchor": {
                    "name": "top",
                    "args": {"rotate": "false", "padding": 0},
                },
                "bottom_outlet_anchor": {
                    "name": "bottom",
                    "args": {"rotate": "false", "padding": 0},
                },
            },
        ),
    ],
)
def test_link_position_mapping(test_input, expected):
    assert link_position_mapping.link_position_mapping[test_input] == expected
