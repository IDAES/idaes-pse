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
link_position_mapping = {
    "cstr": {
        "inlet_anchors": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {
            "name": "bottomRight",
            "args": {"rotate": "false", "padding": 0},
        },
    },
    "equilibrium_reactor": {
        "inlet_anchors": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {
            "name": "bottomRight",
            "args": {"rotate": "false", "padding": 0},
        },
    },
    "gibbs_reactor": {
        "inlet_anchors": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {
            "name": "bottomRight",
            "args": {"rotate": "false", "padding": 0},
        },
    },
    "pfr": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "stoichiometric_reactor": {
        "inlet_anchors": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {
            "name": "bottomRight",
            "args": {"rotate": "false", "padding": 0},
        },
    },
    "flash": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "top_outlet_anchor": {"name": "top", "args": {"rotate": "false", "padding": 0}},
        "bottom_outlet_anchor": {
            "name": "bottom",
            "args": {"rotate": "false", "padding": 0},
        },
    },
    "mixer": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "feed": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "feed_flash": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "product": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "separator": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "heater": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "pressure_changer": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "heat_exchanger": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    # TODO When we get to using different icons we need to figure out how
    #  to deal with mutiple inlets and outlets in different places
    "heat_exchanger1d": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}},
    },
    "packed_column": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "top_outlet_anchor": {"name": "top", "args": {"rotate": "false", "padding": 0}},
        "bottom_outlet_anchor": {
            "name": "bottom",
            "args": {"rotate": "false", "padding": 0},
        },
    },
    "tray_column": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}},
        "top_outlet_anchor": {"name": "top", "args": {"rotate": "false", "padding": 0}},
        "bottom_outlet_anchor": {
            "name": "bottom",
            "args": {"rotate": "false", "padding": 0},
        },
    },
    "default": {
        "inlet_anchors": {"name": "left", "args": {"rotate": "false", "padding": 0}}, 
        "outlet_anchors": {"name": "right", "args": {"rotate": "false", "padding": 0}}
    }
}
