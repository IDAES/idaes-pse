##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

def icon_mapping(unit_model_type):
    mapping = {
        "cstr": "reactor_c.svg",
        "flash": "flash.svg",
        "gibbs_reactor": "reactor_g.svg",
        "heat_exchanger": "heat_exchanger_1.svg",
        "heater": "heater_2.svg",
        "heat_exchanger_1D": "heat_exchanger_1.svg",  # same as HeatExchanger
        "mixer": "mixer.svg",
        "plug_flow_reactor": "reactor_pfr.svg",
        "pressure_changer": "compressor.svg",
        "separator": "splitter.svg",
        "stoichiometric_reactor": "reactor_s.svg",
        "equilibrium_reactor": "reactor_e.svg",
        "feed": "feed.svg",
        "product": "product.svg",
        "feed_flash": "feed.svg",
        "statejunction": "NONE",
        "translator": "NONE",
        "packed_column": "packed_column_1.svg",
        "tray_column": "tray_column_1.svg",
        "default": "default.svg"
    }
    try:
        return mapping[unit_model_type]
    except KeyError:
        return mapping["default"]



# TODO handle multi-icon mappings
#    "Heater": {"cooler.svg", "heater_1.svg", "heater_2.svg"},
#    "PressureChanger": {"compressor.svg", "expander.svg", "pump.svg", "fan.svg"},
# compressor and expander variations are identical when arrows are not present
#    "HeatExchanger": {"heat_exchanger_1.svg", "heat_exchanger_3.svg"},
# 1/2 and 3/4 are identical when arrows not present
#    "HeatExchanger1D": {"heat_exchanger_1.svg", "heat_exchanger_3.svg"},
#  same as HeadExchanger
#    "PackedColumn": {"packed_column_1.svg", "packed_column_2.svg",
#    "packed_column_3.svg", "packed_column_4.svg"},
#    "TrayColumn": {"tray_column_1.svg", "tray_column_2.svg", "tray_column_3.svg",
#    "tray_column_4.svg"},
#


# TODO need images, but the models don't exist yet
#    "StateJunction": "junction.svg", # small black circle "node"
#    "Translator": "translator.svg", # circle with a T in it
#    "Valve": {"valve_1.svg", "valve_2.svg", "valve_3.svg"}
