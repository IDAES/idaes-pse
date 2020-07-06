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
import pytest

from idaes.ui.icon_mapping import icon_mapping


@pytest.mark.unit
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
    assert icon_mapping(test_input) == expected
