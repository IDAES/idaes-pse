#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for the idaes.core.util.intersphinx module
"""
import pytest
from typing import Dict, Union
from collections import namedtuple
from idaes.core.util import intersphinx

Item = namedtuple("Item", ("key", "value", "url"))


def _first_item(m: Dict) -> Union[Item, None]:
    for key, value in m.items():
        return Item(key, value, value[0])
    return None


@pytest.mark.unit
def test_get_intersphinx_mapping():
    mapping = intersphinx.get_intersphinx_mapping()
    assert mapping
    assert mapping.get("idaes", None)


@pytest.mark.unit
def test_get_intersphinx_mapping_args():
    mapping = intersphinx.get_intersphinx_mapping("1.2.3")
    item = _first_item(mapping)
    assert "/1.2.3" in item.url

    mapping = intersphinx.get_intersphinx_mapping(language="jp")
    item = _first_item(mapping)
    assert "/jp/" in item.url

    mapping = intersphinx.get_intersphinx_mapping(language="jp", version="1.2.3")
    item = _first_item(mapping)
    assert "/jp/" in item.url
    assert "/1.2.3" in item.url
