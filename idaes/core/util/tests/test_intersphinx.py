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
from idaes.core.util import intersphinx as ix


@pytest.mark.unit
def test_get_intersphinx_mapping():
    mapping = ix.get_intersphinx_mapping()
    assert mapping
    assert mapping.get("idaes", None)


@pytest.mark.unit
def test_modify_url():
    m = ix.get_intersphinx_mapping()
    key = "idaes"
    assert key in m.keys()

    # modify one, version
    ix.modify_url(m, key, version="1.2.3")
    url = m[key][0]
    assert "/1.2.3" in url

    # modify one, language
    ix.modify_url(m, key, language="jp")
    url = m[key][0]
    assert "/jp" in url

    # modify all, language and version
    ix.modify_url(m, None, version="future", language="klingon")
    url = m[key][0]
    assert "/klingon/future" in url

    # bad key
    with pytest.raises(KeyError):
        ix.modify_url(m, "foobar!", version="x")
