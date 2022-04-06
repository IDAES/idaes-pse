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
"""
Tests get environment info
"""
import json
import pytest
import idaes.ver as ver
from idaes.core.util.env_info import EnvironmentInfo


@pytest.mark.unit
def test_env_info():
    x = EnvironmentInfo()
    d = x.to_dict()
    assert "IDAES" in d
    assert "Pyomo" in d
    assert "OS" in d
    assert hasattr(x, "package_version")
    assert x.version_string == ver.__version__
    s = x.to_json()
    d = json.loads(s)
    assert "IDAES" in d
    assert "Pyomo" in d
    assert "OS" in d
