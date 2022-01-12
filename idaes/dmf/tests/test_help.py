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
Tests for idaes.dmf.help module
"""
import pytest
from .util import tmp_dmf
from pathlib import Path
from idaes.dmf.dmfbase import DMFConfig
from idaes.dmf.help import get_html_docs


@pytest.mark.unit
def test_get_html_docs_noconfig(tmp_dmf):
    # rename config file to create failure
    assert DMFConfig.configuration_exists()
    path = DMFConfig.configuration_path()
    orig_path, new_path = (
        str(path),
        tmp_dmf.workspace_path / "test_get_html_docs_noconfig",
    )
    path = path.rename(new_path)
    assert not DMFConfig.configuration_exists()
    with pytest.raises(ValueError):
        get_html_docs(tmp_dmf, "module", "name")  # parameters don't really matter
    # return file to original name
    path = path.rename(orig_path)


@pytest.mark.unit
def test_get_html_docs_nohtml(tmp_dmf):
    assert DMFConfig.configuration_exists()
    with pytest.raises(ValueError):
        get_html_docs(tmp_dmf, "module", "name")  # parameters don't really matter
