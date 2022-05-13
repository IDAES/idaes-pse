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
Tests for idaes.core.dmf.help module
"""
import pytest
from .util import tmp_dmf
from pathlib import Path
from idaes.core.dmf.dmfbase import DMFConfig
from idaes.core.dmf.help import get_html_docs


@pytest.mark.unit
def test_get_html_docs_noconfig(tmp_dmf):
    # rename config file to create failure
    assert DMFConfig.configuration_exists()
    p = DMFConfig.configuration_path()
    new_p = Path(str(p) + "-moved")
    p.rename(new_p)
    assert not DMFConfig.configuration_exists()
    with pytest.raises(ValueError):
        get_html_docs(tmp_dmf, "module", "name")  # parameters don't really matter
    # return file to original name
    new_p.rename(p)


@pytest.mark.unit
def test_get_html_docs_nohtml(tmp_dmf):
    assert DMFConfig.configuration_exists()
    with pytest.raises(ValueError):
        get_html_docs(tmp_dmf, "module", "name")  # parameters don't really matter
