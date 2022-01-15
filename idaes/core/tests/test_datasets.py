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
Tests for idaes.core.datasets module
"""
import pytest
from idaes.core import datasets


@pytest.mark.unit
def test_pitzer():
    p = datasets.Pitzer()
    assert p
    assert p.list_tables()
    assert p.get_table(p.list_tables()[0])


@pytest.mark.unit
def test_available():
    av = datasets.available()
    assert len(av) >= 1


@pytest.mark.unit
def test_publication_unknown():
    # random, unknown publication
    with pytest.raises(KeyError):
        pub = datasets.Publication("test")


@pytest.mark.unit
def test_publication_known():
    # known publication
    pub = datasets.Publication("Pitzer:1984")
    assert pub
    assert pub.list_tables()
    assert pub.get_table(pub.list_tables()[0])
