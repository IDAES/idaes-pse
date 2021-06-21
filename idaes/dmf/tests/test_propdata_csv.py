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
Test CSV merging for property data
"""
# stdlib
import logging
import sys
from io import StringIO

# third-party
import pytest

# local
from idaes.dmf.propdata import AddedCSVColumnError
from idaes.dmf.propdata import PropertyData as PropData

# for testing
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)

# ------------------
# Setup and teardown
# ------------------

test_data = {
    "main": "Num,State (units 1),Typeof1 Error,Prop (units 2),Typeof2 Error\n"
    "1, 1.0, 0, 100.0, 0.1\n"
    "2, 2.0, 0, 200.0, 0.2\n"
    "3, 3.0, 0, 300.0, 0.3\n",
    "same": [
        "Num,State (units 1),Typeof1 Error,Prop (units 2),Typeof2 Error\n"
        "4, 4.0, 0, 400.0, 0.4\n"
        "5, 5.0, 0, 500.0, 0.5\n",
        "Num,State (units 1),Typeof1 Error,Prop (units 2),Typeof2 Error\n"
        "6, 6.0, 0, 600.0, 0.6\n"
        "7, 7.0, 0, 700.0, 0.7\n",
    ],
    "missing": "Num,State (units 1),Typeof1 Error\n" "4, 4.0, 0\n" "5, 5.0, 0\n",
    "extra": "Num,State (units 1),Typeof1 Error,Prop (units 2),Typeof2 Error,"
    "PropB (units B2), Typeof3 Error\n"
    "4, 4.0, 0, 400.0, 0.3, -400.0, -0.4\n"
    "5, 5.0, 0, 500.0, 0.5, -500.0, -0.5\n",
    "diff": "Num,HELLO (units 1),Typeof1 Error,WORLD (units 2),Typeof2 Error\n"
    "4, 4.0, 0, 400.0, 0.4\n"
    "5, 5.0, 0, 500.0, 0.5\n",
    "alt": "Num,State (units 1),Typeof1 Error,PropB (units 2),Typeof1B Error\n"
    "4, 4.0, 0, 400.0, 0.4\n"
    "5, 5.0, 0, 500.0, 0.5\n",
    "more1": "Num,State (units 1),Typeof1 Error,Prop (units 2),Typeof2 Error,"
    "PropB (units 2),Typeof1B Error\n"
    "4, 4.0, 0, 400.0, 0.4, -400.0, -0.4\n"
    "5, 5.0, 0, 500.0, 0.5, -500.0, -0.5\n",
    "more2": "Num,State (units 1),Typeof1 Error,"
    "PropB (units 2),Typeof1B Error,"
    "PropC (units 2),Typeof2B Error,"
    "PropD (units 2),Typeof3B Error\n"
    "6, 6.0, 0, -600.0, 0.6, 60.0, 0.6, -60.0, -0.6\n"
    "7, 7.0, 0, -700.0, 0.7, 70.0, 0.7, -70.0, -0.7\n",
}


@pytest.fixture
def pd():
    # create main propertydata obj
    return PropData.from_csv(StringIO(test_data["main"]), 1)


# -----
# Tests
# -----


@pytest.mark.unit
def test_same(pd):
    for chunk in test_data["same"]:
        n = pd.num_rows  # previous length

        # add the new data
        num_added = pd.add_csv(StringIO(chunk))

        # Check that rows were appended to PropertyData
        assert pd.num_rows == n + num_added

        # Check that values match
        df = pd.values_dataframe()
        print("@@ got dataframe:\n{}".format(df))
        for name in pd.names():
            pd_values = list(df[name])[n:]  # get new values
            multp = 1.0 if name == "State" else 100.0
            for i, value in enumerate(pd_values):
                assert value == multp * (i + n + 1)


@pytest.mark.unit
def test_extra_strict(pd):
    with pytest.raises(AddedCSVColumnError) as einfo:
        pd.add_csv(StringIO(test_data["extra"]), strict=True)
    print("Extra columns gave expected error:\n(MSG) {}".format(einfo.value))


@pytest.mark.unit
def test_extra(pd):
    pd.add_csv(StringIO(test_data["extra"]))
    assert "PropB" in pd.names()


@pytest.mark.unit
def test_missing_strict(pd):
    with pytest.raises(AddedCSVColumnError) as einfo:
        pd.add_csv(StringIO(test_data["missing"]), strict=True)
    print("Missing columns gave expected error:\n(MSG) {}".format(einfo.value))


@pytest.mark.unit
def test_missing(pd):
    n = pd.add_csv(StringIO(test_data["missing"]))
    assert n == 0, "Adding no new columns should have returned 0"


@pytest.mark.unit
def test_different_strict(pd):
    with pytest.raises(AddedCSVColumnError) as einfo:
        pd.add_csv(StringIO(test_data["diff"]), strict=True)
    print("Different columns gave expected error:\n(MSG) {}".format(einfo.value))


@pytest.mark.unit
def test_different(pd):
    with pytest.raises(AddedCSVColumnError) as einfo:
        pd.add_csv(StringIO(test_data["diff"]))
    print("Different columns gave expected error:\n(MSG) {}".format(einfo.value))


@pytest.mark.unit
def test_alt_strict(pd):
    with pytest.raises(AddedCSVColumnError) as einfo:
        pd.add_csv(StringIO(test_data["alt"]), strict=True)
    print("Different columns gave expected error:\n(MSG) {}".format(einfo.value))


@pytest.mark.unit
def test_alt(pd):
    assert "PropB" not in pd.names()
    n = pd.add_csv(StringIO(test_data["alt"]))
    assert n == 2, "Should have returned 2 rows added"
    assert "PropB" in pd.names()


@pytest.mark.unit
def test_more(pd):
    assert "PropB" not in pd.names()
    pd_rows = pd.num_rows
    n = pd.add_csv(StringIO(test_data["more1"]))
    assert n == 2, "Should have returned 2 rows added"
    assert pd.num_rows == pd_rows + 2
    assert "PropB" in pd.names()
    assert "PropC" not in pd.names()
    n = pd.add_csv(StringIO(test_data["more2"]))
    assert n == 2, "Should have returned 2 rows added"
    assert pd.num_rows == pd_rows + 4
    for ltr in "C", "D":
        label = "Prop" + ltr
        assert label in pd.names()
