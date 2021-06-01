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
Tests for idaes.dmf.propdb.types
"""
# stdlib
from io import StringIO
import json
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Union

# third-party
import pytest

# local
from idaes.dmf.propdata import PropertyTable, PropertyData, PropertyMetadata

# for testing
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)

scratch_dir: Union[str, None] = None
scratch_path: Union[Path, None] = None


def setup_module(module):
    global scratch_dir, scratch_path
    scratch_dir = TemporaryDirectory(prefix="idaes_dmf_")  # easier to remove later
    scratch_path = Path(scratch_dir.name)


def teardown_module(module):
    global scratch_dir
    del scratch_dir


# keep minimal example first for test_Property{Data,Metadata}Parsing
good_data_csv = [
    # minimal
    "Num,State (units 1),Absolute Error,Prop (units 2),Relative Error\n"
    "1,1.0,0,100.0,0.1",
    # two values, "T" and "Density"
    "Data No.,T (K),Absolute Error,Density (g/cm3),Absolute Error\n"
    "1,303.15,0,0.2,0\n"
    "2,304.15,0,0.3,0\n",
    # same but with missing values
    "Data No.,T (K),Absolute Error,Density (g/cm3),Absolute Error\n"
    "1,303.15,,0.2,0\n"
    "2,304.15,0,,0\n",
]

bad_data_csv = ["", "blah, blah", "blah\nblah", "N,X\n1,1,0"]  # missing error in hdr

good_metadata_text = [
    # minimal
    'source , Mr. Author,"A Title",Other Info',
    # Full set of fields
    'Source,Han, J., Jin, J., Eimer, D.A., Melaaen, M.C.,"Density of '
    "Water(1) + Monoethanolamine(2) + CO2(3) from (298.15 to 413.15) "
    "K and Surface Tension of Water(1) + Monethanolamine(2) from "
    '(303.15 to 333.15)K", J. Chem. Eng. Data, 2012, Vol. 57, pg. '
    '1095-1103"\n'
    'Retrieval,"J. Morgan, date unknown"\n'
    "Notes,r is MEA weight fraction in aqueous soln. (CO2-free basis)",
]

good_table_json = [
    {
        "meta": [
            {
                "datatype": "MEA",
                "info": "J. Chem. Eng. Data, 2009, Vol 54, pg. 3096-30100",
                "notes": "r is MEA weight fraction in aqueous soln.",
                "authors": "Amundsen, T.G., Lars, E.O., Eimer, D.A.",
                "title": "Density and Viscosity of Monoethanolamine + .etc.",
                "date": "2009",
            }
        ],
        "data": [
            {
                "name": "Viscosity Value",
                "units": "mPa-s",
                "values": [2.6, 6.2],
                "error_type": "absolute",
                "errors": [0.06, 0.004],
                "type": "property",
            },
            {"name": "r", "units": "", "values": [0.2, 1000], "type": "state"},
        ],
    }
]

dummy_data = [{"name": "dummy", "type": "state", "units": "none", "values": [1, 2, 3]}]


@pytest.mark.unit
def test_property_table_empty():
    with pytest.raises(ValueError):
        PropertyTable(data=[], metadata=[])


@pytest.mark.unit
def test_property_table():
    tbl = PropertyTable(data=dummy_data, metadata=[])
    assert tbl.data is not None
    assert tbl.metadata is not None


@pytest.mark.unit
def test_property_table_objinit():
    d = PropertyData(dummy_data)
    m = PropertyMetadata({"bar": "y"})
    tbl = PropertyTable(data=d, metadata=m)
    assert tbl.data.as_list()[0]["units"] == "none"
    assert tbl.metadata[0].as_dict()["bar"] == "y"


@pytest.mark.unit
def test_property_table_json():
    tbl = PropertyTable(data=dummy_data, metadata={"baz": "2"})
    s = tbl.dumps()
    parsed = json.loads(s)
    assert parsed["data"] == tbl.data.as_list()
    tmp_dir = scratch_path / "property_table_json"
    tmp_dir.mkdir()
    tbl.dump((tmp_dir / "tbl.json").open("w"))
    parsed = json.load((tmp_dir / "tbl.json").open())
    assert parsed["data"] == tbl.data.as_list()


@pytest.mark.unit
def test_property_table_from_json():
    tmp_dir = scratch_path / "property_table_from_json"
    tmp_dir.mkdir()
    prop_path = tmp_dir / "prop.json"
    for gj in good_table_json:
        json.dump(gj, prop_path.open("w"))
        tbl = PropertyTable.load(prop_path.open())
        assert tbl.data.num_columns == 2


@pytest.mark.unit
def test_property_data_good():
    inputs = map(StringIO, good_data_csv)
    states, total = [1], [2]
    for inpfile, state_col, ncol in zip(inputs, states, total):
        obj = PropertyData.from_csv(inpfile, state_col)
        assert obj.num_columns == ncol


@pytest.mark.unit
def test_property_data_column_type():
    csv_input = good_data_csv[0]
    obj = PropertyData.from_csv(StringIO(csv_input), 1)
    csv_header = csv_input.split("\n")[0].split(",")
    for i, col in enumerate((csv_header[1], csv_header[3])):
        units_pos = col.find("(")
        if units_pos >= 0:
            col = col[:units_pos].strip()
        obj.get_column(col)
        index = obj.get_column_index(col)
        if i == 0:
            assert obj.is_state_column(index)
        else:
            assert obj.is_property_column(index)


@pytest.mark.unit
def test_property_data_column_nonexistent():
    csv_input = good_data_csv[0]
    obj = PropertyData.from_csv(StringIO(csv_input), 1)
    try:
        obj.get_column("this is not a column")
        assert False, "Expected failure for bogus column"
    except KeyError:
        pass


@pytest.mark.unit
def test_property_data_bad():
    inputs = map(StringIO, bad_data_csv)
    for inpfile in inputs:
        with pytest.raises(Exception):
            PropertyData.from_csv(inpfile, 0)
