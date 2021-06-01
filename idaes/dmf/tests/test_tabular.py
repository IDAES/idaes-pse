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
Test tabular data module.

Note: This is an attempt a ~100% test coverage of a module.
"""
# stdlib
import json
import logging
import os
from io import StringIO
import sys

# third-party
import pytest

# package-local
from idaes.dmf import tabular
from idaes.dmf import util, errors

# for testing
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)

F = tabular.Fields

# class: Column

data_column = {
    F.DATA_UNITS: "ms",
    F.DATA_VALUES: [1.0, 2.0, 3.0],
    F.DATA_ERRORS: [0.0, -1.0, -2.0],
    F.DATA_ERRTYPE: "whatever",
    "extra": "ignored",
}

column_name = "hello"


@pytest.fixture
def col():
    return tabular.Column(column_name, data_column)


@pytest.mark.unit
def test_column(col):
    assert col.name == column_name


@pytest.mark.unit
def test_column_data(col):
    d = col.data()
    for f in (F.DATA_UNITS, F.DATA_VALUES, F.DATA_ERRORS, F.DATA_ERRTYPE):
        assert d[f] == data_column[f]


@pytest.mark.unit
def test_column_len(col):
    assert len(col) == len(data_column[F.DATA_VALUES])


# class: Metadata


ex_metadata = {
    F.DTYPE: "foobar",
    F.AUTH: "Philip Roth",
    F.DATE: "1969-01-12",
    F.TITLE: "Portnoy's Complaint",
    F.INFO: "Random House",
}

ex_csv = "\n".join(
    [
        "# This is some metadata",
        "#",
        'Source,Philip Roth,"Portnoy\'s Complaint",Random House,1969-01-12',
        "Notes,Confessions of a lust-ridden, mother-addicted Jewish bachelor",
        "Datatype,foobar",
        "",
        "# All done",
    ]
)


@pytest.fixture
def metadata():
    return tabular.Metadata(ex_metadata)


@pytest.fixture
def csvfile():
    fake_file = StringIO(ex_csv)
    return fake_file


@pytest.mark.unit
def test_metadata(metadata):
    assert metadata.author == ex_metadata[F.AUTH]
    assert metadata.date == ex_metadata[F.DATE]


@pytest.mark.unit
def test_metadata_bad():
    inputs = map(StringIO, ["", "Who cares,Some junk", "source,"])
    for inpfile in inputs:
        with pytest.raises(ValueError):
            print('Bad metadata input: "{}"'.format(inpfile.getvalue()))
            tabular.Metadata.from_csv(inpfile)


@pytest.mark.unit
def test_metadata_good():
    inpfile = StringIO(ex_csv)
    obj = tabular.Metadata.from_csv(inpfile)
    assert obj.author != ""
    assert obj.title != ""
    assert obj.info != ""
    whole_len = len(obj.source)
    parts_len = map(len, [obj.author, obj.title, obj.info])
    assert whole_len > sum(parts_len)


@pytest.mark.unit
def test_from_csv(csvfile):
    m = tabular.Metadata.from_csv(csvfile)
    assert m.author == ex_metadata[F.AUTH]
    assert m.date == ex_metadata[F.DATE]


@pytest.mark.unit
def test_md_properties(metadata):
    assert metadata.author is not None
    assert metadata.date is not None
    assert metadata.title is not None
    assert metadata.info is not None
    # setter and getter
    metadata.datatype = "foo"
    assert metadata.datatype == "foo"


# class: TabularData


ex_data = [
    {
        "name": "Density Data",
        "units": "g/cm^3",
        "values": [1.0053, 1.0188, 1.0023],
        "errors": [0.00005, 0.0001, 0.00005],
        "error_type": "absolute",
    }
]

ex_data2 = [
    {
        "name": "Density Data",
        "units": "g/cm^3",
        "values": [1.0053, 1.0188, 1.0023],
        "errors": [0.00005, 0.0001, 0.00005],
        "error_type": "absolute",
    },
    {
        "name": "Density Data2",
        "units": "g/cm^3",
        "values": [1.0053, 1.0188, 1.0023],
        "errors": [0.00005, 0.0001, 0.00005],
        "error_type": "absolute",
    },
]


ex_data2_badlen = [
    {
        "name": "Density Data",
        "units": "g/cm^3",
        "values": [1.0053, 1.0188, 1.0023],
        "errors": [0.00005, 0.0001],
        "error_type": "absolute",
    },
    {
        "name": "Density Data2",
        "units": "g/cm^3",
        "values": [1.0053, 1.0188],
        "errors": [0.00005, 0.0001, 0.00005],
        "error_type": "absolute",
    },
]

ex_data2_csv = [
    "ID,V1 (g/cm^3),Absolute Error,V2 (m/s),Relative Error",
    "1,1.0053, 0, 1.0035, 0.001",
    "2,1.0188, 0, 1.0037, 0.001",
    "3,1.0023, 0, 1.0039, 0.001",
]


ex_data2_csv_short = ["ID", "1"]


@pytest.fixture
def tabdata():
    obj = tabular.TabularData(ex_data)
    return obj


@pytest.fixture
def tabdata2():
    obj = tabular.TabularData(ex_data2)
    return obj


@pytest.mark.unit
def test_init(tabdata, tabdata2):
    pass


@pytest.mark.unit
def test_init_badtype():
    with pytest.raises(TypeError):
        tabular.TabularData("foo")
    with pytest.raises(TypeError):
        tabular.TabularData({"foo": "bar"})


@pytest.mark.unit
def test_init_emptyvalue():
    with pytest.raises(ValueError):
        tabular.TabularData([])


@pytest.mark.unit
def test_init_badvalue():
    with pytest.raises(ValueError):
        tabular.TabularData([{"foo": "bar"}])


@pytest.mark.unit
def test_columns(tabdata):
    c = tabdata.columns
    assert c[0]["units"] == ex_data[0]["units"]


@pytest.mark.unit
def test_len(tabdata, tabdata2):
    assert len(tabdata) == 3
    assert tabdata.num_rows == len(tabdata)
    assert len(tabdata2) == 3
    assert tabdata2.num_rows == len(tabdata2)


@pytest.mark.unit
def test_names(tabdata):
    assert tabdata.names()[0] == ex_data[0]["name"]


@pytest.mark.unit
def test_num_columns(tabdata, tabdata2):
    assert tabdata.num_columns == 1
    assert tabdata2.num_columns == 2


@pytest.mark.unit
def test_get_column(tabdata):
    name = ex_data[0]["name"]
    assert tabdata.get_column(name).values == ex_data[0]["values"]


@pytest.mark.unit
def test_get_column_bad(tabdata):
    name = "NOPE!"
    with pytest.raises(KeyError):
        tabdata.get_column(name)


@pytest.mark.unit
def test_get_column_index(tabdata):
    name = ex_data[0]["name"]
    assert tabdata.get_column_index(name) == 0


@pytest.mark.unit
def test_get_column_index_err(tabdata):
    name = "-no-such-name-"
    with pytest.raises(KeyError):
        tabdata.get_column_index(name)


@pytest.mark.unit
def test_num_rows(tabdata):
    assert tabdata.num_rows == 3


@pytest.mark.unit
def test_as_dict(tabdata):
    assert tabdata.as_list() == ex_data


@pytest.mark.unit
def test_as_arr(tabdata2):
    values, errors = tabdata2.as_arr()
    assert values[0][0] == ex_data2[0][F.DATA_VALUES][0]
    assert values[0][1] == ex_data2[0][F.DATA_VALUES][1]
    assert values[1][0] == ex_data2[1][F.DATA_VALUES][0]
    assert values[1][1] == ex_data2[1][F.DATA_VALUES][1]


@pytest.mark.unit
def test_td_badlen():
    with pytest.raises(ValueError):
        tabular.TabularData(ex_data2_badlen)


@pytest.mark.unit
def test_values_dataframe(tabdata2):
    df = tabdata2.values_dataframe()
    name = tabdata2.names()[0]
    values = ex_data2[0][F.DATA_VALUES]
    for i in range(len(values)):
        assert df[name][i] == values[i]


@pytest.mark.unit
def test_errors_dataframe(tabdata2):
    df = tabdata2.errors_dataframe()
    name = tabdata2.names()[0]
    values = ex_data2[0][F.DATA_ERRORS]
    for i in range(len(values)):
        assert df[name][i] == values[i]


@pytest.mark.unit
def test_td_dataframe_nopandas(tabdata):
    tabular.pd = None
    with pytest.raises(ImportError):
        df = tabdata.values_dataframe()


@pytest.mark.unit
def test_td_from_csv():
    strfile = StringIO("\n".join(ex_data2_csv))
    tabular.TabularData.from_csv(strfile)


@pytest.mark.unit
def test_td_from_csv_badheader():
    csv = ex_data2_csv[:]
    csv[0] += ",Extra"
    strfile = StringIO("\n".join(csv))
    with pytest.raises(ValueError):
        tabular.TabularData.from_csv(strfile)


@pytest.mark.unit
def test_td_from_csv_badrow():
    csv = ex_data2_csv[:]
    csv[1] += ",1"
    strfile = StringIO("\n".join(csv))
    with pytest.raises(ValueError):
        tabular.TabularData.from_csv(strfile)


@pytest.mark.unit
def test_td_csv_short():
    infile = StringIO("\n".join(ex_data2_csv_short))
    with pytest.raises(ValueError):
        tabular.TabularData.from_csv(infile)


# class: Table


@pytest.fixture
def table(tabdata, metadata):
    return tabular.Table(data=tabdata, metadata=metadata)


@pytest.mark.unit
def test_init_table(table):
    pass


@pytest.mark.unit
def test_table_init_variations():
    tabular.Table(data=ex_data2, metadata=ex_metadata)
    tabular.Table(data=ex_data2, metadata=ex_metadata)
    tabular.Table(data=ex_data2, metadata=[ex_metadata])


@pytest.mark.unit
def test_table_init_bad(tabdata, metadata):
    with pytest.raises(TypeError):
        tbl = tabular.Table(data=None, metadata=metadata)
    with pytest.raises(ValueError):
        tbl = tabular.Table(data=[], metadata=metadata)
    with pytest.raises(TypeError):
        tbl = tabular.Table(data=tabdata, metadata="hello")


@pytest.mark.unit
def test_iter_table(table):
    for key, value in table:
        assert key in ["meta", "data"]
        assert isinstance(value, list)


@pytest.mark.unit
def test_table_attrs(table):
    assert table.data is not None
    assert table.metadata is not None


@pytest.mark.unit
def test_dump_table(table):
    fp = StringIO()
    table.dump(fp)
    assert len(fp.getvalue()) > 50


@pytest.mark.unit
def test_dumps_table(table):
    assert len(table.dumps()) > 50
    assert str(table) == table.dumps()


@pytest.mark.unit
def test_roundtrip_table(table):
    fp = StringIO()
    table.dump(fp)
    fp.seek(0)
    table2 = tabular.Table.load(fp)
    for m1, m2 in zip(table.metadata, table2.metadata):
        assert m1.as_dict() == m2.as_dict()


def _table_reload(fp):
    """Table reload and check logic."""
    fp.seek(0)
    with pytest.raises(errors.DataFormatError):
        _ = tabular.Table.load(fp)


@pytest.mark.unit
def test_load_table_missing(tabdata, metadata):
    # ValueError with high-level keys missing
    for delkey in F.DATA, F.META:
        fp = StringIO()
        tbl = tabular.Table(data=tabdata, metadata=metadata)
        tbl_json = tbl.as_dict()
        del tbl_json[delkey]
        json.dump(tbl_json, fp)
        _table_reload(fp)
    # ValueError with metadata keys missing
    for mdelkey in F.DTYPE, F.AUTH, F.INFO:
        fp = StringIO()
        tbl = tabular.Table(data=tabdata, metadata=metadata)
        del tbl.metadata[0]._meta[mdelkey]
        tbl.dump(fp)
        _table_reload(fp)


@pytest.mark.unit
def test_load_table_nometa(tabdata):
    fp = StringIO()
    tbl = tabular.Table(data=tabdata, metadata=[])
    tbl.dump(fp)
    fp.seek(0)
    tabular.Table.load(fp)


# class: TabularObject

# hmmm.. it's entirely abstract..
