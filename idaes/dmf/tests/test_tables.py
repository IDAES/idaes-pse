###############################################################################
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
###############################################################################
"""
Tests for idaes.dmf.tables
"""
__author__ = "Dan Gunter"

from pathlib import Path
import pytest
from idaes.dmf.tables import Table
from idaes.dmf import resource, DMF

# DATA

DATA_DIR = Path(__file__).parent / "data_files"

SIMPLE_CSV_DATA = [
    "C1 (m/s),C2[mol/L],C3 (m^3/kg)",
    "1,1.0,2.0",
    "2,2.0,4.0"
]


@pytest.fixture
def example_csv_file(tmp_path):
    p = tmp_path / "data.csv"
    with p.open("w", encoding="utf-8") as f:
        for line in SIMPLE_CSV_DATA:
            f.write(line)
            f.write("\n")
    return p


@pytest.fixture
def example_excel_file():
    p = DATA_DIR / "tables_data.xlsx"
    return p

@pytest.fixture
def tmp_dmf(tmp_path):
    p = tmp_path / "tables_dmf"
    p.mkdir()
    dmf = DMF(p, create=True)
    return dmf

# TESTS


@pytest.mark.unit
def test_table_create():
    t = Table()
    assert len(t.data) == 0
    assert len(t.units) == 0


@pytest.mark.unit
def test_table_csv(example_csv_file):
    t = Table()
    t.read_csv(example_csv_file)
    validate_example_data(t)


@pytest.mark.unit
def test_table_excel(example_excel_file):
    t = Table()
    t.read_excel(example_excel_file)
    validate_example_data(t)


@pytest.mark.unit
def test_table_excel_multisheet(example_excel_file):
    t = Table()
    with pytest.raises(ValueError):
        t.read_excel(example_excel_file, sheet_name=["Sheet1", "Sheet2"])


@pytest.mark.unit
def test_resource(example_csv_file):
    rsrc = resource.Resource()
    rsrc.v["desc"] = "Example table"
    rsrc.v["version_info"]["version"] = (1, 0, 1)

    t = Table()
    t.read_csv(example_csv_file)
    validate_example_data(t)
    t.add_to_resource(rsrc)

    t2 = Table.get_from_resource(rsrc)
    validate_example_data(t2)


@pytest.mark.component
def test_dmf(example_csv_file, tmp_dmf):
    # create and add to DMF
    rsrc = resource.Resource(type_=resource.ResourceTypes.tabular)
    rsrc.set_field("name", "table_test_resource")
    rsrc.set_field("desc", "Example table")

    t = Table()
    t.read_csv(example_csv_file)
    validate_example_data(t)
    t.add_to_resource(rsrc)

    tmp_dmf.add(rsrc)

    # retrieve from DMF
    r2 = tmp_dmf.find_one(name="table_test_resource")
    t2 = Table.get_from_resource(r2)
    # print(t2.data)
    validate_example_data(t2)


# Helpers

def validate_example_data(t):
    """Shared by test_table_*() functions, since their input data is the same.
    """
    input_data_columns = [f"C{num}" for num in range(1, 4)]
    assert list(t.data.columns) == input_data_columns
    input_data_units = ["m/s", "mol/L", "m^3/kg"]
    assert t.units == {input_data_columns[i]: input_data_units[i]
                       for i in range(len(input_data_columns))}


