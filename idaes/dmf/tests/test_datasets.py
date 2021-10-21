"""
Tests for idaes.dmf.datasets module
"""
import json
import pytest

from idaes.dmf import datasets


# Constants
# ---------

TEST_CONF = {
    "name": "Test",
    "text": {
        "file": "Test.txt",
        "title": "Test",
        "date": "1970",
        "authors": "Dan Gunter",
        "venue": "Nowhere",
        "doi": "https://doi.org/10.1234/5.6789012"
    },
    "tables": [
      {
        "name": "trash",
        "description": "Garbage",
        "datafile": "trash.csv"
      }
    ]
}
TEST_TXT = "This is not a publication."
TRASH_CSV = "Column 1 [T],Column 2 [],Column3\n1,2,3\n"
TEST_DMF_CONF = {
    "_id": "5648dd0b011f42d2801381ea698df9d9",
    "created": '2000-01-01T00:00:00.123456',
    "description": "Test",
    "htmldocs": [
        "https://idaes-pse.readthedocs.io/en/stable",
        "{dmf_root}/docs/build/html"
    ],
    "modified": '2000-01-01T00:00:01.000000',
    "name": "test_data"
}

# Fixtures
# --------


@pytest.fixture(scope="session")
def pub_conf_path(tmp_path_factory):
    p = tmp_path_factory.mktemp("ds_conf")

    # Write configuration file
    with (p / datasets.Dataset.CONF_NAME).open("w") as f:
        json.dump(TEST_CONF, f)

    # Write publication and data file
    (p / "Test.txt").open("w").write(TEST_TXT)
    (p / "trash.csv").open("w").write(TRASH_CSV)

    return p


@pytest.fixture(scope="session")
def workspace_path(tmp_path_factory):
    p = tmp_path_factory.mktemp("ds_workspace")

    with (p / "config.yaml").open("w") as f:
        json.dump(TEST_DMF_CONF, f)

# Tests
# -----


@pytest.mark.unit
def test_get_dataset_workspace():
    w = datasets.get_dataset_workspace()
    assert w
    assert w.exists()


@pytest.mark.unit
def test_dataset_init(workspace_path):
    ds = datasets.Dataset(workspace_path)


@pytest.mark.unit
def test_publication_dataset_init(workspace_path):
    ds = datasets.PublicationDataset(workspace_path)


@pytest.mark.unit
def test_publication_dataset_load_and_retrieve(workspace_path, pub_conf_path):
    ds = datasets.PublicationDataset(workspace_path)
    ds.load(pub_conf_path)
    pub, tables = ds.retrieve("Test")

    assert pub
    assert tables
    assert pub.name == "Test"

    assert len(tables) == 1
    x = TEST_CONF["tables"][0]
    t = tables[x["name"]]
    df = t.data

    # check metadata
    assert t.description == x["description"]

    # check header
    # i. header names
    dx_col = TRASH_CSV.split("\n")[0].split(",")
    for i, col in enumerate(dx_col):
        unit_start = col.find("[")
        col_name = col if unit_start == -1 else col[:unit_start].strip()
        assert df.columns[i] == col_name
    # ii. header units
    for i, col in enumerate(dx_col):
        unit_start = col.find("[")
        if unit_start >= 0:
            unit_end = col.find("]", unit_start)
            assert unit_end >= 0
            unit_str = col[unit_start + 1:unit_end]
        else:
            unit_str = ""
        assert t.units_list[i] == unit_str

    # check data
    dx_val = [int(z) for z in TRASH_CSV.split("\n")[1].split(",")]
    df_val = df.iloc[0].values
    for i, v in enumerate(dx_val):
        assert df_val[i] == v
