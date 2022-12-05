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
Tests for idaes.core.dmf.datasets module
"""
import copy
import json
import pytest
import uuid
from idaes.core.dmf import datasets


# Constants
# ---------

DMF_CONF = {
    "_id": "5648dd0b011f42d2801381ea698df9d9",
    "created": "2000-01-01T00:00:00.123456",
    "description": "Test",
    "htmldocs": [
        "https://idaes-pse.readthedocs.io/en/stable",
        "{dmf_root}/docs/build/html",
    ],
    "modified": "2000-01-01T00:00:01.000000",
    "name": "test_data",
}

DS_CONF = {
    "name": "Test",
    "text": {
        "file": "Test.txt",
        "title": "Test",
        "date": "1970",
        "authors": "Dan Gunter",
        "venue": "Nowhere",
        "doi": "https://doi.org/10.1234/5.6789012",
    },
    "tables": [
        {"name": "example", "description": "An example", "datafile": "example.csv"}
    ],
}
DS_CONF_NO_TABLES = DS_CONF.copy()
DS_CONF_NO_TABLES["name"] = "Test-NT"
DS_CONF_NO_TABLES["tables"] = []
DS_TXT = "This is not a publication."
DS_CSV = "Column 1 [T],Column 2 [],Column3\n1,2,3\n"

# Fixtures
# --------


@pytest.fixture(scope="session")
def dmf_publication_dataset(tmp_path_factory):
    """Create and opulate a temporary directory with an example publication and
       an associated table (CSV).

    Returns:
        (Path) The path to the created directory.
    """
    p = tmp_path_factory.mktemp("ds_conf")
    # Write configuration file
    with (p / datasets.Dataset.CONF_NAME).open("w") as f:
        json.dump(DS_CONF, f)
    # Write publication and data file
    (p / "Test.txt").open("w").write(DS_TXT)
    (p / "example.csv").open("w").write(DS_CSV)
    return p


@pytest.fixture(scope="session")
def dmf_publication_dataset_no_tables(tmp_path_factory):
    p = tmp_path_factory.mktemp("ds_conf_notables")
    # Write configuration file
    with (p / datasets.Dataset.CONF_NAME).open("w") as f:
        json.dump(DS_CONF_NO_TABLES, f)
    # Write publication and data file
    (p / "Test.txt").open("w").write(DS_TXT)

    return p


@pytest.fixture(scope="session")
def dmf_workspace_path(tmp_path_factory):
    """Create a temporary DMF workspace directory and populate it with a
    configuration file.
    """
    p = tmp_path_factory.mktemp("ds_workspace")
    with (p / "config.yaml").open("w") as f:
        json.dump(DMF_CONF, f)
    return p


@pytest.fixture(scope="session")
def pub_datasets(
    dmf_workspace_path, dmf_publication_dataset, dmf_publication_dataset_no_tables
):
    from idaes.core.dmf import datasets

    ds = datasets.PublicationDataset(dmf_workspace_path)
    ds.load(dmf_publication_dataset)
    ds.load(dmf_publication_dataset_no_tables)
    yield ds


@pytest.fixture
def save_restore_config():
    """Use this to avoid changes to global config.yaml"""
    ws = datasets.get_dataset_workspace()
    config = ws / "config.yaml"
    # copy current config value
    with config.open("r") as fp:
        data = fp.read()
    yield config
    # put it back
    with config.open("w") as fp:
        fp.write(data)


# Tests
# -----


@pytest.mark.unit
def test_get_dataset_workspace():
    w = datasets.get_dataset_workspace()
    assert w
    assert w.exists()


@pytest.mark.unit
def test_dataset_init(dmf_workspace_path):
    ds = datasets.Dataset(dmf_workspace_path)


@pytest.mark.unit
def test_publication_dataset_init(dmf_workspace_path):
    ds = datasets.PublicationDataset(dmf_workspace_path)


@pytest.mark.unit
def test_publication_dataset_load_and_retrieve(
    dmf_workspace_path, dmf_publication_dataset
):
    ds = datasets.PublicationDataset(dmf_workspace_path)
    ds.load(dmf_publication_dataset)
    pub, tables = ds.retrieve("Test")

    assert pub
    assert tables
    assert pub.name == "Test"

    assert len(tables) == 1
    x = DS_CONF["tables"][0]
    t = tables[x["name"]]
    df = t.data

    # check metadata
    assert t.description == x["description"]

    # check header
    # i. header names
    dx_col = DS_CSV.split("\n")[0].split(",")
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
            unit_str = col[unit_start + 1 : unit_end]
        else:
            unit_str = ""
        assert t.units_list[i] == unit_str

    # check data
    dx_val = [int(z) for z in DS_CSV.split("\n")[1].split(",")]
    df_val = df.iloc[0].values
    for i, v in enumerate(dx_val):
        assert df_val[i] == v


@pytest.mark.unit
def test_dataset_init_no_workspace(save_restore_config):
    ds = datasets.Dataset()


@pytest.mark.unit
def test_base_load_conf(dmf_workspace_path, tmp_path):
    ds = datasets.Dataset(dmf_workspace_path)
    # impossible directory
    with pytest.raises(datasets.FileMissingError):
        ds._load_conf(f"/{uuid.uuid4()}")
    # missing conf in path
    with pytest.raises(datasets.FileMissingError):
        ds._load_conf(tmp_path)
    # unparseable conf
    conf = tmp_path / datasets.Dataset.CONF_NAME
    conf.open("w").write("Not valid JSON\n")
    with pytest.raises(datasets.ConfigurationError):
        ds._load_conf(tmp_path)


@pytest.mark.unit
def test_retrieve(pub_datasets):
    # does not exist
    with pytest.raises(KeyError):
        pub_datasets.retrieve("foobar")
    # does not have any tables: not an error
    nt_pub = DS_CONF_NO_TABLES["name"]
    pub, tables = pub_datasets.retrieve(nt_pub)
    assert pub
    assert len(tables) == 0


@pytest.mark.unit
def test_load_publication_missing_keys(dmf_workspace_path, tmp_path_factory):
    from idaes.core.dmf import datasets

    # Try to load bad config files
    for key in (
        datasets.PublicationDataset.NAME_KEY,
        datasets.PublicationDataset.PUBLICATION_KEY,
        "file",
    ):
        p = tmp_path_factory.mktemp(f"ds_conf_missing_{key}")
        # Copy config and delete required key
        bad_conf = copy.deepcopy(DS_CONF)
        if key == "file":
            del bad_conf[datasets.PublicationDataset.PUBLICATION_KEY]["file"]
        else:
            del bad_conf[key]
        # Write configuration file
        with (p / datasets.Dataset.CONF_NAME).open("w") as f:
            json.dump(bad_conf, f)
        # load the dataset
        ds = datasets.PublicationDataset(dmf_workspace_path)
        with pytest.raises(datasets.ConfigurationError):
            ds.load(p)


@pytest.mark.unit
def test_load_publication_isbn(dmf_workspace_path, tmp_path_factory):
    from idaes.core.dmf import datasets

    # Try to load config files with isbn
    p = tmp_path_factory.mktemp(f"ds_conf_has_isbn")
    # Copy config and delete required key
    conf = DS_CONF.copy()
    conf[datasets.PublicationDataset.PUBLICATION_KEY]["isbn"] = "ISBN-1"
    conf["tables"] = []
    # Write configuration file
    with (p / datasets.Dataset.CONF_NAME).open("w") as f:
        json.dump(conf, f)
    (p / "Test.txt").open("w").write(DS_TXT)
    # load the dataset
    ds = datasets.PublicationDataset(dmf_workspace_path)
    print(f"@@ load pub from conf: {conf}")
    ds.load(p)
