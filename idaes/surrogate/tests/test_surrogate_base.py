"""
Tests for surrogate base class and functions.
"""
# standard library
import json
import os
from pathlib import Path
import tempfile

# third-party packages
import pytest

# this package
from idaes.surrogate import surrogate_base

ALAMO_VALUE = 1
PYSMO_VALUE = 2.3
CONFIG = {
    surrogate_base.Config.GLOBAL_SECTION: {"toast": "burnt"},
    "alamo": {"value": ALAMO_VALUE},
    "pysmo": {"value": PYSMO_VALUE},
}


@pytest.fixture
def sample_config():
    # make temporary file
    tmpf = tempfile.NamedTemporaryFile(prefix="idaes_", suffix=".json", mode="w+")
    # dump config into file
    json.dump(CONFIG, tmpf)
    tmpf.flush()
    # reset position to start of file
    tmpf.seek(0)
    yield tmpf
    # clean up
    tmpf.close()


def test_config_class(sample_config):
    configs = [
        surrogate_base.Config(input=sample_config),
        surrogate_base.Config(input=sample_config.name),
        surrogate_base.Config(input=Path(sample_config.name)),
        surrogate_base.Config(values=CONFIG),
    ]
    for cfg in configs:
        t = cfg.get_tool("alamo")
        assert t["value"] == ALAMO_VALUE
        t = cfg.get_tool("pysmo")
        assert t["value"] == PYSMO_VALUE
        g = cfg.get_global()
        assert g["toast"] == "burnt"
    # empty
    empty_config = surrogate_base.Config()
    assert empty_config.get_global() == {}
    pytest.raises(KeyError, empty_config.get_tool, "alamo")
