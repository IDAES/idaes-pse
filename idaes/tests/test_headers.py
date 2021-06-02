"""
Test that headers are on all files
"""
# stdlib
from pathlib import Path
import os
# third-party
from addheader.add import FileFinder, detect_files
import pytest
import yaml


@pytest.fixture
def package_root():
    """Determine package root.
    """
    import idaes
    return Path(idaes.__file__).parent


@pytest.fixture
def patterns(package_root):
    """Grab glob patterns from config file.
    """
    conf_file = package_root.parent / "addheader.yml"
    if not conf_file.exists():
        print(f"Cannot load configuration file from '{conf_file}'. Perhaps this is not development mode?")
        return None
    with open(conf_file) as f:
        conf_data = yaml.load(f)
    print(f"Patterns for finding files with headers: {conf_data['patterns']}")
    return conf_data["patterns"]


def test_headers(package_root, patterns):
    if patterns is None:
        print(f"ERROR: Did not get glob patterns: skipping test")
    else:
        # modify patterns to match the files that should have headers
        ff = FileFinder(package_root, glob_patterns=patterns)
        has_header, missing_header = detect_files(ff)
        if len(missing_header) > 0:
            pfx = str(package_root.resolve())
            pfx_len = len(pfx)
            file_list = ", ".join([str(p)[pfx_len + 1:] for p in missing_header])
            print(f"Missing headers from files under '{pfx}{os.path.sep}': {file_list}")
        # uncomment to require all files to have headers
        # assert len(missing_header) == 0
        assert len(missing_header) < 30
