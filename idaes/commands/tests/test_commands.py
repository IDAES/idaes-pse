"""
Tests for idaes.commands
"""
# stdlib
import os
import pathlib
import tempfile

# third-party
import pytest

# package
from idaes.commands import examples

################
# get-examples #
################


def test_examples_list_releases():
    releases = examples.get_releases(True)
    assert len(releases) > 0
    examples.print_releases(releases)


def test_examples_download_bad_version():
    releases = [examples.Release("baddate", "badtag", "info")]
    assert pytest.raises(
        examples.DownloadError, examples.download, releases, "", "1.2.3"
    )


def test_examples_download_target_dir_exists():
    releases = [examples.Release("baddate", "1.2.3", "info")]
    curpath = pathlib.Path(os.curdir)
    assert pytest.raises(
        examples.DownloadError, examples.download, releases, curpath, "1.2.3"
    )


def test_examples_find_python_directories():
    with tempfile.TemporaryDirectory() as tmpd:
        root = pathlib.Path(tmpd)
        # populate a/c/file.py, a/d/file.py, b/c/file.py, b/d/file.py
        for i in ("a", "b"):
            # create parent
            os.mkdir(str(root / i))
            for j in ("c", "d"):
                pydir = root / i / j
                # create child dir
                os.mkdir(str(pydir))
                # put python file in child
                (pydir / "file.py").open("w")
        # find
        found_dirs = examples.find_python_directories(root)
        # check
        rel_found_dirs = [d.relative_to(root) for d in found_dirs]
        for i in ("a", "b"):
            path_i = pathlib.Path(i)
            assert path_i in rel_found_dirs
            for j in ("c", "d"):
                assert path_i / j in rel_found_dirs
