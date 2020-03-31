##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for idaes.commands
"""
# stdlib
import os
import pathlib

# third-party
import pytest

# package
from idaes.commands import examples
from idaes.util.system import TemporaryDirectory

################
# get-examples #
################


def test_examples_list_releases():
    releases = examples.get_releases(True)
    assert len(releases) > 0
    examples.print_releases(releases, True)


def test_examples_download_bad_version():
    releases = [examples.Release("baddate", "badtag", "info")]
    assert pytest.raises(
        examples.DownloadError, examples.download, releases, "", "1.2.3", True
    )


def test_examples_find_python_directories():
    with TemporaryDirectory() as tmpd:
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


def test_illegal_dirs():
    releases = [examples.Release("baddate", "1.2.3", "info")]
    with TemporaryDirectory() as tmpd:
        root = pathlib.Path(tmpd)
        # git
        (root / ".git").mkdir()
        try:
            examples.download(releases, root, "1.2.3", True)
        except examples.DownloadError as err:
            assert ".git" in str(err)
        (root / ".git").rmdir()
