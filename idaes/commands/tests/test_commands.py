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
import json
import os
from pathlib import Path
import subprocess
import uuid

# third-party
import pytest

# package
from idaes.commands import examples
from idaes.util.system import TemporaryDirectory

################
# get-examples #
################


@pytest.mark.nocircleci()  # goes out to network
def test_examples_n():
    target_dir = str(uuid.uuid4())  # pick something that won't exist
    retcode = subprocess.call(["idaes", "get-examples", "-N", "-d", target_dir])
    assert retcode == 255  # result of sys.exit(-1)


@pytest.mark.nocircleci()  # goes out to network
def test_examples_list_releases():
    releases = examples.get_releases(True)
    assert len(releases) > 0
    examples.print_releases(releases, True)


@pytest.mark.nocircleci()  # goes out to network
def test_examples_download_bad_version():
    releases = [examples.Release("baddate", "badtag", "info")]
    assert pytest.raises(examples.DownloadError, examples.download, Path("."), "1.2.3")


def test_examples_find_python_directories():
    with TemporaryDirectory() as tmpd:
        root = Path(tmpd)
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
            path_i = Path(i)
            assert path_i in rel_found_dirs
            for j in ("c", "d"):
                assert path_i / j in rel_found_dirs


def test_illegal_dirs():
    with TemporaryDirectory() as tmpd:
        root = Path(tmpd)
        # git
        (root / ".git").mkdir()
        try:
            examples.download(root, "1.2.3")
        except examples.DownloadError as err:
            assert "exists" in str(err)
        finally:
            (root / ".git").rmdir()


test_cell_nb = {
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Hello, world\n"
     ]
    }
   ],
   "source": [
    "print(\"Hello, world\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {
    "tags": [
     "test",
     "remove_cell"
    ]
   },
   "outputs": [],
   "source": [
    "assert 2 + 2 == 4"
   ]
  }
 ],
 "metadata": {
  "celltoolbar": "Tags",
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}


@pytest.fixture
def remove_cells_notebooks(tmp_path):
    # note: first entry must be in root, one entry per subdir, depth search order
    notebooks = [tmp_path / d for d in (
        "notebook_testing.ipynb",
        Path("a") / "notebook_testing.ipynb",
        Path("a") / "b" / "notebook_testing.ipynb")]
    first = True
    for nb in notebooks:
        if not first:
            nb.parent.mkdir()
        json.dump(test_cell_nb, nb.open("w"))
        first = False
    yield notebooks
    _remove(tmp_path, notebooks)


@pytest.fixture
def remove_cells_notebooks_nosuffix(tmp_path):
    # note: first entry must be in root, one entry per subdir, depth search order
    notebooks = [tmp_path / d for d in (
        "notebook.ipynb",
        Path("a") / "notebook.ipynb",
        Path("a") / "b" / "notebook.ipynb")]
    first = True
    for nb in notebooks:
        if not first:
            nb.parent.mkdir()
        json.dump(test_cell_nb, nb.open("w"))
        first = False
    yield notebooks
    _remove(tmp_path, notebooks)


def _remove(tmp_path, notebooks):
    for f in tmp_path.rglob("*.ipynb"):
        f.unlink()
    removed = set()
    for nb in reversed(notebooks):
        if nb.parent not in removed:
            nb.parent.rmdir()
            removed.add(nb.parent)


def test_find_notebook_files(remove_cells_notebooks):
    root = remove_cells_notebooks[0].parent
    nbfiles = examples.find_notebook_files(root)
    assert(len(nbfiles) == len(remove_cells_notebooks))
    for nbfile in nbfiles:
        assert nbfile in remove_cells_notebooks


def test_strip_test_cells(remove_cells_notebooks):
    root = remove_cells_notebooks[0].parent
    examples.strip_test_cells(root)
    nbfiles = examples.find_notebook_files(root)
    assert(len(nbfiles) == 2 * len(remove_cells_notebooks))
    for nbfile in nbfiles:
        with nbfile.open("r") as f:
            nbdata = json.load(f)
            if not nbfile.stem.endswith("_testing"):
                for cell in nbdata["cells"]:
                    cell_meta = cell["metadata"]
                    tags = cell_meta.get("tags", [])
                    assert examples.REMOVE_CELL_TAG not in tags  # tag is gone
            else:
                n = 0
                for cell in nbdata["cells"]:
                    cell_meta = cell["metadata"]
                    tags = cell_meta.get("tags", [])
                    if examples.REMOVE_CELL_TAG in tags:
                        n += 1
                assert n > 0  # tag still there


def test_strip_test_cells_nosuffix(remove_cells_notebooks_nosuffix):
    root = remove_cells_notebooks_nosuffix[0].parent
    examples.strip_test_cells(root)
    nbfiles = examples.find_notebook_files(root)
    assert(len(nbfiles) == 2 * len(remove_cells_notebooks_nosuffix))
    for nbfile in nbfiles:
        with nbfile.open("r") as f:
            nbdata = json.load(f)
            if nbfile.stem.endswith(examples.STRIPPED_NOTEBOOK_SUFFIX):
                for cell in nbdata["cells"]:
                    cell_meta = cell["metadata"]
                    tags = cell_meta.get("tags", [])
                    assert examples.REMOVE_CELL_TAG not in tags  # tag is gone
            else:
                n = 0
                for cell in nbdata["cells"]:
                    cell_meta = cell["metadata"]
                    tags = cell_meta.get("tags", [])
                    if examples.REMOVE_CELL_TAG in tags:
                        n += 1
                assert n > 0  # tag still there
