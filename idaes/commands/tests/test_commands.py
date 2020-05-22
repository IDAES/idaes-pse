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
import shutil
import subprocess
import uuid

# third-party
from click.testing import CliRunner
import pytest

# package
from idaes.commands import examples
from idaes.util.system import TemporaryDirectory


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


@pytest.fixture
def random_tempdir():
    """Make a completely cross-platform random temporary directory in
    the current working directory, and yield it. As cleanup, recursively
    remove all contents of this temporary directory.
    """
    origdir = os.getcwd()
    random_name = str(uuid.uuid4())
    os.mkdir(random_name)
    tempdir = Path(random_name)
    yield tempdir
    os.chdir(origdir)
    shutil.rmtree(tempdir)


################
# get-examples #
################


def test_examples_cli_noop(runner):
    result = runner.invoke(examples.get_examples, ["--no-install", "--no-download"])
    assert result.exit_code == 0


@pytest.mark.nocircleci()
def test_examples_cli_list(runner):
    result = runner.invoke(examples.get_examples, ["-l"])
    assert result.exit_code == 0


@pytest.mark.nocircleci()
def test_examples_cli_download(runner, random_tempdir):
    # failure with existing dir
    result = runner.invoke(examples.get_examples, ["-d", str(random_tempdir), "-I"])
    assert result.exit_code == -1
    # pick subdir, should be ok; use old version we know exists
    dirname = str(random_tempdir / "examples")
    result = runner.invoke(examples.get_examples, ["-d", dirname, "-I", "-V", "1.5.0"])
    assert result.exit_code == 0


@pytest.mark.nocircleci()
def test_examples_cli_default_version(runner, random_tempdir):
    dirname = str(random_tempdir / "examples")
    result = runner.invoke(examples.get_examples, ["-d", dirname])
    assert result.exit_code == -1


@pytest.mark.nocircleci()
def test_examples_cli_download_unstable(runner, random_tempdir):
    dirname = str(random_tempdir / "examples")
    # unstable version but no --unstable flag
    result = runner.invoke(examples.get_examples, ["-d", dirname, "-V", "1.2.3-beta"])
    assert result.exit_code == -1


@pytest.mark.nocircleci()
def test_examples_cli_copy(runner, random_tempdir):
    dirname = str(random_tempdir / "examples")
    # local path is bad
    badpath = str(random_tempdir / "no-such-dir")
    result = runner.invoke(examples.get_examples, ["-d", dirname, "--local", badpath])
    assert result.exit_code == -1
    # local dir exists, no REPO_DIR in it
    src_dir = random_tempdir / "examples-dev"
    src_dir.mkdir()
    result = runner.invoke(
        examples.get_examples, ["-d", dirname, "--local", str(src_dir), "-I"]
    )
    assert result.exit_code == -1
    # local path ok, target is ok
    # src_dir = random_tempdir / "examples-dev"
    # src_dir.mkdir()
    (src_dir / examples.REPO_DIR).mkdir()
    (src_dir / examples.REPO_DIR / "file.py").open("w")
    result = runner.invoke(
        examples.get_examples, ["-d", dirname, "--local", str(src_dir), "-I"]
    )
    assert result.exit_code == 0
    # make sure the file was copied into the target dir
    assert (Path(dirname) / "file.py").exists()


# non-CLI


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
    assert pytest.raises(examples.DownloadError, examples.download, Path("."), "1.2.3")


@pytest.mark.nocircleci()
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


@pytest.mark.nocircleci()
def test_examples_check_github_response():
    # ok result
    examples.check_github_response([{"result": "ok"}], {})
    # rate limit
    pytest.raises(
        examples.GithubError,
        examples.check_github_response,
        {"message": "API rate limit exceeded,dude"},
        {"X-RateLimit-Reset": "2000000000"},
    )
    # unknown error
    pytest.raises(
        examples.GithubError,
        examples.check_github_response,
        {"message": "Something is seriously wrong"},
        {},
    )
    # unknown, no message - still an error
    pytest.raises(
        examples.GithubError,
        examples.check_github_response,
        {"oscar": "green grouch"},
        {},
    )
    # non-dict: error
    pytest.raises(
        examples.GithubError, examples.check_github_response, 12, "hello",
    )


def test_examples_install_src(random_tempdir):
    # monkey patch a random install package name so as not to have
    # any weird side-effects on other tests
    orig_install_pkg = examples.INSTALL_PKG
    examples.INSTALL_PKG = "i" + str(uuid.uuid4()).replace("-", "_")
    # print(f"1. curdir={os.curdir}")
    # create fake package
    src_dir = random_tempdir / "src"
    src_dir.mkdir()
    m1_dir = src_dir / "module1"
    m1_dir.mkdir()
    (m1_dir / "groot.py").open("w").write("print('I am groot')\n")
    m2_dir = m1_dir / "module1_1"
    m2_dir.mkdir()
    (m2_dir / "groot.py").open("w").write("print('I am groot')\n")
    # install it
    examples.install_src("0.0.0", src_dir)
    # make one of the directories unwritable
    # print(f"2. curdir={os.curdir}")
    m1_dir.chmod(0o400)
    pytest.raises(examples.InstallError, examples.install_src, "0.0.0", src_dir)
    # change it back so we can remove this temp dir
    m1_dir.chmod(0o700)
    # also patch back the proper install package name
    examples.INSTALL_PKG = orig_install_pkg


def test_examples_cleanup(random_tempdir):
    tempdir = random_tempdir
    # put some crap in the temporary dir
    #
    # <tempdir>/
    #   idaes_examples.egg-info/
    #      afile
    #   dist/
    #      idaes_examples-1.0.0-py3.8.egg
    #      idaes_examples-2.0.0-py3.9.egg
    #   <tempdir2>/
    #       afile
    #
    os.chdir(tempdir)
    eggy = Path("idaes_examples.egg-info")
    eggy.mkdir()
    (eggy / "afile").open("w")
    distdir = Path("dist")
    distdir.mkdir()
    (distdir / "idaes_examples-1.0.0-py3.8.egg").open("w")
    (distdir / "idaes_examples-2.0.0-py3.9.egg").open("w")
    name2 = str(uuid.uuid4())
    tempsubdir = Path(name2)
    tempsubdir.mkdir()
    (tempsubdir / "afile").open("w")
    # All of the above should be removed
    # by the cleanup function.
    examples.g_tempdir = tempsubdir
    examples.g_egg = eggy
    examples.clean_up_temporary_files()
    # Check that everything is removed
    assert not distdir.exists()
    assert not eggy.exists()
    assert not tempsubdir.exists()


def test_examples_cleanup_nodist(random_tempdir):
    tempdir = random_tempdir
    # put some crap in the temporary dir
    #
    # <tempdir>/
    #   idaes_examples.egg-info/
    #      afile
    #   <tempdir2>/
    #       afile
    #
    os.chdir(tempdir)
    eggy = Path("idaes_examples.egg-info")
    eggy.mkdir()
    (eggy / "afile").open("w")
    name2 = str(uuid.uuid4())
    tempsubdir = Path(name2)
    tempsubdir.mkdir()
    (tempsubdir / "afile").open("w")
    # All of the above should be removed
    # by the cleanup function.
    examples.g_tempdir = tempsubdir
    examples.g_egg = eggy
    examples.clean_up_temporary_files()
    # Check that everything is removed
    assert not eggy.exists()
    assert not tempsubdir.exists()


def test_examples_cleanup_nodist_noegg(random_tempdir):
    tempdir = random_tempdir
    # put some crap in the temporary dir
    #
    # <tempdir>/
    #   <tempdir2>/
    #       afile
    #
    os.chdir(tempdir)
    name2 = str(uuid.uuid4())
    tempsubdir = Path(name2)
    tempsubdir.mkdir()
    (tempsubdir / "afile").open("w")
    # All of the above should be removed
    # by the cleanup function.
    examples.g_tempdir = tempsubdir
    examples.clean_up_temporary_files()
    # Check that everything is removed
    assert not tempsubdir.exists()


def test_examples_cleanup_nothing(random_tempdir):
    tempdir = random_tempdir
    # nothing to remove, should still be ok
    os.chdir(tempdir)
    examples.clean_up_temporary_files()
    # if we set some globals to bogus values, still OK
    examples.g_egg = Path("no-such-file.egg-info")
    examples.g_tempdir = Path("no-such-file-tempdir")
    examples.clean_up_temporary_files()
    # if we create files and set perms to 000, still OK
    eggy = Path("egg")
    eggy.mkdir()
    eggy.chmod(0)
    examples.g_egg = eggy
    subdir = Path("subdirinho")
    subdir.mkdir()
    subdir.chmod(0)
    examples.g_tempdir = subdir
    dist = Path("dist")
    dist.mkdir()
    dist.chmod(0)
    examples.clean_up_temporary_files()
    # random_tempdir will clean up these files:
    # dist.chmod(700) - removed with .rmdir() which works regardless!
    eggy.chmod(0o700)
    subdir.chmod(0o700)


def test_examples_local(random_tempdir):
    d = random_tempdir
    tgt = d / "examples"
    src = d / "src"
    root = d
    # source doesn't exist yet
    pytest.raises(examples.CopyError, examples.copy_contents, tgt, root)
    # source isn't a directory
    src.open("w")
    pytest.raises(examples.CopyError, examples.copy_contents, tgt, root)
    # source exists = OK
    src.unlink()
    src.mkdir()
    (src / "a").open("w")  # create a file
    examples.copy_contents(tgt, root)
    assert (tgt / "a").exists()
    # target now exists, not ok
    pytest.raises(examples.CopyError, examples.copy_contents, tgt, root)
    # done


@pytest.mark.nocircleci()
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


@pytest.mark.nocircleci()
def test_get_examples_version():
    assert examples.get_examples_version("1.5.0") == "1.5.1"
    assert examples.get_examples_version("foo") == None


test_cell_nb = {
    "cells": [
        {
            "cell_type": "code",
            "execution_count": 2,
            "metadata": {},
            "outputs": [
                {"name": "stdout", "output_type": "stream", "text": ["Hello, world\n"]}
            ],
            "source": ['print("Hello, world")'],
        },
        {
            "cell_type": "code",
            "execution_count": 3,
            "metadata": {"tags": ["test", "remove_cell"]},
            "outputs": [],
            "source": ["assert 2 + 2 == 4"],
        },
    ],
    "metadata": {
        "celltoolbar": "Tags",
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {
            "codemirror_mode": {"name": "ipython", "version": 3},
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "nbconvert_exporter": "python",
            "pygments_lexer": "ipython3",
            "version": "3.8.1",
        },
    },
    "nbformat": 4,
    "nbformat_minor": 2,
}


@pytest.fixture
def remove_cells_notebooks(tmp_path):
    # note: first entry must be in root, one entry per subdir, depth search order
    notebooks = [
        tmp_path / d
        for d in (
            "notebook_testing.ipynb",
            Path("a") / "notebook_testing.ipynb",
            Path("a") / "b" / "notebook_testing.ipynb",
        )
    ]
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
    notebooks = [
        tmp_path / d
        for d in (
            "notebook.ipynb",
            Path("a") / "notebook.ipynb",
            Path("a") / "b" / "notebook.ipynb",
        )
    ]
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
    assert len(nbfiles) == len(remove_cells_notebooks)
    for nbfile in nbfiles:
        assert nbfile in remove_cells_notebooks


def test_strip_test_cells(remove_cells_notebooks):
    root = remove_cells_notebooks[0].parent
    examples.strip_test_cells(root)
    nbfiles = examples.find_notebook_files(root)
    assert len(nbfiles) == 2 * len(remove_cells_notebooks)
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
    assert len(nbfiles) == 2 * len(remove_cells_notebooks_nosuffix)
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
