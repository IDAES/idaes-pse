#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Test doc build.

This does *not* try to build the docs. Instead it looks
for the log of the build and looks for errors in it.
"""
# stdlib
import logging
import os
from pathlib import Path
from subprocess import Popen

# third-party
import pytest


_log = logging.getLogger(__name__)


@pytest.fixture
def docs_path():
    """Find the docs.

    We have to assume that the docs are somewhere
    relative to either this file or the user's current
    working directory.
    Start with the working directory.

    Returns:
        str: Doc path or None if not found.
    """
    cwd = os.path.realpath(os.getcwd())
    # look for "./docs"
    path = os.path.join(cwd, "docs")
    if os.path.exists(path):
        return path
    # look for (last) "docs/" in full current working dir
    cwd_comp = cwd.split(os.path.sep)
    comp_cnt = len(cwd_comp)
    try:
        idx = list(reversed(cwd_comp)).index("docs")
        return os.path.sep.join(cwd_comp[: comp_cnt - idx])
    except ValueError:
        pass
    # look for docs relative to this file
    try:
        path = os.path.join(os.path.dirname(__file__), "..", "..", "docs")
        # look for {this file}/../../docs
        if os.path.exists(path):
            return path
        # look for "docs/" in file's full path
        file_comp = path.split(os.path.sep)
        comp_cnt = len(file_comp)
        try:
            idx = list(reversed(file_comp)).index("docs")
            return os.path.sep.join(file_comp[: comp_cnt - idx])
        except ValueError:
            pass
    except NameError:
        pass  # __file__ not defined(?)


ERRLOG = "sphinx-errors.txt"


@pytest.mark.unit
def test_sphinx_build_log(docs_path):
    """Check the sphinx log for errors or certain warnings."""

    # Find the error log
    _log.info('docs path = "{}"'.format(docs_path))
    if docs_path is None:
        _log.warning('Could not find "docs" directory')
        return
    log_path = os.path.join(docs_path, ERRLOG)
    if not os.path.exists(log_path):
        _log.warning(
            'Could not find "{}" in docs directory: {}'.format(ERRLOG, log_path)
        )
        return

    # Check error log contents; ignore some warnings
    real_errors = []
    with open(log_path) as f:
        for line in f:
            if not (_is_sphinx_warning(line) and _can_ignore_sphinx_warning(line)):
                if not real_errors:  # first error
                    print("+---> Sphinx errors")
                real_errors.append(line)
                print(f"| {line}", end="")
    if real_errors:
        print("+---< Sphinx errors")
        assert False, f"{len(real_errors)} errors found"


def _is_sphinx_warning(line: str) -> bool:
    return "WARNING" in line


def _can_ignore_sphinx_warning(line: str) -> bool:
    if "duplicate object description" in line:
        return True
    return False


def _have_sphinx():
    """Test if a working 'sphinx-build' command exists."""
    have_sphinx = True
    try:
        Popen(["sphinx-build", "--version"]).wait()
    except:
        have_sphinx = False
    return have_sphinx


@pytest.mark.usefixtures("run_in_tmp_path")
@pytest.mark.component
def test_doctests(docs_path):
    if _have_sphinx():
        build_path = os.path.join(docs_path, "build")
        command = ["sphinx-build", "-M", "doctest", docs_path, build_path]
        proc = Popen(command)
        proc.wait(600)
        assert proc.returncode == 0
        _cleanup_generated_files(docs_path)


def _cleanup_generated_files(docs_path):
    docs = Path(docs_path)

    # remove selected autodoc2 Markdown docs in reference
    reference = docs / "reference_guides"
    for pkg in ("core.util.structfs",):
        print(f"Clean up autodoc2 files for {pkg}")
        path = reference
        for subpkg in pkg.split("."):
            path = path / subpkg
        for filename in path.glob("*.md"):
            print(f"- remove: {path / filename}")
            os.unlink(path / filename)
