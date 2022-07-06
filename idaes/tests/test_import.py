# coding: utf-8
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

# stdlib
import importlib
import pathlib
import os
import re
import time

# package
import idaes
from idaes.core.dmf.util import ColorTerm
import pytest

good_modname = re.compile(r"[a-zA-Z][a-zA-Z0-9_]*")


SKIP, OK, BAD = "skipped", "success", "failed"


_term = ColorTerm()


def print_path_status(path, status, msg=""):
    color = {SKIP: "", OK: _term.green, BAD: _term.red}[status]
    print(f"{color}{status.upper():8s}{_term.resetc} {path} {msg}")


def importr(root: pathlib.Path, max_sec=10):
    """Import r_ecursively from the given root path."""
    base, failures, total = root.parent, {}, 0
    # iterate over flattened list of all paths ending a Python source file
    for path in root.rglob("*.py"):
        # check that all path components are valid module names
        # - this is to skip directories like '.ipynb_checkpoints'
        bad = False
        for name in path.parts[:-1]:
            if name == os.path.sep:
                continue
            if (
                name.startswith("Workshop")  # skip workshop tutorials (1)
                or name.startswith("Module_")  # skip workshop tutorials (2)
                or name == "tests"  # skip all tests
                or not good_modname.match(name)  # no _bla
            ):
                bad = True

            if bad:
                break
        # stop if bad directory component or Python filename is invalid
        bad = bad or not good_modname.match(path.parts[-1][:-3])
        if bad:
            print_path_status(path, SKIP)
            continue
        # module is valid, try importing it now
        module_path = path.relative_to(base).with_suffix("")
        module_name = ".".join(module_path.parts)
        try:
            start = time.time()
            importlib.import_module(module_name)
            sec = time.time() - start
            if sec > max_sec:
                print_path_status(path, BAD, msg="time={sec:3.1f}s")
                raise ImportError(f"Import took too long ({sec:.1f}s)")
            print_path_status(path, OK)
        except ImportError as e:
            failures[module_name] = str(e)
            print_path_status(path, BAD, msg="import error")
        total += 1
    return failures, total


@pytest.mark.unit
def test_import():
    root_dir = pathlib.Path(idaes.__file__).parent
    failures, total = importr(root_dir)
    n = len(failures)
    assert n == 0, f"{n:d} failures in {total:d} tests: {failures}"
