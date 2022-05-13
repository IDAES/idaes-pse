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
Tests for Python code style.
"""
import logging
import os
from pathlib import Path
import subprocess
import pytest

_log = logging.getLogger(__name__)


# The most stylish dirs in the project
DIRS = [
    str(p)
    for p in (
        Path("idaes.core.dmf"),
        # Path("apps/ddm-learning/alamo_python/alamopy"),
        # Path("apps/ddm-learning/ripe_python/ripe"),
    )
]


STYLE_CHECK_CMD = "flake8"


@pytest.mark.unit
def test_flake8():
    cwd = os.getcwd()
    for d in DIRS:
        path = os.path.join(cwd, d)
        if not os.path.exists(path):
            _log.warning(
                f"Target path '{d}' not found in current dir, '{cwd}'. " "Skipping test"
            )
            continue
        if not os.path.isdir(path):
            _log.warning(
                f"Target path '{d}' in current dir, '{cwd}', is not a directory. "
                "Skipping test"
            )
            continue
        cmd = [STYLE_CHECK_CMD, d]
        _log.info(f"Test code style with command '{' '.join(cmd)}'")
        try:
            proc = subprocess.Popen(cmd)
        except FileNotFoundError:
            _log.warning(
                f"Style checker {STYLE_CHECK_CMD} not found. Skipping style tests"
            )
            break
        proc.wait()
        status = proc.returncode
        assert status == 0, f"Style checker '{STYLE_CHECK_CMD}' had errors for {path}"
