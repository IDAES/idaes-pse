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
import logging
import os
from pathlib import Path
from shutil import rmtree


def create_module_scratch(module_name):
    """Create a scratch that is easily found (by developers) and removed, and doesn't
    clutter up the working directories.
    """
    path = Path("~").expanduser()
    for path_part in ".idaes", "_scratch", module_name:
        path = path / path_part
        if not path.exists():
            path.mkdir()
    return path


def rmtree_scratch(scratch_path):
    """Do our level best to remove all the temporary files."""
    try:
        rmtree(scratch_path, ignore_errors=True)
    except Exception as err:
        pass
