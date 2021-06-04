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
System-related Python utility commands

- Replace tempfile functions with "safe" versions
  that stop and catch fire if the temporary directory
  is in the home directory. At least for testing, this
  seems like a good safety measure.
"""
import os
import tempfile

# Cache this, since the result of gettempdir() is also cached
_is_curdir = None


def fail_if_tempdir_is_curdir():
    global _is_curdir
    if _is_curdir is None:
        td = tempfile.gettempdir()  # initializes it 1st time
        try:
            curdir = os.getcwd()
        except (AttributeError, OSError):
            curdir = os.curdir
        _is_curdir = td == curdir
    if _is_curdir:
        raise RuntimeError("abort: temporary directory is going to be in "
                           "the current directory. Please set one of the "
                           "following environment variables to a suitable "
                           "directory: TMPDIR, TEMP, or TMP")

# Wrapped classes and functions


class TemporaryDirectory(tempfile.TemporaryDirectory):
    def __init__(self, *args, **kwargs):
        fail_if_tempdir_is_curdir()
        super().__init__(*args, **kwargs)


def NamedTemporaryFile(*args, **kwargs):
    fail_if_tempdir_is_curdir()
    return tempfile.NamedTemporaryFile(*args, **kwargs)


def mkdtemp(*args, **kwargs):
    fail_if_tempdir_is_curdir()
    return tempfile.mkdtemp(*args, **kwargs)


