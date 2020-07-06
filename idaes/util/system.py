##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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


