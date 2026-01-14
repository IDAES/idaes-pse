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
Tests for logutil module
"""

# standard library
import logging
import time

# third-party
import pytest

# package
from .. import logutil

DEBUG = False


def numbytes(p):
    return p.stat().st_size


def dump_console(p):
    if DEBUG:
        with open(p) as f:
            print(f.read())


@pytest.mark.unit
def test_quiet(tmp_path):
    log = logging.getLogger("idaes.test_quiet")
    ni_log = logging.getLogger("ni.test_quiet")
    logfile = tmp_path / "test_quiet.log"
    handler = logging.FileHandler(logfile)
    log.addHandler(handler)
    ni_log.addHandler(handler)
    assert numbytes(logfile) == 0
    log.setLevel(logging.INFO)
    ni_log.setLevel(logging.INFO)

    # messages initially logged
    for i in range(2):
        log.info("this is a log message 1")
        time.sleep(0.1)
    dump_console(logfile)
    sz1 = numbytes(logfile)
    assert sz1 > 0

    logutil.quiet()

    # no more messages from idaes. logger
    for i in range(2):
        log.info("this is a log message 2")
        time.sleep(0.1)
    dump_console(logfile)
    sz2 = numbytes(logfile)
    assert sz2 == sz1

    # still get messages from non-idaes. logger
    for i in range(2):
        ni_log.info("this is a log message 3")
        time.sleep(0.1)
    dump_console(logfile)
    sz3 = numbytes(logfile)
    assert sz3 > sz2

    logutil.unquiet()

    # get messages from idaes. logger again
    for i in range(2):
        log.info("this is a log message 4")
        time.sleep(0.1)
    dump_console(logfile)
    sz4 = numbytes(logfile)
    assert sz4 > sz3
