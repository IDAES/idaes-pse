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
Tests for Python code style.
"""
import logging
import os
import subprocess

_log = logging.getLogger(__name__)

# we want to look at everything:
# DIRS = ['idaes']
# but for now we only do this:
DIRS = ['idaes/dmf']
FLAKE8 = 'flake8'


def test_flake8():
    cwd = os.getcwd()
    for d in DIRS:
        path = os.path.join(cwd, d)
        if not os.path.exists(path):
            _log.warning('Target path "{}" not found in '
                         'current dir, "{}". Skipping test.'.format(d, cwd))
            continue
        if not os.path.isdir(path):
            _log.warning('Target path "{}" in '
                         'current dir, "{}", '
                         'is not a directory. Skipping test.'.format(d, cwd))
            continue
        cmd = [FLAKE8, d]
        _log.info('Test code style with command "{}"'.format(' '.join(cmd)))
        proc = subprocess.Popen(cmd)
        proc.wait()
        status = proc.returncode
        assert status == 0, 'flake8 had errors for {}'.format(path)
