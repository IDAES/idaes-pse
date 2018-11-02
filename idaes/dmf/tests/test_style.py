##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
# 
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for Python code style.
"""
import os
import subprocess

DIRS = ['idaes/dmf']
FLAKE8 = 'flake8'


def test_flake8():
    cwd = os.getcwd()
    for d in DIRS:
        path = os.path.join(cwd, d)
        if not os.path.exists(path):
            raise os.error('Target path "{}" not found in '
                           'current dir, "{}"'.format(path, d))
        if not os.path.isdir(path):
            raise os.error('Target path "{}" in '
                           'current dir, "{}", '
                           'is not a directory'.format(path, d))
        cmd = [FLAKE8, d]
        #print('@@cmd={}'.format(cmd))
        proc = subprocess.Popen(cmd)
        proc.wait()
        status = proc.returncode
        assert status == 0, 'flake8 had errors for {}'.format(path)
