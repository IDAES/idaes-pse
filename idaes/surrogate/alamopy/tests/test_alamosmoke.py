###############################################################################
# ** Copyright Notice **
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
###############################################################################
"""
Smoke tests, to make sure things are working at all.
"""
import pytest

@pytest.mark.unit
def test_doalamo_import():
    from idaes.surrogate.alamopy import alamo


@pytest.mark.integration()
def test_hasalamo():
    from idaes.surrogate import alamopy
    has_alamo_flag = alamopy.multos.has_alamo()
    if has_alamo_flag:
        alamopy.debug['has_alamo'] = True
        version = alamopy.get_alamo_version()
    else:
        alamopy.debug['has_alamo'] = False
        version = "n/a"  # cannot get version w/o alamo present

    print("ALAMO Version: %s"% version)

