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
Smoke tests, to make sure things are working at all.
"""
import pytest


@pytest.mark.unit
def test_doalamo_import():
    from idaes.apps.alamopy_depr import alamo


@pytest.mark.integration()
def test_hasalamo():
    from idaes.apps import alamopy_depr as alamopy

    has_alamo_flag = alamopy.multos.has_alamo()
    if has_alamo_flag:
        alamopy.debug["has_alamo"] = True
        version = alamopy.get_alamo_version()
    else:
        alamopy.debug["has_alamo"] = False
        version = "n/a"  # cannot get version w/o alamo present

    print("ALAMO Version: %s" % version)
