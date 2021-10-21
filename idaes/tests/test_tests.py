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
"Tests used to test the functionality of the test suite itself."

import pytest


# pytestmark = pytest.mark.skip("Tests for the test suite are always skipped under normal circumstances")


@pytest.mark.xfail
def test_function_without_any_required_mark():
    assert True


@pytest.mark.xfail
@pytest.mark.unit
@pytest.mark.component
def test_function_with_too_many_required_marks():
    assert True
