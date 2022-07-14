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
Tests for versioning
"""
# third-party
import pytest

# pkg
import idaes
from idaes import ver


@pytest.mark.unit
def test_idaes_version():
    assert idaes.__version__


@pytest.mark.unit
def test_ver_class():
    v = ver.Version(1, 2, 3)
    assert str(v) == "1.2.3"
    seq = [1, 2, 3]
    for i, n in enumerate(v):
        print("i:", i, "n:", n)
        assert n == seq[i]
    v = ver.Version(1, 2, 3, "beta", 1)
    assert str(v) == "1.2.3.b1"
    v = ver.Version(1, 2, 3, "beta", 2, "test")
    assert str(v) == "1.2.3.b2+test"
    v = ver.Version(1, 2, 3, "beta", label="test")
    assert str(v) == "1.2.3.b+test"
    v = ver.Version(1, 2, 3, "development")
    assert str(v) == "1.2.3.dev"
    pytest.raises(ValueError, ver.Version, 1, 2, 3, "howdy")


class MyVersionedClass(ver.HasVersion):
    def __init__(self):
        super(MyVersionedClass, self).__init__(1, 2, 3)


@pytest.mark.unit
def test_has_version():
    x = MyVersionedClass()
    assert str(x.version) == "1.2.3"


@pytest.mark.unit
def test_bump_version():
    v = ver.Version(1, 2, 3)
    assert tuple(v) == (1, 2, 3)
    v.micro += 1
    assert tuple(v) == (1, 2, 4)
