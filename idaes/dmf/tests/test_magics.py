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
Tests for DMF Jupyter "magics"
"""
# third-party
import pytest
# local
from idaes.dmf import magics
from .util import TempDir


class MockShell(object):
    """Mock object for IPython 'shell'.
    """
    def ev(self, name):
        raise ValueError('Unknown object: {}'.format(name))


@pytest.fixture
def magics_impl():
    shell = MockShell()
    return magics.DmfMagicsImpl(shell)


def test_init(magics_impl):
    assert 1


def test_dmf_init_create(magics_impl):
    with TempDir() as d:
        magics_impl.dmf_init(d, 'create')
        assert magics_impl.initialized
