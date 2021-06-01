# coding: utf-8
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

import logging
import pytest
import sys

from io import StringIO
from pyomo.common.log import LoggingIntercept
from idaes.beta import import_beta

@pytest.mark.unit
def test_beta_module_exception():
    with pytest.raises(
            ImportError, match=r"Module 'idaes.tests.beta_mod' is in beta "
            "and must be imported using idaes.beta.import_beta\(\)."):
        import idaes.tests.beta_mod

@pytest.mark.unit
def test_beta_module_import():
    os = StringIO()
    with LoggingIntercept(os, 'idaes', logging.INFO):
        mod = import_beta('idaes.tests.beta_mod')
    assert mod.__name__ == 'idaes.tests.beta_mod'
    assert mod.reference_value == 42
    assert os.getvalue() == ""

    os = StringIO()
    with LoggingIntercept(os, 'idaes', logging.INFO):
        mod = import_beta('.beta_mod')
    assert mod.__name__ == 'idaes.tests.beta_mod'
    assert mod.reference_value == 42
    assert os.getvalue() == ""

    os = StringIO()
    with LoggingIntercept(os, 'idaes.logger', logging.INFO):
        mod = import_beta('idaes.logger')
    assert mod.__name__ == 'idaes.logger'
    assert os.getvalue().strip() == "Module 'idaes.tests.test_beta' "\
        "imported module 'idaes.logger' as a Beta module.  "\
        "This module is not declared beta and can be "\
        "imported using Python's normal import mechanisms."

@pytest.mark.unit
def test_beta_reimport():
    os = StringIO()
    with LoggingIntercept(os, 'idaes', logging.INFO):
        mod = import_beta('idaes.tests.beta_mod')
    assert mod.__name__ == 'idaes.tests.beta_mod'
    assert mod.reference_value == 42
    assert os.getvalue() == ""

    # Test that subsequent standard imports will still trigger an
    # exception, even if the module was previously imported using
    # import-beta()
    with pytest.raises(
            ImportError, match=r"Module 'idaes.tests.beta_mod' is in beta "
            "and must be imported using idaes.beta.import_beta\(\)."):
        import idaes.tests.beta_mod

    # But subsequent beta_import() modules return the same (original)
    # module object
    os = StringIO()
    with LoggingIntercept(os, 'idaes', logging.INFO):
        mod2 = import_beta('idaes.tests.beta_mod')
    assert mod2 is mod
