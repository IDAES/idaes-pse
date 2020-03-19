# coding: utf-8
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

import pytest
import idaes.beta

def test_beta_module_exception():
    with pytest.raises(
            ImportError, match=r"Module 'idaes.tests.beta_mod' is in beta "
            "and must be imported using idaes.beta.import_beta\(\)."):
        import idaes.tests.beta_mod

def test_beta_module_import():
    mod = idaes.beta.import_beta('idaes.tests.beta_mod')
    assert mod.__name__ == 'idaes.tests.beta_mod'
    assert mod.reference_value == 42

    mod = idaes.beta.import_beta('.beta_mod')
    assert mod.__name__ == 'idaes.tests.beta_mod'
    assert mod.reference_value == 42
