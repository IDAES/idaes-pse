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
Tests for Workshop Modules.

Author: Jaffer Ghouse
"""
# stdlib
import os
from types import ModuleType
# third party
import pytest
# package
from idaes.util.testutil import run_notebook
from idaes.examples.workshops import (
    Module_1_Flash_Unit,
    Module_2_Flowsheet,
    Module_3_Custom_Unit_Model,
)


def module_path(modobj: ModuleType):
    return os.path.dirname(modobj.__file__)


# Test the 'solution' notebooks for all modules
# Note: modules 2 & 3 are skipped in circleci due to failure in free ipopt

def test_module_1():
    assert run_notebook(
        module_path(Module_1_Flash_Unit), "Module_1_Flash_Unit_Solution.ipynb"
    )


@pytest.mark.nocircleci()
def test_module_2():
    assert run_notebook(
        module_path(Module_2_Flowsheet), "Module_2_Flowsheet_Solution.ipynb"
    )


@pytest.mark.nocircleci()
def test_module_3():
    assert run_notebook(
        module_path(Module_3_Custom_Unit_Model), "Module_3_Exercise_1_Solution.ipynb"
    )
