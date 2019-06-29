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
import os
from types import ModuleType
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from idaes.examples.workshops import (
    Module_1_Flash_Unit,
    Module_2_Flowsheet,
    Module_3_Custom_Unit_Model,
)


def module_path(modobj: ModuleType):
    return os.path.dirname(modobj.__file__)


def run_notebook(path: str, name: str):
    """Run a specific jupyter notebook 'name' located at `path`.
    """
    fullpath = os.path.join(path, name)
    with open(fullpath) as f:
        nb = nbformat.read(f, as_version=nbformat.NO_CONVERT)

    ep = ExecutePreprocessor(timeout=100)
    ep.allow_errors = True
    ep.preprocess(nb, {"metadata": {"path": path}})

    failed = False
    for cell in nb.cells:
        if "outputs" in cell:
            for output in cell["outputs"]:
                if output.output_type == "error":
                    failed = True
                    num = cell["execution_count"]
                    exc = f"{output['ename']} = {output['evalue']}"
                    print(f"ERROR in {fullpath} [{num}]: {exc}")
    return not failed


# Test the solution jupyter notebooks for all modules


def test_module_1():
    assert run_notebook(
        module_path(Module_1_Flash_Unit), "Module_1_Flash_Unit_Solution.ipynb"
    )


def test_module_2():
    assert run_notebook(
        module_path(Module_2_Flowsheet), "Module_2_Flowsheet_Solution.ipynb"
    )


def test_module_3():
    assert run_notebook(
        module_path(Module_3_Custom_Unit_Model), "Module_3_Exercise_1_Solution.ipynb"
    )
