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
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

path = os.path.dirname(__file__)


# Method to run a specific jupyter notebook
def run_notebook(notebook_name):
    with open(notebook_name) as f:
        nb = nbformat.read(f, as_version=nbformat.NO_CONVERT)

    ep = ExecutePreprocessor(timeout=100)
    ep.allow_errors = True
    ep.preprocess(nb, {'metadata':
                       {'path': os.path.dirname(notebook_name)}})

    errors = []
    for cell in nb.cells:
        if 'outputs' in cell:
            for output in cell['outputs']:
                if output.output_type == 'error':
                    errors.append(output)
    return errors


# Test the solution jupyter notebooks for all modules
# Module 1
def test_module_1():
    errors = \
        run_notebook(os.path.join(path, '../Module_1_Flash_Unit/'
                                  'Module_1_Flash_Unit_Solution.ipynb'))
    assert errors == []


# Module 2
def test_module_2():
    errors = \
        run_notebook(os.path.join(path, '../Module_2_Flowsheet/'
                                  'Module_2_Flowsheet_Solution.ipynb'))
    assert errors == []


# Module 3
def test_module_3():
    errors = \
        run_notebook(os.path.join(path, '../Module_3_Custom_Unit_Model/'
                                  'Module_3_Exercise_1_Solution.ipynb'))
    assert errors == []
