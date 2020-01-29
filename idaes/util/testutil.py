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
Utility code for testing IDAES code.
"""
import os
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor


__author__ = "Dan Gunter"


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
                    print(f"ERROR in {output['ename']} in {fullpath} [{num}]:")
                    for tb_line in output["traceback"]:
                        print(tb_line)
    return not failed
