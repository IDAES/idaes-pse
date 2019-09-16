"""
Test Python notebooks in property example directories.

Generally, anything that requires a full optimization should be decorated with
`@pytest.mark.nocircleci()`, so our CircleCI tests don't take forever.
"""
# stdlib
import glob
import os
import pathlib
import shutil
# third-party
import pytest
# package
from idaes.util.testutil import run_notebook

__author__ = "Dan Gunter"


root = pathlib.Path(__file__).parent.parent


def move_to_tmp(module_name, tmp_path):
    d = tmp_path / module_name
    d.mkdir()
    src = root / module_name
    # copy source code into temporary dir
    pattern = str(src / "*.py")
    for ext in "py", "ipynb", "yaml", "png", "json":
        for fname in glob.glob(str(src / f"*.{ext}")):
            print(f"COPY {fname} to {str(d)}")
            shutil.copy(fname, str(d))
    # change to temporary dir as working directory
    os.chdir(str(d))
    # return dir
    return d


@pytest.mark.nocircleci()
def test_module_2(tmp_path):
    module_name = "Workshop_Module_2"
    d = move_to_tmp(module_name, tmp_path)
    assert run_notebook(".", "Module_2_Flowsheet_DMF_Solution.ipynb")
