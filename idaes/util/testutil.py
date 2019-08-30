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
                    exc = f"{output['ename']} = {output['evalue']}"
                    print(f"ERROR in {fullpath} [{num}]: {exc}")
    return not failed
