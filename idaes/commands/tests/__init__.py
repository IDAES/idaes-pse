import logging
import os
from pathlib import Path
from shutil import rmtree


def create_module_scratch(module_name):
    """Create a scratch that is easily found (by developers) and removed, and doesn't
    clutter up the working directories.
    """
    path = Path("~").expanduser()
    for path_part in ".idaes", "_scratch", module_name:
        path = path / path_part
        if not path.exists():
            path.mkdir()
    return path


def rmtree_scratch(scratch_path):
    """Do our level best to remove all the temporary files.
    """
    try:
        rmtree(scratch_path, ignore_errors=True)
    except Exception as err:
        pass
