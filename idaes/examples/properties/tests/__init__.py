# stdlib
import os
import pathlib
import shutil
from typing import Tuple

# Utility functions used in tests


def save_config_yaml(nbpath: str) -> Tuple[str, str]:
    """Save `config.yaml` file at given path to a temporary copy.
    Copy will have same metadata as source file.

    Returns:
         (source filename, filename of the temporary file)
    """
    filename = str(pathlib.Path(nbpath) / "config.yaml")
    saved_filename = filename + ".tmp"
    shutil.copy2(filename, saved_filename)
    return filename, saved_filename


def restore_file(src: str, dst: str):
    """Copy `src` to `dst` and remove `src`.
    Preserve all metadata on `src`.
    """
    shutil.copy2(src, dst)
    os.unlink(src)

