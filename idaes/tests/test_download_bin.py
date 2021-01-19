import idaes.solvers
import pytest
import os

def _del_data_file(path):
    path = os.path.join(idaes.data_directory, path)
    try:
        os.remove(path)
    except OSError:
        pass

def _exists(path):
    path = os.path.join(idaes.data_directory, path)
    return os.path.exists(path)

def test_dl_bin():
    _del_data_file("download_tests/version_lib.txt")
    _del_data_file("download_tests/version_solvers.txt")
    idaes.solvers.download_binaries(
        release=idaes.config.default_binary_release,
        to_path="download_tests")
    assert(_exists("download_tests/version_lib.txt"))
    assert(_exists("download_tests/version_solvers.txt"))
