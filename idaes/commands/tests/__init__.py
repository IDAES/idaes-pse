import os
from pathlib import Path
import pytest
import shutil
import uuid


@pytest.fixture(scope="session")
def tempdir_list():
    td_list = []
    yield td_list
    for td in td_list:
        try:
            shutil.rmtree(td)
        except Exception as err:
            # this can happen on Windows due to permissions errors
            print("*" * 40)
            print(f"Failed to remove temporary directory: {td}: {err}")
            print("*" * 40)


@pytest.fixture(scope="function")
def random_tempdir(tempdir_list):
    """Make a completely cross-platform random temporary directory in
    the current working directory, and yield it. As cleanup, recursively
    remove all contents of this temporary directory.
    """
    origdir = os.getcwd()
    random_name = str(uuid.uuid4())
    os.mkdir(random_name)
    tempdir = Path(random_name)
    tempdir_list.append(str(tempdir))
    yield tempdir
    os.chdir(origdir)
