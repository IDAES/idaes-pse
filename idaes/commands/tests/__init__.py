import os
from pathlib import Path
import pytest
import shutil
import uuid


@pytest.fixture(scope="function")
def random_tempdir():
    """Make a completely cross-platform random temporary directory in
    the current working directory, and yield it. As cleanup, recursively
    remove all contents of this temporary directory.
    """
    origdir = os.getcwd()
    random_name = str(uuid.uuid4())
    os.mkdir(random_name)
    tempdir = Path(random_name)
    yield tempdir
    os.chdir(origdir)
    #
    tempdir_name = str(tempdir)
    del tempdir
    try:
        shutil.rmtree(tempdir_name)
    except Exception as err:
        # this can happen on Windows due to permissions errors
        print("*" * 40)
        print(f"Failed to remove temporary directory: {tempdir_name}: {err}")
        print("*" * 40)

