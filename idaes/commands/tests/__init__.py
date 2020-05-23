import os
from pathlib import Path
import pytest
import shutil
import uuid

@pytest.fixture
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
    shutil.rmtree(tempdir)
