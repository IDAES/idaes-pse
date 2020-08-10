import logging
import os
from pathlib import Path
import pytest
import shutil
import uuid


_log = logging.getLogger(__name__)
if os.environ.get("IDAES_TEST_DEBUG", False):
    _log.setLevel(logging.DEBUG)


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
        _log.debug(f"Removing temporary directory '{tempdir_name}'")
        shutil.rmtree(tempdir_name)
    except Exception as err:
        # this can happen on Windows due to permissions errors
        _log.warning(f"Failed to remove temporary directory '{tempdir_name}': {err}")

