import pytest
import sys

output = open("test_files.txt", "w")


def pytest_collection_modifyitems(config, items):
    fspaths = set()
    for item in items:
        fspaths.add(item.fspath)
    for p in fspaths:
        output.write(str(p))
        output.write("\n")
    output.close()

####
# Only run the tests for your particular platform(s), if it is marked with those platform(s)
# e.g. to mark tests for linux:
#
# import pytest
# @pytest.mark.linux
# def test_something():
#    print("this only runs on linux")
#
# The names of the platforms should match what is returned by `sys.platform`, in particular:
#    Linux = 'linux'
#    Windows = 'win32'
#    macOS = 'darwin'


ALL = set("darwin linux win32".split())


def pytest_runtest_setup(item):
    supported_platforms = ALL.intersection(mark.name for mark in item.iter_markers())
    plat = sys.platform
    if supported_platforms and plat not in supported_platforms:
        pytest.skip("cannot run on platform {}".format(plat))

####
