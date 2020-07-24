import pytest
import sys

####
# Uncomment this to collect list of all test files
# into 'test_files.txt'
#
# def pytest_collection_modifyitems(config, items):
#     output = open("test_files.txt", "w")
#     fspaths = set()
#     for item in items:
#         fspaths.add(item.fspath)
#     for p in fspaths:
#         output.write(str(p))
#         output.write("\n")
#     output.close()
####

####
# Only run the tests for your particular platform(s), if it is marked with those platform(s)
# e.g. to mark tests for linux:
#
# import pytest
# @pytest.mark.linux
# def test_something():
#    print("this only runs on linux")
#
#
# In addition, you can use "no<platform>" to exclude
#
# @pytest.mark.nowin32
# def test_something():
#    print("this will not run on windows")
#
# The names of the platforms should match what is returned by `sys.platform`, in particular:
#    Linux = 'linux'
#    Windows = 'win32'
#    macOS = 'darwin'


ALL = {"darwin", "linux", "win32"}
ALL_NO = {"no" + tag for tag in ALL}


def pytest_runtest_setup(item):
    supported_platforms = ALL.intersection(mark.name for mark in item.iter_markers())
    excluded_platforms = ALL_NO.intersection(mark.name for mark in item.iter_markers())
    plat = sys.platform
    if ((excluded_platforms and ("no" + plat) in excluded_platforms) or
        (supported_platforms and plat not in supported_platforms)):
        pytest.skip("cannot run on platform {}".format(plat))

####
