#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
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


REQUIRED_MARKERS = {"unit", "component", "integration", "performance"}
ALL_PLATFORMS = {"darwin", "linux", "win32"}


@pytest.hookimpl
def pytest_runtest_setup(item):
    _validate_required_markers(
        item, required_markers=REQUIRED_MARKERS, expected_count=1
    )
    _skip_for_unsupported_platforms(
        item,
        all_platforms=ALL_PLATFORMS,
        negate_tag=(lambda tag: f"no{tag}"),
    )


def _skip_for_unsupported_platforms(item, all_platforms=None, negate_tag=None):
    """
    Only run the tests for your particular platform(s), if it is marked with those platform(s)
    e.g. to mark tests for linux:

    import pytest
    @pytest.mark.linux
    def test_something():
       print("this only runs on linux")

    In addition, you can use "no<platform>" to exclude

    @pytest.mark.nowin32
    def test_something():
       print("this will not run on windows")

    The names of the platforms should match what is returned by `sys.platform`, in particular:
       Linux = 'linux'
       Windows = 'win32'
       macOS = 'darwin'
    """

    all_platforms = set(all_platforms or [])
    if negate_tag is None:

        def negate_tag(tag):
            return f"no{tag}"

    all_negated_platforms = {negate_tag(tag) for tag in all_platforms}
    item_markers = {marker.name for marker in item.iter_markers()}
    supported_platforms = all_platforms & item_markers
    excluded_platforms = all_negated_platforms & item_markers
    plat = sys.platform
    if (excluded_platforms and negate_tag(plat) in excluded_platforms) or (
        supported_platforms and plat not in supported_platforms
    ):
        pytest.skip("cannot run on platform {}".format(plat))


def _validate_required_markers(item, required_markers=None, expected_count=1):
    required_markers = set(required_markers or [])
    item_markers = {marker.name for marker in item.iter_markers()}
    required_markers_on_item = item_markers & required_markers
    required_count = len(required_markers_on_item)
    reason_to_fail = None
    if required_count < expected_count:
        reason_to_fail = "Too few required markers"
    if required_count > expected_count:
        reason_to_fail = "Too many required markers"
    if reason_to_fail:
        msg = (
            f'{reason_to_fail} for test function "{item.name}". '
            f"Expected: {expected_count} of {required_markers}; "
            f"found: {required_markers_on_item or required_count}"
        )
        pytest.fail(msg)


def pytest_addoption(parser):
    parser.addoption(
        "--performance",
        action="store_true",
        dest="performance",
        default=False,
        help="enable performance decorated tests",
    )


def pytest_configure(config):
    if not config.option.performance:
        if len(config.option.markexpr) > 0:
            setattr(
                config.option,
                "markexpr",
                f"{config.option.markexpr} and not performance",
            )
        else:
            setattr(config.option, "markexpr", "not performance")
    else:
        setattr(config.option, "markexpr", "performance")


####
