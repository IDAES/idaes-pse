#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
# TODO: Missing doc strings
# pylint: disable=missing-module-docstring
# pylint: disable=missing-function-docstring

import importlib.abc
import importlib.machinery
import os
import subprocess
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Union

import pytest

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


def pytest_addoption(parser):
    parser.addoption(
        "--performance",
        action="store_true",
        dest="performance",
        default=False,
        help="enable performance decorated tests",
    )


MARKERS = {
    "build": "test of model build methods",
    "cubic_root": "test requires the compiled cubic root finder",
    "iapws": "test requires the compiled IAPWS95 property package",
    "initialization": "test of initialization methods. These generally require a solver as well",
    "solver": "test requires a solver",
    "ui": "tests of an aspect of the ui",
    "unit": "quick tests that do not require a solver, must run in <2s",
    "component": "quick tests that may require a solver",
    "integration": "long duration tests",
    "performance": "tests for the IDAES performance testing suite",
}


def pytest_configure(config: pytest.Config):
    for name, description in MARKERS.items():
        config.addinivalue_line("markers", f"{name}: {description}")

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
        pytest.skip(f"cannot run on platform {plat}")


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


ModuleName = str


class ImportorskipLoader(importlib.abc.Loader):
    """
    A wrapper class around a concrete Loader instance. If a ModuleNotFoundError is raised
    during module execution and the module name matches one of the registered modules,
    it is replaced with a module-level call to :func:`pyest.skip()`.
    """

    def __init__(
        self, wrapped: importlib.abc.Loader, skip_if_not_found: Iterable[ModuleName]
    ):
        self._wrapped = wrapped
        self.skip_if_not_found = list(skip_if_not_found)

    def module_repr(self, module) -> str:
        return self._wrapped.module_repr(module)

    def create_module(self, spec):
        return self._wrapped.create_module(spec)

    def exec_module(self, module):
        try:
            return self._wrapped.exec_module(module)
        except ModuleNotFoundError as e:
            if e.name in self.skip_if_not_found:
                pytest.skip(allow_module_level=True)
            raise e


class ImportorskipFinder(importlib.abc.MetaPathFinder):
    """
    Custom Finder class to modify import behavior for registered modules.

    If inserted in sys.meta_path before the default finders, it will cause
    a custom Loader to be used for registered modules.
    """

    def __init__(self, registry: dict[ModuleName, list[ModuleName]]):
        self._registry = registry

    def find_spec(self, *args, **kwargs):
        spec = importlib.machinery.PathFinder.find_spec(*args, **kwargs)
        if spec is None:
            return
        registered_for_skipping = self._registry.get(spec.name, None)
        if registered_for_skipping:
            spec.loader = ImportorskipLoader(spec.loader, registered_for_skipping)
        return spec


class Importorskipper:
    """
    A pytest plugin that allows automatically skipping test modules if they attempt to import
    modules distributed in optional dependencies.

    The functionality is similar to pytest.importorskip(), but using a centralized registry
    instead of having to call importorskip() in each test module for each possibly missing module.
    """

    def __init__(self, registry: dict[ModuleName, list[ModuleName]]):
        self._registry = dict(registry)
        self._finder = ImportorskipFinder(self._registry)

    def pytest_configure(self):
        sys.meta_path.insert(0, self._finder)

    def pytest_sessionfinish(self):
        sys.meta_path.remove(self._finder)

    def pytest_report_collectionfinish(self) -> list[str]:
        preamble = [
            "The following modules are registered in the importorskipper plugin",
            " and will cause tests to be skipped if any of the registered modules is not found: ",
        ]
        lines = []
        for importing_mod, mods in self._registry.items():
            lines.append(f"- {importing_mod}:\t{mods}")
        if lines:
            lines = preamble + lines
        return lines


@pytest.fixture(scope="function")
def run_in_tmp_path(tmp_path: Path):
    prev_workdir = os.getcwd()
    try:
        os.chdir(tmp_path)
        yield
    finally:
        os.chdir(prev_workdir)


@pytest.fixture(scope="class")
def run_class_in_tmp_path(
    request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory
):
    prev_workdir = os.getcwd()

    modname = request.module.__name__
    clsname = request.__class__.__name__
    dirname = f"{modname}-{clsname}"

    tmp_path = tmp_path_factory.mktemp(dirname, numbered=True)
    try:
        os.chdir(tmp_path)
        yield
    finally:
        os.chdir(prev_workdir)


@pytest.fixture(scope="module")
def run_module_in_tmp_path(
    request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory
):
    prev_workdir = os.getcwd()

    modname = request.module.__name__

    tmp_path = tmp_path_factory.mktemp(modname, numbered=True)
    try:
        os.chdir(tmp_path)
        yield
    finally:
        os.chdir(prev_workdir)


def _get_repo_root_dir() -> Union[Path, None]:
    try:
        text = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
        ).strip()
    except subprocess.CalledProcessError:
        # no Git directory found
        return None
    return Path(text).resolve()


class VerifyCleanup:
    """
    A pytest plugin to verify that files generated by tests are cleaned up properly.
    """

    def __init__(self, repo_root_dir: Path):
        self._repo_root_dir = Path(repo_root_dir).resolve()
        self._added_by_test = {}

    def _get_files(self) -> list[str]:
        args = [
            "git",
            "-C",
            self._repo_root_dir,
            "ls-files",
            "--others",
            "--exclude-standard",
        ]
        try:
            text = subprocess.check_output(
                args,
                text=True,
            ).strip()
        except subprocess.CalledProcessError as e:
            text = str(e)
        return text.splitlines()

    def pytest_report_collectionfinish(self) -> list[str]:
        return [
            f"Looking for untracked files being created by tests for repo: {self._repo_root_dir}"
        ]

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_protocol(self, item):
        before = set(self._get_files())
        yield
        after = set(self._get_files())
        added = after - before
        if added:
            key = item.nodeid
            self._added_by_test[key] = sorted(added)

    @pytest.hookimpl(trylast=True)
    def pytest_terminal_summary(self, terminalreporter):
        to_report = dict(self._added_by_test)
        if not to_report:
            return

        tr = terminalreporter
        tr.section("Files added (and not cleaned up) by tests")
        for item, files in to_report.items():
            tr.write_line(item)
            for file in files:
                tr.write_line(f"\t{file}")
        if to_report:
            tr.write_line(f"{len(to_report)} test(s) did not clean up after themselves")
            tr.write_line("The exit status of the test session will be set to failed")

    @pytest.hookimpl(trylast=True)
    def pytest_sessionfinish(self, session, exitstatus):
        if self._added_by_test:
            session.exitstatus = pytest.ExitCode.TESTS_FAILED


def pytest_addhooks(pluginmanager: pytest.PytestPluginManager):
    skipper_plugin = Importorskipper(
        {
            "idaes.core.surrogate.keras_surrogate": ["omlt"],
        }
    )

    pluginmanager.register(skipper_plugin, name="importorskipper")
    repo_root_dir = _get_repo_root_dir()
    if repo_root_dir is not None:
        pluginmanager.register(VerifyCleanup(repo_root_dir))


####
