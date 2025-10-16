"""
pytest plugin for testing IDAES "through" IDAES/examples within the IDAES/idaes-pse CI
"""

from contextlib import contextmanager
import fnmatch
import os
from pathlib import Path
import sys
from typing import List

import pytest

import idaes_examples
from idaes_examples import notebooks
from idaes_examples import build


matchmarker = pytest.StashKey()
marked = pytest.StashKey()

# Environment triggers for use either locally
# or in GitHub Action
_ENV_TRUE = {"1", "true", "yes", "on"}
RUN_NOTEBOOKS = str(os.getenv("EXAMPLES_RUN_NOTEBOOKS", "1")).lower() in _ENV_TRUE
RUN_PYFILES = str(os.getenv("EXAMPLES_RUN_PYFILES", "1")).lower() in _ENV_TRUE


def _matches_pattern(item: pytest.Item, pattern: str) -> bool:
    to_match = os.fspath(item.path)
    return fnmatch.fnmatch(to_match, pattern)


def pytest_configure(config: pytest.Config):
    config.stash[matchmarker] = {
        "*/held/*": pytest.mark.xfail(run=False, reason="notebook has 'held' status"),
        "*/archive/*": pytest.mark.skip(reason="notebook is archived"),
        "*/surrogates/sco2/alamo/*": pytest.mark.xfail(
            run=False,
            reason="notebooks require ALAMO to run",
        ),
    }
    config.stash[marked] = []


def pytest_sessionstart(session: pytest.Session):
    build.preprocess(session.config.rootpath)


def pytest_ignore_collect(collection_path: Path, config: pytest.Config):
    """Control what gets collected. By default, notebooks and tests; ignore other files."""
    if "_dev" in collection_path.parts:
        return True
    if not collection_path.is_file():
        return

    suffix = collection_path.suffix
    # Notebook tests (only *_test.ipynb)
    if suffix == ".ipynb":
        if not RUN_NOTEBOOKS:
            return True
        return not collection_path.match("**/*_test.ipynb")

    # Python tests (only test_*.py or *_test.py)
    if suffix == ".py":
        if not RUN_PYFILES:
            return True
        name = collection_path.name
        if fnmatch.fnmatch(name, "test_*.py") or fnmatch.fnmatch(name, "*_test.py"):
            return False
        return True

    # Ignore everything else
    return True


def pytest_collection_modifyitems(config: pytest.Config, items):
    for item in items:
        for pattern, marker in config.stash[matchmarker].items():
            if _matches_pattern(item, pattern):
                item.add_marker(marker)
                config.stash[marked].append((item, pattern, marker))


def pytest_report_collectionfinish(config: pytest.Session):
    lines = []
    for item, pattern, marker in config.stash[marked]:
        lines += [f"\t{item.nodeid=}\t{pattern=}\t{marker.mark}"]
    if lines:
        lines.insert(0, "The following items have been marked:")
    return lines


@contextmanager
def _temp_cwd(path: Path):
    orig_wdir = Path.cwd()
    try:
        os.chdir(path)
        yield path
    finally:
        os.chdir(orig_wdir)


def run_pytest(
    rootdir: Path,
    args: List[str],
    ignore_conftest: bool = True,
    ignore_inifile: bool = True,
    **kwargs,
):
    args = [
        str(rootdir),
        f"--rootdir={rootdir}",
        *args,
    ]

    if ignore_conftest:
        args += ["--noconftest"]
    if ignore_inifile:
        empty_file_for_ignoring = rootdir / "__empty__"
        empty_file_for_ignoring.write_text("")
        args += ["-c", empty_file_for_ignoring]

    with _temp_cwd(rootdir):
        res = pytest.main(
            args,
            **kwargs,
        )

    return res


def main(args):
    rootdir = Path(idaes_examples.__path__[0])

    # Always include --nbmake, but notebook collection can be disabled by env
    res = run_pytest(
        rootdir,
        [
            "--nbmake",
            *args,
        ],
        ignore_conftest=True,
        ignore_inifile=True,
        plugins=[__name__],
    )

    return res


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
