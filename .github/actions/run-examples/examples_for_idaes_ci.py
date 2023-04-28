"""
pytest plugin for testing IDAES "through" IDAES/examples within the IDAES/idaes-pse CI
"""
from contextlib import contextmanager
from dataclasses import dataclass, field
import logging
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


def pytest_configure(config: pytest.Config):
    config.stash[matchmarker] = {
        "*/held/*": pytest.mark.xfail(run=False, reason="notebook has 'held' status"),
        "*/archive/*": pytest.mark.skip(reason="notebook is archived"),
        # TODO: Need to fix this once the Python 3.11 issue is resolved in tensorflow
        "*/surrogates/best_practices_optimization*": pytest.mark.xfail(
            condition=sys.version_info > (3, 11),
            run=True,
            strict=False,
            reason="tensorflow ImportError on 3.11",
        ),
        "*/surrogates/omlt/keras_flowsheet_optimization*": pytest.mark.xfail(
            condition=sys.version_info > (3, 11),
            run=True,
            strict=False,
            reason="tensorflow ImportError on 3.11",
        ),
    }
    config.stash[marked] = []


def pytest_sessionstart(session: pytest.Session):
    build.preprocess(session.config.rootpath)


def pytest_ignore_collect(collection_path: Path, config: pytest.Config):
    if "_dev" in collection_path.parts:
        return True
    if not collection_path.is_file():
        return
    if collection_path.suffix == ".py":
        # specifically ignore python files
        return True
    if not collection_path.match("**/*_test.ipynb"):
        return True


def pytest_collection_modifyitems(config: pytest.Config, items):
    for item in items:
        path = item.path
        for pattern, marker in config.stash[matchmarker].items():
            if path.match(pattern):
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

    with _temp_cwd(rootdir) as p:
        res = pytest.main(
            args,
            **kwargs,
        )

    return res


def main(args):
    rootdir = Path(idaes_examples.__path__[0])

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
