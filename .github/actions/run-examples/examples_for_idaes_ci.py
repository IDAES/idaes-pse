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


def pytest_collection_modifyitems(items):
    for item in items:
        parts = item.path.parts
        if "held" in parts:
            item.add_marker(
                pytest.mark.xfail(run=False, reason="notebook has 'held' status")
            )
        if "archive" in parts:
            item.add_marker(pytest.mark.skip(reason="notebook is archived"))


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
        args += ["-c", os.devnull]

    with _temp_cwd(rootdir):
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
