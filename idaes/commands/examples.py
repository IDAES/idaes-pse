"""
Command to get IDAES examples
"""
# stdlib
import logging
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

# third-party
import click
import requests

# package
from idaes.commands.base import command_base
from idaes.ver import package_version as V

_log = logging.getLogger("idaes.commands")


REPO_ORG = "idaes"
REPO_NAME = "examples-pse"
REPO_DIR = "src"
VERSION = f"{V.major}.{V.minor}.{V.micro}"


class DownloadError(Exception):
    pass


@command_base.command(
    name="get-examples", help="Fetch example scripts and Jupyter Notebooks"
)
@click.option(
    "--dir", "-d", "directory", help="installation target directory", default="examples"
)
@click.option(
    "--version",
    "-V",
    help=f"Version of examples to download (default={VERSION})",
    default=VERSION,
)
def get_examples(directory, version):
    target_dir = Path(directory)
    if target_dir.exists():
        _log.error("target directory {target_dir} already exists")
        sys.exit(-1)
    # download
    try:
        download_contents(version, target_dir)
    except DownloadError as err:
        _log.fatal(f"Download failed: {err}")
        sys.exit(-1)
    print(f"Downloaded example files to {target_dir}")


def download_contents(version, target_dir):
    """Download Github contents.
    """
    url = archive_file_url(version)
    _log.info(f"get examples from: {url}")
    # stream out to a big .zip file
    req = requests.get(url, stream=True)
    if req.status_code != 200:
        if req.status_code in (400, 404):
            raise DownloadError(f"file not found")
        raise DownloadError(f"status={req.status_code}")
    tmpdir = TemporaryDirectory()
    tmpfile = Path(tmpdir.name) / "examples.zip"
    with tmpfile.open("wb") as f:
        for chunk in req.iter_content(chunk_size=65536):
            f.write(chunk)
    _log.info(f"downloaded zipfile to {tmpfile}")
    # open as a zip file, and extract all files into the temporary directory
    _log.debug(f"open zip file: {tmpfile}")
    zipf = ZipFile(str(tmpfile))
    zipf.extractall(path=tmpdir.name)
    # move the REPO_DIR subdirectory into the target dir
    subdir = Path(tmpdir.name) / f"{REPO_NAME}-{version}" / REPO_DIR
    _log.debug(f"move {subdir} -> {target_dir}")
    os.rename(str(subdir), str(target_dir))
    del tmpdir


def archive_file_url(version, org=REPO_ORG, repo=REPO_NAME):
    return f"https://github.com/{org}/{repo}/archive/{version}.zip"
