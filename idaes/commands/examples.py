"""
Command to get IDAES examples
"""
# stdlib
import logging
from operator import itemgetter
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


GITHUB = "https://github.com"
GITHUB_API = "https://api.github.com"
REPO_ORG = "idaes"
REPO_NAME = "examples-pse"
REPO_DIR = "src"
PKG_VERSION = f"{V.major}.{V.minor}.{V.micro}"

# XXX: probably don't *really* want pre-releases, but they are
# XXX: handy for testing this script.
_skip_prereleases = False


class DownloadError(Exception):
    """Used for errors downloading the release files.
    """
    pass


@command_base.command(
    name="get-examples", help="Fetch example scripts and Jupyter Notebooks"
)
@click.option(
    "--dir", "-d", "directory", help="installation target directory", default="examples"
)
@click.option(
    "--list-releases",
    "-l",
    help="List all available released versions, and stop",
    flag_value="list",
)
@click.option(
    "--version",
    "-V",
    help=f"Version of examples to download (default={PKG_VERSION})",
    default=PKG_VERSION,
)
def get_examples(directory, list_releases, version):
    """Get the examples from Github and put them in a local directory.
    """
    # check version
    releases = get_releases()
    if list_releases == "list":
        print_releases(releases)
        sys.exit(0)
    matched_release = None
    for rel in releases:
        if version == rel[1]:  # matches tag
            matched_release = rel
            break
    if matched_release is None:
        # modify message slightly depending on whether they selected a version
        if version == PKG_VERSION:
            how = "installed"
            how_ver = " and -V/--version to choose a desired version"
        else:
            how = "selected"
            how_ver = ""  # they already did this!
        click.echo(
            f"No release found matching {how} IDAES package version '{version}'."
        )
        click.echo(f"Use -l/--list-releases to see all{how_ver}.")
        sys.exit(-1)
    # check target directory
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
    """Download the given version from the Github releases and make
    its `REPO_DIR` subdirectory be the `target_dir`.

    Raises:
        DownloadError: if the GET on the release URL returns non-200 status
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
    _log.debug(f"created temporary directory '{tmpdir.name}'")
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
    _log.debug(f"removing temporary directory '{tmpdir.name}'")
    del tmpdir


def archive_file_url(version, org=REPO_ORG, repo=REPO_NAME):
    """Build & return URL for a given release version.
    """
    return f"{GITHUB}/{org}/{repo}/archive/{version}.zip"


def get_releases():
    """Returns a list of releases, with a tuple for each of (date, tag, info).
    The list is sorted in ascending order by date.
    """
    releases = []
    url = f"{GITHUB_API}/repos/{REPO_ORG}/{REPO_NAME}/releases"
    req = requests.get(url)
    for rel in req.json():
        if _skip_prereleases and rel["prerelease"]:
            continue
        releases.append((rel["published_at"], rel["tag_name"], rel["name"]))
    releases.sort(key=itemgetter(0))  # sort by publication date
    return releases


def print_releases(releases):
    """Print the releases, as returned by `get_releases()`, as a table
    to standard output.
    """
    if len(releases) == 0:
        print("No releases to list")
        return
    # determine column widths
    widths = [4, 7, 7]  # widths of column titles: date,version,details
    widths[0] = len(releases[0][0])  # dates are all the same
    # tags and names can have different widths
    for rel in releases:
        for i in range(1, 3):
            widths[i] = max(widths[i], len(rel[i]))
    # make row format
    pad = "  "
    fmt = f"{{date:{widths[0]}s}}{pad}{{tag:{widths[1]}s}}{pad}{{name:{widths[2]}s}}"
    # print header
    print("")
    print(fmt.format(date="Date", tag="Version", name="Details"))
    print(fmt.format(date="-" * widths[0], tag="-" * widths[1], name="-" * widths[2]))
    # print rows
    for rel in releases:
        print(fmt.format(date=rel[0], tag=rel[1], name=rel[2]))
    # print footer
    print("")
