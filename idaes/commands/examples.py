##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Install IDAES example files locally.

By default, this will download the examples from "examples-pse/src" on Github,
with a release matching the version of the currently installed idaes package,
into a sub-directory of the current directory named "examples". It will
also, unless directed otherwise, install all Python modules from the downloaded
directory into a package called "idaes_examples".

Options let the user choose a different version, directory, and
whether to actually download or install.
"""
# stdlib
from collections import namedtuple
from io import StringIO
import logging
from operator import attrgetter
import os
from pathlib import Path
import re
import shutil
import sys
from typing import List
from zipfile import ZipFile

# third-party
import click
import requests

# package
from idaes.commands.base import command_base
from idaes.ver import package_version as V
from idaes.util.system import TemporaryDirectory

__author__ = "Dan Gunter"


_log = logging.getLogger("idaes.commands.examples")

# Constants

GITHUB = "https://github.com"
GITHUB_API = "https://api.github.com"
REPO_ORG = "idaes"
REPO_NAME = "examples-pse"
REPO_DIR = "src"
PKG_VERSION = f"{V.major}.{V.minor}.{V.micro}"
INSTALL_PKG = "idaes_examples"


class DownloadError(Exception):
    """Used for errors downloading the release files.
    """

    pass


class InstallError(Exception):
    """Used for errors installing the source as a Python package.
    """

    pass


Release = namedtuple("Release", ["date", "tag", "info"])


@command_base.command(
    name="get-examples", help="Fetch example scripts and Jupyter Notebooks."
)
@click.option(
    "--dir", "-d", "directory", help="installation target directory", default="examples",
    type=str,
)
@click.option(
    "--force",
    "-f",
    "force_write",
    help="Do not prompt for remove/overwrite of existing examples",
    is_flag=True
)
@click.option(
    "--no-install", "-I", "no_install", help="Do *not* install examples into 'idaes_examples' package",
    is_flag=True
)
@click.option(
    "--list-releases",
    "-l",
    help="List all available released versions, and stop",
    is_flag=True,
)
@click.option(
    "--no-download",
    "-N",
    "no_download",
    help="Do not download anything",
    is_flag=True
)
@click.option(
    "--unstable",
    "-U",
    help="Allow and list unstable/pre-release versions",
    is_flag=True
)
@click.option(
    "--version",
    "-V",
    help=f"Version of examples to download",
    default=PKG_VERSION,
    show_default=True,
)
def get_examples(directory, force_write, no_install, list_releases, no_download, version,
                 unstable):
    """Get the examples from Github and put them in a local directory.
    """
    # list-releases mode
    if list_releases:
        releases = get_releases(unstable)
        print_releases(releases, unstable)
        sys.exit(0)
    # otherwise..
    target_dir = Path(directory)
    # no-download mode
    if no_download:
        _log.info("skipping download")
    else:
        stable_ver = re.match(r".*-\w+$", version) is None
        if not stable_ver and not unstable:
            click.echo(f"Cannot download unstable version {version} unless you add "
                       f"the -U/--unstable flag")
            sys.exit(-1)
        click.echo("Downloading...")
        try:
            download(get_releases(unstable), target_dir, version, force_write)
        except DownloadError as err:
            _log.fatal(f"abort due to failed download: {err}")
            sys.exit(-1)
        full_dir = os.path.realpath(target_dir)
        click.echo(f"* Downloaded examples to directory '{full_dir}'")
    # install
    if not no_install:
        click.echo("Installing...")
        try:
            install_src(version, target_dir)
        except InstallError as err:
            click.echo(f"Install error: {err}")
            sys.exit(-1)
        click.echo(f"* Installed examples in package {INSTALL_PKG}")


def download(releases, target_dir, version, force):
    """Download `version` into `target_dir`.

    Raises:
        DownloadError
    """
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
            how = f"selected"
            how_ver = ""  # they already did this!
        click.echo(
            f"No release found matching {how} IDAES package version '{version}'."
        )
        click.echo(f"Use -l/--list-releases to see all{how_ver}.")
        raise DownloadError("bad version")
    # check target directory
    if target_dir.exists():
        why_illegal = is_illegal_dir(target_dir)
        if why_illegal:
            raise DownloadError(f"Refuse to remove directory '{target_dir}': "
                                f"{why_illegal}")
        if not force:
            click.confirm(f"Replace existing directory '{target_dir}'", abort=True)
        try:
            shutil.rmtree(target_dir)
        except Exception as err:
            raise DownloadError(f"Cannot remove existing directory: {err}")
    # download
    try:
        download_contents(version, target_dir)
    except DownloadError as err:
        click.echo(f"Download failed: {err}")
        shutil.rmtree(target_dir)  # remove partial download
        raise


def is_illegal_dir(d: Path):
    """Refuse to remove directories for some situations, for safety.
    """
    if (d / ".git").exists():
        return ".git file found"
    if d.absolute() == Path.home().absolute():
        return "cannot replace home directory"
    if d.absolute == Path("/").absolute():
        return "cannot replace root directory"
    return None


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


def get_releases(unstable) -> List[Release]:
    """Returns a list of releases.

    The list is sorted in ascending order by date.
    """
    releases = []
    url = f"{GITHUB_API}/repos/{REPO_ORG}/{REPO_NAME}/releases"
    req = requests.get(url)
    for rel in req.json():
        if not unstable and rel["prerelease"]:
            continue
        releases.append(Release(rel["published_at"], rel["tag_name"], rel["name"]))
    releases.sort(key=attrgetter("date"))  # sort by publication date
    return releases


def print_releases(releases: List[Release], unstable):
    """Print the releases, as returned by `get_releases()`, as a table
    to standard output.
    """
    if len(releases) == 0:
        if unstable:
            print("No releases found")
        else:
            print("No stable releases found. Add -U/--unstable to also look "
                  "for pre-releases.")
        return
    # determine column widths
    widths = [4, 7, 7]  # widths of column titles: date,version,details
    widths[0] = len(releases[0].date)  # dates are all the same
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
        print(fmt.format(date=rel.date, tag=rel.tag, name=rel.info))
    # print footer
    print("")


def install_src(version, target_dir):
    """Install the 'src' subdirectory as a package, given by `INSTALL_PKG`,
    by renaming the directory, adding '__init__.py' files,
    and then running `setuptools.setup()` on the directory tree.
    When done, name the directory back to 'src', and remove '__init__.py' files
    """
    from setuptools import setup, find_packages
    root_dir = target_dir.parent
    examples_dir = root_dir / INSTALL_PKG
    if examples_dir.exists():
        raise InstallError(f"package directory {examples_dir} already exists")
    _log.info(f"install into {INSTALL_PKG} package")
    # set the args to make it look like the 'install' command has been invoked
    saved_args = sys.argv[:]
    sys.argv = ["setup.py", "install"]
    # add some empty __init__.py files
    _log.debug("add temporary __init__.py files")
    pydirs = find_python_directories(target_dir)
    pydirs.append(target_dir)  # include top-level dir
    for d in pydirs:
        init_py = d / "__init__.py"
        init_py.open("w")
    # temporarily rename target directory to the package name
    os.rename(target_dir, examples_dir)
    # if there is a 'build' directory, move it aside
    build_dir = root_dir / 'build'
    if build_dir.exists():
        from uuid import uuid1
        random_letters = str(uuid1())
        moved_build_dir = f"{build_dir}.{random_letters}"
        _log.debug(f"move existing build dir to {moved_build_dir}")
        os.rename(str(build_dir), moved_build_dir)
    else:
        _log.debug("no existing build directory (nothing to do)")
        moved_build_dir = None
    # run setuptools' setup command (in root directory)
    _log.info(f"run setup command in directory {root_dir}")
    orig_dir = os.curdir
    os.chdir(root_dir)
    packages = [d for d in find_packages() if d.startswith(INSTALL_PKG)]
    _log.debug(f"install packages: {packages}")
    # before running, grab stdout
    orig_stdout = sys.stdout
    sys.stdout = setup_out = StringIO()
    # run setup
    setup(
        name=INSTALL_PKG,
        version=version,
        # description='IDAES examples',
        packages=packages,
        python_requires=">=3.5,  <4",
        zip_safe=False
    )
    # restore stdout
    sys.stdout = orig_stdout
    # print/log output
    output_str = setup_out.getvalue()
    if _log.isEnabledFor(logging.DEBUG):
        for line in output_str.split("\n"):
            _log.debug(f"(setup) {line}")
    # name the target directory back to original
    os.rename(examples_dir, target_dir)
    # remove the empty __init__.py files
    # _log.debug("remove temporary __init__.py files")
    for d in pydirs:
        init_py = d / "__init__.py"
        init_py.unlink()
    # remove build dir, and restore any moved build dir
    shutil.rmtree(str(build_dir))
    if moved_build_dir is not None:
        _log.debug(f"restore build dir '{build_dir}' from '{moved_build_dir}'")
        os.rename(moved_build_dir, str(build_dir))
    # restore previous args
    sys.argv = saved_args
    # change back to previous directory
    os.chdir(orig_dir)


def find_python_directories(target_dir: Path) -> List[Path]:
    """Find all directories from target_dir, on down, that contain a
    Python module or sub-package.
    """
    # get directories that contain python files -> pydirs
    pydirs = set((x.parent for x in target_dir.rglob("*.py")))
    # get all directories in the tree leading to the 'pydirs'
    alldirs = set()
    for d in pydirs:
        while d != target_dir:
            alldirs.add(d)
            d = d.parent
    return list(alldirs)
