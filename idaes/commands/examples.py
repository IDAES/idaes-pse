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
from setuptools import setup, find_packages
import shutil
import sys
from tempfile import mkdtemp
from typing import List
from uuid import uuid4
from zipfile import ZipFile

# third-party
import click
from nbconvert.exporters import NotebookExporter
from nbconvert.writers import FilesWriter
from traitlets.config import Config
import nbformat
import requests

# package
from idaes.commands import cb
from idaes.ver import package_version as V

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
JUPYTER_NB_VERSION = 4  # for parsing
REMOVE_CELL_TAG = "remove_cell"  # for stripping Jupyter notebook cells
STRIPPED_NOTEBOOK_SUFFIX = "_s"

# Global vars

g_tempdir, g_egg = None, None

# Exceptions

class DownloadError(Exception):
    """Used for errors downloading the release files.
    """

    pass


class InstallError(Exception):
    """Used for errors installing the source as a Python package.
    """

    pass


Release = namedtuple("Release", ["date", "tag", "info"])


@cb.command(
    name="get-examples", help="Fetch example scripts and Jupyter Notebooks."
)
@click.option(
    "--dir", "-d", "directory", help="installation target directory", default="examples",
    type=str,
)
@click.option(
    "--no-install", "-I", "no_install", help="Do *not* install examples into 'idaes_examples' package",
    is_flag=True
)
@click.option(
    "--list",
    "-l",
    "list_releases",
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
    default=None,
    show_default=False,
)
def get_examples(directory, no_install, list_releases, no_download, version,
                 unstable):
    """Get the examples from Github and put them in a local directory.
    """
    # will always need a list of releases
    releases = get_releases(unstable)
    # list-releases mode
    if list_releases:
        print_releases(releases, unstable)
        sys.exit(0)
    # otherwise..
    target_dir = Path(directory)
    # do nothing?
    if no_download and no_install:
        click.echo("Download and installation disabled. Done.")
        sys.exit(0)
    # download?
    if no_download:  # no download
        _log.info("skipping download step (use existing examples)")
        click.echo("Skip download")
    else:  # download
        # check that unstable flag is given before downloading unstable version
        if version is not None:
            stable_ver = re.match(r".*-\w+$", version) is None
            if not stable_ver and not unstable:
                click.echo(f"Cannot download unstable version {version} unless you add "
                           f"the -U/--unstable flag")
                sys.exit(-1)
        # set version
        if version is None:
            ex_version = get_examples_version(PKG_VERSION)
        else:
            ex_version = version
        # give an error if selected version does not exist
        if ex_version not in [r.tag for r in releases]:
            if version is None:
                click.echo(f"Internal Error: Could not find an examples release matching IDAES {ex_version}.\n"
                           f"As a workaround, you can manually pick a version with -V/--version.\n"
                           f"Use -l/--list-releases to see all available versions.")
            else:
                click.echo(f"Could not find an examples release matching IDAES {version}.\n"
                           f"Use -l/--list-releases to see all available versions.")
            sys.exit(-1)
        click.echo("Downloading...")
        try:
            download(target_dir, ex_version)
        except DownloadError as err:
            _log.warning(f"abort due to failed download: {err}")
            clean_up_temporary_files()
            sys.exit(-1)
        full_dir = os.path.realpath(target_dir)
        _log.info(f"downloaded examples to: '{full_dir}'")
    # install
    if no_install:
        _log.info("skipping installation step")
        click.echo("Skip install")
    else:
        if not target_dir.exists():
            if no_download:
                click.echo(f"Target directory '{target_dir}' does not exist")
                sys.exit(-1)
            else:
                click.echo(f"Internal error: After download, directory '{target_dir}' "
                           f"does not exist")
                sys.exit(-1)
        click.echo("Installing...")
        try:
            install_src(ex_version, target_dir)
        except InstallError as err:
            click.echo(f"Install error: {err}")
            clean_up_temporary_files()
            sys.exit(-1)
        _log.info(f"Installed examples as package {INSTALL_PKG}")

    click.echo("Cleaning up...")
    # strip notebooks
    strip_test_cells(target_dir)
    # temporary files
    clean_up_temporary_files()

    # Done
    print_summary( ex_version, target_dir, not no_install)


def print_summary(version, dirname, installed):
    sep1, sep2 = "-" * 40, "=" * 40
    print(f"{sep1}\nIDAES Examples {version}\n{sep2}")
    print(f"Path   : {dirname}")
    if installed:
        print(f"Package: {INSTALL_PKG}")
    else:
        print(f"Package: not installed")
    print(sep2)


def get_examples_version(idaes_version: str):
    """Given the specified 'idaes-pse' repository release version,
    figure out the matching 'examples-pse' repository release version.

    Args:
        idaes_version: IDAES version, e.g. "1.6.0"

    Returns:
        Examples version, or if there is no match, return None.
    """
    # TODO: fetch mapping from Github and apply mapping
    return idaes_version


def download(target_dir: Path, version: str):
    """Download `version` into `target_dir`.
    Raises:
        DownloadError
    """
    # check target directory
    if target_dir.exists():
        click.echo(f"Directory '{target_dir}' exists. Please move or delete and "
                   f"try this command again")
        raise DownloadError("directory exists")
    # download
    try:
        download_contents(target_dir, version)
    except DownloadError as err:
        click.echo(f"Download failed: {err}")
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


def download_contents(target_dir, version):
    """Download the given version from the Github releases and make
    its `REPO_DIR` subdirectory be the `target_dir`.

    Raises:
        DownloadError: if the GET on the release URL returns non-200 status
    """
    global g_tempdir
    url = archive_file_url(version)
    _log.info(f"get examples from: {url}")
    # stream out to a big .zip file
    req = requests.get(url, stream=True)
    if req.status_code != 200:
        if req.status_code in (400, 404):
            raise DownloadError(f"file not found")
        raise DownloadError(f"status={req.status_code}")
    # note: mkdtemp() creates a directory that seems un-removable on Windows.
    # So, instead, just create the directory yourself in the current directory
    random_name = str(uuid4())
    try:
        os.mkdir(random_name)
    except Exception as err:
        _log.fatal(f"making directory '{random_name}': {err}")
        click.echo("Cannot make temporary directory in current directory. Abort.")
        sys.exit(-1)
    g_tempdir = tempdir = Path(random_name)
    _log.debug(f"created temporary directory '{tempdir.name}'")
    tempfile = tempdir / "examples.zip"
    with tempfile.open("wb") as f:
        for chunk in req.iter_content(chunk_size=65536):
            f.write(chunk)
    _log.info(f"downloaded zipfile to {tempfile}")
    # open as a zip file, and extract all files into the temporary directory
    _log.debug(f"open zip file: {tempfile}")
    zipf = ZipFile(str(tempfile))
    zipf.extractall(path=tempdir.name)
    # move the REPO_DIR subdirectory into the target dir
    subdir = Path(tempdir.name) / f"{REPO_NAME}-{version}" / REPO_DIR
    _log.debug(f"move {subdir} -> {target_dir}")
    os.rename(str(subdir), str(target_dir))
    zipf.close()


def clean_up_temporary_files():
    # temporary directory created for unzipping and renaming
    if g_tempdir:
        d = g_tempdir
        _log.debug(f"remove temporary directory.start name='{d}'")
        try:
            shutil.rmtree(d)
        except Exception as err:  # WTF
            _log.warning(f"remove temporary directory.error name='{d}' msg='{err}'")
        else:
            _log.debug(f"removed temporary directory.end name='{d}'")
    # egg file created by setuptools
    if g_egg:
        _log.debug(f"remove setuptools egg path='{g_egg}'")
        shutil.rmtree(g_egg.name)


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
    When done, name the directory back to 'src', and remove '__init__.py' files.
    Then clean up whatever cruft is left behind..
    """
    global g_egg
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
    # save name of egg file, for later cleanup
    f = Path(".") / (INSTALL_PKG + ".egg-info")
    if f.exists():
        g_egg = f
    else:
        _log.warning(f"build dir not found path='{f}'")


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


def find_notebook_files(target_dir: Path) -> List[Path]:
    """Find all files ending in ".ipynb", at or below target_dir.
    """
    return list(target_dir.rglob("*.ipynb"))


def strip_test_cells(target_dir: Path):
    """Strip test cells from notebooks below `target_dir`.
    See `strip_tags()` for how that is done.
    """
    nb_files = find_notebook_files(target_dir)
    for nb_file in nb_files:
        if has_tagged_cells(nb_file):
            strip_tags(nb_file)


def has_tagged_cells(nb: Path):
    """Quickly check whether this notebook has any cells with the "special" tag.

    Returns:
        True = yes, it does; False = no specially tagged cells
    Raises:
        NotebookFormatError, if notebook at 'entry' is not parseable
    """
    # parse the notebook (assuming this is fast; otherwise should cache it)
    try:
        nb = nbformat.read(str(nb), as_version=JUPYTER_NB_VERSION)
    except nbformat.reader.NotJSONError:
        raise ValueError(f"Notebook '{nb}' is not valid JSON")
    # look for tagged cells; return immediately if one is found
    for i, c in enumerate(nb.cells):
        if "tags" in c.metadata and REMOVE_CELL_TAG in c.metadata.tags:
            _log.debug(f"Found {REMOVE_CELL_TAG} tag in cell {i}")
            return True  # can stop now, one is enough
    # no tagged cells
    return False


def strip_tags(nb: Path):
    """Strip tags from notebook, if there are any.

    Behavior depends on filename:

    - ends with "_test.ipynb" or "_testing.ipynb" -> make a stripped file
      with the suffix removed.
    - otherwise -> make a stripped file with "_s.ipynb" at the end.

    Raises:
        IOError: If the desired target for the stripped notebook already exists
    """
    _log.info(f"stripping tags from notebook {nb}")
    # Set up configuration for removing specially tagged cells
    conf = Config()
    conf.TagRemovePreprocessor.remove_cell_tags = (REMOVE_CELL_TAG,)
    conf.NotebookExporter.preprocessors = [
        "nbconvert.preprocessors.TagRemovePreprocessor"
    ]  # for some reason, this only works with the full module path
    # Convert from Notebook format to Notebook format, stripping tags
    (body, resources) = NotebookExporter(config=conf).from_filename(str(nb))
    # Set up output destination
    wrt = FilesWriter()
    wrt.build_directory = str(nb.parent)
    nb_name = nb.stem
    # Logic: (1) ends in _test.ipynb or _testing.ipynb? -> remove _test/_testing
    # (2) otherwise -> add "_s" (for stripped) before .ipynb extension
    if nb_name.endswith("_test"):
        stripped_name = nb_name[:-5]
    elif nb_name.endswith("_testing"):
        stripped_name = nb_name[:-8]
    else:
        stripped_name = nb_name + STRIPPED_NOTEBOOK_SUFFIX
    new_nb = nb.parent / (stripped_name + ".ipynb")
    # Don't overwrite an existing file
    if new_nb.exists():
        raise IOError(f"Cannot create stripped notebook '{new_nb}': file exists")
    # Write stripped notebook to output file
    wrt.write(body, resources, notebook_name=stripped_name)
