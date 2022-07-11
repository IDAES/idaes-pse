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
# Pyomo utility for delayed import
from pyomo.common.dependencies import attempt_import

# stdlib
from collections import namedtuple
from datetime import datetime
from io import StringIO
import logging
from operator import attrgetter
import os
import re
from pathlib import Path
import shutil
import sys
from typing import List
from uuid import uuid4
from zipfile import ZipFile
import json


# third-party
import click

# third-party slow
nb_exporters = attempt_import("nbconvert.exporters")[0]
nb_writers = attempt_import("nbconvert.writers")[0]
traitlets_config = attempt_import("traitlets.config")[0]
nbformat = attempt_import("nbformat")[0]
requests = attempt_import("requests")[0]

# package
from idaes.commands import cb
from idaes.commands.base import how_to_report_an_error
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
SOLUTION_CELL_TAG = "solution"  # ditto
EXERCISE_CELL_TAG = "exercise"  # ditto

STRIPPED_NOTEBOOK_SUFFIX = "_s"

# Global vars

g_tempdir, g_egg = None, None

# Exceptions


class DownloadError(Exception):
    """Used for errors downloading the release files."""

    pass


class CopyError(Exception):
    """Used for errors copying files."""

    pass


class InstallError(Exception):
    """Used for errors installing the source as a Python package."""

    pass


class GithubError(Exception):
    pass


Release = namedtuple("Release", ["date", "tag", "info"])


@cb.command(name="get-examples", help="Fetch example scripts and Jupyter Notebooks.")
@click.option(
    "--dir",
    "-d",
    "directory",
    help="installation target directory",
    default="examples",
    type=str,
)
@click.option(
    "--local",
    "local_dir",
    help="For developers: instead of downloading, copy from an "
    "idaes-examples repository on local disk",
)
@click.option(
    "--no-install",
    "-I",
    "no_install",
    help="Do *not* install examples into 'idaes_examples' package",
    is_flag=True,
)
@click.option(
    "--list",
    "-l",
    "list_releases",
    help="List all available released versions, and stop",
    is_flag=True,
)
@click.option(
    "--no-download", "-N", "no_download", help="Do not download anything", is_flag=True
)
@click.option(
    "--unstable",
    "-U",
    help="Allow and list unstable/pre-release versions",
    is_flag=True,
)
@click.option(
    "--version",
    "-V",
    help=f"Version of examples to download",
    default=None,
    show_default=False,
)
def get_examples(
    directory, no_install, list_releases, no_download, version, unstable, local_dir
):
    """Get the examples from Github and put them in a local directory."""
    # list-releases mode
    if list_releases:
        try:
            releases = get_releases(unstable)
        except GithubError as err:
            click.echo(f"Error getting data from Github: {err}")
            sys.exit(-1)
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
        if target_dir.exists():
            click.echo(
                f"Target directory '{target_dir}' already exists. Please "
                f"remove it, or choose a different directory."
            )
            sys.exit(-1)
        if local_dir is not None:
            click.echo(f"Copying from local directory: {local_dir}")
            local_path = Path(local_dir)
            if not local_path.exists() or not local_path.is_dir():
                click.echo(
                    f"Cannot copy from local directory '{local_dir}': "
                    f"directory does not exist, or not a directory"
                )
                sys.exit(-1)
            try:
                copy_contents(target_dir, local_path)
            except CopyError as err:
                click.echo(
                    f"Failed to copy from '{local_dir}' to '{target_dir}': " f"{err}"
                )
                sys.exit(-1)
            ex_version = version if version else PKG_VERSION
        else:
            try:
                releases = get_releases(unstable)
            except GithubError as err:
                click.echo(f"Error getting data from Github: {err}")
                sys.exit(-1)
            # check that unstable flag is given before downloading unstable version
            if version is not None:
                stable_ver = re.match(r".*-\w+$", version) is None
                if not stable_ver and not unstable:
                    click.echo(
                        f"Cannot download unstable version {version} unless you add "
                        f"the -U/--unstable flag"
                    )
                    sys.exit(-1)
            # set version
            if version is None:
                ex_version = get_examples_version(PKG_VERSION)
            else:
                ex_version = version
            # give an error if selected version does not exist
            if ex_version not in [r.tag for r in releases]:
                if version is None:
                    click.echo(
                        f"Internal Error: Could not find an examples release\n"
                        f"matching IDAES version {PKG_VERSION}\n"
                        f"You can manually pick  version with '-V/--version'\n"
                        f"or install from a local directory with '--local'.\n"
                        f"Use '-l/--list-releases' to see all versions.\n"
                        f"{how_to_report_an_error()}\n"
                    )
                else:
                    click.echo(
                        f"Could not find an examples release matching IDAES version "
                        f"{version}.\n Use -l/--list-releases to see all "
                        f"available versions."
                    )
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
                click.echo(
                    f"Internal error: After download, directory '{target_dir}'\n"
                    f"does not exist.\n"
                    f"{how_to_report_an_error()}"
                )
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
    _log.info("Stripping test/solution cells from notebooks")
    strip_special_cells(target_dir)
    # temporary files
    _log.info("Removing temporary files")
    clean_up_temporary_files()

    # Done
    print_summary(ex_version, target_dir, not no_install)


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
    identify the matching 'examples-pse' repository release version.

    Args:
        idaes_version: IDAES version, e.g. "1.5.0" or "1.5.0.dev0+e1bbb[...]"

    Returns:
        Examples version, or if there is no match, return None.
    """
    # Fetch the idaes:examples version mapping from Github
    compat_file = "idaes-compatibility.json"
    url = f"{GITHUB_API}/repos/{REPO_ORG}/{REPO_NAME}/contents/{compat_file}"
    headers = {"Accept": "application/vnd.github.v3.raw"}
    _log.debug(f"About to call requests.get({url}, {headers})")
    res = requests.get(url, headers=headers)
    if not res.ok:
        _log.debug(f"Problem getting mapping file: {res.json()}")
        raise DownloadError(res.json())

    try:
        compat_mapping = json.loads(res.text)["mapping"]
    except KeyError:
        # return the latest version instead
        _log.warning("Ill-formed compatibility mapping file for examples repository:")
        _log.debug(f"compat_mapping: {res.text}")
        _log.info("Defaulting to latest released version of examples.")
        return None

    idaes_version_num = idaes_version
    version_numbers = idaes_version.split(".")
    if len(version_numbers) > 3:
        idaes_version_num = ".".join(version_numbers[:3])
        click.echo(
            f"Warning: non-release version of IDAES detected. "
            f"Using IDAES {idaes_version_num} as reference; "
            f"examples version compatibility is not guaranteed."
        )

    try:
        examples_version = compat_mapping[idaes_version_num]
    except KeyError:
        # return the latest version instead, as above
        _log.warning(
            "IDAES version not found in compatibility mapping file. \
                Defaulting to latest released version of examples."
        )
        return None

    _log.debug(f"get_examples_version({idaes_version}: {examples_version}")

    return examples_version


def download(target_dir: Path, version: str):
    """Download `version` into `target_dir`.
    Raises:
        DownloadError
    """
    # check target directory
    if target_dir.exists():
        click.echo(
            f"Directory '{target_dir}' exists. Please move or delete and "
            f"try this command again"
        )
        raise DownloadError("directory exists")
    # download
    try:
        download_contents(target_dir, version)
    except DownloadError as err:
        click.echo(f"Download failed: {err}")
        raise


def is_illegal_dir(d: Path):
    """Refuse to remove directories for some situations, for safety."""
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
    shutil.move(str(subdir), str(target_dir))
    zipf.close()


def copy_contents(target_dir, repo_root):
    subdir = repo_root / REPO_DIR
    if not subdir.is_dir():
        raise CopyError(f"Could not copy from '{subdir}': not a directory")
    _log.info(f"copy.local.start from={subdir} to={target_dir}")
    try:
        shutil.copytree(subdir, target_dir)
    except shutil.Error as err:
        raise CopyError(err)
    except FileNotFoundError:
        raise CopyError(f"Could not find file '{subdir}'")
    except Exception as err:
        raise CopyError(f"Unknown problem copying: {err}")
    _log.info(f"copy.local.end from={subdir} to={target_dir}")


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
    if g_egg and g_egg.exists() and g_egg.is_dir():
        _log.debug(f"remove setuptools egg path='{g_egg}'")
        try:
            shutil.rmtree(g_egg.name)
        except Exception as err:
            _log.warning(f"remove temporary directory.error name='{g_egg}' msg='{err}'")
    # dist directory created by setuptools
    d = Path("dist")
    if d.exists() and d.is_dir():
        for f in d.glob("idaes_examples-*.egg"):
            try:
                f.unlink()
            except Exception as err:
                _log.warning(f"could not remove distribution file {f}: {err}")
        # remove directory, if now empty
        num_files = len(list(d.glob("*")))
        if num_files == 0:
            _log.info(f"removing dist directory '{d.absolute()}'")
            try:
                d.rmdir()
            except Exception as err:
                _log.warning(f"could not remove distribution directory {d}: {err}")


def archive_file_url(version, org=REPO_ORG, repo=REPO_NAME):
    """Build & return URL for a given release version."""
    return f"{GITHUB}/{org}/{repo}/archive/{version}.zip"


def get_releases(unstable) -> List[Release]:
    """Returns a list of releases.

    The list is sorted in ascending order by date.
    """
    releases = []
    url = f"{GITHUB_API}/repos/{REPO_ORG}/{REPO_NAME}/releases"
    resp = requests.get(url)
    data, headers = resp.json(), resp.headers
    check_github_response(data, headers)
    for rel in data:
        if not unstable and rel["prerelease"]:
            continue
        releases.append(Release(rel["published_at"], rel["tag_name"], rel["name"]))
    releases.sort(key=attrgetter("date"))  # sort by publication date
    return releases


def check_github_response(data, headers):
    """Check whether GitHub gave an error message. If so, raise a GithubError
    with a hopefully informative and useful message.
    """
    if isinstance(data, list):
        return  # lists are assumed to be the releases
    if isinstance(data, dict) and "message" in data:
        if "rate limit exceeded" in data["message"]:
            reset_ts = int(headers["X-RateLimit-Reset"])
            now = datetime.now()
            now_ts, tzinfo = now.timestamp(), now.astimezone().tzinfo
            wait_min = int((reset_ts - now_ts) // 60) + 1
            reset_dt = datetime.fromtimestamp(reset_ts)
            datestr = reset_dt.astimezone(tzinfo).isoformat()
            raise GithubError(
                f"API rate limit exceeded.\n"
                f"You will need to wait {wait_min} minutes,"
                f" until {datestr}, to try again from "
                f"this computer."
            )
        else:
            raise GithubError(f"Error connecting to Github: {data['message']}")
    else:
        raise GithubError(f"Invalid result from Github: data={data} headers={headers}")


def print_releases(releases: List[Release], unstable):
    """Print the releases, as returned by `get_releases()`, as a table
    to standard output.
    """
    if len(releases) == 0:
        if unstable:
            print("No releases found")
        else:
            print(
                "No stable releases found. Add -U/--unstable to also look "
                "for pre-releases."
            )
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
    from setuptools import setup, find_packages  # import here due to slowness

    global g_egg
    orig_dir = Path(os.curdir).absolute()
    target_dir = Path(target_dir.absolute())
    root_dir = target_dir.parent
    examples_dir = root_dir.absolute() / INSTALL_PKG
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
    try:
        for d in pydirs:
            init_py = d / "__init__.py"
            init_py.open("w")
    except IOError as err:
        raise InstallError(f"error writing temporary __init__.py files: {err}")
    # temporarily rename target directory to the package name
    _log.info(f"rename {target_dir} -> {examples_dir}")
    shutil.move(target_dir, examples_dir)
    # if there is a 'build' directory, move it aside
    build_dir = root_dir / "build"
    if build_dir.exists():
        from uuid import uuid1

        random_letters = str(uuid1())
        moved_build_dir = f"{build_dir}.{random_letters}"
        _log.debug(f"move existing build dir to {moved_build_dir}")
        shutil.move(str(build_dir), moved_build_dir)
    else:
        _log.debug("no existing build directory (nothing to do)")
        moved_build_dir = None
    # run setuptools' setup command (in root directory)
    _log.info(f"run setup command in directory {root_dir}")
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
        zip_safe=False,
    )
    # restore stdout
    sys.stdout = orig_stdout
    # print/log output
    output_str = setup_out.getvalue()
    if _log.isEnabledFor(logging.DEBUG):
        for line in output_str.split("\n"):
            _log.debug(f"(setup) {line}")
    # name the target directory back to original
    _log.info(f"rename '{examples_dir}' to  '{target_dir}'")
    shutil.move(examples_dir, target_dir)
    # remove the empty __init__.py files
    _log.info("remove temporary __init__.py files")
    for d in pydirs:
        init_py = d / "__init__.py"
        _log.debug(f"remove '{init_py}")
        init_py.unlink()
    # remove build dir, and restore any moved build dir
    _log.info(f"remove build directory '{build_dir}'")
    try:
        shutil.rmtree(build_dir)
    except Exception as err:
        _log.warning(f"failed to remove build directory {build_dir}: {err}")
    if moved_build_dir is not None:
        _log.info(f"restore build dir '{build_dir}' from '{moved_build_dir}'")
        shutil.move(moved_build_dir, build_dir)
    # restore previous args
    sys.argv = saved_args
    # change back to previous directory
    os.chdir(orig_dir)
    # save name of egg file, for later cleanup
    f = root_dir / (INSTALL_PKG + ".egg-info")
    if f.exists():
        g_egg = f
    else:
        _log.warning(f"egg-info file not found path='{f}'")


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
    """Find all files ending in ".ipynb", at or below target_dir."""
    return list(target_dir.rglob("*.ipynb"))


def strip_special_cells(target_dir: Path):
    """Strip 'special' cells from notebooks below `target_dir`.
    See `strip_tags()` for how that is done.
    """
    nb_files = find_notebook_files(target_dir)
    for nb_file in nb_files:
        if strip_tags(nb_file):
            _log.info(f"removing original file '{nb_file}'")
            nb_file.unlink()


def strip_tags(nb: Path) -> bool:
    """Strip tags from notebook, if there are any.

    Behavior depends on filename (case-insensitive).

    - if the notebook file name ends with "_solution_testing", create two files:

       1. chop off "_testing" from name, remove test and exercise cells
       2. chop off "_solution_testing" from name & add "_exercise", remove solution and test cells

    - if the file ends with "_testing" only, create one file:

        1. chop off "_testing" from name, remove test cells

    - otherwise: do nothing

    Returns:
        Was anything stripped (and new files created)?

    Raises:
        IOError: If the desired target for the stripped notebook already exists
    """
    _log.info(f"stripping tags from notebook {nb}")
    nb_name = nb.stem

    # Figure out which tags to strip:
    # - check for case (1), solution_testing.ipynb
    if nb_name.lower().endswith("_solution_testing"):
        pre = nb_name[:-17]
        names = [f"{pre}_solution.ipynb", f"{pre}_exercise.ipynb"]
        tags_to_strip = [
            (REMOVE_CELL_TAG, EXERCISE_CELL_TAG),
            (REMOVE_CELL_TAG, SOLUTION_CELL_TAG),
        ]
    # - check for case (2), _testing.ipynb
    elif nb_name.lower().endswith("_testing"):
        pre = nb_name[:-8]
        names = [f"{pre}.ipynb"]
        tags_to_strip = [(REMOVE_CELL_TAG,)]
    # - if neither, we are done here
    else:
        _log.debug(f"notebook '{nb}' does not need to have tags stripped")
        return False

    # Create all desired tag-stripped copies
    for name, remove_tags in zip(names, tags_to_strip):
        target_nb = nb.parent / name
        _log.info(f"creating new notebook '{target_nb}'")
        # Don't overwrite an existing file
        if target_nb.exists():
            _log.warning(f"cannot create new notebook '{target_nb}': file exists")
            continue
        # Set up configuration for removing specially tagged cells
        conf = traitlets_config.Config()
        conf.TagRemovePreprocessor.remove_cell_tags = remove_tags
        conf.NotebookExporter.preprocessors = [
            # this requires the full module path
            "nbconvert.preprocessors.TagRemovePreprocessor"
        ]
        # Convert from Notebook format to Notebook format, stripping tags
        (body, resources) = nb_exporters.NotebookExporter(config=conf).from_filename(
            str(nb)
        )
        # Set up output destination
        wrt = nb_writers.FilesWriter()
        wrt.build_directory = str(target_nb.parent)
        # Write stripped notebook to output file
        wrt.write(body, resources, notebook_name=target_nb.stem)

    return True


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
