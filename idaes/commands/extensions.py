#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""Commandline Utilities for Managing the IDAES Data Directory"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

# TODO: protected access issues
# pylint: disable=protected-access

__author__ = "John Eslick"

import os
import sys
import logging
import subprocess

import click

import idaes
import idaes.commands.util.download_bin
from idaes.commands import cb
from pyomo.common.dependencies import attempt_import

_log = logging.getLogger("idaes.commands.extensions")


def print_extensions_version(library_only=False):
    click.echo("---------------------------------------------------")
    click.echo("IDAES Extensions Build Versions")
    click.echo("===================================================")
    if not library_only:
        v = os.path.join(idaes.bin_directory, "version_solvers.txt")
        try:
            with open(v, "r") as f:
                v = f.readline().strip()
        except FileNotFoundError:
            v = "no version file found"
        click.echo("Solvers:  v{}".format(v))
    v = os.path.join(idaes.bin_directory, "version_lib.txt")
    try:
        with open(v, "r") as f:
            v = f.readline().strip()
    except FileNotFoundError:
        v = "no version file found"
    click.echo("Library:  v{}".format(v))
    click.echo("===================================================")
    return 0


def print_license():
    click.echo("---------------------------------------------------")
    click.echo("IDAES Extensions License Information")
    click.echo("===================================================")
    fpath = os.path.join(idaes.bin_directory, "license.txt")
    try:
        with open(fpath, "r") as f:
            for line in f.readlines():
                click.echo(line.strip())
    except FileNotFoundError:
        click.echo("no license file found")
    click.echo("")
    click.echo("===================================================")
    return 0


@cb.command(name="get-extensions", help="Get solvers and libraries")
@click.option(
    "--release",
    help="Optional, specify an official binary release to download",
    default=None,
)
@click.option(
    "--url", help="Optional, URL to download solvers/libraries from", default=None
)
@click.option(
    "--distro", help="OS or Linux distribution (default=auto)", default="auto"
)
@click.option("--insecure", is_flag=True, help="Don't verify download location")
@click.option(
    "--cacert",
    help="Specify certificate file to verify download location",
    default=None,
)
@click.option("--nochecksum", is_flag=True, help="Don't verify the file checksum")
@click.option(
    "--library-only",
    is_flag=True,
    help="Only install shared physical property function libraries, not solvers",
)
@click.option(
    "--no-download",
    is_flag=True,
    help="Don't download anything, but report what would be done",
)
@click.option("--extra", multiple=True, help="Install extras")
@click.option("--extras-only", is_flag=True, help="Only install extras")
@click.option("--to", default=None, help="Put extensions in a alternate location")
@click.option("--verbose", help="Show details", is_flag=True)
def get_extensions(
    release,
    url,
    insecure,
    cacert,
    verbose,
    distro,
    nochecksum,
    library_only,
    no_download,
    extras_only,
    extra,
    to,
):
    if url is None and release is None:
        # the default release is only used if neither a release or url is given
        release = idaes.config.default_binary_release
    if url is not None and release is not None:
        click.echo("\n* You must provide either a release or url not both.")
    elif url is not None or release is not None:
        click.echo("Getting files...")
        try:
            d = idaes.commands.util.download_bin.download_binaries(
                release,
                url,
                insecure,
                cacert,
                verbose,
                distro,
                nochecksum,
                library_only,
                no_download,
                extras_only,
                extra,
                alt_path=to,
            )
            click.echo("Done")
        except idaes.commands.util.download_bin.UnsupportedPlatformError as e:
            click.echo("")
            click.echo(e)
            click.echo("")
            click.echo("Specify an os with --distro <os>:")
            return
        if no_download:
            for k, i in d.items():
                click.echo(f"{k:14}: {i}")
        else:
            print_extensions_version(library_only)
    else:
        click.echo("\n* You must provide a download URL for IDAES binary files.")


@cb.command(name="hash-extensions", help="Calculate release hashes")
@click.option(
    "--release",
    help="Optional, specify an official binary release to download",
    default=None,
    required=True,
)
@click.option("--path", help="Directory of release files", default="./")
def hash_extensions(release, path):
    hfile = f"sha256sum_{release}.txt"
    if path is not None:
        hfile = os.path.join(path, hfile)

    def _write_hash(fp, pack, plat):
        f = f"idaes-{pack}-{plat}.tar.gz"
        if path is not None:
            h = idaes.commands.util.download_bin.hash_file_sha256(os.path.join(path, f))
        else:
            h = idaes.commands.util.download_bin.hash_file_sha256(f)
        fp.write(h)
        fp.write("  ")
        fp.write(f)
        fp.write("\n")

    with open(hfile, "w") as f:
        for plat in idaes.config.base_platforms:
            for pack in ["solvers", "lib"]:
                _write_hash(f, pack, plat)
        for plat in idaes.config.base_platforms:
            for pack, sp in idaes.config.extra_binaries.items():
                if plat not in sp:
                    continue
                _write_hash(f, pack, plat)


@cb.command(name="bin-platform", help="Show the compatible binary build.")
@click.option("--distro", default="auto")
def bin_platform(distro):
    fd, _ = idaes.commands.util.download_bin._get_file_downloader(False, None)
    try:
        _, platform = idaes.commands.util.download_bin._get_arch_and_platform(
            fd, distro
        )
        click.echo(idaes.commands.util.download_bin._get_release_platform(platform))
    except idaes.commands.util.download_bin.UnsupportedPlatformError:
        click.echo(f"No supported binaries found for {platform}.")


@cb.command(name="extensions-license", help="show license info for binary extensions")
def extensions_license():
    print_license()


@cb.command(name="extensions-version", help="show license info for binary extensions")
def extensions_version():
    print_extensions_version()


@cb.command(
    name="install-cyipopt",
    help="[BETA] Install CyIpopt and link to IDAES IPOPT library",
)
def install_cyipopt():
    _, cyipopt_available = attempt_import("cyipopt")
    click.echo(
        "WARNING: The install-cyipopt command is beta functionality and may not work"
        " on all platforms."
    )
    if cyipopt_available:
        click.echo(
            "CyIpopt is already available in the current Python environment."
            " Please uninstall CyIpopt before running this command. Note that"
            " you may need to clear your pip cache to re-install with IDAES"
            " binaries."
        )
        return
    pkgconfig_path = os.pathsep.join(
        [
            # Prepend IDAES's pkgconfig directory to this path so we always
            # build cyipopt against IDAES binaries (if they exist) when using this
            # command.
            os.path.join(idaes.bin_directory, "lib", "pkgconfig"),
            os.getenv("PKG_CONFIG_PATH", ""),
        ]
    )
    subprocess_environ = dict(os.environ, PKG_CONFIG_PATH=pkgconfig_path)
    # TODO: Possibly set IPOPTWINDIR on Windows? I haven't been able to get this to
    # work on GHA, so I'd like some advice from somebody who has built CyIpopt
    # on Windows before staring to make assumptions. -RBP
    ret = subprocess.run(["pip", "install", "cyipopt"], env=subprocess_environ)
    if ret.returncode == 1:
        # CyIpopt wheels don't build on Python 3.9 (see
        # https://github.com/mechmotum/cyipopt/issues/225), so we have this workaround
        # in place as an alternative. If #225 gets resolved or we stop supporting 3.9,
        # and no other cases come up where `pip install cyipopt` fails, we can remove
        # this code. -RBP
        click.echo(
            "WARNING: Command `pip install cyipopt` returned 1. Attempting to install"
            " from source"
        )
        orig_cwd = os.getcwd()
        try:
            os.chdir(idaes.bin_directory)
            cyipopt_dir = os.path.join(idaes.bin_directory, "cyipopt")
            if os.path.exists(cyipopt_dir):
                raise RuntimeError(
                    f"{cyipopt_dir} already exists. Please remove and try again."
                )
            subprocess.run(
                [
                    "git",
                    "clone",
                    "https://github.com/mechmotum/cyipopt.git",
                    "--branch=v1.4.1",
                ]
            )
            os.chdir(cyipopt_dir)
            ret = subprocess.run(
                ["python", "setup.py", "develop"], env=subprocess_environ
            )
        finally:
            os.chdir(orig_cwd)
        if ret.returncode == 1:
            raise RuntimeError("Error installing CyIpopt from source.")
    click.echo("Installed CyIpopt in the current Python environment.")
    libdir = os.path.join(idaes.bin_directory, "lib")
    if sys.platform == "nt":
        lib_envvar = "PATH"
    elif sys.platform == "darwin":
        lib_envvar = "DYLD_LIBRARY_PATH"
    else:
        lib_envvar = "LD_LIBRARY_PATH"
    homedir = os.getenv("HOME", "")
    bashrc_file = os.path.join(homedir, ".bashrc")
    zshrc_file = os.path.join(homedir, ".zshrc")
    click.echo(
        f"Note: To use CyIpopt, you may need to add {libdir} to your {lib_envvar}"
        " environment variable. That is, add the following line to your shell"
        f" configuration file (e.g. {bashrc_file} for Bash or {zshrc_file} for Zsh):"
        # Or $PROFILE for PowerShell?
    )
    export_rhs = os.pathsep.join(["$" + lib_envvar, libdir])
    click.echo(f"\n    export {lib_envvar}={export_rhs}\n")
