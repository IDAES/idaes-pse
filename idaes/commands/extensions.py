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
"""Commandline Utilities for Managing the IDAES Data Directory"""

__author__ = "John Eslick"

import os
import click
import logging
import idaes
import idaes.util.download_bin
from idaes.commands import cb

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


@cb.command(name="get-extensions", help="Get solvers and libraries")
@click.option(
    "--release",
    help="Optional, specify an official binary release to download",
    default=None)
@click.option(
    "--url",
    help="Optional, URL to download solvers/libraries from",
    default=None)
@click.option(
    "--platform",
    help="Platform to download binaries for. Use 'idaes get-extensions-platforms'"
         " for a list of available platforms (default=auto)",
    default="auto")
@click.option(
    "--insecure",
    is_flag=True,
    help="Don't verify download location")
@click.option(
    "--cacert",
    help="Specify certificate file to verify download location",
    default=None)
@click.option(
    "--nochecksum",
    is_flag=True,
    help="Don't verify the file checksum")
@click.option(
    "--library-only",
    is_flag=True,
    help="Only install shared physical property function libraries, not solvers")
@click.option(
    "--no-download",
    is_flag=True,
    help="Don't download anything, but report what would be done")
@click.option(
    "--show-current-version",
    is_flag=True,
    help="Show the version information if any for the currently installed")
@click.option(
    "--show-platforms",
    is_flag=True,
    help="Show the platform options")
@click.option(
    "--show-extras",
    is_flag=True,
    help="Show list of binary extras")
@click.option(
    "--extra",
    multiple=True,
    help="Install extras")
@click.option(
    "--extras-only",
    is_flag=True,
    help="Only install extras")
@click.option(
    "--to",
    default=None,
    help="Put extensions in a alternate location")
@click.option("--verbose", help="Show details", is_flag=True)
def get_extensions(
    release,
    url,
    insecure,
    cacert,
    verbose,
    platform,
    nochecksum,
    library_only,
    no_download,
    show_current_version,
    show_platforms,
    show_extras,
    extras_only,
    extra,
    to):

    if show_platforms:
        click.echo("\nSupported platforms for IDAES binary extensions.")
        for key, mes in sorted(idaes.config.known_binary_platform.items()):
            click.echo(f"    {key}: {mes}")
        click.echo("\nBinaries compiled on these platforms:")
        for k in sorted(idaes.config.basic_platforms()):
            click.echo(f"    {k}")
        return
    elif show_extras:
        click.echo("\nBinary Extras")
        for k in sorted(idaes.config.extra_binaries):
            click.echo(f"    {k}")
        return
    elif show_current_version:
        print_extensions_version()
        return
    if url is None and release is None:
        # the default release is only used if neither a release or url is given
        release = idaes.config.default_binary_release
    if url is not None and release is not None:
        click.echo("\n* You must provide either a release or url not both.")
    elif url is not None or release is not None:
        click.echo("Getting files...")
        try:
            d = idaes.util.download_bin.download_binaries(
                release,
                url,
                insecure,
                cacert,
                verbose,
                platform,
                nochecksum,
                library_only,
                no_download,
                extras_only,
                extra,
                alt_path=to)
            click.echo("Done")
        except idaes.util.download_bin.UnsupportedPlatformError as e:
            click.echo("")
            click.echo(e)
            click.echo("")
            click.echo(
                "Specify one of the following platforms with --platform <os>:")
            for i in sorted(idaes.config.known_binary_platform):
                if i == platform:
                    # auto or linux specified, and it didn't work.
                    # will need to be more specific.
                    continue
                click.echo(f"  {i}")
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
    required=True)
@click.option(
    "--path",
    help="Directory of release files",
    default=None,
    required=True)
def hash_extensions(release, path):
    hfile = f"sha256sum_{release}.txt"
    if path is not None:
        hfile = os.path.join(path, hfile)

    def _write_hash(fp, pack, plat):
        f = f"idaes-{pack}-{plat}-64.tar.gz"
        if path is not None:
            h = idaes.util.download_bin.hash_file_sha256(os.path.join(path, f))
        else:
            h = idaes.util.download_bin.hash_file_sha256(f)
        fp.write(h)
        fp.write("  ")
        fp.write(f)
        fp.write("\n")

    with open(hfile, "w") as f:
        for plat in idaes.config.basic_platforms():
            for pack in ["solvers", "lib"]:
                _write_hash(f, pack, plat)
        for plat in idaes.config.basic_platforms():
            for pack, sp in idaes.config.extra_binaries.items():
                if plat not in sp:
                    continue
                _write_hash(f, pack, plat)
