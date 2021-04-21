##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
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
    show_platforms):

    if show_platforms:
        click.echo("\nBuild platforms for IDAES binary Extensions.  Most Linux")
        click.echo("platforms are interchangeable.")
        for key, mes in idaes.config.known_binary_platform.items():
            click.echo("    {}: {}".format(key, mes))
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
                no_download)
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
