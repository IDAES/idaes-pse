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
import idaes.solvers
from idaes.commands import cb

_log = logging.getLogger("idaes.commands.extensions")


def print_extensions_version():
    click.echo("---------------------------------------------------")
    click.echo("IDAES Extensions Build Versions")
    click.echo("===================================================")
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


@cb.command(name="get-extensions-platforms", help="List binary extension platforms")
def get_extensions_platforms():
    click.echo("\nBuild platforms for IDAES binary Extensions.  Most Linux")
    click.echo("platforms are interchangeable.")
    for key, mes in idaes.config.known_binary_platform.items():
        click.echo("    {}: {}".format(key, mes))


@cb.command(name="get-extensions", help="Get solvers and libraries")
@click.option(
    "--url",
    help="URL to download solvers/libraries from",
    default=idaes.config.default_binary_url)
@click.option(
    "--platform",
    help="Platform to download binaries for (default=auto)",
    default="auto")
@click.option("--verbose", help="Show details", is_flag=True)
def get_extensions(url, verbose, platform):
    if url is not None:
        click.echo("Getting files...")
        idaes.solvers.download_binaries(url, verbose, platform)
        click.echo("Done")
        print_extensions_version()
    else:
        click.echo("\n* You must provide a download URL for IDAES binary files.")



@cb.command(name="ver-extensions", help="Get solver and library IDAES package version")
def ver_extensions():
    print_extensions_version()
