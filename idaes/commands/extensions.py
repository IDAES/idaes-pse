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
"""Commandline Utilities for Managing the IDAES Data Directory"""

__author__ = "John Eslick"

import click
import idaes.solvers
from idaes.commands import cb

@cb.command(help="Get IDAES Versions of Solvers")
@click.option("--url", help="URL to download solver", default=None)
@click.option("--verbose", help="Show details", is_flag=True)
def get_extensions(url, verbose):
    if url is not None:
        if verbose:
            click.echo("Getting extensions from: {}".format(url))
        idaes.solvers.download_binaries(url)
        if verbose:
            click.echo("Added IDAES executables to {}".format(idaes.bin_directory))
            click.echo("Added IDAES libraries to {}".format(idaes.lib_directory))
        click.echo("Done")
    else:
        click.echo("\n* You must provide a download URL for IDAES binary files.")
