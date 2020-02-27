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
import logging
import idaes.solvers
from idaes.core.commands import cb

_log = logging.getLogger("idaes.commands.extensions")


@cb.command(name="get-extensions", help="Get solvers and libraries")
@click.option(
    "--url",
    help="URL to download solver",
    default=idaes._config.default_binary_url)
def get_extensions(url):
    verbose = _log.isEnabledFor(logging.INFO)
    if url is not None:
        click.echo("Getting files...")
        idaes.solvers.download_binaries(url, verbose)
        click.echo("Done")
    else:
        click.echo("\n* You must provide a download URL for IDAES binary files.")
