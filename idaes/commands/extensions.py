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
def get_extensions(url):
    if url is not None:
        click.echo("Getting libraries from: {}".format(url))
    else:
        click.echo("Getting libraries from: {}".format(local))
    idaes.solvers.download_binaries(url)
