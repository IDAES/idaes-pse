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
import idaes
import click
from idaes.commands import cb

@cb.command(name="data-directory", help="Show IDAES data directory")
@click.option("--exists", is_flag=True, help="Show if the directory exists")
@click.option("--create", is_flag=True, help="Create the directory")
def data_directory(exists, create):
    if create:
        click.echo("Creating {}".format(idaes.data_directory))
        idaes._create_data_dir()
    elif exists:
        click.echo(os.path.exists(idaes.data_directory))
    else:
        click.echo(idaes.data_directory)

@cb.command(name="bin-directory", help="Show IDAES executable file directory")
@click.option("--exists", is_flag=True, help="Show if the directory exists")
@click.option("--create", is_flag=True, help="Create the directory")
def bin_directory(exists, create):
    if create:
        click.echo("Creating {}".format(idaes.bin_directory))
        idaes._create_bin_dir()
    elif exists:
        click.echo(os.path.exists(idaes.bin_directory))
    else:
        click.echo(idaes.bin_directory)
