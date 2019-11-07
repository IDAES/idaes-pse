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

import os
import idaes
import click
from idaes.commands import cb

@cb.command(help="Show the IDAES data directory path")
# the underscore get turned into a '-' so the command is data-directory
@click.option("--exists", is_flag=True)
def data_directory(exists):
    if exists:
        print(os.path.exists(idaes.data_directory))
    else:
        print(idaes.data_directory)
