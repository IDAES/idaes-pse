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
"""Utility for copying examples to a convienient location"""

__author__ = "John Eslick"

import os
import shutil
import click
from idaes.commands import cb
import pkg_resources


@cb.command(help="Copy example files")
@click.option(
    "--location", default="./idaes-examples", help="Path to copy example files to"
)
# the underscore get turned into a '-' so the command is get-examples
def get_examples(location):
    examples_path = pkg_resources.resource_filename("idaes.examples", "")
    todir = os.path.abspath(location)
    print("Copying '{}' to '{}'...".format(examples_path, todir))
    try:
        shutil.copytree(examples_path, todir)
    except FileExistsError:
        print("New examples directory '{}' already exists.".format(todir))
