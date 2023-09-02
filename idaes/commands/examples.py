#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Install IDAES example files locally.
"""
import sys

# third-party
import click

# package
from idaes.commands import cb

MESSAGE = """
The get-examples command has been removed.

The IDAES examples are now available as the `idaes-examples` package in PyPI.
To install, run the following command: 

     pip install idaes-examples 
     
This also installs the 'idaesx' command that can be used to browse the Jupyter
notebooks. The simplest way to do this is with the embedded desktop UI:

    idaesx gui
    
For more details, see the IDAES documentation at https://idaes-pse.readthedocs.io/
and the idaes-examples page on PyPI at https://pypi.org/project/idaes-examples/
"""


@cb.command(
    name="get-examples", help="(legacy) install IDAES example Jupyter notebooks"
)
def get_examples():
    """Legacy get-examples command."""
    click.echo(MESSAGE)
    sys.exit(0)
