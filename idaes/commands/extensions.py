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

import zipfile
import os
import idaes
import click
from shutil import copyfile
from idaes.commands import cb

_valid_solvers = ["ipopt"]

@cb.command(help="Get IDAES Compiled Libraries (on Windows includes MinGW runtime)")
@click.option("--url", help="URL to download libraries", default=None)
@click.option("--local", help="Use local library zip file", default=None)
def get_libraries(url, local):
    idaes._create_bin_dir()
    zip_file = os.path.join(idaes.lib_directory, "idaes-bin.zip")
    # Here we'll setup the default urls eventally but for now local file
    if local is None:
        # this is temporary I know it has holes (mac...)
        if os.name == 'nt':
            local = "idaes-lib-win.zip"
        else:
            local = "idaes-lib-linux.zip"

    if url is not None:
        print("Getting libraries from: {}".format(url))
        r = requests.get(url)
    else:
        print("Getting libraries from: {}".format(local))
        copyfile(local, zip_file)
    # unzip
    with zipfile.ZipFile(zip_file, 'r') as f:
        f.extractall(idaes.lib_directory)

@cb.command(help="Get IDAES Versions of Solvers")
@click.argument('solver')
@click.option("--url", help="URL to download solver", default=None)
@click.option("--local", help="Use local zip file", default=None)
@click.option("--remove", help="Remove a solver", is_flag=True)
def get_solver(solver, url, local, remove):
    idaes._create_bin_dir()
    # maybe do loop here so you can supply multiple solvers, this colud use
    # improvment espcially if we have more than one solver, but we should
    # be able to fix it up and maintain backward compatablity in the docs.
    if solver not in _valid_solvers:
        raise Exception("The solver '{}' is not available.")
    elif remove:
        if os.name == "nt":
            rmpath = os.path.join(idaes.bin_directory, solver+".exe")
        else:
            rmpath = os.path.join(idaes.bin_directory, solver)
        if click.confirm("Delete {}?".format(rmpath)):
            os.remove(rmpath)
        else:
            print("Not removing {}.".format(solver))
    else:
        zip_file = os.path.join(idaes.bin_directory, "idaes-{}.zip".foramt(solver))
        # Here we'll setup the default urls eventally but for now local file
        if local is None:
            # this is temporary I know it has holes (mac...)
            if os.name == 'nt':
                local = "idaes-{}-win.zip".foramt(solver)
            else:
                local = "idaes-{}-linux.zip".foramt(solver)

        if url is not None:
            print("Getting libraries from: {}".format(url))
            r = requests.get(url)
        else:
            print("Getting libraries from: {}".format(local))
            copyfile(local, zip_file)
        # unzip
        with zipfile.ZipFile(zip_file, 'r') as f:
            f.extractall(idaes.bin_directory)
        if os.name == "nt":
            print("""\n   ***Be sure to install IDAES liraries (idaes get-libraries)
          to ensure you have required MinGW runtime libraires.\n""")
