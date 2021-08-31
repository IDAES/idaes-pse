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
"""Wrapper for running solvers built for Linux on Windows by using WSL2. See
https://docs.microsoft.com/en-us/windows/wsl/install-win10 for information about
setting up the WSL. Combining this utility with a batch file will allow Linux
solver binaries to be used on Windows.
"""

__author__ = "John Eslick"

import os
import re
import idaes
import click
import subprocess
from idaes.commands import cb

@cb.command(
    name="solver-wsl",
    context_settings={"ignore_unknown_options": True},
    help="Run a Linux solver on Windows using WSL")
@click.option("--distribution", required=True, help="Linux distribution")
@click.option("--user", required=True, help="User")
@click.option("--executable", required=True, help="Executable path on WSL")
@click.argument('args', nargs=-1)
def solver_wsl(distribution, user, executable, args):
    al = [None]*len(args)
    for i, a in enumerate(args):
        r = re.match(r"^([A-Za-z]):\\(.*$)", a)
        if r is not None:
            # is a Windows path (well almost certaintly)
            p = r.group(2).replace('\\', '/')
            l = r.group(1).lower()
            a = f"/mnt/{l}/{p}"
        al[i] = a
    r = subprocess.run(["wsl", "-d", distribution, "-u", user, executable] + al)
    return r.returncode
