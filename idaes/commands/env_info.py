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
"""Get info about the environment IDAES is running in and the IDAES version
"""

__author__ = "John Eslick"

import click
from idaes.core.util.env_info import EnvironmentInfo
from idaes.commands import cb


@cb.command(
    name="environment-info",
    help="Print information about idaes, OS, dependencies...")
@click.option(
    "--solver",
    multiple=True,
    help="Add solvers to list of solvers to check"
)
@click.option(
    "--json",
    default=None,
    help="Write output ot a file"
)
def environment_info(solver, json):
    info = EnvironmentInfo(additional_solvers=solver)
    d = info.to_dict()
    if json is None:
        for k, v in d.items():
            click.echo("")
            click.echo(f"{k}")
            for l, q in v.items():
                click.echo(f"    {l}: {q}")
    else:
        info.to_json(json)
