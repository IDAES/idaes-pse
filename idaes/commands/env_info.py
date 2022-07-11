#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""Get info about the environment IDAES is running in and the IDAES version
"""

__author__ = "John Eslick"

import click
from idaes.commands import cb


@cb.command(
    name="environment-info", help="Print information about idaes, OS, dependencies..."
)
@click.option("--solver", multiple=True, help="Add solvers to list of solvers to check")
@click.option("--json", default=None, help="Write output ot a file")
def environment_info(solver, json):
    from idaes.core.util.env_info import EnvironmentInfo

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
