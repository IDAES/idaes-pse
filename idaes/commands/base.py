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
"""
Base command for 'idaes' commandline script
"""

__author__ = "John Eslick"

import click
import logging

# separate command logging from normal IDAES logging
_log = logging.getLogger("idaes.commands")
_h = logging.StreamHandler()
_h.setFormatter(logging.Formatter("%(levelname)-7s %(name)s: %(message)s"))
_log.addHandler(_h)
_log.propagate = False


def level_from_verbosity(vb):
    level = 0  # for pycharm
    if vb >= 3:
        level = logging.DEBUG
    elif vb == 2:
        level = logging.INFO
    elif vb == 1:
        level = logging.WARN
    elif vb == 0:
        level = logging.ERROR
    elif vb == -1:
        level = logging.FATAL
    elif vb <= -2:
        level = logging.FATAL + 1
    return level


def how_to_report_an_error(embed=False):
    msg = [
        "You can report this error by visiting the Github",
        "software page at:",
        "    https://github.com/idaes/idaes-pse/issues",
        "and clicking on 'New issue', or by sending email",
        "to 'idaes-support@idaes.org'. Please include the",
        "command or actions you took, and the resulting",
        "message, in the report.",
    ]
    if not embed:
        bar = "-" * 50
        msg.insert(0, bar)
        msg.append(bar)
    return "\n".join(msg)


@click.group()
@click.version_option(version=None, package_name="idaes-pse")
@click.option(
    "--verbose",
    "-v",
    count=True,
    help="Increase verbosity. Show warnings if given once, "
    "then info, and then debugging messages.",
)
@click.option(
    "--quiet",
    "-q",
    count=True,
    help="Increase quietness. If given once, "
    "only show critical messages. If "
    "given twice, show no messages.",
)
def command_base(verbose, quiet):
    if quiet > 0 and verbose > 0:
        raise click.BadArgumentUsage("Options for verbosity and quietness conflict")
    if verbose > 0:
        _log.setLevel(level_from_verbosity(verbose))
    else:
        _log.setLevel(level_from_verbosity(-quiet))


@command_base.command(help="Show IDAES copyright information")
def copyright():
    click.echo(
        """
================================================================================
 Institute for the Design of Advanced Energy Systems Process Systems
 Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
 software owners: The Regents of the University of California, through
 Lawrence Berkeley National Laboratory,  National Technology & Engineering
 Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
 University Research Corporation, et al. All rights reserved.
 Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
 license information, respectively. Both files are also available online
 at the URL "https://github.com/IDAES/idaes-pse".
================================================================================
"""
    )


@command_base.command(help="Show how long it takes to import command modules")
def import_time(name="import-time"):
    from idaes.commands import _command_import_total_time

    click.echo(f"Time: {_command_import_total_time}")


if __name__ == "__main__":
    # PYLINT-TODO-FIX fix error bypassed by directive
    # pylint: disable=no-value-for-parameter
    command_base()
