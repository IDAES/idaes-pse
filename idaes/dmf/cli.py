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
"""
Command Line Interface for IDAES DMF.

Uses "Click" to handle command-line parsing and dispatch.
"""
# stdlib
from enum import Enum
import logging
import pathlib
import sys

# third-party
from blessings import Terminal
import click
import humanize

# package
from idaes.dmf import DMF, DMFConfig
from idaes.dmf.errors import WorkspaceConfNotFoundError

__author__ = "Dan Gunter"


_log = logging.getLogger(__name__)
_cterm = Terminal()  # color terminal (as per TERM env)
_noterm = Terminal(force_styling=None)  # no styling, regardless of TERM


class Code(Enum):
    """Return codes from the CLI.
    """

    OK = 0
    WORKSPACE_NOT_FOUND = 1
    CONFIGURATION_NOT_FOUND = 2


def level_from_verbosity(vb):
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


class AliasedGroup(click.Group):
    # Command aliases. Key is alias, value is canonical command.
    aliases = {"info": "status"}

    def get_command(self, ctx, cmd_name):
        command = click.Group.get_command(self, ctx, cmd_name)
        if command is None:
            commands = self.list_commands(ctx)
            # if `cmd_name` is an alias, map to canonical name
            matches = [key for key in self.aliases if key.startswith(cmd_name)]
            if len(matches) == 1:
                # check that this alias does not also match a command directly
                alias_matches = [x for x in commands if x.startswith(cmd_name)]
                if len(alias_matches) > 0:
                    ctx.fail(
                        f"Alias prefix {cmd_name}->{matches[0]} also matches command"
                        f"{'s' if len(alias_matches) > 1 else ''}"
                        f": {' '.join(sorted(alias_matches))}"
                    )
                cmd_name = self.aliases[matches[0]]
            elif len(matches) > 1:
                ctx.fail(
                    f"Prefix '{cmd_name}' matches multiple aliases: {' '.join(sorted(matches))}"
                )
            # try to match prefix of command (& alias, for error checking)
            matches = [x for x in commands if x.startswith(cmd_name)]
            if len(matches) == 1:
                command = click.Group.get_command(self, ctx, matches[0])
            elif len(matches) > 1:
                ctx.fail(
                    f"Command '{cmd_name}' could be a prefix for multiple commands: "
                    f"{' '.join(sorted(matches))}"
                )
        # Return found command
        return command  # matched, or None


@click.group(cls=AliasedGroup)
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
def dmfcommand(verbose, quiet):
    if quiet > 0 and verbose > 0:
        raise click.BadArgumentUsage("Options for verbosity and quietness conflict")
    if verbose > 0:
        _log.setLevel(level_from_verbosity(verbose))
    else:
        _log.setLevel(level_from_verbosity(-quiet))


@click.command(help="Initialize")
@click.option(
    "--path", default=".", type=click.Path(), show_default=True, help="Workspace path"
)
@click.option("--create/--no-create", default=False, help="Create new workspace")
@click.option("--name", type=click.STRING, help="Workspace name")
@click.option("--desc", type=click.STRING, help="Workspace description")
@click.option(
    "--html",
    type=click.Path(exists=True, resolve_path=True),
    default=None,
    help="Path to built HTML documentation",
)
def init(path, create, name, desc, html):
    _log.info(f"Initialize with workspace path={path}")
    if create:
        _log.info("Create new workspace")
        if not name:
            name = click.prompt("New workspace name")
        if not desc:
            desc = click.prompt("New workspace description")
        hpath = None if html is None else [html]
        d = DMF(path=path, create=True, name=name, desc=desc, html_paths=hpath)
        click.echo(f"Configuration in '{d.configuration_file}")
    else:
        _log.info("Use existing workspace")
        try:
            _ = DMF(path=path, create=False)
        except WorkspaceConfNotFoundError:
            click.echo(f"Workspace not found at path='{path}'")
            if path == '.':  # probably just the default
                click.echo("Use --path option to set workspace path.")
            return Code.WORKSPACE_NOT_FOUND
    return Code.OK


@click.command(help="Status")
@click.option("--color/--no-color", default=True, help="Use color for output")
@click.option("--show", "-s", type=click.Choice(["files"]), multiple=True)
def status(color, show):
    t = _cterm if color else _noterm
    if not DMFConfig.configuration_exists():
        click.echo(f"No configuration found at '{DMFConfig.configuration_path()}'")
        return Code.CONFIGURATION_NOT_FOUND
    try:
        d = DMF()
    except WorkspaceConfNotFoundError as err:
        _log.fatal(f"Cannot get status: {err}")
        click.echo(str(err))
        return Code.WORKSPACE_NOT_FOUND
    print(f"# {t.bold}Workspace{t.normal}")
    print(f"{t.green}Location:{t.normal} {d.root}")
    print(f"{t.green}Name:{t.normal} {d.name}")
    print(f"{t.green}Description:{t.normal} {d.description}")
    print(f"{t.green}Created:{t.normal} {d.meta[d.CONF_CREATED]}")
    print(f"{t.green}Modified:{t.normal} {d.meta[d.CONF_MODIFIED]}")
    for thing in show:
        if thing == "files":
            print(f"{t.green}Files:{t.normal}")
            fpath = pathlib.Path(d.datafiles_path)
            fdirs = [dr for dr in fpath.glob("*") if dr.is_dir()]
            print(f"  {t.green}Number of files:{t.normal} {len(fdirs)}")
            total_size = sum(
                (sum((fp.stat().st_size for fp in fd.glob("*"))) for fd in fdirs)
            )
            print(
                f"  {t.green}Total size:{t.normal} {humanize.naturalsize(total_size)}"
            )
    return Code.OK


@click.command(help="Add object(s)")
@click.option("--parent", "-p", type=click.STRING, help="Parent object ID")
@click.argument("objects", nargs=-1)
def add(parent, objects):
    if not objects:
        raise click.BadArgumentUsage("Need at least one object to add")
    _log.info(f"Add with parent={parent}, objects={objects}")


dmfcommand.add_command(init)
dmfcommand.add_command(status)

if __name__ == '__main__':
    sys.exit(dmfcommand())

