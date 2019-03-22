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
from urllib.parse import urlparse, ParseResult

# third-party
from blessings import Terminal
import click
import humanize

# package
from . import DMF, DMFConfig, resource
from . import errors

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
    NOT_SUPPORTED = 3
    IMPORT_RESOURCE = 4
    DMF = 5
    DMF_OPER = 6


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
    """Improved click.Group that will accept unique prefixes for the
    commands, as well as a set of aliases.

    For example, the following code will create `mycommand` as a group,
    and alias the subcommand "info" to invoke the subcommand "status".
    Any unique prefix of "info" (not conflicting with other subcommands or
    aliases) or "status" will work, e.g. "inf" or "stat"::

        @click.group(cls=AliasedGroup, aliases={"info": "status"})
        def mycommand():
            pass
    """

    def __init__(self, aliases=None, **attrs):
        super().__init__(**attrs)
        self._aliases = aliases

    def get_command(self, ctx, cmd_name):
        command = click.Group.get_command(self, ctx, cmd_name)
        if command is None:
            commands = self.list_commands(ctx)
            # if `cmd_name` is an alias, map to canonical name
            matches = [key for key in self._aliases if key.startswith(cmd_name)]
            if len(matches) == 1:
                # check that this alias does not also match a command directly
                alias_matches = [x for x in commands if x.startswith(cmd_name)]
                if len(alias_matches) > 0:
                    ctx.fail(
                        f"Alias prefix {cmd_name}->{matches[0]} also matches command"
                        f"{'s' if len(alias_matches) > 1 else ''}"
                        f": {' '.join(sorted(alias_matches))}"
                    )
                # rename to canonical (matched) command
                cmd_name = self._aliases[matches[0]]
            elif len(matches) > 1:
                ctx.fail(
                    f"Prefix '{cmd_name}' matches multiple aliases: "
                    f"{' '.join(sorted(matches))}"
                )
            # try to match prefix of command
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


class URLType(click.ParamType):
    """Click type for URLs.
    """

    name = "URL"
    envvar_list_splitter = ","

    def convert(self, value, param, ctx):
        if isinstance(value, ParseResult):
            result = value
        else:
            result = urlparse(value)
        return result


@click.group(cls=AliasedGroup, aliases={"describe": "status", "add": "register"})
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
def base_command(verbose, quiet):
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
        except errors.WorkspaceConfNotFoundError:
            click.echo(f"Workspace not found at path='{path}'")
            if path == '.':  # probably just the default
                click.echo("Use --path option to set workspace path.")
            return Code.WORKSPACE_NOT_FOUND
    return Code.OK


@click.command(help="Status")
@click.option("--color/--no-color", default=True, help="Use color for output")
@click.option(
    "--show", "-s", type=click.Choice(["files", "htmldocs", "all"]), multiple=True
)
def status(color, show):
    _log.debug(f"Get status. Show items: {' '.join(show) if show else '(basic)'}")
    t = _cterm if color else _noterm
    if not DMFConfig.configuration_exists():
        click.echo(f"No configuration found at '{DMFConfig.configuration_path()}'")
        return Code.CONFIGURATION_NOT_FOUND
    try:
        d = DMF()
    except errors.WorkspaceConfNotFoundError as err:
        _log.fatal(f"Cannot get status: {err}")
        click.echo(str(err))
        return Code.WORKSPACE_NOT_FOUND

    def item(key, value=None, before="", color=t.green):
        after_key = "" if key == "-" else ":"
        if value is None:
            return f"{before}{color}{key}{after_key}{t.normal}"
        return f"{before}{color}{key}{after_key}{t.normal} {value}"

    indent_spc = "  "

    print(item("settings", color=t.blue))
    indent = indent_spc
    conf = DMFConfig()  # note: must exist due to earlier check
    for key, value in conf.c.items():
        print(item(key, value, before=indent))

    print(item("workspace", color=t.blue))
    indent = indent_spc
    for key, value in (
        ("location", d.root),
        ("name", d.name),
        ("description", d.description),
        ("created", d.meta[d.CONF_CREATED]),
        ("modified", d.meta[d.CONF_MODIFIED]),
    ):
        print(item(key, value, before=indent))
    # optional workspace items
    for thing in show:
        if thing == "files" or thing == "all":
            indent = indent_spc
            print(item("files", before=indent))
            fpath = pathlib.Path(d.datafiles_path)
            fdirs = [dr for dr in fpath.glob("*") if dr.is_dir()]
            indent = indent_spc * 2
            print(item("count", len(fdirs), before=indent))
            total_size = sum(
                (sum((fp.stat().st_size for fp in fd.glob("*"))) for fd in fdirs)
            )
            print(item("total_size", humanize.naturalsize(total_size), before=indent))
        if thing == "htmldocs" or thing == "all":
            indent = indent_spc
            doc_paths = d.get_doc_paths()
            print(item("html_documentation_paths", before=indent))
            indent = indent_spc * 2
            for p in doc_paths:
                print(item("-", p, before=indent))

    return Code.OK


@click.command(help="Register a new object in the DMF workspace")
@click.argument("url", type=URLType(), metavar="FILE")
@click.option(
    "--copy/--no-copy",
    help="Whether input file is copied into DMF "
    "workspace, or referred to by location",
)
@click.option(
    "--resource-type",
    "-t",
    type=click.Choice(tuple(resource.RESOURCE_TYPES)),
    help="Resource type (default=determined from file)",
)
@click.option(
    "--strict/--no-strict",
    help="If inferring the type fails, "
    "with --strict report an error, and with --no-strict fall back "
    "to importing as a generic file",
)
@click.option(
    "--contained",
    help="Add 'contained in' relation",
    metavar="OBJECT-ID",
    multiple=True,
)
@click.option(
    "--derived", help="Add 'derived from' relation", metavar="OBJECT-ID", multiple=True
)
@click.option(
    "--used", help="Add 'used by' relation", metavar="OBJECT-ID", multiple=True
)
@click.option(
    "--prev",
    help="Add 'version of previous' relation",
    metavar="OBJECT-ID",
    multiple=True,
)
def register(resource_type, url, copy, strict, contained, derived, used, prev):
    _log.debug(f"Register object type='{resource_type}' url/path='{url.path}'")
    # process url
    if url.scheme in ("file", ""):
        path = url.path
    else:
        click.echo("Only bare or 'file' scheme allowed in URL")
        return Code.NOT_SUPPORTED
    # create the resource
    _log.debug("create resource")
    try:
        rsrc = resource.Resource.from_file(path, as_type=resource_type, strict=strict)
    except resource.Resource.InferResourceTypeError as err:
        click.echo(f"Failed to infer resource: {err}")
        return Code.IMPORT_RESOURCE
    except resource.Resource.LoadResourceError as err:
        click.echo(f"Failed to load resource: {err}")
        return Code.IMPORT_RESOURCE
    # process relations
    _log.debug("add relations")
    # XXX
    # add the resource
    _log.debug("add resource begin")
    try:
        d = DMF()
        new_id = d.add(rsrc)
    except errors.WorkspaceError as err:
        click.echo(f"Failed to connect to DMF: {err}")
        return Code.WORKSPACE_NOT_FOUND
    except errors.DMFError as err:
        click.echo(f"Failed to connect to DMF: {err}")
        return Code.DMF
    except errors.DuplicateResourceError as err:
        click.echo(f"Failed to add resource: {err}")
        return Code.DMF_OPER
    _log.debug(f"added resource: {new_id}")
    click.echo(new_id)


# Register base commands
base_command.add_command(init)
base_command.add_command(register)
base_command.add_command(status)

if __name__ == '__main__':
    sys.exit(base_command())
