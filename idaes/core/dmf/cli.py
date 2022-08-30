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
Command Line Interface for idaes.core.dmf.

Uses "Click" to handle command-line parsing and dispatch.
"""
# stdlib
import time
from collections import namedtuple
from datetime import datetime
from enum import Enum
import json
import logging
from operator import itemgetter
import os
import pathlib
import sys
from typing import List
from urllib.parse import urlparse, ParseResult
import yaml

# third-party
import click

# package
from idaes.core.dmf import DMF, DMFConfig, resource, workspace, create_configuration
from idaes.core.dmf.resource import Predicates
from idaes.core.dmf import errors
from idaes.core.dmf.workspace import Fields
from idaes.core.dmf import util
from idaes.core.dmf.util import (
    ColorTerm,
    yaml_load,
    parse_datetime,
    size_prefix,
)

__author__ = "Dan Gunter"

_log = logging.getLogger(__name__)
_dmf_log = logging.getLogger("idaes.core.dmf")
_timing_log = logging.getLogger("idaes.core.dmf.tm")


class Code(Enum):
    """Return codes from the CLI."""

    OK = 0
    WORKSPACE_NOT_FOUND = 1
    CONFIGURATION_NOT_FOUND = 2
    NOT_SUPPORTED = 3
    IMPORT_RESOURCE = 4
    DMF = 5
    DMF_OPER = 6
    INPUT_VALUE = 7
    CANCELED = 8


# Some common fields we can show
SHOW_FIELDS = "type, desc, creator, created, modified, fikes, codes, version"


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
        t0 = time.time()
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
    """Click type for URLs."""

    name = "URL"
    envvar_list_splitter = ","

    def convert(self, value, param, ctx):
        if isinstance(value, ParseResult):
            result = value
        else:
            result = urlparse(value)
        return result


@click.group(
    cls=AliasedGroup,
    aliases={
        "describe": "status",
        "add": "register",
        "resource": "info",
        "show": "info",
        "list": "ls",
        "delete": "rm",
        "graph": "related",
        "workspace": "init",
    },
    help="Data management framework command wrapper",
)
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
    """Data management framework command wrapper.

    This command does nothing by itself except provide global
    options and list subcommands.
    """
    if quiet > 0 and verbose > 0:
        raise click.BadArgumentUsage("Options for verbosity and quietness conflict")
    if verbose > 0:
        _dmf_log.setLevel(level_from_verbosity(verbose))
    else:
        _dmf_log.setLevel(level_from_verbosity(-quiet))


@click.command(help="Set up the DMF configuration file.")
@click.option(
    "--dir",
    "-d",
    "config_path",
    type=click.Path(),
    default=None,
    help="Specify non-default directory for configuration",
)
def setup(config_path):
    try:
        path = create_configuration(config_path=config_path)
    except ValueError as err:
        click.echo(f"Error creating DMF configuration: {err}")
        sys.exit(Code.DMF_OPER.value)
    click.echo(f"Created configuration in '{path}'")


@click.command(
    help="Initialize the current workspace. Optionally, create a new workspace."
)
@click.argument("path", type=click.Path())
@click.option(
    "--create/--no-create",
    default=False,
    help="Create new workspace. If `--name` and `--desc` are not provided, these will"
    " be prompted for interactively.",
)
@click.option("--name", type=click.STRING, help="Workspace name")
@click.option("--desc", type=click.STRING, help="Workspace description")
@click.option(
    "--html",
    type=click.Path(exists=True, resolve_path=True),
    default=None,
    help="Path to built HTML documentation",
)
def init(path, create, name, desc, html):
    """Initialize the current workspace used for the data management framework commands.
    Optionally, create a new workspace.
    """
    _log.info(f"Initialize with workspace path={path} cwd={os.path.abspath(os.curdir)}")
    if create:
        _log.info("Create new workspace")
        # pre-check that there is no file/dir by this name
        try:
            wspath = pathlib.Path(path)
            if (
                wspath.exists()
                and (wspath / workspace.Workspace.WORKSPACE_CONFIG).exists()
            ):
                click.echo(
                    f"Cannot create workspace: '{path}/{workspace.Workspace.WORKSPACE_CONFIG}' already exists"
                )
                sys.exit(Code.DMF_OPER.value)
        except PermissionError:
            click.echo(f"Cannot create workspace: path '{path}' not accessible")
            sys.exit(Code.DMF_OPER.value)
        if not name:
            name = click.prompt("New workspace name")
        if not desc:
            desc = click.prompt("New workspace description")
        # Note: default HTML paths, and all other default values, are included
        # in the JSON schema at `idaes.core.dmf.workspace.CONFIG_SCHEMA`
        hpath = [html] if html else None
        try:
            d = DMF(
                path=path,
                create=True,
                name=name,
                desc=desc,
                html_paths=hpath,
                add_defaults=True,
            )
        except errors.WorkspaceError as err:
            click.echo(f"Cannot create workspace: {err}")
            sys.exit(Code.DMF_OPER.value)
        click.echo(f"Configuration in '{d.configuration_file}")
    else:
        _log.info("Use existing workspace")
    # In either case, switch to provided config
    try:
        _ = DMF(path=path, create=False, save_path=True)  # noqa: F841
    except errors.WorkspaceConfNotFoundError:
        click.echo(f"Workspace configuration not found at path='{path}'")
        sys.exit(Code.WORKSPACE_NOT_FOUND.value)
    except errors.WorkspaceNotFoundError:
        click.echo(f"Existing workspace not found at path='{path}'")
        click.echo("Add --create flag to create a workspace.")
        sys.exit(Code.WORKSPACE_NOT_FOUND.value)


@click.command(help="Get status of workspace")
@click.option("--color/--no-color", default=True, help="Use color for output")
@click.option(
    "--show",
    "-s",
    type=click.Choice(["files", "htmldocs", "logging", "all"]),
    multiple=True,
)
@click.option(
    "--all", "-a", "show_all", flag_value="yes", help="Synonym for `--show all`"
)
def status(color, show, show_all):
    if show_all == "yes":
        show = ["all"]
    _log.debug(f"Get status. Show items: {' '.join(show) if show else '(basic)'}")
    t = ColorTerm(enabled=color)
    if not DMFConfig.configuration_exists():
        click.echo(f"No configuration found at '{DMFConfig.configuration_path()}'")
        sys.exit(Code.CONFIGURATION_NOT_FOUND.value)
    try:
        d = DMF()
    except errors.WorkspaceConfNotFoundError as err:
        _log.fatal(f"Cannot get status: {err}")
        click.echo(str(err))
        sys.exit(Code.WORKSPACE_NOT_FOUND.value)

    # pretty-display a key/value pair or list value
    def item(key, value=None, before="", color=t.green):
        after_key = "" if key == "" else ":"
        if value is None:
            return f"{before}{color}{key}{after_key}{t.reset}"
        elif key is None:
            return f"{before}{color}{value}{t.reset}"
        return f"{before}{color}{key}{after_key}{t.reset} {value}"

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

    _show_optional_workspace_items(d, show, indent_spc, item, t=t)


def _show_optional_workspace_items(d, items, indent_spc, item_fn, t=None):
    indent = indent_spc
    for thing in items:
        if thing == "files" or thing == "all":
            print(item_fn("files", before=indent))
            fpath = pathlib.Path(d.datafiles_path)
            fdirs = [dr for dr in fpath.glob("*") if dr.is_dir()]
            indent = indent_spc * 2
            print(item_fn("count", len(fdirs), before=indent))
            total_size = sum(
                (sum((fp.stat().st_size for fp in fd.glob("*"))) for fd in fdirs)
            )
            print(item_fn("total_size", size_prefix(total_size), before=indent))
        if thing == "htmldocs" or thing == "all":
            indent = indent_spc
            doc_paths = d.get_doc_paths()
            print(item_fn("html_documentation_paths", before=indent))
            indent = indent_spc * 2
            for p in doc_paths:
                print(item_fn("-", p, before=indent))
        if thing == "logging" or thing == "all":
            indent = indent_spc
            print(item_fn("logging", before=indent))
            log_conf = d.meta.get(Fields.LOG_CONF, None)
            indent = 2 * indent_spc
            if log_conf is None:
                print(item_fn(None, "not configured", before=indent, color=t.yellow))
            else:
                for subconf in sorted(log_conf.keys()):
                    print(item_fn(subconf, before=indent))
                    indent2 = indent + indent_spc
                    for i2, v2 in log_conf[subconf].items():
                        print(item_fn(i2, v2, before=indent2))


@click.command(help="Register a new object in the DMF workspace")
@click.argument("url", type=URLType(), metavar="FILE")
@click.option("--info", help="Show info on created resource", flag_value="yes")
@click.option(
    "--copy/--no-copy",
    default=True,
    help="Whether input file is copied into DMF "
    "workspace, or referred to by location",
)
@click.option(
    "--type",
    "-t",
    "resource_type",
    type=click.Choice(tuple(resource.ResourceTypes.all())),
    help="Resource type (default=determined from file)",
)
@click.option(
    "--strict/--no-strict",
    help="If inferring the type fails, "
    "with `--strict` report an error, and with `--no-strict` fall back "
    "to importing as a generic file",
)
@click.option(
    "--unique/--no-unique",
    default=True,
    help="Stop if another resource has a file matching this file's "
    "name and contents",
)
@click.option(
    "--contained",
    help="Add a 'contained in' relation to the given resource",
    metavar="RESOURCE-ID",
    multiple=True,
)
@click.option(
    "--derived", help="Add 'derived from' relation", metavar="OBJECT-ID", multiple=True
)
@click.option(
    "--used",
    help="Add 'used by' relation to the given resource",
    metavar="RESOURCE-ID",
    multiple=True,
)
@click.option(
    "--prev",
    help="Add 'version of previous' relation to the given resource",
    metavar="RESOURCE-ID",
    multiple=True,
)
@click.option(
    "--is-subject",
    flag_value="yes",
    help="If given, this resource will be the subject rather than the object "
    "of any and all relations added.",
)
@click.option(
    "--version",
    help="Set semantic version for this resource (default=0.0.0)",
    default=None,
)
def register(
    resource_type,
    url,
    info,
    copy,
    strict,
    unique,
    contained,
    derived,
    used,
    prev,
    is_subject,
    version,
):
    _log.debug(f"Register object type='{resource_type}' url/path='{url.path}'")
    # process url
    if url.scheme in ("file", ""):
        path = url.path
    else:
        click.echo("Currently, URL must be a file")
        sys.exit(Code.NOT_SUPPORTED.value)
    # create the resource
    _log.debug("create resource")
    try:
        rsrc = resource.Resource.from_file(
            path, as_type=resource_type, strict=strict, do_copy=copy
        )
    except resource.Resource.InferResourceTypeError as err:
        click.echo(f"Failed to infer resource: {err}")
        sys.exit(Code.IMPORT_RESOURCE.value)
    except resource.Resource.LoadResourceError as err:
        click.echo(f"Failed to load resource: {err}")
        sys.exit(Code.IMPORT_RESOURCE.value)
    # connect to DMF
    try:
        dmf = DMF()
    except errors.WorkspaceError as err:
        click.echo(f"Failed to connect to DMF: {err}")
        sys.exit(Code.WORKSPACE_NOT_FOUND.value)
    except errors.DMFError as err:
        click.echo(f"Failed to connect to DMF: {err}")
        sys.exit(Code.DMF.value)
    # check uniqueness
    if unique:
        df = rsrc.v["datafiles"][0]  # file info for this upload
        query = {"datafiles": [{"sha1": df["sha1"]}]}
        query_result, dup_ids = dmf.find(query), []
        for dup in query_result:
            dup_df = dup.v["datafiles"][0]
            if dup_df["path"] in df["path"]:
                dup_ids.append(dup.id)
        n_dup = len(dup_ids)
        if n_dup > 0:
            click.echo(
                f"This file is already in {n_dup} resource(s): " f"{' '.join(dup_ids)}"
            )
            sys.exit(Code.DMF_OPER.value)
    # process relations
    _log.debug("add relations")
    rel_to_add = {  # translate into standard relation names
        Predicates.contains: contained,
        Predicates.derived: derived,
        Predicates.uses: used,
        Predicates.version: prev,
    }
    target_resources = {}  # keep target resources in dict, update at end
    for rel_name, rel_ids in rel_to_add.items():
        for rel_id in rel_ids:
            if rel_id in target_resources:
                rel_subj = target_resources[rel_id]
            else:
                rel_subj = dmf.fetch_one(rel_id)
                target_resources[rel_id] = rel_subj
            if rel_subj is None:
                click.echo(f"Relation {rel_name} target not found: {rel_id}")
                sys.exit(Code.DMF_OPER.value)
            if is_subject == "yes":
                resource.create_relation(rsrc, rel_name, rel_subj)
            else:
                resource.create_relation(rel_subj, rel_name, rsrc)
            _log.debug(f"added relation {rsrc.id} <-- {rel_name} -- {rel_id}")
    _log.debug("update resource relations")
    for rel_rsrc in target_resources.values():
        dmf.update(rel_rsrc)
    # add metadata
    if version:
        try:
            vlist = resource.version_list(version)
        except ValueError:
            click.echo(f"Invalid version `{version}`")
            sys.exit(Code.INPUT_VALUE.value)
        else:
            rsrc.v["version_info"]["version"] = vlist
    # add the resource
    _log.debug("add resource begin")
    try:
        new_id = dmf.add(rsrc)
    except errors.DuplicateResourceError as err:
        click.echo(f"Failed to add resource: {err}")
        sys.exit(Code.DMF_OPER.value)
    _log.debug(f"added resource: {new_id}")
    if info == "yes":
        pfxlen = len(new_id)
        si = _ShowInfo("term", pfxlen)
        for rsrc in dmf.find_by_id(new_id):
            si.show(rsrc)
    else:
        click.echo(new_id)


@click.command(help="List resources in the workspace")
@click.option("--color/--no-color", default=True, help="Use color for output")
@click.option(
    "--show",
    "-s",
    multiple=True,
    help=f"Show given field(s) in listing. Common fields: {SHOW_FIELDS}",
)
@click.option(
    "--sort",
    "-S",
    "sort_by",
    type=click.Choice(["id", "type", "desc", "created", "modified", "version", "name"]),
    multiple=True,
    help="Sort by given field; if repeated, combine to make a compound sort key",
)
@click.option(
    "--prefix/--no-prefix",
    "prefix",
    help="By default, shown 'id' is the shortest unique prefix; "
    "`--no-prefix` shows full id",
    default=True,
)
@click.option("--reverse", "-r", "reverse", flag_value="yes", help="Reverse sort order")
def ls(color, show, sort_by, reverse, prefix):
    try:
        d = DMF()
    except errors.WorkspaceError as e:
        print(f"Workspace error: {e}")
        sys.exit(Code.WORKSPACE_NOT_FOUND.value)
    if not show:
        show = ["type", "desc", "modified"]  # note: 'id' is always first
    else:
        try:
            show = _split_and_validate_fields(show)
        except ValueError as err:
            click.echo(f"Bad fields for --show option: {err}")
            sys.exit(Code.INPUT_VALUE.value)
    reverse = bool(reverse == "yes")
    if not sort_by:
        sort_by = ["id"]
    resources = list(d.find())
    _print_resource_table(resources, show, sort_by, reverse, prefix, color)


def _split_and_validate_fields(fields: List[str]) -> List[str]:
    """Split comma-separated 'show' fields, validate that they are
    part of a resource, and return a list of tuples.

    Raises:
        ValueError: if the fields do not validate
    Returns:

    """
    dummy = resource.DUMMY_RESOURCE
    # Find every field in the resource
    result = []
    for grp in fields:
        for f in grp.split(","):
            f = f.strip()
            if f in _show_fields:
                result.append(f)
                continue
            if "." in f:
                keys, v = f.split("."), dummy
                k = None
                try:
                    for k in keys:
                        try:
                            v = v[k]
                        except TypeError:
                            raise ValueError(f"bad field {f}")
                except KeyError:
                    raise ValueError(f"bad field {f} (error at " f"sub-field {k})")
                result.append(tuple(keys))
            elif f not in dummy:
                raise ValueError(f"bad field {f}")
            else:
                result.append(f)
    return result


def _print_resource_table(resources, show_fields, sort_by, reverse, prefix, color):
    """Text-mode `ls`."""
    t = ColorTerm(enabled=color)
    if len(resources) == 0:
        print("no resources to display")
        return
    uuid_len = util.uuid_prefix_len([r.id for r in resources])
    full_len = max((len(r.id) for r in resources)) if not prefix else uuid_len
    fields = ["id"] + list(show_fields)
    nfields = len(fields)
    # calculate table body. do this first to get widths.
    hdr_fields = [".".join(f) if isinstance(f, tuple) else f for f in fields]
    rows, maxwid, widths = [], 60, [len(f) for f in hdr_fields]
    for r in resources:
        row = []
        for i, fld in enumerate(fields):
            # transform field for display
            try:
                transformer = _show_fields[fld]
            except KeyError:
                transformer = _IdentityField(fld)
            transformer.set_value(r)
            # if it's a UUID field, add info about unique prefix length
            is_id_field = isinstance(transformer, _IdField)
            if is_id_field:
                transformer.term = t
                transformer.pfxlen = uuid_len
                transformer.flen = full_len
            # extract display string and set field width
            s = str(transformer)
            slen = full_len if is_id_field else len(s)
            if slen > widths[i]:
                if slen > maxwid:
                    s, widths[i] = s[:maxwid], maxwid
                else:
                    widths[i] = slen
            row.append(s)
        # append sort keys (not displayed)
        sort_obj = [_show_fields[fld] for fld in sort_by]
        for o in sort_obj:
            o.set_value(r)
        row.extend([o.value for o in sort_obj])
        # add row to table body
        rows.append(row)
    # sort rows
    nsort = len(sort_by)
    if nsort > 0:
        sort_key_idx = tuple(range(nfields, nfields + nsort))
        rows.sort(key=itemgetter(*sort_key_idx), reverse=reverse)
    # print table header
    hdr_columns = [
        t.bold + "{{f:{w}s}}".format(w=w).format(f=f)
        for f, w in zip(hdr_fields, widths)
    ]
    print(" ".join(hdr_columns) + t.normal)
    # print table body
    for row in rows:
        col = []
        for i, fld in enumerate(fields):
            if i == 0:
                col.append(t.red)
            elif fld == resource.Resource.TYPE_FIELD:
                col.append(t.yellow)
            elif fld != "desc":
                col.append(t.green)
            else:
                col.append("")
        row_columns = [
            f"{col[i]}{f:{w}}{t.reset}"
            for i, f, w in zip(range(len(widths)), row[:nfields], widths)
        ]
        print(" ".join(row_columns))


@click.command(help="Find resources in the workspace")
@click.option("--color/--no-color", default=True, help="Use color for output")
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["list", "info", "json"]),
    default="list",
    help="Output format",
)
@click.option(
    "--show",
    "-s",
    multiple=True,
    help=f"Show given field(s) in listing. Common fields: {SHOW_FIELDS}",
)
@click.option(
    "--sort",
    "-S",
    "sort_by",
    type=click.Choice(["id", "type", "desc", "created", "modified", "version"]),
    multiple=True,
    help="Sort by given field; if repeated, combine to make a compound sort key",
)
@click.option("--reverse", "-r", "reverse", flag_value="yes", help="Reverse sort order")
@click.option(
    "--prefix/--no-prefix",
    "prefix",
    help="By default, shown 'id' is the shortest unique prefix; "
    "`--no-prefix` shows full id",
    default=True,
)
@click.option("--by", default="", help="Creator name")
@click.option("--created", default="", help="Creation date or date range")
@click.option("--file", "filedesc", default="", help="File desc(ription)")
@click.option("--modified", default="", help="Modification date or date range")
@click.option("--name", default="", help="Matches any of the aliases")
@click.option("--type", "datatype", default="", help="The resource type")
def find(
    output_format,
    color,
    show,
    sort_by,
    prefix,
    reverse,
    by,
    created,
    filedesc,
    modified,
    name,
    datatype,
):
    d = DMF()
    if output_format == "list":
        if not show:
            show = ["type", "desc", "modified"]  # note: 'id' is always first
        else:
            try:
                show = _split_and_validate_fields(show)
            except ValueError as err:
                click.echo(f"Bad fields for --show option: {err}")
                sys.exit(Code.INPUT_VALUE.value)
    reverse = bool(reverse == "yes")
    if not sort_by:
        sort_by = ["id"]

    # Build query
    query = {}
    if by:
        query["creator.name"] = by
    if created:
        try:
            query["created"] = _date_query(created)
        except ValueError as err:
            click.echo(f"bad date for 'created': {err}")
            sys.exit(Code.INPUT_VALUE.value)
    if filedesc:
        query["datafiles"] = [{"desc": filedesc}]
    if modified:
        try:
            query["modified"] = _date_query(modified)
        except ValueError as err:
            click.echo(f"bad date for 'modified': {err}")
            sys.exit(Code.INPUT_VALUE.value)
    if name:
        query["aliases"] = [name]
    if datatype:
        query["type"] = datatype

    # Execute query
    _log.info(f"find: query = '{query}'")
    _log.debug("find.begin")
    resources = list(d.find(query))
    _log.debug("find.end")

    # Print result
    if output_format == "list":
        # print resources like `ls`
        _print_resource_table(resources, show, sort_by, reverse, prefix, color)
    elif output_format == "info":
        # print resources one by one
        si = _ShowInfo("term", 32, color=color)
        for rsrc in resources:
            si.show(rsrc)
    elif output_format == "json":
        # print resources as JSON
        for r in resources:
            print(json.dumps(r.v, indent=2))


def _date_query(datestr):
    """Transform date(s) into a query.

    Raises:
        ValueError: If date parsing fails
    """

    def _ts(x):  # get timestamp of date or time
        try:
            return x.timestamp()
        except AttributeError:
            return datetime(*x.timetuple()[:6]).timestamp()

    # Parse "<Date1>..<Date2>" or just "<Date>"
    dates = datestr.split("..", 1)
    _log.debug(f"date '{datestr}', split: {dates}")
    parsed_dates = [parse_datetime(d) if d else None for d in dates]
    # Range <Date1>..<Date2>
    if len(parsed_dates) == 2:
        begin_date, end_date = parsed_dates
        query = {}
        if begin_date:
            query["$ge"] = _ts(begin_date)
        if end_date:
            query["$le"] = _ts(end_date)
    # Single date(time)
    else:
        pd = parsed_dates[0]
        if pd is None:
            raise ValueError(f"Empty date: {dates[0]}")
        if isinstance(pd, datetime):
            _log.warning("datetime given, must match to the second")
            query = _ts(pd)
        else:
            query = {"$ge": _ts(pd), "$le": _ts(pd) + 60 * 60 * 24}  # 1 day
    return query


@click.command(help="Show detailed information about a resource")  # aliases: resource
@click.argument("identifier")
@click.option("--multiple/--no-multiple", default=False)
@click.option(
    "--format",
    "-f",
    "output_format",
    type=click.Choice(["term", "json", "jsonc"]),
    multiple=False,
    default="term",
    help="Output format: `term` (the default) shows terminal-friendly output, with "
    "some colors; `json` outputs JSON text; `jsonc` outputs compact JSON (no indents"
    "or line breaks)",
)
@click.option(
    "--color/--no-color",
    default=True,
    help="In term mode, use (or do not use)" " color for output",
)
def info(identifier, multiple, output_format, color):
    _log.debug(f"info for resource id='{identifier}'")
    try:
        resource.identifier_str(identifier, allow_prefix=True)
    except ValueError as err:
        click.echo(f"{err}")
        sys.exit(Code.INPUT_VALUE.value)
    rsrc_list = list(find_by_id(identifier))
    n = len(rsrc_list)
    if n > 1 and not multiple:
        click.echo(
            f"Too many ({n}) resources match prefix '{identifier}'. "
            "Add option --multiple to allow multiple matches."
        )
        sys.exit(Code.DMF_OPER.value)
    elif n == 0:
        click.echo("Resource not found")
        sys.exit(Code.DMF_OPER.value)
    pfxlen = len(identifier)
    si = _ShowInfo(output_format, pfxlen, color=color)
    for rsrc in rsrc_list:
        si.show(rsrc)


@click.command(help="Show resources related (connected) to a given resource")
@click.argument("identifier")
@click.option(
    "-d",
    "--direction",
    type=click.Choice(["out", "in"]),
    default="out",
    help="Direction of relationship to show (default=out[going])",
)
@click.option(
    "--color/--no-color",
    default=True,
    help="In term mode, use (or do not use)" " color for output",
)
@click.option(
    "--unicode/--no-unicode",
    default=True,
    help="With `--unicode`, allow unicode characters for connecting "
    "output items with 'lines'. With --no-unicode, use plain ASCII "
    "characters for this",
)
def related(identifier, direction, color, unicode):
    _log.info(f"related to resource id='{identifier}'")
    t = ColorTerm(enabled=color)
    dmf = DMF()
    try:
        resource.identifier_str(identifier, allow_prefix=True)
    except ValueError as err:
        click.echo(f"{err}")
        sys.exit(Code.INPUT_VALUE.value)
    _log.debug(f"begin: finding root resource {identifier}")
    rsrc_list = list(find_by_id(identifier, dmf=dmf))
    n = len(rsrc_list)
    if n > 1:
        click.echo(f"Too many resources matching `{identifier}`")
        sys.exit(Code.INPUT_VALUE)
    rsrc = rsrc_list[0]
    _log.debug(f"end: finding root resource {identifier}")
    # get related resources
    _log.debug(f"begin: finding related resources for {identifier}")
    outgoing = direction == "out"
    rr = list(dmf.find_related(rsrc, meta=["aliases", "type"], outgoing=outgoing))
    _log.debug(f"end: finding related resources for {identifier}")
    # stop if no relations
    if not rr:
        _log.warning(f"no resource related to {identifier}")
        click.echo(f"No relations for resource `{identifier}`")
        sys.exit(0)
    _log.info(f"got {len(rr)} related resources")
    # debugging
    if _log.isEnabledFor(logging.DEBUG):
        dbgtree = "\n".join(["  " + str(x) for x in rr])
        _log.debug(f"related resources:\n{dbgtree}")
    # extract uuids & determine common UUID prefix length
    uuids = [item[2][resource.Resource.ID_FIELD] for item in rr]
    pfx = util.uuid_prefix_len(uuids)
    # initialize queue with depth=1 items
    q = [item for item in rr if item[0] == 1]
    # print root resource
    print(_related_item(rsrc.id, rsrc.name, rsrc.type, pfx, t, unicode))
    # set up printing style
    if unicode:
        # connector chars
        vrt, vrd, relbow, relbow2, rarr = (
            "\u2502",
            "\u2506",
            "\u2514",
            "\u251C",
            "\u2500\u2500",
        )
        # relation prefix and arrow
        relpre, relarr = (
            ["\u2500\u25C0\u2500\u2524", "\u2524"][outgoing],
            ["\u2502", "\u251C\u2500\u25B6"][outgoing],
        )
    else:
        # connector chars
        vrt, vrd, relbow, relbow2, rarr = "|", ".", "+", "+", "--"
        # relation prefix and arrow
        relpre, relarr = ["<-[", "-["][outgoing], ["]-", "]->"][outgoing]
    # create map of #items at each level, so we can easily
    # know when there are more at a given level, for drawing
    n_at_level = {0: 0}
    for item in rr:
        depth = item[0]
        if depth in n_at_level:
            n_at_level[depth] += 1
        else:
            n_at_level[depth] = 1
    # print tree
    while q:
        depth, rel, meta = q.pop()
        n_at_level[depth] -= 1
        indent = "".join(
            [
                f" {t.blue}{vrd if n_at_level[i - 1] else ' '}{t.resetc} "
                for i in range(1, depth + 1)
            ]
        )
        print(f"{indent} {t.blue}{vrt}")
        rstr = f"{t.blue}{relpre}{t.yellow}{rel.predicate}" f"{t.blue}{relarr}{t.reset}"
        if meta["aliases"]:
            item_name = meta["aliases"][0]
        else:
            item_name = meta.get("desc", "-")
        istr = _related_item(
            meta[resource.Resource.ID_FIELD], item_name, meta["type"], pfx, t, unicode
        )
        # determine correct connector (whether there is another one down the stack)
        elbow = relbow if (not q or q[-1][0] != depth) else relbow2
        print(f"{indent} {t.blue}{elbow}{rarr}{t.resetc}{rstr} {istr}")
        new_rr = []
        for d2, rel2, _ in rr:
            if outgoing:
                is_same = rel2.subject == rel.object
            else:
                is_same = rel2.object == rel.subject
            if d2 == depth + 1 and is_same:
                q.append((d2, rel2, _))
            else:
                new_rr.append((d2, rel2, _))
        rr = new_rr


def _related_item(id_, name, type_, pfx, t, unicode):
    return f"{t.red}{id_[:pfx]} {t.green}{type_} {t.resetc}{name}"


@click.command(help="Remove a resource")  # aliases: delete
@click.argument("identifier", nargs=-1)
@click.option(
    "-y",
    "--yes",
    flag_value="yes",
    help="No interactive confirmations; assume 'yes' answer to all",
)
@click.option("--list/--no-list", "list_resources", default=True)
@click.option("--multiple/--no-multiple", default=False)
def rm(identifier, yes, multiple, list_resources):
    for ident in identifier:
        _log.info(f"remove resource '{ident}'")
        try:
            resource.identifier_str(ident, allow_prefix=True)
        except ValueError as errmsg:
            click.echo(f"Invalid identifier. Details: {errmsg}")
            sys.exit(Code.INPUT_VALUE.value)
        rsrc_list = list(find_by_id(ident))
        found_multiple = len(rsrc_list) > 1
        if found_multiple and not multiple:
            click.echo(
                f"Too many ({len(rsrc_list)}) resources match prefix '{ident}'. "
                "Add option --multiple to allow multiple matches."
            )
            sys.exit(Code.DMF_OPER.value)
        fields = ["type", "desc", "modified"]  # "id" is prepended by _ls_basic()
        if list_resources:
            _print_resource_table(rsrc_list, fields, ["id"], False, False, True)
        if yes != "yes":
            if found_multiple:
                s = f"these {len(rsrc_list)} resources"
            else:
                s = "this resource"
            do_remove = click.confirm(f"Remove {s}", prompt_suffix="? ", default=False)
            if not do_remove:
                click.echo("aborted")
                sys.exit(Code.CANCELED.value)
        d = DMF()
        for r in rsrc_list:
            _log.debug(f"begin remove-resource id={r.id}")
            d.remove(identifier=r.id)
            _log.debug(f"end remove-resource id={r.id}")
        if found_multiple:
            s = f"{len(rsrc_list)} resources removed"
        else:
            s = "resource removed"
        click.echo(s)


@click.command(help="Load a directory of data and associated reference")
@click.option(
    "-d",
    "--datadir",
    default=None,
    help="Data & configuration directory (default is current working directory)",
)
@click.option(
    "--global/--no-global",
    "global_",
    default=False,
    help="Load into global IDAES data workspace (by default, use current workspace)",
)
def load_data(datadir, global_):
    from idaes.core.dmf import datasets

    def echolog(msg, error=False):
        if error:
            _log.error(f"Error: {msg}")
        else:
            _log.info(msg)
        click.echo(msg)

    # Set DMF workspace
    if global_:
        # Setting workspace to None forces the default, global, one
        data_workspace = None
        # logging
        echolog(f"Load data into data workspace '{datasets.get_dataset_workspace()}'")
    else:
        # Query DMF to get current workspace
        try:
            data_workspace = DMFConfig().workspace
        except (IOError, ValueError) as err:
            echolog(f"Cannot retrieve current workspace: {err}", error=True)
            sys.exit(-1)
        echolog(f"Load data into current DMF workspace at '{data_workspace}'")

    # Set input data directory
    if datadir is None:
        data_directory = pathlib.Path(".").absolute()
    else:
        data_directory = pathlib.Path(datadir).absolute()
    echolog(f"Load data from {data_directory}")

    # Load the dataset
    pd = datasets.PublicationDataset(workspace=data_workspace)
    try:
        pd.load(data_directory)
    except datasets.ConfigurationError as err:
        echolog(f"Bad configuration. Details: {err}", error=True)
        sys.exit(1)
    except datasets.FileMissingError as err:
        echolog(f"Bad input. Details: {err}", error=True)
        sys.exit(2)

    echolog(f"Loaded data in '{data_directory}' into the DMF")


######################################################################################


def find_by_id(identifier, dmf=None):
    if dmf is None:
        dmf = DMF()
    return dmf.find_by_id(identifier)


class _ShowInfo:
    """Container for methods, etc. to show info about a resource."""

    contents_indent, json_indent = 4, 2  # for `term` output

    def __init__(self, output_format, pfxlen, color=None, unicode=True):
        self._terminal = ColorTerm(enabled=color)
        self._pfxlen = pfxlen
        self._fmt = output_format
        self._resource = None
        C = namedtuple("Corners", ["nw", "ne", "se", "sw"])
        if unicode:
            self._corners = C._make("\u250C\u2510\u2518\u2514")
            self._hz, self._vt = "\u2500", "\u2502"
        else:
            self._corners = C._make("++++")
            self._hz, self._vt = "-", "|"

    def show(self, rsrc):
        self._resource = rsrc
        getattr(self, f"show_{self._fmt}")()

    def show_json(self):
        json.dump(self._resource.v, sys.stdout, indent=self.json_indent)
        print()

    def show_jsonc(self):
        json.dump(self._resource.v, sys.stdout)
        print()

    def show_term(self):
        t = self._terminal
        rval = self._human_readable_values()
        term_width = t.width
        width = min(term_width, self._longest_line(rval) + 3 + self.contents_indent)
        self._print_info_term_header(width)
        top_keys = sorted(rval.keys())
        for tk in top_keys:
            val = rval[tk]
            if self._has_values(val):
                contents_str = yaml.dump(
                    yaml_load(json.dumps(val)),
                    default_flow_style=False,
                    explicit_end=False,
                )
                if contents_str.endswith("...\n"):
                    contents_str = contents_str[:-4]
                print(
                    f"{t.cyan}{self._vt}{t.resetc} {t.bold}{t.cyan}{tk}{t.reset}"
                    f"{' ' * (width - len(tk) - 3)}{t.cyan}{self._vt}{t.resetc}"
                )
                self._print_info_contents_term(contents_str, width)
        print(
            f"{t.cyan}{self._corners.sw}{self._hz * (width - 2)}"
            f"{self._corners.se}{t.resetc}"
        )

    def _longest_line(self, formatted_resource):
        longest = 0
        for v in formatted_resource.values():
            v_longest = max(
                (len(s) for s in json.dumps(v, indent=self.json_indent).split("\n"))
            )
            longest = max((longest, v_longest))
        return longest

    def _print_info_term_header(self, width):
        t, r = self._terminal, self._resource
        padding = width - len(r.id) - 10
        if padding % 2 == 1:
            lpad = padding // 2 - 1
            rpad = padding // 2 - 1
        else:
            lpad, rpad = padding // 2 - 1, padding // 2 - 2
        print(
            f"{t.cyan}{self._corners.nw}{lpad * self._hz}{t.resetc} Resource "
            f"{t.red}{r.id[:self._pfxlen]}"
            f"{t.green}{r.id[self._pfxlen:]} "
            f"{t.cyan}{rpad * self._hz}{self._corners.ne}{t.resetc}"
        )

    def _print_info_contents_term(self, s, width):
        t, n = self._terminal, self.contents_indent
        contents_width = width - 2 - n
        indent = " " * n
        for line in s.split("\n"):
            p = 0
            while p < len(line):
                print_line = line[p : p + contents_width]
                print(
                    f"{t.cyan}{self._vt}{t.resetc}{indent}{print_line}"
                    f"{' ' * (contents_width - len(print_line))}{t.cyan}{self._vt}"
                    f"{t.resetc}"
                )
                p += contents_width

    def _human_readable_values(self):
        val, result = self._resource.v, {}

        def dateize(v):
            return datetime.isoformat(datetime.fromtimestamp(v))

        for k, v in val.items():
            if k in ("created", "modified"):
                result[k] = dateize(v)
            elif k == "version_info":
                result[
                    "version"
                ] = f"{resource.format_version(v['version'])} @ {dateize(v['created'])}"
            elif k == "relations":
                relations = []
                for rel in v:
                    is_object = rel["role"] == "object"
                    predicate = f"--[{rel['predicate']}]--"
                    if is_object:
                        s = f"{rel['identifier']} {predicate}> ME"
                    else:
                        s = f"{rel['identifier']} <{predicate} ME"
                    relations.append(s)
                if relations:
                    result[k] = relations
            else:
                if isinstance(v, list):
                    result[k] = v[:]
                elif isinstance(v, dict):
                    result[k] = v.copy()
                else:
                    result[k] = v
        return result

    @staticmethod
    def _has_values(data):
        if isinstance(data, list) or isinstance(data, dict):
            return bool(data)
        return True


class _LsField:
    """Subclasses define how to retrieve and transform the
    user-specified field into a string for display by 'ls'.
    """

    def __init__(self, key):
        self._key = key
        self.value = ""

    def set_value(self, rsrc):
        if isinstance(self._key, list) or isinstance(self._key, tuple):
            v = rsrc.v
            for k in self._key:
                v = v[k]
            self.value = v
        else:
            self.value = rsrc.v[self._key]


class _IdentityField(_LsField):
    def __str__(self):
        return str(self.value)


class _DateField(_LsField):
    def __str__(self):
        return datetime.isoformat(datetime.fromtimestamp(float(self.value)))


class _FilesField(_LsField):
    def __str__(self):
        if len(self.value) > 1:
            return f"{self.value[0]['path']} ..."
        elif len(self.value) == 1:
            return f"{self.value[0]['path']}"
        return "<none>"


class _CodesField(_LsField):
    def __str__(self):
        if len(self.value) > 1:
            return f"{self.value[0]['name']} ..."
        elif len(self.value) == 1:
            return f"{self.value[0]['name']}"
        return "<none>"


class _IdField(_LsField):
    def __init__(self, key):
        super().__init__(key)
        self.pfxlen = 0
        self.flen = 0
        self.term = None

    def __str__(self):
        t = self.term
        if self.pfxlen < len(self.value):
            if self.flen > self.pfxlen:
                return (
                    f"{t.cyan}{self.value[: self.pfxlen]}{t.resetc}"
                    f"{self.value[self.pfxlen: self.flen]}"
                )
            else:
                return self.value[: self.pfxlen]
        return self.value


class _VersionField(_LsField):
    def __str__(self):
        return resource.format_version(self.value)


class _NameField(_LsField):
    def __init__(self):
        super().__init__("aliases")
        self.value = []

    def __str__(self):
        if len(self.value) == 0:
            return ""
        return self.value[0]


# Mapping from the field name to an instance of a subclass
# of _Field that can extract sortable and formatted values.


_show_fields = {
    "id": _IdField(resource.Resource.ID_FIELD),
    "type": _IdentityField("type"),
    "desc": _IdentityField("desc"),
    "creator": _IdentityField(("creator", "name")),
    "created": _DateField("created"),
    "modified": _DateField("modified"),
    "files": _FilesField("datafiles"),
    "codes": _CodesField("codes"),
    "version": _VersionField(("version_info", "version")),
    "name": _NameField(),
}


######################################################################################

# Register base commands
base_command.add_command(setup)
base_command.add_command(init)
base_command.add_command(register)
base_command.add_command(status)
base_command.add_command(ls)
base_command.add_command(find)
base_command.add_command(info)
base_command.add_command(related)
base_command.add_command(rm)
base_command.add_command(load_data)

# if __name__ == '__main__':
#     base_command()
