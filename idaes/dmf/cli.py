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
import json
import logging
from operator import itemgetter
import pathlib
import sys
from urllib.parse import urlparse, ParseResult
import yaml

# third-party
from blessings import Terminal
import click
import humanize
import pendulum

# package
from idaes.dmf import DMF, DMFConfig, resource
from idaes.dmf import errors
from idaes.dmf.workspace import Fields

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
    INPUT_VALUE = 7


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


@click.group(
    cls=AliasedGroup,
    aliases={
        "describe": "status",
        "add": "register",
        "resource": "info",
        "show": "info",
        "list": "ls",
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
        _log.setLevel(level_from_verbosity(verbose))
    else:
        _log.setLevel(level_from_verbosity(-quiet))


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
    _log.info(f"Initialize with workspace path={path}")
    if create:
        _log.info("Create new workspace")
        # pre-check that there is no file/dir by this name
        if pathlib.Path(path).exists():
            click.echo(f"Cannot create workspace: path '{path}' already exists")
            sys.exit(Code.DMF_OPER.value)
        if not name:
            name = click.prompt("New workspace name")
        if not desc:
            desc = click.prompt("New workspace description")
        if html is None:
            # guess html path
            # XXX: don't try to verify the guess
            errfile = pathlib.Path(errors.__file__)
            docsdir = errfile.parent.parent.parent / 'docs'
            hpath = [str(docsdir / 'build')]
        else:
            hpath = [html]
        try:
            d = DMF(path=path, create=True, name=name, desc=desc, html_paths=hpath)
        except errors.WorkspaceError as err:
            click.echo(f"Cannot create workspace: {err}")
            sys.exit(Code.DMF_OPER.value)
        click.echo(f"Configuration in '{d.configuration_file}")
    else:
        _log.info("Use existing workspace")
        try:
            _ = DMF(path=path, create=False, save_path=True)
        except errors.WorkspaceConfNotFoundError:
            click.echo(f"Workspace configuration not found at path='{path}'")
            if path == '.':  # probably just the default
                click.echo("Use --path option to set workspace path.")
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
        if show:
            click.echo(f"note: option `--all` overrides `--show`")
        show = ["all"]
    _log.debug(f"Get status. Show items: {' '.join(show) if show else '(basic)'}")
    t = _cterm if color else _noterm
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
            return f"{before}{color}{key}{after_key}{t.normal}"
        elif key is None:
            return f"{before}{color}{value}{t.normal}"
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
            print(
                item_fn("total_size", humanize.naturalsize(total_size), before=indent)
            )
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
def register(
    resource_type, url, info, copy, strict, unique, contained, derived, used, prev
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
                f"This file is already in {n_dup} resources: " f"{' '.join(dup_ids)}"
            )
            sys.exit(Code.DMF_OPER.value)
    # process relations
    _log.debug("add relations")
    rel_to_add = {  # translate into standard relation names
        resource.PR_CONTAINS: contained,
        resource.PR_DERIVED: derived,
        resource.PR_USES: used,
        resource.PR_VERSION: prev,
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
            resource.create_relation_args(rel_subj, rel_name, rsrc)
            _log.debug(f"added relation {rsrc.id} <-- {rel_name} -- {rel_id}")
    _log.debug("update resource relations")
    for rel_rsrc in target_resources.values():
        dmf.update(rel_rsrc)
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
@click.option(
    "--show",
    "-s",
    type=click.Choice(
        ["type", "desc", "created", "modified", "version", "files", "codes"]
    ),
    multiple=True,
    help="Show given field in listing",
)
@click.option(
    "--sort",
    "-S",
    "sort_by",
    type=click.Choice(["id", "type", "desc", "created", "modified", "version"]),
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
def ls(show, sort_by, reverse, prefix):
    d = DMF()
    if not show:
        show = ["type", "desc", "modified"]  # note: 'id' is always first
    reverse = bool(reverse == "yes")
    if not sort_by:
        sort_by = ["id"]
    _ls_basic(d, show, sort_by, reverse, prefix)


def _ls_basic(d, show_fields, sort_by, reverse, prefix):
    """Text-mode `ls`.
    """
    t = Terminal()
    resources = list(d.find())
    if len(resources) == 0:
        print("no resources to display")
        return
    uuid_len = _uuid_prefix_len([r.id for r in resources])
    full_len = max((len(r.id) for r in resources)) if not prefix else uuid_len
    fields = ["id"] + list(show_fields)
    nfields = len(fields)
    # calculate table body. do this first to get widths.
    rows, maxwid, widths = [], 60, [0] * nfields
    for r in resources:
        row = []
        for i, fld in enumerate(fields):
            # transform field for display
            transformer = _show_fields[fld]
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
    hdr_columns = [t.bold + f"{f:{w}}" for f, w in zip(fields, widths)]
    print(" ".join(hdr_columns) + t.normal)
    # print table body
    for row in rows:
        row_columns = [f"{f:{w}}" for f, w in zip(row[:nfields], widths)]
        print(" ".join(row_columns))


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
def info(identifier, multiple, output_format):
    _log.debug(f"info for resource id='{identifier}'")
    try:
        resource.identifier_str(identifier, allow_prefix=True)
    except ValueError as err:
        click.echo(f"{err}")
        sys.exit(Code.INPUT_VALUE.value)
    d = DMF()
    rsrc_list = list(d.find_by_id(identifier))
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
    si = _ShowInfo(output_format, pfxlen)
    for rsrc in rsrc_list:
        si.show(rsrc)


class _ShowInfo:
    """Container for methods, etc. to show info about a resource.
    """

    contents_indent, json_indent = 4, 2  # for `term` output

    def __init__(self, output_format, pfxlen):
        self._terminal = Terminal()
        self._pfxlen = pfxlen
        self._fmt = output_format
        self._resource = None

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
        t, r = self._terminal, self._resource  # aliases
        rval = self._human_readable_values()
        width = min(t.width, self._longest_line(rval) + 3 + self.contents_indent)
        self._print_info_term_header(width)
        top_keys = sorted(rval.keys())
        for rownum, tk in enumerate(top_keys):
            val = rval[tk]
            if self._has_values(val):
                contents_str = yaml.dump(
                    yaml.load(json.dumps(val)),
                    default_flow_style=False,
                    explicit_end=False,
                )
                if contents_str.endswith('...\n'):
                    contents_str = contents_str[:-4]
                print(
                    f"{t.on_cyan} {t.normal} {t.bold}{t.cyan}{tk}{t.normal}"
                    f"{' ' * (width - len(tk) - 3)}{t.on_cyan} {t.normal}"
                )
                self._print_info_contents_term(contents_str, width)
        print(f"{t.on_cyan}{' ' * width}{t.normal}")

    def _longest_line(self, formatted_resource):
        longest = 0
        for k, v in formatted_resource.items():
            v_longest = max(
                (len(s) for s in json.dumps(v, indent=self.json_indent).split('\n'))
            )
            longest = max((longest, v_longest))
        return longest

    def _print_info_term_header(self, width):
        t, r = self._terminal, self._resource
        padding = width - len(r.id) - 10
        if padding % 2 == 1:
            lpad = padding // 2
            rpad = padding // 2 + 1
        else:
            lpad = rpad = padding // 2
        print(
            f"{t.bold_white_on_cyan}{lpad * ' '} Resource {t.normal}"
            f"{t.black_on_cyan}{r.id[:self._pfxlen]}"
            f"{t.white_on_cyan}{r.id[self._pfxlen:]}"
            f"{rpad * ' '}{t.normal}"
        )

    def _print_info_contents_term(self, s, width):
        t, n = self._terminal, self.contents_indent
        contents_width = width - 2 - n
        indent = " " * n
        for line in s.split('\n'):
            p = 0
            while p < len(line):
                print_line = line[p : p + contents_width]
                print(
                    f"{t.on_cyan} {t.normal}{indent}{print_line}"
                    f"{' ' * (contents_width - len(print_line))}{t.on_cyan} {t.normal}"
                )
                p += contents_width

    def _human_readable_values(self):
        val, result = self._resource.v, {}

        def dateize(v):
            return pendulum.from_timestamp(v).to_datetime_string()

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
        return pendulum.from_timestamp(self.value).to_datetime_string()


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
                    f"{t.cyan}{self.value[: self.pfxlen]}{t.normal}"
                    f"{self.value[self.pfxlen: self.flen]}"
                )
            else:
                return self.value[: self.pfxlen]
        return self.value


class _VersionField(_LsField):
    def __str__(self):
        return resource.format_version(self.value)


# Mapping from the field name to an instance of a subclass
# of _Field that can extract sortable and formatted values.

_show_fields = {
    "id": _IdField(resource.Resource.ID_FIELD),
    "type": _IdentityField("type"),
    "desc": _IdentityField("desc"),
    "created": _DateField("created"),
    "modified": _DateField("modified"),
    "files": _FilesField("datafiles"),
    "codes": _CodesField("codes"),
    "version": _VersionField(("version_info", "version")),
}


def _uuid_prefix_len(uuids, step=4, maxlen=32):
    """Get smallest multiple of `step` len prefix that gives unique values.

    The algorithm is not fancy, but good enough: build *sets* of
    the ids at increasing prefix lengths until the set has all ids (no duplicates).
    Experimentally this takes ~.1ms for 1000 duplicate ids (the worst case).
    """
    full = set(uuids)
    all_of_them = len(full)
    for n in range(step, maxlen, step):
        prefixes = {u[:n] for u in uuids}
        if len(prefixes) == all_of_them:
            return n
    return maxlen


# Register base commands
base_command.add_command(init)
base_command.add_command(register)
base_command.add_command(status)
base_command.add_command(ls)
base_command.add_command(info)

if __name__ == '__main__':
    base_command()
