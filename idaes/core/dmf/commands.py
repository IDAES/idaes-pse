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
Perform all logic, input, output of commands that is
particular to the CLI.

Call functions defined in 'api' module to handle logic
that is common to the API and CLI.
"""
# stdlib
from datetime import datetime
import glob
import json
import logging
import math
import os
import re
import sys

# Third-party
import jsonschema

# Local
from .dmfbase import DMF, DMFConfig
from .util import strlist
from .util import is_jupyter_notebook, is_python, is_resource_json
from .util import ColorTerm
from .errors import (
    ParseError,
    CommandError,
    WorkspaceNotFoundError,
    WorkspaceConfNotFoundError,
    WorkspaceError,
    BadResourceError,
)
from . import resource
from .workspace import Workspace, find_workspaces

__author__ = "Dan Gunter"

_log = logging.getLogger(__name__)


def workspace_init(dirname, metadata):
    # type: (str, dict) -> None
    """Initialize from root at `dirname`, set environment variable
    for other commands, and parse config file.
    """
    try:
        ws = Workspace(dirname, create=True, add_defaults=True)
    except OSError as err:
        raise CommandError("init", "initialize workspace", str(err))
    except ParseError as err:
        raise CommandError("init", "parse config", str(err))
    except WorkspaceError as err:
        raise CommandError("init", "initialize workspace", str(err))
    _log.info("Created new workspace in: {}".format(dirname))
    if metadata:
        ws.set_meta(metadata)
        _log.info("Set metadata for: {}".format(strlist(list(metadata))))


def workspace_info(dirname):
    # type: (str) -> None
    t = ColorTerm()
    try:
        ws = Workspace(dirname, create=False)
    except WorkspaceNotFoundError:
        print("Workspace not found at path: {}".format(dirname))
        raise CommandError("info", "find workspace", "not found at: {}".format(dirname))
    except WorkspaceConfNotFoundError:
        print("No configuration found for workspace for path: {}".format(dirname))
        raise CommandError(
            "info", "find workspace configuration", "not found at: {}".format(dirname)
        )
    num_obj = DMF(path=ws.root).count()
    bullet = " - "
    print(f"\n{t.blue}Workspace")
    if ws.name and (ws.name != "none"):
        if ws.description and (ws.description != "none"):
            print(f"  {t.blue}[{ws.name}] - {ws.description}")
        else:
            print("  {t.blue}{ws.name} - (no description)")
    elif ws.description and (ws.description != "none"):
        print(f"  {t.blue}(no name) - {ws.description}")
    else:
        print(f"  {t.blue}(no name or description)")
    print("\nGeneral information")
    print(f"{bullet}{t.blue}Location = {ws.root}")
    info = ws.meta.copy()
    if "_id" in info:
        print(f"{bullet}{t.blue}Workspace identifier (_id) = {info['_id']}")
        del info["_id"]
    else:
        print(f"{bullet}{t.blue}Workspace identifier (_id) = unknown")
    if "created" in info:
        print(f"{bullet}{t.blue}Created = {info[ws.CONF_CREATED]}")
    else:
        print(f"{bullet}{t.blue}Created = unknown")
    if "modified" in info:
        print(f"{bullet}{t.blue}Modified = {info[ws.CONF_MODIFIED]}")
    else:
        print(f"{bullet}{t.blue}Modified = unknown")
    print(f"{bullet}{t.blue}Num. resources = {num_obj}")
    print(f"\n{t.magenta}{t.bold}Configuration")
    already_shown = (ws.CONF_MODIFIED, ws.CONF_CREATED, ws.CONF_NAME, ws.CONF_DESC)
    for k in info.keys():
        if k in already_shown:
            continue
        v = info[k]
        print(f"{bullet}{t.blue}{k} = {v}")
    print("")


def init_conf(workspace):
    # type: (str) -> int
    """Initialize the workspace."""
    t = ColorTerm()
    # Open/create configuration file
    try:
        conf = DMFConfig()
    except IOError as err:
        print(f"Failed to open global configuration: {err}")
        try:
            open(DMFConfig._filename, "w")
        except IOError:
            print("Failed to create new configuration file")
            return -1
        print("Created new configuration file")
        conf = DMFConfig()
    # If a workspace argument is given, save this value,
    # as the default workspace, in the configuration file
    if workspace:
        fullpath = os.path.abspath(workspace)
        conf.c[conf.WORKSPACE] = fullpath
        conf.save()
    # Print contents of configuration file to standard output
    print(
        f"{t.magenta}{t.bold}DMF global configuration{t.reset} "
        f"<{t.green}{conf._filename}>"
    )
    keys = conf.c.keys()
    if keys:
        for k in sorted(keys):
            print(f" > {t.blue}{k}{t.reset} = {t.bold}{conf.c[k]}]")
    else:
        print(f"{t.blue}(empty)")
    return 0


def workspace_import(path, patterns, exit_on_error):
    # type: (str, list, bool) -> int
    """Import files into workspace.

    Args:
        path (str): Target workspace directory
        patterns (list): List of Unix-style glob for files to import.
                       Files are expected to be resource JSON or a
                       Jupyter Notebook.
        exit_on_error (bool): If False, continue trying to import resources
                              even if one or more fail.

    Returns:
        int: Number of things imported
    Raises:
        BadResourceError, if there is a problem
    """
    d = DMF(path)
    count = 0
    for pattern in patterns:
        for filename in glob.glob(pattern):
            # Skip directories
            if os.path.isdir(filename):
                _log.warning('Not importing directory "{}"'.format(filename))
                continue
            # For Jupyter Notebooks, first create a (temporary)
            # JSON resource from the original data.
            if is_jupyter_notebook(filename):
                try:
                    rsrc = _import_jupyternb(filename)
                except ValueError as e:
                    msg = (
                        "Cannot create resource from Jupyter Notebook "
                        '"{}": {}'.format(filename, e)
                    )
                    if exit_on_error:
                        raise BadResourceError(msg)
                    _log.error(msg)
                    continue
            # For Python files, first create a (temporary)
            # JSON resource from the original data.
            elif is_python(filename):
                try:
                    rsrc = _import_python(filename)
                except ValueError as e:
                    msg = "Cannot create resource from Python file " '"{}": {}'.format(
                        filename, e
                    )
                    if exit_on_error:
                        raise BadResourceError(msg)
                    _log.error(msg)
                    continue
            # JSON resource file
            elif is_resource_json(filename):
                try:
                    rsrc = _import_resource(filename)
                except ValueError as e:
                    msg = 'Bad resource from file "{}": {}'.format(filename, e)
                    if exit_on_error:
                        raise BadResourceError(msg)
                    _log.error(msg)
                    continue
            # Generic file import
            else:
                try:
                    rsrc = _import_file(filename)
                except ValueError as e:
                    msg = "Cannot create resource from file " '"{}": {}'.format(
                        filename, e
                    )
                    if exit_on_error:
                        raise BadResourceError(msg)
                    _log.error(msg)
                    continue
            # Resource in hand. Now add it.
            d.add(rsrc)
            count += 1
        return count


def list_workspaces(root, stream=None):
    """List workspaces found from a given root path.

    Args:
        root: root path
        stream: Output stream (must have .write() method)
    """
    workspaces = find_workspaces(root)
    if stream is None or stream == sys.stdout:
        colors = True
    else:
        colors = False
    t = ColorTerm(enabled=colors)
    if colors:
        output_table = [("Path", "Name")]
    else:
        output_table = [("Path", "Name"), ("----", "----")]
    widths = [4, 4]
    any_good_workspaces = False
    for w in sorted(workspaces):
        try:
            ws = Workspace(w)
            output_table.append((w, ws.name))
            widths = [max(len(w), widths[0]), max(len(ws.name), widths[1])]
            any_good_workspaces = True
        except WorkspaceError:
            pass  # XXX: Should we print a warning?
    if not any_good_workspaces:
        # either no paths, or all paths raised an error
        stream.write("ERROR: No valid workspaces found\n")
    else:
        colfmts = ["{{:{:d}s}}".format(width) for width in widths]
        first_row = True
        for row in output_table:
            for i in (0, 1):
                if colors:
                    if first_row:
                        fmt = f"{t.bold}{colfmts[i]}"
                    else:
                        fmt = f"{[t.blue, t.white][i]}{colfmts[i]}"
                    fmt += t.reset
                else:
                    fmt = colfmts[i]
                stream.write(fmt.format(row[i]))
                stream.write("\n" if i == 1 else " ")
            first_row = False


def list_resources(path, long_format=None, relations=False):
    """List resources in a given DMF workspace.

    Args:
        path (str): Path to the workspace
        long_format (bool): List in long format flag
        relations (bool): Show relationships, in long format

    Returns:
        None
    """
    t = ColorTerm()
    d = DMF(path)
    if long_format:
        resources = list(d.find())
        uuid_pfx = _uuid_prefix([r.uuid for r in resources])
        fields = ("uuid", "name", "type", "modified", "created")
        widths = (uuid_pfx, 30, 20, 19, 19)
        colors = (t.green, t.white, t.yellow, t.white, t.white)
        fmts = [f"{{:{w}s}}" for w in widths]
        left_gutter = "| " if relations else ""
        # table header
        print(
            " " * len(left_gutter)
            + t.bold
            + "  ".join([f.format(v) for f, v in zip(fmts, fields)])
            + t.reset
        )

        def datestr(t):
            return datetime.isoformat(datetime.fromtimestamp(t))

        # table body
        for r in resources:
            values = list(getattr(r, k) for k in fields[:-2])
            values.append(datestr(r.modified))
            values.append(datestr(r.created))
            if not values[1] and r.desc:
                values[1] = r.desc[: widths[1]]
            else:
                values[1] = values[1][: widths[1]]
            if uuid_pfx < 32:
                values[0] = values[0][:uuid_pfx]
            print(
                left_gutter
                + "  ".join([c + f.format(v) for c, f, v in zip(colors, fmts, values)])
                + t.reset
            )
            if relations and len(r.relations) > 0:
                relitems = []
                for rel in r.relations:
                    if rel.subject == r.uuid:
                        fmt = f"{t.white}{{p}}->{t.blue}{{o}}"
                    else:
                        fmt = f"{t.blue}{{s}}->{t.white}{{p}}"
                    item = fmt.format(
                        s=rel.subject[:uuid_pfx],
                        p=rel.predicate,
                        o=rel.object[:uuid_pfx],
                    )
                    relitems.append(item)
                print(f"+-- {' / '.join(relitems)}")
    else:
        items = []
        for r in d.find():
            name_color = "w"
            if r.name:
                name = r.name
            elif r.desc:
                name = r.desc[:40]
                name_color = t.blue
            else:
                name = r.uuid
                name_color = t.green
            item = f"{name_color}{name}{t.yellow}:{r.type}"
            items.append(item)
        if items:
            columnized = _display_in_columns(items, max_line=t.width)
            print(columnized + t.reset)


def _uuid_prefix(uuids, step=4, maxlen=32):
    """Get smallest multiple of `step` len prefix that gives unique values."""
    full = set(uuids)
    for n in range(step, maxlen, step):
        prefixes = {u[:n] for u in uuids}
        if len(prefixes) == len(full):
            return n
    return maxlen


def cat_resources(path, objects=(), color=True):
    d = DMF(path=path)
    t = ColorTerm(enabled=color)
    unmatched = set(objects)
    first = True
    # get all resources,
    # display any that match an object as a prefix
    for r in d.find():
        for oid in unmatched:
            if r.uuid.startswith(oid):
                unmatched.remove(oid)  # don't show twice
                if not first:
                    _cat_resource_sep(t)
                _cat_resource_show(t, r)
                first = False
                break


def _cat_resource_sep(t):
    print(f"{t.blue}{'-' * 60}")


def _cat_resource_show(cp, r):
    d = r.as_dict()
    json.dump(d, cp, indent=2)
    print()


# regular expression to find VT100 color escape codes
_noprint_re = re.compile(r"\033\[[0-9]+m")


def _display_in_columns(items, max_line=80, col_sep="  ", row_sep="\n"):
    """Take a list of items and max line width, and calculate  display
    of the items in columns.

    The algorithm is simple, just trying increasing numbers of columns and
    picking the largest number that did not result in a row that was too wide.

    The input items are not re-ordered.

    Args:
        items (List[str]): String items
        max_line (int): Maximum width for any displayed line (row)
        col_sep (str): Separator between columns, after each item
        row_sep (str): Separator between rows, at the end of each line

    Returns:
        str:
    """
    if not items:
        return ""
    # Calculate item lengths, stripping terminal escapes
    lengths, nplengths = [], []
    for item in items:
        clean = _noprint_re.sub("", item)
        lengths.append(len(clean))
        nplengths.append(len(item) - len(clean))
    col_sep_len = len(col_sep)  # useful later
    # Give up immediately, putting everything in one column,
    # if any single item doesn't fit
    if max_line <= max(lengths) + col_sep_len:
        return row_sep.join(items)
    # Determine maximum number of columns
    max_columns = 1  # number of columns
    max_widths = [max_line]  # width of each column
    n = len(lengths)
    # Determine number of columns.
    # Start at 2 columns, stop when cannot fit items side-by-side
    for i in range(2, n):
        # initialize calculated widths of each column
        widths = [0] * i
        # for display where all columns are same length except last,
        # number of items per column is ceiling of items/#col
        nitems = int(math.ceil(n / i))
        # put items in each column
        for col in range(i):
            pad = 0 if col == (i - 1) else col_sep_len  # sep. between columns
            # put items in the current column, adjusting the column
            # max width to widest item
            maxj = min(n, (col + 1) * nitems)  # don't overshoot on last col
            for j in range(col * nitems, maxj):
                widths[col] = max(widths[col], lengths[j] + pad)
        # total width is sum of column widths
        line_len = sum(widths)
        # if we went over, then stop
        if line_len > max_line:
            break
        # otherwise, this is valid -- save and continue
        max_columns, max_widths = i, widths[:]
    # Put items into rows of max. number of columns determined above
    nrows, rows = int(math.ceil(n / max_columns)), []
    for row in range(nrows):
        col, row_items = 0, []
        # skip through items by nrows at a time, to move across the columns,
        # and fill in the items for the current row (which acts as an offset)
        for i in range(row, len(items), nrows):
            # create format string with width = col. max including esc chars,
            # but without padding since we will add that when we join
            # the row items together
            pad = 0 if col == (max_columns - 1) else col_sep_len
            fmt = "{{:{n}s}}".format(n=max_widths[col] + nplengths[i] - pad)
            # format row item for column
            row_items.append(fmt.format(items[i]))
            col += 1  # move to next column
        # append the row items as big string
        rows.append(col_sep.join(row_items))
    # Final result is a big string of the rows joined together
    return row_sep.join(rows)


def _import_resource(filename):
    """Import a resource from 'filename'. Raises a ValueError if that
    fails. Most of the code is simply generating error messages.
    """
    if not os.path.exists(filename):
        raise ValueError('File "{}" not found'.format(filename))
    try:
        f = open(filename)
    except Exception as e:
        raise ValueError('Cannot open file "{}": {}'.format(filename, e))
    try:
        j = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError('Cannot parse JSON file "{}": {}'.format(filename, e))
    try:
        r = resource.Resource(value=j)
        r.validate()
    except (ValueError, jsonschema.ValidationError) as err:
        raise ValueError("Invalid resource: {}".format(err))
    return r


def _import_jupyternb(path):
    """Create & import a resource from a Jupyter Notebook file at `path`.
    Assume that `path` exists and is a Jupyter Notebook.

    Args:
        path (str): Jupyter Notebook file.
    Returns:
          (Resource) DMF Resource representing the notebook.
    """
    r = resource.Resource(type_=resource.ResourceTypes.notebook)
    filename = os.path.splitext(os.path.split(path)[1])[0]
    # XXX: add notebook 'metadata' as FilePath metadata attr
    r.v["datafiles"].append({"desc": filename, "path": path})
    r.v["desc"] = filename
    r.validate()
    return r


def _import_python(path):
    """Create & import a resource from a Python file at `path`.
    Assume that `path` exists and is a valid Python file.

    Args:
        path (str): Python file name.
    Returns:
          (Resource) DMF Resource representing the notebook.
    """
    r = resource.Resource(type_=resource.ResourceTypes.code)
    filename = os.path.splitext(os.path.split(path)[1])[0]
    r.v["codes"].append({"name": filename, "language": "python", "type": "module"})
    r.v["datafiles"].append({"desc": filename, "path": path})
    r.validate()
    return r


def _import_file(path):
    """Create & import a resource from a generic file at `path`.
    Assume that `path` exists.

    Args:
        path (str): File name.
    Returns:
          (Resource) DMF Resource representing the notebook.
    """
    r = resource.Resource(type_=resource.ResourceTypes.data)
    filename = os.path.split(path)[1]
    r.v["datafiles"].append({"desc": filename, "path": path})
    r.v["desc"] = filename
    r.validate()
    return r
