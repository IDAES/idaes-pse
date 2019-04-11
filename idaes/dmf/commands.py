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
Perform all logic, input, output of commands that is
particular to the CLI.

Call functions defined in 'api' module to handle logic
that is common to the API and CLI.
"""
# stdlib
import glob
import json
import logging
import math
import os
import re
import sys

# Third-party
from backports.shutil_get_terminal_size import get_terminal_size
import jsonschema
import pendulum

# Local
from .dmfbase import DMF, DMFConfig
from .util import strlist
from .util import is_jupyter_notebook, is_python, is_resource_json
from .util import CPrint
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

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'

_log = logging.getLogger(__name__)


def workspace_init(dirname, metadata):
    # type: (str, dict) -> None
    """Initialize from root at `dirname`, set environment variable
    for other commands, and parse config file.
    """
    try:
        ws = Workspace(dirname, create=True, add_defaults=True)
    except OSError as err:
        raise CommandError('init', 'initialize workspace', str(err))
    except ParseError as err:
        raise CommandError('init', 'parse config', str(err))
    except WorkspaceError as err:
        raise CommandError('init', 'initialize workspace', str(err))
    _log.info('Created new workspace in: {}'.format(dirname))
    if metadata:
        ws.set_meta(metadata)
        _log.info('Set metadata for: {}'.format(strlist(list(metadata))))


def workspace_info(dirname):
    # type: (str) -> None
    cp = CPrint()
    try:
        ws = Workspace(dirname, create=False)
    except WorkspaceNotFoundError:
        print('Workspace not found at path: {}'.format(dirname))
        raise CommandError('info', 'find workspace', 'not found at: {}'.format(dirname))
    except WorkspaceConfNotFoundError:
        print('No configuration found for workspace for path: {}'.format(dirname))
        raise CommandError(
            'info', 'find workspace configuration', 'not found at: {}'.format(dirname)
        )
    num_obj = DMF(path=ws.root).count()
    bullet = ' - '
    cp('\n@h[Workspace]')
    if ws.name and (ws.name != 'none'):
        if ws.description and (ws.description != 'none'):
            cp('  @b[{}] - {}'.format(ws.name, ws.description))
        else:
            cp('  @b[{}] - (no description)'.format(ws.name))
    elif ws.description and (ws.description != 'none'):
        cp('  @b[(no name)] - {}'.format(ws.description))
    else:
        cp('  @b[(no name or description)]')
    cp('\n@h[General information]')
    cp('{}@b[Location] = {}'.format(bullet, ws.root))
    info = ws.meta.copy()
    if '_id' in info:
        cp('{}@b[Workspace identifier (_id)] = {}'.format(bullet, info['_id']))
        del info['_id']
    else:
        cp('{}@b[Workspace identifier (_id)] = unknown'.format(bullet))
    if 'created' in info:
        cp('{}@b[Created] = {}'.format(bullet, info[ws.CONF_CREATED]))
    else:
        cp('{}@b[Created] = unknown'.format(bullet))
    if 'modified' in info:
        cp('{}@b[Modified] = {}'.format(bullet, info[ws.CONF_MODIFIED]))
    else:
        cp('{}@b[Modified] = unknown'.format(bullet))
    cp('{}@b[Num. resources] = {:d}'.format(bullet, num_obj))
    cp('\n@h[Configuration]')
    already_shown = (ws.CONF_MODIFIED, ws.CONF_CREATED, ws.CONF_NAME, ws.CONF_DESC)
    for k in info.keys():
        if k in already_shown:
            continue
        v = info[k]
        cp('{}@b[{}] = {}'.format(bullet, k, v))
    cp('')


def init_conf(workspace):
    # type: (str) -> int
    """Initialize the workspace.
    """
    # Open/create configuration file
    try:
        conf = DMFConfig()
    except IOError as err:
        print('Failed to open global configuration: {}'.format(err))
        try:
            open(DMFConfig.filename, 'w')
        except IOError:
            print('Failed to create new configuration file')
            return -1
        print('Created new configuration file')
        conf = DMFConfig()
    # If a workspace argument is given, save this value,
    # as the default workspace, in the configuration file
    if workspace:
        fullpath = os.path.abspath(workspace)
        conf.c[conf.WORKSPACE] = fullpath
        conf.save()
    # Print contents of configuration file to standard output
    cp = CPrint()
    cp('@h[DMF global configuration] <@g[{}]>'.format(conf.filename))
    keys = conf.c.keys()
    if keys:
        for k in sorted(keys):
            cp(' > @b[{}] = @*[{}]'.format(k, conf.c[k]))
    else:
        print('@b[(empty)]')
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
                        'Cannot create resource from Jupyter Notebook '
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
                    msg = 'Cannot create resource from Python file ' '"{}": {}'.format(
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
                    msg = 'Cannot create resource from file ' '"{}": {}'.format(
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
        stream = CPrint()
        colors = True
    else:
        colors = False
    if colors:
        output_table = [('Path', 'Name')]
    else:
        output_table = [('Path', 'Name'), ('----', '----')]
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
        stream.write('ERROR: No valid workspaces found\n')
    else:
        colfmts = ['{{:{:d}s}}'.format(width) for width in widths]
        first_row = True
        for row in output_table:
            for i in (0, 1):
                if colors:
                    if first_row:
                        fmt = '@_h[{t}]'.format(t=colfmts[i])
                    else:
                        fmt = '@{c}[{t}]'.format(c=('b', 'w')[i], t=colfmts[i])
                else:
                    fmt = colfmts[i]
                stream.write(fmt.format(row[i]))
                stream.write('\n' if i == 1 else ' ')
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
    d, cp = DMF(path), CPrint()
    if long_format:
        resources = list(d.find())
        uuid_pfx = _uuid_prefix([r.uuid for r in resources])
        fields = ('uuid', 'name', 'type', 'modified', 'created')
        widths = (uuid_pfx, 30, 20, 19, 19)
        colors = ('g', 'w', 'y', 'w', 'w')
        fmts = ['{{:{:d}s}}'.format(w) for w in widths]
        left_gutter = '| ' if relations else ''
        # table header
        cp.println(
            ' ' * len(left_gutter)
            + '  '.join(
                [
                    cp.colorize('@_h[{}]'.format(f).format(v))
                    for f, v in zip(fmts, fields)
                ]
            )
        )

        def datestr(t):
            p = pendulum.from_timestamp(t)
            return p.to_datetime_string()

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
            cp.println(
                left_gutter
                + '  '.join(
                    [
                        cp.colorize('@{}[{}]'.format(c, f.format(v)))
                        for c, f, v in zip(colors, fmts, values)
                    ]
                )
            )
            if relations and len(r.relations) > 0:
                relitems = []
                for rel in r.relations:
                    if rel.subject == r.uuid:
                        fmt = '@-w[{p}]->@b[{o}]'
                    else:
                        fmt = '@b[{s}]->@-w[{p}]'
                    item = fmt.format(
                        s=rel.subject[:uuid_pfx],
                        p=rel.predicate,
                        o=rel.object[:uuid_pfx],
                    )
                    relitems.append(item)
                cp.println('+-- {}'.format(' / '.join(relitems)))
    else:
        items = []
        for r in d.find():
            name_color = 'w'
            if r.name:
                name = r.name
            elif r.desc:
                name = r.desc[:40]
                name_color = 'b'
            else:
                name = r.uuid
                name_color = 'g'
            item = cp.colorize('@{}[{}]@y[:{}]'.format(name_color, name, r.type))
            items.append(item)
        if items:
            tsz = get_terminal_size((80, 20))
            term_width = max(tsz.columns, 1)
            columnized = _display_in_columns(items, max_line=term_width)
            print(columnized)


def _uuid_prefix(uuids, step=4, maxlen=32):
    """Get smallest multiple of `step` len prefix that gives unique values.
    """
    full = set(uuids)
    for n in range(step, maxlen, step):
        prefixes = {u[:n] for u in uuids}
        if len(prefixes) == len(full):
            return n
    return maxlen


def cat_resources(path, objects=(), color=True):
    d = DMF(path=path)
    cp = CPrint(color=color)
    unmatched = set(objects)
    first = True
    # get all resources,
    # display any that match an object as a prefix
    for r in d.find():
        for oid in unmatched:
            if r.uuid.startswith(oid):
                unmatched.remove(oid)  # don't show twice
                if not first:
                    _cat_resource_sep(cp)
                _cat_resource_show(cp, r)
                first = False
                break


def _cat_resource_sep(cp):
    cp.println('@b[{}]'.format('-' * 60))


def _cat_resource_show(cp, r):
    d = r.as_dict()
    json.dump(d, cp, indent=2)
    print()


# regular expression to find VT100 color escape codes
_noprint_re = re.compile(r'\033\[[0-9]+m')


def _display_in_columns(items, max_line=80, col_sep='  ', row_sep='\n'):
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
        return ''
    # Calculate item lengths, stripping terminal escapes
    lengths, nplengths = [], []
    for item in items:
        clean = _noprint_re.sub('', item)
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
            fmt = '{{:{n}s}}'.format(n=max_widths[col] + nplengths[i] - pad)
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
        raise ValueError('Invalid resource: {}'.format(err))
    return r


def _import_jupyternb(path):
    """Create & import a resource from a Jupyter Notebook file at `path`.
    Assume that `path` exists and is a Jupyter Notebook.

    Args:
        path (str): Jupyter Notebook file.
    Returns:
          (Resource) DMF Resource representing the notebook.
    """
    r = resource.Resource(type_=resource.TY_NOTEBOOK)
    filename = os.path.splitext(os.path.split(path)[1])[0]
    # XXX: add notebook 'metadata' as FilePath metadata attr
    r.v['datafiles'].append({'desc': filename, 'path': path})
    r.v['desc'] = filename
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
    r = resource.Resource(type_=resource.TY_CODE)
    filename = os.path.splitext(os.path.split(path)[1])[0]
    r.v['codes'].append({'name': filename, 'language': 'python', 'type': 'module'})
    r.v['datafiles'].append({'desc': filename, 'path': path})
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
    r = resource.Resource(type_=resource.TY_DATA)
    filename = os.path.split(path)[1]
    r.v['datafiles'].append({'desc': filename, 'path': path})
    r.v['desc'] = filename
    r.validate()
    return r
