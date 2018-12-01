##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Utility functions.
"""
# stdlib
import importlib
import logging
import os
import re
import shutil
import signal
import sys
import tempfile
import time
# third party
import psutil

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'

_log = logging.getLogger(__name__)


def strlist(x, sep=', '):
    # type: (list, str) -> str
    return sep.join([str(item) for item in x])


# def pluck(obj, key):
#     """Remove and return obj[key].
#     """
#     value = obj[key]
#     del obj[key]
#     return value


def get_file(file_or_path, mode='r'):
    """Open a file for reading, or simply return the file object.
    """
    if hasattr(file_or_path, 'read'):
        return file_or_path
    return open(file_or_path, mode)


def import_module(name):
    mod = importlib.import_module(name)
    return mod


def get_module_version(mod):
    """Find and return the module version.

    Version must look like a semantic version with
    `<a>.<b>.<c>` parts; there can be arbitrary extra
    stuff after the `<c>`. For example::

        1.0.12
        0.3.6
        1.2.3-alpha-rel0

    Args:
        mod (module): Python module
    Returns:
        (str) Version string or None if not found
    Raises:
        ValueError if version is found but not valid
    """
    v = getattr(mod, '__version__', None)
    if v is None:
        return None
    pat = r'\d+\.\d+\.\d+.*'
    if not re.match(pat, v):
        raise ValueError('Version "{}" does not match regular expression '
                         'pattern "{}"'.format(v, pat))
    return v


def get_module_author(mod):
    """Find and return the module author.

    Args:
        mod (module): Python module
    Returns:
        (str) Author string or None if not found
    Raises:
        nothing
    """
    return getattr(mod, '__author__', None)


class TempDir(object):
    """Simple context manager for mkdtemp().
    """
    def __init__(self, *args):
        self._d = None
        self._a = args

    def __enter__(self):
        self._d = tempfile.mkdtemp(*self._a)
        return self._d

    def __exit__(self, *args):
        if self._d:
            shutil.rmtree(self._d)
            self._d = None


def is_jupyter_notebook(filename):
    # type: (str) -> bool
    """See if this is a Jupyter notebook.
    """
    if not filename.endswith('.ipynb'):
        return False
    return True  # XXX: look inside?


def is_python(filename):
    # type: (str) -> bool
    """See if this is a Python file.
    Do *not* import the source code.
    """
    if not filename.endswith('.py'):
        return False
    return True  # XXX: look inside?


def is_resource_json(filename):
    import json

    if not filename.endswith('.json'):
        return False
    # if it's under 1MB, try parsing it:
    st = os.stat(filename)
    if st.st_size < 1e6:
        try:
            d = json.load(open(filename))
        except (UnicodeDecodeError, json.JSONDecodeError):
            return False
        # look for a couple distinctive keys
        if '_id' in d and 'sources' in d:
            return True
    else:
        return False  # screw it


def find_process_byname(name, uid=None):
    # type: (str, uid: int) -> Generator[int, None, None]
    """Generate zero or more PIDs where 'name' is part of either the
    first or second token in the command line.
    Optionally also filter the returned PIDs to only those with
    a 'real' user UID (UID) equal to the provided uid. If None,
    the default, is given, then use the current process UID. Providing
    a value of < 0 will skip the filter.
    """
    try:
        proc_entries = os.listdir('/proc')
    except OSError:
        _log.error('Attempt to find process "{}" in /proc, but '
                   'cannot open /proc for reading'.format(name))
        return
    if uid is None:
        uid = os.getuid()
    for pid in proc_entries:
        try:
            pid = int(pid)
        except ValueError:
            continue
        found = None
        try:
            with open('/proc/{}/cmdline'.format(pid), mode='rb') as fd:
                content = fd.read().decode().split('\x00')
            if len(content) > 0:
                if content[0].endswith(name):
                    found = pid
                elif len(content) > 1 and content[1].endswith(name):
                    found = pid
        except (OSError, IOError):
            continue
        if found:
            if uid is None:
                yield found
            else:
                found_uid = psutil.Process(pid).uids().real
                if found_uid == uid:
                    yield found
    return


def terminate_pid(pid, waitfor=1):
    result = True
    os.kill(pid, signal.SIGTERM)
    time.sleep(waitfor)
    if psutil.pid_exists(pid):
        os.kill(pid, signal.SIGKILL)
        time.sleep(waitfor)
        if psutil.pid_exists(pid):
            result = False
    return result


def datetime_timestamp(v):
    """Get numeric timestamp.
    This will work under both Python 2 and 3.
    """
    if hasattr(v, 'timestamp'):  # Python 2/3 test
        # Python 2
        result = v.timestamp()
    else:
        # Python 3
        result = time.mktime(v.timetuple()) + v.microsecond / 1e6
    return result


class CPrint(object):
    """Colorized terminal printing.

    Codes are below. To use:

        cprint = CPrint()
        cprint('This has no colors')  # just like print()
        cprint('This is @b[blue] and @_r[red underlined]')

    You can use the same class as a no-op by just passing `color=False` to
    the constructor.
    """
    COLORS = {'h': '\033[1m\033[95m',
              'r': '\033[91m',
              'g': '\033[92m',
              'y': '\033[93m',
              'b': '\033[94m',
              'm': '\033[95m',
              'c': '\033[96m',
              'w': '\033[97m',
              '.': '\033[0m',
              '*': '\033[1m',
              '-': '\033[2m',
              '_': '\033[4m'}

    _styled = re.compile(r'@([*_-]?[hbgyrwcm*_-])\[([^]]*)\]')

    def __init__(self, color=True):
        self._c = color

    def println(self, s):
        print(self.colorize(s))

    def __call__(self, *args):
        return self.println(args[0])

    def write(self, s):
        sys.stdout.write(self.colorize(s))

    def colorize(self, s):
        chunks = []
        last = 0
        c, stop = '', ''
        for m in re.finditer(self._styled, s):
            code, text = m.groups()
            clen = len(code)
            if self._c:
                if clen == 2:
                    if code[0] == code[1]:
                        c = self.COLORS[code]
                    else:
                        c = self.COLORS[code[0]] + self.COLORS[code[1]]
                else:
                    c = self.COLORS[code]
                stop = self.COLORS['.']
            x, y = m.span()
            chunks.append(s[last:x])  # text since last found piece
            chunks.append(c + s[x + 2 + clen:y - 1] + stop)  # colorized
            last = y
        chunks.append(s[last:])  # to end of string
        return ''.join(chunks)


def mkdir_p(path, *args):
    """Try to create all non-existent components of a path.

    Args:
        path (str): Path to create
        args: Other arguments for `os.mkdir()`.
    Returns:
        None
    Raises:
        os.error: Raised from `os.mkdir()`
    """
    plist = path.split(os.path.sep)
    if plist[0] == '':
        # do not try to create filesystem root
        dir_name = os.path.sep
        plist = plist[1:]
    else:
        dir_name = ''
    for p in plist:
        dir_name = os.path.join(dir_name, p)
        if not os.path.exists(dir_name):
            os.mkdir(dir_name, *args)
