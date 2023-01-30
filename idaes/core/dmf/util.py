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
Utility functions.
"""
# stdlib
from datetime import datetime, timedelta, timezone
import importlib
import json
from json import JSONDecodeError
import logging
from math import floor, log
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
from typing import Union
import yaml

# third-party
import colorama

__author__ = "Dan Gunter"

_log = logging.getLogger(__name__)


def yaml_load(*args):
    return yaml.load(*args, Loader=yaml.SafeLoader)


def strlist(x, sep=", "):
    # type: (list, str) -> str
    return sep.join([str(item) for item in x])


def get_file(file_or_path, mode="r"):
    """Open a file for reading, or simply return the file object."""
    if hasattr(file_or_path, "read"):
        return file_or_path
    return open(file_or_path, mode)


def as_path(
    f: Union[Path, str, None],
    must_exist: bool = False,
    must_be_dir: bool = False,
    must_be_file: bool = False,
) -> Path:
    """Simply coerce input to a `Path`.

    Args:
        f: Path, string path, or None
        must_exist: Path must exist

    Returns:
        Input if it's a Path or None, else the Path(input).

    Raises:
        ValueError: if one of the 'must' conditions fails
    """
    if f is not None:
        f = Path(f)
        must_exist = must_exist or (must_be_dir or must_be_file)
        if must_exist and not f.exists():
            raise ValueError(f"Path '{f}' must exist")
        if must_be_file and not f.is_file():
            raise ValueError(f"Path '{f}' must be a file")
        if must_be_dir and not f.is_dir():
            raise ValueError(f"Path '{f}' must be a directory")
    return f


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
        (str): Version string or None if not found
    Raises:
        ValueError if version is found but not valid
    """
    v = getattr(mod, "__version__", None)
    if v is None:
        return None
    pat = r"\d+\.\d+\.\d+.*"
    if not re.match(pat, v):
        raise ValueError(
            'Version "{}" does not match regular expression '
            'pattern "{}"'.format(v, pat)
        )
    return v


def get_module_author(mod):
    """Find and return the module author.

    Args:
        mod (module): Python module
    Returns:
        (str): Author string or None if not found
    Raises:
        nothing
    """
    return getattr(mod, "__author__", None)


class TempDir(object):
    """Simple context manager for mkdtemp()."""

    def __init__(self, *args):
        self._d = None
        self._a = args

    def __enter__(self):
        self._d = mkdtemp(*self._a)
        return self._d

    def __exit__(self, *args):
        if self._d:
            shutil.rmtree(self._d)
            self._d = None


def is_jupyter_notebook(filename, check_contents=True):
    # type: (str, bool) -> bool
    """See if this is a Jupyter notebook.

    Args:
        filename (str): Full path to file
        check_contents (bool): Check contents of filename in addition to
            if filename extension is `.ipynb`

    Returns:
        (bool): If filename is a jupyter notebook
    """
    if not filename.endswith(".ipynb"):
        return False
    if check_contents:
        try:
            nb = json.load(open(filename))
        except (UnicodeDecodeError, JSONDecodeError):
            return False
        for key in "cells", "metadata", "nbformat":
            if key not in nb:
                return False
    return True


def is_python(filename):
    # type: (str) -> bool
    """See if this is a Python file.
    Do *not* import the source code.
    """
    if not filename.endswith(".py"):
        return False
    return True  # XXX: look inside?


def is_resource_json(filename, max_bytes=1e6):
    """Is this file a JSON Resource?

    Args:
        filename (str): Full path to file
        max_bytes (int): Max. allowable size. Since we try to parse
             the file, this saves potential DoS issues. Large files
             are a bad idea anyways, since this is metadata and may
             be stored somewhere with a record size limit (like MongoDB).

    Returns:
        (bool): Whether it's a resource JSON file.
    """
    if not filename.endswith(".json"):
        return False
    # get size
    st = os.stat(filename)
    # if it's under max_bytes, parse it
    if st.st_size <= max_bytes:
        try:
            d = json.load(open(filename))
        except (UnicodeDecodeError, JSONDecodeError):
            return False
        # look for a couple distinctive keys
        for key in "id_", "type":
            if key not in d:
                return False
        return True
    else:
        # if it's over max_bytes, it's "bad"
        return False


def datetime_timestamp(v):
    """Get numeric timestamp.
    This will work under both Python 2 and 3.

    Args:
        v (datetime.datetime): Date/time value

    Returns:
        (float): Floating point timestamp
    """
    if hasattr(v, "timestamp"):  # Python 2/3 test
        # Python 2
        result = v.timestamp()
    else:
        # Python 3
        result = time.mktime(v.timetuple()) + v.microsecond / 1e6
    return result


def parse_datetime(datetime_string: str) -> datetime:
    """Wrapper around datetime parser to allow for easy migration to other
    (better) implementations.

    Args:
        datetime_string: String in expected format

    Returns:
        A `datetime.datetime` object

    Raises:
        ValueError: Cannot parse this string
        TypeError: argument must be a string
    """
    if not isinstance(datetime_string, str):
        raise TypeError("parse_datetime: argument must be str")
    if hasattr(datetime, "fromisoformat"):
        # new in 3.7, and we support back to 3.6
        return datetime.fromisoformat(datetime_string)
    else:
        # The following code is from the datetime module in Python 3.8
        dtstr = datetime_string[0:10]
        if len(dtstr) != 10:
            raise ValueError(f"Invalid isoformat string: '{dtstr}'")
        year = int(dtstr[0:4])
        if dtstr[4] != "-":
            raise ValueError(f"Invalid date separator: '{dtstr[4]}'")
        if len(dtstr) < 7:
            raise ValueError(f"Invalid month: '{dtstr[5:]}'")
        month = int(dtstr[5:7])
        if dtstr[7] != "-":
            raise ValueError("Invalid date separator")
        day = int(dtstr[8:10])
        date_components = [year, month, day]

        tstr = datetime_string[11:]
        if tstr:
            try:
                time_components = _parse_isoformat_time(tstr)
            except ValueError:
                raise ValueError(f"Invalid isoformat string: {datetime_string!r}")
        else:
            time_components = [0, 0, 0, 0, None]

        return datetime(*(date_components + time_components))


#################
#
# The following two functions are from 3.8 Python standard library
# They are not present in 3.6, but are present in 3.7+
#


def _parse_isoformat_time(tstr):
    # Format supported is HH[:MM[:SS[.fff[fff]]]][+HH:MM[:SS[.ffffff]]]
    len_str = len(tstr)
    if len_str < 2:
        raise ValueError("Isoformat time too short")

    # This is equivalent to re.search('[+-]', tstr), but faster
    tz_pos = tstr.find("-") + 1 or tstr.find("+") + 1
    timestr = tstr[: tz_pos - 1] if tz_pos > 0 else tstr

    time_comps = _parse_hh_mm_ss_ff(timestr)

    tzi = None
    if tz_pos > 0:
        tzstr = tstr[tz_pos:]

        # Valid time zone strings are:
        # HH:MM               len: 5
        # HH:MM:SS            len: 8
        # HH:MM:SS.ffffff     len: 15

        if len(tzstr) not in (5, 8, 15):
            raise ValueError("Malformed time zone string")

        tz_comps = _parse_hh_mm_ss_ff(tzstr)
        if all(x == 0 for x in tz_comps):
            tzi = timezone.utc
        else:
            tzsign = -1 if tstr[tz_pos - 1] == "-" else 1

            td = timedelta(
                hours=tz_comps[0],
                minutes=tz_comps[1],
                seconds=tz_comps[2],
                microseconds=tz_comps[3],
            )

            tzi = timezone(tzsign * td)

    time_comps.append(tzi)

    return time_comps


def _parse_hh_mm_ss_ff(tstr):
    # Parses things of the form HH[:MM[:SS[.fff[fff]]]]
    len_str = len(tstr)

    time_comps = [0, 0, 0, 0]
    pos = 0
    for comp in range(0, 3):
        if (len_str - pos) < 2:
            raise ValueError("Incomplete time component")

        time_comps[comp] = int(tstr[pos : pos + 2])

        pos += 2
        next_char = tstr[pos : pos + 1]

        if not next_char or comp >= 2:
            break

        if next_char != ":":
            raise ValueError("Invalid time separator: %c" % next_char)

        pos += 1

    if pos < len_str:
        if tstr[pos] != ".":
            raise ValueError("Invalid microsecond component")
        else:
            pos += 1

            len_remainder = len_str - pos
            if len_remainder not in (3, 6):
                raise ValueError("Invalid microsecond component")

            time_comps[3] = int(tstr[pos:])
            if len_remainder == 3:
                time_comps[3] *= 1000

    return time_comps


#
################


class ColorTerm:
    """For colorized printing, a very simple wrapper that
    allows colorama objects, or nothing, to be used.
    """

    class EmptyStr:
        """Return an empty string on any attribute requested."""

        def __getattr__(self, a):
            return ""

    def __init__(self, enabled=True):
        self._width = None
        if enabled:
            colorama.init(autoreset=True)
            # Colorama colors and styles
            F = self.Fore = colorama.Fore
            B = self.Back = colorama.Back
            S = self.Style = colorama.Style
            # Aliases for colors and styles
            (
                self.blue,
                self.green,
                self.yellow,
                self.red,
                self.cyan,
                self.magenta,
                self.bluebg,
                self.greenbg,
                self.yellowbg,
                self.redbg,
                self.cyanbg,
                self.magentabg,
                self.white,
                self.black,
                self.whitebg,
                self.blackbg,
                self.bold,
                self.dim,
                self.normal,
                self.reset,
                self.resetc,
            ) = (
                F.BLUE,
                F.GREEN,
                F.YELLOW,
                F.RED,
                F.CYAN,
                F.MAGENTA,
                B.BLUE,
                B.GREEN,
                B.YELLOW,
                B.RED,
                B.CYAN,
                B.MAGENTA,
                F.WHITE,
                F.BLACK,
                B.WHITE,
                B.BLACK,
                S.BRIGHT,
                S.DIM,
                S.NORMAL,
                S.RESET_ALL,
                F.RESET,
            )
        else:
            # Make any attribute (e.g. `Fore.BLUE`) return an empty
            # string, thus disabling all the color codes.
            self.Fore, self.Back, self.Style = (
                self.EmptyStr(),
                self.EmptyStr(),
                self.EmptyStr(),
            )
            (
                self.blue,
                self.green,
                self.yellow,
                self.red,
                self.cyan,
                self.magenta,
                self.bluebg,
                self.greenbg,
                self.yellowbg,
                self.redbg,
                self.cyanbg,
                self.magentabg,
                self.white,
                self.black,
                self.whitebg,
                self.blackbg,
                self.bold,
                self.dim,
                self.normal,
                self.reset,
                self.resetc,
            ) = [""] * 21

    @property
    def width(self):
        if self._width is None:
            self._width = shutil.get_terminal_size()[0]
        return self._width


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
    if plist[0] == "":
        # do not try to create filesystem root
        dir_name = os.path.sep
        plist = plist[1:]
    elif plist[0].endswith(":"):  # windows root
        dir_name = plist[0] + os.path.sep
        plist = plist[1:]
    else:
        dir_name = ""
    for p in plist:
        dir_name = os.path.join(dir_name, p)
        if not os.path.exists(dir_name):
            _log.debug(f"mkdir_p: create directory '{dir_name}' args={args}")
            os.mkdir(dir_name, *args)
        else:
            _log.debug(f"mkdir_p: directory '{dir_name}' exists, do not create")


def uuid_prefix_len(uuids, step=4, maxlen=32):
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


_pow_abbr = "KMGTPE"


def size_prefix(number, base2=False):
    """Return an abbreviation for the size, up to E (exa)"""
    fudge = 1e-12
    if number == 0:
        return "0"
    elif number < 0:
        number = -number
        negative = "-"
    else:
        negative = ""
    if base2:
        nlog = int(floor(log(number, 2) / 10 + fudge))
        print(nlog)
        npow = number / (2 ** (10 * nlog))
        extra = "i"
    else:
        nlog = int(floor(log(number, 10) / 3 + fudge))
        npow = number / (10 ** (3 * nlog))
        extra = ""
    abbr = _pow_abbr[int(floor(nlog))]
    if abbr == "-":
        return f"{negative}{number}"
    if floor(npow) == (floor(npow * 10) / 10):  # don't print ".0"
        num_str = f"{negative}{npow:.0f}{abbr}{extra}"
    else:
        num_str = f"{negative}{npow:.1f}{abbr}{extra}"
    return num_str


# Edge-case temporary dir handling

# Cache this, since the result of gettempdir() is also cached
_is_curdir = None


def fail_if_tempdir_is_curdir():
    global _is_curdir
    if _is_curdir is None:
        td = tempfile.gettempdir()  # initializes it 1st time
        try:
            curdir = os.getcwd()
        except (AttributeError, OSError):
            curdir = os.curdir
        _is_curdir = td == curdir
    if _is_curdir:
        raise RuntimeError(
            "abort: temporary directory is going to be in "
            "the current directory. Please set one of the "
            "following environment variables to a suitable "
            "directory: TMPDIR, TEMP, or TMP"
        )


def mkdtemp(*args, **kwargs):
    fail_if_tempdir_is_curdir()
    return tempfile.mkdtemp(*args, **kwargs)


def NamedTemporaryFile(*args, **kwargs):
    fail_if_tempdir_is_curdir()
    return tempfile.NamedTemporaryFile(*args, **kwargs)
