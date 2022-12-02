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
"""The API in this module is mostly for internal use, e.g. from 'setup.py' to get the version of
the package. But :class:`Version` has been written to be usable as a general
versioning interface.

Example of using the class directly:

.. doctest::

    >>> from idaes.ver import Version
    >>> my_version = Version(1, 2, 3)
    >>> print(my_version)
    1.2.3
    >>> tuple(my_version)
    (1, 2, 3)
    >>> my_version = Version(1, 2, 3, 'alpha')
    >>> print(my_version)
    1.2.3.a
    >>> tuple(my_version)
    (1, 2, 3, 'alpha')
    >>> my_version = Version(1, 2, 3, 'candidate', 1)
    >>> print(my_version)
    1.2.3.rc1
    >>> tuple(my_version)
    (1, 2, 3, 'candidate', 1)

If you want to add a version to a class, e.g. a model, then
simply inherit from ``HasVersion`` and initialize it with the
same arguments you would give the :class:`Version` constructor:

.. doctest::

    >>> from idaes.ver import HasVersion
    >>> class MyClass(HasVersion):
    ...     def __init__(self):
    ...         super(MyClass, self).__init__(1, 2, 3, 'alpha')
    ...
    >>> obj = MyClass()
    >>> print(obj.version)
    1.2.3.a

"""
import os
import re
import sys

__author__ = "Dan Gunter"


class Version(object):
    """This class attempts to be compliant with a subset of
    `PEP 440 <https://www.python.org/dev/peps/pep-0440/>`_.

    Note: If you actually happen to read the PEP, you will notice
    that pre- and post- releases, as well as "release epochs", are not
    supported.
    """

    _specifiers = {
        "alpha": "a",
        "beta": "b",
        "candidate": "rc",
        "development": "dev",
        "final": "",  # this is the default
    }

    def __init__(
        self, major, minor, micro, releaselevel="final", serial=None, label=None
    ):
        """Create new version object.

        Provided arguments are stored in public class
        attributes by the same name.

        Args:
            major (int): Major version
            minor (int):  Minor version
            micro (int):  Micro (aka patchlevel) version
            releaselevel (str): Optional PEP 440 specifier
            serial (int): Optional number associated with releaselevel
            label (str): Optional local version label
        """
        if releaselevel not in self._specifiers:
            raise ValueError(
                'Value "{}" for releaselevel not in ({})'.format(
                    releaselevel, ",".join(sorted(self._specifiers.keys()))
                )
            )
        self.major, self.minor, self.micro = major, minor, micro
        self.releaselevel, self.serial, self.label = releaselevel, serial, label

    def __iter__(self):
        """Return version information as a sequence."""
        items = [self.major, self.minor, self.micro]
        if self.releaselevel != "final":
            items.append(self.releaselevel)
            if self.serial is not None:
                items.append(self.serial)
                if self.label is not None:
                    items.append(self.label)
            elif self.label is not None:
                items.append(0)  # placeholder for serial
                items.append(self.label)
        for it in items:
            yield it

    def __str__(self):
        """Return version information as a string."""
        return "{}.{}.{}{}".format(
            self.major,
            self.minor,
            self.micro,
            (
                ""
                if self.releaselevel == "final"
                else "."
                + self._specifiers[self.releaselevel]
                + ("" if self.serial is None else str(self.serial))
                + ("" if self.label is None else "+" + self.label)
            ),
        )


class HasVersion(object):
    """Interface for a versioned class."""

    def __init__(self, *args):
        """Constructor creates a `version` attribute that is
        an instance of :class:`Version` initialized with the provided args.

        Args:
            *args: Arguments to be passed to Version constructor.
        """
        self.version = Version(*args)


def git_hash():
    """Get current git hash, with no dependencies on external packages."""
    # find git root (in grandparent dir to this file, if anywhere)
    git_root = os.path.realpath(os.path.join(__file__, "..", "..", ".git"))
    if not os.path.exists(git_root) or not os.path.isdir(git_root):
        raise ValueError(f"git root '{git_root}' not found")
    # get HEAD ref's file
    try:
        head = open(os.path.join(git_root, "HEAD"))
    except FileNotFoundError as err:
        raise ValueError(f"cannot open HEAD: {err}")
    # parse file looking for 'ref: <path>'
    head_ref = None
    for line in head:
        ref_match = re.match(r"ref:\s+(\S+)", line)
        if ref_match:
            head_ref = ref_match.group(1)
            break
    if head_ref is None:
        raise ValueError(f"no ref found in HEAD '{head}'")
    # read value of ref in <path> found previously
    ref_file = os.path.join(git_root, head_ref)
    try:
        ref = open(ref_file).read().strip()
    except FileNotFoundError:
        raise ValueError(f"ref file '{ref_file}' not found")
    return ref


# Get git hash. No output unless IDAES_DEBUG is set in env
gh = None
try:
    try:
        gh = git_hash()
        if os.environ.get("IDAES_DEBUG", None):
            print(f"git hash = {gh}", file=sys.stderr)
    except ValueError as err:
        if os.environ.get("IDAES_DEBUG", None):
            print(f"git_hash() error: {err}", file=sys.stderr)
except NameError:  # eg, if invoked from setup.py
    pass

#: Package's version as an object
package_version = Version(2, 0, 0, "beta", 2, gh)

#: Package's version as a simple string
__version__ = str(package_version)
