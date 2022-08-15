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
Get project version and/or git hash
"""
# stdlib
from collections import namedtuple
import importlib
import json
import logging
from pathlib import Path
import pkg_resources
import subprocess
from types import ModuleType
from typing import Union

_log = logging.getLogger(__name__)


# Exceptions


class GetVersionError(Exception):
    pass


class ModuleImportError(GetVersionError):
    pass


class PipError(GetVersionError):
    pass


class GitHashError(GetVersionError):
    pass


# Classes and functions


#: Return from `Versioned.get_info()`
VersionInfo = namedtuple("VersionInfo", "version git_hash")


class Versioned:
    """Get version information for a module.

    Sample usage::

        from idaes.core.dmf.getver import Versioned
        import some_module

        vrs = Versioned(some_module)
        print(vrs.version)
        print(vrs.git_hash)

    """

    #: PIP command
    PIP = "pip"

    def __init__(self, module: Union[ModuleType, str], package: str = None):
        """Constructor.

        Args:
            module: Imported module or module namespace.
            package: Package name, overrides 'module' if present
        """
        if package:
            self._name = package
            self._mod = None
        else:
            if hasattr(module, "__package__"):
                mod = module
            else:
                try:
                    mod = importlib.import_module(module)
                except ModuleNotFoundError as err:
                    raise ModuleImportError(
                        f"Could not import module '{module}': {err}"
                    )
            self._mod = mod
            self._name = mod.__package__
        try:
            self._version = self._get_version()
        except PipError:
            raise GetVersionError(
                f"Could not get version for '{self._name}' from pip command '{self.PIP}'"
            )
        if self._version is None:
            raise GetVersionError(f"Could not get version for '{self._name}'")
        self._hash = None

    def get_info(self) -> VersionInfo:
        """Get hash and version without throwing any exceptions.
        If there was an exception, the hash may be None.

        Returns:
            The Git hash (which may be None) and the version (a string)
        """
        hash = None
        try:
            hash = self.git_hash
        except GetVersionError:
            pass
        return VersionInfo(git_hash=hash, version=self.version)

    @property
    def name(self) -> str:
        """Get the name of the module.

        Returns:
            Module name.
        """
        return self._name

    @property
    def version(self) -> str:
        """Get the version.

        Return:
            The version as a string, e.g. "2.1.0"
        """
        return self._version

    @property
    def git_hash(self) -> str:
        """Get the git hash.

        Returns:
            Git hash as a string of hexadecimal digits, e.g. "a5327b889b4606792cd31194f939c0d358e68e43"

        Raises:
            GitHashError
        """

        if self._hash is None:
            self._hash = self._get_git_hash()
        return self._hash

    # -- Protected methods -- #

    def _get_version(self):
        version = None
        try:
            version = pkg_resources.get_distribution(self._name).version
        except Exception as err:
            _log.warning(
                f"Could not get version from pkg_resources.get_distribution; trying pip ({self.PIP}). "
                f"Error message was: {err}"
            )
            status, output = subprocess.getstatusoutput(f"{self.PIP} list")
            if status != 0:
                raise PipError(f"Shell command '{self.PIP} list' failed")
            for s in output.split("\n"):
                if s.startswith(self._name):
                    version = s.split()[1]
                    break
            if version is None:
                raise PipError("Could not find version in pip output")
        return version

    def _get_git_hash(self):
        module_path = self._get_module_path()
        repo, git_hash = None, None
        # (1) try .git in module directory
        _log.debug(f"Looking for .git/ORIG_HEAD in module path '{module_path.parent}'")
        for git_head in module_path.parent.glob(r"**/.git/ORIG_HEAD"):
            _log.debug(f"Found ORIG_HEAD at {git_head}")
            repo = git_head.parent
            git_hash = git_head.open().read().strip()
            break
        if repo:
            _log.info(f"Get Git hash from repository rooted at: {repo.parent}")
        else:
            # (2) try installed dir
            _log.info(
                f"Get Git hash from installed module ({module_path}) dist-info directory"
            )
            # XXX: Find a less fragile way to do this (e.g. with setuptools)
            path = (
                module_path.parent
                / f"{self._name}-{self.version}.dist-info"
                / "direct_url.json"
            )
            if not path.exists():
                raise GitHashError(f"Cannot find magic file with Git hash: {path}")
            data = json.load(path.open())
            try:
                git_hash = data["vcs_info"]["commit_id"]
            except KeyError:
                raise GitHashError(f"Cannot find commit_id in metadata: {data}")
        _log.debug(f"Got git hash: {git_hash}")
        return git_hash

    def _get_module_path(self) -> Path:
        if self._mod is None:  # package name, use pkg_resources
            return Path(pkg_resources.get_distribution(self._name).module_path)
        paths = getattr(self._mod, "__path__", [])
        if paths:
            return Path(paths[0])
        loader = self._mod.__loader__
        if hasattr(loader, "path"):
            return Path(loader.path)
        raise ModuleImportError(f"Cannot get path for module: {self._mod}")


def get_version_info(
    module: Union[ModuleType, str] = None, package: str = None
) -> VersionInfo:
    """Get version info for a module.

    This is simply syntactic sugar around the `Versioned` object.

    Args:
        module: The same as accepted by the :class:`Versioned` constructor.
        package: If you have a package name instead of a module, use this instead (ignore 'module' if given)

    Returns:
        The same as :meth:`Versioned.get_info()` for the given module.
    """
    if package and module:
        _log.warning(
            "Since 'package' argument is given to get_version_info(), 'module' argument will be ignored"
        )
    vv = Versioned(module, package=package)
    return vv.get_info()


# If run from command-line, provide hash of a module
if __name__ == "__main__":
    import argparse

    logging.basicConfig()
    prs = argparse.ArgumentParser()
    prs.add_argument("module")
    args = prs.parse_args()
    _log.setLevel(logging.DEBUG)
    try:
        info = get_version_info(args.module)
    except ModuleImportError:
        info = get_version_info(package=args.module)
    print(f"{args.module}: version={info.version} hash={info.git_hash}")
