"""
Get project version and/or git hash
"""
# stdlib
from collections import namedtuple
import importlib
import json
import logging
from pathlib import Path
import subprocess
from types import ModuleType
from typing import Union

# ext
import git

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

        from idaes.dmf.getver import Versioned
        import some_module

        vrs = Versioned(some_module)
        print(vrs.version)
        print(vrs.git_hash)

    """

    #: PIP command
    PIP = "pip"

    def __init__(self, obj: Union[ModuleType, str]):
        """Constructor.

        Args:
            obj: Imported package or package name.
        """
        if hasattr(obj, "__package__"):
            mod = obj
        else:
            try:
                mod = importlib.import_module(obj)
            except ModuleNotFoundError as err:
                raise ModuleImportError(f"Could not import module '{obj}': {err}")
        self._mod = mod
        self._name = mod.__package__
        self._version = self._get_version()
        if self._version is None:
            raise PipError(f"Could not get version from pip for '{self._name}'")
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
        status, output = subprocess.getstatusoutput(f"{self.PIP} list")
        if status != 0:
            raise PipError(f"Shell command '{self.PIP} list' failed")
        version = None
        for s in output.split("\n"):
            if s.startswith(self._name):
                version = s.split()[1]
                break
        return version

    def _get_git_hash(self):
        module_path = self._get_module_path()
        repo, git_hash = None, None
        try:
            repo = git.Repo(module_path, search_parent_directories=True)
        except git.InvalidGitRepositoryError:
            pass
        if repo:
            _log.info(f"Get Git hash from repository rooted at: {repo.working_dir}")
            git_hash = repo.head.commit.hexsha
        else:
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
        paths = getattr(self._mod, "__path__", [])
        if paths:
            return Path(paths[0])
        loader = self._mod.__loader__
        if hasattr(loader, "path"):
            return Path(loader.path)
        raise ModuleImportError(f"Cannot get path for module: {self._mod}")


def get_version_info(module: Union[ModuleType, str]) -> VersionInfo:
    """Get version info for a module.

    This is simply syntactic sugar around the `Versioned` object.

    Args:
        module: The same as accepted by the :class:`Versioned` constructor.

    Returns:
        The same as :meth:`Versioned.get_info()` for the given module.
    """
    vv = Versioned(module)
    return vv.get_info()


# if __name__ == "__main__":
#     import sys
#
#     if len(sys.argv) != 2:
#         print(f"usage: getver.py <module>")
#         sys.exit(1)
#
#     module_name = sys.argv[1].strip()
#
#     hnd = logging.StreamHandler()
#     hnd.setFormatter(logging.Formatter("[%(levelname)s] %(asctime)s - %(message)s"))
#     _log.addHandler(hnd)
#     _log.setLevel(logging.DEBUG)
#
#     if False:
#         try:
#             vp = Versioned(module_name)
#         except RuntimeError as err:
#             _log.error(err)
#             print(f"Error getting version/hash: {err}")
#             sys.exit(2)
#
#         info = vp.get_info()
#     else:
#         info = get_version_info(module_name)
#
#     print(f"{module_name} version={info.version} hash={info.git_hash}")
#
#     sys.exit(0)
