#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
# TODO: Missing doc strings
# pylint: disable=missing-module-docstring
# pylint: disable=missing-class-docstring

# TODO: protected access issues
# pylint: disable=protected-access

import os
import hashlib
from platform import machine
import tarfile
import urllib

from pyomo.common.download import FileDownloader

import idaes
import idaes.logger as idaeslog
from idaes.config import release_major, base_platforms

_log = idaeslog.getLogger(__name__)
_release_base_url = idaes.config.release_base_url


class UnsupportedPlatformError(RuntimeError):
    pass


def hash_file_sha256(fname):
    """Calculate sha256 hash of potentially large files.

    Args:
        fname: path to file to hash

    Returns:
        hash as a string
    """
    with open(fname, "rb") as f:
        fh = hashlib.sha256()
        while True:
            fb = f.read(10000)
            if len(fb) <= 0:
                break
            fh.update(fb)
    return str(fh.hexdigest())


_hash = hash_file_sha256


def _get_file_downloader(insecure, cacert):
    fd = FileDownloader(insecure=insecure, cacert=cacert)
    arch = fd.get_sysinfo()
    return fd, arch


def _get_arch_and_platform(fd, platform):
    arch = fd.get_sysinfo()

    if platform == "auto":
        platform = arch[0]
    if platform == "linux":
        platform = fd.get_os_version().replace(".", "")
        _log.debug(f"Detected Linux distribution: {platform}")

    return arch, platform


def _map_legacy_distro(platform, release):
    """
    Map canonical distro names to the closest supported pre-4.x build target.
    """
    if release_major(release) >= 4:
        return platform

    legacy_distro_map = {
        "ubuntu2404": "ubuntu2204",
        "el9": "ubuntu2204",
    }
    return legacy_distro_map.get(platform, platform)


def _get_release_platform(platform, release):
    # Check if platform (OS) maps to another platform
    platform = idaes.config.canonical_distro(platform)

    # Apply release-specific fallback for older binary releases
    platform = _map_legacy_distro(platform, release)

    # Get machine type (e.g. x86_64, ...)
    mach = idaes.config.canonical_arch(machine(), release=release, platform=platform)

    # full platform name
    platform = f"{platform}-{mach}"
    # See if machine is supported
    if platform not in idaes.config.base_platforms:
        raise UnsupportedPlatformError(f"Unsupported platform: {platform}.")
    _log.debug(f"Downloading binaries for {platform}")

    return platform


def _get_release_url(release, url):
    if release is not None:
        url = "/".join([_release_base_url, release])
    else:
        _log.debug("Release not specified.")

    if url is None:
        _log.debug("No release or URL was provided.")
        # PYLINT-TODO
        # pylint: disable-next=broad-exception-raised
        raise Exception("Must provide a location to download binaries")

    if url.endswith("/"):
        url = url[0:-1]  # if url ends with "/" remove it for proper join later

    _log.debug(f"Downloading binaries from {url}")
    return url


def _get_checksum_paths(to_path, release):
    check_to = os.path.join(to_path, f"sha256sum_{release}.txt")
    check_from = idaes.config.release_checksum_url.format(release)
    _log.debug(f"Getting checksum file {check_from}")

    return check_to, check_from


def _download_checksum(fd, check_to, check_from):
    fd.set_destination_filename(check_to)
    fd.get_binary_file(check_from)


def _read_checksum_file(check_to):
    # read the hashes file and store them in checksum dict
    checksum = {}  # storage for hashes if install from release

    with open(check_to, "r") as f:
        for i in range(1000):  # pylint: disable=unused-variable
            line = f.readline(1000)
            if line == "":
                break
            line = line.split(sep="  ")
            checksum[line[1].strip()] = line[0].strip()

    return checksum


def _get_checksums(fd, to_path, release, nochecksum):
    if nochecksum:
        _log.debug("Skip release checksum verification at user request.")
        return False  # default if not checking checksums

    check_to, check_from = _get_checksum_paths(to_path, release)
    _download_checksum(fd, check_to, check_from)

    return _read_checksum_file(check_to)


def _supported_extra_binaries(release):
    major = release_major(release)
    if major >= 4:
        return {}
    return {
        "petsc": base_platforms,
    }


def _normalize_extras(extra, release, platform=None):
    """
    Filter requested extras to those supported for the given release,
    optionally also filtering by platform.
    This is mostly to support the 4.x series.
    """
    extra = tuple(extra)

    if release_major(release) >= 4 and "petsc" in extra:
        _log.info(
            "Ignoring extra 'petsc' for IDAES %s: PETSc is included in the base package.",
            release,
        )

    supported = _supported_extra_binaries(release)
    kept = []

    for e in extra:
        if e not in supported:
            _log.debug("Ignoring unknown or unsupported extra binary package '%s'", e)
            continue
        if platform is not None and platform not in supported[e]:
            _log.debug(
                "Ignoring extra binary package '%s' for unsupported platform '%s'",
                e,
                platform,
            )
            continue
        kept.append(e)

    return tuple(kept)


def _create_download_package(
    release, platform, to_path, url, extra, extras_only, library_only
):
    major = release_major(release)
    extra = _normalize_extras(extra, release)

    pname = []
    ptar = []
    ftar = []
    furl = []

    def _add_pack(name):
        f = f"idaes-{name}-{platform}.tar.gz"
        ftar.append(f)
        ptar.append(os.path.join(to_path, f))
        pname.append(name)
        furl.append("/".join([url, f]))

    if major >= 4:
        if not extras_only:
            if not library_only:
                _add_pack("solvers")
            _add_pack("functions")
    else:
        if not extras_only:
            if not library_only:
                _add_pack("solvers")
            _add_pack("lib")
        for e in extra:
            _add_pack(e)

    return pname, ptar, ftar, furl


def _download_package(fd, name, frm, to, platform):
    _log.debug(f"Getting {name} from: {frm}")
    _log.debug(f"Saving solvers to: {to}")
    fd.set_destination_filename(to)

    try:
        fd.get_binary_file(frm)
    except urllib.error.HTTPError as e:
        # PYLINT-TODO
        # pylint: disable-next=broad-exception-raised
        raise Exception(f"{name} binaries are unavailable for {platform}. \nError: {e}")


def _verify_checksums(checksum, pname, ptar, ftar):
    # If release checksum is not False, nochecksum opt allows hash to be ignored
    if checksum:
        for n, p, f in zip(pname, ptar, ftar):
            hash_l = _hash(p)
            _log.debug(f"{n} Hash {hash_l}")
            if checksum.get(f, "") != hash_l:
                # PYLINT-TODO
                # pylint: disable-next=broad-exception-raised
                raise Exception(f"{n} hash does not match expected")


def _splitpath(path):
    """Split path into elements

    This routine is similar to str.split(), except when the pieces are
    re-joined by os.path.join(), the result will still be an absolute
    path (if ``path`` was an absolute path). Note that the path
    elements are returned in reverse order (which is convenient for
    _resolve_path below)

    """
    ans = []
    while 1:
        a, b = os.path.split(path)
        if b:
            ans.append(b)
        if a == path:
            if a:
                ans.append(a)
            return ans
        path = a


def _resolve_path(path, links):
    """Resolve path to the actual filesystem location

    This routine is similar to os.realpath(), but with the added feature
    of recognizing not only hard / symbolic links on the filesystem, but
    also links in the ``links`` dict that have not yet been committed to
    the filesystem.

    """
    dirs = _splitpath(path)
    path = ""
    target = None
    while dirs:
        path = os.path.realpath(os.path.join(path, dirs.pop()))
        if path in links:
            target = path
        while target in links:
            if links[target][0] == tarfile.SYMTYPE:
                if os.path.isabs(links[target][1]):
                    path = ""
                else:
                    path = os.path.dirname(path)
                dirs.extend(_splitpath(links[target][1]))
                target = None
            else:  # links[target][0] == tarfile.LNKTYPE:
                target_dir, target_name = os.path.split(links[target][1])
                target = os.path.join(_resolve_path(target_dir, links), target_name)
        if target is not None:
            path = target
            target = None
    return path


def _verify_tar_member_targets(tar, to_path, links=None):
    """Check tar for unsafe members (references outside to_path)

    Nominally, this would be a simple task: verify that every member of
    the tarfile (when converted to an absolute path and normalized)
    starts with the to_path. Unfortunately, symbolic links make this
    more difficult. It is not sufficient to just check the absolute
    path of the link target, as sequences of links could walk out of the
    target space. In this routine, we attempt to identify symbolic
    links and resolve the actual location where each tar member would be
    stored on the filesystem.

    Example::

        [DIRTYPE] dir1
        [SYMTYPE] dir1/dir2 -> ..
        [SYMTYPE] dir1/dir2/dir3 -> ..
        [REGTYPE] dir1/dir2/dir3/file  # <-- this is ../file

    Parameters
    ----------
    tar : TarFile
        tarfile to verify

    to_path : str
        target directory to extract tarfile into

    links : dict[str, tuple[bytes, str]]
        dictionary that records "virtual" symbolic links; that is,
        symbolic links in the tarfile that have not yet been extracted,
        but will exist during / after the extraction process. This maps
        the actual symbolic link location to a (type, target) tuple.
        The target is a resolved (absolute) path for hard links and a
        potentially relative path for soft links.

    """
    to_path = os.path.realpath(to_path)

    if links is None:
        # A dict mapping {from: (type, to)}
        links = {}

    for member in tar.getmembers():
        member_path = _resolve_path(os.path.join(to_path, member.name), links)
        if not member_path.startswith(to_path):
            # PYLINT-TODO
            # pylint: disable-next=broad-exception-raised
            raise Exception(
                f"Tarball {tar.name} contained potentially unsafe member "
                f"{member.name} that would extract to {member_path} outside "
                f"target directory {to_path}"
            )
        if member.issym():
            links[member_path] = (member.type, member.linkname)
        elif member.islnk() and member.name != member.linkname:
            links[member_path] = (member.type, os.path.join(to_path, member.linkname))


def _warn_if_libraries_exist(release):
    """
    Warn when installing 4.x+ extensions over an environment that appears to
    contain an older flat-layout binary installation.

    In older releases, shared libraries lived directly in idaes.bin_directory.
    In 4.x+, libraries move to idaes.data_directory/lib while executables remain
    in idaes.data_directory/bin. If old files are still present in the old bin
    directory, users may encounter overwritten binaries or runtime linking issues.
    """
    if release_major(release) < 4:
        return

    old_lib = os.path.join(idaes.bin_directory, "libpynumero_ASL.dylib")
    if os.path.exists(old_lib):
        # Setting this as error so it shows up by default. It's really a warning
        # but I don't think it's worth it to change the logging level of the
        # whole system just to print this message as a warning to then immediately
        # change it back.
        _log.error(
            "Existing pre-4.x extension files were detected in %s. "
            "IDAES 4.x+ extensions use a different install layout and may "
            "overwrite existing binaries or conflict with older libraries. "
            "Consider removing old extension files from %s before installing.",
            idaes.bin_directory,
            idaes.bin_directory,
        )


def download_binaries(
    release=None,
    url=None,
    insecure=False,
    cacert=None,
    verbose=False,
    platform="auto",
    nochecksum=False,
    library_only=False,
    no_download=False,
    extras_only=False,
    extra=(),
    to_path=None,
    alt_path=None,
):
    """
    Download IDAES solvers and libraries and put them in the right location. Need
    to supply either local or url argument.

    Args:
        url (str): a url to download binary files to install files
        to_path: for testing only, an alternate subdirectory of the idaes data
                 directory for a test binary install.
    Returns:
        None
    """
    if verbose:
        _log.setLevel(idaeslog.DEBUG)
    if release is None:
        release = idaes.config.default_binary_release
    if no_download:
        nochecksum = True
    else:
        _warn_if_libraries_exist(release)

    # set the locations to download files to, to_path is an alternate
    # subdirectory of idaes.data_directory that can optionally be used to test
    # this function without interfering with anything else. It's a subdirectory
    # of the data directory because that should be a safe place to store some
    # test files.
    if alt_path is not None:
        to_path = os.path.abspath(alt_path)
    if to_path is None:
        if release_major(release) >= 4:
            to_path = idaes.data_directory
        else:
            to_path = idaes.bin_directory
    else:
        to_path = os.path.join(idaes.data_directory, to_path)
    idaes._create_bin_dir(to_path)

    # Create FileDownloader object
    fd = FileDownloader(insecure=insecure, cacert=cacert)

    # Get platform information
    arch, platform = _get_arch_and_platform(  # pylint: disable=unused-variable
        fd, platform
    )
    platform = _get_release_platform(platform, release=release)

    # Get release url
    url = _get_release_url(release, url)

    # Get checksum
    checksum = _get_checksums(fd, to_path, release, nochecksum)

    # Set the binary file destinations
    pname, ptar, ftar, furl = _create_download_package(
        release, platform, to_path, url, extra, extras_only, library_only
    )

    if no_download:
        d = {
            "release": release,
            "platform": platform,
        }
        for n, p, u in zip(pname, ptar, furl):
            d[u] = p
        return d

    # download packages
    for n, p, u in zip(pname, ptar, furl):
        _download_package(fd, n, frm=u, to=p, platform=platform)

    # Verify checksums
    _verify_checksums(checksum, pname, ptar, ftar)

    # Extract solvers
    links = {}
    for n, p in zip(pname, ptar):
        _log.debug(f"Extracting files in {p} to {to_path}")
        with tarfile.open(p, "r") as f:
            _verify_tar_member_targets(f, to_path, links)
            f.extractall(to_path)
