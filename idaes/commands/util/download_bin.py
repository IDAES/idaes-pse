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
import os
import hashlib
from platform import machine
import idaes.logger as idaeslog
import tarfile
import idaes
from idaes.config import extra_binaries
from pyomo.common.download import FileDownloader
import urllib

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


def _get_release_platform(platform):
    # Check if platform (OS) maps to another platform
    platform = idaes.config.canonical_distro(platform)

    # Get machine type (e.g. x86_64, ...)
    mach = idaes.config.canonical_arch(machine())

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
        for i in range(1000):
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


def _create_download_package(platform, to_path, url, extra, extras_only, library_only):
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

    for e in extra:
        if e not in extra_binaries:
            _log.warning(f"Unknown extra package {e}, not installed.")
            continue
        if platform not in extra_binaries[e]:
            _log.warning(
                f"Extra package {e} not available for {platform}, not installed."
            )
            continue
        _add_pack(e)  # you have to explicitly ask for extras so assume you want
    if not extras_only:
        _add_pack("lib")
    if not library_only and not extras_only:
        _add_pack("solvers")

    return pname, ptar, ftar, furl


def _download_package(fd, name, frm, to, platform):
    _log.debug(f"Getting {name} from: {frm}")
    _log.debug(f"Saving solvers to: {to}")
    fd.set_destination_filename(to)

    try:
        fd.get_binary_file(frm)
    except urllib.error.HTTPError:
        raise Exception(f"{name} binaries are unavailable for {platform}")


def _verfiy_checksums(checksum, pname, ptar, ftar):
    # If release checksum is not False, nochecksum opt allows hash to be ignored
    if checksum:
        for n, p, f in zip(pname, ptar, ftar):
            hash_l = _hash(p)
            _log.debug(f"{n} Hash {hash_l}")
            if checksum.get(f, "") != hash_l:
                raise Exception(f"{n} hash does not match expected")


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
    if no_download:
        nochecksum = True

    # set the locations to download files to, to_path is an alternate
    # subdirectory of idaes.data_directory that can optionally be used to test
    # this function without interfering with anything else.  It's a subdirectory
    # of the data directory because that should be a safe place to store some
    # test files.
    if alt_path is not None:
        to_path = os.path.abspath(alt_path)
    if to_path is None:
        to_path = idaes.bin_directory
    else:
        to_path = os.path.join(idaes.data_directory, to_path)
    idaes._create_bin_dir(to_path)

    # Create FileDownloader object
    fd = FileDownloader(insecure=insecure, cacert=cacert)

    # Get platform information
    arch, platform = _get_arch_and_platform(fd, platform)
    platform = _get_release_platform(platform)

    # Get release url
    url = _get_release_url(release, url)

    # Get checksum
    checksum = _get_checksums(fd, to_path, release, nochecksum)

    # Set the binary file destinations
    pname, ptar, ftar, furl = _create_download_package(
        platform, to_path, url, extra, extras_only, library_only
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
    _verfiy_checksums(checksum, pname, ptar, ftar)

    # Extract solvers
    for n, p in zip(pname, ptar):
        _log.debug(f"Extracting files in {p} to {to_path}")
        with tarfile.open(p, "r") as f:
            f.extractall(to_path)
