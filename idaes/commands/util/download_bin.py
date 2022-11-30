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


def _get_platform(fd, platform, arch):
    if platform == "auto":
        platform = arch[0]
    if platform == "linux":
        platform = fd.get_os_version().replace(".", "")
        _log.debug(f"Detected Linux distribution: {platform}")
    # Check if platform (OS) maps to another platform
    platform = idaes.config.canonical_distro(platform)
    # Get machine type (e.g. x86_64, ...)
    mach = idaes.config.canonical_arch(machine())
    # full platform name
    platform = f"{platform}-{mach}"
    # See if machine is supported
    if platform not in idaes.config.base_platforms:
        raise UnsupportedPlatformError(f"Unsupported platfrom: {platform}.")
    _log.debug(f"Downloading binaries for {platform}")
    return platform


def _get_checksums(fd, to_path, release):
    checksum = {}  # storage for hashes if install from release
    check_to = os.path.join(to_path, f"sha256sum_{release}.txt")
    check_from = idaes.config.release_checksum_url.format(release)
    _log.debug(f"Getting checksum file {check_from}")
    fd.set_destination_filename(check_to)
    fd.get_binary_file(check_from)
    # read the hashes file and store then in checksum dict
    with open(check_to, "r") as f:
        for i in range(1000):
            line = f.readline(1000)
            if line == "":
                break
            line = line.split(sep="  ")
            checksum[line[1].strip()] = line[0].strip()
    return checksum


def _get_release_url_and_checksum(fd, to_path, release, url, nochecksum):
    checksum = False  # default if not checking checksums
    if release is not None:
        url = "/".join([_release_base_url, release])
        if not nochecksum:
            checksum = _get_checksums(fd, to_path, release)
        else:
            _log.debug("Skip release checksum verification at user request.")
    else:
        _log.debug("Release not specified.")
    if url is None:
        _log.debug("No release or URL was provided.")
        raise Exception("Must provide a location to download binaries")
    if url.endswith("/"):
        url = url[0:-1]  # if url ends with "/" remove it for proper join later
    _log.debug(f"Downloading binaries from {url}")
    return url, checksum


def _download_package(fd, name, frm, to, platform):
    _log.debug(f"Getting {name} from: {frm}")
    _log.debug(f"Saving solvers to: {to}")
    fd.set_destination_filename(to)
    try:
        fd.get_binary_file(frm)
    except urllib.error.HTTPError:
        raise Exception(f"{name} binaries are unavailable for {platform}")


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
    # this function without interfereing with anything else.  It's a subdirectory
    # of the data directory because that should be a safe place to store some
    # test files.
    if alt_path is not None:
        to_path = os.path.abspath(alt_path)
    if to_path is None:
        to_path = idaes.bin_directory
    else:
        to_path = os.path.join(idaes.data_directory, to_path)
    idaes._create_bin_dir(to_path)
    fd, arch = _get_file_downloader(insecure, cacert)
    platform = _get_platform(fd, platform, arch)
    url, checksum = _get_release_url_and_checksum(fd, to_path, release, url, nochecksum)
    # Set the binary file destinations
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

    # If release checksum is not empty, nochecksum opt allows hash to be ignored
    if checksum:
        for n, p, f in zip(pname, ptar, ftar):
            hash_l = _hash(p)
            _log.debug(f"{n} Hash {hash_l}")
            if checksum.get(f, "") != hash_l:
                raise Exception(f"{n} hash does not match expected")

    # Extract solvers
    for n, p in zip(pname, ptar):
        _log.debug(f"Extracting files in {p} to {to_path}")
        with tarfile.open(p, "r") as f:
            f.extractall(to_path)
