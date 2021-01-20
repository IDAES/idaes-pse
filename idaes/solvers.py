##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import os
import hashlib
import idaes.logger as idaeslog
import tarfile
import idaes
from shutil import copyfile
from pyomo.common.download import FileDownloader
import urllib
import idaes.config

_log = idaeslog.getLogger(__name__)


_release_base_url = idaes.config.release_base_url

def _hash(fname):
    with open(fname, 'rb') as f:
        fh = hashlib.sha256()
        while True:
            fb = f.read(10000)
            if len(fb) <= 0:
                break
            fh.update(fb)
    return str(fh.hexdigest())


def download_binaries(
    release=None,
    url=None,
    insecure=False,
    cacert=None,
    verbose=False,
    platform="auto",
    nochecksum=False,
    library_only=False,
    to_path=None):
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
    # set the locations to download files to, to_path is an alternate
    # subdirectory of idaes.data_directory that can optionally be used to test
    # this function without interfereing with anything else.  It's a subdirectory
    # of the data directory because that should be a safe place to store some
    # test files.
    if to_path is None:
        to_path = idaes.bin_directory
    else:
        to_path = os.path.join(idaes.data_directory, to_path)
    idaes._create_bin_dir(to_path)
    solvers_tar = os.path.join(to_path, "idaes-solvers.tar.gz")
    libs_tar = os.path.join(to_path, "idaes-lib.tar.gz")
    # Get a pyomo file downloader object and check Arch
    fd = FileDownloader(insecure=insecure, cacert=cacert)
    arch = fd.get_sysinfo()
    if arch[1] != 64:
        _log.error("IDAES Extensions currently only supports 64bit Python.")
        raise RuntimeError("IDAES Extensions currently only supports 64bit Python.")
    if platform == "auto":
        platform = arch[0]
        if platform == "linux":
            linux_dist = fd.get_os_version().replace(".", "")
            if linux_dist in idaes.config.known_binary_platform:
                platform = linux_dist
    if platform not in idaes.config.known_binary_platform:
        raise Exception("Unknow platform {}".format(platform))
    if platform in idaes.config.binary_platform_map:
        platform = idaes.config.binary_platform_map[platform]

    checksum = {} # storage for hashes if install from release
    # If a release is specified, use that to set the URL, also get hash info
    if release is not None:
        url = "/".join([_release_base_url, release])
        # if we're downloading an official release get hashes
        # check_to = place to store hash file
        check_to = os.path.join(to_path, f"sha256sum_{release}.txt")
        # check_from = release hashes url
        check_from = f"https://raw.githubusercontent.com/IDAES/idaes-ext/main/releases/sha256sum_{release}.txt"
        _log.debug("Getting release {}\n  checksum file {}".format(release, check_from))
        fd.set_destination_filename(check_to)
        fd.get_binary_file(check_from)
        # read the hashes file and store then in checksum dict
        with open(check_to, 'r') as f:
            for i in range(1000):
                line = f.readline(1000)
                if line == "":
                    break
                line = line.split(sep="  ")
                checksum[line[1].strip()] = line[0].strip()
    # Either URL should have been specified or obtained from specified release
    if url is not None:
        if not url.endswith("/"):
            c = "/"
        else:
            c = ""
        # Get specific package file URLs
        solvers_from = c.join([url, "idaes-solvers-{}-{}.tar.gz".format(platform, arch[1])])
        libs_from = c.join([url, "idaes-lib-{}-{}.tar.gz".format(platform, arch[1])])
        _log.debug("URLs \n  {}\n  {}\n  {}".format(url, solvers_from, libs_from))
        _log.debug("Destinations \n  {}\n  {}".format(solvers_tar, libs_tar))
        # Download solvers
        if not library_only:
            fd.set_destination_filename(solvers_tar)
            try:
                fd.get_binary_file(solvers_from)
            except urllib.error.HTTPError:
                 raise Exception(f"{platform} solver binaries are unavailable")
        # Download Libraries
        fd.set_destination_filename(libs_tar)
        try:
            fd.get_binary_file(libs_from)
        except urllib.error.HTTPError:
             raise Exception(f"{platform} library binaries are unavailable")
    else:
        raise Exception("Must provide a location to download binaries")
    # If release checksum is not empty, nochecksum opt allows hash to be ignored
    if checksum and not nochecksum:
        fn_s = "idaes-solvers-{}-{}.tar.gz".format(platform, arch[1])
        fn_l = "idaes-lib-{}-{}.tar.gz".format(platform, arch[1])
        # Check solvers package hash
        if not library_only:
            hash_s = _hash(solvers_tar)
            _log.debug("Solvers Hash {}".format(hash_s))
            if checksum.get(fn_s, "") != hash_s:
                raise Exception("Solver files hash does not match expected")
        # Check libraries package hash
        hash_l = _hash(libs_tar)
        _log.debug("Libs Hash {}".format(hash_l))
        if checksum.get(fn_l, "") != hash_l:
            raise Exception("Library files hash does not match expected")

    # Extract solvers
    if not library_only:
        _log.debug("Extracting files in {}".format(idaes.bin_directory))
        with tarfile.open(solvers_tar, 'r') as f:
            f.extractall(to_path)
    # Extract libraries
    _log.debug("Extracting files in {}".format(idaes.bin_directory))
    with tarfile.open(libs_tar, 'r') as f:
        f.extractall(to_path)
