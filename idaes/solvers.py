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


_log = idaeslog.getLogger(__name__)


_release_base_url = "https://github.com/IDAES/idaes-ext/releases/download"

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
    platform="auto"):
    """
    Download IDAES solvers and libraries and put them in the right location. Need
    to supply either local or url argument.

    Args:
        url (str): a url to download binary files to install files

    Returns:
        None
    """
    if verbose:
        _log.setLevel(idaeslog.DEBUG)
    idaes._create_bin_dir()
    solvers_tar = os.path.join(idaes.bin_directory, "idaes-solvers.tar.gz")
    libs_tar = os.path.join(idaes.bin_directory, "idaes-lib.tar.gz")
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

    checksum = {}
    if release is not None:
        # if a release is specified it takes precedence over a url
        url = "/".join([_release_base_url, release])
        # if we're downloading an official release check checksum
        check_to = os.path.join(idaes.bin_directory, f"sha256sum_{release}.txt")
        check_from = f"https://raw.githubusercontent.com/IDAES/idaes-ext/main/releases/sha256sum_{release}.txt"
        _log.debug("Getting release {}\n  checksum file {}".format(release, check_from))
        fd.set_destination_filename(check_to)
        fd.get_binary_file(check_from)
        with open(check_to, 'r') as f:
            for i in range(1000):
                line = f.readline(1000)
                if line == "":
                    break
                line = line.split(sep="  ")
                checksum[line[1].strip()] = line[0].strip()
    if url is not None:
        if not url.endswith("/"):
            c = "/"
        else:
            c = ""
        solvers_from = c.join([url, "idaes-solvers-{}-{}.tar.gz".format(platform, arch[1])])
        libs_from = c.join([url, "idaes-lib-{}-{}.tar.gz".format(platform, arch[1])])
        _log.debug("URLs \n  {}\n  {}\n  {}".format(url, solvers_from, libs_from))
        _log.debug("Destinations \n  {}\n  {}".format(solvers_tar, libs_tar))
        if platform == 'darwin':
            raise Exception('Mac OSX currently unsupported')
        fd.set_destination_filename(solvers_tar)
        fd.get_binary_file(solvers_from)
        fd.set_destination_filename(libs_tar)
        fd.get_binary_file(libs_from)
    else:
        raise Exception("Must provide a location to download binaries")

    if checksum:
        # if you are downloading a release and not a specific URL verify checksum
        fn_s = "idaes-solvers-{}-{}.tar.gz".format(platform, arch[1])
        fn_l = "idaes-lib-{}-{}.tar.gz".format(platform, arch[1])
        hash_s = _hash(solvers_tar)
        hash_l = _hash(libs_tar)
        _log.debug("Solvers Hash {}".format(hash_s))
        _log.debug("Libs Hash {}".format(hash_l))
        if checksum.get(fn_s, "") != hash_s:
            raise Exception("Solver files hash does not match expected")
        if checksum.get(fn_l, "") != hash_l:
            raise Exception("Library files hash does not match expected")

    _log.debug("Extracting files in {}".format(idaes.bin_directory))
    with tarfile.open(solvers_tar, 'r') as f:
        f.extractall(idaes.bin_directory)
    _log.debug("Extracting files in {}".format(idaes.bin_directory))
    with tarfile.open(libs_tar, 'r') as f:
        f.extractall(idaes.bin_directory)
