#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
import pytest
import base64
import io
import os
import tarfile
import tempfile
import shutil
from platform import machine

from pyomo.common.download import FileDownloader

import idaes
import idaes.commands.util.download_bin as dlb
import idaes.logger as idaeslog
import idaes.config as idaes_config


_log = idaeslog.getLogger(__name__)


def _del_data_file(path):
    try:
        os.remove(path)
    except OSError:
        pass


@pytest.mark.unit
def test_get_file_downloader():
    fd, arch = dlb._get_file_downloader(insecure=True, cacert=None)

    assert isinstance(fd, FileDownloader)
    assert fd.insecure
    assert fd.cacert is None


@pytest.mark.unit
def test_get_file_downloader_cacert_error():
    with pytest.raises(
        RuntimeError, match="cacert='foo' does not refer to a valid file."
    ):
        dlb._get_file_downloader(insecure=False, cacert="foo")


@pytest.mark.unit
def test_get_arch_and_platform_auto():
    fd = FileDownloader()

    arch, platform = dlb._get_arch_and_platform(fd, platform="auto")

    assert arch == fd.get_sysinfo()

    if arch[0] == "linux":
        assert platform == fd.get_os_version().replace(".", "")
    else:
        assert platform == arch[0]


@pytest.mark.unit
def test_get_release_platform_mapping():
    mach = idaes.config.canonical_arch(machine())
    for p, m in idaes_config.binary_distro_map.items():
        if mach == "aarch64" and m == "el7":
            continue
        if mach == "x86_64" and m == "darwin":
            continue
        output = dlb._get_release_platform(p)
        assert output.startswith(m)
        assert output in idaes_config.base_platforms


@pytest.mark.unit
def test_get_release_url_release():
    output = dlb._get_release_url("foo/", None)

    assert output == "https://github.com/IDAES/idaes-ext/releases/download/foo"


@pytest.mark.unit
def test_get_release_url_release_and_url():
    output = dlb._get_release_url("foo/", "bar/")

    assert output == "https://github.com/IDAES/idaes-ext/releases/download/foo"


@pytest.mark.unit
def test_get_release_url():
    output = dlb._get_release_url(None, "bar/")

    assert output == "bar"


@pytest.mark.unit
def test_get_release_url_none():
    with pytest.raises(Exception, match="Must provide a location to download binaries"):
        dlb._get_release_url(None, None)


@pytest.mark.unit
def test_get_checksum_paths():
    cto, cfrom = dlb._get_checksum_paths("foo", "bar")

    assert cto == os.path.join("foo", f"sha256sum_bar.txt")
    assert (
        cfrom
        == "https://raw.githubusercontent.com/IDAES/idaes-ext/main/releases/sha256sum_bar.txt"
    )


@pytest.mark.unit
def test_download_checksum():
    # Mock paths to avoid downloads
    fd = FileDownloader()
    fd.retrieve_url = lambda url: bytes("\n", encoding="utf-8")

    tmpdir = tempfile.mkdtemp()
    target = os.path.join(tmpdir, "bin.txt")

    output = dlb._download_checksum(fd, target, None)

    assert output is None

    assert os.path.getsize(target) == 1

    shutil.rmtree(tmpdir)


@pytest.mark.unit
def test_read_checksum():
    path = os.path.dirname(__file__)
    testfile = os.path.join(path, "checksum.txt")

    checksum = dlb._read_checksum_file(testfile)

    assert checksum == {
        "idaes-lib-centos6-64.tar.gz": "449457bc91fac41e8c443b4d9cf2ba9181201e84318e890833f9685c01dcd531",
        "idaes-solvers-centos6-64.tar.gz": "2efda2ce61f9fc1d7180b5b141077ce52c2682f6532e9c8d873164af13162948",
        "idaes-lib-centos7-64.tar.gz": "dfae8a6e1bb4afe27220062d8b47abd780207f7fd85cf27126cf2e2c41ef6d36",
        "idaes-solvers-centos7-64.tar.gz": "ae465b064b7e74c91127558c6f0ebbae1b01756e5c4d59d5681a70088a70ffcd",
        "idaes-lib-centos8-64.tar.gz": "5301a517e098262be2e4536207114c855eb9c0cf17aad085a4964ba64674e457",
        "idaes-solvers-centos8-64.tar.gz": "a597b72d71b12e019a0d629e5d8eb089382d8051488f6a51f597511fc2872b3f",
        "idaes-lib-ubuntu1804-64.tar.gz": "bc6bbf8e2ede3e5afc80fed8656eddb575a2464cebb593afdab9171c2b5723e6",
        "idaes-solvers-ubuntu1804-64.tar.gz": "25299b15e270a2e7224e5aa4d6665bec074ceb727d7a4cf87b11d01fc96acdef",
        "idaes-lib-ubuntu1910-64.tar.gz": "6fbc3034d8c94763c8db9798d6349eb31a1980af92829cef5fecdfca3bc0454a",
        "idaes-solvers-ubuntu1910-64.tar.gz": "f635ceccb71c0e5c6c0b8bcb26b507d5531ec34630aa897d794d0e7d3dc0d401",
        "idaes-lib-ubuntu2004-64.tar.gz": "272ee8fc338074e056fed2c6115cc504526c389e5880c8819af7e6c7b9f1fc61",
        "idaes-solvers-ubuntu2004-64.tar.gz": "586fe646436e27f3e3e0f4c7a3c983c500f5048f2e65eb28f9303d790d600986",
        "idaes-lib-windows-64.tar.gz": "3e17985dfc12c7dea8da47b6d15c45d2d27e68accaacc521d708ef31c08d001e",
        "idaes-solvers-windows-64.tar.gz": "c9c141a38b208025bd16a3f3af8d6deb2dcd03777b63c04c5f8842f2743b7e7b",
    }


@pytest.mark.parametrize("platform", idaes_config.base_platforms)
@pytest.mark.unit
def test_create_download_package(platform):
    pname, ptar, ftar, furl = dlb._create_download_package(
        platform, "foo", "bar", [], False, False
    )

    assert len(ptar) == 2
    assert (os.path.join("foo", f"idaes-lib-{platform}.tar.gz")) in ptar
    assert (os.path.join("foo", f"idaes-solvers-{platform}.tar.gz")) in ptar

    assert ftar == [f"idaes-lib-{platform}.tar.gz", f"idaes-solvers-{platform}.tar.gz"]
    assert furl == [
        f"bar/idaes-lib-{platform}.tar.gz",
        f"bar/idaes-solvers-{platform}.tar.gz",
    ]


@pytest.mark.parametrize("platform", idaes_config.base_platforms)
@pytest.mark.unit
def test_create_download_package_extras(platform):
    pname, ptar, ftar, furl = dlb._create_download_package(
        platform, "foo", "bar", ["petsc", "baz"], False, False
    )

    # baz should be ignored as an unknown extra
    assert len(ptar) == 3
    assert (os.path.join("foo", f"idaes-lib-{platform}.tar.gz")) in ptar
    assert (os.path.join("foo", f"idaes-solvers-{platform}.tar.gz")) in ptar
    assert (os.path.join("foo", f"idaes-petsc-{platform}.tar.gz")) in ptar

    assert ftar == [
        f"idaes-petsc-{platform}.tar.gz",
        f"idaes-lib-{platform}.tar.gz",
        f"idaes-solvers-{platform}.tar.gz",
    ]
    assert furl == [
        f"bar/idaes-petsc-{platform}.tar.gz",
        f"bar/idaes-lib-{platform}.tar.gz",
        f"bar/idaes-solvers-{platform}.tar.gz",
    ]


@pytest.mark.parametrize("platform", idaes_config.base_platforms)
@pytest.mark.unit
def test_create_download_package_extras_only(platform):
    pname, ptar, ftar, furl = dlb._create_download_package(
        platform, "foo", "bar", ["petsc", "baz"], True, False
    )

    # baz should be ignored as an unknown extra
    assert len(ptar) == 1
    assert (os.path.join("foo", f"idaes-petsc-{platform}.tar.gz")) in ptar

    assert ftar == [f"idaes-petsc-{platform}.tar.gz"]
    assert furl == [f"bar/idaes-petsc-{platform}.tar.gz"]


@pytest.mark.parametrize("platform", idaes_config.base_platforms)
@pytest.mark.unit
def test_create_download_package_library_only(platform):
    pname, ptar, ftar, furl = dlb._create_download_package(
        platform, "foo", "bar", [], False, True
    )

    assert len(ptar) == 1
    assert (os.path.join("foo", f"idaes-lib-{platform}.tar.gz")) in ptar

    assert ftar == [f"idaes-lib-{platform}.tar.gz"]
    assert furl == [f"bar/idaes-lib-{platform}.tar.gz"]


@pytest.mark.unit
def test_download_package():
    # Mock paths to avoid downloads
    fd = FileDownloader()
    fd.retrieve_url = lambda url: bytes("\n", encoding="utf-8")

    tmpdir = tempfile.mkdtemp()
    target = os.path.join(tmpdir, "bin.txt")

    output = dlb._download_package(fd, "foo", None, target, "bar")

    assert output is None

    assert os.path.getsize(target) == 1

    shutil.rmtree(tmpdir)


@pytest.mark.unit
def test_dl_bin_unknown():
    _del_data_file(os.path.join(idaes.testing_directory, "version_lib.txt"))
    _del_data_file(os.path.join(idaes.testing_directory, "version_solvers.txt"))
    with pytest.raises(Exception):
        dlb.download_binaries(
            platform="unknown platform",
            release=idaes.config.default_binary_release,
            to_path="testing",
        )


@pytest.mark.unit
def test_verify_tarfile_members():
    # Tarball with:
    #
    # [DIRTYPE] dir1
    # [SYMTYPE] dir1/dir2 -> ..
    # [LNKTYPE] dir1/dir2 -> dir1/dir2
    # [SYMTYPE] dir1/dir2/dir3 -> ..
    # [REGTYPE] dir1/dir2/dir3/file
    #
    softlink_tar_gz = (
        b"H4sICLbv92MAA3NvZnRsaW5rLnRhcgDt0+0KgjAUgOFdyq5Ad7a5XU9gwSIItO4/l0EFR"
        b"QlO+nifHwoeYQdfbFMntSrLDGJslBEfw+39Sol3sclPox/mEo1Vuim819mxP6w6rdW27V"
        b"Pq9k/fezX/Um3uP1xswTPG/nFCfxskKG2rquBWF/T/wP5ejCgtSyxH//ET54src8aU/sG"
        b"4Ye58Y/j/l3Dfv96k3Xr2M3LgEPx7/UXy3FqntJl9kwf+vD8AAAAAAAAAAAAA4DecAEnw"
        b"Tw0AKAAA"
    )
    data = io.BytesIO(base64.b64decode(softlink_tar_gz))
    tar = tarfile.open(None, "r:gz", data)
    # verify tarfile contents:
    # for f in tar.getmembers():
    #     print(f.type, f.name, f.linkname)
    with pytest.raises(
        Exception,
        match="Tarball None contained potentially unsafe member dir1/dir2/dir3/file",
    ):
        dlb._verify_tar_member_targets(tar, os.path.abspath(""))

    # Tarball with:
    #
    # [DIRTYPE] dir1
    # [SYMTYPE] dir1/dir2 -> ..
    # [LNKTYPE] dir1/dir2 -> dir1/dir2
    # [LNKTYPE] dir1/dir2/dir3 -> dir1/dir2
    # [REGTYPE] dir3/file
    #
    hardlink_tar_gz = (
        b"H4sICHg++GMAA2hhcmRsaW5rLnRhcgDt1FkKwjAUQNG3lKzA5qUZ1iNUISIIre7fxor6I"
        b"w6Q4nDPRwpNIY9caJd7baQuO0opiFWf4u3zStS3KZS3yY/7mqwTEyrPdXIY9sveGNl0Q8"
        b"797u53j/a/VFf6j4ureMbUP73Q30WNYtxiUXGqM/p/YH+vVsXoHMPRf7risrR1znijf/K"
        b"O/nMo1Zt13q4qnlECx+if6h9Vy/8/OC/GVpzp4s/7AwAAAAAAAAAAAAB+wxHqgNJYACgA"
        b"AA=="
    )
    data = io.BytesIO(base64.b64decode(hardlink_tar_gz))
    tar = tarfile.open(None, "r:gz", data)
    # verify tarfile contents:
    # for f in tar.getmembers():
    #     print(f.type, f.name, f.linkname)
    with pytest.raises(
        Exception, match="Tarball None contained potentially unsafe member dir3/file"
    ):
        dlb._verify_tar_member_targets(tar, os.path.abspath(""))
