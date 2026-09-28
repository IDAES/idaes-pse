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

import os
import tempfile
import pytest
import idaes
import logging


class TestIdaesConfigure(object):
    @pytest.mark.unit
    def test_write_config(self):
        idaes.write_config(
            os.path.join(idaes.testing_directory, "config_testing_orig.json")
        )
        idaes.write_config(
            os.path.join(idaes.testing_directory, "config_testing_default.json")
        )

    @pytest.mark.unit
    def test_read_default_config(self):
        # set some stuff to not the default
        idaes.cfg.use_idaes_solvers = False
        assert idaes.cfg.use_idaes_solvers == False
        # read the default
        idaes.read_config(
            os.path.join(idaes.testing_directory, "config_testing_default.json")
        )
        assert idaes.cfg.ipopt.options.nlp_scaling_method == "gradient-based"
        assert idaes.cfg.use_idaes_solvers == True

    @pytest.mark.unit
    def test_change_env(self):
        with idaes.temporary_config_ctx():
            # default config uses idaes solvers
            idaes.cfg.use_idaes_solvers = False
            idaes.reconfig()
            pf = os.environ["PATH"]
        pt = os.environ["PATH"]
        assert pf != pt

    @pytest.mark.unit
    def test_setup_environment_old_layout(self):
        old_path = os.environ.get("PATH", "")
        old_ld = os.environ.get("LD_LIBRARY_PATH", "")
        old_dyld = os.environ.get("DYLD_LIBRARY_PATH", "")

        with tempfile.TemporaryDirectory() as tmpdir:
            bin_dir = os.path.join(tmpdir, "bin")
            os.mkdir(bin_dir)

            idaes.config.setup_environment(bin_dir, use_idaes_solvers=True)

            assert bin_dir in os.environ["PATH"]
            if os.name != "nt":
                assert bin_dir in os.environ["LD_LIBRARY_PATH"]
                assert bin_dir in os.environ["DYLD_LIBRARY_PATH"]

        os.environ["PATH"] = old_path
        if os.name != "nt":
            os.environ["LD_LIBRARY_PATH"] = old_ld
            os.environ["DYLD_LIBRARY_PATH"] = old_dyld

    @pytest.mark.unit
    @pytest.mark.skipif(os.name == "nt", reason="non-Windows only")
    def test_setup_environment_new_layout_unix(self):
        old_path = os.environ.get("PATH", "")
        old_ld = os.environ.get("LD_LIBRARY_PATH", "")
        old_dyld = os.environ.get("DYLD_LIBRARY_PATH", "")

        with tempfile.TemporaryDirectory() as tmpdir:
            bin_dir = os.path.join(tmpdir, "bin")
            lib_dir = os.path.join(tmpdir, "lib")
            os.mkdir(bin_dir)
            os.mkdir(lib_dir)

            idaes.config.setup_environment(bin_dir, use_idaes_solvers=True)

            assert bin_dir in os.environ["PATH"]
            assert lib_dir in os.environ["LD_LIBRARY_PATH"]
            assert lib_dir in os.environ["DYLD_LIBRARY_PATH"]

        os.environ["PATH"] = old_path
        os.environ["LD_LIBRARY_PATH"] = old_ld
        os.environ["DYLD_LIBRARY_PATH"] = old_dyld

    @pytest.mark.unit
    @pytest.mark.skipif(os.name != "nt", reason="Windows only")
    def test_setup_environment_new_layout_windows(self):
        old_path = os.environ.get("PATH", "")

        with tempfile.TemporaryDirectory() as tmpdir:
            bin_dir = os.path.join(tmpdir, "bin")
            lib_dir = os.path.join(tmpdir, "lib")
            os.mkdir(bin_dir)
            os.mkdir(lib_dir)

            idaes.config.setup_environment(bin_dir, use_idaes_solvers=True)

            path_parts = os.environ["PATH"].split(os.pathsep)
            assert bin_dir in path_parts
            assert lib_dir in path_parts

        os.environ["PATH"] = old_path


@pytest.mark.unit
def test_canonical_arch():
    assert idaes.config.canonical_arch("Intel64") == "x86_64"
    assert idaes.config.canonical_arch("AMD64") == "x86_64"
    assert idaes.config.canonical_arch("ARM64") == "aarch64"
    # If we provide a release, it should kick in some extra logic
    assert (
        idaes.config.canonical_arch("ARM64", release="4.0.1", platform="darwin")
        == "arm64"
    )


@pytest.mark.unit
def test_canonical_distro():
    assert idaes.config.canonical_distro("kubUntu1804") == "ubuntu1804"
    assert idaes.config.canonical_distro("Ubuntu1804") == "ubuntu1804"
    assert idaes.config.canonical_distro("Windows") == "windows"
    assert idaes.config.canonical_distro("daRwin") == "darwin"


@pytest.mark.unit
def test_warning_to_except():
    with idaes.temporary_config_ctx():
        with pytest.raises(RuntimeError):
            _log = logging.getLogger("idaes")
            idaes.cfg.warning_to_exception = True
            idaes.reconfig()
            _log.warning("Hey! Don't do that.")


@pytest.mark.unit
def test_deprecate_to_except():
    with idaes.temporary_config_ctx():
        with pytest.raises(RuntimeError):
            _log = logging.getLogger("idaes")
            idaes.cfg.deprecation_to_exception = True
            idaes.reconfig()
            _log.warning("DEPRECATED: Hey! Don't use that.")


@pytest.mark.unit
def test_get_data_directory():
    # this tests the function that give that data, binary and testing directories
    # it doesn't actually change the idaes configuration.
    if "IDAES_DATA" in os.environ:
        odat = os.environ["IDAES_DATA"]
    else:
        odat = None
    # Test set by env variable on a directory that should exist
    os.environ["IDAES_DATA"] = idaes.config.data_directory
    dd, bd, td = idaes.config.get_data_directory()
    assert dd == idaes.config.data_directory
    # Test set by env variable on a directory that should exist
    os.environ["IDAES_DATA"] = "/this_directory_doesnt_exist/for_real/not_here"
    dd, bd, td = idaes.config.get_data_directory()
    assert dd is None
    if odat is not None:
        os.environ["IDAES_DATA"] = odat
    else:
        del os.environ["IDAES_DATA"]
