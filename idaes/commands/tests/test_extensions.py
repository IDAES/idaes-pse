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
"""
Tests for idaes.commands.extensions
"""

import os
import tempfile

from click.testing import CliRunner
import pytest

import idaes
from idaes.commands import extensions

test_release_3 = "3.4.2"
test_release_4 = "4.0.1"


@pytest.fixture
def runner():
    return CliRunner()


##################
# print helpers  #
##################
@pytest.mark.unit
def test_print_header_and_footer_use_echo_callback():
    lines = []
    extensions.print_header("My Title", echo=lines.append)
    extensions.print_footer(echo=lines.append)

    assert any("IDAES Extensions My Title" in line for line in lines)


###########################
# _get_extensions_root    #
###########################
@pytest.mark.unit
def test_get_extensions_root_explicit_bin_directory():
    with tempfile.TemporaryDirectory() as tmpdir:
        assert extensions._get_extensions_root(bin_directory=tmpdir) == tmpdir


@pytest.mark.unit
def test_get_extensions_root_default_is_a_directory():
    root = extensions._get_extensions_root()
    assert root in (idaes.data_directory, idaes.bin_directory)


###########################
# _find_extensions_file   #
###########################
@pytest.mark.unit
def test_find_extensions_file_in_root():
    with tempfile.TemporaryDirectory() as tmpdir:
        target = os.path.join(tmpdir, "VERSION_SOLVERS.md")
        with open(target, "w") as f:
            f.write("solver v2\n")

        found = extensions._find_extensions_file(
            "VERSION_SOLVERS.md", bin_directory=tmpdir
        )
        assert found == target


@pytest.mark.unit
def test_find_extensions_file_in_bin_subdir():
    with tempfile.TemporaryDirectory() as tmpdir:
        bin_subdir = os.path.join(tmpdir, "bin")
        os.mkdir(bin_subdir)
        target = os.path.join(bin_subdir, "version_solvers.txt")
        with open(target, "w") as f:
            f.write("solver v1\n")

        found = extensions._find_extensions_file(
            "version_solvers.txt", bin_directory=tmpdir
        )
        assert found == target


@pytest.mark.unit
def test_find_extensions_file_missing_returns_none():
    with tempfile.TemporaryDirectory() as tmpdir:
        found = extensions._find_extensions_file("fake_file.md", bin_directory=tmpdir)
        assert found is None


##############################
# print_extensions_version   #
##############################
@pytest.mark.unit
def test_print_extensions_version_new_layout(capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(os.path.join(tmpdir, "VERSION_SOLVERS.md"), "w") as f:
            f.write(f"solvers-{test_release_4}\n")
        with open(os.path.join(tmpdir, "VERSION_FUNCTIONS.md"), "w") as f:
            f.write(f"functions-{test_release_4}\n")

        rc = extensions.print_extensions_version(bin_directory=tmpdir)
        out = capsys.readouterr().out

    assert rc == 0
    assert "solvers" in out
    assert "functions" in out


@pytest.mark.unit
def test_print_extensions_version_old_layout(capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(os.path.join(tmpdir, "version_solvers.txt"), "w") as f:
            f.write(f"solvers-{test_release_3}\n")
        with open(os.path.join(tmpdir, "version_lib.txt"), "w") as f:
            f.write(f"lib-{test_release_3}\n")

        rc = extensions.print_extensions_version(bin_directory=tmpdir)
        out = capsys.readouterr().out

    assert rc == 0
    assert "solvers" in out
    assert "lib" in out


@pytest.mark.unit
def test_print_extensions_version_library_only_skips_solvers(capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(os.path.join(tmpdir, "VERSION_FUNCTIONS.md"), "w") as f:
            f.write("functions-4.0.1\n")

        rc = extensions.print_extensions_version(
            library_only=True, bin_directory=tmpdir
        )
        out = capsys.readouterr().out

    assert rc == 0
    assert "functions-4.0.1" in out
    assert "Solvers:" not in out


@pytest.mark.unit
def test_print_extensions_version_missing_files(capsys):
    with tempfile.TemporaryDirectory() as tmpdir:
        rc = extensions.print_extensions_version(bin_directory=tmpdir)
        out = capsys.readouterr().out

    assert rc == 0
    assert "no version file found" in out


######################
# print_build_info   #
######################
@pytest.mark.unit
@pytest.mark.parametrize("release", [test_release_3, test_release_4])
def test_print_build_info(capsys, release):
    extensions.print_build_info(release)
    out = capsys.readouterr().out

    assert "Build Information" in out
    assert "Platform:" in out
    assert "Architecture:" in out


########################
# get-extensions CLI   #
########################
@pytest.mark.unit
def test_get_extensions_url_and_release_conflict(runner):
    result = runner.invoke(
        extensions.get_extensions, ["--url", "http://example", "--release", "4.0.1"]
    )
    assert result.exit_code == 0
    assert "either a release or url not both" in result.output


@pytest.mark.unit
@pytest.mark.parametrize("release", [test_release_3, test_release_4])
def test_get_extensions_info(runner, release):
    result = runner.invoke(extensions.get_extensions, ["--info", "--release", release])
    assert result.exit_code == 0
    assert "Build Information" in result.output
