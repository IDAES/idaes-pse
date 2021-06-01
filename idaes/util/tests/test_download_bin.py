###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
#
###############################################################################
import idaes
import idaes.util.download_bin
import pytest
import os
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

def _del_data_file(path):
    try:
        os.remove(path)
    except OSError:
        pass

@pytest.mark.integration
def test_dl_bin():
    idaes._create_testing_dir()
    _del_data_file(os.path.join(idaes.testing_directory, "version_lib.txt"))
    _del_data_file(os.path.join(idaes.testing_directory, "version_solvers.txt"))
    ll = _log.getEffectiveLevel() # verbose will set level to DEBUG
    idaes.util.download_bin.download_binaries(
        release=idaes.config.default_binary_release,
        verbose=True,
        to_path="testing")
    _log.setLevel(ll) # set logger level bakc to whatever it was
    assert os.path.exists(os.path.join(idaes.testing_directory, "version_lib.txt"))
    assert os.path.exists(os.path.join(idaes.testing_directory, "version_solvers.txt"))

@pytest.mark.unit
def test_dl_bin_unknown():
    _del_data_file(os.path.join(idaes.testing_directory, "version_lib.txt"))
    _del_data_file(os.path.join(idaes.testing_directory, "version_solvers.txt"))
    with pytest.raises(Exception):
        idaes.util.download_bin.download_binaries(
            platform="unknown platform",
            release=idaes.config.default_binary_release,
            to_path="testing")
