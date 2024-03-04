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
import os

import pytest
from pyomo.environ import SolverFactory

from idaes.core.util.testing import _enable_scip_solver_for_testing


ampl_module_scip = pytest.importorskip(
    "ampl_module_scip", reason="'ampl_module_scip' not available"
)


@pytest.mark.unit
def test_path_manipulation():
    path_before_enabling = str(os.environ["PATH"])

    func_to_revert_changes = _enable_scip_solver_for_testing()
    path_after_enabling = str(os.environ["PATH"])
    sf = SolverFactory("scip")
    assert len(path_after_enabling) > len(path_before_enabling)
    assert str(ampl_module_scip.bin_dir) in str(sf.executable())
    assert sf.available()

    func_to_revert_changes()
    path_after_reverting = str(os.environ["PATH"])
    assert path_after_reverting == path_before_enabling
