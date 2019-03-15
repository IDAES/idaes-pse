##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for feed with flash.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.unit_models.feed_flash import FeedFlash
from idaes.property_models.ideal.BTX_ideal_VLE import (
                        BTXParameterBlock)
from idaes.ui.report import degrees_of_freedom


# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


# -----------------------------------------------------------------------------
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = BTXParameterBlock()

    m.fs.ff = FeedFlash(
            default={"property_package": m.fs.properties})

    assert hasattr(m.fs.ff, "flow_mol")
    assert hasattr(m.fs.ff, "mole_frac")
    assert hasattr(m.fs.ff, "temperature")
    assert hasattr(m.fs.ff, "pressure")

    assert hasattr(m.fs.ff, "outlet")
    assert len(m.fs.ff.outlet.vars) == 4
    assert hasattr(m.fs.ff.outlet, "flow_mol")
    assert hasattr(m.fs.ff.outlet, "mole_frac")
    assert hasattr(m.fs.ff.outlet, "temperature")
    assert hasattr(m.fs.ff.outlet, "pressure")

    assert hasattr(m.fs.ff, "isothermal")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = BTXParameterBlock()

    m.fs.ff = FeedFlash(
            default={"property_package": m.fs.properties})

    m.fs.ff.flow_mol.fix(1)
    m.fs.ff.temperature.fix(368)
    m.fs.ff.pressure.fix(101325)
    m.fs.ff.mole_frac[0, "benzene"].fix(0.5)
    m.fs.ff.mole_frac[0, "toluene"].fix(0.5)

    assert degrees_of_freedom(m) == 0

    m.fs.ff.initialize(outlvl=5,
                         optarg={'tol': 1e-6})

    assert (pytest.approx(101325.0, abs=1e3) ==
            m.fs.ff.outlet.pressure[0].value)
    assert (pytest.approx(368.00, abs=1e-0) ==
            m.fs.ff.outlet.temperature[0].value)
    assert (pytest.approx(1.0, abs=1e-2) ==
            m.fs.ff.outlet.flow_mol[0].value)
    assert (pytest.approx(0.3961, abs=1e-3) ==
            m.fs.ff.control_volume.
                properties_out[0].flow_mol_phase["Vap"].value)
