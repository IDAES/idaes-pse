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
Tests for StateJunction unit model.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.unit_models.statejunction import StateJunction
from idaes.property_models.examples.saponification_thermo import (
    SaponificationParameterBlock)
from idaes.core.util.model_statistics import degrees_of_freedom


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
@pytest.mark.build
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()

    m.fs.sj = StateJunction(default={"property_package": m.fs.properties})

    assert hasattr(m.fs.sj, "inlet")
    assert len(m.fs.sj.inlet.vars) == 4
    assert hasattr(m.fs.sj.inlet, "flow_vol")
    assert hasattr(m.fs.sj.inlet, "conc_mol_comp")
    assert hasattr(m.fs.sj.inlet, "temperature")
    assert hasattr(m.fs.sj.inlet, "pressure")

    assert hasattr(m.fs.sj, "outlet")
    assert len(m.fs.sj.outlet.vars) == 4
    assert hasattr(m.fs.sj.outlet, "flow_vol")
    assert hasattr(m.fs.sj.outlet, "conc_mol_comp")
    assert hasattr(m.fs.sj.outlet, "temperature")
    assert hasattr(m.fs.sj.outlet, "pressure")

    assert degrees_of_freedom(m) == 0


@pytest.mark.initialization
@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()

    m.fs.sj = StateJunction(default={"property_package": m.fs.properties})

    m.fs.sj.inlet.flow_vol.fix(1.0e-03)
    m.fs.sj.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.sj.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.sj.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.sj.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.sj.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.sj.inlet.temperature.fix(303.15)
    m.fs.sj.inlet.pressure.fix(101325.0)

    m.fs.sj.initialize(outlvl=5,
                       optarg={'tol': 1e-6})

    assert (pytest.approx(101325.0, abs=1e-2) ==
            m.fs.sj.outlet.pressure[0].value)
    assert (pytest.approx(303.15, abs=1e-2) ==
            m.fs.sj.outlet.temperature[0].value)
    assert (pytest.approx(100, abs=1e-2) ==
            m.fs.sj.outlet.conc_mol_comp[0, "EthylAcetate"].value)


@pytest.mark.ui
def test_report():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()

    m.fs.sj = StateJunction(default={"property_package": m.fs.properties})

    m.fs.sj.report()
