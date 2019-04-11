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
Tests for Node unit model.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.unit_models.node import Node
from idaes.property_models.examples.saponification_thermo import (
    SaponificationParameterBlock)
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

    m.fs.properties = SaponificationParameterBlock()

    m.fs.node = Node(default={"property_package": m.fs.properties})

    assert hasattr(m.fs.node, "inlet")
    assert len(m.fs.node.inlet.vars) == 4
    assert hasattr(m.fs.node.inlet, "flow_vol")
    assert hasattr(m.fs.node.inlet, "conc_mol_comp")
    assert hasattr(m.fs.node.inlet, "temperature")
    assert hasattr(m.fs.node.inlet, "pressure")

    assert hasattr(m.fs.node, "outlet")
    assert len(m.fs.node.outlet.vars) == 4
    assert hasattr(m.fs.node.outlet, "flow_vol")
    assert hasattr(m.fs.node.outlet, "conc_mol_comp")
    assert hasattr(m.fs.node.outlet, "temperature")
    assert hasattr(m.fs.node.outlet, "pressure")

    assert degrees_of_freedom(m) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()

    m.fs.node = Node(default={"property_package": m.fs.properties})

    m.fs.node.inlet.flow_vol.fix(1.0e-03)
    m.fs.node.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.node.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.node.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.node.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.node.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.node.inlet.temperature.fix(303.15)
    m.fs.node.inlet.pressure.fix(101325.0)

    m.fs.node.initialize(outlvl=5,
                         optarg={'tol': 1e-6})

    assert (pytest.approx(101325.0, abs=1e-2) ==
            m.fs.node.outlet.pressure[0].value)
    assert (pytest.approx(303.15, abs=1e-2) ==
            m.fs.node.outlet.temperature[0].value)
    assert (pytest.approx(100, abs=1e-2) ==
            m.fs.node.outlet.conc_mol_comp[0, "EthylAcetate"].value)
