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
Tests for feed block.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.unit_models.feed import Feed
from idaes.property_models.examples.saponification_thermo import (
                        SaponificationParameterBlock)
from idaes.core.util.model_statistics import calculate_degrees_of_freedom


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

    m.fs.feed = Feed(default={"property_package": m.fs.properties})

    assert hasattr(m.fs.feed, "flow_vol")
    assert hasattr(m.fs.feed, "conc_mol_comp")
    assert hasattr(m.fs.feed, "temperature")
    assert hasattr(m.fs.feed, "pressure")

    assert hasattr(m.fs.feed, "outlet")
    assert len(m.fs.feed.outlet.vars) == 4
    assert hasattr(m.fs.feed.outlet, "flow_vol")
    assert hasattr(m.fs.feed.outlet, "conc_mol_comp")
    assert hasattr(m.fs.feed.outlet, "temperature")
    assert hasattr(m.fs.feed.outlet, "pressure")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()

    m.fs.feed = Feed(default={"property_package": m.fs.properties})

    m.fs.feed.flow_vol.fix(1.0e-03)
    m.fs.feed.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.feed.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.feed.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.feed.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.feed.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.feed.temperature.fix(303.15)
    m.fs.feed.pressure.fix(101325.0)

    assert calculate_degrees_of_freedom(m) == 0

    m.fs.feed.initialize(outlvl=5,
                         optarg={'tol': 1e-6})

    assert m.fs.feed.outlet.flow_vol[0].value == 1.0e-03
    assert m.fs.feed.outlet.conc_mol_comp[0, "H2O"].value == 55388.0
    assert m.fs.feed.outlet.conc_mol_comp[0, "NaOH"].value == 100.0
    assert m.fs.feed.outlet.conc_mol_comp[0, "EthylAcetate"].value == 100.0
    assert m.fs.feed.outlet.conc_mol_comp[0, "SodiumAcetate"].value == 0.0
    assert m.fs.feed.outlet.conc_mol_comp[0, "Ethanol"].value == 0.0
