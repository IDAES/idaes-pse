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
Tests for translator block.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.unit_models.translator import Translator
from idaes.property_models.ideal.BTX_ideal_VLE import (
                        BTXParameterBlock)
from idaes.property_models.examples.saponification_thermo import (
                        SaponificationParameterBlock)


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

    m.fs.ideal = BTXParameterBlock()
    m.fs.sap = SaponificationParameterBlock()

    m.fs.trans = Translator(
                default={"inlet_property_package": m.fs.ideal,
                         "outlet_property_package": m.fs.sap})

    assert hasattr(m.fs.trans, "inlet")
    assert len(m.fs.trans.inlet.vars) == 4
    assert hasattr(m.fs.trans.inlet, "flow_mol")
    assert hasattr(m.fs.trans.inlet, "mole_frac")
    assert hasattr(m.fs.trans.inlet, "temperature")
    assert hasattr(m.fs.trans.inlet, "pressure")

    assert hasattr(m.fs.trans, "outlet")
    assert len(m.fs.trans.outlet.vars) == 4
    assert hasattr(m.fs.trans.outlet, "flow_vol")
    assert hasattr(m.fs.trans.outlet, "conc_mol_comp")
    assert hasattr(m.fs.trans.outlet, "temperature")
    assert hasattr(m.fs.trans.outlet, "pressure")


@pytest.mark.initialization
@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    # Very basic test of initialization routine - only checks data flow
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.sap = SaponificationParameterBlock()

    m.fs.trans = Translator(
                default={"inlet_property_package": m.fs.sap,
                         "outlet_property_package": m.fs.sap})

    m.fs.trans.initialize(outlvl=5,
                          optarg={'tol': 1e-6})


@pytest.mark.ui
def test_report():
    # Very basic test of initialization routine - only checks data flow
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.sap = SaponificationParameterBlock()

    m.fs.trans = Translator(
                default={"inlet_property_package": m.fs.sap,
                         "outlet_property_package": m.fs.sap})

    m.fs.trans.report()
