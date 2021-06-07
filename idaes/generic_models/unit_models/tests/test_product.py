#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for product block.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models.product import Product

from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.generic_models.properties import iapws95
from idaes.generic_models.properties.examples.saponification_thermo import \
    SaponificationParameterBlock

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import (PhysicalParameterTestBlock,
                                     initialization_tester)
from idaes.core.util import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = Product(default={"property_package": m.fs.properties})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 4

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.property_package is m.fs.properties


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = Product(default={"property_package": m.fs.properties})

        m.fs.unit.flow_vol.fix(1.0e-03)
        m.fs.unit.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.temperature.fix(303.15)
        m.fs.unit.pressure.fix(101325.0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):

        assert hasattr(sapon.fs.unit, "inlet")
        assert len(sapon.fs.unit.inlet.vars) == 4
        assert hasattr(sapon.fs.unit.inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.inlet, "temperature")
        assert hasattr(sapon.fs.unit.inlet, "pressure")

        assert hasattr(sapon.fs.unit, "flow_vol")
        assert hasattr(sapon.fs.unit, "conc_mol_comp")
        assert hasattr(sapon.fs.unit, "temperature")
        assert hasattr(sapon.fs.unit, "pressure")

        assert number_variables(sapon) == 8
        assert number_total_constraints(sapon) == 0
        assert number_unused_variables(sapon) == 8

    @pytest.mark.component
    def test_units(self, sapon):
        assert_units_consistent(sapon)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    # No solve tests, as Product block has nothing to solve

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, sapon):
        sapon.fs.unit.report()


# -----------------------------------------------------------------------------
class TestBTX(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = BTXParameterBlock(default={"valid_phase": 'Liq'})

        m.fs.unit = Product(default={"property_package": m.fs.properties})

        m.fs.unit.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.temperature[0].fix(365)  # K
        m.fs.unit.pressure[0].fix(101325)  # Pa
        m.fs.unit.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.mole_frac_comp[0, "toluene"].fix(0.5)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "inlet")
        assert len(btx.fs.unit.inlet.vars) == 4
        assert hasattr(btx.fs.unit.inlet, "flow_mol")
        assert hasattr(btx.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.inlet, "temperature")
        assert hasattr(btx.fs.unit.inlet, "pressure")

        assert hasattr(btx.fs.unit, "flow_mol")
        assert hasattr(btx.fs.unit, "mole_frac_comp")
        assert hasattr(btx.fs.unit, "temperature")
        assert hasattr(btx.fs.unit, "pressure")

        assert number_variables(btx) == 8
        assert number_total_constraints(btx) == 3
        assert number_unused_variables(btx) == 2

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert (pytest.approx(5, abs=1e-3) ==
                value(btx.fs.unit.properties[0].flow_mol_phase["Liq"]))
        assert (pytest.approx(0.5, abs=1e-3) == value(
            btx.fs.unit.properties[0].mole_frac_phase_comp["Liq", "benzene"]))
        assert (pytest.approx(0.5, abs=1e-3) == value(
            btx.fs.unit.properties[0].mole_frac_phase_comp["Liq", "toluene"]))

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, btx):
        btx.fs.unit.report()


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(),
                    reason="IAPWS not available")
class TestIAPWS(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = Product(default={"property_package": m.fs.properties})

        m.fs.unit.flow_mol[0].fix(100)
        m.fs.unit.enth_mol[0].fix(4000)
        m.fs.unit.pressure[0].fix(101325)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iapws):
        assert len(iapws.fs.unit.inlet.vars) == 3
        assert hasattr(iapws.fs.unit.inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.inlet, "pressure")

        assert hasattr(iapws.fs.unit, "flow_mol")
        assert hasattr(iapws.fs.unit, "enth_mol")
        assert hasattr(iapws.fs.unit, "pressure")

        assert number_variables(iapws) == 3
        assert number_total_constraints(iapws) == 0
        assert number_unused_variables(iapws) == 3

    @pytest.mark.component
    def test_units(self, iapws):
        assert_units_consistent(iapws)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    # No solve as there is nothing to solve for

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, iapws):
        iapws.fs.unit.report()
