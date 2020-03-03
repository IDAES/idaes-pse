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
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value)
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models.feed import Feed
from idaes.generic_models.properties.examples.saponification_thermo import (
                        SaponificationParameterBlock)
from idaes.generic_models.properties import iapws95
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from idaes.core.util.testing import (get_default_solver,
                                     PhysicalParameterTestBlock)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = Feed(default={"property_package": m.fs.properties})

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

        m.fs.unit = Feed(default={"property_package": m.fs.properties})

        return m

    @pytest.mark.build
    def test_build(self, sapon):

        assert hasattr(sapon.fs.unit, "flow_vol")
        assert hasattr(sapon.fs.unit, "conc_mol_comp")
        assert hasattr(sapon.fs.unit, "temperature")
        assert hasattr(sapon.fs.unit, "pressure")

        assert hasattr(sapon.fs.unit, "outlet")
        assert len(sapon.fs.unit.outlet.vars) == 4
        assert hasattr(sapon.fs.unit.outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.outlet, "temperature")
        assert hasattr(sapon.fs.unit.outlet, "pressure")

        assert number_variables(sapon) == 8
        assert number_total_constraints(sapon) == 0
        assert number_unused_variables(sapon) == 8

    def test_dof(self, sapon):
        sapon.fs.unit.flow_vol.fix(1.0e-03)
        sapon.fs.unit.conc_mol_comp[0, "H2O"].fix(55388.0)
        sapon.fs.unit.conc_mol_comp[0, "NaOH"].fix(100.0)
        sapon.fs.unit.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        sapon.fs.unit.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        sapon.fs.unit.conc_mol_comp[0, "Ethanol"].fix(0.0)

        sapon.fs.unit.temperature.fix(303.15)
        sapon.fs.unit.pressure.fix(101325.0)

        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, sapon):
        orig_fixed_vars = fixed_variables_set(sapon)
        orig_act_consts = activated_constraints_set(sapon)

        sapon.fs.unit.initialize(optarg={'tol': 1e-6})

        assert degrees_of_freedom(sapon) == 0

        fin_fixed_vars = fixed_variables_set(sapon)
        fin_act_consts = activated_constraints_set(sapon)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    # No solve, as nothing to solve for

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, sapon):
        assert (pytest.approx(101325.0, abs=1e-2) ==
                value(sapon.fs.unit.outlet.pressure[0]))
        assert (pytest.approx(303.15, abs=1e-2) ==
                value(sapon.fs.unit.outlet.temperature[0]))
        assert (pytest.approx(1e-3, abs=1e-5) ==
                value(sapon.fs.unit.outlet.flow_vol[0]))
        assert (pytest.approx(55388, abs=1e0) ==
                value(sapon.fs.unit.outlet.conc_mol_comp[0, "H2O"]))
        assert (pytest.approx(100.0, abs=1e-2) ==
                value(sapon.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]))
        assert (pytest.approx(100.0, abs=1e-2) ==
                value(sapon.fs.unit.outlet.conc_mol_comp[0, "NaOH"]))
        assert (pytest.approx(0.00, abs=1e-2) ==
                value(sapon.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]))
        assert (pytest.approx(0.00, abs=1e-2) ==
                value(sapon.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]))

    @pytest.mark.ui
    def test_report(self, sapon):
        sapon.fs.unit.report()


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(),
                    reason="IAPWS not available")
class TestIAPWS(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = iapws95.Iapws95ParameterBlock(default={
                "phase_presentation": iapws95.PhaseType.LG})

        m.fs.unit = Feed(default={"property_package": m.fs.properties})

        return m

    @pytest.mark.build
    def test_build(self, iapws):
        assert hasattr(iapws.fs.unit, "flow_mol")
        assert hasattr(iapws.fs.unit, "enth_mol")
        assert hasattr(iapws.fs.unit, "pressure")

        assert hasattr(iapws.fs.unit, "outlet")
        assert len(iapws.fs.unit.outlet.vars) == 3
        assert hasattr(iapws.fs.unit.outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.outlet, "pressure")

        assert number_variables(iapws) == 3
        assert number_total_constraints(iapws) == 0
        assert number_unused_variables(iapws) == 3

    def test_dof(self, iapws):
        iapws.fs.unit.flow_mol.fix(100)
        iapws.fs.unit.enth_mol.fix(24000)
        iapws.fs.unit.pressure.fix(101325)

        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.initialization
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, iapws):
        orig_fixed_vars = fixed_variables_set(iapws)
        orig_act_consts = activated_constraints_set(iapws)

        iapws.fs.unit.initialize(optarg={'tol': 1e-6})

        assert degrees_of_freedom(iapws) == 0

        fin_fixed_vars = fixed_variables_set(iapws)
        fin_act_consts = activated_constraints_set(iapws)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    # No solve test, as nothing to solve for

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, iapws):
        assert (pytest.approx(101325.0, abs=1e3) ==
                value(iapws.fs.unit.outlet.pressure[0]))
        assert (pytest.approx(24000, abs=1e3) ==
                value(iapws.fs.unit.outlet.enth_mol[0]))
        assert (pytest.approx(100.0, abs=1e-2) ==
                value(iapws.fs.unit.outlet.flow_mol[0]))

        assert (pytest.approx(373.12, abs=1e-2) == value(
            iapws.fs.unit.properties[0].temperature))
        assert (pytest.approx(0.5953, abs=1e-4) == value(
            iapws.fs.unit.properties[0].phase_frac["Liq"]))

    @pytest.mark.ui
    def test_report(self, iapws):
        iapws.fs.unit.report()
