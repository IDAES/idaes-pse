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
from pyomo.environ import ConcreteModel, SolverFactory, value
from idaes.core import FlowsheetBlock
from idaes.unit_models.feed_flash import FeedFlash, FlashType
from idaes.property_models import iapws95
from idaes.property_models.ideal.BTX_ideal_VLE import BTXParameterBlock
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)


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
class TestBTX(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = BTXParameterBlock()

        m.fs.unit = FeedFlash(default={"property_package": m.fs.properties})

        return m

    @pytest.mark.build
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "flow_mol")
        assert hasattr(btx.fs.unit, "mole_frac")
        assert hasattr(btx.fs.unit, "temperature")
        assert hasattr(btx.fs.unit, "pressure")

        assert hasattr(btx.fs.unit, "outlet")
        assert len(btx.fs.unit.outlet.vars) == 4
        assert hasattr(btx.fs.unit.outlet, "flow_mol")
        assert hasattr(btx.fs.unit.outlet, "mole_frac")
        assert hasattr(btx.fs.unit.outlet, "temperature")
        assert hasattr(btx.fs.unit.outlet, "pressure")

        assert hasattr(btx.fs.unit, "isothermal")

        assert number_variables(btx) == 36
        assert number_total_constraints(btx) == 31
        assert number_unused_variables(btx) == 0

    def test_dof(self, btx):
        btx.fs.unit.flow_mol.fix(1)
        btx.fs.unit.temperature.fix(368)
        btx.fs.unit.pressure.fix(101325)
        btx.fs.unit.mole_frac[0, "benzene"].fix(0.5)
        btx.fs.unit.mole_frac[0, "toluene"].fix(0.5)

        assert degrees_of_freedom(btx) == 0

    @pytest.mark.initialization
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, btx):
        orig_fixed_vars = fixed_variables_set(btx)
        orig_act_consts = activated_constraints_set(btx)

        btx.fs.unit.initialize(optarg={'tol': 1e-6})

        assert degrees_of_freedom(btx) == 0

        fin_fixed_vars = fixed_variables_set(btx)
        fin_act_consts = activated_constraints_set(btx)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, btx):
        assert (pytest.approx(101325.0, abs=1e3) ==
                btx.fs.unit.outlet.pressure[0].value)
        assert (pytest.approx(368.00, abs=1e-0) ==
                btx.fs.unit.outlet.temperature[0].value)
        assert (pytest.approx(1.0, abs=1e-2) ==
                btx.fs.unit.outlet.flow_mol[0].value)
        assert (pytest.approx(0.355, abs=1e-3) ==
                btx.fs.unit.control_volume.
                properties_out[0].flow_mol_phase["Vap"].value)

    @pytest.mark.ui
    def test_report(self, btx):
        btx.fs.unit.report()


# -----------------------------------------------------------------------------
@pytest.mark.iapws
class TestIAPWS(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = iapws95.Iapws95ParameterBlock(default={
                "phase_presentation": iapws95.PhaseType.LG})

        m.fs.unit = FeedFlash(default={"property_package": m.fs.properties,
                                       "flash_type": FlashType.isenthalpic})

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

        assert hasattr(iapws.fs.unit, "isenthalpic")

        assert number_variables(iapws) == 7
        assert number_total_constraints(iapws) == 4
        assert number_unused_variables(iapws) == 0

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

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, iapws):
        assert (pytest.approx(101325.0, abs=1e3) ==
                iapws.fs.unit.outlet.pressure[0].value)
        assert (pytest.approx(24000, abs=1e3) ==
                iapws.fs.unit.outlet.enth_mol[0].value)
        assert (pytest.approx(100.0, abs=1e-2) ==
                iapws.fs.unit.outlet.flow_mol[0].value)

        assert (pytest.approx(373.12, abs=1e-2) == value(
            iapws.fs.unit.control_volume.properties_out[0].temperature))
        assert (pytest.approx(0.5953, abs=1e-4) == value(
            iapws.fs.unit.control_volume.properties_out[0].phase_frac["Liq"]))

    @pytest.mark.ui
    def test_report(self, iapws):
        iapws.fs.unit.report()
