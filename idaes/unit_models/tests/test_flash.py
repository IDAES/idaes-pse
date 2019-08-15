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
Tests for Flash unit model.
Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value)

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType)
from idaes.unit_models.flash import Flash, EnergySplittingType
from idaes.property_models.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.property_models import iapws95
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

    m.fs.unit = Flash(default={"property_package": m.fs.properties})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 11

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == \
        MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == \
        EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == \
        MomentumBalanceType.useDefault
    assert m.fs.unit.config.ideal_separation
    assert m.fs.unit.config.has_heat_transfer
    assert m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.energy_split_basis == \
        EnergySplittingType.equal_temperature


# -----------------------------------------------------------------------------
class TestBTXIdeal(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = BTXParameterBlock(default={"valid_phase":
                                                     ('Liq', 'Vap'),
                                                     "activity_coeff_model":
                                                     "Ideal"})

        m.fs.unit = Flash(default={"property_package": m.fs.properties})

        return m

    @pytest.mark.build
    def test_build(self, btx):
        assert hasattr(btx.fs.unit.inlet, "flow_mol")
        assert hasattr(btx.fs.unit.inlet, "mole_frac")
        assert hasattr(btx.fs.unit.inlet, "temperature")
        assert hasattr(btx.fs.unit.inlet, "pressure")

        assert hasattr(btx.fs.unit, "vap_outlet")
        assert len(btx.fs.unit.vap_outlet.vars) == 4
        assert hasattr(btx.fs.unit.vap_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.vap_outlet, "mole_frac")
        assert hasattr(btx.fs.unit.vap_outlet, "temperature")
        assert hasattr(btx.fs.unit.vap_outlet, "pressure")

        assert hasattr(btx.fs.unit, "liq_outlet")
        assert len(btx.fs.unit.liq_outlet.vars) == 4
        assert hasattr(btx.fs.unit.liq_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.liq_outlet, "mole_frac")
        assert hasattr(btx.fs.unit.liq_outlet, "temperature")
        assert hasattr(btx.fs.unit.liq_outlet, "pressure")

        assert hasattr(btx.fs.unit, "split")
        assert hasattr(btx.fs.unit, "heat_duty")
        assert hasattr(btx.fs.unit, "deltaP")

        assert number_variables(btx) == 48
        assert number_total_constraints(btx) == 41
        assert number_unused_variables(btx) == 0

    def test_dof(self, btx):
        btx.fs.unit.inlet.flow_mol.fix(1)
        btx.fs.unit.inlet.temperature.fix(368)
        btx.fs.unit.inlet.pressure.fix(101325)
        btx.fs.unit.inlet.mole_frac[0, "benzene"].fix(0.5)
        btx.fs.unit.inlet.mole_frac[0, "toluene"].fix(0.5)

        btx.fs.unit.heat_duty.fix(0)
        btx.fs.unit.deltaP.fix(0)

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

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, btx):
        assert (pytest.approx(0.603, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.flow_mol[0]))
        assert (pytest.approx(0.396, abs=1e-3) ==
                value(btx.fs.unit.vap_outlet.flow_mol[0]))
        assert (pytest.approx(368, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.temperature[0]))
        assert (pytest.approx(101325, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.pressure[0]))

        assert (pytest.approx(0.412, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.mole_frac[0, "benzene"]))
        assert (pytest.approx(0.588, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.mole_frac[0, "toluene"]))
        assert (pytest.approx(0.634, abs=1e-3) ==
                value(btx.fs.unit.vap_outlet.mole_frac[0, "benzene"]))
        assert (pytest.approx(0.366, abs=1e-3) ==
                value(btx.fs.unit.vap_outlet.mole_frac[0, "toluene"]))

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_conservation(self, btx):
        assert abs(value(btx.fs.unit.inlet.flow_mol[0] -
                         (btx.fs.unit.vap_outlet.flow_mol[0] +
                          btx.fs.unit.liq_outlet.flow_mol[0]))) <= 1e-6

    @pytest.mark.ui
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

        m.fs.properties = iapws95.Iapws95ParameterBlock(default={
                "phase_presentation": iapws95.PhaseType.LG})

        m.fs.unit = Flash(default={
                "property_package": m.fs.properties,
                "ideal_separation": False,
                "energy_split_basis": EnergySplittingType.enthalpy_split})

        return m

    @pytest.mark.build
    def test_build(self, iapws):
        assert len(iapws.fs.unit.inlet.vars) == 3
        assert hasattr(iapws.fs.unit.inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.inlet, "pressure")

        assert hasattr(iapws.fs.unit, "vap_outlet")
        assert len(iapws.fs.unit.vap_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.vap_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.vap_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.vap_outlet, "pressure")

        assert hasattr(iapws.fs.unit, "liq_outlet")
        assert len(iapws.fs.unit.liq_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.liq_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.liq_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.liq_outlet, "pressure")

        assert hasattr(iapws.fs.unit, "split")
        assert hasattr(iapws.fs.unit, "heat_duty")
        assert hasattr(iapws.fs.unit, "deltaP")

        assert number_variables(iapws) == 18
        assert number_total_constraints(iapws) == 13
        assert number_unused_variables(iapws) == 0

    def test_dof(self, iapws):
        iapws.fs.unit.inlet.flow_mol.fix(100)
        iapws.fs.unit.inlet.enth_mol.fix(24000)
        iapws.fs.unit.inlet.pressure.fix(101325)

        iapws.fs.unit.heat_duty.fix(0)
        iapws.fs.unit.deltaP.fix(0)

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

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, iapws):
        t_in = value(iapws.fs.unit.control_volume.properties_in[0].temperature)
        assert (pytest.approx(373.12, abs=1e-2) == t_in)
        assert (pytest.approx(t_in, abs=1e-2) == value(
                iapws.fs.unit.control_volume.properties_out[0].temperature))
        assert (pytest.approx(t_in, abs=1e-2) == value(
                iapws.fs.unit.split.Vap_state[0].temperature))
        assert (pytest.approx(t_in, abs=1e-2) == value(
                iapws.fs.unit.split.Liq_state[0].temperature))

        assert (pytest.approx(101325.0, abs=1e3) ==
                value(iapws.fs.unit.vap_outlet.pressure[0]))
        assert (pytest.approx(48200, abs=1e3) ==
                value(iapws.fs.unit.vap_outlet.enth_mol[0]))
        assert (pytest.approx(40.467, abs=1e-3) ==
                value(iapws.fs.unit.vap_outlet.flow_mol[0]))

        assert (pytest.approx(101325.0, abs=1e3) ==
                value(iapws.fs.unit.liq_outlet.pressure[0]))
        assert (pytest.approx(7549.4, abs=1e-1) ==
                value(iapws.fs.unit.liq_outlet.enth_mol[0]))
        assert (pytest.approx(59.532, abs=1e-3) ==
                value(iapws.fs.unit.liq_outlet.flow_mol[0]))

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_conservation(self, iapws):
        assert abs(value(iapws.fs.unit.inlet.flow_mol[0] -
                         (iapws.fs.unit.vap_outlet.flow_mol[0] +
                          iapws.fs.unit.liq_outlet.flow_mol[0]))) <= 1e-6
        assert abs(value(iapws.fs.unit.inlet.flow_mol[0] *
                         iapws.fs.unit.inlet.enth_mol[0] -
                         (iapws.fs.unit.vap_outlet.flow_mol[0] *
                          iapws.fs.unit.vap_outlet.enth_mol[0] +
                          iapws.fs.unit.liq_outlet.flow_mol[0] *
                          iapws.fs.unit.liq_outlet.enth_mol[0]))) <= 1e-6

    @pytest.mark.ui
    def test_report(self, iapws):
        iapws.fs.unit.report()
