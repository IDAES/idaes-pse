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
Tests for Flash unit model.
Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           units,
                           Var)
from pyomo.util.check_units import (assert_units_consistent,
                                    assert_units_equivalent)

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
                        MomentumBalanceType)
from idaes.generic_models.unit_models.flash import Flash, EnergySplittingType
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.generic_models.properties import iapws95
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import (PhysicalParameterTestBlock,
                                     initialization_tester)
from idaes.core.util import get_solver
import idaes.core.util.scaling as iscale


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
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
        MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.ideal_separation
    assert m.fs.unit.config.has_heat_transfer
    assert m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.energy_split_basis == \
        EnergySplittingType.equal_temperature


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_calc_scale():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.unit = Flash(default={"property_package": m.fs.properties})
    iscale.calculate_scaling_factors(m)


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

        m.fs.unit.inlet.flow_mol.fix(1)
        m.fs.unit.inlet.temperature.fix(368)
        m.fs.unit.inlet.pressure.fix(101325)
        m.fs.unit.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit.inlet, "flow_mol")
        assert hasattr(btx.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.inlet, "temperature")
        assert hasattr(btx.fs.unit.inlet, "pressure")

        assert hasattr(btx.fs.unit, "vap_outlet")
        assert len(btx.fs.unit.vap_outlet.vars) == 4
        assert hasattr(btx.fs.unit.vap_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.vap_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.vap_outlet, "temperature")
        assert hasattr(btx.fs.unit.vap_outlet, "pressure")

        assert hasattr(btx.fs.unit, "liq_outlet")
        assert len(btx.fs.unit.liq_outlet.vars) == 4
        assert hasattr(btx.fs.unit.liq_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.liq_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.liq_outlet, "temperature")
        assert hasattr(btx.fs.unit.liq_outlet, "pressure")

        assert hasattr(btx.fs.unit, "split")
        assert hasattr(btx.fs.unit, "heat_duty")
        assert hasattr(btx.fs.unit, "deltaP")

        assert number_variables(btx) == 48
        assert number_total_constraints(btx) == 41
        assert number_unused_variables(btx) == 0

    @pytest.mark.component
    def test_units(self, btx):
        assert_units_consistent(btx)
        assert_units_equivalent(btx.fs.unit.heat_duty[0], units.W)
        assert_units_equivalent(btx.fs.unit.deltaP[0], units.Pa)

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
        assert (pytest.approx(0.603, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.flow_mol[0]))
        assert (pytest.approx(0.396, abs=1e-3) ==
                value(btx.fs.unit.vap_outlet.flow_mol[0]))
        assert (pytest.approx(368, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.temperature[0]))
        assert (pytest.approx(101325, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.pressure[0]))

        assert (pytest.approx(0.412, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.mole_frac_comp[0, "benzene"]))
        assert (pytest.approx(0.588, abs=1e-3) ==
                value(btx.fs.unit.liq_outlet.mole_frac_comp[0, "toluene"]))
        assert (pytest.approx(0.634, abs=1e-3) ==
                value(btx.fs.unit.vap_outlet.mole_frac_comp[0, "benzene"]))
        assert (pytest.approx(0.366, abs=1e-3) ==
                value(btx.fs.unit.vap_outlet.mole_frac_comp[0, "toluene"]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert abs(value(btx.fs.unit.inlet.flow_mol[0] -
                         (btx.fs.unit.vap_outlet.flow_mol[0] +
                          btx.fs.unit.liq_outlet.flow_mol[0]))) <= 1e-6

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

        m.fs.properties = iapws95.Iapws95ParameterBlock(default={
                "phase_presentation": iapws95.PhaseType.LG})

        m.fs.unit = Flash(default={
                "property_package": m.fs.properties,
                "ideal_separation": False,
                "energy_split_basis": EnergySplittingType.enthalpy_split})

        m.fs.unit.inlet.flow_mol.fix(100)
        m.fs.unit.inlet.enth_mol.fix(24000)
        m.fs.unit.inlet.pressure.fix(101325)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
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

    @pytest.mark.component
    def test_units(self, iapws):
        assert_units_consistent(iapws)
        assert_units_equivalent(iapws.fs.unit.heat_duty[0], units.W)
        assert_units_equivalent(iapws.fs.unit.deltaP[0], units.Pa)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
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
        assert (pytest.approx(7549.4, abs=1) ==
                value(iapws.fs.unit.liq_outlet.enth_mol[0]))
        assert (pytest.approx(59.532, abs=1e-3) ==
                value(iapws.fs.unit.liq_outlet.flow_mol[0]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
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
    @pytest.mark.component
    def test_report(self, iapws):
        iapws.fs.unit.report()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_costing(self, iapws):
        iapws.fs.unit.get_costing()
        assert isinstance(iapws.fs.unit.costing.purchase_cost, Var)
        iapws.fs.unit.diameter.fix(2)
        iapws.fs.unit.length.fix(4)
        # initialize unit with costing block
        iapws.fs.unit.initialize()
        # check costing initialized correct
        assert (pytest.approx(86957.195, abs=1e-3) ==
                value(iapws.fs.unit.costing.purchase_cost))

        results = solver.solve(iapws)
        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok
        assert (pytest.approx(63787.06525, abs=1e3) ==
                value(iapws.fs.unit.costing.base_cost))
        assert (pytest.approx(97660.6169, abs=1e3) ==
                value(iapws.fs.unit.costing.purchase_cost))

        assert_units_consistent(iapws.fs.unit.costing)
