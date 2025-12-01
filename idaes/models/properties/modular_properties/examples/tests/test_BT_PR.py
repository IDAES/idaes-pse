#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Authors: Andrew Lee, Douglas Allan
"""

import pytest

from numpy import logspace

from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import (
    assert_optimal_termination,
    check_optimal_termination,
    ConcreteModel,
    Objective,
    value,
)

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from idaes.models.properties.modular_properties.examples.BT_PR import configuration
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.models.properties.tests.test_harness import PropertyTestHarness

import idaes.logger as idaeslog

SOUT = idaeslog.INFO

# Set module level pyest marker
pytestmark = pytest.mark.cubic_root


# -----------------------------------------------------------------------------
# Get default solver for testing
# Limit iterations to make sure sweeps aren't getting out of hand
solver = get_solver(solver="ipopt_v2", solver_options={"max_iter": 100})


@pytest.mark.skipif(not cubic_roots_available(), reason="Cubic functions not available")
class TestBTPR(PropertyTestHarness):
    def configure(self):
        self.prop_pack = GenericParameterBlock
        self.param_args = configuration
        self.prop_args = {}
        self.has_density_terms = False


# -----------------------------------------------------------------------------
# Test robustness and some outputs
@pytest.mark.skipif(not cubic_roots_available(), reason="Cubic functions not available")
class TestBTExampleLegacyScaling(object):
    @pytest.fixture()
    def m(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(**configuration)

        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        # Set small values of epsilon to get accurate results
        # Initialization will handle finding the correct region
        m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
        m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)

        iscale.calculate_scaling_factors(m.fs.props)
        iscale.calculate_scaling_factors(m.fs.state[1])
        return m

    @pytest.mark.integration
    def test_T_sweep(self, m):
        assert_units_consistent(m)

        m.fs.obj = Objective(expr=(m.fs.state[1].temperature - 510) ** 2)
        m.fs.state[1].temperature.setub(600)

        for P in logspace(4.8, 5.9, 8):

            m.fs.state[1].flow_mol.fix(100)
            m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
            m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
            m.fs.state[1].temperature.fix(300)
            m.fs.state[1].pressure.fix(P)

            # For optimization sweep, use a large eps to avoid getting stuck at
            # bubble and dew points
            m.fs.state[1].eps_t_Vap_Liq.set_value(10)
            m.fs.state[1].eps_z_Vap_Liq.set_value(10)

            m.fs.state.initialize()

            m.fs.state[1].temperature.unfix()
            m.fs.obj.activate()

            results = solver.solve(m)
            assert_optimal_termination(results)

            # Switch to small eps and re-solve to refine result
            m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
            m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)

            results = solver.solve(m)

            assert_optimal_termination(results)
            assert m.fs.state[1].flow_mol_phase["Liq"].value <= 1e-5

    @pytest.mark.integration
    def test_P_sweep(self, m):
        for T in range(370, 500, 25):
            m.fs.state[1].flow_mol.fix(100)
            m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
            m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
            m.fs.state[1].temperature.fix(T)
            m.fs.state[1].pressure.fix(1e5)

            m.fs.state.initialize()

            results = solver.solve(m)

            assert_optimal_termination(results)

            while m.fs.state[1].pressure.value <= 1e6:

                results = solver.solve(m)
                assert_optimal_termination(results)

                m.fs.state[1].pressure.value = m.fs.state[1].pressure.value + 1e5

    @pytest.mark.component
    def test_T350_P1_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(350)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), abs=1e-1) == 365
        assert 0.0035346 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.966749 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.894676
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.347566
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.971072
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.959791
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.70584
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.29416
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38942.8
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 78048.7
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -361.794
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -264.0181
        )

    @pytest.mark.component
    def test_T350_P5_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(350)
        m.fs.state[1].pressure.fix(5e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 431.47
        assert (
            pytest.approx(value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5)
            == 0.01766
        )
        assert (
            pytest.approx(value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5)
            == 0.80245
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.181229
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.070601
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.856523
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.799237
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.65415
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.34585
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38966.9
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 75150.7
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -361.8433
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -281.9703
        )

    @pytest.mark.component
    def test_T450_P1_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(450)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 371.4
        assert 0.0033583 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.9821368 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 8.069323
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 4.304955
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.985365
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.979457
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.29861
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.70139
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.5
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 49441.2
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 84175.1
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -328.766
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -241.622
        )

    @pytest.mark.component
    def test_T450_P5_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(450)
        m.fs.state[1].pressure.fix(5e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 436.93
        assert 0.0166181 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.9053766 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 1.63308
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.873213
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.927534
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.898324
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.3488737
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.6511263
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.5
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 51095.2
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 83362.3
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -326.299
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -256.198
        )

    @pytest.mark.component
    def test_T368_P1_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(368)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 368
        assert 0.003504 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.97 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 1.492049
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.621563
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.97469
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.964642
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.4012128
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.5987872
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.6141738
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.3858262
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38235.1
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 77155.4
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -359.256
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -262.348
        )

    @pytest.mark.component
    def test_T376_P1_x2(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.2)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.8)
        m.fs.state[1].temperature.fix(376)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 376
        assert 0.00361333 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.968749 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 1.8394188
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.7871415
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.9763608
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.9663611
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.17342
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.82658
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.3267155
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.6732845
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 31535.8
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 69175.3
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -369.033
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -273.513
        )


# -----------------------------------------------------------------------------
# Test robustness and some outputs
@pytest.mark.skipif(not cubic_roots_available(), reason="Cubic functions not available")
class TestBTExampleScalerObject(object):
    @pytest.fixture()
    def m(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(**configuration)

        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        # Set small values of epsilon to get accurate results
        # Initialization will handle finding the correct region
        m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
        m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)

        # Try to prevent the equilibrium constraints from
        # being satisfied by using an extremely supercritical
        # temperature
        m.fs.state[1]._teq[("Vap", "Liq")].setub(600)

        scaler_obj = m.fs.state[1].default_scaler()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 0.01
        scaler_obj.default_scaling_factors["mole_frac_phase_comp"] = 2
        scaler_obj.default_scaling_factors["temperature"] = 1e-2
        scaler_obj.default_scaling_factors["pressure"] = 1e-5
        scaler_obj.default_scaling_factors["enth_mol_phase"] = 1e-4
        scaler_obj.scale_model(m.fs.state[1])
        return m

    @pytest.mark.integration
    def test_T_sweep(self, m):
        assert_units_consistent(m)

        m.fs.obj = Objective(expr=(m.fs.state[1].temperature - 510) ** 2)
        m.fs.state[1].temperature.setub(600)

        for P in logspace(4.8, 5.9, 8):

            m.fs.state[1].flow_mol.fix(100)
            m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
            m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
            m.fs.state[1].temperature.fix(300)
            m.fs.state[1].pressure.fix(P)

            # For optimization sweep, use a large eps to avoid getting stuck at
            # bubble and dew points
            m.fs.state[1].eps_t_Vap_Liq.set_value(10)
            m.fs.state[1].eps_z_Vap_Liq.set_value(10)

            m.fs.state.initialize()

            m.fs.state[1].temperature.unfix()
            m.fs.obj.activate()

            results = solver.solve(m)
            assert_optimal_termination(results)

            # Switch to small eps and re-solve to refine result
            m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
            m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)

            results = solver.solve(m)

            assert_optimal_termination(results)
            assert m.fs.state[1].flow_mol_phase["Liq"].value <= 1e-5

    @pytest.mark.integration
    def test_P_sweep(self, m):
        for T in range(370, 500, 25):
            m.fs.state[1].flow_mol.fix(100)
            m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
            m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
            m.fs.state[1].temperature.fix(T)
            m.fs.state[1].pressure.fix(1e5)

            m.fs.state.initialize()

            # Use a less strict complementarity condition
            # to encourage convergence.
            m.fs.state[1].eps_t_Vap_Liq.set_value(1e-2)
            m.fs.state[1].eps_z_Vap_Liq.set_value(1e-2)

            results = solver.solve(m)

            assert_optimal_termination(results)

            while m.fs.state[1].pressure.value <= 1e6:

                results = solver.solve(m)

                assert_optimal_termination(results)

                m.fs.state[1].pressure.value = m.fs.state[1].pressure.value + 1e5

            # Clean up for later tests, which may require
            # greater precision in the flash calculations
            m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
            m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)

    @pytest.mark.component
    def test_T350_P1_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(350)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), abs=1e-1) == 365
        assert 0.0035346 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.966749 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.894676
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.347566
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.971072
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.959791
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.70584
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.29416
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38942.8
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 78048.7
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -361.794
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -264.0181
        )

    @pytest.mark.component
    def test_T350_P5_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(350)
        m.fs.state[1].pressure.fix(5e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 431.47
        assert (
            pytest.approx(value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5)
            == 0.01766
        )
        assert (
            pytest.approx(value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5)
            == 0.80245
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.181229
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.070601
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.856523
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.799237
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.65415
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.34585
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38966.9
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 75150.7
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -361.8433
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -281.9703
        )

    @pytest.mark.component
    def test_T450_P1_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(450)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 371.4
        assert 0.0033583 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.9821368 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 8.069323
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 4.304955
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.985365
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.979457
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.29861
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.70139
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.5
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 49441.2
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 84175.1
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -328.766
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -241.622
        )

    @pytest.mark.component
    def test_T450_P5_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(450)
        m.fs.state[1].pressure.fix(5e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 436.93
        assert 0.0166181 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.9053766 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 1.63308
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.873213
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.927534
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.898324
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.3488737
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.6511263
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.5
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 51095.2
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 83362.3
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -326.299
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -256.198
        )

    @pytest.mark.component
    def test_T368_P1_x5(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
        m.fs.state[1].temperature.fix(368)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 368
        assert 0.003504 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.97 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 1.492049
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.621563
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.97469
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.964642
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.4012128
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.5987872
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.6141738
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.3858262
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 38235.1
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 77155.4
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -359.256
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -262.348
        )

    @pytest.mark.component
    def test_T376_P1_x2(self, m):
        m.fs.state[1].flow_mol.fix(100)
        m.fs.state[1].mole_frac_comp["benzene"].fix(0.2)
        m.fs.state[1].mole_frac_comp["toluene"].fix(0.8)
        m.fs.state[1].temperature.fix(376)
        m.fs.state[1].pressure.fix(1e5)

        # Trigger build of enthalpy and entropy
        m.fs.state[1].enth_mol_phase
        m.fs.state[1].entr_mol_phase

        m.fs.state.initialize(outlvl=SOUT)

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        assert pytest.approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-5) == 376
        assert 0.00361333 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-5
        )
        assert 0.968749 == pytest.approx(
            value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-5
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 1.8394188
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.7871415
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.9763608
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.9663611
        )

        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-5
            )
            == 0.17342
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-5
            )
            == 0.82658
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-5
            )
            == 0.3267155
        )
        assert (
            pytest.approx(
                value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-5
            )
            == 0.6732845
        )

        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-5) == 31535.8
        )
        assert (
            pytest.approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-5) == 69175.3
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-5) == -369.033
        )
        assert (
            pytest.approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-5) == -273.513
        )
