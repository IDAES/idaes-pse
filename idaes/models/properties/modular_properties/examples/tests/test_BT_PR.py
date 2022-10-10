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
Author: Andrew Lee
"""


import pytest

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from idaes.models.properties.modular_properties.examples.BT_PR import configuration
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from pyomo.util.check_units import assert_units_consistent

from pyomo.environ import check_optimal_termination, ConcreteModel, Objective, value

from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.models.properties.tests.test_harness import PropertyTestHarness

import idaes.logger as idaeslog

SOUT = idaeslog.INFO

# Set module level pyest marker
pytestmark = pytest.mark.cubic_root
prop_available = cubic_roots_available()


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
# Limit iterations to make sure sweeps aren;t getting out of hand
solver.options["max_iter"] = 50


class TestBTPR(PropertyTestHarness):
    def configure(self):
        self.prop_pack = GenericParameterBlock
        self.param_args = configuration
        self.prop_args = {}
        self.has_density_terms = False


# -----------------------------------------------------------------------------
# Test robustness and some outputs
class TestBTExample(object):
    @pytest.fixture()
    def m(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(**configuration)

        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        iscale.calculate_scaling_factors(m.fs.props)
        iscale.calculate_scaling_factors(m.fs.state[1])
        return m

    @pytest.mark.integration
    def test_T_sweep(self, m):
        assert_units_consistent(m)

        m.fs.obj = Objective(expr=(m.fs.state[1].temperature - 510) ** 2)
        m.fs.state[1].temperature.setub(600)

        for logP in range(8, 13, 1):
            m.fs.obj.deactivate()

            m.fs.state[1].flow_mol.fix(100)
            m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
            m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
            m.fs.state[1].temperature.fix(300)
            m.fs.state[1].pressure.fix(10 ** (0.5 * logP))

            m.fs.state.initialize()

            m.fs.state[1].temperature.unfix()
            m.fs.obj.activate()

            results = solver.solve(m)

            assert check_optimal_termination(results)
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

            assert check_optimal_termination(results)

            while m.fs.state[1].pressure.value <= 1e6:
                m.fs.state[1].pressure.value = m.fs.state[1].pressure.value + 1e5

                results = solver.solve(m)
                assert check_optimal_termination(results)
                print(T, m.fs.state[1].pressure.value)

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
        assert check_optimal_termination(results)

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
        assert check_optimal_termination(results)

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
        assert check_optimal_termination(results)

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
        assert check_optimal_termination(results)

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
        assert check_optimal_termination(results)

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

        m.fs.state[1].mole_frac_phase_comp.display()
        m.fs.state[1].enth_mol_phase_comp.display()

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
        assert check_optimal_termination(results)

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

    @pytest.mark.unit
    def test_basic_scaling(self, m):

        assert len(m.fs.state[1].scaling_factor) == 23
        assert m.fs.state[1].scaling_factor[m.fs.state[1].flow_mol] == 1e-2
        assert m.fs.state[1].scaling_factor[m.fs.state[1].flow_mol_phase["Liq"]] == 1e-2
        assert m.fs.state[1].scaling_factor[m.fs.state[1].flow_mol_phase["Vap"]] == 1e-2
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1].flow_mol_phase_comp["Liq", "benzene"]
            ]
            == 1e-2
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1].flow_mol_phase_comp["Liq", "toluene"]
            ]
            == 1e-2
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1].flow_mol_phase_comp["Vap", "benzene"]
            ]
            == 1e-2
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1].flow_mol_phase_comp["Vap", "toluene"]
            ]
            == 1e-2
        )
        assert (
            m.fs.state[1].scaling_factor[m.fs.state[1].mole_frac_comp["benzene"]]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[m.fs.state[1].mole_frac_comp["toluene"]]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]
            ]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]
            ]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]
            ]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]
            ]
            == 1000
        )
        assert m.fs.state[1].scaling_factor[m.fs.state[1].pressure] == 1e-5
        assert m.fs.state[1].scaling_factor[m.fs.state[1].temperature] == 1e-2
        assert m.fs.state[1].scaling_factor[m.fs.state[1]._teq["Vap", "Liq"]] == 1e-2
        assert m.fs.state[1].scaling_factor[m.fs.state[1]._t1_Vap_Liq] == 1e-2

        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1]._mole_frac_tbub["Vap", "Liq", "benzene"]
            ]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1]._mole_frac_tbub["Vap", "Liq", "toluene"]
            ]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1]._mole_frac_tdew["Vap", "Liq", "benzene"]
            ]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[
                m.fs.state[1]._mole_frac_tdew["Vap", "Liq", "toluene"]
            ]
            == 1000
        )
        assert (
            m.fs.state[1].scaling_factor[m.fs.state[1].temperature_bubble["Vap", "Liq"]]
            == 1e-2
        )
        assert (
            m.fs.state[1].scaling_factor[m.fs.state[1].temperature_dew["Vap", "Liq"]]
            == 1e-2
        )
