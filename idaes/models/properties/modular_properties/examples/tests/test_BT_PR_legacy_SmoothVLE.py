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
    check_optimal_termination,
    ConcreteModel,
    Objective,
    TransformationFactory,
    units as pyunits,
    value,
)

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.core.solvers import get_solver
from idaes.core.scaling import get_scaling_factor

from idaes.models.properties.tests.test_harness import PropertyTestHarness
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.pure import RPP4

import idaes.logger as idaeslog

SOUT = idaeslog.INFO

# Set module level pyest marker
pytestmark = pytest.mark.cubic_root


# -----------------------------------------------------------------------------
# Get default solver for testing
# Limit iterations to make sure sweeps aren't getting out of hand
solver = get_solver(
    solver="ipopt_v2",
    solver_options={"max_iter": 50},
    writer_config={
        "scale_model": True,
        "linear_presolve": True,
    },
)

# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019

configuration = {
    # Specifying components
    "components": {
        "benzene": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (78.1136e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                "temperature_crit": (562.2, pyunits.K),  # [1]
                "omega": 0.212,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (-3.392e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (4.739e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (-3.017e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (7.130e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (82.9e3, pyunits.J / pyunits.mol),  # [3]
                "entr_mol_form_vap_comp_ref": (
                    -269,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [3]
                "pressure_sat_comp_coeff": {
                    "A": (-6.98273, None),  # [1]
                    "B": (1.33213, None),
                    "C": (-2.62863, None),
                    "D": (-3.33399, None),
                },
            },
        },
        "toluene": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (92.1405e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (41e5, pyunits.Pa),  # [1]
                "temperature_crit": (591.8, pyunits.K),  # [1]
                "omega": 0.263,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (-2.435e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (5.125e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (-2.765e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (4.911e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (50.1e3, pyunits.J / pyunits.mol),  # [3]
                "entr_mol_form_vap_comp_ref": (
                    -321,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [3]
                "pressure_sat_comp_coeff": {
                    "A": (-7.28607, None),  # [1]
                    "B": (1.38091, None),
                    "C": (-2.83433, None),
                    "D": (-2.79168, None),
                },
            },
        },
    },
    # Specifying phases
    "phases": {
        "Liq": {
            "type": LiquidPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        },
        "Vap": {
            "type": VaporPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        },
    },
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 500, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {
        "PR_kappa": {
            ("benzene", "benzene"): 0.000,
            ("benzene", "toluene"): 0.000,
            ("toluene", "benzene"): 0.000,
            ("toluene", "toluene"): 0.000,
        }
    },
}


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
class TestBTExample(object):
    @pytest.fixture()
    def m(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = GenericParameterBlock(**configuration)

        m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

        scaler_obj = m.fs.state.default_scaler()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 0.01
        scaler_obj.scale_model(m.fs.state[1])

        return m

    @pytest.mark.component
    def test_scaling(self, m):
        assert len(m.fs.state[1].scaling_factor) == 57
        assert len(m.fs.state[1].scaling_hint) == 6

        # Variables
        assert get_scaling_factor(m.fs.state[1].flow_mol) == 1e-2
        for vdata in m.fs.state[1].mole_frac_comp.values():
            assert get_scaling_factor(vdata) == 10
        assert get_scaling_factor(m.fs.state[1].pressure) == 1e-5
        assert get_scaling_factor(m.fs.state[1].temperature) == 1 / 300
        for vdata in m.fs.state[1].flow_mol_phase.values():
            assert get_scaling_factor(vdata) == 1e-2
        for vdata in m.fs.state[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vdata) == 10
        for vdata in m.fs.state[1].phase_frac.values():
            assert get_scaling_factor(vdata) == 1
        assert get_scaling_factor(m.fs.state[1]._teq["Vap", "Liq"]) == 1 / 300
        assert get_scaling_factor(m.fs.state[1]._t1_Vap_Liq) == 1 / 300
        assert (
            get_scaling_factor(m.fs.state[1].temperature_bubble["Vap", "Liq"])
            == 1 / 300
        )
        for vdata in m.fs.state[1]._mole_frac_tbub.values():
            assert get_scaling_factor(vdata) == 10
        for vdata in m.fs.state[1].log_mole_frac_comp.values():
            assert get_scaling_factor(vdata) == 1
        assert (
            get_scaling_factor(m.fs.state[1].temperature_dew["Vap", "Liq"]) == 1 / 300
        )
        for vdata in m.fs.state[1]._mole_frac_tdew.values():
            assert get_scaling_factor(vdata) == 10
        for vdata in m.fs.state[1].log_mole_frac_tdew.values():
            assert get_scaling_factor(vdata) == 1
        for vdata in m.fs.state[1].log_mole_frac_phase_comp.values():
            assert get_scaling_factor(vdata) == 1

        # Constraints
        assert get_scaling_factor(m.fs.state[1].total_flow_balance) == 1e-2
        for cdata in m.fs.state[1].component_flow_balances.values():
            assert get_scaling_factor(cdata) == 1e-1
        assert get_scaling_factor(m.fs.state[1].sum_mole_frac) == 10
        for cdata in m.fs.state[1].phase_fraction_constraint.values():
            assert get_scaling_factor(cdata) == 1e-2
        assert get_scaling_factor(m.fs.state[1]._t1_constraint_Vap_Liq) == 1 / 300
        for cdata in m.fs.state[1].eq_temperature_bubble.values():
            assert get_scaling_factor(cdata) == 1
        for cdata in m.fs.state[1].log_mole_frac_comp_eqn.values():
            assert get_scaling_factor(cdata) == 10
        for cdata in m.fs.state[1].log_mole_frac_tbub_eqn.values():
            assert get_scaling_factor(cdata) == 10
        assert get_scaling_factor(m.fs.state[1].eq_mole_frac_tbub["Vap", "Liq"]) == 1
        assert get_scaling_factor(m.fs.state[1]._teq_constraint_Vap_Liq) == 1 / 300
        for cdata in m.fs.state[1].eq_temperature_dew.values():
            assert get_scaling_factor(cdata) == 1
        for cdata in m.fs.state[1].log_mole_frac_tdew_eqn.values():
            assert get_scaling_factor(cdata) == 10
        assert get_scaling_factor(m.fs.state[1].eq_mole_frac_tdew["Vap", "Liq"]) == 1
        for cdata in m.fs.state[1].equilibrium_constraint.values():
            assert get_scaling_factor(cdata) == 1
        for cdata in m.fs.state[1].log_mole_frac_phase_comp_eqn.values():
            assert get_scaling_factor(cdata) == 10

        # Expressions
        for edata in m.fs.state[1].flow_mol_phase_comp.values():
            assert get_scaling_factor(edata) == 0.1

    @pytest.mark.integration
    def test_T_sweep(self, m):
        assert_units_consistent(m)

        m.fs.obj = Objective(expr=((m.fs.state[1].temperature - 510) / 100) ** 2)
        m.fs.state[1].temperature.setub(600)

        for P in logspace(4.8, 5.9, 8):
            m.fs.obj.deactivate()

            m.fs.state[1].flow_mol.fix(100)
            m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
            m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
            m.fs.state[1].temperature.fix(300)
            m.fs.state[1].pressure.fix(float(P))

            m.fs.state.initialize()

            m.fs.state[1].temperature.unfix()
            m.fs.obj.activate()
            results = solver.solve(m)
            assert check_optimal_termination(results)

            assert m.fs.state[1].flow_mol_phase["Liq"].value <= 1e-4

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
                results = solver.solve(m)
                assert check_optimal_termination(results)

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
