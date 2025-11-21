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
Tests for feed with flash.
Authors: Andrew Lee, Daison Caballero
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.unittest import assertStructuredAlmostEqual

from idaes.core import FlowsheetBlock, Component
from idaes.models.unit_models.feed_flash import FeedFlash
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties.modular_properties import GenericParameterBlock
from idaes.models.properties.modular_properties.state_definitions import FpTPxpc
from idaes.models.properties.modular_properties.examples.BT_ideal_FpTPxpc import configuration

from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.tests.test_harness import PropertyTestHarness

from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver

from idaes.core.util import DiagnosticsToolbox



# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver("ipopt_v2")


# -----------------------------------------------------------------------------

def _as_quantity(x):
    unit = pyunits.get_units(x)
    if unit is None:
        unit = pyunits.dimensionless
    return value(x) * unit._get_pint_unit()


class TestBTIdeal(PropertyTestHarness):
    def configure(self):
        self.prop_pack = GenericParameterBlock
        self.param_args = configuration
        self.prop_args = {}
        self.has_density_terms = True


class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ["benzene", "toluene"]
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 4
        for i in model.params._phase_component_set:
            assert i in [
                ("Liq", "benzene"),
                ("Liq", "toluene"),
                ("Vap", "benzene"),
                ("Vap", "toluene"),
            ]

        assert model.params.config.state_definition == FpTPxpc

        assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            {
                "flow_mol_phase": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 450, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 2
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2"]

        assert model.params.phase_equilibrium_list == {
            "PE1": ["benzene", ("Vap", "Liq")],
            "PE2": ["toluene", ("Vap", "Liq")],
        }

        assert model.params.pressure_ref.value == 1e5
        assert model.params.temperature_ref.value == 300

        assert model.params.benzene.mw.value == 78.1136e-3
        assert model.params.benzene.pressure_crit.value == 48.9e5
        assert model.params.benzene.temperature_crit.value == 562.2

        assert model.params.toluene.mw.value == 92.1405e-3
        assert model.params.toluene.pressure_crit.value == 41e5
        assert model.params.toluene.temperature_crit.value == 591.8

        assert_units_consistent(model)

# -----------------------------------------------------------------------------
class TestBTIdeal_FpTPxpc(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = GenericParameterBlock(**configuration)

        m.fs.unit = FeedFlash(property_package=m.fs.properties)

        m.fs.unit.temperature.fix(368.0)
        m.fs.unit.pressure.fix(101325.0)

        m.fs.unit.flow_mol_phase[0, "Liq"].fix(0.5)
        m.fs.unit.flow_mol_phase[0, "Vap"].fix(0.5)
        m.fs.unit.mole_frac_phase_comp[0, "Liq", "benzene"].fix(0.5)
        m.fs.unit.mole_frac_phase_comp[0, "Liq", "toluene"].fix(0.5)
        m.fs.unit.mole_frac_phase_comp[0, "Vap", "benzene"].fix(0.5)
        m.fs.unit.mole_frac_phase_comp[0, "Vap", "toluene"].fix(0.5)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):

        # Check state variable values and bounds
        assert isinstance(model.fs.unit.flow_mol_phase, Var)
        for p in model.fs.unit.flow_mol_phase:
            assert value(model.fs.unit.flow_mol_phase[p]) == 0.5
            assert model.fs.unit.flow_mol_phase[p].ub == 1000
            assert model.fs.unit.flow_mol_phase[p].lb == 0

        assert isinstance(model.fs.unit.pressure, Var)
        assert value(model.fs.unit.pressure[0]) == 101325
        assert model.fs.unit.pressure[0].ub == 1e6
        assert model.fs.unit.pressure[0].lb == 5e4

        assert isinstance(model.fs.unit.temperature, Var)
        assert value(model.fs.unit.temperature[0]) == 368
        assert model.fs.unit.temperature[0].ub == 450
        assert model.fs.unit.temperature[0].lb == 273.15

        assert isinstance(model.fs.unit.mole_frac_phase_comp, Var)
        for (j, p, i) in model.fs.unit.mole_frac_phase_comp:
            assert value(model.fs.unit.mole_frac_phase_comp[j, p, i]) == 0.5

        assert_units_consistent(model)

    @pytest.mark.component
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_structural_warnings()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, model):
        perf_dict = model.fs.unit._get_performance_contents()

        assert perf_dict is None

    # TODO the formatting for modular properties is broken (see #1684).
    # This test can be fixed once that issue is fixed
    # @pytest.mark.ui
    # @pytest.mark.unit
    # def test_get_stream_table_contents(self, model):
    #     stable = model.fs.unit._get_stream_table_contents()

        # expected = {
        #     "Units": {
        #         "flow_mol_phase": getattr(pyunits.pint_registry, "mole/second"),
        #         "mole_frac_phase_comp Liq benzene": getattr(
        #             pyunits.pint_registry, "dimensionless"
        #         ),
        #         "mole_frac_phase_comp Liq toluene": getattr(
        #             pyunits.pint_registry, "dimensionless"
        #         ),
        #         "temperature": getattr(pyunits.pint_registry, "K"),
        #         "pressure": getattr(pyunits.pint_registry, "Pa"),
        #     },
        #     "Outlet": {
        #         "flow_mol": pytest.approx(1.0, rel=1e-4),
        #         "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
        #         "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
        #         "temperature": pytest.approx(368.0, rel=1e-4),
        #         "pressure": pytest.approx(101325.0, rel=1e-4),
        #     },
        # }
        #
        # assert stable.to_dict() == expected

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(101325.0, abs=1e-3) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(368.0, abs=1e-3) == value(
            model.fs.unit.outlet.temperature[0]
        )
        
        assert pytest.approx(1.0, abs=1e-3) == value(
            model.fs.unit.control_volume.properties_out[0].flow_mol
        )

        assert pytest.approx(0.603, abs=1e-3) == value(
            model.fs.unit.outlet.flow_mol_phase[0, "Liq"]
        )
        assert pytest.approx(0.396, abs=1e-3) == value(
            model.fs.unit.outlet.flow_mol_phase[0, "Vap"]
        )

        assert pytest.approx(0.412, abs=1e-3) == value(
            model.fs.unit.outlet.mole_frac_phase_comp[0,
                "Liq", "benzene"
            ]
        )
        assert pytest.approx(0.588, abs=1e-3) == value(
            model.fs.unit.outlet.mole_frac_phase_comp[0,
                "Liq", "toluene"
            ]
        )
        assert pytest.approx(0.634, abs=1e-3) == value(
            model.fs.unit.outlet.mole_frac_phase_comp[0,
                "Vap", "benzene"
            ]
        )
        assert pytest.approx(0.366, abs=1e-3) == value(
            model.fs.unit.outlet.mole_frac_phase_comp[0,
                "Vap", "toluene"
            ]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()