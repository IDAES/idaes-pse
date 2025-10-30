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
Tests for Membrane 1D model
"""
__author__ = "Maojian Wang"

# pylint: disable=unused-import
import pytest

from pyomo.environ import (
    check_optimal_termination,
    assert_optimal_termination,
    ConcreteModel,
    value,
)
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.solvers import get_solver
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.co2_capture_and_utilization.unit_models import (
    Membrane1D,
    MembraneFlowPattern,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config_countercurrent():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = GenericParameterBlock(
        **get_prop(["CO2", "H2O", "N2"], ["Vap"], eos=EosType.IDEAL),
        doc="Key flue gas property parameters",
    )

    m.fs.unit = Membrane1D(
        finite_elements=3,
        dynamic=False,
        sweep_flow=True,
        flow_type=MembraneFlowPattern.COUNTERCURRENT,
        feed_side={"property_package": m.fs.properties},
        sweep_side={"property_package": m.fs.properties},
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 7
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup


@pytest.mark.unit
def test_config_cocurrent():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = GenericParameterBlock(
        **get_prop(["CO2", "H2O", "N2"], ["Vap"], eos=EosType.IDEAL),
        doc="Key flue gas property parameters",
    )

    m.fs.unit = Membrane1D(
        finite_elements=3,
        dynamic=False,
        sweep_flow=True,
        flow_type=MembraneFlowPattern.COCURRENT,
        feed_side={"property_package": m.fs.properties},
        sweep_side={"property_package": m.fs.properties},
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 7
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup


class TestMembrane:
    @pytest.fixture(scope="class")
    def membrane(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = GenericParameterBlock(
            **get_prop(["CO2", "H2O", "N2"], ["Vap"], eos=EosType.IDEAL),
            doc="Key flue gas property parameters",
        )
        m.fs.unit = Membrane1D(
            finite_elements=3,
            dynamic=False,
            sweep_flow=True,
            flow_type=MembraneFlowPattern.COUNTERCURRENT,
            feed_side={"property_package": m.fs.properties},
            sweep_side={"property_package": m.fs.properties},
        )

        m.fs.unit.permeance[:, :, "CO2"].fix(1500)
        m.fs.unit.permeance[:, :, "H2O"].fix(1500 / 25)
        m.fs.unit.permeance[:, :, "N2"].fix(1500 / 25)
        m.fs.unit.area.fix(100)
        m.fs.unit.length.fix(10)

        m.fs.unit.feed_side_inlet.flow_mol[0].fix(100)
        m.fs.unit.feed_side_inlet.temperature[0].fix(365)
        m.fs.unit.feed_side_inlet.pressure[0].fix(120000)
        m.fs.unit.feed_side_inlet.mole_frac_comp[0, "N2"].fix(0.76)
        m.fs.unit.feed_side_inlet.mole_frac_comp[0, "CO2"].fix(0.13)
        m.fs.unit.feed_side_inlet.mole_frac_comp[0, "H2O"].fix(0.11)

        m.fs.unit.sweep_side_inlet.flow_mol[0].fix(0.01)
        m.fs.unit.sweep_side_inlet.temperature[0].fix(300)
        m.fs.unit.sweep_side_inlet.pressure[0].fix(51325)
        m.fs.unit.sweep_side_inlet.mole_frac_comp[0, "H2O"].fix(0.9986)
        m.fs.unit.sweep_side_inlet.mole_frac_comp[0, "CO2"].fix(0.0003)
        m.fs.unit.sweep_side_inlet.mole_frac_comp[0, "N2"].fix(0.0001)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, membrane):
        assert hasattr(membrane.fs.unit, "feed_side_inlet")
        assert len(membrane.fs.unit.feed_side_inlet.vars) == 4
        assert hasattr(membrane.fs.unit.feed_side_inlet, "flow_mol")
        assert hasattr(membrane.fs.unit.feed_side_inlet, "mole_frac_comp")
        assert hasattr(membrane.fs.unit.feed_side_inlet, "temperature")
        assert hasattr(membrane.fs.unit.feed_side_inlet, "pressure")

        assert hasattr(membrane.fs.unit, "sweep_side_inlet")
        assert len(membrane.fs.unit.sweep_side_inlet.vars) == 4
        assert hasattr(membrane.fs.unit.sweep_side_inlet, "flow_mol")
        assert hasattr(membrane.fs.unit.sweep_side_inlet, "mole_frac_comp")
        assert hasattr(membrane.fs.unit.sweep_side_inlet, "temperature")
        assert hasattr(membrane.fs.unit.sweep_side_inlet, "pressure")

        assert hasattr(membrane.fs.unit, "feed_side_outlet")
        assert len(membrane.fs.unit.feed_side_outlet.vars) == 4
        assert hasattr(membrane.fs.unit.feed_side_outlet, "flow_mol")
        assert hasattr(membrane.fs.unit.feed_side_outlet, "mole_frac_comp")
        assert hasattr(membrane.fs.unit.feed_side_outlet, "temperature")
        assert hasattr(membrane.fs.unit.feed_side_outlet, "pressure")

        assert hasattr(membrane.fs.unit, "sweep_side_outlet")
        assert len(membrane.fs.unit.sweep_side_outlet.vars) == 4
        assert hasattr(membrane.fs.unit.sweep_side_outlet, "flow_mol")
        assert hasattr(membrane.fs.unit.sweep_side_outlet, "mole_frac_comp")
        assert hasattr(membrane.fs.unit.sweep_side_outlet, "temperature")
        assert hasattr(membrane.fs.unit.sweep_side_outlet, "pressure")

        assert hasattr(membrane.fs.unit, "mscontactor")
        assert hasattr(membrane.fs.unit, "permeability_calculation")
        assert hasattr(membrane.fs.unit, "isothermal_constraint")

        assert number_variables(membrane) == 157
        assert number_total_constraints(membrane) == 89
        assert number_unused_variables(membrane) == 28

    @pytest.mark.component
    def test_structural_issues(self, membrane):
        dt = DiagnosticsToolbox(membrane)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, membrane):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(membrane.fs.unit)
        results = solver.solve(membrane)
        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, membrane):
        dt = DiagnosticsToolbox(membrane)
        dt.assert_no_numerical_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_feed_solution(self, membrane):
        assert pytest.approx(99.99, abs=1e-2) == value(
            membrane.fs.unit.feed_side_outlet.flow_mol[0]
        )
        assert pytest.approx(0.1299, abs=1e-4) == value(
            membrane.fs.unit.feed_side_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.11, abs=1e-4) == value(
            membrane.fs.unit.feed_side_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.76, abs=1e-4) == value(
            membrane.fs.unit.feed_side_outlet.mole_frac_comp[0, "N2"]
        )

        assert pytest.approx(365, abs=1e-2) == value(
            membrane.fs.unit.feed_side_outlet.temperature[0]
        )
        assert pytest.approx(120000, abs=1e-2) == value(
            membrane.fs.unit.feed_side_outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_sweep_side_solution(self, membrane):
        assert pytest.approx(0.01006, abs=1e-4) == value(
            membrane.fs.unit.sweep_side_outlet.flow_mol[0]
        )
        assert pytest.approx(0.0070, abs=1e-4) == value(
            membrane.fs.unit.sweep_side_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.99120, abs=1e-4) == value(
            membrane.fs.unit.sweep_side_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.001710, abs=1e-4) == value(
            membrane.fs.unit.sweep_side_outlet.mole_frac_comp[0, "N2"]
        )

        assert pytest.approx(365, abs=1e-2) == value(
            membrane.fs.unit.sweep_side_outlet.temperature[0]
        )
        assert pytest.approx(51325.0, abs=1e-2) == value(
            membrane.fs.unit.sweep_side_outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_enthalpy_balance(self, membrane):

        assert (
            abs(
                value(
                    (
                        membrane.fs.unit.feed_side_inlet.flow_mol[0]
                        * membrane.fs.unit.mscontactor.feed_side_inlet_state[
                            0
                        ].enth_mol_phase["Vap"]
                        + membrane.fs.unit.sweep_side_inlet.flow_mol[0]
                        * membrane.fs.unit.mscontactor.sweep_side_inlet_state[
                            0
                        ].enth_mol_phase["Vap"]
                        - membrane.fs.unit.feed_side_outlet.flow_mol[0]
                        * membrane.fs.unit.mscontactor.feed_side[0, 3].enth_mol_phase[
                            "Vap"
                        ]
                        - membrane.fs.unit.sweep_side_outlet.flow_mol[0]
                        * membrane.fs.unit.mscontactor.sweep_side[0, 1].enth_mol_phase[
                            "Vap"
                        ]
                    )
                )
            )
            <= 1e-6
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_material_balance(self, membrane):

        assert (
            abs(
                value(
                    (
                        membrane.fs.unit.feed_side_inlet.flow_mol[0]
                        + membrane.fs.unit.sweep_side_inlet.flow_mol[0]
                        - membrane.fs.unit.feed_side_outlet.flow_mol[0]
                        - membrane.fs.unit.sweep_side_outlet.flow_mol[0]
                    )
                )
            )
            <= 1e-3
        )
