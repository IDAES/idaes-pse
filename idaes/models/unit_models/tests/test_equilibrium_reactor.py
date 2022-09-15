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
Tests for equilibrium reactor unit model.
Authors: Andrew Lee
"""

import pytest
from io import StringIO

from pyomo.environ import check_optimal_termination, ConcreteModel, value, units
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models.unit_models.equilibrium_reactor import EquilibriumReactor
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
    initialization_tester,
)
from idaes.core.solvers import get_solver
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from idaes.core.util.exceptions import InitializationError


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.reactions = ReactionParameterTestBlock(property_package=m.fs.properties)

    m.fs.unit = EquilibriumReactor(
        property_package=m.fs.properties, reaction_package=m.fs.reactions
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 15

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_rate_reactions
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.has_equilibrium_reactions
    assert not m.fs.unit.config.has_phase_equilibrium
    assert not m.fs.unit.config.has_heat_of_reaction
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.reaction_package is m.fs.reactions


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = EquilibriumReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol.fix(1.0e-03)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(101325.0)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

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

        assert hasattr(sapon.fs.unit, "outlet")
        assert len(sapon.fs.unit.outlet.vars) == 4
        assert hasattr(sapon.fs.unit.outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.outlet, "temperature")
        assert hasattr(sapon.fs.unit.outlet, "pressure")

        assert hasattr(sapon.fs.unit, "rate_reaction_constraint")
        assert hasattr(sapon.fs.unit, "heat_duty")
        assert hasattr(sapon.fs.unit, "deltaP")

        assert number_variables(sapon) == 26
        assert number_total_constraints(sapon) == 16
        assert number_unused_variables(sapon) == 0

    @pytest.mark.component
    def test_units(self, sapon):
        assert_units_consistent(sapon)
        assert_units_equivalent(sapon.fs.unit.heat_duty[0], units.W)
        assert_units_equivalent(sapon.fs.unit.deltaP[0], units.Pa)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(101325.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(304.32, abs=1e-2) == value(
            sapon.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(0.01, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0] - sapon.fs.unit.outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.inlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                    - sapon.fs.unit.outlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.outlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0]
                    * sapon.fs.properties.cp_mol
                    * (
                        sapon.fs.unit.inlet.temperature[0]
                        - sapon.fs.properties.temperature_ref
                    )
                    - sapon.fs.unit.outlet.flow_vol[0]
                    * sapon.fs.properties.cp_mol
                    * (
                        sapon.fs.unit.outlet.temperature[0]
                        - sapon.fs.properties.temperature_ref
                    )
                    + sapon.fs.unit.heat_duty[0]
                )
            )
            <= 1e-1
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Heat Duty": sapon.fs.unit.heat_duty[0],
                "Pressure Change": sapon.fs.unit.deltaP[0],
            }
        }

    @pytest.mark.component
    def test_initialization_error(self, sapon):
        sapon.fs.unit.outlet.pressure[0].fix(1)

        with pytest.raises(InitializationError):
            sapon.fs.unit.initialize()
