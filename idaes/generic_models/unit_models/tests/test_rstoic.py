##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Tests for IDAES Stoichiometric reactor.

Author: Chinedu Okoli, Andrew Lee
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

from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)

from idaes.generic_models.unit_models.stoichiometric_reactor import \
    StoichiometricReactor

from idaes.generic_models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock)
from idaes.generic_models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock)

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import (PhysicalParameterTestBlock,
                                     ReactionParameterTestBlock,
                                     initialization_tester)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.reactions = ReactionParameterTestBlock(default={
                            "property_package": m.fs.properties})

    m.fs.unit = StoichiometricReactor(default={
            "property_package": m.fs.properties,
            "reaction_package": m.fs.reactions})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 12

    assert m.fs.unit.config.material_balance_type == \
        MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == \
        EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert not m.fs.unit.config.has_heat_of_reaction
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.reaction_package is m.fs.reactions


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(default={
                                "property_package": m.fs.properties})

        m.fs.unit = StoichiometricReactor(default={
                "property_package": m.fs.properties,
                "reaction_package": m.fs.reactions,
                "has_heat_transfer": True,
                "has_heat_of_reaction": True,
                "has_pressure_change": True})

        m.fs.unit.inlet.flow_vol.fix(1)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(101325.0)

        m.fs.unit.rate_reaction_extent[0, 'R1'].fix(90)
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

        assert hasattr(sapon.fs.unit, "rate_reaction_extent")

        assert number_variables(sapon) == 24
        assert number_total_constraints(sapon) == 13
        assert number_unused_variables(sapon) == 0

    @pytest.mark.component
    def test_units(self, sapon):
        assert_units_consistent(sapon)
        assert_units_equivalent(sapon.fs.unit.heat_duty[0], units.W)
        assert_units_equivalent(sapon.fs.unit.deltaP[0], units.Pa)
        assert_units_equivalent(sapon.fs.unit.rate_reaction_extent[0, "R1"],
                                units.mol/units.s)

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
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert (pytest.approx(101325.0, abs=1e-2) ==
                value(sapon.fs.unit.outlet.pressure[0]))
        assert (pytest.approx(304.21, abs=1e-2) ==
                value(sapon.fs.unit.outlet.temperature[0]))
        assert (pytest.approx(90, abs=1e-2) ==
                value(sapon.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert abs(value(sapon.fs.unit.inlet.flow_vol[0] -
                         sapon.fs.unit.outlet.flow_vol[0])) <= 1e-6
        assert (abs(value(sapon.fs.unit.inlet.flow_vol[0] *
                          sum(sapon.fs.unit.inlet.conc_mol_comp[0, j]
                              for j in sapon.fs.properties.component_list) -
                          sapon.fs.unit.outlet.flow_vol[0] *
                          sum(sapon.fs.unit.outlet.conc_mol_comp[0, j]
                              for j in sapon.fs.properties.component_list)))
                <= 1e-6)

        assert (pytest.approx(4410000, abs=1e3) == value(
                sapon.fs.unit.control_volume.heat_of_reaction[0]))
        assert abs(value(
                (sapon.fs.unit.inlet.flow_vol[0] *
                 sapon.fs.properties.dens_mol *
                 sapon.fs.properties.cp_mol *
                 (sapon.fs.unit.inlet.temperature[0] -
                    sapon.fs.properties.temperature_ref)) -
                (sapon.fs.unit.outlet.flow_vol[0] *
                 sapon.fs.properties.dens_mol *
                 sapon.fs.properties.cp_mol *
                 (sapon.fs.unit.outlet.temperature[0] -
                  sapon.fs.properties.temperature_ref)) +
                sapon.fs.unit.control_volume.heat_of_reaction[0])) <= 1e-3

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, sapon):
        sapon.fs.unit.report()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_costing(self, sapon):
        sapon.fs.unit.get_costing()
        assert isinstance(sapon.fs.unit.costing.purchase_cost, Var)
        sapon.fs.unit.diameter.fix(2)
        sapon.fs.unit.length.fix(3)
        results = solver.solve(sapon)
        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok
        assert (pytest.approx(56327.5803, abs=1e3) ==
                value(sapon.fs.unit.costing.base_cost))
        assert (pytest.approx(85432.06008, abs=1e3) ==
                value(sapon.fs.unit.costing.purchase_cost))

        assert_units_consistent(sapon.fs.unit.costing)
